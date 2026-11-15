//! Gene-set modes of probuild 2.52 built on the gene parser object
//! ([`GeneParser`], models.md §0.5): models.md §5-10.
//!
//! - `--parse` (0x4100ec): writes the coding sequences and the non-coding
//!   segments of a genome + gene list as two sequence files.
//! - `--mkchain` (0x4109c2): trains coding and non-coding chains from a
//!   genome + gene list and writes their tables to two files.
//! - `--intergenic` (0x40fb32): writes the sequence outside a list of
//!   `a b` regions.
//! - `--pre_start` (0x40eda4): writes the pre-start sequences (upstream of
//!   the start codons) as FASTA, the input of the Gibbs3 motif finder.
//! - `--STAT_NON_VS_ALL` (0x40f654): non-coding fraction and GC of the
//!   non-coding part vs the whole sequence.
//! - `--GENE_SET_STAT` (0x409d16): counts of neighbouring gene pairs by
//!   strand arrangement and overlap.
//!
//! Except for --GENE_SET_STAT, every handler constructs the gene parser
//! first (as main does in the original), so the parser's option messages
//! come before the handler's own.

use super::geneset::{self, Gene};
use super::gp::{out_endl, seqstr, GeneParser, SeqObj};
use super::io::{self, OFile};
use super::markov::Chain;
use super::params::{atoi, Params};
use super::pm::{fixed, g, P};

/// Writes `text` to a new file (truncate; an open failure is silent).
fn write(name: &[u8], text: &str) {
    let mut f = OFile::open(name, false);
    f.w(text.as_bytes());
    f.close();
}

/// String check 0x41fd58 of `k` (prints the value error on failure),
/// followed by the handler's own message `msg + k` to stderr.
fn need(p: &mut P, k: &str, msg: &str) -> bool {
    if !p.chk_str(k) {
        io::err(format!("{msg}{k}\n").as_bytes());
        return false;
    }
    true
}

/// The `--parse` coding file (0x40d4d4): `>cod1`, the set-1 sequences one
/// per line, then `>cod2` and set 2 if set 2 is not empty.
fn cod_text(gp: &GeneParser) -> String {
    let mut o = String::from(">cod1\n");
    for x in &gp.cod1 {
        o += &seqstr(x);
    }
    if !gp.cod2.is_empty() {
        o += ">cod2\n";
        for x in &gp.cod2 {
            o += &seqstr(x);
        }
    }
    o
}

/// The non-coding file of --parse / --intergenic (0x40eed8): one segment
/// per line, no header.
fn non_text(gp: &GeneParser) -> String {
    gp.non.iter().map(|x| seqstr(x)).collect()
}

/// `--parse` (0x4100ec), models.md §5.
///
/// Checks `--cod` and `--non` (output names), reads `--COD_L_MARGIN` and
/// the signed `--COD_R_MARGIN`, then the sequence and gene list, extracts
/// the coding set with these margins and the non-coding segments, and
/// writes both files.
pub fn parse(p: &Params) -> bool {
    let mut p = P::new(p);
    let p = &mut p;
    let mut gp = GeneParser::new(p);
    if !p.has("--cod") {
        io::err(b"error, key missing --cod\n");
        return false;
    }
    if !need(p, "--cod", "error, file name missing ") {
        return false;
    }
    // The same check twice: the second message is unreachable (kept from
    // the original).  A missing --non is inserted by the check, so it
    // reports as a value error `...--non: ` first.
    if !need(p, "--non", "error, key missing ") || !need(p, "--non", "error, file name missing ") {
        return false;
    }
    let (ok, lm) = p.int_opt("--COD_L_MARGIN", 0, false);
    if !ok {
        return false;
    }
    let (ok, rm) = p.int_opt("--COD_R_MARGIN", 0, true);
    if !ok {
        return false;
    }
    let mut so = SeqObj::new();
    if !so.load(p) {
        return false;
    }
    let Some(genes) = geneset::load_param(p, so.seq.len() as i32) else { return false };
    gp.collect_coding(&so, &genes, lm, rm);
    gp.collect_noncoding(&so, &genes);
    write(&p.val("--cod"), &cod_text(&gp));
    write(&p.val("--non"), &non_text(&gp));
    true
}

/// `--mkchain` (0x4109c2), models.md §6.
///
/// Like the gene-set path of --mkmod (coding set without gene margins,
/// chain margins from MKCHAIN_*, non-coding chain with L = R = 0, no
/// RANDOM_NONCOD) but the order comes from `--order` and the tables go to
/// `--chain_cod` (with `COD1` / `COD2` headers if SEPARATE_COD) and
/// `--chain_non`.  The spelling `erorr` is the original's.
pub fn mkchain(p: &Params) -> bool {
    let mut p = P::new(p);
    let p = &mut p;
    let mut gp = GeneParser::new(p);
    let e = |m: &str| io::err(m.as_bytes());
    if !p.has("--chain_cod") || !p.chk_str("--chain_cod") {
        e("error, key missing --chain_cod\n");
        return false;
    }
    if !p.chk_str("--chain_non") {
        e("error, key missing --chain_non\n");
        return false;
    }
    // unreachable second check, kept from the original
    if !p.chk_str("--chain_non") {
        e("erorr, file name missing --chain_non\n");
        return false;
    }
    if !p.has("--order") {
        e("error, key missing --order\n");
        return false;
    }
    if !p.chk_int("--order") {
        e("erorr, key value missing --order\n");
        return false;
    }
    let order = atoi(&p.val("--order"));
    let mut so = SeqObj::new();
    if !so.load(p) {
        return false;
    }
    let Some(genes) = geneset::load_param(p, so.seq.len() as i32) else { return false };
    gp.collect_coding(&so, &genes, 0, 0);
    gp.collect_noncoding(&so, &genes);
    let mut cod1_chain = Chain::from_params(p, order, 3);
    let mut cod2_chain = Chain::from_params(p, order, 3);
    let (l, r) = (cod1_chain.l_margin, cod1_chain.r_margin);
    cod1_chain.train_set(&gp.cod1, l, r);
    if gp.separate_cod {
        let (l, r) = (cod2_chain.l_margin, cod2_chain.r_margin);
        cod2_chain.train_set(&gp.cod2, l, r);
    }
    let mut non_chain = Chain::from_params(p, order, 1);
    non_chain.l_margin = 0;
    non_chain.r_margin = 0;
    non_chain.train_set(&gp.non, 0, 0);
    let cod_out = if !gp.separate_cod {
        cod1_chain.table()
    } else {
        format!("COD1\n{}COD2\n{}", cod1_chain.table(), cod2_chain.table())
    };
    write(&p.val("--chain_cod"), &cod_out);
    write(&p.val("--chain_non"), &non_chain.table());
    true
}

/// `--intergenic` (0x40fb32), models.md §7.
///
/// The `--geneset` file is read in the pair format (0x409c4a): lines
/// `a b`, 1-based inclusive, shifted by `--shift_zero`.  The sequence
/// outside all listed regions is split into segments of at least MIN_NON
/// bases and written to `--non`.  `--COD_L_MARGIN` / `--COD_R_MARGIN` are
/// read (a bad value fails) but not used.
pub fn intergenic(p: &Params) -> bool {
    let mut p = P::new(p);
    let p = &mut p;
    let mut gp = GeneParser::new(p);
    if !need(p, "--non", "error, key missing ") || !need(p, "--non", "error, file name missing ") {
        return false;
    }
    if !p.int_opt("--COD_L_MARGIN", 0, false).0 || !p.int_opt("--COD_R_MARGIN", 0, true).0 {
        return false;
    }
    let mut so = SeqObj::new();
    if !so.load(p) {
        return false;
    }
    // 0x409c4a: pair format, --shift_zero offset
    let (ok, offset) = p.int_opt("--shift_zero", 0, true);
    if !ok {
        return false;
    }
    let name = p.str_req("--geneset").unwrap_or_default();
    let Some(d) = io::load_file(&name) else { return false };
    let genes = geneset::parse(&d, so.seq.len() as i32, true, offset);
    gp.collect_noncoding(&so, &genes);
    write(&p.val("--non"), &non_text(&gp));
    true
}

/// Pre-start output (0x40e4be, writer 0x40ebd2), models.md §8.  Used by
/// `--pre_start` and by `--mkmod` with `--pre_start` but without
/// `--fixmotif`.
///
/// Options: `--pre_start` (output file), `--PRE_START_WIDTH` (total
/// length of a sequence), `--PRE_START_RELATIVE_POS` (rel: how many of
/// these bases lie downstream, from the first base of the start codon on),
/// `--PRE_START_INTERGENIC`, `--PRE_START_LINES` (maximum number of
/// records, negative = all).  The first failure returns `false`.  The
/// sequences have `WIDTH - rel` upstream and `rel` downstream bases.
/// (Inside --mkmod's RBS block, by contrast, PRE_START_WIDTH is the
/// upstream part and `--width` the total.)
///
/// Output: `>i L R s` (1-based gene ends) and the sequence, for i = 1 ..
/// LINES.  Not reproduced: the original's `reserve(LINES)` aborts for a
/// huge LINES value.
pub fn pre_start_out(p: &mut P, gp: &mut GeneParser, so: &SeqObj, genes: &[Gene]) -> bool {
    if p.str_opt("--pre_start", Vec::new()).is_none() {
        return false;
    }
    let (ok, width) = p.int_opt("--PRE_START_WIDTH", 0, false);
    if !ok {
        return false;
    }
    let (ok, rel) = p.int_opt("--PRE_START_RELATIVE_POS", 0, false);
    if !ok {
        return false;
    }
    let (ok, intergenic) = p.bool_opt("--PRE_START_INTERGENIC", true);
    if !ok {
        return false;
    }
    let (ok, lines) = p.int_opt("--PRE_START_LINES", -1, true);
    if !ok {
        return false;
    }
    let names = gp.collect_pre_start(so, genes, so.seq.len() as i32, width.wrapping_sub(rel), rel, intergenic);
    let fname = p.val("--pre_start");
    // negative LINES: one more than there are sequences, i.e. no limit
    let limit: i64 = if lines < 0 { gp.pre_start.len() as i64 + 1 } else { lines as i64 };
    let mut o = String::new();
    let mut i: i64 = 1;
    while ((i - 1) as usize) < gp.pre_start.len() && i <= limit {
        o += &format!(">{}{}{}", i, names[(i - 1) as usize], seqstr(&gp.pre_start[(i - 1) as usize]));
        i += 1;
    }
    write(&fname, &o);
    true
}

/// `--pre_start` (0x40eda4): sequence, gene list (maxpos = length), then
/// [`pre_start_out`].
pub fn pre_start(p: &Params) -> bool {
    let mut p = P::new(p);
    let p = &mut p;
    let mut gp = GeneParser::new(p);
    let mut so = SeqObj::new();
    if !so.load(p) {
        return false;
    }
    let Some(genes) = geneset::load_param(p, so.seq.len() as i32) else { return false };
    pre_start_out(p, &mut gp, &so, &genes)
}

/// `--STAT_NON_VS_ALL` (0x40f654), models.md §9.
///
/// Prints `frac  gc_all gc_non` + endl to stdout: `frac` = total length of
/// the non-coding segments (MIN_NON applies) / sequence length (`%g`),
/// `gc_all` = GC percentage of the whole sequence (`%.1f`, 0x43d29e),
/// `gc_non` = GC percentage of the segments (`%g`).  Note the two spaces
/// after `frac`.  GC percentages use int counts: `(C + G) * 100.0 /
/// (T + A + C + G)`, 0.0 without bases.  An empty sequence fails with
/// `error: sequence size is <0>`.
pub fn stat_non_vs_all(p: &Params) -> bool {
    let mut p = P::new(p);
    let p = &mut p;
    let mut gp = GeneParser::new(p);
    let mut so = SeqObj::new();
    if !so.load(p) {
        return false;
    }
    let seq_len = so.seq.len() as i32;
    let Some(genes) = geneset::load_param(p, seq_len) else { return false };
    gp.collect_noncoding(&so, &genes);
    let mut non_counts = [0i32; 256];
    let mut non_total: i32 = 0;
    for x in &gp.non {
        non_total = non_total.wrapping_add(x.len() as i32);
        for &v in x {
            non_counts[v as usize] += 1;
        }
    }
    let gc_percent = |c: &[i32; 256]| -> f64 {
        let num = c[2].wrapping_add(c[1]) as f64 * 100.0;
        let t = c[3].wrapping_add(c[0]).wrapping_add(c[1]).wrapping_add(c[2]);
        if t == 0 {
            0.0
        } else {
            num / t as f64
        }
    };
    let gc_non = gc_percent(&non_counts);
    if seq_len == 0 {
        io::err(b"error: sequence size is <0>\n");
        return false;
    }
    let frac = non_total as f64 / seq_len as f64;
    let mut all_counts = [0i32; 256];
    for &v in &so.seq {
        all_counts[v as usize] += 1;
    }
    out_endl(&format!("{}  {} {}", g(frac), fixed(1, gc_percent(&all_counts), false), g(gc_non)));
    true
}

/// `--GENE_SET_STAT` (0x409afe + 0x409d16), models.md §10.
///
/// For every pair of consecutive list entries (a, b), including header
/// lines: if b starts or ends before a, print `a.left b.left a.right
/// b.right check input data` (0-based); otherwise classify the pair by the
/// strands and by whether they are separate (`a.right < b.left`) or
/// overlap.  Entries without strand print `no direction set for gene`.
///
/// Final line `same_sep same_ovl +-_sep +-_ovl -+_sep -+_ovl ` (trailing
/// space), where same = same strand (both '+' and both '-' summed), +- =
/// '+' followed by '-' (tail to tail), -+ = '-' followed by '+' (head to
/// head).
pub fn gene_set_stat(p: &Params) -> bool {
    let mut p = P::new(p);
    let p = &mut p;
    let Some(gs) = geneset::load_param(p, -1) else { return false };
    // pair_counts[0..4]: same strand  ++ sep, ++ overlap, -- sep, -- overlap
    // pair_counts[4..8]: other strand +- sep, +- overlap, -+ sep, -+ overlap
    let mut pair_counts = [0i32; 8];
    let mut i = 0;
    while i < gs.len() && i + 1 < gs.len() {
        let (a, b) = (&gs[i], &gs[i + 1]);
        if a.left > b.left || a.right > b.right {
            out_endl(&format!("{} {} {} {} check input data", a.left, b.left, a.right, b.right));
        } else {
            let base = if a.strand == b.strand { 0 } else { 4 };
            match a.strand {
                1 => pair_counts[base + if a.right < b.left { 0 } else { 1 }] += 1,
                2 => pair_counts[base + if a.right < b.left { 2 } else { 3 }] += 1,
                _ => out_endl("no direction set for gene"),
            }
        }
        i += 1;
    }
    let k = &pair_counts;
    out_endl(&format!("{} {} {} {} {} {} ", k[2] + k[0], k[1] + k[3], k[4], k[5], k[6], k[7]));
    true
}
