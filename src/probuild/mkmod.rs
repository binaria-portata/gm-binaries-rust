//! Model file writers of probuild 2.52: `--mkmod_euk` (0x440730),
//! `--mkmod_non` (0x442010), `--mkmod` (0x442432) and `--combine`
//! (0x442c58).  models.md §2-4, §15a; core.md section 9.
//!
//! - `--mkmod_euk` (GeneMark-ES training, gmes_petap.pl): a period-3 coding
//!   chain and a period-1 non-coding chain trained on two plain sequence
//!   files, written one word per line as
//!   `WORD cod0 cod1 cod2 non non` (the non-coding column twice: gmhmme3
//!   uses the second copy for introns / intergenic).  build_mod.pl stores
//!   the file verbatim as the `$MARKOV` value of the gmhmme3 model.
//! - `--mkmod_non` (GeneMark-ES `--fungus`): only the non-coding chain,
//!   `WORD v` per line, trained on `spacer.seq` (stored by build_mod.pl as
//!   `$MARKOV_BP_SPACER`).
//! - `--mkmod` (GeneMarkS-T, gmst.pl): a complete gmhmmp model.  The
//!   training sets come either from a genome and gene list (`--seq`,
//!   `--geneset`: coding sets, non-coding segments, start codon
//!   frequencies, optionally an RBS matrix) or from two plain sequence
//!   files (`--cod`, `--non`).  The scalar model keys (PHMM, NAME, ORDM,
//!   stop codons, length limits, ...) are copied verbatim from the
//!   parameters (normally the `--par` file); COD1, COD2 and NONC are the
//!   chain tables; ATG_ .. CTG_ the start probabilities; RBSS/RBSM/RBSL/RBSD
//!   the RBS motif block.
//! - `--combine`: merges two model files (COD1 from the first, COD2 from
//!   the second, NONC from either).
//!
//! All model files go through the report object [`Report`]: a map of
//! `key -> (label, value)` plus a request list naming which keys to write
//! in which order.

use super::geneset;
use super::gp::{GeneParser, SeqObj};
use super::io::{self, OFile};
use super::markov::{read_stream, word, Chain};
use super::params::{Bytes, Params};
use super::pm::P;
use std::collections::BTreeMap;

/// Report object (0x405a08, models.md §0.4).
///
/// The original's append flag (+0x10) is passed to [`Report::write`]
/// instead of being stored.
#[derive(Default)]
pub struct Report {
    /// +0x18: `std::map<key, (label, value)>`, byte-wise key order.
    pub entries: BTreeMap<Bytes, (Bytes, Bytes)>,
    /// +0: requested keys in output order (a key may appear twice and is
    /// then written twice).  Empty = write every entry in key order.
    pub requests: Vec<Bytes>,
}

impl Report {
    /// Sets an entry (0x420f68).
    pub fn put(&mut self, key: &str, label: &[u8], value: &[u8]) {
        self.entries.insert(key.as_bytes().to_vec(), (label.to_vec(), value.to_vec()));
    }

    /// Writes `label + value + "\n"` per entry (0x421806 truncating /
    /// 0x421782 with `append`).
    ///
    /// Without requests every entry is written in key order (0x421290).
    /// With requests (0x4210ee, 0x4214e0) all of them are checked first:
    /// an empty key, a key without entry or an entry with an empty value
    /// prints its error to stderr and nothing is written (`false`).  An
    /// empty file name writes to stdout; a file that cannot be opened is
    /// silently ignored (still `true`).
    pub fn write(&self, name: &[u8], append: bool) -> bool {
        let mut text = Vec::new();
        if self.requests.is_empty() {
            for (l, v) in self.entries.values() {
                text.extend_from_slice(l);
                text.extend_from_slice(v);
                text.push(b'\n');
            }
        } else {
            for k in &self.requests {
                let ks = String::from_utf8_lossy(k);
                if k.is_empty() {
                    io::err(b"error: no value in request\n");
                    return false;
                }
                match self.entries.get(k) {
                    None => {
                        io::err(format!("error: no entry to report, {ks}\n").as_bytes());
                        return false;
                    }
                    Some((_, v)) if v.is_empty() => {
                        io::err(format!("error: no value to report, {ks}\n").as_bytes());
                        return false;
                    }
                    _ => {}
                }
            }
            for k in &self.requests {
                let (l, v) = &self.entries[k];
                text.extend_from_slice(l);
                text.extend_from_slice(v);
                text.push(b'\n');
            }
        }
        if name.is_empty() {
            io::out(&text);
            io::flush_stdout();
        } else {
            let mut f = OFile::open(name, append);
            f.w(&text);
            f.close();
        }
        true
    }
}

/// `--mkmod_euk` (0x441fd8 -> 0x440730), models.md §2.
///
/// Reads `--ORDM`, `--cod`, `--non` (all required) and both streams; the
/// GeneMark-ES files have one sequence per line without deflines, so all
/// lines are concatenated into one stream (the padding `n`s break the
/// windows at the joins).  Trains a coding chain (ORDM, period 3) and a
/// non-coding chain (ORDM, period 1, `--order_non`, `--revcomp_non`) on
/// the whole streams, then writes `--mkmod_euk`.  Coding phase j counts
/// windows whose last base is at stream position `ORDM + j (mod 3)`.
pub fn mkmod_euk(p: &Params) -> bool {
    let mut p = P::new(p);
    let p = &mut p;
    let order = match p.int_req("--ORDM") {
        Some(o) => o,
        None => return false,
    };
    let Some(cod_name) = p.str_req("--cod") else { return false };
    let Some(non_name) = p.str_req("--non") else { return false };
    let Some(cod_seq) = read_stream(&cod_name) else { return false };
    let Some(non_seq) = read_stream(&non_name) else { return false };
    // The MKCHAIN_*_MARGING options are read by the constructor (and their
    // errors printed) but overwritten with 0 here; stream training ignores
    // them anyway.
    let mut cod_chain = Chain::from_params(p, order, 3);
    cod_chain.l_margin = 0;
    cod_chain.r_margin = 0;
    cod_chain.train_stream(&cod_seq);
    let mut non_chain = Chain::from_params(p, order, 1);
    non_chain.l_margin = 0;
    non_chain.r_margin = 0;
    non_chain.train_stream(&non_seq);
    let Some(out) = p.str_req("--mkmod_euk") else { return false };
    let mut f = OFile::open(&out, false);
    for i in 0..cod_chain.size() {
        // non-coding first, as in the original (only matters for which
        // chain reports a bad precision; both have the same one)
        let non_col = non_chain.word_values(i);
        let cod_cols = cod_chain.word_values(i);
        f.w(format!("{} {} {} {}\n", word(i, order), cod_cols, non_col, non_col).as_bytes());
    }
    f.close();
    true
}

/// `--mkmod_non` (0x442010), models.md §3: one non-coding chain (ORDM,
/// period 1, `--order_non`, `--revcomp_non`) from `--non`, written as
/// `WORD v` per line to the `--mkmod_non` file.
pub fn mkmod_non(p: &Params) -> bool {
    let mut p = P::new(p);
    let p = &mut p;
    let order = match p.int_req("--ORDM") {
        Some(o) => o,
        None => return false,
    };
    let Some(non_name) = p.str_req("--non") else { return false };
    let Some(non_seq) = read_stream(&non_name) else { return false };
    let mut non_chain = Chain::from_params(p, order, 1);
    non_chain.l_margin = 0;
    non_chain.r_margin = 0;
    non_chain.train_stream(&non_seq);
    let Some(out) = p.str_req("--mkmod_non") else { return false };
    let mut f = OFile::open(&out, false);
    for i in 0..non_chain.size() {
        f.w(format!("{} {}\n", word(i, order), non_chain.word_values(i)).as_bytes());
    }
    f.close();
    true
}

/// The keys of a GeneMarkS-T model in file order (request list of the
/// model writer 0x445560).
const MODEL_KEYS: &[&str] = &[
    "--PHMM", "--NAME", "--ORDM", "--ATG_", "--GTG_", "--TTG_", "--CTG_", "--TAA_", "--TAG_",
    "--TGA_", "--MINC", "--MAXC", "--MAXN", "--NDEC", "--CDEC", "--CDCD", "--CD1P", "--CD2P",
    "--COD1", "--COD2", "--NONC",
];

/// Model key validator (0x444850, models.md §4.1).  Checks the scalar model
/// keys in the original order and stops at the first failure (which prints
/// its value error).  The checks insert missing keys with an empty value
/// (see [`P::val`]); an empty value passes the digit and float checks but
/// not the string check, so without `--par` the first message is
/// `error in value of parameter --PHMM: `.  COD1/COD2/NONC only need to be
/// present (their values are never used); missing ones fail WITHOUT a
/// message.
pub fn validate(p: &mut P) -> bool {
    for k in ["--PHMM", "--NAME"] {
        if !p.chk_str(k) {
            return false;
        }
    }
    if !p.chk_int("--ORDM") {
        return false;
    }
    for k in ["--ATG_", "--GTG_", "--TTG_", "--CTG_", "--TAA_", "--TAG_", "--TGA_"] {
        if !p.chk_float_nonneg(k) {
            return false;
        }
    }
    for k in ["--MINC", "--MAXC", "--MAXN"] {
        if !p.chk_int(k) {
            return false;
        }
    }
    for k in ["--NDEC", "--CDEC", "--CDCD", "--CD1P", "--CD2P"] {
        if !p.chk_float_nonneg(k) {
            return false;
        }
    }
    ["--COD1", "--COD2", "--NONC"].iter().all(|k| p.has(k))
}

/// Model writer (0x445560, models.md §4.4); `p` is the parameter map whose
/// scalar keys are copied (the command line for --mkmod, the first model
/// for --combine).
///
/// Puts the raw parameter strings of PHMM .. CD2P (so `ORDM 4` keeps the
/// text as given), requests all [`MODEL_KEYS`], then RBSS/RBSM/RBSL/RBSD
/// for each of them that is a PARAMETER, then all four again if
/// `--fixmotif` is present.  As in the original: `--RBSS 1 --fixmotif` writes the
/// RBSS line twice.  Written truncating to `name`.
pub fn write_model(p: &mut P, rep: &mut Report, name: &[u8]) -> bool {
    for k in [
        "--PHMM", "--NAME", "--ORDM", "--TAA_", "--TAG_", "--TGA_", "--MINC", "--MAXC", "--MAXN",
        "--NDEC", "--CDEC", "--CDCD", "--CD1P", "--CD2P",
    ] {
        let v = p.val(k);
        // label = key without "--" plus a space, e.g. "PHMM "
        rep.put(k, format!("{} ", &k[2..]).as_bytes(), &v);
    }
    rep.requests = MODEL_KEYS.iter().map(|k| k.as_bytes().to_vec()).collect();
    for k in ["--RBSS", "--RBSM", "--RBSL", "--RBSD"] {
        if p.has(k) {
            rep.requests.push(k.as_bytes().to_vec());
        }
    }
    if p.has("--fixmotif") {
        for k in ["--RBSS", "--RBSM", "--RBSL", "--RBSD"] {
            rep.requests.push(k.as_bytes().to_vec());
        }
    }
    rep.write(name, false)
}

/// The `--RANDOM_NONCOD` test of 0x43dbc6 / 0x440dba: present, passes the
/// float-character check (0x41fd10, which prints its error twice and then
/// counts as absent) and `strtod(value) != 0`.
fn random_noncod(p: &mut P) -> bool {
    if !p.has("--RANDOM_NONCOD") {
        return false;
    }
    if !p.chk_float2("--RANDOM_NONCOD") {
        return false;
    }
    p.bool_val("--RANDOM_NONCOD")
}

/// `--RANDOM_NONCOD_GC` through the double getter 0x41fff2.  A value that
/// fails the character check is reported by the getter and leaves the
/// original's local uninitialised (a stack slot that holds 0.0 on both
/// paths in practice, verified with the original).  The port uses 0.0, which
/// reproduces that output, and warns.
fn random_noncod_gc(p: &mut P) -> f64 {
    let (ok, gc) = p.f64_opt("--RANDOM_NONCOD_GC", 0.0);
    if !ok {
        io::err(
            b"probuild-rs: warning: --RANDOM_NONCOD_GC is not a number; using a GC fraction of 0 \
              (not handled by the original)\n",
        );
    }
    gc
}

/// Puts the ATG_ GTG_ TTG_ CTG_ entries.  With a gene parser (gene-set
/// path) an enabled codon (`strtod(value) != 0`) gets its trained
/// probability, `%#.Pf` with START_PRECISION (e.g. `1.000`); otherwise the
/// raw parameter string is copied (e.g. `0`).
fn start_entries(p: &mut P, rep: &mut Report, gp: Option<&GeneParser>) {
    for (k, codon) in [("--ATG_", 14), ("--GTG_", 46), ("--TTG_", 62), ("--CTG_", 30)] {
        let label = format!("{} ", &k[2..]);
        let v = match gp {
            Some(gp) if p.bool_val(k) => gp.start_prob_str(codon).into_bytes(),
            _ => p.val(k),
        };
        rep.put(k, label.as_bytes(), &v);
    }
}

/// --mkmod with `--geneset` (0x43dbc6, models.md §4.2): trains everything
/// from a genome and a gene list.  Puts COD1, COD2, NONC, the start
/// entries and (with `--PRE_START_WIDTH` and `--fixmotif`) the RBS block
/// into `rep`.
fn geneset_path(p: &mut P, rep: &mut Report) -> bool {
    // steps 1-4: order, sequence, gene list (maxpos = sequence length),
    // gene parser (its option messages are printed here), training sets
    let Some(order) = p.int_req("--ORDM") else { return false };
    let mut so = SeqObj::new();
    if !so.load(p) {
        return false;
    }
    let Some(genes) = geneset::load_param(p, so.seq.len() as i32) else { return false };
    let mut gp = GeneParser::new(p);
    gp.collect_coding(&so, &genes, 0, 0);
    gp.collect_noncoding(&so, &genes);
    // step 5: coding chains with the MKCHAIN margins (default 3: the start
    // and stop codons are not counted); both are built (and their option
    // errors printed) even without SEPARATE_COD
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
    // step 6: non-coding chain, trained or a random GC model
    if random_noncod(p) {
        let gc = if !p.has("--RANDOM_NONCOD_GC") {
            gp.gc_non()
        } else {
            random_noncod_gc(p)
        };
        non_chain.random_gc(gc);
    } else {
        non_chain.train_set(&gp.non, 0, 0);
    }
    // step 7: the chain tables (COD2 repeats COD1 without SEPARATE_COD)
    rep.put("--COD1", b"COD1 ", cod1_chain.table().as_bytes());
    rep.put("--COD2", b"COD2 ", if gp.separate_cod { cod2_chain.table() } else { cod1_chain.table() }.as_bytes());
    rep.put("--NONC", b"NONC ", non_chain.table().as_bytes());
    // step 8: start codon probabilities
    gp.start_stats();
    start_entries(p, rep, Some(&gp));
    // step 9: the --pre_start file (input of the Gibbs3 route, fixmotif off)
    if p.has("--pre_start") && !p.has("--fixmotif") && !super::genemodes::pre_start_out(p, &mut gp, &so, &genes) {
        return false;
    }
    // step 10: RBS block, only with --PRE_START_WIDTH and --fixmotif
    if !(p.has("--PRE_START_WIDTH") && p.has("--fixmotif")) {
        return true;
    }
    // width: total length of a pre-start sequence (upstream + start region)
    let (ok, width) = p.int_opt("--width", 0, false);
    if !ok {
        return false;
    }
    // upstream: bases before the start codon
    let (ok, upstream) = p.int_opt("--PRE_START_WIDTH", 0, true);
    if !ok {
        return false;
    }
    let (ok, intergenic) = p.bool_opt("--PRE_START_INTERGENIC", true);
    if !ok {
        return false;
    }
    let (ok, pc) = p.int_opt("--GIBBS_PSEUDOCOUNTS_MAT", 1, true);
    if !ok {
        return false;
    }
    let (ok, prec) = p.int_opt("--GIBBS_PRECISION_MAT", 3, true);
    if !ok {
        return false;
    }
    gp.collect_pre_start(&so, &genes, so.seq.len() as i32, upstream, width.wrapping_sub(upstream), intergenic);
    // a zero row total fails silently (no pre-start sequence: error, see pre_start_matrix)
    if !gp.pre_start_matrix(pc) {
        return false;
    }
    // RBSS = width - 3 (matrix rows without the start codon),
    // RBSL = upstream + 3 - width, RBSD = 1 (models.md §4.2 step 10)
    rep.put("--RBSS", b"RBSS ", width.wrapping_sub(3).to_string().as_bytes());
    rep.put("--RBSM", b"RBSM \n", gp.matrix_text(1, upstream, prec).as_bytes());
    rep.put("--RBSL", b"RBSL ", upstream.wrapping_add(3).wrapping_sub(width).to_string().as_bytes());
    rep.put("--RBSD", b"RBSD ", b"1");
    if !p.has("--log") {
        return true;
    }
    // Log: the matrix in mode 0, appended; the report adds the final
    // newline, so the log block ends with an empty line.  The request list
    // is replaced again by the model writer afterwards.
    let Some(log_name) = p.str_opt("--log", Vec::new()) else { return false };
    rep.put("--log", b"", gp.matrix_text(0, upstream, prec).as_bytes());
    rep.requests = vec![b"--log".to_vec()];
    rep.write(&log_name, true)
}

/// GC of a whole stream as a PERCENTAGE (0x43d5dc; int counts,
/// `(C + G) * 100.0 / (T + A + C + G)`, 0.0 without bases).  Used as the
/// default `--RANDOM_NONCOD` GC on the --cod/--non path, where the random
/// chain expects a fraction: this does not give a valid chain.
fn stream_gc_percent(s: &[u8]) -> f64 {
    let mut c = [0i32; 256];
    for &x in s {
        c[x as usize] += 1;
    }
    let t = c[3].wrapping_add(c[0]).wrapping_add(c[1]).wrapping_add(c[2]);
    if t == 0 {
        return 0.0;
    }
    (c[2].wrapping_add(c[1]) as f64 * 100.0) / t as f64
}

/// --mkmod without `--geneset` (0x440dba, models.md §4.3): chains from
/// two plain sequence files as in --mkmod_euk.  COD2 repeats COD1; the
/// start entries are the raw parameter strings; no RBS block (so
/// `--fixmotif` later fails with `error: no entry to report, --RBSS`).
fn cod_non_path(p: &mut P, rep: &mut Report) -> bool {
    let Some(order) = p.int_req("--ORDM") else { return false };
    let Some(cod_name) = p.str_req("--cod") else { return false };
    let Some(non_name) = p.str_req("--non") else { return false };
    let Some(cod_seq) = read_stream(&cod_name) else { return false };
    let Some(non_seq) = read_stream(&non_name) else { return false };
    let mut cod_chain = Chain::from_params(p, order, 3);
    cod_chain.l_margin = 0;
    cod_chain.r_margin = 0;
    cod_chain.train_stream(&cod_seq);
    let mut non_chain = Chain::from_params(p, order, 1);
    non_chain.l_margin = 0;
    non_chain.r_margin = 0;
    if random_noncod(p) {
        let gc = if !p.has("--RANDOM_NONCOD_GC") {
            stream_gc_percent(&non_seq)
        } else {
            random_noncod_gc(p)
        };
        non_chain.random_gc(gc);
    } else {
        non_chain.train_stream(&non_seq);
    }
    rep.put("--COD1", b"COD1 ", cod_chain.table().as_bytes());
    rep.put("--COD2", b"COD2 ", cod_chain.table().as_bytes());
    rep.put("--NONC", b"NONC ", non_chain.table().as_bytes());
    start_entries(p, rep, None);
    true
}

/// `--mkmod` (0x442432), models.md §4: validate the model keys, read the
/// output name, build the trained parts on one of the two paths, then
/// write the model (truncating).
pub fn mkmod(p: &Params) -> bool {
    let mut p = P::new(p);
    let p = &mut p;
    let mut rep = Report::default();
    if !validate(p) {
        return false;
    }
    let Some(name) = p.str_req("--mkmod") else { return false };
    let ok = if p.has("--geneset") { geneset_path(p, &mut rep) } else { cod_non_path(p, &mut rep) };
    if !ok {
        return false;
    }
    write_model(p, &mut rep, &name)
}

/// Makes a model file parseable as parameter text (0x44263c): inserts "--"
/// before the FIRST occurrence of each model keyword, anywhere in the text
/// (also inside other words), in this keyword order.
fn model_text(d: &[u8]) -> Vec<u8> {
    let mut s = d.to_vec();
    for k in [
        "PHMM", "NAME", "ORDM", "ATG_", "GTG_", "TTG_", "CTG_", "TAA_", "TAG_", "TGA_", "MINC",
        "MAXC", "MAXN", "NDEC", "CDEC", "CDCD", "CD1P", "CD2P", "COD1", "COD2", "NONC", "RBSS",
        "RBSM", "RBSL", "RBSD",
    ] {
        let k = k.as_bytes();
        if let Some(pos) = s.windows(k.len()).position(|w| w == k) {
            s.splice(pos..pos, b"--".iter().copied());
        }
    }
    s
}

/// Reads a model file into a parameter map (loader 0x40a306, keyword
/// marking 0x44263c, par-text parser 0x4480f2) and validates it (0x444850).
///
/// `Err(false)`: the loader failed (its message only).  `Err(true)`: parse
/// or validation failed, `error in fomat of model NAME` (sic) was printed.
/// Because a multi-word unquoted value is a parse error, every real
/// --mkmod model (multi-line COD1 ...) fails here (models.md §15a).
fn read_model(name: &[u8]) -> Result<P, bool> {
    let d = match io::load_file(name) {
        Some(d) => d,
        None => return Err(false),
    };
    let t = model_text(&d);
    let mut q = Params::default();
    let ok = q.parse_par_text(&t);
    let mut m = P::from_map(q.map);
    if ok && validate(&mut m) {
        Ok(m)
    } else {
        let mut e = b"error in fomat of model ".to_vec();
        e.extend_from_slice(name);
        e.push(b'\n');
        io::err(&e);
        Err(true)
    }
}

/// `--combine` (0x442c58), models.md §15a: COD1 of `--first`, COD2 of
/// `--second`, NONC of the second (of the first with `--first_NON`), start
/// entries and the RBS block of the first; the scalar keys come from the
/// first model via [`write_model`].
pub fn combine(p: &Params) -> bool {
    let mut p = P::new(p);
    let p = &mut p;
    let Some(out) = p.str_req("--combine") else { return false };
    let Some(first) = p.str_req("--first") else { return false };
    let Some(second) = p.str_req("--second") else { return false };
    // result not checked: a bad value prints and keeps false
    let first_non = p.bool_opt("--first_NON", false).1;
    let Ok(mut a) = read_model(&first) else { return false };
    let Ok(mut b) = read_model(&second) else { return false };
    let mut rep = Report::default();
    let v = a.val("--COD1");
    rep.put("--COD1", b"COD1 \n", &v);
    let v = b.val("--COD2");
    rep.put("--COD2", b"COD2 \n", &v);
    let v = if first_non { a.val("--NONC") } else { b.val("--NONC") };
    rep.put("--NONC", b"NONC \n", &v);
    for k in ["--ATG_", "--GTG_", "--TTG_", "--CTG_"] {
        let v = a.val(k);
        rep.put(k, format!("{} ", &k[2..]).as_bytes(), &v);
    }
    if a.has("--RBSS") {
        // missing RBSM/RBSL/RBSD become empty values (val inserts them)
        for (k, l) in [("--RBSS", "RBSS "), ("--RBSM", "RBSM \n"), ("--RBSL", "RBSL "), ("--RBSD", "RBSD \n")] {
            let v = a.val(k);
            rep.put(k, l.as_bytes(), &v);
        }
    }
    write_model(&mut a, &mut rep, &out)
}
