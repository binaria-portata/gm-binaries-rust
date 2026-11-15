//! Traceback to state records and all output writers.
//! Port of 0x4154b0 (traceback), 0x414e00 (record typing), 0x414660
//! (numbering), 0x414880 / 0x4147d0 (gene objects), 0x414a70 / 0x414b20
//! (reverse complement, translation), 0x414c50 / 0x414d30 (`-p` / `-n`),
//! 0x4157f0 (lst), 0x417b60 (gtf), 0x417340 (gff3), 0x415db0 (tr).
//! See re-notes/output.md.
//!
//! Last stage of the pipeline:
//! 1. [`traceback`] follows the Viterbi back pointers from the best tail
//!    sentinel and turns every segment `j -> i` into a [`StateRec`]
//!    (Intergenic, Intron, Initial, Internal, Terminal or Single, with
//!    strand, 0-based inclusive coordinates and frames), then numbers genes
//!    and exons.
//! 2. [`genes`] groups the exons of each gene ([`Gene`]) and collects their
//!    bases (for `-p`/`-n`, the `-w` filter and the tr gene table).
//! 3. One writer formats the result. The formats are reproduced byte for
//!    byte, including all details (fixed `intron00001`/`cds00001` ids in
//!    GFF3, the spelling "defilne" in tr, the 10-character padded type names in
//!    lst, ...).
//!
//! Coordinates are printed 1-based (`0-based + 1`), plus the `-z` offset in
//! lst/gtf/gff3.

use crate::model::{Model, NEG};
use crate::sequence::Sequence;
use crate::sitemap::Rec;
use std::fmt::Write as _;

/// Version string printed in every header.
pub const VERSION: &str = "3.68";

/// Standard genetic code (0x5fe180), index = 16*c0 + 4*c1 + c2 with
/// A0 C1 G2 T3. Disabled stops are patched to 'x' (see `Ctx::code_table`).
const CODE_TABLE: &[u8; 64] = b"KNKNTTTTRSRSIIMIQHQHPPPPRRRRLLLLEDEDAAAAGGGGVVVV*Y*YSSSS*CWCLFLF";
/// Base code -> letter, forward strand.
const FWD: &[u8; 5] = b"ACGTN";
/// Base code -> letter of the complement.
const REVC: &[u8; 5] = b"TGCAN";

/// Segment type of a state record.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Kind {
    Initial,
    Internal,
    Terminal,
    Single,
    Intron,
    Intergenic,
}

impl Kind {
    /// 10-character name as stored in the record (+0x13)
    pub fn padded(self) -> &'static str {
        match self {
            Kind::Initial => "Initial   ",
            Kind::Internal => "Internal  ",
            Kind::Terminal => "Terminal  ",
            Kind::Single => "Single    ",
            Kind::Intron => "Intron    ",
            Kind::Intergenic => "Intergenic",
        }
    }
    /// Name without the padding (gtf `cds_type`, tr state names).
    pub fn name(self) -> &'static str {
        self.padded().trim_end()
    }
}

/// State record (0x34 bytes in the binary): one segment of the best parse.
///
/// Exon type names follow the gene's transcript direction: on the '-'
/// strand the "Initial" exon is the rightmost one. A partial exon touching
/// a sequence end is typed by its other end.
#[allow(dead_code)]
#[derive(Debug, Clone)]
pub struct StateRec {
    /// first base, 0-based inclusive (+0x00)
    pub start: i32,
    /// last base, 0-based inclusive (+0x04)
    pub end: i32,
    /// type of the left map record (+0x08)
    pub sflags: i32,
    /// type of the right map record (+0x0c)
    pub eflags: i32,
    /// "start frame" (+0x10). Exons: codon position 1..3 of the leftmost
    /// base. Introns / intergenic: raw `sub` of the left record.
    pub f0: i32,
    /// "end frame" (+0x11). Exons: codon position 1..3 of the rightmost
    /// base. Introns / intergenic: raw `sub` of the right record.
    pub f1: i32,
    /// '+', '-', or '.' (intergenic) (+0x12)
    pub strand: u8,
    /// segment type (+0x13)
    pub kind: Kind,
    /// gene number, 0 for intergenic (+0x20)
    pub gene: i32,
    /// exon number in transcript order, 0 for introns (+0x24)
    pub exon: i32,
    /// '+' if the left signal has matching -d evidence, else '-' (+0x28)
    pub ev0: u8,
    /// '+' if the right signal has matching -d evidence, else '-' (+0x29)
    pub ev1: u8,
    /// map index of the left record (+0x2c, used by -b)
    pub sidx: usize,
    /// map index of the right record (+0x30, used by -b)
    pub eidx: usize,
}

/// 0x414e00: type the segment between map records `si` and `ei`.
///
/// The tests run in the original's order, first match wins: intergenic,
/// intron+, intron-, then by the left record type (ATG+, acceptor+,
/// donor-, stop-) with the right record type choosing the exon kind.
/// Coordinates start as the two record positions and are moved inwards by
/// one where the record position is a site base that belongs to the
/// neighbouring segment (`fix_start` / `fix_end`; never for sentinels).
/// Exons include their start and stop codons; introns are exactly GT..AG.
/// Frames: '+' exons get `f0 + 1` and `f1 == 0 -> 3`; '-' exons get
/// `f0 == 0 -> 3` and `f1 + 1`. Any other type pair is fatal
/// ("error in code", exit 1).
fn make_record(map: &[Rec], si: usize, ei: usize) -> StateRec {
    let left = &map[si];
    let right = &map[ei];
    let s = left.typ;
    let e = right.typ;
    let mut r = StateRec {
        start: left.pos,
        end: right.pos,
        sflags: s,
        eflags: e,
        f0: left.sub as i32,
        f1: right.sub as i32,
        strand: b'.',
        kind: Kind::Intergenic,
        gene: 0,
        exon: 0,
        ev0: if (s as u32 & left.mark) != 0 { b'+' } else { b'-' },
        ev1: if (e as u32 & right.mark) != 0 { b'+' } else { b'-' },
        sidx: si,
        eidx: ei,
    };
    let fix_start = |r: &mut StateRec| {
        if s & 0x100 == 0 {
            r.start += 1;
        }
    };
    let fix_end = |r: &mut StateRec| {
        if e & 0x200 == 0 {
            r.end -= 1;
        }
    };
    let plus_frames = |r: &mut StateRec| {
        r.f0 += 1;
        if r.f1 == 0 {
            r.f1 = 3;
        }
    };
    let minus_frames = |r: &mut StateRec| {
        if r.f0 == 0 {
            r.f0 = 3;
        }
        r.f1 += 1;
    };
    if (s & 2 != 0 && e & 0x21 != 0) || (s & 0x10 != 0 && e & 0x21 != 0) {
        r.strand = b'.';
        r.kind = Kind::Intergenic;
        fix_start(&mut r);
        fix_end(&mut r);
    } else if s & 4 != 0 && e & 8 != 0 {
        r.strand = b'+';
        r.kind = Kind::Intron;
    } else if s & 0x80 != 0 && e & 0x40 != 0 {
        r.strand = b'-';
        r.kind = Kind::Intron;
    } else if s & 1 != 0 && (e & 4 != 0 || e & 2 != 0) {
        r.strand = b'+';
        plus_frames(&mut r);
        if e & 4 != 0 {
            r.kind = Kind::Initial;
            fix_end(&mut r);
        } else {
            r.kind = Kind::Single;
        }
    } else if s & 8 != 0 && (e & 4 != 0 || e & 2 != 0) {
        r.strand = b'+';
        plus_frames(&mut r);
        fix_start(&mut r);
        if e & 4 != 0 {
            r.kind = Kind::Internal;
            fix_end(&mut r);
        } else {
            r.kind = Kind::Terminal;
        }
    } else if s & 0x40 != 0 && (e & 0x10 != 0 || e & 0x80 != 0) {
        r.strand = b'-';
        minus_frames(&mut r);
        fix_start(&mut r);
        if e & 0x10 != 0 {
            r.kind = Kind::Initial;
        } else {
            r.kind = Kind::Internal;
            fix_end(&mut r);
        }
    } else if s & 0x20 != 0 && (e & 0x80 != 0 || e & 0x10 != 0) {
        r.strand = b'-';
        minus_frames(&mut r);
        if e & 0x80 != 0 {
            r.kind = Kind::Terminal;
            fix_end(&mut r);
        } else {
            r.kind = Kind::Single;
        }
    } else {
        crate::out("error in code\n");
        let _ = std::io::Write::flush(&mut std::io::stdout());
        std::process::exit(1);
    }
    r
}

/// 0x4154b0 + 0x414660: best parse as state records, left to right.
///
/// The end cell is the best of the 13 tail sentinels: start with N-13
/// (treated as NEG / index 0 if its score is <= NEG), then take a later one
/// only if strictly better, so the lowest index wins ties. If no sentinel
/// is reachable (index < 13) the result is empty. Otherwise follow back
/// pointers until a head sentinel (index <= 12); the segment from the head
/// sentinel is included.
pub fn traceback(map: &[Rec], back: &[i32], score: &[f64]) -> Vec<StateRec> {
    let n = score.len();
    let mut e = n - 13;
    if score[n - 13] <= NEG {
        e = 0;
    }
    let mut best = if score[n - 13] <= NEG { NEG } else { score[n - 13] };
    for k in n - 12..n {
        if best < score[k] {
            best = score[k];
            e = k;
        }
    }
    let mut recs = Vec::new();
    if e < 13 {
        return recs;
    }
    loop {
        let s = back[e] as usize;
        recs.push(make_record(map, s, e));
        e = s;
        if s <= 12 {
            break;
        }
    }
    recs.reverse();
    number(&mut recs);
    recs
}

/// 0x414660: gene and exon numbers.
///
/// Pass 1, left to right: an Intergenic record closes the current gene; the
/// next non-intergenic record opens gene g+1; introns get the gene number,
/// exons also a running exon number. Pass 2, right to left: exons on '-'
/// are renumbered from the right, so exon 1 is always the 5' exon of the
/// transcript (the counter resets at every Intergenic).
fn number(recs: &mut [StateRec]) {
    let mut g = 0;
    let mut in_gene = false;
    let mut cnt = 0;
    for r in recs.iter_mut() {
        if r.kind == Kind::Intergenic {
            in_gene = false;
            continue;
        }
        if !in_gene {
            g += 1;
            cnt = 0;
            in_gene = true;
        }
        r.gene = g;
        if r.kind != Kind::Intron {
            cnt += 1;
            r.exon = cnt;
        }
    }
    let mut c = 0;
    for r in recs.iter_mut().rev() {
        if r.kind == Kind::Intergenic {
            c = 0;
            continue;
        }
        if r.kind != Kind::Intron && r.strand == b'-' {
            c += 1;
            r.exon = c;
        }
    }
}

/// Gene object (0x414880 / 0x4147d0, 0x40 bytes): one per run of
/// non-intergenic records.
#[derive(Debug, Clone)]
pub struct Gene {
    /// start of the first (leftmost) exon, 0-based
    pub start: i32,
    /// end of the last exon, 0-based inclusive
    pub end: i32,
    /// gene number (from the first exon)
    pub number: i32,
    /// total exon length in bases (stop codon included); used by `-w`
    pub len: i32,
    /// concatenated exon base codes, left to right (forward strand)
    pub bases: Vec<u8>,
    /// number of exons
    pub exons: i32,
    /// f0 of the first exon
    pub f_first: i32,
    /// f1 of the last exon
    pub f_last: i32,
    /// strand of the first exon
    pub strand: u8,
    /// any exon with evidence on either end ('+' in ev0 or ev1); such genes
    /// are never removed by `-w`
    pub evidence: bool,
}

/// 0x414880 + 0x4147d0: build the gene objects from the state records.
pub fn genes(recs: &[StateRec], codes: &[u8]) -> Vec<Gene> {
    let mut out: Vec<Gene> = Vec::new();
    let mut open = false;
    for r in recs {
        if r.kind == Kind::Intergenic {
            open = false;
            continue;
        }
        if r.kind == Kind::Intron {
            continue;
        }
        if !open {
            out.push(Gene {
                start: r.start,
                end: r.end,
                number: r.gene,
                len: 0,
                bases: Vec::new(),
                exons: 0,
                f_first: r.f0,
                f_last: r.f1,
                strand: r.strand,
                evidence: false,
            });
            open = true;
        }
        let g = out.last_mut().unwrap();
        g.end = r.end;
        g.len += r.end - r.start + 1;
        g.bases.extend_from_slice(&codes[r.start as usize..=r.end as usize]);
        g.exons += 1;
        g.f_last = r.f1;
        if r.ev0 == b'+' || r.ev1 == b'+' {
            g.evidence = true;
        }
    }
    out
}

/// 0x414a70 + 0x414b20: translate a gene. Returns (protein, transcript
/// bases); the protein may end with '*'.
///
/// '-' genes are reverse complemented first (code 4 stays 4) and their
/// frames swapped, so the result is in transcript orientation (the original
/// mutates the gene object; here a copy is returned, which `-n` then
/// prints). Translation starts at the first full codon given by the codon
/// position of the first transcript base (1 -> 0, 2 -> 2, 3 -> 1, other
/// -> 0); partial codons at either end are dropped; a codon with an N
/// gives 'X'.
fn translate(g: &Gene, table: &[u8; 64]) -> (Vec<u8>, Vec<u8>) {
    let mut nt = g.bases.clone();
    let mut f0 = g.f_first;
    if g.strand == b'-' {
        nt.reverse();
        for c in nt.iter_mut() {
            if *c < 4 {
                *c = 3 - *c;
            }
        }
        f0 = g.f_last;
    }
    let mut k = match f0 {
        2 => 2,
        3 => 1,
        _ => 0,
    };
    let mut aa = Vec::new();
    while k + 2 < nt.len() {
        let c = &nt[k..k + 3];
        if c.iter().any(|&x| x >= 4) {
            aa.push(b'X');
        } else {
            aa.push(table[(c[0] as usize) * 16 + (c[1] as usize) * 4 + c[2] as usize]);
        }
        k += 3;
    }
    (aa, nt)
}

/// 0x414c50 / 0x414d30 line wrapping: 60 characters per line, no newline
/// after the last character, then a blank line ("\n\n").
fn wrap60(o: &mut String, s: &[u8]) {
    for (i, &c) in s.iter().enumerate() {
        o.push(c as char);
        if i + 1 != s.len() && (i + 1) % 60 == 0 {
            o.push('\n');
        }
    }
    o.push_str("\n\n");
}

/// Everything the writers need besides the records and genes.
pub struct Ctx<'a> {
    /// sequence file name as typed ("Sequence name")
    pub seqfile: &'a str,
    /// model file name as typed ("Matrices file")
    pub modfile: &'a str,
    pub seq: &'a Sequence,
    pub model: &'a Model,
    /// `-w`
    pub min_gene: i32,
    /// `-z`
    pub trace_z: bool,
    /// `-p`
    pub protein: bool,
    /// `-n`
    pub nucleotide: bool,
    /// ctime() line, including its newline
    pub date: String,
}

impl Ctx<'_> {
    /// Sequence id and coordinate offset for lst/gtf/gff3: with `-z`,
    /// defline token 2 and token 3 - 1 (a missing or zero token 3 gives an
    /// empty id and offset -1, as in the original); else the sequence
    /// id and 0.
    fn id_off(&self) -> (&str, i32) {
        if self.trace_z {
            // int arithmetic of the original (wraps for a token 3 near the
            // int limits, as do the coordinate sums with this offset)
            (self.seq.token2.as_str(), self.seq.token3.wrapping_sub(1))
        } else {
            (self.seq.seqid.as_str(), 0)
        }
    }
    /// FASTA defline; "(null)" without one (printf of a NULL `%s` in glibc).
    fn defline(&self) -> &str {
        self.seq.defline.as_deref().unwrap_or("(null)")
    }
    /// Genetic code with disabled stops (`$TAx_ON 0`) replaced by 'x'.
    fn code_table(&self) -> [u8; 64] {
        let mut t = *CODE_TABLE;
        if !self.model.taa_on {
            t[0x30] = b'x';
        }
        if !self.model.tag_on {
            t[0x32] = b'x';
        }
        if !self.model.tga_on {
            t[0x38] = b'x';
        }
        t
    }
    /// `-w` filter (lst, gtf, gff3 only): drop genes shorter than the
    /// minimum unless they have evidence. Numbering is not changed.
    fn skip(&self, g: &Gene) -> bool {
        g.len < self.min_gene && !g.evidence
    }
}

/// lst writer 0x4157f0.
///
/// Header, then one line per exon (introns print nothing, every Intergenic
/// record prints an empty line, and one more empty line follows if the
/// parse does not end intergenic). Then the optional `-p` protein and `-n`
/// nucleotide sections, which list all genes (no `-w` filter).
pub fn lst(c: &Ctx, recs: &[StateRec], gs: &[Gene]) -> String {
    let mut o = String::new();
    let (_, off) = c.id_off();
    let _ = write!(o, "Eukaryotic GeneMark.hmm version {VERSION}\n");
    let _ = write!(o, "Sequence name: {}\n", c.seqfile);
    let _ = write!(o, "FASTA defline: {}\n", c.defline());
    let _ = write!(o, "Sequence length: {} bp\n", c.seq.len());
    let _ = write!(o, "G+C content: {:.2}%\n", c.seq.gc);
    let _ = write!(o, "Matrices file: {}\n", c.modfile);
    let _ = write!(o, "{}\n\n", c.date);
    o.push_str("Predicted genes/exons\n\n");
    o.push_str("Gene Exon Strand Exon           Exon Range         Exon      Start/End\n");
    o.push_str("  #    #         Type                             Length       Frame\n\n");
    // gi: next gene object; cur: gene number of the current exon group
    let mut gi = 0usize;
    let mut cur: i32 = -1;
    let mut skip = false;
    for r in recs {
        match r.kind {
            Kind::Intergenic => o.push('\n'),
            Kind::Intron => {}
            _ => {
                if r.gene != cur {
                    cur = r.gene;
                    match gs.get(gi) {
                        Some(g) if g.number == r.gene => {
                            skip = c.skip(g);
                            gi += 1;
                        }
                        _ => {
                            // unreachable: the genes are built from these
                            // records. The original puts() this to stdout;
                            // here it only goes into the discarded buffer.
                            o.push_str("fatal error, check code here\n");
                            std::process::exit(2);
                        }
                    }
                }
                if !skip {
                    let _ = write!(
                        o,
                        "{:3} {:3}   {}  {}  {:10} {:10} {:10} {:10} {} {} {}\n",
                        r.gene,
                        r.exon,
                        r.strand as char,
                        r.kind.padded(),
                        r.start.wrapping_add(1).wrapping_add(off),
                        r.end.wrapping_add(1).wrapping_add(off),
                        r.end - r.start + 1,
                        r.f0,
                        r.f1,
                        r.ev0 as char,
                        r.ev1 as char
                    );
                }
            }
        }
    }
    if recs.last().map(|r| r.kind != Kind::Intergenic).unwrap_or(true) {
        o.push('\n');
    }
    let table = c.code_table();
    let tr: Vec<(Vec<u8>, Vec<u8>)> = if c.protein || c.nucleotide { gs.iter().map(|g| translate(g, &table)).collect() } else { Vec::new() };
    if c.protein {
        o.push_str("# protein sequence of predicted genes\n\n");
        for (g, (aa, _)) in gs.iter().zip(&tr) {
            // a final '*' is neither counted nor printed
            let n = if aa.last() == Some(&b'*') { aa.len() - 1 } else { aa.len() };
            let _ = write!(o, ">gene_{}|GeneMark.hmm|{}_aa\n", g.number, n);
            // an internal stop is printed, plus this message on stdout
            if aa[..n].contains(&b'*') {
                crate::out("Internal error! Please report to code developers!");
            }
            wrap60(&mut o, &aa[..n]);
        }
        o.push_str("# end protein sequence\n\n");
    }
    if c.nucleotide {
        o.push_str("# nucleotide sequence of predicted genes\n\n");
        for (g, (_, nt)) in gs.iter().zip(&tr) {
            // length in the header = exon bases, sequence in transcript
            // orientation
            let _ = write!(o, ">gene_{}|GeneMark.hmm|{}_nt\n", g.number, g.len);
            let s: Vec<u8> = nt.iter().map(|&x| FWD[x as usize]).collect();
            wrap60(&mut o, &s);
        }
        o.push_str("# end nucleotide sequence\n\n");
    }
    o
}

/// Header lines shared by gtf and gff3 (gff3 prepends `##gff-version 3`).
fn gff_header(c: &Ctx, o: &mut String) {
    let _ = write!(o, "# Eukaryotic GeneMark.hmm version {VERSION}\n");
    let _ = write!(o, "# Sequence name: {}\n", c.seqfile);
    let _ = write!(o, "# FASTA defline: {}\n", c.defline());
    let _ = write!(o, "# Sequence length: {} bp\n", c.seq.len());
    let _ = write!(o, "# G+C content: {:.2}%\n", c.seq.gc);
    let _ = write!(o, "# Matrices file: {}\n", c.modfile);
    let _ = write!(o, "# {}", c.date);
    let _ = write!(o, "# FASTA definition line: {}\n", c.defline());
}

/// gtf writer 0x417b60.
///
/// Per gene: `gene` and `mRNA` lines, then per exon start_codon / stop_codon
/// (in coordinate order around the CDS) and the CDS line with GTF frame
/// `(4 - f) % 3` of the transcript-start codon position, plus evidence and
/// `count "exon_total"` attributes. Introns are printed only if they lie
/// inside the sequence with at least 3 bases to spare (`start0 > 3 &&
/// end0 < L-3`), which drops partial introns at the borders.
pub fn gtf(c: &Ctx, recs: &[StateRec], gs: &[Gene]) -> String {
    let mut o = String::new();
    let (id, off) = c.id_off();
    let l = c.seq.len() as i32;
    gff_header(c, &mut o);
    let _ = write!(o, "# sequence-region {} 1 {}\n", id, l);
    let mut gi = 0usize;
    let mut cur: i32 = -1;
    let mut skip = false;
    let mut cur_gene: Option<&Gene> = None;
    for r in recs {
        let a = format!("gene_id \"{}_g\"; transcript_id \"{}_t\";", r.gene, r.gene);
        match r.kind {
            Kind::Intergenic => {}
            Kind::Intron => {
                if !skip && r.start > 3 && r.end < l - 3 {
                    if r.ev0 == b'+' && r.ev1 == b'+' {
                        let _ = write!(
                            o,
                            "{id}\tGeneMark.hmm3\tintron\t{}\t{}\t.\t{}\t{}\t{a} evidence \"1_1\";\n",
                            r.start.wrapping_add(1).wrapping_add(off),
                            r.end.wrapping_add(1).wrapping_add(off),
                            r.strand as char,
                            r.f0
                        );
                    } else {
                        let _ = write!(
                            o,
                            "{id}\tGeneMark.hmm3\tintron\t{}\t{}\t.\t{}\t{}\t{a}\n",
                            r.start.wrapping_add(1).wrapping_add(off),
                            r.end.wrapping_add(1).wrapping_add(off),
                            r.strand as char,
                            r.f0
                        );
                    }
                }
            }
            _ => {
                if r.gene != cur {
                    cur = r.gene;
                    let g = &gs[gi];
                    gi += 1;
                    skip = c.skip(g);
                    cur_gene = Some(g);
                    if !skip {
                        for t in ["gene", "mRNA"] {
                            let _ = write!(
                                o,
                                "{id}\tGeneMark.hmm3\t{t}\t{}\t{}\t.\t{}\t.\t{a}\n",
                                g.start.wrapping_add(1).wrapping_add(off),
                                g.end.wrapping_add(1).wrapping_add(off),
                                r.strand as char
                            );
                        }
                    }
                }
                if skip {
                    continue;
                }
                let st = r.strand as char;
                let k = r.kind;
                let ev = |e: u8| if e == b'+' { " evidence \"1_1\";" } else { "" };
                if r.strand == b'+' && matches!(k, Kind::Initial | Kind::Single) {
                    let _ = write!(o, "{id}\tGeneMark.hmm3\tstart_codon\t{}\t{}\t.\t{st}\t0\t{a}{}\n", r.start.wrapping_add(1).wrapping_add(off), r.start.wrapping_add(3).wrapping_add(off), ev(r.ev0));
                }
                if r.strand == b'-' && matches!(k, Kind::Terminal | Kind::Single) {
                    let _ = write!(o, "{id}\tGeneMark.hmm3\tstop_codon\t{}\t{}\t.\t{st}\t0\t{a}{}\n", r.start.wrapping_add(1).wrapping_add(off), r.start.wrapping_add(3).wrapping_add(off), ev(r.ev0));
                }
                let fr = if r.strand == b'+' { (4 - r.f0) % 3 } else { (4 - r.f1) % 3 };
                let evs = format!("{}{}", if r.ev0 == b'+' { "1_" } else { "0_" }, if r.ev1 == b'+' { "1" } else { "0" });
                let nex = cur_gene.map(|g| g.exons).unwrap_or(0);
                let tail = if evs == "0_0" {
                    format!(" cds_type \"{}\"; count \"{}_{}\";", k.name(), r.exon, nex)
                } else {
                    format!(" evidence \"{evs}\"; cds_type \"{}\"; count \"{}_{}\";", k.name(), r.exon, nex)
                };
                let _ = write!(o, "{id}\tGeneMark.hmm3\tCDS\t{}\t{}\t.\t{st}\t{fr}\t{a}{tail}\n", r.start.wrapping_add(1).wrapping_add(off), r.end.wrapping_add(1).wrapping_add(off));
                if r.strand == b'+' && matches!(k, Kind::Terminal | Kind::Single) {
                    let _ = write!(o, "{id}\tGeneMark.hmm3\tstop_codon\t{}\t{}\t.\t{st}\t0\t{a}{}\n", r.end.wrapping_sub(1).wrapping_add(off), r.end.wrapping_add(1).wrapping_add(off), ev(r.ev1));
                }
                if r.strand == b'-' && matches!(k, Kind::Initial | Kind::Single) {
                    let _ = write!(o, "{id}\tGeneMark.hmm3\tstart_codon\t{}\t{}\t.\t{st}\t0\t{a}{}\n", r.end.wrapping_sub(1).wrapping_add(off), r.end.wrapping_add(1).wrapping_add(off), ev(r.ev1));
                }
            }
        }
    }
    o
}

/// gff3 writer 0x417340.
///
/// gene, mRNA, then exon + CDS per exon (CDS phase `(4 - f) % 3`). Kept as in
/// the original: every intron is `ID=intron00001` and every CDS
/// `ID=cds00001` (literal 1), exon ids use the exon number. Same intron
/// filter as gtf; no codon lines, no evidence.
pub fn gff3(c: &Ctx, recs: &[StateRec], gs: &[Gene]) -> String {
    let mut o = String::from("##gff-version 3\n");
    let (id, off) = c.id_off();
    let l = c.seq.len() as i32;
    gff_header(c, &mut o);
    let _ = write!(o, "##sequence-region {} 1 {}\n", id, l);
    let mut gi = 0usize;
    let mut cur: i32 = -1;
    let mut skip = false;
    for r in recs {
        let st = r.strand as char;
        match r.kind {
            Kind::Intergenic => {}
            Kind::Intron => {
                if !skip && r.start > 3 && r.end < l - 3 {
                    let _ = write!(
                        o,
                        "{id}\tGeneMark.hmm3\tintron\t{}\t{}\t.\t{st}\t{}\tID=intron{:05};Parent=mRNA{:05};Name={};\n",
                        r.start.wrapping_add(1).wrapping_add(off),
                        r.end.wrapping_add(1).wrapping_add(off),
                        r.f0,
                        1,
                        r.gene,
                        r.gene
                    );
                }
            }
            _ => {
                if r.gene != cur {
                    cur = r.gene;
                    let g = &gs[gi];
                    gi += 1;
                    skip = c.skip(g);
                    if !skip {
                        let _ = write!(o, "{id}\tGeneMark.hmm3\tgene\t{}\t{}\t.\t{st}\t.\tID=gene{:05};Name={};\n", g.start.wrapping_add(1).wrapping_add(off), g.end.wrapping_add(1).wrapping_add(off), g.number, g.number);
                        let _ = write!(
                            o,
                            "{id}\tGeneMark.hmm3\tmRNA\t{}\t{}\t.\t{st}\t.\tID=mRNA{:05};Parent=gene{:05};Name={};\n",
                            g.start.wrapping_add(1).wrapping_add(off),
                            g.end.wrapping_add(1).wrapping_add(off),
                            g.number,
                            g.number,
                            g.number
                        );
                    }
                }
                if skip {
                    continue;
                }
                let _ = write!(
                    o,
                    "{id}\tGeneMark.hmm3\texon\t{}\t{}\t.\t{st}\t.\tID=exon{:05};Parent=mRNA{:05};Name={};\n",
                    r.start.wrapping_add(1).wrapping_add(off),
                    r.end.wrapping_add(1).wrapping_add(off),
                    r.exon,
                    r.gene,
                    r.gene
                );
                let f = if r.strand == b'+' { r.f0 } else { r.f1 };
                let _ = write!(
                    o,
                    "{id}\tGeneMark.hmm3\tCDS\t{}\t{}\t.\t{st}\t{}\tID=cds{:05};Parent=mRNA{:05};Name={};\n",
                    r.start.wrapping_add(1).wrapping_add(off),
                    r.end.wrapping_add(1).wrapping_add(off),
                    (4 - f) % 3,
                    1,
                    r.gene,
                    r.gene
                );
            }
        }
    }
    o
}

/// Failure of the tr writer.
pub enum TrError {
    /// exon frame outside 1..3: the original calls exit(1) (no message)
    Frame,
    /// model values on which the original aborts: the port's error message
    Port(String),
}

/// tr writer 0x415db0: training data for GeneMark-ES/ET/EP self-training.
///
/// Header, a gene table, then one line per state record with its sequence
/// (reverse complemented on '-'), preceded for exons by START / STOP_xxx
/// site windows and for full introns by DON / ACC windows. Window sizes
/// come from the `_MARGIN_P`/`_WIDTH_P` values of INI, TERM_TAA (used for
/// all stops), DONOR_0 and ACCEPTOR_0. No date, no `-w`, no `-z`.
///
/// As in the original: every line carries the *current* gene object, which starts
/// as the last gene (so a leading Intergenic line shows the last gene) and
/// advances at the first record of each new gene; exons are padded with
/// lowercase `n` to whole codons; an exon frame outside 1..3 makes the
/// original exit(1) (returned as `Err`). Windows reaching outside the
/// sequence (not handled by the original) print 'N'.
pub fn tr(c: &Ctx, recs: &[StateRec], gs: &[Gene]) -> Result<String, TrError> {
    let codes = &c.seq.codes;
    let l = codes.len() as i64;
    let at = |i: i64, table: &[u8; 5]| -> char {
        if i < 0 || i >= l {
            'N'
        } else {
            table[codes[i as usize] as usize] as char
        }
    };
    // A negative _WIDTH_P makes the original's std::string::resize throw
    // std::length_error, which is not caught (abort).
    let check_width = |n: i64, key: &str| -> Result<(), TrError> {
        if n < 0 {
            return Err(TrError::Port(format!(
                "model file {}: {key}_WIDTH_P is {n}; a site window width cannot be negative for -f tr \
                 (not handled by the original)",
                c.modfile
            )));
        }
        Ok(())
    };
    // Site windows of genes near the sequence ends can reach outside the
    // sequence; the original then prints letters looked up with heap bytes
    // (undefined memory). The port prints 'N' there and warns once. A window
    // reaching 2^30 or more bases outside makes the original run into
    // unmapped memory (SIGSEGV); the port stops with an error instead of
    // building such a window.
    let check_range = |a: i64, n: i64| -> Result<(), TrError> {
        if a < 0 || a + n > l {
            const FAR: i64 = 1 << 30;
            if a <= -FAR || a + n - l >= FAR {
                return Err(TrError::Port(format!(
                    "model file {}: a site window of -f tr (_MARGIN_P/_WIDTH_P) reaches {} bases outside the \
                     sequence; use window sizes that fit the sequence (not handled by the original)",
                    c.modfile,
                    if a < 0 { -a } else { a + n - l }
                )));
            }
            crate::warn(&format!(
                "-f tr: a site window (_MARGIN_P/_WIDTH_P of model {}) reaches outside the sequence; bases \
                 outside the sequence are printed as N (not handled by the original)",
                c.modfile
            ));
        }
        Ok(())
    };
    // n bases from a, forward / reverse complement (read right to left)
    let fwd = |a: i64, n: i64| -> Result<String, TrError> {
        check_range(a, n)?;
        Ok((0..n).map(|i| at(a + i, FWD)).collect())
    };
    let rc = |a: i64, n: i64| -> Result<String, TrError> {
        check_range(a, n)?;
        Ok((0..n).rev().map(|i| at(a + i, REVC)).collect())
    };
    let mut o = String::new();
    o.push_str("## TR formated data for training \n");
    let _ = write!(o, "## Eukaryotic GeneMark.hmm version {VERSION}\n");
    let _ = write!(o, "## mod_file: {}\n", c.modfile);
    let _ = write!(o, "## seq_from_file: {}\n", c.seqfile);
    let _ = write!(o, "## defilne: {}\n", c.defline());
    let _ = write!(o, "## seq_length: {} bp\n", l);
    let _ = write!(o, "## seq_gc: {:.0}\n", c.seq.gc);
    let _ = write!(o, "## seq_id: {}\n", c.seq.seqid);
    let _ = write!(o, "## seq_L: {}\n", 1);
    let _ = write!(o, "## seq_R: {}\n", l);
    if gs.is_empty() {
        o.push_str("# seq_id_end\n---\n");
        return Ok(o);
    }
    o.push_str("## gene_id\tcds_length\texons_in_gene\n");
    for g in gs {
        let _ = write!(o, "##\t{}\t{}\t{}\n", g.number, g.len, g.exons);
    }
    o.push_str("## state\tgene_id\tgene_length\tleft\tright\tlength\tstrand\tph_left\tph_right\tseq\n");
    let m = c.model;
    let (mg_start, w_start) = (m.start.margin_p as i64, m.start.width_p as i64);
    let (mg_stop, w_stop) = (m.stop[0].margin_p as i64, m.stop[0].width_p as i64);
    let (mg_don, w_don) = (m.donor[0].margin_p as i64, m.donor[0].width_p as i64);
    let (mg_acc, w_acc) = (m.acceptor[0].margin_p as i64, m.acceptor[0].width_p as i64);
    let mut cur = gs.last().unwrap();
    let mut gi = 0usize;
    let mut cur_no = -1;
    for r in recs {
        if r.kind != Kind::Intergenic && r.gene != cur_no {
            cur_no = r.gene;
            if let Some(g) = gs.get(gi) {
                cur = g;
                gi += 1;
            }
        }
        let (s, e) = (r.start as i64, r.end as i64);
        let n = e - s + 1;
        let mut sq = if r.strand == b'-' { rc(s, n)? } else { fwd(s, n)? };
        let line = |o: &mut String, name: &str, sq: &str| {
            let _ = write!(
                o,
                "# {}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n",
                name,
                cur.number,
                cur.len,
                s + 1,
                e + 1,
                n,
                r.strand as char,
                r.f0,
                r.f1,
                sq
            );
        };
        match r.kind {
            Kind::Intron => {
                // no site windows for partial introns at the sequence ends
                if s != 0 && e != l - 1 {
                    check_width(w_don, "DONOR_0")?;
                    let d = if r.strand == b'+' { fwd(s - mg_don, w_don)? } else { rc(mg_don + 1 - w_don + e, w_don)? };
                    line(&mut o, "DON", &d);
                    check_width(w_acc, "ACCEPTOR_0")?;
                    let a = if r.strand == b'+' { fwd(e - mg_acc - 1, w_acc)? } else { rc(mg_acc + 2 - w_acc + s, w_acc)? };
                    line(&mut o, "ACC", &a);
                }
            }
            Kind::Intergenic => {}
            k => {
                // pad to whole codons: by the codon position of the
                // transcript's first base (5') and of its last base (3')
                let (p0, p1) = if r.strand == b'+' { (r.f0, r.f1) } else { (r.f1, r.f0) };
                let pre = match p0 {
                    1 => "",
                    2 => "n",
                    3 => "nn",
                    _ => return Err(TrError::Frame),
                };
                let post = match p1 {
                    1 => "nn",
                    2 => "n",
                    3 => "",
                    _ => return Err(TrError::Frame),
                };
                sq = format!("{pre}{sq}{post}");
                if matches!(k, Kind::Initial | Kind::Single) {
                    check_width(w_start, "INI")?;
                    let w = if r.strand == b'+' { fwd(s - mg_start, w_start)? } else { rc(mg_start + 1 - w_start + e, w_start)? };
                    line(&mut o, "START", &w);
                }
                if matches!(k, Kind::Terminal | Kind::Single) {
                    check_width(w_stop, "TERM_TAA")?;
                    let w = if r.strand == b'+' { fwd(e - 2 - mg_stop, w_stop)? } else { rc(mg_stop + 3 - w_stop + s, w_stop)? };
                    // The stop codon sits at window offset MARGIN_P
                    // (`window.compare(MARGIN_P, 3, "TAA")`); any other
                    // triplet, also a shorter one at the window end, gives
                    // an empty name. An offset outside 0..=WIDTH_P makes
                    // compare throw std::out_of_range, not caught (abort).
                    if mg_stop < 0 || mg_stop > w_stop {
                        return Err(TrError::Port(format!(
                            "model file {}: TERM_TAA_MARGIN_P {mg_stop} is outside the window of TERM_TAA_WIDTH_P \
                             {w_stop} (needed for -f tr; not handled by the original)",
                            c.modfile
                        )));
                    }
                    let m0 = mg_stop as usize;
                    let cod = w.get(m0..(m0 + 3).min(w.len())).unwrap_or("");
                    let name = match cod {
                        "TAA" => "STOP_TAA",
                        "TAG" => "STOP_TAG",
                        "TGA" => "STOP_TGA",
                        _ => "",
                    };
                    line(&mut o, name, &w);
                }
            }
        }
        line(&mut o, r.kind.name(), &sq);
    }
    o.push_str("# seq_id_end\n---\n");
    Ok(o)
}
