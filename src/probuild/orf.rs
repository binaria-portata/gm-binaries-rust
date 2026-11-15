//! Open reading frame scans of probuild 2.52: `--orf` (0x44a182) and
//! `--ORF_PER_SEQ` (0x44a734).  models.md §12, §13.
//!
//! An ORF here is the longest start-to-stop stretch in one frame: for
//! every stop codon (or sequence end) the scan walks upstream codon by
//! codon to the previous in-frame stop and keeps the most upstream start
//! codon that passes the `--MIN_ORF_SIZE` test (`(f - g) - 3 >= MIN`
//! with f, g the last bases of the stop and start codon).  Both strands are
//! scanned on the forward sequence; the reverse strand uses the reverse
//! complements of the start/stop triplets (e.g. TTA for TAA, CAT for ATG).
//!
//! - `--orf`: prints the start/stop counts, then the ORFs of `--seq` as
//!   coordinate lines or as FASTA-like records with the ORF sequence.
//! - `--ORF_PER_SEQ`: counts ORFs in the intergenic segments (outside the
//!   genes of `--geneset`) and prints ORFs per megabase, a measure of how
//!   much random ORF signal the non-coding sequence has.
//!
//! Positions in the ORF list are the LAST bases of the start and stop
//! codons (codon arrays are indexed by the codon's last base), 0-based.
//!
//! Setup 0x4483dc, codon array 0x448cec, direct ORFs 0x448fb2, reverse
//! ORFs 0x449440, output 0x4498f6 / 0x449de4.

use super::geneset;
use super::gp::{GeneParser, SeqObj};
use super::io::{self, OFile};
use super::params::Params;
use super::pm::{g, P};

/// Codon index meaning "no codon" (first two positions, or an N inside).
const NO_CODON: usize = 64;

/// Stop options: (key, direct codon index, reverse-strand codon index).
/// TAA 48 / TTA 60, TGA 56 / TCA 52, TAG 50 / CTA 28.
const STOPS: [(&str, usize, usize); 3] = [("--TAA_", 48, 60), ("--TGA_", 56, 52), ("--TAG_", 50, 28)];
/// Start options: ATG 14 / CAT 19, GTG 46 / CAC 17, TTG 62 / CAA 16,
/// CTG 30 / CAG 18.
const STARTS: [(&str, usize, usize); 4] =
    [("--ATG_", 14, 19), ("--GTG_", 46, 17), ("--TTG_", 62, 16), ("--CTG_", 30, 18)];

/// ORF object (0x406582).  The tables have 65 entries so that
/// [`NO_CODON`] can be looked up (always false).
struct Orf {
    /// +0x98: direct-strand start codons.
    dir_start: [bool; 65],
    /// +0xc0: reverse-strand start codons (reverse complements).
    rev_start: [bool; 65],
    /// +0xe8: direct-strand stop codons.
    dir_stop: [bool; 65],
    /// +0x110: reverse-strand stop codons.
    rev_stop: [bool; 65],
    /// +0x60: histogram of the codon array of the last scanned sequence.
    codon_count: [i32; 65],
    /// +0x48: found ORFs as (strand 1|2, left, right); all direct ORFs of a
    /// sequence, then all reverse ones, appended over all sequences.  The
    /// original's open/closed end flags are not stored (nothing reads
    /// them).
    orfs: Vec<(i32, usize, usize)>,
    /// The empty-sequence warning of `reverse` was printed (once per run).
    warned_empty: bool,
}

impl Orf {
    fn new() -> Orf {
        Orf {
            dir_start: [false; 65],
            rev_start: [false; 65],
            dir_stop: [false; 65],
            rev_stop: [false; 65],
            codon_count: [0; 65],
            orfs: Vec::new(),
            warned_empty: false,
        }
    }

    /// Reads the start/stop options (0x4483dc; bools, default false) in the
    /// order TAA_ TGA_ TAG_ ATG_ GTG_ TTG_ CTG_.  The first failure (value
    /// error printed twice) returns `false`.
    fn setup(&mut self, p: &mut P) -> bool {
        for (k, d, r) in STOPS {
            let (ok, v) = p.bool_opt(k, false);
            if !ok {
                return false;
            }
            self.dir_stop[d] = v;
            self.rev_stop[r] = v;
        }
        for (k, d, r) in STARTS {
            let (ok, v) = p.bool_opt(k, false);
            if !ok {
                return false;
            }
            self.dir_start[d] = v;
            self.rev_start[r] = v;
        }
        true
    }

    /// Codon array (0x448cec): entry i = index of the codon `s[i-2..=i]`, or
    /// [`NO_CODON`] for i < 2 or a codon with an N.  Also rebuilds the
    /// histogram, which counts positions i >= 2 only.
    fn codons(&mut self, s: &[u8]) -> Vec<usize> {
        self.codon_count = [0; 65];
        let mut c = vec![NO_CODON; s.len()];
        // codon_idx: last three bases in base 4; n_mask: 3 per N slot
        let (mut codon_idx, mut n_mask) = (0usize, 0usize);
        for (i, &x) in s.iter().enumerate() {
            codon_idx <<= 2;
            n_mask <<= 2;
            if x < 4 {
                codon_idx += x as usize;
            } else {
                n_mask += 3;
            }
            if i < 2 {
                continue;
            }
            codon_idx &= 63;
            n_mask &= 63;
            c[i] = if n_mask == 0 { codon_idx } else { NO_CODON };
            self.codon_count[c[i]] += 1;
        }
        c
    }

    /// Direct-strand ORFs (0x448fb2), skipped if `min_len > len`.
    ///
    /// Every position f that ends a stop codon, or lies in the last three
    /// positions (open ORF end, one per frame), closes a candidate.  The
    /// scan goes upstream g = f-3, f-6, ... until the previous stop; every
    /// start with `(f - g) - 3 >= min_len` overwrites the candidate, so the
    /// most upstream qualifying start wins.  Appends (1, g, f).
    fn direct(&mut self, c: &[usize], min_len: i32) {
        let l = c.len();
        if min_len as i64 > l as i64 {
            return;
        }
        for f in 0..l {
            if !self.dir_stop[c[f]] && (f as i64) < l as i64 - 3 {
                continue;
            }
            let mut best = None;
            let mut g = f as i64 - 3;
            while g >= 0 {
                let gu = g as usize;
                if self.dir_stop[c[gu]] {
                    break;
                }
                if self.dir_start[c[gu]] && (f as i64 - g) - 3 >= min_len as i64 {
                    best = Some(gu);
                }
                g -= 3;
            }
            if let Some(b) = best {
                self.orfs.push((1, b, f));
            }
        }
    }

    /// Reverse-strand ORFs (0x449440), the mirror image of [`Orf::direct`]:
    /// f runs from the end down, a candidate closes at every reverse stop
    /// (the stop codon is at the LEFT end on this strand) or at f <= 2, and
    /// the scan goes downstream g = f+3, f+6, ... keeping the most
    /// downstream qualifying reverse start.  Appends (2, f, g).
    ///
    /// The original iterates with pointers from `end - 1` down to `begin`;
    /// on an empty sequence (with min_len 0) `end - 1` wraps and is read, a
    /// SIGSEGV after the two `#` lines of --orf were printed. An empty
    /// sequence has no ORFs, so the port warns once, finds none and goes on.
    fn reverse(&mut self, c: &[usize], min_len: i32) {
        let l = c.len();
        if min_len as i64 > l as i64 {
            return;
        }
        if l == 0 {
            if !self.warned_empty {
                self.warned_empty = true;
                io::err(
                    b"probuild-rs: warning: ORF search on an empty sequence (no sequence letters \
                      in --seq, or an empty segment); reporting no ORFs for it \
                      (not handled by the original)\n",
                );
            }
            return;
        }
        for f in (0..l).rev() {
            if !self.rev_stop[c[f]] && f > 2 {
                continue;
            }
            let mut best = None;
            let mut g = f + 3;
            while g < l {
                if self.rev_stop[c[g]] {
                    break;
                }
                if self.rev_start[c[g]] && (g - f) as i64 - 3 >= min_len as i64 {
                    best = Some(g);
                }
                g += 3;
            }
            if let Some(b) = best {
                self.orfs.push((2, f, b));
            }
        }
    }
}

/// Reads `--seq` through the sequence object (0x43c64c).  Its stop getters
/// re-read TAA_/TGA_/TAG_ (printing their errors again); the result is
/// unused here.
fn read_seq(p: &mut P) -> Option<Vec<u8>> {
    let mut so = SeqObj::new();
    if so.load(p) {
        Some(so.seq)
    } else {
        None
    }
}

/// Rounds a margin toward zero to a multiple of 3 (C int division).
fn trunc3(x: i32) -> i32 {
    x / 3 * 3
}

/// Bases `[a, b]` (0-based inclusive) as letters, reverse complemented for
/// strand 2; N for anything not ACGT (0x43d058).
///
/// Negative margins can move `[a, b]` past either end of the sequence. The
/// original then reads neighbouring heap memory (undefined result: mostly
/// zero bytes printed as `A`, other bytes as `N` or raw characters). The port
/// prints `N` for every position outside the sequence. A range reaching 2^30
/// or more positions past the sequence makes the original run into unmapped
/// memory and terminate; the port stops with an error message instead.
fn orf_letters(s: &[u8], a: i32, b: i32, strand: i32) -> String {
    const FAR_OUTSIDE: i64 = 1 << 30;
    let len = s.len() as i64;
    if a as i64 <= -FAR_OUTSIDE || b as i64 - len >= FAR_OUTSIDE {
        io::fatal(&format!(
            "--orf: --ORF_MARGIN_LEFT/--ORF_MARGIN_RIGHT move an ORF to [{a}, {b}], \
             far outside the sequence of length {len}; use margins of a sensible size"
        ));
    }
    let code_at = |i: i32| -> u8 {
        if i >= 0 && (i as i64) < len {
            s[i as usize]
        } else {
            4
        }
    };
    let mut sub: Vec<u8> = (a..=b).map(code_at).collect();
    if strand == 2 {
        sub = sub.iter().rev().map(|&c| if c < 4 { 3 - c } else { 4 }).collect();
    }
    sub.iter().map(|&c| if c < 4 { b"ACGT"[c as usize] as char } else { 'N' }).collect()
}

/// `--orf` (0x44a182), models.md §12.
///
/// Prints `# Starts S` / `# Stops T` (occurrences on both strands, all
/// frames), finds the ORFs and reports them for `--ORF_OUT` 0 as
/// `a b length strand` lines, otherwise as `> POS_L a POS_R b STRAND(1/2)
/// s` + sequence.  `a = left + 2`, `b = right - 2`, each further moved
/// inward by the 3-rounded margins `--ORF_MARGIN_LEFT/RIGHT` (see the
/// comment in the loop); ORFs with `a >= b` are dropped.  The text plus a
/// final "\n" goes to the `--orf` file, or to stdout if its value is
/// empty.
pub fn orf(p: &Params) -> bool {
    let mut p = P::new(p);
    let p = &mut p;
    let mut o = Orf::new();
    if !o.setup(p) {
        return false;
    }
    let Some(s) = read_seq(p) else { return false };
    let c = o.codons(&s);
    // start / stop totals; NO_CODON is never in a table, so including
    // index 64 adds nothing
    let (mut n_starts, mut n_stops) = (0i32, 0i32);
    for k in 0..65 {
        if o.dir_start[k] {
            n_starts = n_starts.wrapping_add(o.codon_count[k]);
        }
        if o.rev_start[k] {
            n_starts = n_starts.wrapping_add(o.codon_count[k]);
        }
        if o.dir_stop[k] {
            n_stops = n_stops.wrapping_add(o.codon_count[k]);
        }
        if o.rev_stop[k] {
            n_stops = n_stops.wrapping_add(o.codon_count[k]);
        }
    }
    io::out(format!("# Starts {n_starts}\n# Stops {n_stops}\n").as_bytes());
    io::flush_stdout();
    let (ok, min_len) = p.int_opt("--MIN_ORF_SIZE", 0, false);
    if !ok {
        return false;
    }
    o.direct(&c, min_len);
    o.reverse(&c, min_len);
    // output getters: signed, failures only print
    let ml = trunc3(p.int_opt("--ORF_MARGIN_LEFT", 0, true).1);
    let mr = trunc3(p.int_opt("--ORF_MARGIN_RIGHT", 0, true).1);
    let out = p.int_opt("--ORF_OUT", 0, true).1;
    let mut res = String::new();
    for &(strand, left, right) in &o.orfs {
        let (left, right) = (left as i32, right as i32);
        // left/right are the last bases of the boundary codons (0-based),
        // so right - 2 is the first base of the right codon.  The +2 / -2
        // are the original's arithmetic; the results are not 1-based gene
        // ends and are printed as computed.  The left margin is applied
        // at the 5' end of the ORF's strand (left end on '+', right end
        // on '-').
        // C int arithmetic: extreme margins wrap around (and may then give
        // a >= b, i.e. no output for this ORF, as in the original).
        let (a, b) = if strand == 1 {
            (ml.wrapping_add(left).wrapping_add(2), right.wrapping_sub(2).wrapping_sub(mr))
        } else {
            (mr.wrapping_add(left).wrapping_add(2), right.wrapping_sub(2).wrapping_sub(ml))
        };
        if a < b {
            if out == 0 {
                res += &format!("{} {} {} {}\n", a, b, b.wrapping_sub(a).wrapping_add(1), strand);
            } else {
                // negative margins can put a/b outside the sequence, see
                // `orf_letters`
                res += &format!("> POS_L {a} POS_R {b} STRAND(1/2) {strand}\n");
                res += &orf_letters(&s, a, b, strand);
                res.push('\n');
            }
        }
    }
    // the report object adds endl: output ends with an extra "\n"
    res.push('\n');
    let target = p.val("--orf");
    if target.is_empty() {
        io::out(res.as_bytes());
        io::flush_stdout();
    } else {
        let mut f = OFile::open(&target, false);
        f.w(res.as_bytes());
        f.close();
    }
    true
}

/// `--ORF_PER_SEQ` (0x44a734), models.md §13.
///
/// Builds the non-coding segments of `--seq` outside the `--geneset`
/// genes (MIN_NON applies), scans each segment for ORFs (positions
/// relative to the segment) and prints, for `--ORF_OUT` 0, ORFs per
/// megabase: `n_orfs * 1000000.0 / total_length` (`%g`).  Otherwise prints
/// `a b length strand` lines with a = left + 2, b = right - 2 (no margins)
/// and `error orf` on stderr for ORFs with a >= b.
///
/// As in the original: a failure to read the sequence or the gene list (0x40f518)
/// returns SUCCESS (rc 0) with no further output.
pub fn orf_per_seq(p: &Params) -> bool {
    let mut p = P::new(p);
    let p = &mut p;
    let mut o = Orf::new();
    if !o.setup(p) {
        return false;
    }
    // only MIN_NON of the gene parser is used here
    let mut gp = GeneParser::new(p);
    // 0x40f518
    let mut so = SeqObj::new();
    if !so.load(p) {
        return true;
    }
    let Some(genes) = geneset::load_param(p, so.seq.len() as i32) else { return true };
    gp.collect_noncoding(&so, &genes);
    let (ok, min_len) = p.int_opt("--MIN_ORF_SIZE", 0, false);
    if !ok {
        return false;
    }
    let mut total_len: i32 = 0;
    for seg in &gp.non {
        let c = o.codons(seg);
        o.direct(&c, min_len);
        o.reverse(&c, min_len);
        total_len = total_len.wrapping_add(seg.len() as i32);
    }
    let out = p.int_opt("--ORF_OUT", 0, true).1;
    if out == 0 {
        let num = o.orfs.len() as f64 * 1000000.0;
        // x/0.0 in the original: inf, or for 0/0 the x86 default NaN,
        // which has the sign bit set and prints as "-nan"
        let s = if total_len == 0 {
            if num != 0.0 {
                "inf".to_string()
            } else {
                "-nan".to_string()
            }
        } else {
            g(num / total_len as f64)
        };
        io::out(s.as_bytes());
        io::out(b"\n");
        io::flush_stdout();
    } else {
        // 0x449de4
        let mut res = String::new();
        for &(strand, left, right) in &o.orfs {
            let (a, b) = (left as i32 + 2, right as i32 - 2);
            if a < b {
                res += &format!("{} {} {} {}\n", a, b, b - a + 1, strand);
            } else {
                io::err(b"error orf\n");
            }
        }
        // no final endl here
        io::out(res.as_bytes());
    }
    true
}
