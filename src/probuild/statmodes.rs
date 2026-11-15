//! Modes working on the FastaX object (0x405dce / 0x405e8a, see
//! `fasta::FastaX`): per-record statistics and simple record filters
//! (re-notes/probuild/core.md section 7).
//!
//! - `--clean_fasta`: renumber records as `LABEL.i` and write an index of
//!   new versus old definition lines.
//! - `--stat_fasta`: per-record length, GC and nucleotide composition table.
//! - `--count_substring_n`: lengths of all N runs per record.
//! - `--substring_n_distr`: histogram of N (gap) and X (mask) run lengths.
//! - `--count_substring_atcg`: coordinates of the N-free segments.
//! - `--CLEAN_CODING`: keep coding sequences without in-frame stop codons.
//!
//! FastaX keeps every letter (unknown letters become `N`), keeps lower case
//! `acgt`, and maps `x`/`X` to `X` only for FastaX(1) (main passes
//! `--allow_x` presence for `--stat_fasta`, `--count_substring_n` and
//! `--substring_n_distr`).

use super::fasta::{runs_of, FastaX};
use super::io::{self, OFile};
use super::params::{Bytes, Params};
use std::collections::BTreeMap;

/// `ostream << fixed << setprecision(0) << v`, i.e. printf `%.0f`.  Rust's
/// `{:.0}` rounds the exact binary value half to even, like glibc, so the
/// result is identical for the finite values used here.
fn fixed0(v: f64) -> String {
    format!("{v:.0}")
}

/// --clean_fasta (0x44b7a6): write record `i` (0-based) as `>LABEL.i`, its
/// sequence and an empty line to `--clean_fasta_out`, and
/// `>LABEL.i<TAB>OLD_DEFINITION` to `--clean_fasta_index`.  The output files
/// are opened without checking (a failed open just writes nothing).
pub fn clean_fasta(p: &Params) -> bool {
    let (mut input, mut index_name, mut out_name, mut label) = (Bytes::new(), Bytes::new(), Bytes::new(), Bytes::new());
    if !(p.req_str(&mut input, "--clean_fasta")
        && p.req_str(&mut index_name, "--clean_fasta_index")
        && p.req_str(&mut out_name, "--clean_fasta_out")
        && p.req_str(&mut label, "--clean_fasta_label"))
    {
        return false;
    }
    let mut f = FastaX::new(false);
    if !f.load(&input) {
        return false;
    }
    // Opened in this order (out first, then index), both truncated.
    let mut out_file = OFile::open(&out_name, false);
    let mut index_file = OFile::open(&index_name, false);
    for (i, r) in f.recs.iter().enumerate() {
        let mut a = b">".to_vec();
        a.extend_from_slice(&label);
        a.extend_from_slice(format!(".{i}\n").as_bytes());
        a.extend_from_slice(&r.seq);
        a.extend_from_slice(b"\n\n");
        out_file.w(&a);
        let mut b = b">".to_vec();
        b.extend_from_slice(&label);
        b.extend_from_slice(format!(".{i}\t").as_bytes());
        b.extend_from_slice(&r.def);
        b.push(b'\n');
        index_file.w(&b);
    }
    out_file.close();
    index_file.close();
    true
}

/// --stat_fasta (0x44c824, report text 0x44beda): one tab-separated line per
/// record, written to the file named by `--stat_fasta` (stdout if empty):
///
/// ```text
/// DEF  LEN  GC%  A%  C%  G%  T%  N%  [X%]  [lower%]  A  C  G  T  N  [X]  [lower]
/// ```
///
/// Percentages of A/C/G/T/GC are relative to the unambiguous length
/// `LEN - N - X`; N%, X% and lower% relative to `LEN`.  All are printed with
/// `%.0f`.  Counts merge both cases (A = `a` + `A`).  The X columns appear
/// with `--allow_x`, the lower case columns with `--soft` (presence tests
/// only).  With no unambiguous letter the five composition columns are `-`.
///
/// `allow_x_obj` selects FastaX(1) and is the same presence test main does,
/// so it always equals the local `allow_x`.
pub fn stat_fasta(p: &Params, allow_x_obj: bool) -> bool {
    let mut seq = Bytes::new();
    if !p.req_str(&mut seq, "--seq") {
        return false;
    }
    let allow_x = p.has("--allow_x");
    let soft = p.has("--soft");
    let mut f = FastaX::new(allow_x_obj);
    if !f.load(&seq) {
        return false;
    }
    let mut t = Vec::new();
    for r in &f.recs {
        let mut counts = [0i64; 256];
        for &b in &r.seq {
            counts[b as usize] += 1;
        }
        let k = |x: u8| counts[x as usize];
        let len = r.seq.len() as i64;
        t.extend_from_slice(&r.def);
        t.extend_from_slice(format!("\t{}\t", len).as_bytes());
        // As in the original: an empty record ends after "DEF\t0\t" without a newline, so
        // the next record continues on the same line.
        if len == 0 {
            continue;
        }
        let acgt = len - k(b'N') - k(b'X');
        if acgt < 1 {
            t.extend_from_slice(b"-\t-\t-\t-\t-\t");
        } else {
            // Operand order of the sums as in the original (integer sums, so
            // the order does not affect the value).
            let d = acgt as f64;
            let gc = (k(b'G') + k(b'c') + k(b'C') + k(b'g')) as f64 * 100.0 / d;
            let a = (k(b'A') + k(b'a')) as f64 * 100.0 / d;
            let c = (k(b'C') + k(b'c')) as f64 * 100.0 / d;
            let g = (k(b'G') + k(b'g')) as f64 * 100.0 / d;
            let tt = (k(b'T') + k(b't')) as f64 * 100.0 / d;
            for v in [gc, a, c, g, tt] {
                t.extend_from_slice(fixed0(v).as_bytes());
                t.push(b'\t');
            }
        }
        // The original divides by `(double)(unsigned long)len` here.
        let len_f = len as u64 as f64;
        t.extend_from_slice(fixed0(k(b'N') as f64 * 100.0 / len_f).as_bytes());
        t.push(b'\t');
        if allow_x {
            t.extend_from_slice(fixed0(k(b'X') as f64 * 100.0 / len_f).as_bytes());
            t.push(b'\t');
        }
        let lower = k(b't') + k(b'a') + k(b'c') + k(b'g');
        if soft {
            t.extend_from_slice(fixed0(lower as f64 * 100.0 / len_f).as_bytes());
            t.push(b'\t');
        }
        t.extend_from_slice(
            format!(
                "{}\t{}\t{}\t{}\t{}",
                k(b'a') + k(b'A'),
                k(b'c') + k(b'C'),
                k(b'g') + k(b'G'),
                k(b't') + k(b'T'),
                k(b'N')
            )
            .as_bytes(),
        );
        if allow_x {
            t.extend_from_slice(format!("\t{}", k(b'X')).as_bytes());
        }
        if soft {
            t.extend_from_slice(format!("\t{}", lower).as_bytes());
        }
        t.push(b'\n');
    }
    // One report section with an empty head; the request list is empty, so
    // the writer adds one final "\n" after the body.
    f.report.set("--stat_fasta", b"", &t);
    let name = p.val("--stat_fasta");
    f.report.write(&name);
    true
}

/// --count_substring_n (0x44d9b4): to stdout, per record
/// `DEF<` followed by ` LEN` for each run of `N`, then a newline.
pub fn count_substring_n(p: &Params, allow_x: bool) -> bool {
    let mut seq = Bytes::new();
    if !p.req_str(&mut seq, "--seq") {
        return false;
    }
    let mut f = FastaX::new(allow_x);
    if !f.load(&seq) {
        return false;
    }
    let mut o = Vec::new();
    for r in &f.recs {
        o.extend_from_slice(&r.def);
        o.push(b'<');
        for (s, e) in runs_of(&r.seq, b'N') {
            o.extend_from_slice(format!(" {}", e - s).as_bytes());
        }
        o.push(b'\n');
    }
    io::out(&o);
    true
}

/// 0x44dc16: histogram of the lengths of runs of letter `c` over all
/// records, as `LEN<TAB>COUNT` lines in increasing length (`map<int,int>`).
/// The length is computed as `(int)end - (int)start`.  Without any run the
/// body is a single "\n", so that the report writer does not reject it as
/// empty.
fn run_length_histogram(f: &FastaX, c: u8) -> Bytes {
    let mut hist: BTreeMap<i32, i32> = BTreeMap::new();
    for r in &f.recs {
        for (s, e) in runs_of(&r.seq, c) {
            let l = (e as i32).wrapping_sub(s as i32);
            *hist.entry(l).or_insert(0) += 1;
        }
    }
    let mut o = Vec::new();
    for (l, n) in &hist {
        o.extend_from_slice(format!("{l}\t{n}\n").as_bytes());
    }
    if hist.is_empty() {
        o.push(b'\n');
    }
    o
}

/// --substring_n_distr (0x44dedc): gap length distribution (runs of `N`)
/// and, with `--allow_x`, mask length distribution (runs of `X`), written as
/// requested report sections to the `--substring_n_distr` file (stdout if
/// empty).  The writer's result is ignored, so the mode always succeeds once
/// the input was read.
pub fn substring_n_distr(p: &Params, allow_x: bool) -> bool {
    let mut seq = Bytes::new();
    if !p.req_str(&mut seq, "--seq") {
        return false;
    }
    let mut f = FastaX::new(allow_x);
    if !f.load(&seq) {
        return false;
    }
    let gaps = run_length_histogram(&f, b'N');
    f.report.set("gap_distribution", b"# gap distribution\n", &gaps);
    f.report.request("gap_distribution");
    if p.has("--allow_x") {
        let masks = run_length_histogram(&f, b'X');
        f.report.set("mask_distribution", b"# mask distribution\n", &masks);
        f.report.request("mask_distribution");
    }
    let name = p.val("--substring_n_distr");
    f.report.write(&name);
    true
}

/// --count_substring_atcg (0x44e478): to stdout, for each record the
/// segments between N runs longer than `--split_at_n`, one line per segment:
/// `START<TAB>END<TAB>LENGTH<TAB>DEF` with 1-based inclusive coordinates.
/// `--split_at_x` is read (and validated) but unused.
pub fn count_substring_atcg(p: &Params) -> bool {
    let (mut split_at_n, mut split_at_x) = (0i32, 0i32);
    if !p.opt_int(&mut split_at_n, "--split_at_n") || !p.opt_int(&mut split_at_x, "--split_at_x") {
        return false;
    }
    let mut seq = Bytes::new();
    if !p.req_str(&mut seq, "--seq") {
        return false;
    }
    let mut f = FastaX::new(false);
    if !f.load(&seq) {
        return false;
    }
    // First pass (segment table, never printed): a letter outside
    // N/X/ACGT (lower case acgt) makes find_first_not_of("ATCG", pos)
    // return pos itself, the loop never advances and the segment vector
    // grows until allocation fails ("Out of memory!", exit 1).  The original
    // only gets there after exhausting memory; the port throws the same
    // exception at once (core.md section 10).
    for r in &f.recs {
        if r.seq.iter().any(|&c| !matches!(c, b'N' | b'X' | b'A' | b'C' | b'G' | b'T')) {
            io::throw("std::bad_alloc", "std::bad_alloc");
        }
    }
    let mut o = Vec::new();
    // `(u64)split_at_n`: a negative value becomes huge, so no N run is long
    // enough and each record is one segment.
    let threshold = split_at_n as i64 as u64;
    for r in &f.recs {
        let s = &r.seq;
        let len = s.len() as u64;
        // End (0-based, exclusive) of the last splitting N run.
        let mut last = 0u64;
        let runs = runs_of(s, b'N');
        let tail = len;
        for (b, e) in runs {
            if threshold < e - b {
                if last < b {
                    o.extend_from_slice(format!("{}\t{}\t{}\t", last + 1, b, b - last).as_bytes());
                    o.extend_from_slice(&r.def);
                    o.push(b'\n');
                }
                last = e;
            }
        }
        if last < tail {
            o.extend_from_slice(format!("{}\t{}\t{}\t", last + 1, tail, tail - last).as_bytes());
            o.extend_from_slice(&r.def);
            o.push(b'\n');
        }
    }
    io::out(&o);
    true
}

/// --CLEAN_CODING (0x44d07e, record test 0x44cf66 / stop scan 0x44ce64):
/// filter a FASTA of coding sequences.
///
/// A record is "bad" (stdout `error in record DEF`) unless its length is a
/// multiple of 3, the inner region `[L, len-R)` is not empty
/// (`len - L != R`), and that region contains no enabled stop codon in frame
/// with L (L, R = `--MKCHAIN_L_MARGING`, `--MKCHAIN_R_MARGING`).  Stop codons
/// are enabled by `--TAA_`, `--TGA_`, `--TAG_`.  A good record longer than
/// `--MIN_COD` is written to the `--CLEAN_CODING` file without its margins
/// (preceded by its definition line with `--PRINT_LABEL 1`); a shorter one
/// counts as "short".  A summary line goes to stdout.
pub fn clean_coding(p: &Params) -> bool {
    let (mut seq, mut out_name) = (Bytes::new(), Bytes::new());
    let mut print_label = false;
    let (mut left_margin, mut right_margin, mut min_cod) = (0i32, 0i32, 0i32);
    let ok = p.req_str(&mut seq, "--seq")
        && p.req_str(&mut out_name, "--CLEAN_CODING")
        && p.req_bool(&mut print_label, "--PRINT_LABEL")
        && p.req_uint(&mut left_margin, "--MKCHAIN_L_MARGING")
        && p.req_uint(&mut right_margin, "--MKCHAIN_R_MARGING")
        && p.req_uint(&mut min_cod, "--MIN_COD");
    if !ok {
        return false;
    }
    let mut f = FastaX::new(false);
    if !f.load(&seq) {
        return false;
    }
    // 0x44cc30: stop codons enabled by --TAA_ --TGA_ --TAG_ (getter results
    // ignored: a bad value is only reported and leaves the codon disabled).
    // The original joins them into one '|'-separated string and tests
    // `find(codon) != npos`; as codons never contain '|', that is the same as
    // comparing with each enabled codon.
    let mut stops: Vec<&[u8]> = Vec::new();
    for (k, c) in [("--TAA_", b"TAA"), ("--TGA_", b"TGA"), ("--TAG_", b"TAG")] {
        let mut on = false;
        p.opt_bool(&mut on, k);
        if on {
            stops.push(c);
        }
    }
    // Truncated without an open check.
    let mut out_file = OFile::open(&out_name, false);
    let (mut n_records, mut good, mut bad, mut short) = (0i32, 0i32, 0i32, 0i32);
    for rec in &f.recs {
        n_records += 1;
        let len = rec.seq.len() as i64;
        let (left, right) = (left_margin as i64, right_margin as i64);
        // 0x44cf66: length a multiple of 3 and region [L, len-R) not empty;
        // 0x44ce64 walks `p = L; while p < len-R { codon = seq[p, p+3); p += 3 }`
        // (iterator `<`): a negative region is not scanned, a last partial
        // codon includes the string's NUL terminator and never matches.
        //
        // The original scans with raw pointers (`begin + L < end - R`,
        // unsigned) and copies the 3 bytes at `begin + L`. Margins above
        // INT_MAX given as unsigned wrap to negative ints, so codons can start
        // before or end after the sequence. Bytes in front of the data are the
        // std::string header (length, capacity, refcount), bytes after it the
        // NUL terminator and then other heap data. This is undefined
        // behaviour; the port never reads outside the record and resolves it
        // as follows:
        // - header and terminator bytes are never letters (exact). Other bytes
        //   outside the record are treated the same way, although the
        //   original may find a stop codon in neighbouring heap data there.
        // - L <= -2^30: `begin + L` wraps below address 0 (this static,
        //   non-PIE binary keeps its heap far below 1 GiB), the unsigned loop
        //   test fails at once and nothing is scanned.
        // - R <= -2^30: the original scans gigabytes past the record until it
        //   meets a stop codon in heap data or terminates; the port counts the
        //   record as containing a stop.
        // Otherwise only codons fully inside [0, len) can match: start at the
        // first in-frame position >= 0 and stop before len - 2.
        const FAR_OUTSIDE: i64 = -(1 << 30);
        let passes = len % 3 == 0 && len - left != right && {
            let has_stop = if left <= FAR_OUTSIDE {
                false
            } else if right <= FAR_OUTSIDE {
                true
            } else {
                let mut codon_pos = left;
                if codon_pos < 0 {
                    codon_pos += 3 * ((2 - codon_pos) / 3);
                }
                let scan_end = (len - right).min(len - 2);
                let mut found = false;
                while codon_pos < scan_end {
                    let cp = codon_pos as usize;
                    if stops.iter().any(|s| *s == &rec.seq[cp..cp + 3]) {
                        found = true;
                        break;
                    }
                    codon_pos += 3;
                }
                found
            };
            !has_stop
        };
        if !passes {
            io::out(&[b"error in record ".as_slice(), &rec.def, b"\n"].concat());
            bad += 1;
            continue;
        }
        good += 1;
        // `(u64)MIN_COD < len`: strictly longer than the minimum.
        if (min_cod as i64 as u64) < rec.seq.len() as u64 {
            if print_label {
                out_file.w(&[rec.def.as_slice(), b"\n"].concat());
            }
            // substr(L, len-R-L): a negative count wraps to a huge size_t and
            // reads to the end; L > len throws out_of_range ("error").
            let count = (len - right - left) as u64 as usize;
            let sub = if left as u64 > len as u64 {
                io::throw_substr()
            } else {
                &rec.seq[left as usize..(left as usize).saturating_add(count).min(len as usize)]
            };
            out_file.w(&[sub, b"\n"].concat());
        } else {
            short += 1;
        }
    }
    out_file.close();
    let mut m = seq.clone();
    m.extend_from_slice(format!(" records {n_records} good {good} bad {bad} short {short}\n").as_bytes());
    io::out(&m);
    true
}
