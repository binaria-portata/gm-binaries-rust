//! Modes working on the SeqUtil object (0x405aa0, see `fasta::SeqUtil`):
//! FASTA reformatting, splitting, cutting, joining and composition
//! statistics (re-notes/probuild/core.md section 8).
//!
//! | mode | what it does |
//! |------|--------------|
//! | `--reformat_fasta` | rename records, filter by length or id, soft-mask handling, re-wrap lines |
//! | `--split` | split records at long N (and X / lower case) runs into contig files |
//! | `--cut` | write the regions of a GFF file (or their complement) to files |
//! | `--split_fasta` | distribute records over K files |
//! | `--select_fasta` | write listed records to one file each |
//! | `--stat`, `--gc` | nucleotide counts and GC content of a whole file |
//! | `--clean_join` | join records into one sequence with a gap filler |
//! | `--clean` | join, replace long N runs by fillers, write clean and masked versions |
//! | `--mask` | replace listed intervals of a single-record file by N |
//!
//! Conventions: sequence coordinates are 0-based half-open `(start, end)`
//! pairs of `u64` (`Interval`); files and trace lines show them 1-based
//! inclusive (`start + 1`, `end`).  Integer options are `i32` from `atoi`
//! and are converted to `u64` with sign extension (`as i64 as u64`) wherever
//! the original compares them with `size_t` values, so a negative option
//! becomes a huge number.  Interval arithmetic wraps like the original's
//! unsigned arithmetic.
//!
//! The piece-writing modes (`--split`, `--cut`) name their outputs
//! `BASE_1`, `BASE_2`, ... (BASE = mode value, or the `--seq` value if that
//! is empty) and can write a trace table mapping new to old names and
//! coordinates.

use super::fasta::{
    first_word_after_gt, runs_of, runs_of_set, table_allow_x, table_default, table_details,
    table_protein, table_soft, write_lines, Record, SeqUtil,
};
use super::io::{self, err, OFile};
use super::params::{atoi, Bytes, Params};
use std::collections::{BTreeMap, BTreeSet};

/// Half-open interval `[start, end)` of sequence positions (0-based).
type Interval = (u64, u64);

/// `--reformat_fasta`: the value `--mask_soft` gets when only `--low2gff` is
/// given: longer than any real run, so lower case runs are reported in the
/// GFF but never replaced.
const MASK_SOFT_GFF_ONLY: i32 = 1_000_000_000;

/// Largest size `std::string::reserve` accepts in the original's (COW)
/// libstdc++ (`max_size()`); larger requests throw `std::length_error`.
const STRING_MAX_SIZE: u64 = 0x3fff_ffff_ffff_fff9;

/// Concatenate byte slices.
fn cat(parts: &[&[u8]]) -> Bytes {
    let mut v = Vec::new();
    for p in parts {
        v.extend_from_slice(p);
    }
    v
}

/// `ostream << fixed << setprecision(prec) << v` (printf `%.precf`; a
/// negative precision is 6, as in libstdc++), formatted by the C library
/// (`pm::fixed`).  Rust's `format!` cannot be used: it panics for a
/// precision above 65535, which `--GC_PRECISION` accepts.  (The original
/// formats into a stack buffer of the precision's size and ends with
/// SIGSEGV for precisions of several million; the port prints them.)
fn fixed(v: f64, prec: i32) -> String {
    super::pm::fixed(prec, v, false)
}

/// `s.find(c, from)`.
fn find_byte(s: &[u8], c: u8, from: usize) -> Option<usize> {
    (from..s.len()).find(|&i| s[i] == c)
}
/// `s.find_first_of(set, from)`.
fn first_of(s: &[u8], set: &[u8], from: usize) -> Option<usize> {
    (from..s.len()).find(|&i| set.contains(&s[i]))
}
/// `s.find_first_not_of(set, from)`.
fn first_not_of(s: &[u8], set: &[u8], from: usize) -> Option<usize> {
    (from..s.len()).find(|&i| !set.contains(&s[i]))
}
/// `s.find(c, from) + 1` as the original computes it: not found
/// (`npos + 1`) wraps to 0, so a following search restarts at the beginning.
fn find_plus1(s: &[u8], c: u8, from: usize) -> usize {
    find_byte(s, c, from).map(|x| x + 1).unwrap_or(0)
}
/// `s.substr(p, n)`: throws `std::out_of_range` if `p > len`; `n` is clamped
/// to the end (so a "negative" `n` cast to `usize` means "to the end").
fn substr(s: &[u8], p: usize, n: usize) -> Bytes {
    if p > s.len() {
        io::throw_substr();
    }
    let e = p.saturating_add(n).min(s.len());
    s[p..e].to_vec()
}

// ---------------------------------------------------------------------------
// Interval helpers of --split / --cut (0x42ca2c, 0x42cc3a, 0x42cd4a, 0x42cfba,
// 0x42d1ce).

/// 0x42cc3a: append to `out` the runs at least `min` long, shrunk by
/// `margin` on both sides if `margin < min` (so a gap keeps `margin` letters
/// of N on each side of the cut).  `min == 0` disables this kind of gap
/// entirely.  Comparisons are unsigned (`(u64)min`), the shrink wraps.
fn filter_runs(runs: &[Interval], out: &mut Vec<Interval>, min: i32, margin: i32) {
    if min == 0 {
        return;
    }
    for &(s, e) in runs {
        if (min as i64 as u64) <= e.wrapping_sub(s) {
            let mut r = (s, e);
            if margin < min {
                r.0 = r.0.wrapping_add(margin as i64 as u64);
                r.1 = r.1.wrapping_sub(margin as i64 as u64);
            }
            out.push(r);
        }
    }
}

/// 0x42cd4a: the pieces of `[lo, hi)` not covered by `gaps` (sorted by
/// start).  Exact transcription of the original loop, including all details:
///
/// - the cursor starts at 0, not at `lo` (all callers pass `lo == 0`);
/// - a gap covering the whole range stops the loop with the cursor still 0,
///   so the whole range is output as one piece;
/// - adjacent gaps produce an empty piece `(cur, cur)` (later filters or
///   `write_pieces` may drop it).
fn complement(gaps: &[Interval], lo: u64, hi: u64) -> Vec<Interval> {
    let mut out = Vec::new();
    if gaps.is_empty() {
        out.push((lo, hi));
        return out;
    }
    let mut cur = 0u64;
    for &(s, e) in gaps {
        if e <= lo {
            continue;
        }
        if s >= hi {
            break;
        }
        if s <= lo && e >= hi {
            break;
        }
        if s <= lo {
            cur = e;
            continue;
        }
        // Overlapping the previous gap: extend the cursor.
        if s < cur {
            if e > cur {
                cur = e;
            }
            continue;
        }
        out.push((cur, s));
        cur = e;
    }
    if cur < hi {
        out.push((cur, hi));
    }
    out
}

/// 0x42cfba: cut every piece into consecutive chunks of at most `max`
/// letters (`max == 0`: unchanged).  Wrapping arithmetic as in the original.
fn split_max(v: Vec<Interval>, max: i32) -> Vec<Interval> {
    if max == 0 {
        return v;
    }
    let max_len = max as i64 as u64;
    let mut out = Vec::new();
    for (s, e) in v {
        let mut p = s;
        while p < e {
            if p.wrapping_add(max_len) < e {
                out.push((p, p.wrapping_add(max_len)));
            } else {
                out.push((p, e));
            }
            p = p.wrapping_add(max_len);
        }
    }
    out
}

/// 0x42d1ce: keep pieces at least `min` long (`min == 0`: unchanged).
fn filter_min(v: Vec<Interval>, min: i32) -> Vec<Interval> {
    if min == 0 {
        return v;
    }
    v.into_iter().filter(|&(s, e)| (min as i64 as u64) <= e.wrapping_sub(s)).collect()
}

/// 0x42d304: write each piece of `rec` at least `minlen` long to its own
/// file `BASE_n` (n counts on from `n`, across records) and append a trace
/// line for it.  Returns the next free `n`.
///
/// File content: header `>BASE_n<TAB>ID<TAB>START<TAB>END` (1-based
/// inclusive; without `ID<TAB>` if the record has no definition line), then
/// the piece through the line writer.  ID is the first word of the
/// definition (`first_word_after_gt`); a definition without one (e.g. `>`
/// only) throws `out_of_range` ("error", exit 1) after the file was created.
/// An open failure prints `error on open file: NAME` to stdout and calls
/// `exit(1)` (open output files of earlier pieces are already closed).
fn write_pieces(
    rec: &Record,
    pieces: &[Interval],
    base: &[u8],
    mut n: i32,
    lpl: i32,
    trace: &mut Bytes,
    minlen: i32,
) -> i32 {
    for &(s, e) in pieces {
        if (minlen as i64 as u64) > e.wrapping_sub(s) {
            continue;
        }
        let fname = cat(&[base, b"_", n.to_string().as_bytes()]);
        let mut f = OFile::open(&fname, false);
        if !f.is_open() {
            io::out(&cat(&[b"error on open file: ", &fname, b"\n"]));
            io::c_exit(1);
        }
        let coords = format!("{}\t{}", s.wrapping_add(1), e);
        if !rec.def.is_empty() {
            let id = match first_word_after_gt(&rec.def) {
                Some(x) => x,
                None => io::throw_substr(),
            };
            f.w(&cat(&[b">", &fname, b"\t", &id, b"\t", coords.as_bytes(), b"\n"]));
            trace.extend_from_slice(&cat(&[&fname, b"\t", &id, b"\t", coords.as_bytes(), b"\n"]));
        } else {
            f.w(&cat(&[b">", &fname, b"\t", coords.as_bytes(), b"\n"]));
            trace.extend_from_slice(&cat(&[&fname, b"\t", coords.as_bytes(), b"\n"]));
        }
        // `seq.substr(s, e - s)`: a start beyond the sequence (possible with
        // --cut regions) throws out_of_range; the header is already written
        // and is flushed when `f` is dropped during unwinding.
        if s > rec.seq.len() as u64 {
            io::throw_substr();
        }
        let sub = substr(&rec.seq, s as usize, e.wrapping_sub(s) as usize);
        let mut o = Vec::new();
        write_lines(&sub, lpl, &mut o);
        f.w(&o);
        f.close();
        n = n.wrapping_add(1);
    }
    n
}

/// Trace file of --split / --cut: header line plus the collected trace
/// lines (`--cut` adds one empty line at the end).  Open failure: stdout
/// `error on open file: NAME` and `exit(1)`.
fn write_trace(name: &[u8], body: &[u8], extra_nl: bool) {
    let mut f = OFile::open(name, false);
    if !f.is_open() {
        io::out(&cat(&[b"error on open file: ", name, b"\n"]));
        io::c_exit(1);
    }
    f.w(b"#new_label\told_label\tstart\tend\n");
    f.w(body);
    if extra_nl {
        f.w(b"\n");
    }
    f.close();
}

// ---------------------------------------------------------------------------

/// --split (0x42e992): split every record into contigs at gaps and write each
/// contig to `BASE_n`.
///
/// Gaps are runs of `N` at least `--split_at_n` long, plus runs of `X`
/// (with `--allow_x`) and of lower case `acgt` (with `--soft`) at least
/// `--split_at_x` long, each shrunk by `--nnn_margin` (`filter_runs`).  The
/// contigs between gaps (`complement`) are cut into chunks of at most
/// `--max_contig`, chunks shorter than `--min_contig` are dropped, and
/// `--x_to_n` turns every `X` into `N` before writing.  All size options
/// default to 0 = off.
pub fn split(p: &Params) -> bool {
    let (mut seq, mut base, mut trace) = (Bytes::new(), Bytes::new(), Bytes::new());
    let (mut split_at_n, mut split_at_x, mut max_contig, mut min_contig, mut lpl, mut margin) =
        (0i32, 0i32, 0i32, 0i32, 0i32, 0i32);
    let (mut soft, mut allow_x, mut x_to_n) = (false, false, false);
    let ok = p.req_raw(&mut seq, "--seq")
        && p.req_raw(&mut base, "--split")
        && p.opt_int(&mut split_at_n, "--split_at_n")
        && p.opt_int(&mut split_at_x, "--split_at_x")
        && p.opt_int(&mut max_contig, "--max_contig")
        && p.opt_int(&mut min_contig, "--min_contig")
        && p.opt_int(&mut lpl, "--letters_per_line")
        && p.opt_int(&mut margin, "--nnn_margin")
        && p.opt_bool(&mut soft, "--soft")
        && p.opt_bool(&mut allow_x, "--allow_x")
        && p.opt_bool(&mut x_to_n, "--x_to_n")
        && p.opt_raw(&mut trace, "--trace");
    if !ok {
        return false;
    }
    if base.is_empty() {
        base = seq.clone();
    }
    // Lower case is only preserved with --soft; X only with --allow_x.
    let mut su = SeqUtil::new();
    if soft {
        table_soft(&mut su.table);
    }
    if allow_x {
        table_allow_x(&mut su.table, soft);
    }
    if !su.load(&seq) {
        return false;
    }
    if su.recs.is_empty() {
        err(b"error: no sequence to split\n");
        return false;
    }
    // Output file counter, continues across records.
    let mut n = 1i32;
    let mut trace_lines = Vec::new();
    let recs = std::mem::take(&mut su.recs);
    for mut rec in recs {
        let n_runs = runs_of(&rec.seq, b'N');
        let mut x_runs = Vec::new();
        if allow_x {
            x_runs = runs_of(&rec.seq, b'X');
        }
        if soft {
            x_runs.extend(runs_of_set(&rec.seq, b"atcg"));
        }
        let mut gaps = Vec::new();
        filter_runs(&n_runs, &mut gaps, split_at_n, margin);
        filter_runs(&x_runs, &mut gaps, split_at_x, margin);
        // The original uses std::sort with a start-only comparator; the
        // starts of different runs are distinct, so a stable sort by start
        // gives the same order.
        gaps.sort_by_key(|g| g.0);
        let pieces = complement(&gaps, 0, rec.seq.len() as u64);
        let pieces = split_max(pieces, max_contig);
        let pieces = filter_min(pieces, min_contig);
        if x_to_n {
            for c in rec.seq.iter_mut() {
                if *c == b'X' {
                    *c = b'N';
                }
            }
        }
        n = write_pieces(&rec, &pieces, &base, n, lpl, &mut trace_lines, 0);
    }
    if !trace.is_empty() {
        write_trace(&trace, &trace_lines, false);
    }
    true
}

/// 0x430a7a: read the `--cut` regions file (GFF-like) into
/// `sequence id -> list of [start, end)`.
///
/// Read with the `while (!eof)` getline loop.  Lines that are blank or start
/// with `#` (after leading blanks) are skipped.  Fields are separated by runs
/// of TAB (empty fields vanish) and at most 9 are collected, so lines with
/// more than 9 fields are accepted.  Fewer than 9: stdout
/// `Data format error, 9 columns is expected in GFF<k> was found in: LINE`
/// and false (the message lacks a blank before the count, as in the
/// original).  Columns 4 and 5 (1-based inclusive) become
/// `[atoi(col4) - 1, atoi(col5))`; a start of 0 wraps to a huge value, which
/// later throws in `write_pieces`.
fn read_regions(name: &[u8], m: &mut BTreeMap<Bytes, Vec<Interval>>) -> bool {
    let lines = match io::read_lines_eof(name) {
        Some(l) => l,
        None => {
            io::out(&cat(&[b"error on open file ", name, b"\n"]));
            return false;
        }
    };
    for l in lines {
        let line_start = match first_not_of(&l, b" \t", 0) {
            Some(x) if l[x] != b'#' => x,
            _ => continue,
        };
        let mut cols: Vec<Bytes> = Vec::new();
        let mut field_start = Some(line_start);
        let mut field_end = find_byte(&l, b'\t', line_start);
        while cols.len() < 9 {
            let fs = match field_start {
                Some(x) => x,
                None => break,
            };
            let fe = field_end.unwrap_or(l.len());
            cols.push(l[fs..fe].to_vec());
            field_start = match field_end {
                Some(x) => first_not_of(&l, b"\t", x),
                None => None,
            };
            field_end = match field_start {
                Some(x) => find_byte(&l, b'\t', x),
                None => None,
            };
        }
        if cols.len() == 9 {
            let s = atoi(&cols[3]).wrapping_sub(1) as i64 as u64;
            let en = atoi(&cols[4]) as i64 as u64;
            m.entry(cols[0].clone()).or_default().push((s, en));
        } else {
            let mut msg = format!("Data format error, 9 columns is expected in GFF{} was found in: ", cols.len()).into_bytes();
            msg.extend_from_slice(&l);
            msg.push(b'\n');
            io::out(&msg);
            return false;
        }
    }
    true
}

/// --cut (0x42fbce): for each record whose id (first word of the definition)
/// appears in the `--regions` file, write its regions to `BASE_n` files.
/// With `--reverse` the record id is printed to stdout and the complement of
/// the (sorted) regions is written instead, i.e. everything outside them.
/// Pieces shorter than `--min_contig` are skipped; lines have
/// `--letters_per_line` letters (default 60).
pub fn cut(p: &Params) -> bool {
    let (mut seq, mut base, mut regions, mut trace) = (Bytes::new(), Bytes::new(), Bytes::new(), Bytes::new());
    let (mut allow_x, mut soft, mut reverse) = (false, false, false);
    let (mut lpl, mut min_contig) = (60i32, 0i32);
    let ok = p.req_str(&mut seq, "--seq")
        && p.req_raw(&mut base, "--cut")
        && p.opt_bool(&mut allow_x, "--allow_x")
        && p.opt_bool(&mut soft, "--soft")
        && p.req_raw(&mut regions, "--regions")
        && p.opt_raw(&mut trace, "--trace")
        && p.opt_int(&mut lpl, "--letters_per_line")
        && p.opt_bool(&mut reverse, "--reverse")
        && p.opt_int(&mut min_contig, "--min_contig");
    if !ok {
        return false;
    }
    if base.is_empty() {
        base = seq.clone();
    }
    let mut su = SeqUtil::new();
    if soft {
        table_soft(&mut su.table);
    }
    if allow_x {
        table_allow_x(&mut su.table, soft);
    }
    if !su.load(&seq) {
        return false;
    }
    // 0x430e8a: definition line -> its first word (the sequence id).  This
    // happens before the empty check.  A definition without a word (`>`
    // only) throws out_of_range, and so does a record without a definition
    // (text before the first '>').
    for r in su.recs.iter_mut() {
        r.def = match first_word_after_gt(&r.def) {
            Some(x) => x,
            None => io::throw_substr(),
        };
    }
    if su.recs.is_empty() {
        err(b"error: no sequence to split\n");
        return false;
    }
    let mut regions_by_id: BTreeMap<Bytes, Vec<Interval>> = BTreeMap::new();
    if !read_regions(&regions, &mut regions_by_id) {
        return false;
    }
    let mut n = 1i32;
    let mut trace_lines = Vec::new();
    for rec in &su.recs {
        if let Some(list) = regions_by_id.get_mut(&rec.def) {
            if reverse {
                io::out(&cat(&[&rec.def, b"\n"]));
                // Sorted in place (pair order: start, then end), as the
                // original sorts the vector stored in its map.
                list.sort();
                let pieces = complement(list, 0, rec.seq.len() as u64);
                n = write_pieces(rec, &pieces, &base, n, lpl, &mut trace_lines, min_contig);
            } else {
                // Regions in file order, unsorted.
                n = write_pieces(rec, list, &base, n, lpl, &mut trace_lines, min_contig);
            }
        }
    }
    if !trace.is_empty() {
        write_trace(&trace, &trace_lines, true);
    }
    true
}

/// 0x4280fe: read a two-column sequence id table (`--include_sid`,
/// `--exclude_sid`) with the `while (getline)` loop.  Blank lines are
/// skipped; key = first word, value = second word (words split on blanks
/// and tabs).  `map::insert` keeps the first entry of a repeated key.  A
/// line with a single word prints `error in file format: NAME` and stops
/// reading (entries read so far stay); an unreadable file prints
/// `error on open file for reading: NAME`.  The caller ignores both errors.
fn read_sid(name: &[u8], m: &mut BTreeMap<Bytes, Bytes>) {
    let lines = match io::read_lines(name) {
        Some(l) => l,
        None => {
            err(&cat(&[b"error on open file for reading: ", name, b"\n"]));
            return;
        }
    };
    for l in lines {
        let key_start = match first_not_of(&l, b" \t", 0) {
            Some(x) => x,
            None => continue,
        };
        let key_end = match first_of(&l, b" \t", key_start) {
            Some(x) => x,
            None => {
                err(&cat(&[b"error in file format: ", name, b"\n"]));
                return;
            }
        };
        let value_start = match first_not_of(&l, b" \t", key_end) {
            Some(x) => x,
            None => {
                err(&cat(&[b"error in file format: ", name, b"\n"]));
                return;
            }
        };
        let value_end = first_of(&l, b" \t", value_start).unwrap_or(l.len());
        m.entry(l[key_start..key_end].to_vec()).or_insert_with(|| l[value_start..value_end].to_vec());
    }
}

/// --reformat_fasta (0x42845a): the GeneMark pipelines' genome preparation
/// step.  Reads one or more FASTA files (`--in`, TAB-separated list) and
/// writes all records to `--out`:
///
/// - records shorter than `--min_contig` are dropped;
/// - each record is renamed (`ID + --label` with a running number ID over
///   all input files, or `--original` / `--ncbi_gi` / `--first_w` names),
///   then filtered by the `--include_sid` / `--exclude_sid` tables;
/// - with `--mask_soft L`, runs of lower case letters at least L long are
///   replaced by `N` (or `X` with `--allow_x`), keeping `--mask_margin`
///   lower case letters at each end; `--low2gff` writes every lower case run
///   as a GFF line;
/// - `--uppercase` (default on) converts the rest to upper case;
/// - lines have `--letters_per_line` letters (default 60).
///
/// `--trace` writes `>NAME<TAB>FILE<TAB>OLD_DEFINITION` per written record.
pub fn reformat_fasta(p: &Params) -> bool {
    let mut in_files: Vec<Bytes> = Vec::new();
    let (mut out_name, mut trace, mut label, mut low2gff, mut include_sid, mut exclude_sid) =
        (Bytes::new(), Bytes::new(), Bytes::new(), Bytes::new(), Bytes::new(), Bytes::new());
    let (mut min_contig, mut lpl, mut mask_soft, mut margin) = (0i32, 60i32, 0i32, 0i32);
    let (mut uppercase, mut allow_x, mut original, mut ncbi_gi, mut first_w, mut swap_sid) =
        (true, false, false, false, false, false);
    let ok = p.req_list(&mut in_files, b"\t", "--in")
        && p.req_str(&mut out_name, "--out")
        && p.opt_uint(&mut min_contig, "--min_contig")
        && p.opt_str(&mut trace, "--trace")
        && p.opt_str(&mut label, "--label")
        && p.opt_bool(&mut uppercase, "--uppercase")
        && p.opt_uint(&mut lpl, "--letters_per_line")
        && p.opt_bool(&mut allow_x, "--allow_x")
        && p.opt_uint(&mut mask_soft, "--mask_soft")
        && p.opt_bool(&mut original, "--original")
        && p.opt_bool(&mut ncbi_gi, "--ncbi_gi")
        && p.opt_bool(&mut first_w, "--first_w")
        && p.opt_raw(&mut low2gff, "--low2gff")
        && p.opt_uint(&mut margin, "--mask_margin")
        && p.opt_str(&mut include_sid, "--include_sid")
        && p.opt_str(&mut exclude_sid, "--exclude_sid")
        && p.opt_bool(&mut swap_sid, "--swap_sid");
    if !ok {
        // On any getter failure the original dumps the current option values
        // to stdout.  It starts with `in_files[0]` without a size check: an
        // empty list (missing or bad --in, or a value of TABs only) terminates
        // the original (SIGSEGV). The port shows an empty --in value and
        // takes the normal error exit.
        let flag = |x: bool| if x { "1" } else { "0" };
        let mut o = Vec::new();
        let first_in: &[u8] = in_files.first().map(|f| f.as_slice()).unwrap_or(b"");
        o.extend_from_slice(&cat(&[b"--in ", first_in, b"\n"]));
        o.extend_from_slice(&cat(&[b"--out ", &out_name, b"\n"]));
        o.extend_from_slice(format!("--min_contig {min_contig}\n").as_bytes());
        o.extend_from_slice(&cat(&[b"--trace ", &trace, b"\n"]));
        o.extend_from_slice(format!("--uppercase {}\n", flag(uppercase)).as_bytes());
        o.extend_from_slice(&cat(&[b"--label ", &label, b"\n"]));
        o.extend_from_slice(format!("--letters_per_line {lpl}\n").as_bytes());
        o.extend_from_slice(format!("--allow_x {}\n", flag(allow_x)).as_bytes());
        o.extend_from_slice(format!("--mask_soft {mask_soft}\n").as_bytes());
        o.extend_from_slice(format!("--original {}\n", flag(original)).as_bytes());
        o.extend_from_slice(format!("--ncbi_gi {}\n", flag(ncbi_gi)).as_bytes());
        o.extend_from_slice(format!("--first_w {}\n", flag(first_w)).as_bytes());
        o.extend_from_slice(&cat(&[b"--low2gff ", &low2gff, b"\n"]));
        o.extend_from_slice(format!("--mask_margin {margin}\n").as_bytes());
        o.extend_from_slice(&cat(&[b"--include_sid ", &include_sid, b"\n"]));
        o.extend_from_slice(&cat(&[b"--exclude_sid ", &exclude_sid, b"\n"]));
        o.extend_from_slice(format!("--swap_sid {}\n", flag(swap_sid)).as_bytes());
        io::out(&o);
        return false;
    }
    // A masked run must be long enough to keep both margins (message text
    // verbatim, with its double blank).
    if mask_soft != 0 && mask_soft < margin.wrapping_mul(2) {
        io::out(b"error, check  mask margin and mask soft values\n");
        return false;
    }
    let mut su = SeqUtil::new();
    su.table = table_default();
    // --low2gff alone: report lower case runs without masking any of them.
    if !low2gff.is_empty() && mask_soft == 0 {
        mask_soft = MASK_SOFT_GFF_ONLY;
        margin = 0;
    }
    // Lower case must survive reading whenever it is needed later: for
    // masking / GFF, or for output without --uppercase.
    if mask_soft != 0 || !uppercase {
        table_soft(&mut su.table);
    }
    if allow_x {
        table_allow_x(&mut su.table, !(mask_soft == 0 && uppercase));
    }
    let mut include_map = BTreeMap::new();
    let mut exclude_map = BTreeMap::new();
    if !include_sid.is_empty() {
        read_sid(&include_sid, &mut include_map);
    }
    if !exclude_sid.is_empty() {
        read_sid(&exclude_sid, &mut exclude_map);
    }
    let mut out_file = OFile::open(&out_name, false);
    if !out_file.is_open() {
        err(&cat(&[b"error on open file for writing: ", &out_name, b"\n"]));
        return false;
    }
    // Note: these two messages have no ':' after "writing" (unlike the one
    // for --out), as in the original.
    let mut trace_file = OFile::closed();
    if !trace.is_empty() {
        trace_file.reopen(&trace, false);
        if !trace_file.is_open() {
            err(&cat(&[b"error on open file for writing ", &trace, b"\n"]));
            return false;
        }
    }
    let mut gff_file = OFile::closed();
    if !low2gff.is_empty() {
        gff_file.reopen(&low2gff, false);
        if !gff_file.is_open() {
            err(&cat(&[b"error on open file for writing ", &low2gff, b"\n"]));
            return false;
        }
    }
    // Running record number over all input files, counted before the length
    // filter (skipped records use up a number).
    let mut id = 0i32;
    for file in &in_files {
        su.recs.clear();
        // A read failure ends the mode; outputs written so far are kept
        // (flushed by the drops).
        if !su.load(file) {
            return false;
        }
        let recs = std::mem::take(&mut su.recs);
        for mut rec in recs {
            id = id.wrapping_add(1);
            if (rec.seq.len() as u64) < (min_contig as i64 as u64) {
                continue;
            }
            let def = &rec.def;
            let mut name = if original {
                // Definition without the '>' (throws for an empty definition).
                substr(def, 1, usize::MAX)
            } else if ncbi_gi {
                // ">gi|123|ref|NC_000001.1|..." -> "NC_000001": text after the
                // third '|' up to the next '.'.  A missing '|' restarts the
                // search at 0 (npos + 1 == 0).
                let mut pos = find_plus1(def, b'|', 0);
                pos = find_plus1(def, b'|', pos);
                pos = find_plus1(def, b'|', pos);
                let dot = find_byte(def, b'.', pos);
                substr(def, pos, dot.map(|e| e.wrapping_sub(pos)).unwrap_or(usize::MAX))
            } else if first_w {
                // First word after the first '>' (or from the start if there
                // is no '>'); no word throws out_of_range.
                let pos = find_plus1(def, b'>', 0);
                let pos = match first_not_of(def, b" \t", pos) {
                    Some(x) => x,
                    None => io::throw_substr(),
                };
                let word_end = first_of(def, b" \t", pos);
                substr(def, pos, word_end.map(|e| e - pos).unwrap_or(usize::MAX))
            } else {
                cat(&[id.to_string().as_bytes(), &label])
            };
            // The filters test the option value, not the table: an unreadable
            // include file therefore drops every record.
            if !include_sid.is_empty() {
                match include_map.get(&name) {
                    None => continue,
                    Some(v) => {
                        if swap_sid {
                            name = v.clone();
                        }
                    }
                }
            }
            if !exclude_sid.is_empty() && exclude_map.contains_key(&name) {
                continue;
            }
            // Both lines end with std::endl in the original, i.e. they are
            // flushed; visible when the masking below terminates.
            out_file.w(&cat(&[b">", &name, b"\n"]));
            out_file.flush();
            if !trace.is_empty() {
                trace_file.w(&cat(&[b">", &name, b"\t", file, b"\t", &rec.def, b"\n"]));
                trace_file.flush();
            }
            if mask_soft != 0 {
                let min_run = mask_soft as i64 as u64;
                let margin_u = margin as i64 as u64;
                let fill = if allow_x { b'X' } else { b'N' };
                for (s, e) in runs_of_set(&rec.seq, b"acgtnx") {
                    if min_run <= e - s {
                        let fill_start = s.wrapping_add(margin_u) as i64;
                        let fill_end = e.wrapping_sub(margin_u) as i64;
                        // memset(begin + start, fill, end - start). With the
                        // margin check above, the cases below need wrapped
                        // (out of range) option values above INT_MAX that wrap:
                        // - negative length or a fill of a gigabyte or more
                        //   (margin near +-2^31): the original terminates
                        //   (SIGSEGV); the port stops with an error message;
                        // - a smaller negative margin widens the run past the
                        //   sequence. The original then overwrites the
                        //   std::string header in front of the data (harmless
                        //   up to 16 bytes, after that it overwrites the length
                        //   and aborts later) and the heap behind it (harmless
                        //   within the allocation's slack, which depends on
                        //   how the string grew, otherwise glibc aborts).
                        //   Undefined behaviour; the port masks only inside
                        //   the sequence, which is what the original produces
                        //   whenever it survives.
                        if fill_end < fill_start || fill_end - fill_start >= 1 << 30 {
                            io::fatal(&format!(
                                "--reformat_fasta: --mask_margin {} is out of range \
                                 (it must be a small non-negative number)",
                                margin as u32
                            ));
                        }
                        let len = rec.seq.len() as i64;
                        let (lo, hi) = (fill_start.clamp(0, len), fill_end.clamp(0, len));
                        for c in &mut rec.seq[lo as usize..hi as usize] {
                            *c = fill;
                        }
                    }
                    // GFF line for every lower case run, masked or not; 1-based
                    // inclusive coordinates.
                    if !low2gff.is_empty() {
                        gff_file.w(&cat(&[
                            &name,
                            format!("\tlowcase\tmask\t{}\t{}\t0\t.\t.\tsoftmask\n", s + 1, e).as_bytes(),
                        ]));
                    }
                }
            }
            if uppercase {
                rec.seq.make_ascii_uppercase();
            }
            let mut o = Vec::new();
            write_lines(&rec.seq, lpl, &mut o);
            out_file.w(&o);
        }
    }
    out_file.close();
    trace_file.close();
    gff_file.close();
    true
}

/// --split_fasta (0x42e1c6): distribute the records over `--split_numfile`
/// (K) files `BASE_1`, `BASE_2`, ..., `ceil(records / K)` records per file,
/// in input order.  Files are opened only when needed, so fewer than K files
/// can result.  Uses the protein table (works for protein FASTA); lines have
/// `--letters_per_line` letters (default 0 = one line per record).
pub fn split_fasta(p: &Params) -> bool {
    let (mut seq, mut base) = (Bytes::new(), Bytes::new());
    let (mut n_files, mut lpl) = (0i32, 0i32);
    let ok = p.req_str(&mut seq, "--seq")
        && p.opt_raw(&mut base, "--split_fasta")
        && p.req_uint(&mut n_files, "--split_numfile")
        && p.opt_int(&mut lpl, "--letters_per_line");
    if !ok {
        return false;
    }
    if base.is_empty() {
        base = seq.clone();
    }
    if n_files == 0 {
        err(b"error, parameter number_of_file should be non zero\n");
        return false;
    }
    let mut su = SeqUtil::new();
    su.table = table_protein();
    if !su.load(&seq) {
        return false;
    }
    if su.recs.is_empty() {
        err(b"error: no sequence to split\n");
        return false;
    }
    // Records per file: ceil(n / K) with unsigned 64-bit division, stored in
    // an int.
    let n_recs = su.recs.len() as u64;
    let n_files_u = n_files as i64 as u64;
    let mut per_file = (n_recs / n_files_u) as i32;
    if n_recs % n_files_u != 0 {
        per_file = per_file.wrapping_add(1);
    }
    let mut f = OFile::closed();
    let mut in_file = 0i32;
    let mut file_no = 1i32;
    for r in &su.recs {
        if !f.is_open() {
            f.reopen(&cat(&[&base, b"_", file_no.to_string().as_bytes()]), false);
            if !f.is_open() {
                err(b"error on open file for writing\n");
                return false;
            }
        }
        let mut o = cat(&[&r.def, b"\n"]);
        write_lines(&r.seq, lpl, &mut o);
        f.w(&o);
        in_file = in_file.wrapping_add(1);
        if in_file == per_file {
            f.close();
            file_no = file_no.wrapping_add(1);
            in_file = 0;
        }
    }
    f.close();
    true
}

/// --select_fasta (0x42d97a): write every record whose id is listed in
/// `--fasta_list` to a file named after the id (in the working directory),
/// as `DEF` plus the sequence lines.
///
/// List file (`while (getline)` loop): only lines containing a blank or tab
/// count, the id is the text before the first one (a line starting with a
/// blank gives the empty id).  Record id: the definition without one
/// leading '>', then its first word if it has any non-blank.
pub fn select_fasta(p: &Params) -> bool {
    let (mut seq, mut list) = (Bytes::new(), Bytes::new());
    let mut lpl = 0i32;
    let ok = p.req_str(&mut seq, "--seq") && p.req_str(&mut list, "--fasta_list") && p.opt_int(&mut lpl, "--letters_per_line");
    if !ok {
        return false;
    }
    let mut su = SeqUtil::new();
    su.table = table_protein();
    // The list is read before the FASTA file.
    let lines = match io::read_lines(&list) {
        Some(l) => l,
        None => {
            err(&cat(&[b"error on open file for reading: ", &list, b"\n"]));
            return false;
        }
    };
    let mut wanted = BTreeSet::new();
    for l in lines {
        if let Some(i) = first_of(&l, b" \t", 0) {
            wanted.insert(l[..i].to_vec());
        }
    }
    if !su.load(&seq) {
        return false;
    }
    for r in &su.recs {
        let mut id = r.def.clone();
        if !id.is_empty() && id[0] == b'>' {
            id.remove(0);
        }
        if let Some(b) = first_not_of(&id, b" \t", 0) {
            let e = first_of(&id, b" \t", b).unwrap_or(id.len());
            id = id[b..e].to_vec();
        }
        if wanted.contains(&id) {
            let mut f = OFile::open(&id, false);
            if !f.is_open() {
                err(b"error on open file for writing\n");
                return false;
            }
            let mut o = cat(&[&r.def, b"\n"]);
            write_lines(&r.seq, lpl, &mut o);
            f.w(&o);
            f.close();
        }
    }
    true
}

/// Letter counts over all records (0x423c38), indexed by byte value.
fn counts(recs: &[Record]) -> [i64; 256] {
    let mut c = [0i64; 256];
    for r in recs {
        for &b in &r.seq {
            c[b as usize] += 1;
        }
    }
    c
}

/// GC content in percent, 0x4252d4: `(g + C + G + c) * 100.0 / (unambiguous
/// letters of both cases)`, 0 if there are none.  N and X are excluded.
/// The integer sums keep the original's operand order.
fn gc_percent(c: &[i64; 256]) -> f64 {
    let k = |x: u8| c[x as usize];
    let tot = k(b't') + k(b'A') + k(b'C') + k(b'G') + k(b'T') + k(b'a') + k(b'c') + k(b'g');
    if tot == 0 {
        0.0
    } else {
        (k(b'g') + k(b'C') + k(b'G') + k(b'c')) as f64 * 100.0 / tot as f64
    }
}

/// --gc (0x42565e): GC content of the whole file (all records together) as
/// `LABEL` + value with `--GC_PRECISION` decimals (default 1), e.g.
/// `GC% = 41.3`.  Written as report section `--gc` to the file named by
/// `--gc`, or stdout if that is empty.  Uses the default (upper case) table,
/// so lower case letters count too.
pub fn gc(p: &Params) -> bool {
    let mut seq = Bytes::new();
    let mut prec = 1i32;
    let mut label = b"GC% = ".to_vec();
    if !(p.req_str(&mut seq, "--seq") && p.opt_uint(&mut prec, "--GC_PRECISION") && p.opt_raw(&mut label, "--GC_LABEL")) {
        return false;
    }
    let mut su = SeqUtil::new();
    if !su.load(&seq) {
        return false;
    }
    let c = counts(&su.recs);
    let t = cat(&[&label, fixed(gc_percent(&c), prec).as_bytes()]);
    su.report.set("--gc", b"", &t);
    su.report.write(&p.val("--gc"));
    true
}

/// --stat (0x4271d2, report text 0x425bb6): composition summary of the whole
/// file as `NAME` (padded to 15 columns) + value lines:
///
/// - default: `LABEL`, `SEQUENCE_SIZE`, `SEQUENCE_ACGT`, `NT_A` .. `NT_N`
///   (both cases merged; the default table already upper-cases), `NT_X` with
///   `--allow_x`;
/// - `--soft`: upper and lower case counted separately (`SEQUENCE_ACGT` and
///   `NT_A`.. for upper case, `SEQUENCE_atcg` and `NT_a`.. for lower case);
/// - `--details`: like `--soft`, but `SEQUENCE_SIZE` counts every letter,
///   followed by `SEQUENCE_other` and one `NT_<c>` line per other IUPAC
///   letter that occurs;
/// - always `GC` (`--GC_PRECISION` decimals) and `RECORDS`.
///
/// `LABEL` is the `--seq` value, plus the definition line if the file has
/// exactly one record with a definition.  Written as report section `--stat`
/// to the file named by `--stat` (stdout if empty).
pub fn stat(p: &Params) -> bool {
    let mut seq = Bytes::new();
    let mut prec = 1i32;
    let (mut allow_x, mut soft, mut details) = (false, false, false);
    let ok = p.req_str(&mut seq, "--seq")
        && p.opt_uint(&mut prec, "--GC_PRECISION")
        && p.opt_bool(&mut allow_x, "--allow_x")
        && p.opt_bool(&mut soft, "--soft")
        && p.opt_bool(&mut details, "--details");
    if !ok {
        io::out(b"error, required parameters are missing\n");
        return false;
    }
    let mut su = SeqUtil::new();
    if soft {
        table_soft(&mut su.table);
    }
    if allow_x {
        table_allow_x(&mut su.table, soft);
    }
    // --details replaces the table entirely (identity on all IUPAC letters).
    if details {
        su.table = table_details();
    }
    if !su.load(&seq) {
        return false;
    }
    let c = counts(&su.recs);
    let k = |x: u8| c[x as usize];
    let mut label = seq.clone();
    if su.recs.len() == 1 && !su.recs[0].def.is_empty() {
        label.push(b' ');
        label.extend_from_slice(&su.recs[0].def);
    }
    let mut t: Bytes = Vec::new();
    let line = |t: &mut Bytes, name: &str, v: i64| {
        t.extend_from_slice(format!("{name:<15}{v}\n").as_bytes());
    };
    let lab = |t: &mut Bytes| {
        t.extend_from_slice(b"LABEL          ");
        t.extend_from_slice(&label);
        t.push(b'\n');
    };
    if !soft && !details {
        lab(&mut t);
        let size = k(b'X') + k(b'N') + k(b'T') + k(b'G') + k(b'C') + k(b'A') + k(b't') + k(b'g') + k(b'c') + k(b'a');
        line(&mut t, "SEQUENCE_SIZE", size);
        line(&mut t, "SEQUENCE_ACGT", k(b'T') + k(b'G') + k(b'C') + k(b'A') + k(b't') + k(b'g') + k(b'c') + k(b'a'));
        line(&mut t, "NT_A", k(b'a') + k(b'A'));
        line(&mut t, "NT_C", k(b'c') + k(b'C'));
        line(&mut t, "NT_G", k(b'g') + k(b'G'));
        line(&mut t, "NT_T", k(b't') + k(b'T'));
        line(&mut t, "NT_N", k(b'n') + k(b'N'));
        if allow_x {
            line(&mut t, "NT_X", k(b'x') + k(b'X'));
        }
    }
    // Upper and lower case block shared by --soft and --details; only
    // SEQUENCE_SIZE differs between them.
    let upper_lower = |t: &mut Bytes, size: i64| {
        let mut line = |name: &str, v: i64| t.extend_from_slice(format!("{name:<15}{v}\n").as_bytes());
        line("SEQUENCE_SIZE", size);
        line("SEQUENCE_ACGT", k(b'T') + k(b'G') + k(b'C') + k(b'A'));
        line("NT_A", k(b'A'));
        line("NT_C", k(b'C'));
        line("NT_G", k(b'G'));
        line("NT_T", k(b'T'));
        line("NT_N", k(b'N'));
        if allow_x {
            line("NT_X", k(b'X'));
        }
        line("SEQUENCE_atcg", k(b't') + k(b'g') + k(b'c') + k(b'a'));
        line("NT_a", k(b'a'));
        line("NT_c", k(b'c'));
        line("NT_g", k(b'g'));
        line("NT_t", k(b't'));
        line("NT_n", k(b'n'));
        if allow_x {
            line("NT_x", k(b'x'));
        }
    };
    if soft && !details {
        lab(&mut t);
        // Same sum as the default block: lower case n and x are not counted.
        let size = k(b'X') + k(b'N') + k(b'T') + k(b'G') + k(b'C') + k(b'A') + k(b't') + k(b'g') + k(b'c') + k(b'a');
        upper_lower(&mut t, size);
    }
    if details {
        let all: i64 = c.iter().sum();
        lab(&mut t);
        upper_lower(&mut t, all);
        // "Other" letters: everything counted except a c g t A C G T n N.
        let skip = |i: usize| c[i] == 0 || b"acgtACGTnN".contains(&(i as u8));
        let other: i64 = (0..256).filter(|&i| !skip(i)).map(|i| c[i]).sum();
        line(&mut t, "SEQUENCE_other", other);
        // One line per other letter in byte order; "NT_" + letter is 4
        // columns, so 11 blanks pad it to 15.
        for i in 0..256 {
            if !skip(i) {
                t.extend_from_slice(b"NT_");
                t.push(su.table[i]);
                t.extend_from_slice(format!("           {}\n", c[i]).as_bytes());
            }
        }
    }
    t.extend_from_slice(format!("GC             {}\n", fixed(gc_percent(&c), prec)).as_bytes());
    t.extend_from_slice(format!("RECORDS        {}\n", su.recs.len()).as_bytes());
    su.report.set("--stat", b"", &t);
    su.report.write(&p.val("--stat"));
    true
}

/// Join records into one sequence (0x423ca4 for --clean_join, 0x427996 for
/// --clean).  Records with at least `min` unambiguous letters (A+C+G+T after
/// the default table) are appended in order, separated by `filler`.
///
/// Returns (joined sequence, number of joined records, counts A,C,G,T,N over
/// all records, the same over joined records, filler intervals in the
/// joined sequence).
///
/// The original first calls `reserve((nrec - 1) * len(filler) + total)` in
/// unsigned 64-bit arithmetic: with no records and a non-empty filler this
/// wraps to a huge size, `reserve` throws `std::length_error` and the
/// program ends with "error" (exit 1).
fn join(recs: &[Record], filler: &[u8], min: i32) -> (Bytes, i32, [i64; 5], [i64; 5], Vec<Interval>) {
    let total: u64 = recs.iter().map(|r| r.seq.len() as u64).sum();
    let cap = (recs.len() as u64).wrapping_sub(1).wrapping_mul(filler.len() as u64).wrapping_add(total);
    if cap > STRING_MAX_SIZE {
        io::throw("std::length_error", "basic_string::_S_create");
    }
    let mut joined = Vec::new();
    let mut n_joined = 0i32;
    let (mut counts_in, mut counts_out) = ([0i64; 5], [0i64; 5]);
    let mut filler_intervals = Vec::new();
    for r in recs {
        let c = counts(std::slice::from_ref(r));
        let acgtn = [c[b'A' as usize], c[b'C' as usize], c[b'G' as usize], c[b'T' as usize], c[b'N' as usize]];
        let acgt = acgtn[3] + acgtn[0] + acgtn[1] + acgtn[2];
        if (min as i64) <= acgt {
            if !joined.is_empty() {
                let s = joined.len() as u64;
                joined.extend_from_slice(filler);
                filler_intervals.push((s, joined.len() as u64));
            }
            joined.extend_from_slice(&r.seq);
            for i in 0..5 {
                counts_out[i] += acgtn[i];
            }
            n_joined += 1;
        }
        for i in 0..5 {
            counts_in[i] += acgtn[i];
        }
    }
    (joined, n_joined, counts_in, counts_out, filler_intervals)
}

/// --clean_join (0x424220): join all records with at least
/// `--MIN_CONTIG_SIZE` unambiguous letters (default 1) into one sequence,
/// separated by `--GAP_FILLER` (default empty), optionally trimming every N
/// run to at most `--MAX_N_SUBSTRING` letters (default -1 = no trimming).
/// The result goes without a newline to the file named by `--clean_join`,
/// or to stdout with a newline if that value is empty.  With `--log` a
/// summary table is appended to the log file (stdout if its value is empty).
/// `--MIN_SEQ_SIZE` is only printed in the log.
pub fn clean_join(p: &Params) -> bool {
    // Output name (empty = stdout), read before the other options.
    let out_name = p.val("--clean_join");
    let mut seq = Bytes::new();
    let mut min_contig = 1i32;
    let mut filler = Bytes::new();
    let mut max_n = -1i32;
    let mut min_seq = 1i32;
    let ok = p.req_str(&mut seq, "--seq")
        && p.opt_uint(&mut min_contig, "--MIN_CONTIG_SIZE")
        && p.opt_raw(&mut filler, "--GAP_FILLER")
        && p.opt_int(&mut max_n, "--MAX_N_SUBSTRING")
        && p.opt_int(&mut min_seq, "--MIN_SEQ_SIZE");
    if !ok {
        return false;
    }
    let mut su = SeqUtil::new();
    if !su.load(&seq) {
        return false;
    }
    let (mut joined, n_joined, counts_in, counts_out, _) = join(&su.recs, &filler, min_contig);
    if max_n >= 0 {
        // 0x424144: an N run longer than max_n is shortened to max_n by
        // erasing the excess at the start of the run.  The filler is part of
        // the joined text, so N in the filler counts too.
        let mut pos = find_byte(&joined, b'N', 0);
        while let Some(s) = pos {
            let mut e = first_not_of(&joined, b"N", s).unwrap_or(joined.len());
            let len = (e as i32).wrapping_sub(s as i32);
            if max_n < len {
                let excess = (e - s) - max_n as usize;
                joined.drain(s..s + excess);
                e = s + max_n as usize;
            }
            pos = find_byte(&joined, b'N', e);
        }
    }
    if p.has("--log") {
        let log = p.val("--log");
        let mut t = String::new();
        t += "Clean-Join      IN   OUT\n";
        t += &format!("Contigs - {} {}\n", su.recs.len(), n_joined);
        t += &format!("All letters - {}   {}\n", counts_in.iter().sum::<i64>(), joined.len());
        t += &format!(
            "Unambiguous letters - {} {}\n",
            counts_in[0] + counts_in[1] + counts_in[2] + counts_in[3],
            counts_out[0] + counts_out[1] + counts_out[2] + counts_out[3]
        );
        for (i, l) in ["A", "C", "G", "T", "N"].iter().enumerate() {
            t += &format!("{l} - {} {}\n", counts_in[i], counts_out[i]);
        }
        t += "Parameters:\n";
        t += &format!("  Minimum of unambiguous letters per record: {min_contig}\n");
        t += &format!("  Trim down substrings \"N...N\" to maximum length: {max_n}\n");
        // The filler is raw bytes, so the rest is assembled as bytes.
        let mut tb = t.into_bytes();
        tb.extend_from_slice(b"  Connect records using gap filler: ");
        tb.extend_from_slice(&filler);
        tb.extend_from_slice(format!("\n  Minimum of unambiguous letters in output: {min_seq}\n").as_bytes());
        if !log.is_empty() {
            let mut f = OFile::open(&log, true);
            f.w(&tb);
            f.close();
        } else {
            io::out(&tb);
        }
    }
    if !out_name.is_empty() {
        let mut f = OFile::open(&out_name, false);
        f.w(&joined);
        f.close();
    } else {
        io::out(&cat(&[&joined, b"\n"]));
    }
    true
}

/// --clean (0x42aa42): join records like `--clean_join` (no N trimming) and
/// write two versions, both without newline:
///
/// - `--out_clean`: the joined sequence; with `--REPLACE_N_STRING`, every N
///   run at least `len(left) + len(right)` long gets `--LEFT_GAP_FILLER`
///   written over its first letters and `--RIGHT_GAP_FILLER` over its last;
/// - `--out_masked`: the same with every gap filler interval and every
///   replaced N run turned back into N.
///
/// Fails (exit 1) with a message if fewer than `--MIN_SEQ_SIZE` unambiguous
/// letters were joined; both files are written anyway.
pub fn clean(p: &Params) -> bool {
    let (mut seq, mut out_clean, mut out_masked) = (Bytes::new(), Bytes::new(), Bytes::new());
    let mut min_contig = 1i32;
    let mut filler = Bytes::new();
    let mut min_seq = 1i32;
    let mut replace_n = false;
    let (mut left_filler, mut right_filler) = (Bytes::new(), Bytes::new());
    let ok = p.req_str(&mut seq, "--seq")
        && p.req_str(&mut out_clean, "--out_clean")
        && p.req_str(&mut out_masked, "--out_masked")
        && p.opt_uint(&mut min_contig, "--MIN_CONTIG_SIZE")
        && p.opt_raw(&mut filler, "--GAP_FILLER")
        && p.opt_int(&mut min_seq, "--MIN_SEQ_SIZE")
        && p.opt_bool(&mut replace_n, "--REPLACE_N_STRING")
        && p.opt_raw(&mut left_filler, "--LEFT_GAP_FILLER")
        && p.opt_raw(&mut right_filler, "--RIGHT_GAP_FILLER");
    if !ok {
        return false;
    }
    let mut su = SeqUtil::new();
    if !su.load(&seq) {
        return false;
    }
    let (mut joined, _n_joined, _counts_in, counts_out, mut masked_intervals) = join(&su.recs, &filler, min_contig);
    if replace_n {
        // 0x427e9e
        let needed = (left_filler.len() + right_filler.len()) as u64;
        let mut pos = find_byte(&joined, b'N', 0);
        while let Some(b) = pos {
            let e = first_not_of(&joined, b"N", b).unwrap_or(joined.len());
            if needed <= (e - b) as u64 {
                joined[b..b + left_filler.len()].copy_from_slice(&left_filler);
                let right_start = e - right_filler.len();
                joined[right_start..e].copy_from_slice(&right_filler);
                masked_intervals.push((b as u64, e as u64));
            }
            pos = find_byte(&joined, b'N', e);
        }
    }
    let mut f = OFile::open(&out_clean, false);
    f.w(&joined);
    f.close();
    for &(a, b) in &masked_intervals {
        for c in &mut joined[a as usize..b as usize] {
            *c = b'N';
        }
    }
    f.reopen(&out_masked, false);
    f.w(&joined);
    f.close();
    let acgt = counts_out[0] + counts_out[1] + counts_out[2] + counts_out[3];
    if (min_seq as i64) <= acgt {
        true
    } else {
        err(format!("error: number of unambiguous letters in output sequence {acgt} less than allowed {min_seq}\n").as_bytes());
        false
    }
}

/// --mask (0x42b7ce): replace intervals of a single-record FASTA by `N` and
/// write `DEF\nSEQ\n` (one sequence line) to `--out`.
///
/// The `--mask` file is read with the `while (!eof)` getline loop; a line
/// with exactly two blank-separated words `A B` masks the 1-based inclusive
/// interval `[A, B]`; lines with fewer or more words are ignored.  Values go
/// through `(unsigned)atoi`.  A file with not exactly one record fails
/// silently.  An interval reaching beyond the sequence fails before anything
/// is written.  `B < A - 1` is a `memset` with negative length in the
/// original (SIGSEGV); the port stops with an error naming the interval.
pub fn mask(p: &Params) -> bool {
    let (mut seq, mut mask_name, mut out_name) = (Bytes::new(), Bytes::new(), Bytes::new());
    if !(p.req_str(&mut seq, "--seq") && p.req_str(&mut mask_name, "--mask") && p.req_str(&mut out_name, "--out")) {
        return false;
    }
    let mut su = SeqUtil::new();
    if !su.load(&seq) {
        return false;
    }
    if su.recs.len() != 1 {
        return false;
    }
    let lines = match io::read_lines_eof(&mask_name) {
        Some(l) => l,
        None => {
            err(&cat(&[b" : error, opening file ", &mask_name, b"\n"]));
            return false;
        }
    };
    let mut intervals: Vec<(u32, u32)> = Vec::new();
    for l in lines {
        // Up to 3 words are collected; only lines with exactly 2 are used.
        let mut words: Vec<Bytes> = Vec::new();
        let mut word_start = first_not_of(&l, b" \t", 0);
        while words.len() < 3 {
            let ws = match word_start {
                Some(x) => x,
                None => break,
            };
            let we = first_of(&l, b" \t", ws).unwrap_or(l.len());
            words.push(l[ws..we].to_vec());
            word_start = first_not_of(&l, b" \t", we);
        }
        if words.len() == 2 {
            intervals.push((atoi(&words[0]) as u32, atoi(&words[1]) as u32));
        }
    }
    let rec = &mut su.recs[0];
    let len = rec.seq.len() as u64;
    for (a, b) in intervals {
        if a as u64 > len || b as u64 > len {
            return false;
        }
        // memset(seq + a - 1, 'N', b - (a - 1)).  For a == 0 the original
        // starts one byte before the sequence; the visible effect is a fill
        // from position 0, which is what the port does.
        let s = (a as i64) - 1;
        let e = b as i64;
        if e < s {
            // memset with a negative length: the original terminates (SIGSEGV).
            io::fatal(&format!(
                "--mask: interval \"{a} {b}\" in {} has its end before its start \
                 (expected \"start end\", 1-based, start <= end)",
                String::from_utf8_lossy(&mask_name)
            ));
        }
        for c in &mut rec.seq[s.max(0) as usize..e as usize] {
            *c = b'N';
        }
    }
    let mut f = OFile::open(&out_name, false);
    if !f.is_open() {
        // As in the original: the message names the mask file, not the output file.
        err(&cat(&[b" : error, opening file ", &mask_name, b"\n"]));
        return false;
    }
    f.w(&cat(&[&rec.def, b"\n", &rec.seq, b"\n"]));
    f.close();
    true
}
