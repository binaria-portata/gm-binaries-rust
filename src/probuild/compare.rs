//! --compare (handler 0x452264; matching 0x451ce2, value 0x451f98;
//! re-notes/probuild/models.md section 11).
//!
//! Compares two gene lists (`--source`, `--target`, read with
//! `geneset::load`) and prints an agreement score: twice the number of
//! matched source genes divided by the total number of entries in both
//! lists (1.0 = identical lists).  The GeneMark pipelines use it to compare
//! the gene sets of successive training iterations.
//!
//! - START score (`--COMPARE_START`, default on): genes with identical
//!   start and stop (both ends equal).
//! - END score (`--COMPARE_END`, default on): genes with the same stop
//!   (exact matches included).
//!
//! The stop end is `right` on the direct strand and `left` on the reverse
//! strand.  Every entry the list reader returns counts in the denominator,
//! including header and unparsable lines (strand 0), which also match each
//! other.

use super::geneset;
use super::io::{self, OFile};
use super::params::Params;
use super::pm::{fixed, P};

/// Run --compare.  Output: `START<TAB>END` (or just one of them) with
/// `fixed` precision `--COMPARE_PRECISION` (default 1) and a newline, to the
/// file named by `--compare`, or to stdout if that value is empty.
pub fn run(p: &Params) -> bool {
    let mut p = P::new(p);
    let p = &mut p;
    let Some(source_name) = p.str_req("--source") else { return false };
    let Some(target_name) = p.str_req("--target") else { return false };
    let (ok, precision) = p.int_opt("--COMPARE_PRECISION", 1, false);
    if !ok {
        return false;
    }
    let (ok, compare_start) = p.bool_opt("--COMPARE_START", true);
    if !ok {
        return false;
    }
    let (ok, compare_end) = p.bool_opt("--COMPARE_END", true);
    if !ok {
        return false;
    }
    let Some(source) = geneset::load(&source_name, -1) else { return false };
    let Some(target) = geneset::load(&target_name, -1) else { return false };
    if source.is_empty() || target.is_empty() {
        // Message text verbatim from the original (including its spelling).
        // The original ends it with std::endl, which flushes stdout.
        io::out(b"warning, gene set in comparision is empty\n");
        io::flush_stdout();
    }
    // Counters of 0x451ce2 (object +0xa8, +0xac, +0xb0, +0xb4): exact
    // matches and stop-only matches on the direct strand (strand 1) and on
    // every other strand (2 = reverse, and 0 = non-gene lines).
    let (mut exact_direct, mut exact_reverse, mut stop_direct, mut stop_reverse) = (0u32, 0u32, 0u32, 0u32);
    for x in &source {
        // First target entry on the same strand that matches ends the scan.
        // As in the original: a stop-only match also ends it, so a later exact match in
        // the target is never seen.
        for y in &target {
            if x.strand != y.strand {
                continue;
            }
            if x.strand == 1 {
                if x.left == y.left && x.right == y.right {
                    exact_direct += 1;
                    break;
                }
                if x.right == y.right {
                    stop_direct += 1;
                    break;
                }
            } else {
                if x.left == y.left && x.right == y.right {
                    exact_reverse += 1;
                    break;
                }
                if x.left == y.left {
                    stop_reverse += 1;
                    break;
                }
            }
        }
    }
    // 0x451f98.  n = (#source + #target) as `(double)(unsigned long)`.  The
    // int counters are summed in the original's order, converted to double,
    // and the score is computed as `(x + x) / n`.
    let n = (source.len() + target.len()) as f64;
    let mut out = String::new();
    if !source.is_empty() || !target.is_empty() {
        if compare_start {
            let x = (exact_reverse.wrapping_add(exact_direct) as i32) as f64;
            out += &fixed(precision, (x + x) / n, false);
            if compare_end {
                out.push('\t');
            }
        }
        if compare_end {
            let x = (exact_reverse.wrapping_add(stop_direct).wrapping_add(stop_reverse).wrapping_add(exact_direct)
                as i32) as f64;
            out += &fixed(precision, (x + x) / n, false);
        }
    } else {
        // Both lists empty: 0 in the requested precision (e.g. "0.0").
        out = fixed(precision, 0.0, false);
    }
    out.push('\n');
    // Written through the report object (key "--compare", empty head, no
    // request list): head + value + endl to stdout or to the named file
    // (truncated; an open failure is silent and still gives exit code 0).
    let target = p.val("--compare");
    if target.is_empty() {
        io::out(out.as_bytes());
        io::flush_stdout();
    } else {
        let mut f = OFile::open(&target, false);
        f.w(out.as_bytes());
        f.close();
    }
    true
}
