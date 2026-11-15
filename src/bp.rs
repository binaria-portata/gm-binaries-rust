//! Branch point model (`$BP 1`, GeneMark-ES --fungus): score data
//! (0x41f6d0 with 0x41f260 / 0x409b20 / 0x41f600, then 0x41f9a0 / 0x41b230 /
//! 0x41fa40), the extra acceptor candidates of the intron helpers
//! (0x41fad0 forward, 0x41fec0 reverse; the border-intron variants 0x41fea0 /
//! 0x41feb0 return NEG and are not modelled) and the `-b` intron statistics
//! writer 0x418900 with 0x41fd30 / 0x420070. See re-notes/bp.md.
//!
//! # The model
//!
//! In fungal genomes the branch point sits a short, fairly constant distance
//! upstream of the acceptor. A BP model lets a full intron be scored a
//! second way, as donor .. branch point .. spacer .. acceptor:
//!
//! ```text
//! donor        branch point (BRANCH, width W)      acceptor
//!   |------ DON_BP ------|#########|---- spacer s (BP_ACC) ----|
//! ```
//!
//! The Viterbi takes the better of the normal intron candidate and the BP
//! candidate, which maximises over the spacer length `s`. The BP candidate
//! reuses the normal intron content (`run`) and corrects it to the spacer
//! Markov chain over the spacer bases via prefix sums of the per-base
//! difference (spacer - intron); the intron length distribution is replaced
//! by DON_BP and BP_ACC; the acceptor site score is replaced by the ACC_BP
//! score + log TO_BP (the per-record term `arr`).
//!
//! # Pipeline position
//!
//! [`build`] runs after the k-mer index (main 0x41f6d0), [`acceptor_terms`]
//! after the site track, the candidates inside `viterbi::Dp::intron_end`,
//! and [`intron_stats`] after the main output.

use crate::model::{BranchPoint, LenDist, Model, NEG};
use crate::output::{Kind, StateRec};
use crate::sequence::kmer_index;
use crate::sitemap::Rec;
use crate::tracks::{site_fwd, site_fwd_x, site_rev, site_rev_x};
use std::fmt::Write;

/// The BP object (in main). Per-base arrays have L entries and
/// start as NEG.
pub struct BpData {
    /// +0x00: per map record, ACC_BP score + log TO_BP - SITE (0x41b230, 0x41fa40)
    pub arr: Vec<f64>,
    /// +0x08: sequence length
    pub l: i32,
    /// +0x28: forward BRANCH score of the window ending at base i (0x41c1f0 sh = -1)
    pub fwd: Vec<f64>,
    /// +0x30: reverse BRANCH score of the window starting at base i (0x41c300 sh = 0)
    pub rev: Vec<f64>,
    /// +0x38: prefix sums of spacer nf - intron nf (forward)
    pub csum_f: Vec<f64>,
    /// +0x40: prefix sums of spacer rc_nf - intron rc_nf
    pub csum_r: Vec<f64>,
    /// +0x4c: BP_ACC MAX
    pub acc_max: i32,
    /// +0x50: BP_ACC MIN - 1
    pub acc_min1: i32,
}

/// Out-of-range reads are uninitialised heap data in the binary; only reachable with
/// unusual site margins. Return NEG instead of panicking.
#[inline]
fn at(v: &[f64], i: i32) -> f64 {
    if i < 0 {
        NEG
    } else {
        v.get(i as usize).copied().unwrap_or(NEG)
    }
}

/// `lp(D, x)` as in the Viterbi: `tab_base[x]` if `lo < x < hi`, else NEG.
#[inline]
fn lp(d: &LenDist, x: i32) -> f64 {
    if d.lo() < x && x < d.hi() {
        d.base(x)
    } else {
        NEG
    }
}

/// 0x41f6d0 (base-indexed part). `kmer` is the main k-mer index (seq+0x18).
///
/// * `fwd[i]` / `rev[i]`: BRANCH site scores per base (0x41f260);
/// * spacer k-mer index: the main one if BP_SPACER_ORDER == ORDER, else a
///   new one (0x409b20, same checks and messages as the main index);
/// * `csum_f[i]` / `csum_r[i]`: running sums from 0.0, left to right, of
///   `spacer.nf[sk[i]] - intron.nf[k[i]]` (the difference is formed first,
///   then added), for `i < L - max(BP_SPACER_ORDER, ORDER)`; the rest stays
///   NEG (0x41f600).
pub fn build(m: &Model, b: &BranchPoint, codes: &[u8], kmer: &[i32]) -> Result<BpData, String> {
    let l = codes.len() as i32;
    let w = b.branch.width;
    let mut fwd = vec![NEG; l as usize];
    let mut rev = vec![NEG; l as usize];
    // 0x41f260. With MARGIN != WIDTH the windows of the first/last positions
    // reach outside the sequence buffer (uninitialised heap data in the binary); those
    // values are never used by the candidates unless |WIDTH-MARGIN| > 30.
    // Bases outside the sequence are read as A (code 0), without a warning
    // (see tracks::site_fwd_x).
    for i in (w - 1).max(0)..l {
        fwd[i as usize] = site_fwd_x(&b.branch, codes, i as i64, -1, false);
    }
    let mut i = 0;
    while i <= l - w {
        rev[i as usize] = site_rev_x(&b.branch, codes, i as i64, 0, false);
        i += 1;
    }
    // spacer k-mer index (0x409b20, only if the orders differ)
    let own;
    let spacer_kmer: &[i32] = if m.order == b.spacer_order {
        kmer
    } else {
        if l <= b.spacer_order {
            return Err(format!("gmhmme3 : error, size of sequence {l}, order {}\n", b.spacer_order));
        }
        if 2 * b.spacer_order + 2 >= 32 {
            return Err(format!("gmhmme3 : error, can't convert address for order {} to 'int'\n", b.spacer_order));
        }
        own = kmer_index(codes, b.spacer_order);
        &own
    };
    // 0x41f600
    let mut csum_f = vec![NEG; l as usize];
    let mut csum_r = vec![NEG; l as usize];
    let n = l - b.spacer_order.max(m.order);
    if n > 0 {
        let mut d = 0.0;
        for i in 0..n as usize {
            d += b.spacer.fwd_nf[spacer_kmer[i] as usize] - m.nc3.fwd_nf[kmer[i] as usize];
            csum_f[i] = d;
        }
        let mut d = 0.0;
        for i in 0..n as usize {
            d += b.spacer.rc_nf[spacer_kmer[i] as usize] - m.nc3.rc_nf[kmer[i] as usize];
            csum_r[i] = d;
        }
    }
    Ok(BpData {
        arr: Vec::new(),
        l,
        fwd,
        rev,
        csum_f,
        csum_r,
        acc_max: b.bp_acc.max,
        acc_min1: b.bp_acc.min.wrapping_sub(1),
    })
}

/// 0x41f9a0 + 0x41b230 + 0x41fa40: per map record acceptor term.
///
/// For acceptor records (`T & 0x88`) with `30 < P < L-30`: ACC_BP score of
/// the record's phase + log TO_BP; everything else NEG. Then the normal
/// site score is subtracted (`arr - SITE`, NEG if either side is below
/// -1e117). The path through the acceptor still collects `SITE[r]` the
/// normal way (forward: from the following exon segment; reverse: it is
/// already in `cell[j]` of the acceptor- record), so `arr[r] + SITE[r]`
/// replaces the normal acceptor score by the BP acceptor score; the
/// AROUND_BP and +/-100 evidence terms of SITE cancel. The SITE track
/// itself is not changed.
pub fn acceptor_terms(d: &mut BpData, b: &BranchPoint, map: &[Rec], codes: &[u8], site: &[f64]) {
    let l = codes.len() as i64;
    let mut arr = vec![NEG; map.len()];
    for (r, rec) in map.iter().enumerate() {
        let p = rec.pos as i64;
        if !(0x1e < p && p < l - 0x1e) || rec.typ & 0x88 == 0 || rec.sub > 2 {
            continue;
        }
        let s = &b.acc_bp[rec.sub as usize];
        let v = if rec.typ & 8 != 0 { site_fwd(s, codes, p, 1) } else { site_rev(s, codes, p, 2) };
        arr[r] = v + b.log_to_bp;
    }
    // -9.999999999999999e116 is the original's threshold constant
    for (a, &s) in arr.iter_mut().zip(site) {
        *a = if -9.999999999999999e+116 < *a {
            if s <= -9.999999999999999e+116 {
                NEG
            } else {
                *a - s
            }
        } else {
            NEG
        };
    }
    d.arr = arr;
}

/// Spacer range shared by the four candidate functions: `(s0, s_end)` with
/// s running over s0 .. s_end-1, or None if the intron is too short.
///
/// With `D = pi - pj`: the intron must hold DON_BP.MIN + W + BP_ACC.MIN
/// bases (`D + 1` counted as in the original); s starts at BP_ACC.MIN - 1
/// and ends before `min(D + 1 - W - DON_BP.MIN, BP_ACC.MAX)`.
#[inline]
fn spacer_range(d: &BpData, b: &BranchPoint, pj: i32, pi: i32) -> Option<(i32, i32)> {
    let d1 = pi - pj + 1;
    let w = b.branch.width;
    // 32-bit int arithmetic of the original (wraps for MIN values near the
    // int limits)
    if d1 < b.don_bp.min.wrapping_add(w).wrapping_add(b.bp_acc.min) {
        return None;
    }
    let hi = d1.wrapping_sub(w).wrapping_sub(b.don_bp.min).min(d.acc_max);
    let s0 = d.acc_min1;
    if hi <= s0 {
        return None;
    }
    Some((s0, hi))
}

/// 0x41fad0 (Viterbi) / 0x41fd30 (-b): forward intron donor `pj`, acceptor
/// `pi` (record `ri`). Returns (best, spacer index of the first maximum, -1 if none).
///
/// Per spacer `s` (branch site window ending at `k = pi - s - 1`):
/// `(((fwd[k] + lp(BP_ACC, s)) + lp(DON_BP, D-1-W-s)) + run) +
/// (csum_f[pi-2-M0] - csum_f[k])`, then `+ arr[ri]`; strict maximum. `M0`
/// is the MARGIN of ACC_BP_0 for every phase (as in the original).
pub fn cand_fwd(d: &BpData, b: &BranchPoint, run: f64, pj: i32, pi: i32, ri: usize) -> (f64, i32) {
    let Some((s0, hi)) = spacer_range(d, b, pj, pi) else {
        return (NEG, -1);
    };
    let w = b.branch.width;
    let dist = pi - pj;
    // prefix sum up to the base before the BP acceptor window
    let csum_at_acc = at(&d.csum_f, pi.wrapping_sub(2).wrapping_sub(b.acc_bp[0].margin));
    let acc_term = d.arr[ri];
    let mut best = NEG;
    let mut arg = -1;
    for s in s0..hi {
        // (int arithmetic of the original; the spacer range comes from the
        // model's BP_ACC MIN/MAX and wraps only for absurd values)
        let k = pi.wrapping_sub(s).wrapping_sub(1);
        let v = (((at(&d.fwd, k) + lp(&b.bp_acc, s)) + lp(&b.don_bp, dist.wrapping_sub(1).wrapping_sub(w).wrapping_sub(s))) + run) + (csum_at_acc - at(&d.csum_f, k));
        let v = v + acc_term;
        if best < v {
            best = v;
            arg = s;
        }
    }
    (best, arg)
}

/// 0x41fec0 (Viterbi) / 0x420070 (-b): reverse intron, acceptor- at `pj`
/// (record `rj`), donor- at `pi`.
///
/// Per spacer `s` (`k = pj + s`; the scan stops once `k` reaches L):
/// `(csum_r[k] - csum_r[M0+pj+1]) + (((rev[k+1] + lp(BP_ACC, s)) +
/// lp(DON_BP, D-1-W-s)) + run)`, then `+ arr[rj]` (the acceptor- record,
/// not the donor). The bracketing is that of the original.
pub fn cand_rev(d: &BpData, b: &BranchPoint, run: f64, pj: i32, pi: i32, rj: usize) -> (f64, i32) {
    let Some((s0, hi)) = spacer_range(d, b, pj, pi) else {
        return (NEG, -1);
    };
    let w = b.branch.width;
    let dist = pi - pj;
    let csum_at_acc = at(&d.csum_r, b.acc_bp[0].margin.wrapping_add(pj).wrapping_add(1));
    let acc_term = d.arr[rj];
    let mut best = NEG;
    let mut arg = -1;
    for s in s0..hi {
        if d.l <= pj.wrapping_add(s) {
            break;
        }
        let k = pj.wrapping_add(s);
        let a = ((at(&d.rev, k.wrapping_add(1)) + lp(&b.bp_acc, s)) + lp(&b.don_bp, dist.wrapping_sub(1).wrapping_sub(w).wrapping_sub(s))) + run;
        let v = ((at(&d.csum_r, k) - csum_at_acc) + a) + acc_term;
        if best < v {
            best = v;
            arg = s;
        }
    }
    (best, arg)
}

/// -b writer 0x418900. `intron` = IN track, `tr` = output format is tr.
/// Returns the text to write; the caller opens the file with "w", or "a"
/// for tr (where the text is always empty, see the comment in the loop).
/// bp_off / bp_on / spacer persist across introns like the binary's locals.
///
/// Per intron of the best parse: `bp_off` = the normal intron score
/// (content + length), recomputed by scanning back from the acceptor-side
/// record to the donor-side record like the Viterbi does (no cut-off and no
/// stop check, which cannot matter on a traced path); `bp_on` = the BP
/// candidate; status 'U' if the BP score is higher, else 'D'. A border
/// intron gets `bp_on = NEG` but keeps the spacer of the previous intron
/// (stale local in the original). `bp_on = NEG` prints as the full 121-digit
/// `%.2f` expansion of -1e120, as in C.
pub fn intron_stats(
    d: &BpData,
    m: &Model,
    map: &[Rec],
    intron: &[f64],
    recs: &[StateRec],
    tr: bool,
) -> String {
    let b = m.bp.as_ref().unwrap();
    let k = m.order;
    let mut o = String::new();
    if !tr {
        o.push_str("status\t#\tgeneId\tstrand\tleft\tright\tlength\tbp_off\tbp_on\tspacer\n");
    }
    let mut count = 0;
    let mut bp_off = 0.0f64;
    let mut bp_on = 0.0f64;
    let mut spacer = -1i32;
    for r in recs {
        // the tr writer has already stripped the type names in place, so
        // the comparison with "Intron    " never matches there
        if r.kind != Kind::Intron || tr {
            continue;
        }
        let e = r.eidx;
        let j0 = r.sidx;
        let init = intron[e];
        let mut run = 0.0;
        count += 1;
        let h = map[e].sub;
        let pe = map[e].pos;
        let (full, border) = if r.strand == b'+' { (4, 0x104) } else { (0x80, 0x180) };
        // rebuild `run` exactly as the Viterbi did, down to j0
        for j in (j0..e).rev() {
            if map[j].sub != h {
                continue;
            }
            let l = pe - map[j].pos;
            if map[j].typ == full {
                bp_off = NEG;
                if k < l - 3 {
                    run = if run == 0.0 { init } else { run + intron[j] };
                    bp_off = lp(&m.intron, l) + run;
                }
            } else if map[j].typ == border {
                bp_off = NEG;
                if k < l - 1 {
                    run = if run == 0.0 { init } else { run + intron[j] };
                    let t = if m.intron.tail_lo() < l && l < m.intron.tail_hi() { m.intron.tail_base(l) } else { NEG };
                    bp_off = run + t;
                }
            }
        }
        let l = pe - map[j0].pos;
        if map[j0].typ == full {
            if k < l - 3 {
                let (v, a) = if r.strand == b'+' {
                    cand_fwd(d, b, run, map[j0].pos, pe, e)
                } else {
                    cand_rev(d, b, run, map[j0].pos, pe, j0)
                };
                bp_on = v;
                spacer = a;
            } else {
                bp_on = NEG;
            }
        } else if map[j0].typ == border {
            bp_on = NEG;
        } else {
            crate::out("error in code\n");
        }
        let status = if bp_off < bp_on { 'U' } else { 'D' };
        let end1 = r.end + 1;
        let _ = write!(
            o,
            "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{:.2}\t{:.2}\t{}\n",
            status,
            count,
            r.gene,
            r.strand as char,
            r.start + 1,
            end1,
            end1 - r.start,
            bp_off,
            bp_on,
            spacer + 1
        );
        let _ = write!(o, "{}\t{}\n", j0, e);
    }
    o
}

/// Upper-bound pruning of the BP candidates in the Viterbi (speed only).
///
/// The Viterbi only uses a BP candidate if `cand + cell[j]` beats the
/// running best (strict `>`). These bounds let it skip candidates that
/// provably cannot: every term of the candidate sum is bounded by a
/// precomputed maximum, and a safety margin far above the possible rounding
/// error is added. If any term is huge (|x| > 1e100, i.e. NEG-like
/// sentinels), no bound is given and the candidate is always evaluated.
///
/// Why this is exact: a skipped candidate has computed value `v <= ub <
/// best` (the margin covers the rounding of the few additions in both `v`
/// and the bound), so the strict update `best < v` would have rejected it
/// anyway and no cell or back pointer can differ. `GMHMM_RS_NO_PRUNE=1`
/// disables the pruning for checking this claim.
pub struct Bounds {
    /// per map record (acceptor+ i): max_s fwd[k] - csum_f[k], k = P[i]-s-1
    pub fw: Vec<f64>,
    /// per map record (acceptor- j): max_s rev[k+1] + csum_r[k], k = P[j]+s
    pub rw: Vec<f64>,
    /// max of lp(BP_ACC, s) over the spacer range
    pub bmax: f64,
    /// max of the DON_BP log table
    pub dmax: f64,
}

/// Terms with a magnitude above this are NEG-like; no bound is attempted.
const HUGE: f64 = 1e100;

/// Precompute the per-record maxima over the widest possible spacer range
/// `acc_min1 .. acc_max` (each candidate's actual range is a subset).
pub fn bounds(d: &BpData, b: &BranchPoint, map: &[Rec]) -> Bounds {
    let s0 = d.acc_min1;
    let s1 = d.acc_max;
    let mut bmax = f64::NEG_INFINITY;
    for s in s0..s1 {
        bmax = bmax.max(lp(&b.bp_acc, s));
    }
    let dmax = b.don_bp.tab.iter().copied().fold(f64::NEG_INFINITY, f64::max);
    let mut fw = vec![f64::INFINITY; map.len()];
    let mut rw = vec![f64::INFINITY; map.len()];
    for (r, rec) in map.iter().enumerate() {
        let p = rec.pos;
        if rec.typ == 8 {
            let mut m = f64::NEG_INFINITY;
            for s in s0..s1 {
                let k = p.wrapping_sub(s).wrapping_sub(1);
                m = m.max(at(&d.fwd, k) - at(&d.csum_f, k));
            }
            fw[r] = m;
        } else if rec.typ == 0x80 {
            let mut m = f64::NEG_INFINITY;
            for s in s0..s1 {
                if d.l <= p.wrapping_add(s) {
                    break;
                }
                let k = p.wrapping_add(s);
                m = m.max(at(&d.rev, k.wrapping_add(1)) + at(&d.csum_r, k));
            }
            rw[r] = m;
        }
    }
    Bounds { fw, rw, bmax, dmax }
}

/// Upper bound of a sum of terms, or +inf if any term is huge / not finite.
#[inline]
fn safe_sum(terms: &[f64]) -> f64 {
    let mut s = 0.0;
    let mut a = 0.0;
    for &t in terms {
        if !(t.abs() < HUGE) {
            return f64::INFINITY;
        }
        s += t;
        a += t.abs();
    }
    // rounding of a handful of additions is < 1e-15 * a; use a wide margin
    s + a * 1e-9 + 1e-6
}

/// Upper bound of `cand_fwd(..).0 + cell_j` (acceptor+ record `ri`).
#[inline]
pub fn bound_fwd(bd: &Bounds, d: &BpData, b: &BranchPoint, run: f64, pi: i32, ri: usize, cell_j: f64) -> f64 {
    let d1 = at(&d.csum_f, pi.wrapping_sub(2).wrapping_sub(b.acc_bp[0].margin));
    safe_sum(&[bd.fw[ri], bd.bmax, bd.dmax, run, d1, d.arr[ri], cell_j])
}

/// Upper bound of `cand_rev(..).0 + cell_j` (acceptor- record `rj` at `pj`).
#[inline]
pub fn bound_rev(bd: &Bounds, d: &BpData, b: &BranchPoint, run: f64, pj: i32, rj: usize, cell_j: f64) -> f64 {
    let base = at(&d.csum_r, b.acc_bp[0].margin.wrapping_add(pj).wrapping_add(1));
    safe_sum(&[bd.rw[rj], -base, bd.bmax, bd.dmax, run, d.arr[rj], cell_j])
}
