//! Segment Viterbi over map records (0x415220 init, 0x41ef70 driver,
//! helpers 0x41d180..0x41edf0). See re-notes/viterbi.md sections 2 and 4.
//! Addition order follows the binary exactly (bit-exact cells).
//!
//! # The DP
//!
//! The HMM is a generalized (explicit-duration) HMM, and the original runs
//! it as a segment DP over the site map: one cell per map record, records in
//! map order. `cell[i]` is the best log score of a parse of the sequence
//! start up to the signal `i`, where the last segment ends at `i`. The state
//! entered after a signal is implied by its type (after donor+ k: intron+
//! of phase k; after ATG+: exon+; after stop+: intergenic, ...), so a
//! segment `j -> i` is fully determined by the two record types.
//!
//! For each record `i` (13 .. N-1) a type-specific helper scans the earlier
//! records `j` backwards (nearest first) and evaluates
//! `segment score(j -> i) + cell[j]`:
//! * content: from the score tracks, accumulated in `run` (below);
//! * duration: `lp(D, L)` of the segment's length distribution, L = P[i] -
//!   P[j] (map distance; each segment type maps it to the true length with
//!   its own offset, see re-notes/viterbi.md section 2);
//! * signals: site scores `SITE[j]` (and `SITE[i]` through `init`);
//! * gene structure: exon phase terms and the `To*` transition logs;
//! * for branch point models, a second intron candidate ([`crate::bp`]).
//!
//! Exon helpers only look at records of the same frame (so the content is
//! always in frame) and stop at the first in-frame stop codon (ORF rule).
//! Intron helpers only look at records of the same intron phase and skip
//! introns that would create a stop codon across the splice junction.
//! Every scan ends when the segment gets longer than the model's maximum.
//!
//! # Exactness rules
//!
//! * Updates are strict (`best < s`): the first (nearest) `j` wins ties.
//! * `run` is the running content sum: `run = (run == 0.0) ? init : run +
//!   track[j]`. The exact `== 0.0` test is the original's way to detect the
//!   first admissible `j`; it is kept as is.
//! * The bracketing of every sum is the original's (some helpers write
//!   `x + run` instead of `run + x`; a single addition is commutative in
//!   IEEE, only the tree shape matters, and that is preserved literally).
//! * Exons at the sequence borders score `cell[j] - 1e120`, which rounds to
//!   exactly NEG and can never win: genes cannot be cut by the sequence
//!   ends except through a partial intron.

use crate::bp::{self, BpData};
use crate::model::{BranchPoint, LenDist, Model, NEG};
use crate::sitemap::Rec;
use crate::tracks::Tracks;

/// -ln 4: per-base content score of an exon too short for the Markov chain
/// (uniform base distribution), used as `(L-3) * MINUS_LN4`.
const MINUS_LN4: f64 = -1.3862943611198906;

/// Result of the DP: back pointer and score per map record (the 0x18-byte
/// cells of the original, whose third field is never read).
pub struct Cells {
    /// index of the left record `j` of the best last segment (0 if none)
    pub bp: Vec<i32>,
    /// best log score; NEG if the record is unreachable
    pub score: Vec<f64>,
}

/// Duration log-probability `lp(D, L)`: `tab_base[L]` if `lo < L < hi`
/// (strict on both sides), else NEG.
#[inline]
fn lp(d: &LenDist, l: i32) -> f64 {
    if d.lo() < l && l < d.hi() {
        d.base(l)
    } else {
        NEG
    }
}

/// Duration log-probability of a border (truncated) intron, from the tail
/// table: `tail_base[L]` if `tail_lo < L < tail_hi`, else NEG.
#[inline]
fn lp2(d: &LenDist, l: i32) -> f64 {
    if d.tail_lo() < l && l < d.tail_hi() {
        d.tail_base(l)
    } else {
        NEG
    }
}

/// DP state shared by all helpers.
struct Dp<'a> {
    map: &'a [Rec],
    tracks: &'a Tracks,
    m: &'a Model,
    codes: &'a [u8],
    /// Markov order K: content needs more than K bases
    k: i32,
    /// back pointers (output)
    back: Vec<i32>,
    /// cell scores (output)
    score: Vec<f64>,
    /// running maximum of the current helper call; reset to NEG per record
    best: f64,
    /// TAA_ON, TAG_ON, TGA_ON (for the stop-across-splice check)
    stop_on: [bool; 3],
    /// branch point data and model, only for $BP models
    bpd: Option<(&'a BpData, &'a BranchPoint)>,
    /// pruning bounds for the BP candidates (None: no pruning)
    bounds: Option<bp::Bounds>,
}

impl Dp<'_> {
    /// `upd(i, j, s)`: strict update of cell `i` with candidate `s` from `j`.
    #[inline]
    fn upd(&mut self, i: usize, j: usize, s: f64) {
        if self.best < s {
            self.best = s;
            self.back[i] = j as i32;
            self.score[i] = s;
        }
    }
    /// Is codon `16*c0 + 4*c1 + c2` (forward reading) an enabled stop?
    /// 0x30 TAA, 0x32 TAG, 0x38 TGA.
    fn is_stop_fwd(&self, c: usize) -> bool {
        (c == 0x30 && self.stop_on[0]) || (c == 0x32 && self.stop_on[1]) || (c == 0x38 && self.stop_on[2])
    }
    /// Same for a codon read on the forward strand that is a stop on the
    /// reverse strand: 0x3c TTA (TAA), 0x1c CTA (TAG), 0x34 TCA (TGA).
    fn is_stop_rev(&self, c: usize) -> bool {
        (c == 0x3c && self.stop_on[0]) || (c == 0x1c && self.stop_on[1]) || (c == 0x34 && self.stop_on[2])
    }
    /// Base code; reads past the sequence end return 0 (A), which is what
    /// the binary's heap buffer holds there (out-of-bounds read in 0x41d490).
    #[inline]
    fn code_at(&self, i: i32) -> usize {
        self.codes.get(i as usize).copied().unwrap_or(0) as usize
    }

    /// 0x41d180: donor+ record `i` ends an exon+. Left records of the same
    /// frame: acceptor+ (internal exon), ATG+ (initial exon), 0x108 (exon
    /// from the sequence start, impossible). An in-frame stop+ ends the scan.
    /// Cut-off: L >= internal-exon MAX.
    fn donor_fwd(&mut self, i: usize, js: impl Iterator<Item = usize>) {
        let m = self.m;
        let (frame, phase_i) = (self.map[i].frame, self.map[i].sub as usize);
        let init = self.tracks.cod[i] + self.tracks.site[i];
        let mut run = 0.0;
        self.best = NEG;
        for j in js {
            if self.map[j].frame != frame {
                continue;
            }
            let l = self.map[i].pos - self.map[j].pos;
            if m.internal.max <= l {
                return;
            }
            match self.map[j].typ {
                8 => {
                    // internal exon acceptor+ -> donor+
                    let mut s = NEG;
                    if self.k < l - 1 {
                        run = if run == 0.0 { init } else { run + self.tracks.cod[j] };
                        s = ((run + lp(&m.internal, l)) + self.tracks.site[j]) + m.phase_internal[phase_i + 3 * self.map[j].sub as usize];
                    }
                    let s = (s + m.log_to_internal) + self.score[j];
                    self.upd(i, j, s);
                }
                2 => return,
                1 => {
                    // initial exon ATG+ -> donor+; phase term indexed by the
                    // left record's sub (always 0 for a codon record)
                    let s;
                    if self.k < l - 3 {
                        run = if run == 0.0 { init } else { run + self.tracks.cod[j] };
                        s = ((run + lp(&m.initial, l)) + self.tracks.site[j]) + m.phase_initial[self.map[j].sub as usize];
                    } else {
                        // too short for the chain: uniform content, and
                        // SITE[i] added explicitly (it is not in a run)
                        let mut x = NEG;
                        if l > 2 {
                            x = (((lp(&m.initial, l) + self.tracks.site[j]) + m.phase_initial[self.map[j].sub as usize]) + (l - 3) as f64 * MINUS_LN4) + self.tracks.site[i];
                        }
                        s = x;
                    }
                    let s = (s + m.log_to_multi) + self.score[j];
                    self.upd(i, j, s);
                }
                0x108 => {
                    let s = self.score[j] - 1e120;
                    self.upd(i, j, s);
                }
                _ => {}
            }
        }
    }

    /// 0x41d490 (fwd = true, acceptor+ `i`) / 0x41d870 (fwd = false,
    /// donor- `i`): `i` ends an intron. Left records of the same phase: the
    /// opening splice site (full intron) or the head sentinel 0x104/0x180
    /// (intron from the sequence start, tail-table duration). Cut-off:
    /// L >= INTRON MAX.
    ///
    /// For a full intron of phase 1 or 2 the split codon (bases before the
    /// donor + bases after the acceptor) must not be an enabled stop; if it
    /// is, `j` is skipped after `run` has been updated (so later `j` still
    /// see the right content sum).
    fn intron_end(&mut self, i: usize, fwd: bool, js: impl Iterator<Item = usize>) {
        let m = self.m;
        let phase = self.map[i].sub;
        let init = self.tracks.intron[i];
        let mut run = 0.0;
        let (t_full, t_border) = if fwd { (4, 0x104) } else { (0x80, 0x180) };
        self.best = NEG;
        for j in js {
            if self.map[j].sub != phase {
                continue;
            }
            let (pj, pi) = (self.map[j].pos, self.map[i].pos);
            let l = pi - pj;
            if m.intron.max <= l {
                return;
            }
            let t = self.map[j].typ;
            if t == t_full {
                let mut s = NEG;
                if self.k < l - 3 {
                    run = if run == 0.0 { init } else { run + self.tracks.intron[j] };
                    s = run + lp(&m.intron, l);
                }
                // stop codon across the splice junction (forward reading
                // c0*16 + c1*4 + c2); bounds tests as in the original
                if fwd {
                    if phase == 1 && pj >= 1 {
                        if self.is_stop_fwd(self.code_at(pj - 1) * 16 + self.code_at(pi + 1) * 4 + self.code_at(pi + 2)) {
                            continue;
                        }
                    } else if phase == 2 && pj > 1 && self.is_stop_fwd(self.code_at(pj - 2) * 16 + self.code_at(pj - 1) * 4 + self.code_at(pi + 1)) {
                        continue;
                    }
                } else if phase == 1 && pj >= 2 {
                    if self.is_stop_rev(self.code_at(pj - 2) * 16 + self.code_at(pj - 1) * 4 + self.code_at(pi + 1)) {
                        continue;
                    }
                } else if phase == 2 && pj > 0 && self.is_stop_rev(self.code_at(pj - 1) * 16 + self.code_at(pi + 1) * 4 + self.code_at(pi + 2)) {
                    continue;
                }
                let s = s + self.score[j];
                self.upd(i, j, s);
                // BP candidate (0x41fad0 / 0x41fec0); NEG + cell[j] can never win
                if let Some((d, b)) = self.bpd {
                    if self.k < l - 3 {
                        // skip candidates that provably cannot beat `best`
                        // (speed only; see bp::Bounds for why this is exact)
                        if let Some(bd) = &self.bounds {
                            let ub = if fwd { bp::bound_fwd(bd, d, b, run, pi, i, self.score[j]) } else { bp::bound_rev(bd, d, b, run, pj, j, self.score[j]) };
                            if ub < self.best {
                                continue;
                            }
                        }
                        // the reverse candidate takes record j (the
                        // acceptor-), the forward one record i (acceptor+)
                        let v = if fwd { bp::cand_fwd(d, b, run, pj, pi, i).0 } else { bp::cand_rev(d, b, run, pj, pi, j).0 };
                        let v = v + self.score[j];
                        self.upd(i, j, v);
                    }
                }
            } else if t == t_border {
                // intron from the sequence start: no stop check, no BP
                // candidate (0x41fea0 / 0x41feb0 return NEG)
                let mut s = NEG;
                if self.k < l - 1 {
                    run = if run == 0.0 { init } else { run + self.tracks.intron[j] };
                    s = run + lp2(&m.intron, l);
                }
                let s = s + self.score[j];
                self.upd(i, j, s);
            }
        }
    }

    /// 0x41dcb0: acceptor- record `i` ends (in forward coordinates) an
    /// exon-. Left records of the same frame: donor- (internal exon),
    /// stop- (terminal exon of the reverse gene; evaluated, then the scan
    /// ends), 0x140 (impossible). Cut-off: L >= internal-exon MAX.
    fn acceptor_rev(&mut self, i: usize, js: impl Iterator<Item = usize>) {
        let m = self.m;
        let (frame, phase_i) = (self.map[i].frame, self.map[i].sub as usize);
        let init = self.tracks.cod[i] + self.tracks.site[i];
        let mut run = 0.0;
        self.best = NEG;
        for j in js {
            if self.map[j].frame != frame {
                continue;
            }
            let l = self.map[i].pos - self.map[j].pos;
            if m.internal.max <= l {
                return;
            }
            match self.map[j].typ {
                0x40 => {
                    let mut s = NEG;
                    if self.k < l - 1 {
                        run = if run == 0.0 { init } else { run + self.tracks.cod[j] };
                        s = ((run + lp(&m.internal, l)) + self.tracks.site[j]) + m.phase_internal[phase_i + 3 * self.map[j].sub as usize];
                    }
                    let s = (s + m.log_to_internal) + self.score[j];
                    self.upd(i, j, s);
                }
                0x20 => {
                    // terminal exon of a reverse gene; leftmost exon of a
                    // multi-exon gene in forward order, hence ToMultiGene
                    let ph = self.map[j].sub as usize;
                    let mut s = NEG;
                    if self.k < l - 3 {
                        run = if run == 0.0 { init } else { run + self.tracks.cod[j] };
                        s = ((run + lp(&m.terminal, l)) + self.tracks.site[j]) + m.phase_terminal[ph];
                    } else if l > 2 {
                        // note the bracketing differs from the forward helpers
                        s = ((m.phase_terminal[ph] + (lp(&m.terminal, l) + self.tracks.site[j])) + (l - 3) as f64 * MINUS_LN4) + self.tracks.site[i];
                    }
                    let s = (s + m.log_to_multi) + self.score[j];
                    self.upd(i, j, s);
                    return;
                }
                0x140 => {
                    let s = self.score[j] - 1e120;
                    self.upd(i, j, s);
                }
                _ => {}
            }
        }
    }

    /// 0x41dfc0: stop+ record `i` ends an exon+. Left records of the same
    /// frame: acceptor+ (terminal exon), ATG+ (single-exon gene), 0x108
    /// (impossible). An in-frame stop+ ends the scan.
    fn stop_fwd(&mut self, i: usize, js: impl Iterator<Item = usize>) {
        let m = self.m;
        let frame = self.map[i].frame;
        let init = self.tracks.cod[i] + self.tracks.site[i];
        let mut run = 0.0;
        self.best = NEG;
        for j in js {
            if self.map[j].frame != frame {
                continue;
            }
            let l = self.map[i].pos - self.map[j].pos;
            if m.internal.max <= l {
                return;
            }
            match self.map[j].typ {
                8 => {
                    let ph = self.map[j].sub as usize;
                    let s;
                    if self.k < l - 3 {
                        run = if run == 0.0 { init } else { run + self.tracks.cod[j] };
                        s = ((run + lp(&m.terminal, l)) + self.tracks.site[j]) + m.phase_terminal[ph];
                    } else {
                        let mut x = NEG;
                        if l > 2 {
                            x = (((lp(&m.terminal, l) + self.tracks.site[j]) + m.phase_terminal[ph]) + (l - 3) as f64 * MINUS_LN4) + self.tracks.site[i];
                        }
                        s = x;
                    }
                    let s = (s + m.log_to_terminal) + self.score[j];
                    self.upd(i, j, s);
                }
                2 => return,
                1 => {
                    // single-exon gene: needs more than K+5 bases, no
                    // short-exon fallback, no phase term
                    let mut s = NEG;
                    if self.k < l - 5 {
                        run = if run == 0.0 { init } else { run + self.tracks.cod[j] };
                        s = (run + lp(&m.single, l)) + self.tracks.site[j];
                    }
                    let s = (s + m.log_to_single) + self.score[j];
                    self.upd(i, j, s);
                }
                0x108 => {
                    let s = self.score[j] - 1e120;
                    self.upd(i, j, s);
                }
                _ => {}
            }
        }
    }

    /// 0x41e470: start- record `i` (CAT) ends an exon-. Left records of the
    /// same frame: donor- (initial exon of the reverse gene; rightmost in
    /// forward order, hence ToTerminalExon), stop- (single-exon gene;
    /// evaluated, then the scan ends), 0x140 (impossible).
    fn start_rev(&mut self, i: usize, js: impl Iterator<Item = usize>) {
        let m = self.m;
        let frame = self.map[i].frame;
        let init = self.tracks.cod[i] + self.tracks.site[i];
        let mut run = 0.0;
        self.best = NEG;
        for j in js {
            if self.map[j].frame != frame {
                continue;
            }
            let l = self.map[i].pos - self.map[j].pos;
            if m.internal.max <= l {
                return;
            }
            match self.map[j].typ {
                0x40 => {
                    let ph = self.map[j].sub as usize;
                    let s;
                    if self.k < l - 3 {
                        run = if run == 0.0 { init } else { run + self.tracks.cod[j] };
                        s = ((run + lp(&m.initial, l)) + self.tracks.site[j]) + m.phase_initial[ph];
                    } else {
                        let mut x = NEG;
                        if l > 2 {
                            x = (((lp(&m.initial, l) + self.tracks.site[j]) + m.phase_initial[ph]) + (l - 3) as f64 * MINUS_LN4) + self.tracks.site[i];
                        }
                        s = x;
                    }
                    let s = (s + m.log_to_terminal) + self.score[j];
                    self.upd(i, j, s);
                }
                0x20 => {
                    let mut s = NEG;
                    if self.k < l - 5 {
                        run = if run == 0.0 { init } else { run + self.tracks.cod[j] };
                        // bracketing of the original (differs from stop_fwd)
                        s = self.tracks.site[j] + (lp(&m.single, l) + run);
                    }
                    let s = (s + m.log_to_single) + self.score[j];
                    self.upd(i, j, s);
                    return;
                }
                0x140 => {
                    let s = self.score[j] - 1e120;
                    self.upd(i, j, s);
                }
                _ => {}
            }
        }
    }

    /// 0x41e2b0 (at_end = false: `i` is ATG+ or stop-, a gene starts) /
    /// 0x41e760 (at_end = true: `i` is the 0x220 tail sentinel): `i` ends an
    /// intergenic segment. Left records (no filter): stop+ / start- (a gene
    /// ended) or the 0x102 head sentinel. Cut-off: `L - 1 > INTERGENIC MAX`.
    ///
    /// Segments touching a sequence border (from 0x102 or to 0x220) get a
    /// +100 bonus (once even if both apply). Both border cases use
    /// P(len = L-1) although the true length is L (as in the original).
    fn ig_end(&mut self, i: usize, at_end: bool, js: impl Iterator<Item = usize>) {
        let m = self.m;
        let init = self.tracks.ig[i];
        let mut run = 0.0;
        let pi = self.map[i].pos;
        self.best = NEG;
        for j in js {
            let l = pi - self.map[j].pos;
            if m.intergenic.max < l - 1 {
                return;
            }
            let t = self.map[j].typ;
            let s;
            if t == 2 || t == 0x10 {
                let mut x = NEG;
                if self.k < l - 1 {
                    run = if run == 0.0 { init } else { self.tracks.ig[j] + run };
                    x = lp(&m.intergenic, l) + run;
                    if at_end {
                        x += 100.0;
                    }
                }
                s = x;
            } else if t == 0x102 {
                let mut x = NEG;
                if self.k < l - 1 {
                    run = if run == 0.0 { init } else { run + self.tracks.ig[j] };
                    x = (run + lp(&m.intergenic, l)) + 100.0;
                }
                s = x;
            } else {
                continue;
            }
            let s = s + self.score[j];
            self.upd(i, j, s);
        }
    }

    /// 0x41e940 (0x208, t1 = 4, t2 = 0x104) / 0x41eae0 (0x240, t1 = 0x80,
    /// t2 = 0x180): an intron of `i`'s phase runs into the sequence end.
    /// Tail-table duration; no stop check, no BP candidate.
    fn intron_to_end(&mut self, i: usize, t1: i32, t2: i32, js: impl Iterator<Item = usize>) {
        let m = self.m;
        let phase = self.map[i].sub;
        let init = self.tracks.intron[i];
        let mut run = 0.0;
        self.best = NEG;
        for j in js {
            if self.map[j].sub != phase {
                continue;
            }
            let l = self.map[i].pos - self.map[j].pos;
            if m.intron.max <= l {
                return;
            }
            let t = self.map[j].typ;
            if t != t1 && t != t2 {
                continue;
            }
            let mut s = NEG;
            if self.k < l - 1 {
                // same values in both branches; only the bracketing of the
                // original differs
                if t == t1 {
                    run = if run == 0.0 { init } else { self.tracks.intron[j] + run };
                    s = lp2(&m.intron, l) + run;
                } else {
                    run = if run == 0.0 { init } else { run + self.tracks.intron[j] };
                    s = run + lp2(&m.intron, l);
                }
            }
            let s = s + self.score[j];
            self.upd(i, j, s);
        }
    }

    /// 0x41ec80 (fwd, 0x204) / 0x41edf0 (reverse, 0x280): an exon running
    /// into the sequence end. Every candidate scores `cell[j] - 1e120`,
    /// which rounds to NEG for any ordinary `cell[j]` and so never passes
    /// the strict update: these cells stay NEG with back pointer 0 in
    /// practice. Kept literally for the unusual cases. Forward: an in-frame
    /// stop+ ends the scan; reverse: an in-frame stop- is the last
    /// candidate.
    fn exon_to_end(&mut self, i: usize, fwd: bool, js: impl Iterator<Item = usize>) {
        let m = self.m;
        let frame = self.map[i].frame;
        self.best = NEG;
        for j in js {
            if self.map[j].frame != frame {
                continue;
            }
            let l = self.map[i].pos - self.map[j].pos;
            if m.internal.max <= l {
                return;
            }
            let t = self.map[j].typ;
            if fwd {
                if t == 8 || t == 1 || t == 0x108 {
                    let s = self.score[j] - 1e120;
                    self.upd(i, j, s);
                } else if t == 2 {
                    return;
                }
            } else if t == 0x40 || t == 0x140 {
                let s = self.score[j] - 1e120;
                self.upd(i, j, s);
            } else if t == 0x20 {
                let s = self.score[j] - 1e120;
                self.upd(i, j, s);
                return;
            }
        }
    }
}

/// 0x415220 + 0x41ef70: run the DP over the whole map.
///
/// `bpd` is the branch point data (only for `$BP` models). With `verbose`
/// (`-v`) progress lines are printed to stdout. Returns all cells; the
/// traceback starts from the 13 tail-sentinel cells.
pub fn run(map: &[Rec], t: &Tracks, m: &Model, codes: &[u8], bpd: Option<&BpData>, verbose: bool) -> Cells {
    let n = map.len();
    let mut dp = Dp {
        map,
        tracks: t,
        m,
        codes,
        k: m.order,
        back: vec![0; n],
        score: vec![NEG; n],
        best: NEG,
        stop_on: [m.taa_on, m.tag_on, m.tga_on],
        bpd: bpd.zip(m.bp.as_ref()),
        bounds: match (bpd, m.bp.as_ref()) {
            (Some(d), Some(b)) if std::env::var_os("GMHMM_RS_NO_PRUNE").is_none() => Some(bp::bounds(d, b, map)),
            _ => None,
        },
    };
    // 0x415220: initial state log-probabilities of the head sentinels
    // (log_initial = [direct intron, reverse intron, intergenic, direct
    // coding, reverse coding]). There is no end probability.
    dp.score[0] = m.log_initial[2];
    for k in 1..=3 {
        dp.score[k] = m.log_initial[0];
    }
    for k in 4..=6 {
        dp.score[k] = m.log_initial[1];
    }
    for k in 7..=9 {
        dp.score[k] = m.log_initial[3];
    }
    for k in 10..=12 {
        dp.score[k] = m.log_initial[4];
    }
    if verbose {
        crate::out(&format!("map_points total: {n}\n"));
    }
    // Candidate lists: each helper only acts on a few record types (with a
    // frame or phase filter). Scanning only those records, nearest first,
    // performs the same floating-point operations in the same order as the
    // full backward scan. The early cut-offs (segment longer than MAX) are
    // also reached at the same point, because the map is sorted by position
    // and the distance only grows going backwards. Unsorted maps (never
    // produced by sitemap::build) fall back to the full scan.
    //
    // (In detail: the original tests the cut-off on every record that
    // passes the frame/phase filter, also on types the helper ignores. If
    // such a skipped record would have triggered it, every later candidate
    // is at least as far away and triggers it before doing anything else,
    // so the helper's observable effects are the same.)
    let sorted = map[13..n.saturating_sub(13).max(13)].windows(2).all(|w| w[0].pos <= w[1].pos);
    // exon+ / exon- left records per frame, intron+ / intron- left records
    // per phase, intergenic left records
    let mut ex_p: [Vec<usize>; 3] = Default::default();
    let mut ex_m: [Vec<usize>; 3] = Default::default();
    let mut in_p: [Vec<usize>; 3] = Default::default();
    let mut in_m: [Vec<usize>; 3] = Default::default();
    let mut ig: Vec<usize> = Vec::new();
    for (j, r) in map.iter().enumerate() {
        let (f, h) = (r.frame as usize % 3, r.sub as usize % 3);
        match r.typ {
            8 | 1 | 0x108 | 2 => ex_p[f].push(j),
            0x40 | 0x20 | 0x140 => ex_m[f].push(j),
            _ => {}
        }
        match r.typ {
            4 | 0x104 => in_p[h].push(j),
            0x80 | 0x180 => in_m[h].push(j),
            _ => {}
        }
        if matches!(r.typ, 2 | 0x10 | 0x102) {
            ig.push(j);
        }
    }
    // indices < i of a list, nearest first
    fn before(list: &[usize], i: usize) -> std::iter::Rev<std::iter::Copied<std::slice::Iter<'_, usize>>> {
        let k = list.partition_point(|&x| x < i);
        list[..k].iter().copied().rev()
    }
    // call a helper with the candidate list, or with all j < i (nearest
    // first) if the map is not sorted
    macro_rules! scan {
        ($list:expr, $i:expr, $call:expr) => {
            if sorted {
                $call(before($list, $i))
            } else {
                $call(Box::new((0..$i).rev()) as Box<dyn Iterator<Item = usize>>)
            }
        };
    }
    // 0x41ef70: dispatch on the exact record type
    for i in 13..n {
        if verbose && i % 10000 == 0 {
            crate::out(&format!("working: {} {:.1}\n", i, 100.0 * i as f64 / n as f64));
        }
        let f = map[i].frame as usize % 3;
        let h = map[i].sub as usize % 3;
        match map[i].typ {
            4 => scan!(&ex_p[f], i, |js| dp.donor_fwd(i, js)),
            8 => scan!(&in_p[h], i, |js| dp.intron_end(i, true, js)),
            0x40 => scan!(&in_m[h], i, |js| dp.intron_end(i, false, js)),
            0x80 => scan!(&ex_m[f], i, |js| dp.acceptor_rev(i, js)),
            2 => scan!(&ex_p[f], i, |js| dp.stop_fwd(i, js)),
            0x20 | 1 => scan!(&ig, i, |js| dp.ig_end(i, false, js)),
            0x10 => scan!(&ex_m[f], i, |js| dp.start_rev(i, js)),
            0x220 => scan!(&ig, i, |js| dp.ig_end(i, true, js)),
            0x208 => scan!(&in_p[h], i, |js| dp.intron_to_end(i, 4, 0x104, js)),
            0x240 => scan!(&in_m[h], i, |js| dp.intron_to_end(i, 0x80, 0x180, js)),
            0x204 => scan!(&ex_p[f], i, |js| dp.exon_to_end(i, true, js)),
            0x280 => scan!(&ex_m[f], i, |js| dp.exon_to_end(i, false, js)),
            _ => crate::out(&format!("error in map, position {}\n", map[i].pos)),
        }
    }
    Cells { bp: dp.back, score: dp.score }
}

/// 0x4153c0 (`-r`): the best path log-probability, the maximum over the 13
/// tail-sentinel cells, starting from `max(cell[N-13], NEG)`. Non-strict
/// `<=` here (unlike the traceback), which only matters for which of equal
/// values is kept, not for the value itself.
pub fn best_path(score: &[f64]) -> f64 {
    let n = score.len();
    let mut b = score[n - 13];
    if b <= NEG {
        b = NEG;
    }
    for &s in &score[n - 12..] {
        if b <= s {
            b = s;
        }
    }
    b
}
