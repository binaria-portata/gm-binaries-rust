//! Score tracks indexed by map record: intergenic (0x41aae0), intron
//! (0x41aa30, with 0x41a030 / 0x41a390), coding (0x419f60, with 0x4193e0,
//! 0x4197e0, 0x419110, 0x419230) and site (0x41abe0, with 0x41c1f0 /
//! 0x41c300). See re-notes/viterbi.md section 3.
//!
//! Stage between the site map and the Viterbi. Every array here has one f64
//! (log-probability) per map record, not per base.
//!
//! # Content scores and why they are split into "begin" and "end" entries
//!
//! The Viterbi needs, for a segment from record j to record i, the Markov
//! log-probability of the bases between them. Each segment class
//! (intergenic; intron per strand and phase; coding per strand and frame)
//! scores a half-open base range `[b, e)`: `b` from the left record, `e`
//! from the right record, with type-specific offsets (so that start/stop
//! codons and GT/AG are scored by the site models, not by the content
//! models; table in re-notes/viterbi.md 1.1).
//!
//! Instead of summing per segment, the original precomputes per class:
//! * end entry `track[E]` for every record E that can end a segment: the
//!   score of `[b(r), e(E))` where r is the nearest admissible begin record
//!   before E (the first one leaving more than K bases), including the
//!   final K+1 bases scored with the joint table;
//! * begin entry `track[r]` for every begin record r: the increment
//!   `INC(b(r), b(r'))` up to the next begin record r' of the same class
//!   above it.
//!
//! The Viterbi then rebuilds the content of `j -> i` as `track[i]` plus the
//! begin increments of the admissible begins it passes while scanning back
//! (its `run` accumulator). This grouping of the additions is part of the
//! bit-exact result: summing `[b(j), e(i))` directly differs in the last
//! bits, so the tracks must be built exactly as here.
//!
//! The non-coding chains read the sequence right to left: with `A` = the
//! "nf" table (P(first base | next K bases)) and `B` = the joint table,
//! `SEG(b, e) = B[kmer[e-K-1]] + A[kmer[b]] + ... + A[kmer[e-K-2]]`, summed
//! in exactly this order. The coding chain does the same with a ring of
//! three codon-position tables that advances by one per base.
//!
//! # Site scores
//!
//! [`site_track`] scores every candidate signal with its site matrix, plus
//! the GT/GC donor prior, the AROUND_BP term of branch point models and the
//! +/-100 evidence terms of `-d`.

use crate::model::{Model, Site, NEG};
use crate::sitemap::Rec;

/// The four per-record score arrays (main locals e68, e28, c28[0], c28[1]).
/// Records not written by a builder keep 0.0.
pub struct Tracks {
    /// intergenic content (MARKOV column 4)
    pub ig: Vec<f64>,
    /// intron content, both strands, all phases (MARKOV column 3 and its
    /// reverse complement)
    pub intron: Vec<f64>,
    /// coding content, both strands, all frames (MARKOV columns 0..2),
    /// including the `-k` mask penalty
    pub cod: Vec<f64>,
    /// site scores (log-odds or log-likelihood of the signal, priors and
    /// evidence terms)
    pub site: Vec<f64>,
}

/// Read-only context of the content-track builders (score context in
/// main).
struct Ctx<'a> {
    map: &'a [Rec],
    /// per-base k-mer index (sequence::kmer_index)
    kmer: &'a [i32],
    /// Markov order K
    k: i64,
}

impl Ctx<'_> {
    fn n(&self) -> usize {
        self.map.len()
    }
    /// Class key of record r: frame for coding classes, intron phase
    /// (`sub`) for intron and intergenic classes (intergenic only uses 0).
    #[inline]
    fn key(&self, r: usize, by_frame: bool) -> u8 {
        if by_frame {
            self.map[r].frame
        } else {
            self.map[r].sub
        }
    }
    /// Does record r belong to the class (type bit in `mask`, same key)?
    /// Sentinels carry site bits too (0x104 has bit 4, ...), so they match
    /// by mask.
    #[inline]
    fn matches(&self, r: usize, mask: i32, key: u8, by_frame: bool) -> bool {
        (self.map[r].typ & mask) != 0 && self.key(r, by_frame) == key
    }
    /// The head sentinel of the class: the matching record with the highest
    /// index among 12..=0.
    fn sentinel(&self, mask: i32, key: u8, by_frame: bool) -> Option<usize> {
        (0..=12).rev().find(|&r| self.matches(r, mask, key, by_frame))
    }
    #[inline]
    fn pos(&self, r: usize) -> i64 {
        self.map[r].pos as i64
    }
    /// k-mer index at base i (table row for the window c[i..=i+K]).
    #[inline]
    fn kmer_at(&self, i: i64) -> usize {
        self.kmer[i as usize] as usize
    }
}

/// 0x41a030: begin increments for one intergenic / intron class.
///
/// `nf` is the class's A table, `mask` its begin types, `begin_off` the
/// offset from a begin record's position to its first content base.
/// Walking the normal records top-down, each matching record r (except the
/// topmost one, which stays 0) gets `sum A[kmer[i]]` over
/// `[b(r), b(previous match))`, summed from 0.0 left to right. Finally the
/// head sentinel gets the increment from base 0 (its b is its position,
/// without offset) to the lowest matching record. Nothing at all is
/// written if no normal record matches.
fn begin_noncoding(c: &Ctx, trk: &mut [f64], nf: &[f64], key: u8, mask: i32, begin_off: i64) {
    let n = c.n();
    let mut prev: Option<usize> = None;
    for r in (13..=n - 14).rev() {
        if c.matches(r, mask, key, false) {
            if let Some(p) = prev {
                let mut acc = 0.0;
                for i in c.pos(r) + begin_off..c.pos(p) + begin_off {
                    acc += nf[c.kmer_at(i)];
                }
                trk[r] = acc;
            }
            prev = Some(r);
        }
    }
    let Some(p) = prev else { return };
    let rs = c.sentinel(mask, key, false).expect("sentinel");
    let mut acc = 0.0;
    for i in c.pos(rs)..c.pos(p) + begin_off {
        acc += nf[c.kmer_at(i)];
    }
    trk[rs] = acc;
}

/// `SEG(nf, joint, beg, endx)`: content score of bases `[beg, endx)` read
/// right to left: start with the joint probability of the last K+1 bases,
/// then add the conditional terms from left to right. The caller ensures
/// `endx - beg > K`.
fn seg_score(c: &Ctx, nf: &[f64], joint: &[f64], beg: i64, endx: i64) -> f64 {
    let e1 = endx - c.k - 1;
    let mut acc = joint[c.kmer_at(e1)];
    for i in beg..e1 {
        acc += nf[c.kmer_at(i)];
    }
    acc
}

/// 0x41a390: end entries for one intergenic / intron class.
///
/// `begin_mask`/`end_mask` are the class's begin and end types,
/// `begin_off`/`end_off` the offsets from record position to `b` / `e`.
/// End records: the class's tail sentinel (highest matching index among the
/// last 13, `e = P + 1 = L`) and every matching normal record. For each end,
/// the nearest matching begin below it that leaves more than K bases is
/// used; failing that, the head sentinel (b = 0). If even that is too short
/// the entry is NEG.
#[allow(clippy::too_many_arguments)]
fn end_noncoding(c: &Ctx, trk: &mut [f64], nf: &[f64], joint: &[f64], key: u8, begin_mask: i32, end_mask: i32, begin_off: i64, end_off: i64) {
    let n = c.n();
    let mut ends: Vec<(usize, i64)> = Vec::new();
    if let Some(top) = (n - 13..n).rev().find(|&r| c.matches(r, end_mask, key, false)) {
        ends.push((top, c.pos(top) + 1));
    }
    for r in (13..=n - 14).rev() {
        if c.matches(r, end_mask, key, false) {
            ends.push((r, c.pos(r) + end_off));
        }
    }
    // each end writes only its own entry and reads no track value, so the
    // processing order is irrelevant
    for (ei, endx) in ends {
        let mut val = None;
        for r in (13..ei.min(n - 13)).rev() {
            if c.matches(r, begin_mask, key, false) {
                let bb = c.pos(r) + begin_off;
                if endx - bb > c.k {
                    val = Some(seg_score(c, nf, joint, bb, endx));
                    break;
                }
            }
        }
        let v = match val {
            Some(v) => v,
            None => {
                let rs = c.sentinel(begin_mask, key, false).expect("sentinel");
                let bb = c.pos(rs);
                if endx - bb > c.k {
                    seg_score(c, nf, joint, bb, endx)
                } else {
                    NEG
                }
            }
        };
        trk[ei] = v;
    }
}

/// Coding ring node (records F0..F2, Ra..Rc of model.md 3.3):
/// `{A = nf table, B = joint table, next}`.
struct Node<'a> {
    a: &'a [f64],
    b: &'a [f64],
    next: usize,
}

/// The two three-periodic rings of coding tables. Scoring a coding stretch
/// starts at a head chosen by the begin record's `sub` byte and moves to
/// `next` after every base.
struct Ring<'a> {
    nodes: Vec<Node<'a>>,
    head_fwd: [usize; 3],
    head_rev: [usize; 3],
}

/// Build the rings from the coding block (pointer layout of 0x40f960).
/// Forward ring F0 -> F1 -> F2 -> F0 over columns 0, 1, 2; reverse ring
/// Ra -> Rc -> Rb -> Ra over the reverse-complement columns given by
/// `Coding::rc_ring`. Forward heads by phase: F0, F1, F2; reverse heads:
/// Rc, Ra, Rb.
fn ring(m: &Model) -> Ring<'_> {
    let c = &m.coding;
    let [a, b, cc] = c.rc_ring();
    // 0..2: F0,F1,F2 ; 3: Ra, 4: Rb, 5: Rc
    let nodes = vec![
        Node { a: &c.fwd_nf[0], b: &c.fwd[0], next: 1 },
        Node { a: &c.fwd_nf[1], b: &c.fwd[1], next: 2 },
        Node { a: &c.fwd_nf[2], b: &c.fwd[2], next: 0 },
        Node { a: &c.rc_nf[a], b: &c.rc[a], next: 5 },
        Node { a: &c.rc_nf[b], b: &c.rc[b], next: 3 },
        Node { a: &c.rc_nf[cc], b: &c.rc[cc], next: 4 },
    ];
    Ring { nodes, head_fwd: [0, 1, 2], head_rev: [5, 3, 4] }
}

/// Context of the coding-track builders: rings plus the `-k` mask penalty.
struct CodingCtx<'a> {
    ring: Ring<'a>,
    /// evidence marks (with -d); the penalty applies to bases with MaskMark
    marks: Option<&'a [u32]>,
    /// `-k` value, added as given per masked base
    pen: f64,
}

impl CodingCtx<'_> {
    fn head(&self, fwd: bool, ph: u8) -> usize {
        if fwd {
            self.ring.head_fwd[ph as usize]
        } else {
            self.ring.head_rev[ph as usize]
        }
    }
    /// 0x419110: begin increment, `sum A` over `[st, en)` walking the ring
    /// from the head of phase `ph`, then `+ pen` per masked base in the same
    /// range (added after the Markov sum, in a second loop, as in the
    /// original). 0.0 for an empty range.
    fn inc(&self, c: &Ctx, st: i64, en: i64, ph: u8, fwd: bool) -> f64 {
        if en <= st {
            return 0.0;
        }
        let mut n = self.head(fwd, ph);
        let mut acc = 0.0;
        for i in st..en {
            let nd = &self.ring.nodes[n];
            acc += nd.a[c.kmer_at(i)];
            n = nd.next;
        }
        if let Some(mk) = self.marks {
            for i in st..en {
                if mk[i as usize] & 0x20000 != 0 {
                    acc += self.pen;
                }
            }
        }
        acc
    }
    /// 0x419230: end entry, the coding analogue of `seg_score`: A terms over
    /// `[st, en-K-1)` walking the ring, then the B (joint) term of the node
    /// reached, at `en-K-1`; then the mask penalty over `[st, en-K-1)`
    /// only (the last K+1 bases get none). Note the different summation
    /// order from the non-coding version: here the B term is added last.
    fn end(&self, c: &Ctx, st: i64, en: i64, ph: u8, fwd: bool) -> f64 {
        let e1 = en - (c.k + 1);
        let mut n = self.head(fwd, ph);
        let mut acc = 0.0;
        let last;
        if st < e1 {
            for i in st..e1 {
                let nd = &self.ring.nodes[n];
                acc += nd.a[c.kmer_at(i)];
                n = nd.next;
            }
            last = e1;
        } else {
            // callers guarantee e1 >= st, so this is e1 == st
            last = st;
        }
        acc += self.ring.nodes[n].b[c.kmer_at(last)];
        if let Some(mk) = self.marks {
            if st < e1 {
                for i in st..e1 {
                    if mk[i as usize] & 0x20000 != 0 {
                        acc += self.pen;
                    }
                }
            }
        }
        acc
    }
}

/// 0x4193e0: coding begin increments for one strand and frame `f`.
///
/// Begin types are `m1 | m2` with offsets `o1` (record type in `m1`) or
/// `o2` from the record position to the first coding base (forward: ATG
/// +3, acceptor +1; reverse: stop +3, donor +1). Same scheme as
/// [`begin_noncoding`]; the ring head is chosen by the begin record's `sub`.
#[allow(clippy::too_many_arguments)]
fn begin_coding(c: &Ctx, cc: &CodingCtx, trk: &mut [f64], f: u8, m1: i32, m2: i32, o1: i64, o2: i64, fwd: bool) {
    let n = c.n();
    let m = m1 | m2;
    let mut prevb: Option<i64> = None;
    for r in (13..=n - 14).rev() {
        if c.matches(r, m, f, true) {
            let b = c.pos(r) + if c.map[r].typ & m1 != 0 { o1 } else { o2 };
            if let Some(pb) = prevb {
                trk[r] = cc.inc(c, b, pb, c.map[r].sub, fwd);
            }
            prevb = Some(b);
        }
    }
    let Some(pb) = prevb else { return };
    let rs = c.sentinel(m, f, true).expect("sentinel");
    trk[rs] = cc.inc(c, c.pos(rs), pb, c.map[rs].sub, fwd);
}

/// 0x4197e0: coding end entries for one strand and frame `f`.
///
/// Begin types `bm1 | bm2` with offsets `bo1`/`bo2`, end types `em1 | em2`
/// with offsets `eo1`/`eo2` (forward: stop -2, donor 0; reverse: start -2,
/// acceptor 0). As [`end_noncoding`], except that an end without any
/// admissible begin (or with `e < 0`) is left untouched (stays 0.0)
/// instead of being set to NEG.
#[allow(clippy::too_many_arguments)]
fn end_coding(c: &Ctx, cc: &CodingCtx, trk: &mut [f64], f: u8, bm1: i32, bm2: i32, bo1: i64, bo2: i64, em1: i32, em2: i32, eo1: i64, eo2: i64, fwd: bool) {
    let n = c.n();
    let bm = bm1 | bm2;
    let em = em1 | em2;
    let mut ends: Vec<(usize, i64)> = Vec::new();
    if let Some(top) = (n - 13..n).rev().find(|&r| c.matches(r, em, f, true)) {
        ends.push((top, c.pos(top) + 1));
    }
    for r in (13..=n - 14).rev() {
        if c.matches(r, em, f, true) {
            ends.push((r, c.pos(r) + if c.map[r].typ & em1 != 0 { eo1 } else { eo2 }));
        }
    }
    for (ei, endx) in ends {
        let valid = endx >= 0;
        let mut done = false;
        for r in (13..ei.min(n - 13)).rev() {
            if c.matches(r, bm, f, true) {
                let b = c.pos(r) + if c.map[r].typ & bm1 != 0 { bo1 } else { bo2 };
                // with !valid no begin is accepted, and the scan runs on to
                // the sentinel, which is rejected as well
                if endx - b > c.k && valid {
                    trk[ei] = cc.end(c, b, endx, c.map[r].sub, fwd);
                    done = true;
                    break;
                }
            }
        }
        if !done {
            if let Some(rs) = c.sentinel(bm, f, true) {
                let b = c.pos(rs);
                if endx - b > c.k && valid {
                    trk[ei] = cc.end(c, b, endx, c.map[rs].sub, fwd);
                }
            }
        }
    }
}

/// 0x41c1f0: forward-strand site score of `s` for the site record at
/// position `p`.
///
/// `sh` moves `p` to the first base of the consensus (0 for donor G and
/// start A, 1 for the acceptor G of AG, 2 for the last base of a stop, -1
/// for the BP "window ending at p"); the window then starts MARGIN bases
/// before it: `st = p - (sh + MARGIN)`. The sum over the WIDTH columns
/// uses a rolling k-mer index `u` of the site order. As in the original: `u` is not
/// reset at the window start, so the first `order` columns see a context
/// padded with A (digit 0), and an N contributes digit 0 as well.
pub fn site_fwd(s: &Site, codes: &[u8], p: i64, sh: i64) -> f64 {
    site_fwd_x(s, codes, p, sh, true)
}

/// Sites whose windows reached outside the sequence (names, in order of
/// first occurrence), for one warning each.
static OUTSIDE: std::sync::Mutex<Vec<String>> = std::sync::Mutex::new(Vec::new());

/// Base code at `i`, or N (4) outside the sequence.
///
/// Site windows are placed by the model's MARGIN and WIDTH, and the site
/// track only keeps them inside the sequence for margins up to about 30
/// (records with `31 <= P < L - 30`). With larger or negative margins the
/// original reads the heap around its sequence buffer (or terminates with
/// SIGSEGV when that is far away). The port reads such bases as N, which
/// contributes k-mer digit 0 like any N, and (with `note`) records the site
/// for a warning (see [`warn_outside`]).
///
/// Without `note` (the branch point windows of [`crate::bp`]) a base outside
/// the sequence reads as A (code 0), as from the zero-padded copy those
/// windows were scored on before; their values outside the sequence are
/// never used by the candidates of normal models.
#[inline]
fn code_or_n(codes: &[u8], i: i64, s: &Site, note: bool) -> u8 {
    if i >= 0 && (i as u64) < codes.len() as u64 {
        return codes[i as usize];
    }
    if !note {
        return 0;
    }
    let mut v = OUTSIDE.lock().unwrap_or_else(|e| e.into_inner());
    if !v.contains(&s.name) {
        v.push(s.name.clone());
    }
    4
}

/// Print one warning per site model whose window reached outside the
/// sequence (see [`code_or_n`]).
pub fn warn_outside(model_file: &str) {
    let v = OUTSIDE.lock().unwrap_or_else(|e| e.into_inner());
    for name in v.iter() {
        crate::warn(&format!(
            "model file {model_file}: the window of site {name} (MARGIN/WIDTH) reaches outside the sequence for \
             signals near its ends; bases outside the sequence are scored as N (not handled by the original)"
        ));
    }
}

/// [`site_fwd`]; `note` = record a window outside the sequence for the
/// warning (the branch point windows of [`crate::bp`] reach outside for the
/// first and last bases in normal models, so they are not reported).
pub fn site_fwd_x(s: &Site, codes: &[u8], p: i64, sh: i64, note: bool) -> f64 {
    if s.width < 1 {
        return 0.0;
    }
    let mask = s.n.wrapping_sub(1);
    let st = p - (sh + s.margin as i64);
    let mut u = 0usize;
    let mut acc = 0.0;
    for k in 0..s.width as usize {
        let c = code_or_n(codes, st + k as i64, s, note);
        u = if c == 4 { (u * 4) & mask } else { (u * 4 + c as usize) & mask };
        acc += s.mat[k * s.n + u];
    }
    acc
}

/// 0x41c300: reverse-strand site score with the same (forward) matrix.
///
/// The window `[base, base + WIDTH)` with `base = p - (WIDTH - MARGIN - sh)`
/// is read right to left with complemented bases, so window column k is
/// `3 - codes[base + WIDTH-1-k]`; the consensus lands at column MARGIN as on
/// the forward strand. N contributes digit 0 (not its complement), and `u`
/// is not reset at the window start, as in [`site_fwd`].
pub fn site_rev(s: &Site, codes: &[u8], p: i64, sh: i64) -> f64 {
    site_rev_x(s, codes, p, sh, true)
}

/// [`site_rev`] with the `note` flag of [`site_fwd_x`].
pub fn site_rev_x(s: &Site, codes: &[u8], p: i64, sh: i64, note: bool) -> f64 {
    if s.width < 1 {
        return 0.0;
    }
    let mask = s.n.wrapping_sub(1);
    let w = s.width as i64;
    let base = p - (w - s.margin as i64 - sh);
    let mut u = 0usize;
    let mut acc = 0.0;
    for k in 0..w {
        let c = code_or_n(codes, base + w - 1 - k, s, note);
        u *= 4;
        if c != 4 {
            u += 3 - c as usize;
        }
        u &= mask;
        acc += s.mat[k as usize * s.n + u];
    }
    acc
}

/// 0x41abe0: site track.
///
/// Only records with `31 <= P < L - 30` are scored (all others keep 0.0;
/// in particular every sentinel), so site windows never leave the sequence
/// for site margins up to about 30. Per record:
/// 1. evidence (with -d and a non-zero mark): +100 if the mark has this
///    record's own type bit, else -100 if the record is not a stop and its
///    base is marked IntronMark (0x8000);
/// 2. donors: the GT or GC donor prior (log GT_DONOR_PROB / GC_DONOR_PROB)
///    is added first, chosen by the second consensus base;
/// 3. the site matrix value `v` (per phase `sub` for splice sites, by codon
///    for stops; acceptors get AROUND_BP added in BP models) is added as
///    `v + track` (this operand order is the original's).
fn site_track(map: &[Rec], codes: &[u8], has_marks: bool, m: &Model) -> Vec<f64> {
    let l = codes.len() as i64;
    let mut site = vec![0.0; map.len()];
    let (gtp, gcp) = (m.log_gt_donor, m.log_gc_donor);
    // AROUND_BP is added to acceptor scores of BP models (0x41ac05)
    let around = m.bp.as_ref().map(|b| b.log_around_bp);
    for (r, rec) in map.iter().enumerate() {
        let p = rec.pos as i64;
        let typ = rec.typ;
        let phase = rec.sub as usize;
        if p < 0x1f || p >= l - 0x1e {
            continue;
        }
        if has_marks && rec.mark != 0 {
            let mark = rec.mark;
            if mark & typ as u32 != 0 {
                site[r] += 100.0;
            } else if (typ & 0x22) == 0 && mark & 0x8000 != 0 {
                site[r] -= 100.0;
            }
        }
        let code_at = |i: i64| codes[i as usize];
        let mut v = None;
        // the type tests run in this order; a record has exactly one site bit
        if typ & 4 != 0 {
            // donor+: G at p, T or C at p+1
            if code_at(p + 1) == 3 {
                site[r] += gtp;
                v = Some(site_fwd(&m.donor[phase], codes, p, 0));
            } else if code_at(p + 1) == 1 {
                site[r] += gcp;
                if let Some(g) = &m.donor_gc {
                    v = Some(site_fwd(&g[phase], codes, p, 0));
                }
            }
        } else if typ & 8 != 0 {
            let x = site_fwd(&m.acceptor[phase], codes, p, 1);
            v = Some(around.map_or(x, |a| x + a));
        } else if typ & 1 != 0 {
            v = Some(site_fwd(&m.start, codes, p, 0));
        } else if typ & 2 != 0 {
            // stop+ ends at p: T, a, b = codes[p-2..=p]
            let (a, b) = (code_at(p - 1), code_at(p));
            if a == 0 && b == 0 {
                v = Some(site_fwd(&m.stop[0], codes, p, 2));
            } else if a == 0 && b == 2 {
                v = Some(site_fwd(&m.stop[1], codes, p, 2));
            } else if a == 2 && b == 0 {
                v = Some(site_fwd(&m.stop[2], codes, p, 2));
            }
        } else if typ & 0x40 != 0 {
            // donor-: C at p, A (GT donor) or G (GC donor) at p-1
            if code_at(p - 1) == 0 {
                site[r] += gtp;
                v = Some(site_rev(&m.donor[phase], codes, p, 1));
            } else if code_at(p - 1) == 2 {
                site[r] += gcp;
                if let Some(g) = &m.donor_gc {
                    v = Some(site_rev(&g[phase], codes, p, 1));
                }
            }
        } else if typ & 0x80 != 0 {
            let x = site_rev(&m.acceptor[phase], codes, p, 2);
            v = Some(around.map_or(x, |a| x + a));
        } else if typ & 0x10 != 0 {
            v = Some(site_rev(&m.start, codes, p, 1));
        } else if typ & 0x20 != 0 {
            // stop- starts at p: TTA -> TAA, TCA -> TGA, CTA -> TAG
            let (a, b) = (code_at(p), code_at(p + 1));
            if a == 3 && b == 3 {
                v = Some(site_rev(&m.stop[0], codes, p, 3));
            } else if a == 3 && b == 1 {
                v = Some(site_rev(&m.stop[2], codes, p, 3));
            } else if a == 1 && b == 3 {
                v = Some(site_rev(&m.stop[1], codes, p, 3));
            }
        }
        if let Some(v) = v {
            site[r] = v + site[r];
        }
    }
    site
}

/// Build all four tracks. `marks` are the -d evidence marks, `pen` the `-k`
/// penalty. The builders write disjoint records of each array (a record is
/// never both a begin and an end of the same class), so their order does
/// not matter; it follows main.
///
/// Classes and offsets (record position -> first / one-past-last content
/// base), see re-notes/viterbi.md 1.1:
/// * intergenic (sub 0): begins stop+ / start- (`P+1`), ends ATG+ / stop-
///   (`P`);
/// * intron+ per phase: begins donor+ (`P+2`), ends acceptor+ (`P-1`);
///   intron- per phase: begins acceptor- (`P+2`), ends donor- (`P-1`);
/// * coding+ per frame: begins ATG+ (`P+3`) / acceptor+ (`P+1`), ends stop+
///   (`P-2`) / donor+ (`P`); coding- per frame: begins stop- (`P+3`) /
///   donor- (`P+1`), ends start- (`P-2`) / acceptor- (`P`).
pub fn build(map: &[Rec], codes: &[u8], kmer: &[i32], marks: Option<&[u32]>, m: &Model, pen: f64) -> Tracks {
    let n = map.len();
    let c = Ctx { map, kmer, k: m.order as i64 };
    // 0x41aae0: intergenic, MARKOV column 4
    let mut ig = vec![0.0; n];
    begin_noncoding(&c, &mut ig, &m.nc4.fwd_nf, 0, 0x12, 1);
    end_noncoding(&c, &mut ig, &m.nc4.fwd_nf, &m.nc4.fwd, 0, 0x12, 0x21, 1, 0);
    // 0x41aa30: introns, MARKOV column 3 (forward / reverse complement)
    let mut intron = vec![0.0; n];
    for s in 0..3u8 {
        begin_noncoding(&c, &mut intron, &m.nc3.fwd_nf, s, 4, 2);
        end_noncoding(&c, &mut intron, &m.nc3.fwd_nf, &m.nc3.fwd, s, 4, 8, 2, -1);
        begin_noncoding(&c, &mut intron, &m.nc3.rc_nf, s, 0x80, 2);
        end_noncoding(&c, &mut intron, &m.nc3.rc_nf, &m.nc3.rc, s, 0x80, 0x40, 2, -1);
    }
    // 0x419f60: coding
    let cc = CodingCtx { ring: ring(m), marks, pen };
    let mut cod = vec![0.0; n];
    for f in 0..3u8 {
        begin_coding(&c, &cc, &mut cod, f, 1, 8, 3, 1, true);
        end_coding(&c, &cc, &mut cod, f, 1, 8, 3, 1, 2, 4, -2, 0, true);
        begin_coding(&c, &cc, &mut cod, f, 0x20, 0x40, 3, 1, false);
        end_coding(&c, &cc, &mut cod, f, 0x20, 0x40, 3, 1, 0x10, 0x80, -2, 0, false);
    }
    // 0x41abe0
    let site = site_track(map, codes, marks.is_some(), m);
    Tracks { ig, intron, cod, site }
}
