//! Site map (0x40ed90 and helpers 0x40db60, 0x40e3c0, 0x40ebe0, 0x40e830):
//! the list of candidate signal sites, sorted by position.
//! See re-notes/input.md section 4 and re-notes/viterbi.md section 1.
//!
//! The binary calls this step "encoding", but it does not encode bases: it
//! enumerates every place where a state boundary could be. Each entry
//! ("map record", [`Rec`]) is one candidate signal; the Viterbi
//! ([`crate::viterbi`]) runs over this list, not over bases, with one DP
//! cell per record. The state after a record is fixed by its type (after a
//! donor+ the path is in an intron+, after an ATG in an exon+, ...).
//!
//! Splice sites get three records (`sub` = intron phase 0, 1, 2), because
//! the intron phase is part of the HMM state; codon sites get one. `frame`
//! is the residue mod 3 of the codon starts of the adjacent coding segment,
//! so exon segments only connect records with the same `frame`.
//!
//! Layout of the finished list: 13 head sentinels (position 0), the body
//! records sorted by position, 13 tail sentinels (position L-1). The
//! sentinels let a path start or end in any state at the sequence borders.

use crate::sequence::{A, C, G, T};

/// ATG on the forward strand; pos = the A.
pub const DIR_START: i32 = 0x01;
/// Stop codon TAA/TAG/TGA, forward strand; pos = last base of the codon.
pub const DIR_STOP: i32 = 0x02;
/// Donor GT (or GC with GC_DONOR_ON), forward; pos = the G (first intron base).
pub const DIR_DONOR: i32 = 0x04;
/// Acceptor AG, forward; pos = the G (last intron base).
pub const DIR_ACC: i32 = 0x08;
/// Reverse-strand start (CAT); pos = the T (highest base).
pub const REV_START: i32 = 0x10;
/// Reverse-strand stop (TTA/CTA/TCA); pos = first base.
pub const REV_STOP: i32 = 0x20;
/// Reverse-strand donor (AC, or GC with GC_DONOR_ON); pos = the C (last
/// intron base in forward coordinates).
pub const REV_DONOR: i32 = 0x40;
/// Reverse-strand acceptor (CT); pos = the C (first intron base).
pub const REV_ACC: i32 = 0x80;

/// Number of sentinel records at each end of the map.
pub const SENTINELS: usize = 13;

/// One map record (16 bytes in the binary).
///
/// `typ` is one of the eight site bits above, or a sentinel: `0x100` (head,
/// sequence begin) or `0x200` (tail, sequence end) or'ed with the site bit
/// whose state the sentinel opens or closes. The low byte uses the same bit
/// meanings as the evidence point marks, so `mark & typ` tests "evidence for
/// this very signal".
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct Rec {
    /// site type bitmask (+0)
    pub typ: i32,
    /// 0-based sequence position (+4), see the type constants
    pub pos: i32,
    /// codon-start residue mod 3 of the adjacent coding segment (+8)
    pub frame: u8,
    /// intron phase 0..2 for splice sites (which of the three copies), 0 for
    /// codon sites (+9); called `H` in re-notes/viterbi.md
    pub sub: u8,
    /// evidence marks at `pos` (+12), 0 without -d
    pub mark: u32,
}

impl Rec {
    fn new(typ: i32, pos: i32, frame: u8, sub: u8) -> Rec {
        Rec { typ, pos, frame, sub, mark: 0 }
    }
}

/// Model switches that decide which sites exist.
pub struct Flags {
    /// `$TAA_ON`: TAA (and TTA reverse) are stops
    pub taa: bool,
    /// `$TAG_ON`: TAG (and CTA reverse) are stops
    pub tag: bool,
    /// `$TGA_ON`: TGA (and TCA reverse) are stops
    pub tga: bool,
    /// `$GC_DONOR_ON`: GC donors in addition to GT
    pub gc_donor: bool,
}

/// 0x40db60: body records in emission order (not yet sorted).
///
/// Walks every dinucleotide `(c[p], c[p+1])`; code-4 bases never match.
/// The three copies of a splice site get the frames in the rotation order
/// of the original (`u, u+2, u+1` or `u, u+1, u+2`, mod 3), which is what
/// makes `frame` the codon-start residue of the adjacent exon for each
/// intron phase `sub`.
fn build_body(c: &[u8], f: &Flags) -> Vec<Rec> {
    let l = c.len();
    let mut r = Vec::with_capacity(l + l / 2);
    let push_three = |r: &mut Vec<Rec>, t: i32, pos: usize, frames: [usize; 3]| {
        for (k, &x) in frames.iter().enumerate() {
            r.push(Rec::new(t, pos as i32, x as u8, k as u8));
        }
    };
    for p in 0..l.saturating_sub(1) {
        let (a, b) = (c[p], c[p + 1]);
        // dinucleotide sites; each emits three records (sub = 0, 1, 2)
        if a == G && (b == T || (f.gc_donor && b == C)) {
            let u = p % 3;
            push_three(&mut r, DIR_DONOR, p, [u, (u + 2) % 3, (u + 1) % 3]);
        }
        // GC is a reverse donor too when GC donors are on (GC is its own
        // reverse complement)
        if (a == A && b == C) || (f.gc_donor && a == G && b == C) {
            let u = (p + 2) % 3;
            push_three(&mut r, REV_DONOR, p + 1, [u, (u + 1) % 3, (u + 2) % 3]);
        }
        if a == A && b == G {
            let u = (p + 2) % 3;
            push_three(&mut r, DIR_ACC, p + 1, [u, (u + 2) % 3, (u + 1) % 3]);
        }
        if a == C && b == T {
            let u = p % 3;
            push_three(&mut r, REV_ACC, p, [u, (u + 1) % 3, (u + 2) % 3]);
        }
        // codon sites: codon c[p-1] c[p] c[p+1], one record, frame = residue
        // of the codon's first base
        if p >= 1 {
            let x = c[p - 1];
            let fr = ((p + 2) % 3) as u8;
            let cod = (x, a, b);
            if cod == (A, T, G) {
                r.push(Rec::new(DIR_START, (p - 1) as i32, fr, 0));
            } else if (cod == (T, A, A) && f.taa) || (cod == (T, A, G) && f.tag) || (cod == (T, G, A) && f.tga) {
                r.push(Rec::new(DIR_STOP, (p + 1) as i32, fr, 0));
            } else if cod == (C, A, T) {
                r.push(Rec::new(REV_START, (p + 1) as i32, fr, 0));
            } else if (cod == (T, T, A) && f.taa) || (cod == (C, T, A) && f.tag) || (cod == (T, C, A) && f.tga) {
                r.push(Rec::new(REV_STOP, (p - 1) as i32, fr, 0));
            }
        }
    }
    r
}

/// 0x40e3c0: the 13 head and 13 tail sentinels for a sequence of length `l`.
///
/// Head (pos 0): 0x102 intergenic starts; 0x104 x3 intron+ starts;
/// 0x180 x3 intron- starts; 0x108 x3 exon+ starts; 0x140 x3 exon- starts.
/// Tail (pos L-1): 0x208 x3 intron+ ends; 0x240 x3 intron- ends; 0x204 x3
/// exon+ ends; 0x280 x3 exon- ends; 0x220 intergenic ends. The fixed index
/// ranges of these groups (0, 1-3, 4-6, ...) are relied on by the Viterbi
/// initialisation and the traceback.
fn sentinels(l: usize) -> (Vec<Rec>, Vec<Rec>) {
    let mut head = vec![Rec::new(0x102, 0, 0, 0)];
    for (t, fr) in [(0x104, [0, 0, 0]), (0x180, [0, 0, 0]), (0x108, [0, 2, 1]), (0x140, [0, 1, 2])] {
        for k in 0..3 {
            head.push(Rec::new(t, 0, fr[k], k as u8));
        }
    }
    let u = (l % 3) as u8;
    let last = l as i32 - 1;
    let mut tail = Vec::new();
    for (t, fr) in [
        (0x208, [0, 0, 0]),
        (0x240, [0, 0, 0]),
        (0x204, [u, (u + 2) % 3, (u + 1) % 3]),
        (0x280, [u, (u + 1) % 3, (u + 2) % 3]),
    ] {
        for k in 0..3 {
            tail.push(Rec::new(t, last, fr[k], k as u8));
        }
    }
    tail.push(Rec::new(0x220, last, 0, 0));
    (head, tail)
}

/// 0x40e830: order the body records (indices 13 .. len-13) by position.
///
/// The emission order of [`build_body`] is already sorted by position except
/// for two cases, which this fixes exactly as the original does:
/// * runs of equal `pos` whose length is not 1 or 3 (for example a codon
///   site at the position of a splice triple) are rewritten
///   grouped by type in the order 1, 4, 2, 8, 0x80, 0x40, 0x20, 0x10,
///   keeping the input order within a type; runs of length 1 or 3 are left
///   as they are;
/// * `CATG`: the 0x10 record at p+1 is emitted before the ATG record at p;
///   one forward pass swaps such adjacent pairs.
///
/// The resulting tie order is what the Viterbi's "first maximum wins" rule
/// sees, so it must not be replaced by a generic sort.
fn sort_map(m: &mut [Rec]) {
    const ORDER: [i32; 8] = [1, 4, 2, 8, 0x80, 0x40, 0x20, 0x10];
    let lo = SENTINELS;
    let hi = m.len() - SENTINELS;
    let mut i = lo;
    while i < hi {
        let mut j = i;
        while j < hi && m[j].pos == m[i].pos {
            j += 1;
        }
        let n = j - i;
        if n != 1 && n != 3 {
            let group: Vec<Rec> = m[i..j].to_vec();
            let mut k = i;
            for t in ORDER {
                for r in group.iter().filter(|r| r.typ == t) {
                    m[k] = *r;
                    k += 1;
                }
            }
            // records of other types would be lost by the binary's rewrite;
            // none exist in the body.
        }
        i = j;
    }
    // starts at 14 so that i-1 is the first body record
    for i in 14..hi {
        if m[i].typ == 1 && m[i - 1].typ == 0x10 && m[i].pos < m[i - 1].pos {
            m.swap(i, i - 1);
        }
    }
}

/// 0x40ed90(&map, codes, marks, L, model): build the complete, sorted map.
///
/// `c` are the base codes, `marks` the per-base evidence masks (with -d).
/// Each record gets `mark = marks[pos]`, sentinels included.
pub fn build(c: &[u8], marks: Option<&[u32]>, f: &Flags) -> Vec<Rec> {
    // the original also writes no body for L < 2 (count = 26)
    let body = if c.len() < 2 { Vec::new() } else { build_body(c, f) };
    let (head, tail) = sentinels(c.len());
    let mut m = Vec::with_capacity(body.len() + 2 * SENTINELS);
    m.extend(head);
    m.extend(body);
    m.extend(tail);
    if let Some(mk) = marks {
        for r in m.iter_mut() {
            r.mark = mk[r.pos as usize];
        }
        m.retain(|r| r.typ != 0); // 0x40ebe0, a no-op in practice
    }
    sort_map(&mut m);
    m
}
