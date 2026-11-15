//! Gene prediction for one FASTA record: from the encoded sequence to the
//! list of predicted genes.
//!
//! # Place in the pipeline
//!
//! `mod.rs` drives, per record: command line (`cli.rs`) -> model
//! (`model.rs`) -> FASTA record (`seq.rs`) -> **this module** -> output
//! (`output.rs`). The functions here are called in this order:
//!
//! 1. [`markov_image`] (0x40c020/0x40c040): for every position the index of
//!    the (k+1)-mer starting there, so that every Markov table lookup below
//!    is a single array access.
//! 2. [`build_map`] (0x40ea20 and helpers): the *site map*. gmhmmp is not a
//!    per-nucleotide HMM. It works on a sorted list of candidate gene
//!    boundaries ("sites"): the start and stop codons of both strands plus
//!    14 boundary pseudo sites, linked into ORF chains.
//! 3. [`rbs_prepare`] and [`site_scores`] (0x418b00, 0x418bf0): log score of
//!    the start codon type plus, with `-r`, an RBS (spacer model) or
//!    Kozak-like context term, for every start site.
//! 4. [`run_dp`] (0x4164a0 setup, 0x417660 = 0x415e70 + 0x415d70 + 0x4165c0
//!    + 0x415f40): a dynamic programme over the sites, called "Viterbi" in
//!    the original. Sites are visited once, left to right; each gets the best
//!    log score of a parse of the sequence up to it and a back pointer.
//! 5. [`traceback`] (0x4123b0), [`n_ends`] (0x4120f0), [`logodd`]
//!    (0x429d70/0x429e10): the gene list from the back pointers, partial
//!    marks for genes ending on N, and the per-gene log-odds score.
//!
//! The specification is `re-notes/gmhmmp/viterbi.md` (section numbers in
//! the comments below refer to it). The order of every floating-point
//! operation follows a reference model that reproduced the binary bit
//! for bit. Do not "simplify" sums:
//! `(a + b) + c` and `a + (b + c)` give different bits.
//!
//! # Conventions
//!
//! * Positions are 0-based sequence indices (`i64`). Boundary pseudo sites
//!   lie outside the sequence (-4 and n+3).
//! * Every score is a natural log (`f64`); [`NEG`] = -1e100 stands for
//!   log 0. A score of exactly 0.0 doubles as "not set yet" (see [`better`]).
//! * Codon flags (one byte per site, from `seq::codon_image` and the model's
//!   codon table): `0x20` forward start, `0x10` forward stop, `0x08` reverse
//!   start (CAT, CAC, CAA, CAG), `0x04` reverse stop (TTA, CTA, TCA); bits
//!   0..1 = start type (0 ATG, 1 GTG, 2 TTG, 3 CTG); bits 6..7 = frame tag
//!   `(i % 3) + 1` of the codon ending at index i (0 only for the two
//!   sentinels).
//! * Sites are described by where a gene sits relative to them in forward
//!   coordinates, not by biology. **GB** ("gene begins to the right",
//!   `flags & 0x24`) = forward start or reverse stop; **GE** ("gene ends to
//!   the left", `flags & 0x18`) = forward stop or reverse start.
//! * `pos` of a site is a coding base next to its codon: for GB sites the
//!   first base AFTER the codon (`codon end + 1`), for GE sites the last
//!   base BEFORE it (`codon end - 3`). A gene from GB site `b` to GE site
//!   `e` has the coding region `[b.pos, e.pos]` (codons excluded) and the
//!   length `L = e.pos - b.pos + 7` (codons included). Length tables are
//!   indexed with `L - 1 = e.pos - b.pos + 6`.
//! * Coding and noncoding content are k-th order Markov chains evaluated
//!   RIGHT TO LEFT: `ln P(base i | bases i+1..i+k)`, plus one joint
//!   (k+1)-mer term at the right end of the region (model.md 6.5).
//!
//! # Speed tricks (all exact)
//!
//! The original scans the site array linearly in several places. The port
//! replaces these scans by precomputed index lists (`Dp::ge_sites`,
//! `Dp::next_fstop`, `Dp::rstarts`) that visit exactly the sites the scan
//! acts on, in the same order, so every floating-point operation is the
//! same. The Markov image turns each table lookup into one load. Coding
//! sums of the starts of one ORF are extended incrementally (as in the
//! original), never recomputed.

use super::evidence::Maps;
use super::model::{Markov, Model, NEG};

/// ln 4 as the double literal of the binary (0x40e0b0): the price of one
/// base that no model explains (probability 1/4).
const LN4: f64 = 1.3862943611198906;
/// ln 16 = 2 * ln 4 (two unexplained bases), double literal of the binary.
const LN16: f64 = 2.772588722239781;

/// Markov image (0x40c040): for every position i the index of the
/// (k+1)-mer `codes[i..=i+k]` with base i as the most significant base-4
/// digit, so `table[idx[i]]` is `ln P(base i | bases i+1..i+k)` for the
/// conditional tables and `ln P(bases i..i+k)` for the joint ones.
///
/// * `codes`: base codes 0..4 (A, C, G, T, other) as the image sees them:
///   after gap splitting but BEFORE R/Y stop resolution (main calls
///   0x40c020 before 0x407260), so R/Y count as 4 here.
/// * Windows containing a code 4 get the fill index `4^(k+1)`, the extra
///   last slot of every Markov table (ln 0.48 coding, ln 0.52 noncoding).
/// * The last k entries (no complete window) are 0, as in the binary.
///
/// Literal port of the binary's right-to-left loop with its 16-bit index
/// field: `idx = (u16)(code << 2k) + u9` (or `(u16)(4 << 2k)` with an N
/// in the window), then `u9 = ((code << 2k) + u9) & 0xffff) >> 2`. For
/// k <= 7 this is the exact (k+1)-mer index. For k = 8 (the largest order
/// main accepts) the first base's digit (bits 16..17) is cut off before
/// the rolling sum is shifted, so no new base ever enters the sum: the
/// index of position i is the last-8-bases index of the sequence shifted
/// right by 2 per position (0 from 17 positions before the end on), and N
/// windows get `(u16)4^9 = 0`. Reproduced as is (deterministic, no read
/// outside the buffers) [V order-8 model].
pub fn markov_image(codes: &[u8], order: i32) -> Image {
    let n = codes.len();
    let k = order as u32;
    let mut idx = vec![0usize; n];
    let mut u8n: u32 = 0; // N tracker (3 << 2t per N, shifted out after k+1 steps)
    let mut u9: u32 = 0; // rolling sum
    // the last k positions: no complete window, index 0 (n > k always:
    // records have at least 20 letters and k <= 8)
    let mut p = n as i64 - 1;
    for t in 0..k {
        let b = codes[p as usize] as u32;
        if b == 4 {
            u8n += 3 << (2 * t);
        } else {
            u9 += b << (2 * t);
        }
        p -= 1;
    }
    let fill16 = ((4u32 << (2 * k)) & 0xffff) as usize;
    while p >= 0 {
        let b = codes[p as usize] as u32;
        let add = if b == 4 {
            u8n += 3 << (2 * k);
            0
        } else {
            b << (2 * k)
        };
        idx[p as usize] = if u8n == 0 { ((add + u9) & 0xffff) as usize } else { fill16 };
        u8n >>= 2;
        u9 = ((add + u9) & 0xffff) >> 2;
        p -= 1;
    }
    // fill16 is the fill slot 4^(k+1) for k <= 7 (0 for k = 8, see above);
    // reads outside the image use the real fill slot
    Image { v: idx, fill: 4usize << (2 * k), oob: std::cell::Cell::new(false) }
}

/// The Markov image of one record with checked reads.
///
/// The DP and the logodd read `idx[x]` for windows next to gene ends. With
/// MINC < ORDM + 9 (only small hand-made MINC values) a gene can be so
/// short that a window lies partly outside the sequence: the binary then
/// reads the 4-byte image cells in front of (or behind) its array, i.e.
/// heap metadata or neighbouring data (undefined). Intended difference:
/// the port reads such windows as the fill slot `4^(k+1)` (a window with
/// an unknown base, like an N) and sets `oob`, so that main can warn once
/// per record.
pub struct Image {
    pub v: Vec<usize>,
    /// table slot used for reads outside the image
    fill: usize,
    /// set when a read fell outside the image
    pub oob: std::cell::Cell<bool>,
}

impl Image {
    /// `idx[i]`, or the fill slot (and `oob` set) outside 0..n.
    #[inline]
    pub fn at(&self, i: i64) -> usize {
        match self.v.get(i as usize) {
            Some(&x) => x,
            None => {
                self.oob.set(true);
                self.fill
            }
        }
    }
}

/// One site of the site map (0x60-byte record in the binary, viterbi.md
/// 4.2): a candidate gene boundary, and at the same time a DP cell.
#[derive(Debug, Clone)]
pub struct Site {
    /// Codon flags (module docs): direction, start type and frame tag.
    /// 0 marks a deleted site while the map is built.
    pub flags: u8,
    /// 0-based position of the coding base next to the codon (module docs);
    /// -4 and n+3 for the two sentinels.
    pub pos: i64,
    /// ORF chain link, always to a site further LEFT. Forward ORF: stop ->
    /// nearest start -> next start further left -> ... -> None. Reverse
    /// ORF: rightmost start -> ... -> leftmost start -> stop (None).
    pub link: Option<usize>,
    /// DP back pointer: the predecessor of this site on the best parse.
    pub back: Option<usize>,
    /// DP value: best log score of a parse ending at this site. Real sites
    /// start at NEG, boundary sites at their fixed initial scores.
    pub score: f64,
    /// Noncoding sum from this GE site to the next GE site
    /// (`sum NC[idx[E]]` for E in [pos+4, next.pos+3], see `Dp::nc_sums`).
    pub ncs: f64,
    /// Start term: ln P(start codon type) + RBS/context log ratio
    /// ([`site_scores`]). Set for every site, used for starts only.
    pub st: f64,
    /// Class of the best gene ending at this GE site: b'1' (COD1 model won
    /// or tie) or b'2' (COD2 strictly better). Written into GB sites too,
    /// where it is never read.
    pub cls: u8,
    /// RBS spacer chosen in spacer mode (0..RBSL-1), -1 = none.
    pub spacer: i8,
    /// Evidence value; after compaction `ev + E` (`E` = `-E`, default 100),
    /// so exactly E for every site without `-e` evidence.
    pub ev: f64,
    /// Evidence flag (only with `-e`).
    pub evf: u8,
    /// Gene score: the winning class's coding log likelihood (incl. the
    /// CD1P/CD2P prior, without length, start and path terms) of the gene
    /// whose start is this site. -1.0 for boundary sites, 0.0 otherwise
    /// until set.
    pub gs: f64,
    /// `link` is a stale pointer of the binary (see [`compact`]): it points
    /// to the freed temporary copy of the linked site, whose link leads to
    /// the final copy (or it is a NULL that replaced a real chain link).
    /// Only possible with MINC < 7.
    pub link_stale: bool,
}

impl Site {
    fn new(flags: u8, pos: i64) -> Site {
        Site { flags, pos, link: None, back: None, score: 0.0, ncs: 0.0, st: 0.0, cls: 0, spacer: -1, ev: 0.0, evf: 0, gs: 0.0, link_stale: false }
    }
}

/// Link of a temporary site. During chain building it names another
/// temporary site; after compaction the binary reuses the link field of a
/// copied site as a forward pointer to its copy in the final array.
#[derive(Clone, Copy)]
enum TLink {
    /// index into the temporary array
    Temp(usize),
    /// index into the compacted (final) array
    Final(usize),
}

/// Site of the temporary array (0x40ddb0/0x40e0b0) before compaction.
struct TempSite {
    s: Site,
    tlink: Option<TLink>,
}

/// Frame tag bits (`((j % 3) + 1) << 6`) of a codon ending at index j.
fn frame_tag_of_codon_end(j: i64) -> u8 {
    (((j % 3) + 1) << 6) as u8
}

/// Temporary site array (0x40ddb0 real sites + 0x40e0b0 boundaries,
/// viterbi.md 4.3): 7 left boundary sites, one site per nonzero codon flag
/// in sequence order, 7 right boundary sites.
///
/// * `codon`: codon flag image (one byte per codon end, `seq::codon_image`).
/// * `p`: the `-i` probability ("initiation/termination in noncoding
///   state", possibly replaced per record by `-c`). `lp = ln p` starts a
///   parse in noncoding state, `l1p = ln(1-p)` starts it inside a gene.
///
/// Left boundary (all with gs = -1.0):
/// * 0: GE sentinel (flags 0x18) at pos -4, score lp: "a gene ended just
///   before the sequence", so noncoding may start at base 0.
/// * 1..6: pseudo forward starts and reverse stops at pos 0, 1, 2 (one per
///   frame) for genes running off the left end (`<`). Starting at pos 1 or
///   2 leaves 1 or 2 bases unexplained, charged ln 4 each.
///
/// Right boundary: for each of pos n-3, n-2, n-1 a pseudo forward stop and
/// a pseudo reverse start (genes running off the right end, `>`), frame tag
/// of a codon ending at n, n+1, n+2; and the GB end sentinel (0x24) at
/// pos n+3. All with score 0.0.
fn build_temp_sites(codon: &[u8], n: i64, p: f64, ev: Option<&Maps>) -> Vec<TempSite> {
    // Both tests are `ucomisd` + `ja/jp` (0x40e0c4) and `ucomisd` + `jb`
    // (0x40e0f0): an unordered compare takes the log branch, so a NaN p
    // (`-i nan`) gives ln(NaN) = NaN for lp AND l1p (Ghidra's pseudo-C
    // shows the second test as a plain `0 < 1-p`). Same in Dp::terminate.
    let lp = if p > 0.0 || p.is_nan() { p.ln() } else { NEG };
    let l1p = if (1.0 - p) > 0.0 || (1.0 - p).is_nan() { (1.0 - p).ln() } else { NEG };
    let boundary_site = |flags: u8, pos: i64, score: f64| {
        let mut s = Site::new(flags, pos);
        s.score = score;
        s.gs = -1.0;
        TempSite { s, tlink: None }
    };
    let mut t = vec![
        boundary_site(0x18, -4, lp),
        boundary_site(0xe0, 0, l1p),
        boundary_site(0xc4, 0, l1p),
        boundary_site(0x60, 1, l1p - LN4),
        boundary_site(0x44, 1, l1p - LN4),
        boundary_site(0xa0, 2, l1p - LN16),
        boundary_site(0x84, 2, l1p - LN16),
    ];
    for (i, &c) in codon.iter().enumerate() {
        if c != 0 {
            // GB (fwd start, rev stop): first base after the codon;
            // GE (fwd stop, rev start): last base before the codon
            let pos = if c & 0x18 == 0 { i as i64 + 1 } else { i as i64 - 3 };
            let mut s = Site::new(c, pos);
            s.score = NEG;
            s.gs = 0.0;
            if let Some(m) = ev {
                // 0x40ddb0: exact key lookups (key = codon end index) in
                // tree 0x58 (fwd stops) / 0x60 (rev stops)
                let key = i as i64;
                if c & 0x10 != 0 && m.rec_plus.contains_key(&key) {
                    s.evf = 1;
                    s.ev = m.rec_plus[&key].score;
                } else if c & 0x04 != 0 && m.rec_minus.contains_key(&key) {
                    s.evf = 1;
                    s.ev = m.rec_minus[&key].score;
                }
            }
            t.push(TempSite { s, tlink: None });
        }
    }
    if let Some(m) = ev {
        // left pseudo reverse stops (temp 2, 4, 6 = frame tags 3, 1, 2) look
        // up the first codon end of their frame (2, 0, 1); only the flag is set
        for (i, key) in [(2usize, 2i64), (4, 0), (6, 1)] {
            if m.rec_minus.contains_key(&key) {
                t[i].s.evf = 1;
            }
        }
    }
    for pp in [n - 3, n - 2, n - 1] {
        let tg = frame_tag_of_codon_end(pp + 3);
        let mut a = boundary_site(tg | 0x10, pp, 0.0);
        if let Some(m) = ev {
            if m.rec_plus.contains_key(&pp) {
                a.s.evf = 1;
            }
        }
        t.push(a);
        t.push(boundary_site(tg | 0x08, pp, 0.0));
    }
    t.push(boundary_site(0x24, n + 3, 0.0));
    t
}

/// Reverse ORF chains of one frame (0x40d290, viterbi.md 4.4), ascending.
///
/// On the reverse strand the stop (0x04) is LEFT of its starts (0x08).
/// Every start is linked to the previous chain element on its left (start
/// or the stop); the stop's link is None. Deletions (flags := 0):
/// * reverse starts of this frame before the first reverse stop (cannot
///   happen in practice: the left boundary has a pseudo stop per frame);
/// * a stop followed by no start before the next stop;
/// * starts too close to the stop: while no start has been accepted yet,
///   a start giving a gene shorter than `minlen` (MINC) is dropped; the
///   first one that passes ends the test, further starts are longer.
fn rev_chain(t: &mut [TempSite], frame: u8, minlen: i64) {
    let cnt = t.len();
    let mut k = 0;
    while !((t[k].s.flags & 4 != 0) && (t[k].s.flags & 0xc0) == frame) {
        if (t[k].s.flags & 8 != 0) && (t[k].s.flags & 0xc0) == frame {
            t[k].s.flags = 0;
        }
        k += 1;
    }
    t[k].tlink = None;
    // prev = last element of the current chain (stop or start)
    let mut prev = k;
    for j in k + 1..cnt {
        let b = t[j].s.flags;
        if b & 0xc0 != frame {
            continue;
        }
        if b & 8 == 0 {
            if b & 4 != 0 {
                // next stop: the previous stop without any start is dropped
                if t[prev].s.flags & 8 == 0 {
                    t[prev].s.flags = 0;
                }
                t[j].tlink = None;
                prev = j;
            }
        } else if t[prev].s.flags & 8 == 0 && (t[j].s.pos - t[prev].s.pos) + 7 < minlen {
            // first start after the stop and gene length < MINC
            t[j].s.flags = 0;
        } else {
            t[j].tlink = Some(TLink::Temp(prev));
            prev = j;
        }
    }
    if (t[prev].s.flags & 4 != 0) && (t[prev].s.flags & 0xc0) == frame {
        // trailing stop without start
        t[prev].s.flags = 0;
    }
}

/// Forward ORF chains of one frame (0x40d3b0, viterbi.md 4.4): mirror of
/// [`rev_chain`], descending from the right end. The forward stop (0x10)
/// is RIGHT of its starts (0x20): stop -> nearest start -> next start
/// further left -> ... -> None. The same deletions apply.
///
/// The control flow follows the binary's do-while loop closely (including
/// the case where the first stop is at index 0, which deletes it).
fn fwd_chain(t: &mut [TempSite], frame: u8, minlen: i64) {
    let cnt = t.len();
    let mut k = cnt as i64 - 1;
    while !((t[k as usize].s.flags & 0x10 != 0) && (t[k as usize].s.flags & 0xc0) == frame) {
        if (t[k as usize].s.flags & 0x20 != 0) && (t[k as usize].s.flags & 0xc0) == frame {
            t[k as usize].s.flags = 0;
        }
        k -= 1;
    }
    // cur = last element of the current chain; the first stop's link is not
    // cleared here (it is still None from construction)
    let mut cur = k as usize;
    let mut j = k - 1;
    if j >= 0 {
        loop {
            let ju = j as usize;
            let b = t[ju].s.flags;
            if b & 0xc0 == frame {
                if b & 0x20 == 0 {
                    if b & 0x10 != 0 {
                        // next stop: close the previous chain
                        if t[cur].s.flags & 0x20 == 0 {
                            t[cur].s.flags = 0; // stop without start
                        } else {
                            t[cur].tlink = None; // leftmost start of the chain
                        }
                        cur = ju;
                    }
                } else if t[cur].s.flags & 0x20 == 0 && (t[cur].s.pos - t[ju].s.pos) + 7 < minlen {
                    // first start before the stop and gene length < MINC
                    t[ju].s.flags = 0;
                } else {
                    t[cur].tlink = Some(TLink::Temp(ju));
                    cur = ju;
                }
            }
            if j - 1 <= -1 {
                break;
            }
            j -= 1;
        }
        if t[cur].s.flags & 0x10 == 0 || (t[cur].s.flags & 0xc0) != frame {
            t[cur].tlink = None;
            return;
        }
    }
    // trailing stop without start (or the loop never ran)
    t[cur].s.flags = 0;
}

/// Compaction (0x40e030 + 0x40d4f0, viterbi.md 4.5): copy the surviving
/// temporary sites into the final array, approximately sorted by `pos`.
///
/// Real sites come in codon-end order, but GB and GE sites use different
/// pos conventions (end+1 vs end-3), so the order by pos differs locally.
/// The binary fixes this with a window sort: among the current site and
/// the next 3 temporary entries (deleted ones count towards the 3 but are
/// never picked) it copies the one with the strictly smallest pos (first
/// minimum wins), then repeats. The result is not guaranteed to be sorted;
/// the DP relies on exactly this order.
///
/// Links are translated through the copy's forward pointer (0x40d4f0: the
/// copy gets `target->link`, and the temporary site's link then points to
/// its copy). With MINC >= 7 a link target is always copied before the
/// site linking to it [V]. With MINC < 7 a gene can have its stop before
/// its start in pos order, the window sort then copies the linking site
/// first, and the binary copies the target's own temporary link: a pointer
/// to the next chain element's TEMPORARY site (freed after compaction; its
/// link field then points to that element's copy). The port resolves such
/// a link to the final copy and marks it `link_stale` ([`gc_error`] needs
/// the extra node). `ev` becomes `ev + e` (`e` = `-E`), or exactly `e` for
/// a flagged site with ev == 0.0.
fn compact(t: &mut [TempSite], e: f64) -> Vec<Site> {
    let cnt = t.len();
    let mut out: Vec<Site> = Vec::new();
    // (final index, temporary index of the stale target)
    let mut stale: Vec<(usize, usize)> = Vec::new();
    let mut k = 0;
    while k < cnt && t[k].s.flags == 0 {
        k += 1;
    }
    loop {
        while k < cnt && t[k].s.flags == 0 {
            k += 1;
        }
        if k >= cnt {
            break;
        }
        let mut best = k;
        let mut bp = t[k].s.pos;
        for j in k + 1..k + 4 {
            if j >= cnt {
                break;
            }
            if t[j].s.flags != 0 && t[j].s.pos < bp {
                best = j;
                bp = t[j].s.pos;
            }
        }
        let s = &t[best].s;
        let mut c = s.clone();
        c.back = None;
        c.ev = if s.evf == 0 || s.ev != 0.0 { s.ev + e } else { e };
        c.link = match t[best].tlink {
            None => None,
            Some(TLink::Temp(x)) => match t[x].tlink {
                Some(TLink::Final(ci)) => Some(ci),
                Some(TLink::Temp(y)) => {
                    // target not copied yet: stale pointer to temp site y
                    stale.push((out.len(), y));
                    None
                }
                None => {
                    // target not copied yet and the end of its chain: the
                    // binary copies NULL, and its DP later dereferences it
                    // (SIGSEGV [V]); marked stale so that main stops
                    c.link_stale = true;
                    None
                }
            },
            // unreachable: a copied site has flags 0 and is never picked
            Some(TLink::Final(ci)) => Some(ci),
        };
        t[best].tlink = Some(TLink::Final(out.len()));
        t[best].s.flags = 0;
        out.push(c);
    }
    for (ci, y) in stale {
        // every surviving site is copied by now; a deleted target (never
        // linked in practice) stays None
        if let Some(TLink::Final(fi)) = t[y].tlink {
            out[ci].link = Some(fi);
        }
        out[ci].link_stale = true;
    }
    out
}

/// GC-of-ORFs pass (0x40dcc0 with 0x4286a0, viterbi.md 9), run by main
/// after the length tables and before the DP. Its GC values are never
/// read, so only its error exit is observable: `GCofInterval(a, b)` with
/// `b < a` writes `error in code GCofInterval` to stderr and exits 1.
/// Calls: for every forward stop with a link, `GC(x.pos, stop.pos)` for
/// each chain element x except the last (leftmost) one; for every reverse
/// start with a link, `GC(root.pos, start.pos)` (root = chain end). An
/// error needs a gene shorter than 7 bp, i.e. MINC < 7 (MINC == ORDM == 4
/// [V]). A stale link (see [`compact`]) first visits the temporary copy,
/// whose link is never NULL, so that element is checked even if it is the
/// last one. Returns true if the binary would exit there.
pub fn gc_error(s: &[Site]) -> bool {
    for st in s.iter() {
        if st.flags & 0x10 == 0 {
            continue;
        }
        let mut cur = st.link;
        let mut via_stale = st.link_stale;
        while let Some(x) = cur {
            if via_stale {
                if st.pos < s[x].pos {
                    return true;
                }
                via_stale = false;
                continue;
            }
            if s[x].link.is_none() {
                break;
            }
            if st.pos < s[x].pos {
                return true;
            }
            via_stale = s[x].link_stale;
            cur = s[x].link;
        }
    }
    for st in s.iter() {
        if st.flags & 8 == 0 || st.link.is_none() {
            continue;
        }
        let mut r = st.link.unwrap();
        while let Some(l) = s[r].link {
            r = l;
        }
        if st.pos < s[r].pos {
            return true;
        }
    }
    false
}

/// Size requests for the length tables (0x40d6a0, 0x40d700).
///
/// * `max_gene_len` (map+0x50, "m50"): 7 + the longest span `pos -
///   pos(chain root)` over all linked sites = the longest possible gene
///   length L. 7 for an empty map.
/// * `max_nc_span` (map+0x54, "m54"): the largest `pos(GB) - pos(last GE
///   before it)` = longest noncoding length + 7; if smaller than
///   `max_gene_len` it becomes `max_gene_len + 1000`.
///
/// The binary updates "last GE" to j+1 before testing j+1 and reads one
/// record past the array at the last step; the `j + 1 < len` guard skips
/// that read, whose value only affected a local.
fn max_lengths(c: &[Site]) -> (i64, i64) {
    let mut max_gene_len = 7;
    if !c.is_empty() {
        let mut mx = 0;
        for s in c {
            if s.link.is_some() {
                let mut r = s;
                while let Some(l) = r.link {
                    r = &c[l];
                }
                mx = mx.max(s.pos - r.pos);
            }
        }
        max_gene_len = mx + 7;
    }
    let mut max_nc_span = 0;
    let mut last = 0;
    for j in 1..c.len() {
        if c[j].flags & 0x24 != 0 {
            max_nc_span = max_nc_span.max(c[j].pos - c[last].pos);
        }
        if j + 1 < c.len() && c[j + 1].flags & 0x18 != 0 {
            last = j + 1;
        }
    }
    if max_nc_span < max_gene_len {
        max_nc_span = max_gene_len + 1000;
    }
    (max_gene_len, max_nc_span)
}

/// The site map of one record (in main).
pub struct SiteMap {
    /// Final site array S, approximately sorted by pos; S[0] is the left GE
    /// sentinel, the last site the GB end sentinel.
    pub sites: Vec<Site>,
    /// Longest possible gene length (map+0x50); request for
    /// `Model::coding_lengths`.
    pub max_gene_len: i64,
    /// Longest noncoding span + 7, or `max_gene_len + 1000` (map+0x54); main
    /// passes `max_nc_span + 8` to `Model::noncoding_lengths`.
    pub max_nc_span: i64,
}

/// Build the site map of one record (0x40d1f0, 0x40d7a0, 0x40ea20).
///
/// * `codon`: codon flag image (after R/Y resolution and `-e` edits).
/// * `p`: `-i` probability (per record with `-c`), see [`build_temp_sites`].
/// * `minlen`: MINC, minimum gene length incl. codons.
/// * `strand`: `-s` value; b'd' forward only, b'r' reverse only, anything
///   else both. Sites of the unwanted strand are deleted outright (the two
///   sentinels at index 0 and len-1 are kept).
/// * `e`: `-E` evidence score (added to every site's `ev`).
pub fn build_map(codon: &[u8], p: f64, minlen: i64, strand: u8, e: f64, ev: Option<&Maps>) -> SiteMap {
    let n = codon.len() as i64;
    let mut t = build_temp_sites(codon, n, p, ev);
    let len = t.len();
    if strand == b'd' {
        for frame in [0x40, 0x80, 0xc0] {
            fwd_chain(&mut t, frame, minlen);
        }
        for k in 1..len - 1 {
            if t[k].s.flags & 0x0c != 0 {
                t[k].s.flags = 0;
            }
        }
    } else {
        for frame in [0x40, 0x80, 0xc0] {
            rev_chain(&mut t, frame, minlen);
        }
        if strand == b'r' {
            for k in 1..len - 1 {
                if t[k].s.flags & 0x30 != 0 {
                    t[k].s.flags = 0;
                }
            }
        } else {
            for frame in [0x40, 0x80, 0xc0] {
                fwd_chain(&mut t, frame, minlen);
            }
        }
    }
    let sites = compact(&mut t, e);
    let (max_gene_len, max_nc_span) = max_lengths(&sites);
    SiteMap { sites, max_gene_len, max_nc_span }
}

// ---------------------------------------------------------------- RBS

/// Precomputed RBS spacer-mode tables of one record (0x418b00, model.md 8),
/// W = RBSS (motif width), RBSL = number of spacer lengths. All values are
/// LINEAR probabilities (products), not logs.
pub struct RbsPrep {
    /// `prod_t RBSM[t][digit_t(u)]` for every W-mer u (R[0]): motif
    /// probability. digit t counts from the least significant digit, so
    /// matrix row t scores motif base W-1-t (rows are stored right to left).
    motif_fwd: Vec<f64>,
    /// `motif_fwd` of the reverse complement W-mer (R[1]).
    motif_rev: Vec<f64>,
    /// `prod_t nonc0[digit_t(u)]` (R[2]): background from the NONC order-0
    /// sums (linear, not normalised, as stored in the Markov object).
    bg_fwd: Vec<f64>,
    /// `bg_fwd` of the reverse complement W-mer (R[3]).
    bg_rev: Vec<f64>,
    /// `exp(-(W + j) / NDEC)` for spacer j (R[4]): background of W + j
    /// noncoding bases under the geometric noncoding length model.
    spacer_bg: Vec<f64>,
    /// W-mer index starting at each position (R[10]), first base most
    /// significant; N counts as C; the last W-1 entries hold the shorter
    /// suffix (right-aligned).
    wmer: Vec<usize>,
}

/// 0x418b00, only called for RBSS > 0 and RBSL >= 0 (spacer mode).
/// `nuc` = codes after gap split, before R/Y resolution (image byte 0).
///
/// Side effect kept from the binary: the model's spacer distribution RBSD
/// is renormalised IN PLACE on every record (sum in index order, then
/// divide), so its bits can drift from record to record [V gdb].
pub fn rbs_prepare(m: &mut Model, nuc: &[u8]) -> RbsPrep {
    let w = m.rbss as usize;
    let l = m.rbsl as usize;
    let n_wmers = 1usize << (2 * w);
    let bg = m.markov.nonc0;
    let mx = &m.rbsm;
    let mut motif_fwd = vec![1.0f64; n_wmers];
    let mut bg_fwd = vec![1.0f64; n_wmers];
    let n = nuc.len();
    let mut wmer = vec![0usize; n];
    // last W positions: right-aligned partial W-mers (0x418700)
    let mut v = 0usize;
    for t in 0..w.min(n) {
        let b = nuc[n - 1 - t];
        let b = if b == 4 { 1 } else { b as usize };
        v += b << (2 * t);
        wmer[n - 1 - t] = v;
    }
    if n > w {
        // roll leftwards: new base enters as the most significant digit,
        // the base W to the right drops out
        for i in (0..n - w).rev() {
            let b = nuc[i];
            let b = if b == 4 { 1 } else { b as usize };
            v = (b << (2 * w - 2)) + (v >> 2);
            wmer[i] = v;
        }
    }
    // products in the order t = 0..W-1 (0x4187a0)
    for t in 0..w {
        for u in 0..n_wmers {
            let dg = (u >> (2 * t)) & 3;
            bg_fwd[u] *= bg[dg];
            motif_fwd[u] *= mx[dg + t * 4];
        }
    }
    let mut motif_rev = vec![0.0; n_wmers];
    let mut bg_rev = vec![0.0; n_wmers];
    for u in 0..n_wmers {
        // r = W-digit reverse complement of u
        let mut r = 0usize;
        let mut x = u;
        for _ in 0..w {
            r = (r * 4 + 3) - (x & 3);
            x >>= 2;
        }
        motif_rev[u] = motif_fwd[r];
        bg_rev[u] = bg_fwd[r];
    }
    let spacer_bg: Vec<f64> = (0..l).map(|j| (-((w + j) as f64) / m.ndec).exp()).collect();
    // 0x4188d0: in-place renormalisation of RBSD
    let mut tot = 0.0;
    for j in 0..l {
        tot += m.rbsd[j];
    }
    for j in 0..l {
        m.rbsd[j] /= tot;
    }
    RbsPrep { motif_fwd, motif_rev, bg_fwd, bg_rev, spacer_bg, wmer }
}

/// Flat views of the binary's memory for the context-mode RBS score
/// (RBSL < 0), so that out-of-row reads caused by an N (code 4, or -1
/// after complementing) in the window give the same values as in the
/// binary (viterbi.md sec. 10):
///
/// * the motif matrix `mx` is indexed `row * 4 + base`; base 4 reads the
///   next row's A, base -1 the previous row's T;
/// * `bgph` = the 17 doubles of the Markov object from K+0xd0: NONC
///   order-0 sums [0..4), COD1 phase 0..2 compositions [4..16), and at
///   [16] the two i32 at K+0x150 (order, N) read as one double. So
///   `bg[4]` = `ph[0][0]`, `ph[j][4]` = `ph[j+1][0]`, `ph[0][-1]` = `bg[3]`.
struct CtxMemory<'a> {
    mx: &'a [f64],
    bgph: [f64; 17],
}

impl CtxMemory<'_> {
    /// `Mx[i]` (flat index). `Mx[-1]` is the glibc chunk-size field in
    /// front of the matrix (`new double[4W]`, no array cookie):
    /// `((4W*8 + 8 + 15) & ~15) | PREV_INUSE`, minimum chunk 32, read as a
    /// double (a subnormal). This undefined read is reproduced exactly.
    /// Reads past the end hit uninitialised padding in the original; 0.0
    /// is used there (documented divergence, README).
    #[inline]
    fn motif_at(&self, i: i64) -> f64 {
        if i == -1 {
            let req = self.mx.len() as u64 * 8;
            let chunk = ((req + 8 + 15) & !15).max(32) | 1;
            f64::from_bits(chunk)
        } else if i < 0 {
            0.0
        } else {
            self.mx.get(i as usize).copied().unwrap_or(0.0)
        }
    }
    /// Background `bg[b]` (NONC order-0). `bg[-1]` is K+0xc8, an unused
    /// table pointer slot that is NULL, i.e. 0.0 (the ratio then becomes
    /// inf or NaN, as in the binary).
    #[inline]
    fn bg_at(&self, b: i64) -> f64 {
        if b < 0 {
            0.0
        } else {
            self.bgph[(b as usize).min(16)]
        }
    }
    /// COD1 phase-j composition `ph[j][b]`, with the neighbouring doubles
    /// for b = -1 and b = 4.
    #[inline]
    fn phase_at(&self, j: usize, b: i64) -> f64 {
        if b < 0 {
            self.bgph[4 + 4 * j - 1]
        } else {
            self.bgph[(4 + 4 * j + b as usize).min(16)]
        }
    }
}

/// Context score from the four products: `ln((pm/pb) * (pm2/pb2))`, NEG
/// for an exact 0.
fn ctx_log_ratio(pm: f64, pb: f64, pm2: f64, pb2: f64) -> f64 {
    let r = (pm / pb) * (pm2 / pb2);
    if r == 0.0 {
        NEG
    } else {
        r.ln()
    }
}

/// Context-mode score of a forward start (0x418bf0, RBSL = -g < 0).
///
/// The W-row matrix is a Kozak-like profile around the start codon, which
/// itself is not scored: rows W-1..g score the W-g bases just before the
/// codon against the NONC background, rows g-1..0 the g bases right after
/// it against the COD1 composition of codon positions 0, 1, 2, 0, ...
/// `pos` = first base after the codon, `i3` = first base of the codon.
/// The caller's range test keeps the window inside the sequence if
/// W >= 2g + 3 (forward) and W >= 2g + 1 (reverse); the shipped model has
/// W = 9, g = 3.
fn ctx_fwd(cx: &CtxMemory, nuc: &[u8], pos: i64, w: i64, l: i64) -> f64 {
    let i3 = pos - 3;
    let upstream_first = (-l - w) + i3; // i3 - (W - g)
    let (mut pm, mut pb) = (1.0f64, 1.0f64);
    let mut row = w - 1;
    for e in upstream_first..i3 {
        let b = nuc[e as usize] as i64;
        pm *= cx.motif_at(b + row * 4);
        pb *= cx.bg_at(b);
        row -= 1;
    }
    // rows continue downwards after the upstream part (g-1 for W-g bases)
    let mut row = if upstream_first < i3 { w - 1 - (i3 - upstream_first) } else { w - 1 };
    let (mut pm2, mut pb2) = (1.0f64, 1.0f64);
    let mut phi = 0usize;
    for e in pos..pos - l {
        let b = nuc[e as usize] as i64;
        pm2 *= cx.motif_at(b + row * 4);
        pb2 *= cx.phase_at(phi, b);
        phi = if phi != 2 { phi + 1 } else { 0 };
        row -= 1;
    }
    ctx_log_ratio(pm, pb, pm2, pb2)
}

/// Context-mode score of a reverse start: mirror of [`ctx_fwd`] on the
/// complemented bases (`3 - code`; an N gives -1). `pos` = last coding
/// base left of the CAT codon, `i4` = last base of CAT. The upstream part
/// runs from `i4 + W - g` down to `i4 + 1`, then the g coding bases from
/// `pos` leftwards. Note `pb` is multiplied before `pm` here (independent
/// products, order irrelevant).
fn ctx_rev(cx: &CtxMemory, nuc: &[u8], pos: i64, w: i64, l: i64) -> f64 {
    let i3 = pos;
    let i4 = i3 + 3;
    let upstream_last = w + l + i4;
    let (mut pm, mut pb) = (1.0f64, 1.0f64);
    let mut row = w - 1;
    let mut e = upstream_last;
    while e > i4 {
        let b = 3 - nuc[e as usize] as i64;
        pb *= cx.bg_at(b);
        pm *= cx.motif_at(b + row * 4);
        row -= 1;
        e -= 1;
    }
    let mut row = if i4 < upstream_last { w - 1 - (upstream_last - i4) } else { w - 1 };
    let (mut pm2, mut pb2) = (1.0f64, 1.0f64);
    let mut phi = 0usize;
    let mut e = i3;
    while e != i3 + l {
        let b = 3 - nuc[e as usize] as i64;
        pm2 *= cx.motif_at(b + row * 4);
        pb2 *= cx.phase_at(phi, b);
        phi = if phi != 2 { phi + 1 } else { 0 };
        row -= 1;
        e -= 1;
    }
    ctx_log_ratio(pm, pb, pm2, pb2)
}

/// Spacer-mode RBS score: the best motif position among the RBSL spacer
/// lengths j (motif ending j bases before the start codon), scanned from
/// j = RBSL-1 down to 0. Ratio per spacer:
/// `(motif(u) * RBSD[j]) / (bg(u) * exp(-(W+j)/NDEC))`; a strictly larger
/// ratio wins (so ties keep the longer spacer), `best` starts at 0.0.
/// Sets `s.spacer` (-1 if no ratio > 0) and returns `ln best` (NEG if 0).
///
/// `first_wmer` = W-mer start for j = RBSL-1, `step` = +1 (forward: the
/// motif moves right as j shrinks) or -1 (reverse).
fn spacer_best(s: &mut Site, prep: &RbsPrep, first_wmer: i64, step: i64, l: i64, spacer_dist: &[f64], fwd: bool) -> f64 {
    s.spacer = -1;
    let mut best = 0.0f64;
    let mut wmer_pos = first_wmer;
    for j in (0..l).rev() {
        let u = prep.wmer[wmer_pos as usize];
        let j = j as usize;
        let v = if fwd {
            (prep.motif_fwd[u] * spacer_dist[j]) / (prep.bg_fwd[u] * prep.spacer_bg[j])
        } else {
            (prep.motif_rev[u] * spacer_dist[j]) / (prep.bg_rev[u] * prep.spacer_bg[j])
        };
        if best < v {
            s.spacer = j as i8;
            best = v;
        }
        wmer_pos += step;
    }
    if best != 0.0 {
        best.ln()
    } else {
        NEG
    }
}

/// Start term of every site (0x418bf0, viterbi.md 8): `st = ln P(start
/// type)` (`log_start[flags & 3]`; stops get the ATG value, never used),
/// then for starts with RBSS = W > 0 (only with `-r`) the RBS log ratio is
/// added:
///
/// * context mode (RBSL < 0, e.g. the GeneMarkS-T trained model W = 9,
///   RBSL = -3): [`ctx_fwd`] / [`ctx_rev`]; spacer stays -1;
/// * spacer mode (RBSL >= 0, prokaryotic RBS): [`spacer_best`].
///
/// Starts whose window would leave the sequence get +0.0. The window
/// tests are inclusive for forward starts and exclusive for reverse
/// starts, as in the binary.
///
/// `nuc` = codes after gap split, before R/Y resolution (image byte 0).
///
/// Intended differences (context mode, g = -RBSL; the binary reads memory
/// it does not own, the port adds 0.0 instead and reports it in the
/// returned [`RbsIssues`] so that main can warn):
///
/// * W < g: the rows for the bases after the start codon run below row 0,
///   i.e. in front of the motif matrix, for every start site;
/// * the range test keeps the window inside the sequence only if
///   W >= 2g + 3 (forward) and W >= 2g + 1 (reverse); otherwise starts
///   next to the sequence ends read image cells past the end / in front
///   of the start of the sequence.
///
/// Spacer mode with W > [`RBS_MAX_SPACER_W`]: no tables (see main), +0.0.
pub fn site_scores(sites: &mut [Site], nuc: &[u8], m: &Model, prep: Option<&RbsPrep>) -> RbsIssues {
    let mut issues = RbsIssues::default();
    for s in sites.iter_mut() {
        s.st = m.log_start[(s.flags & 3) as usize];
    }
    let w = m.rbss as i64;
    if w <= 0 {
        return issues;
    }
    let l = m.rbsl as i64;
    let g = -l;
    if l < 0 && w < g {
        issues.rows_missing = true;
        return issues;
    }
    if l >= 0 && prep.is_none() {
        issues.too_wide = true;
        return issues;
    }
    let n = nuc.len() as i64;
    let k: &Markov = &m.markov;
    let mut bgph = [0.0f64; 17];
    bgph[..4].copy_from_slice(&k.nonc0);
    for j in 0..3 {
        bgph[4 + 4 * j..8 + 4 * j].copy_from_slice(&k.cod1mono[j]);
    }
    // K+0x150: int order, int N, read as one double
    bgph[16] = f64::from_bits((k.order as u32 as u64) | ((k.n as u32 as u64) << 32));
    let cx = CtxMemory { mx: &m.rbsm, bgph };
    for s in sites.iter_mut() {
        let f = s.flags;
        if f & 0x20 != 0 {
            let base = s.st;
            let i3 = s.pos - 3; // first base of the start codon
            let r = if l < 0 {
                if i3 < w + l || (n - w) - l < i3 {
                    0.0
                } else if s.pos + g - 1 > n - 1 {
                    // window past the sequence end (W < 2g + 3)
                    issues.window_outside = true;
                    0.0
                } else {
                    ctx_fwd(&cx, nuc, s.pos, w, l)
                }
            } else if !(i3 < w + l || (n - w) - l < i3) {
                // 0x418bf0 forward spacer branch: inclusive bounds (unlike reverse)
                spacer_best(s, prep.unwrap(), i3 - l - w + 1, 1, l, &m.rbsd, true)
            } else {
                0.0
            };
            s.st = base + r;
        } else if f & 8 != 0 {
            let base = s.st;
            let i4 = s.pos + 3; // last base of the CAT codon
            let r = if l < 0 {
                if (i4 < (n - l) - w) && (w + l < i4) {
                    if s.pos - g + 1 < 0 {
                        // window before the sequence start (W < 2g + 1)
                        issues.window_outside = true;
                        0.0
                    } else {
                        ctx_rev(&cx, nuc, s.pos, w, l)
                    }
                } else {
                    0.0
                }
            } else if (i4 < (n - w) - l) && (w + l < i4) {
                spacer_best(s, prep.unwrap(), i4 + l, -1, l, &m.rbsd, false)
            } else {
                0.0
            };
            s.st = base + r;
        }
    }
    issues
}

/// RBS model problems met by [`site_scores`] (the binary reads memory it
/// does not own there; the port adds 0.0 to the affected start scores).
#[derive(Default, Clone, Copy)]
pub struct RbsIssues {
    /// context mode with RBSS < -RBSL (model level, all start sites)
    pub rows_missing: bool,
    /// context window outside the sequence for some start sites
    pub window_outside: bool,
    /// spacer mode with RBSS > RBS_MAX_SPACER_W (model level)
    pub too_wide: bool,
}

/// Largest motif width for the spacer-mode tables (4 tables of 4^W
/// doubles: 4^13 * 8 bytes = 512 MiB each). From W = 16 on the binary's
/// int size `4 << (2W - 2)` wraps (to 0 for W = 16) and its tables are
/// then read far out of bounds; W = 14, 15 need 2 and 8 GiB per table.
pub const RBS_MAX_SPACER_W: i32 = 13;

/// `-K` file text (0x418bf0 tail, output.md 6): one line per start site in
/// array order, `1-based codon position \t strand \t st`. Forward: first
/// base of the start codon (pos - 2); reverse: last base of CAT (pos + 4).
/// `pos >= 3` is the binary's test "byte offset into the image > 11".
pub fn k_file(sites: &[Site]) -> String {
    let mut o = String::new();
    for s in sites {
        let idx = s.pos;
        if s.flags & 0x20 != 0 && idx >= 3 {
            o.push_str(&format!("{}\t+\t{:.6}\n", idx - 2, s.st));
        }
        if s.flags & 0x08 != 0 && idx >= 3 {
            o.push_str(&format!("{}\t-\t{:.6}\n", idx + 4, s.st));
        }
    }
    o
}

// ---------------------------------------------------------------- DP

/// Command-line parameters of the DP (DP object slots, viterbi.md 6.1).
pub struct DpParams {
    /// `-h`: genes with `L - 1 < hmin` get `ev_score * -0.5` added to the
    /// COD2 (class 2) sum (default 0 = off).
    pub hmin: i64,
    /// `-E` evidence score (default 100.0); also the hmin penalty scale.
    pub ev_score: f64,
    /// `-p`: gene overlaps allowed (default true).
    pub allow_overlap: bool,
    /// `-x` overlap parameter q (default 0.26); `ln q` rewards codons that
    /// share bases, q also enters the overlap penalty.
    pub q: f64,
    /// `-i` probability p (default 0.5, per record with `-c`): ln p to end
    /// (and start) in noncoding state, ln(1-p) inside a gene.
    pub p: f64,
}

/// Per reverse frame state of [`Dp::rev_start`] (DP slots [0xd + 12f ..],
/// viterbi.md 6.6). The reverse ORF of frame f is scanned from its stop
/// (left) through its starts (rightwards); this state carries what can be
/// reused from one start to the next of the same ORF.
#[derive(Clone, Copy, Default)]
struct FrameState {
    /// Running coding sums of the ORF for set A (COD1), set B (COD2) and
    /// the third set (C, `-u` only; here always about NEG and never part of
    /// a score), covering
    /// the conditional entries [stop.pos, start.pos - k - 1] of the latest
    /// start.
    sum_a: f64,
    sum_b: f64,
    sum_c: f64,
    /// First forward stop of frame tag 0x40 / 0x80 / 0xc0 met in the ORF
    /// region (candidates for a forward gene overlapping this one).
    fwd_stops: [Option<usize>; 3],
    /// Overlap penalties of `fwd_stops`, computed once when the stop is
    /// found, with the gene length of the FIRST start of the ORF (they are
    /// not recomputed for later starts: literal).
    fwd_stop_pen: [f64; 3],
    /// Last reverse start of each other frame found right of the stop
    /// (candidates for an overlapping reverse gene).
    rev_starts: [Option<usize>; 3],
}

/// DP state (in main, viterbi.md 6.1).
struct Dp<'a> {
    /// The site array S, updated in place (score, back, cls, gs, ncs).
    s: &'a mut Vec<Site>,
    /// Markov image (see [`markov_image`]).
    idx: &'a Image,
    /// Codes after gap split, before R/Y resolution (image byte 0); used
    /// only to index the short-gap tables.
    nuc: &'a [u8],
    m: &'a Model,
    /// Markov order k (ORDM).
    k: i64,
    fs: [FrameState; 3],
    /// Termination mode (+0x1a1): true only for the final noncoding
    /// candidates into the end sentinel (0x415f40).
    term: bool,
    /// Class byte of the last gene evaluated (+0x60, uninitialised stack
    /// uninitialised data in the binary; 0 here). Copied into `cls` by every update
    /// except the forward overlap ones; only GE sites' values are ever read.
    cls: u8,
    /// "Other site" slot [9]: written as in the binary, read only as the
    /// reverse ORF root inside [`Dp::rev_start`].
    other_site: Option<usize>,
    /// Current site index (slot [8]); after the main loop the end sentinel.
    cur: usize,
    /// Position 2 (`image+8`): a start with pos <= p0 is a left boundary
    /// pseudo start, i.e. a partial gene.
    p0: i64,
    /// Position n-3: an end with pos >= p1 is a right boundary pseudo
    /// stop/start. Genes with `p0 < begin.pos && end.pos < p1` are complete
    /// and use the complete-gene length table, others the partial one.
    p1: i64,
    /// ln q (`-x`).
    logq: f64,
    par: &'a DpParams,
    /// Speed-up: indices of gene-end sites (flags & 0x18), ascending.
    ge_sites: Vec<usize>,
    /// Speed-up: next index >= i with a forward stop (flags & 0x10), or len.
    next_fstop: Vec<usize>,
    /// Speed-up: indices of reverse starts (flags & 8), ascending.
    rstarts: Vec<usize>,
}

/// DP update rule used almost everywhere: take the candidate if the cell
/// is still 0.0 ("unset") or strictly smaller. Ties keep the earlier
/// candidate. A real score of exactly 0.0 would be overwritten by any
/// candidate (never observed, kept literally).
#[inline]
fn better(cur: f64, s: f64) -> bool {
    cur == 0.0 || cur < s
}

/// Pointer comparison of the binary on site indices, with None (NULL)
/// smaller than every site.
#[inline]
fn ptr_lt(a: Option<usize>, b: Option<usize>) -> bool {
    a.map(|x| x as i64).unwrap_or(-1) < b.map(|x| x as i64).unwrap_or(-1)
}

impl Dp<'_> {
    /// Markov coding table j: 0..7 = set A (from COD1, `K+0x00..0x38`),
    /// 8..15 = set B (from COD2, `K+0x40..0x78`). Within a set (model.md
    /// 6.5): 0..2 forward conditional tables per codon position, 3 forward
    /// joint (k+1)-mer, 4..6 reverse-strand conditional tables (reverse
    /// complemented, permuted by k mod 3), 7 reverse joint.
    #[inline]
    fn t(&self, j: usize) -> &[f64] {
        let mk = &self.m.markov;
        if j < 8 {
            &mk.c1[j]
        } else {
            &mk.c2[j - 8]
        }
    }
    /// Noncoding length term `nl_len[d]` (model table M+0xa8) for
    /// `d = pos(GB) - pos(GE)` = noncoding length + 7: [0] = NEG,
    /// [1..=8] = ln(1-e), [8+i] = i ln e + ln(1-e) with e = exp(-1/NDEC).
    #[inline]
    fn nc_dur(&self, d: i64) -> f64 {
        self.m.dur_n[d as usize]
    }

    /// Overlap penalty (0x417690 and the inline copies in 0x4165c0): our
    /// gene runs from GB site `our_begin` to GE site `our_end`; the other
    /// gene from `other_begin` to `other_end` and overlaps our left part.
    /// With D = `(other_end.pos - our_begin.pos) + 7` (overlap length):
    ///
    /// `ln((1-q) * 0.2 * (D / L_ours + D / L_other) + q) * D - nl_len[7]`
    ///
    /// (L = gene lengths incl. codons). The caller subtracts it. A missing
    /// `other_begin` (NULL back pointer) gives 1e100.
    fn overlap_pen(&self, other_end: usize, our_begin: usize, our_end: usize, other_begin: Option<usize>) -> f64 {
        let s = &self.s;
        let q = self.par.q;
        let other_begin = match other_begin {
            None => return 1e100,
            Some(other_begin) => other_begin,
        };
        let overlap_len = ((s[other_end].pos - s[our_begin].pos) + 7) as f64;
        let x = ((1.0 - q) * 0.2 * (overlap_len / ((s[our_end].pos - s[our_begin].pos) + 7) as f64 + overlap_len / ((s[other_end].pos - s[other_begin].pos) + 7) as f64) + q).ln();
        x * overlap_len - self.nc_dur(7)
    }

    /// Noncoding log likelihood of a long gap (d > k+7) from GE site `left`
    /// to GB site `right`: the gap is [left.pos+4, right.pos-4]; first the
    /// joint (k+1)-mer `nl[k]` ending at right.pos-4, then the conditional
    /// NONC terms right to left down to left.pos+4 (this order of
    /// additions).
    fn nc_sum_long(&self, left: usize, right: usize) -> f64 {
        let s = &self.s;
        let nc = &self.m.markov.nc;
        let last_nc_kmer = s[right].pos - self.k - 4;
        let mut acc = self.m.markov.nl[self.k as usize][self.idx.at(last_nc_kmer)];
        // the guard is redundant with the while condition (kept from the binary)
        if s[left].pos + 3 < last_nc_kmer - 1 {
            let mut e = last_nc_kmer - 1;
            while e >= s[left].pos + 4 {
                acc += nc[self.idx.at(e)];
                e -= 1;
            }
        }
        acc
    }

    /// Index of a short gap (1..k bases, [left.pos+4, right.pos-4]) into
    /// the lower-order joint table `nl[len-1]`: base-4 number, first base
    /// most significant, built from the raw codes (an N is digit 4, so the
    /// index may run past the table, see [`Dp::nc_candidates`]).
    fn short_gap_index(&self, left: usize, right: usize) -> usize {
        let mut v = 0usize;
        for e in self.s[left].pos + 4..self.s[right].pos - 3 {
            v = self.nuc[e as usize] as usize + v * 4;
        }
        v
    }

    /// Noncoding transitions into GB site `c` (0x4165c0 GB branch, also
    /// used by 0x415f40 with `term`), viterbi.md 6.4: candidates are GE
    /// sites g at index <= `start`, scanned right to left;
    /// `d = c.pos - g.pos` (noncoding length + 7):
    ///
    /// * `d < 1`: skipped (g lies right of c; the array is only roughly
    ///   sorted).
    /// * `1 <= d <= 7`: the two codons overlap or touch:
    ///   `(g.score + nl_len[d]) - (7-d) ln q` (ln q < 0, so shared bases
    ///   earn a bonus).
    /// * `8 <= d <= k+7` (1..k noncoding bases): `(g.score + nl_len[d]) +
    ///   nl[d-8][gap]`. Scanning continues to the left after these.
    /// * `d > k+7`: the first such g gets the full sum ([`Dp::nc_sum_long`]);
    ///   then the scan extends the same sum leftwards through the remaining
    ///   GE sites by adding their precomputed `ncs` segments, until
    ///   `d > MAXN`, and stops (`return`).
    ///
    /// With `term` (end sentinel) every non-overlap candidate gets
    /// `(s + nl_len[9]) + nl_len[8] * -2.0`.
    ///
    /// Speed trick: the binary walks every site from `start` down and
    /// skips non-GE sites; walking the GE list visits the same sites in
    /// the same order.
    fn nc_candidates(&mut self, c: usize, start: i64) {
        let k = self.k;
        let mut gi = self.ge_sites.partition_point(|&x| (x as i64) <= start);
        while gi > 0 {
            gi -= 1;
            let g_end = self.ge_sites[gi];
            let d = self.s[c].pos - self.s[g_end].pos;
            if d <= k + 7 {
                let cand;
                if d < 8 {
                    if d < 1 {
                        continue;
                    }
                    cand = (self.s[g_end].score + self.nc_dur(d)) - (7 - d) as f64 * self.logq;
                } else {
                    let tab = &self.m.markov.nl[(d - 8) as usize];
                    let w = self.short_gap_index(g_end, c);
                    // N in the gap reads past the table in the binary (undefined read); use the fill slot.
                    // (An N that keeps the index inside the table reads the same
                    // entry as the binary, so that case is exact.)
                    if w >= tab.len() && std::env::var_os("GMHMMP_RS_DEBUG").is_some() {
                        eprintln!("short-gap read past table: d={d} w={w} len={} g.pos={} r.pos={}", tab.len(), self.s[g_end].pos, self.s[c].pos);
                    }
                    let v = tab.get(w).copied().unwrap_or(*tab.last().unwrap());
                    let mut x = (self.s[g_end].score + self.nc_dur(d)) + v;
                    if self.term {
                        x = (x + self.nc_dur(9)) + self.nc_dur(8) * -2.0;
                    }
                    cand = x;
                }
                // 0x417367/0x417371: `cur == 0` takes, else `ucomisd` + `jbe`
                // skips `cand <= cur` AND unordered, i.e. better() (a NaN
                // candidate is not taken; `-x nan`)
                if better(self.s[c].score, cand) {
                    self.s[c].back = Some(g_end);
                    self.s[c].score = cand;
                    self.s[c].cls = self.cls;
                }
                continue;
            }
            // first long gap: explicit sum
            let mut acc = self.nc_sum_long(g_end, c);
            let mut cand = (self.s[g_end].score + self.nc_dur(d)) + acc;
            if self.term {
                cand = (cand + self.nc_dur(9)) + self.nc_dur(8) * -2.0;
            }
            if better(self.s[c].score, cand) {
                self.s[c].back = Some(g_end);
                self.s[c].score = cand;
                self.s[c].cls = self.cls;
            }
            // extend the same noncoding stretch leftwards: g_end2.ncs covers
            // [g_end2.pos+4, next GE pos+3], so acc grows by one segment per site
            let maxn = self.m.maxn as i64;
            let cpos = self.s[c].pos;
            while gi > 0 {
                gi -= 1;
                let g_end2 = self.ge_sites[gi];
                let d = cpos - self.s[g_end2].pos;
                if maxn < d {
                    break;
                }
                acc += self.s[g_end2].ncs;
                let mut cand = (self.s[g_end2].score + self.nc_dur(d)) + acc;
                if self.term {
                    cand = (cand + self.nc_dur(9)) + self.nc_dur(8) * -2.0;
                }
                if better(self.s[c].score, cand) {
                    self.s[c].back = Some(g_end2);
                    self.s[c].score = cand;
                    self.s[c].cls = self.cls;
                }
            }
            return;
        }
    }

    /// 0x415e70 (viterbi.md 6.3): clear `back` of the sites with
    /// `pos - S[0].pos < 7` (the left boundary pseudo starts at pos 0..2,
    /// which keep their initial scores and end every traceback), give a GB
    /// site at pos 3 (a codon at bases 0..2) the score of the left sentinel
    /// directly, and set `cur` to the first site the main loop visits.
    fn init(&mut self) {
        let s = &mut self.s;
        let base = s[0].pos;
        let mut c = 0usize;
        let mut d = s[0].pos - base;
        if d < 7 {
            loop {
                s[c].back = None;
                c += 1;
                d = s[c].pos - base;
                if d >= 7 {
                    break;
                }
            }
        }
        if s[c].flags & 0x24 != 0 && d == 7 {
            // `+ 0.0` is literal (only turns -0.0 into +0.0)
            let v = s[0].score + 0.0;
            if better(s[c].score, v) {
                s[c].score = v;
                s[c].back = Some(0);
                s[c].cls = self.cls;
            }
            c += 1;
        }
        self.cur = c;
    }

    /// 0x415d70 (viterbi.md 6.3): noncoding segment sums `ncs`. For
    /// consecutive GE sites g < g', `g.ncs` = sum of NC over
    /// [g.pos+4, g'.pos+3] (left to right). A running cursor makes sure no
    /// base is counted twice when sites are out of order (then the segment
    /// is empty, 0.0). Non-GE sites, the last GE site and every site within
    /// 7 of the end sentinel get 0.0. Used by [`Dp::nc_candidates`] to
    /// extend long noncoding sums site by site.
    fn nc_sums(&mut self) {
        let nc = &self.m.markov.nc;
        let idx = self.idx;
        let s = &mut self.s;
        let mut k = 0;
        while s[k].flags & 0x18 == 0 {
            s[k].ncs = 0.0;
            k += 1;
        }
        let mut prev_ge = k;
        let last = s.len() - 1;
        let lastpos = s[last].pos;
        // j = last site with pos <= lastpos - 7; sites after it get 0.0
        let mut j = last;
        loop {
            s[j].ncs = 0.0;
            j -= 1;
            if lastpos - s[j].pos >= 7 {
                break;
            }
        }
        let mut next_base = s[prev_ge].pos + 4;
        for m in prev_ge + 1..=j {
            if s[m].flags & 0x18 == 0 {
                s[m].ncs = 0.0;
            } else {
                let mut acc = 0.0;
                let e = s[m].pos + 3;
                if next_base <= e {
                    for x in next_base..=e {
                        acc += nc[idx.at(x)];
                    }
                    next_base = e + 1;
                }
                s[prev_ge].ncs = acc;
                prev_ge = m;
            }
        }
    }

    /// Evaluate one forward gene (start `start_site`, stop `c`) given its
    /// coding sums `a` (set A) and `b` (set B) and `d = L - 1`
    /// (0x4165c0, EVAL in viterbi.md 6.5):
    ///
    /// * `-h`: `b += E * -0.5` if `d < hmin`. The penalised `b` is returned
    ///   and the caller keeps extending it, so the penalty stays for all
    ///   starts further left (literal);
    /// * class: '2' if `!(a >= b)` (also for NaN), else '1' (ties go to class 1); `s` = the
    ///   winner, stored as the start's gene score `gs`;
    /// * `t = (s + len) + st` with the complete or partial gene length
    ///   table and the start term of `start_site`;
    /// * candidate `start.score + t` for the stop (back = start), and with
    ///   `-p 1` the overlap candidates ([`Dp::fwd_overlap`]).
    fn fwd_eval(&mut self, c: usize, start_site: usize, a: f64, mut b: f64, d: i64) -> (f64, f64) {
        if d < self.par.hmin {
            b += self.par.ev_score * -0.5;
        }
        let s;
        // `ucomisd` + `jb` (0x416829, 0x416d19): a NaN a or b (e.g. `-E nan`
        // with `-h`) gives class 2, so the test is `!(a >= b)`
        if !(a >= b) {
            self.cls = 0x32;
            s = b;
        } else {
            s = a;
            self.cls = 0x31;
        }
        self.s[start_site].gs = s;
        let len = if self.p0 < self.s[start_site].pos && self.s[c].pos < self.p1 { self.m.dur_c[d as usize] } else { self.m.dur_p[d as usize] };
        let t = (s + len) + self.s[start_site].st;
        let tot = self.s[start_site].score + t;
        if better(self.s[c].score, tot) {
            self.s[c].back = Some(start_site);
            self.s[c].score = tot;
            self.s[c].cls = self.cls;
        }
        if self.par.allow_overlap {
            self.fwd_overlap(c, start_site, t);
        }
        (a, b)
    }

    /// Overlap candidates for the forward gene `start_site`..`c` (viterbi.md
    /// 6.5 fwd_overlap). Instead of entering the gene from the noncoding
    /// path at the start, the path comes from another gene's right end j
    /// lying inside this gene: candidate `t + (j.score - penalty)`. When
    /// one wins, `c.back = start` and the START's back pointer is re-routed
    /// to j (the traceback follows it); `cls` is not changed here.
    ///
    /// (a) forward stops j in (start, c): the first one, then the first
    ///     one of a different frame after it, i.e. the first stop of each
    ///     of the two other frames (the binary looks at no others). Taken
    ///     only if `j.back < start` (a NULL back counts as smaller and gives
    ///     penalty 1e100, which still enters the accumulated `t`).
    /// (b) every reverse start j in (start, c) with `j.back < start`.
    ///
    /// Kept literally as in the original [V]: `t` is NOT reset between candidates, each
    /// candidate adds its term to the previous one. Both comparisons are
    /// `t > score` (`ucomisd` + `jbe`: a NaN never wins, checked with
    /// `-x nan`).
    ///
    /// Speed trick: `next_fstop` and `rstarts` visit exactly the sites of
    /// the binary's linear scans, in the same order.
    fn fwd_overlap(&mut self, c: usize, start_site: usize, mut t: f64) {
        let mut hit = start_site;
        let mut hit_frame: i64 = -1;
        loop {
            let prev_frame = hit_frame;
            let mut j = hit + 1;
            let mut found = false;
            while j < c {
                j = self.next_fstop[j];
                if j >= c {
                    break;
                }
                let fr = (self.s[j].flags >> 6) as i64 - 1;
                if fr != prev_frame {
                    hit_frame = fr;
                    found = true;
                    break;
                }
                j += 1;
            }
            if !found {
                break;
            }
            hit = j;
            let o = self.s[hit].back;
            if !ptr_lt(o, Some(start_site)) {
                if prev_frame != -1 {
                    break;
                }
                continue;
            }
            let other_score = self.s[hit].score;
            let pen = self.overlap_pen(hit, start_site, c, o);
            t += other_score - pen;
            // `ucomisd` + `jbe` (0x416a68, 0x416f58): skipped if `t <= score`
            // or unordered, so a NaN t (or score) is never taken (`-x nan`)
            if !(t > self.s[c].score) {
                if prev_frame != -1 {
                    break;
                }
                continue;
            }
            self.s[c].back = Some(start_site);
            self.s[c].score = t;
            self.s[start_site].back = Some(hit);
            if prev_frame != -1 {
                break;
            }
        }
        let lo = self.rstarts.partition_point(|&x| x <= start_site);
        let hi = self.rstarts.partition_point(|&x| x < c);
        for ri in lo..hi {
            let p = self.rstarts[ri];
            if ptr_lt(self.s[p].back, Some(start_site)) {
                let other_score = self.s[p].score;
                let pen = self.overlap_pen(p, start_site, c, self.s[p].back);
                t += other_score - pen;
                // same as `score < t`
                if self.s[c].score <= t && t != self.s[c].score {
                    self.s[c].back = Some(start_site);
                    self.s[c].score = t;
                    self.s[start_site].back = Some(p);
                }
            }
        }
    }

    /// Forward genes ending at forward stop `c` (0x4165c0 branch for flags
    /// & 0x10, viterbi.md 6.5). Walks the ORF chain from the nearest start
    /// to the farthest, extending the coding sums incrementally, and
    /// evaluates every start with [`Dp::fwd_eval`].
    ///
    /// Coding sum for model X (A or B), right to left over [b, e]
    /// (e = c.pos, b = start pos), u = e - k:
    /// `prior + X[3][idx[u]]` (joint (k+1)-mer ending at e), then 0, 1 or 2
    /// single conditional terms depending on k mod 3 so that the remaining
    /// part is whole codons, then per codon (m = first base of the codon
    /// to the right): `x += X[2][m-1]; x = (x + X[1][m-2]) + X[0][m-3]`.
    /// Table X[p] scores a base at codon position p+1. The loop ends exactly
    /// at b because b and e are whole codons apart.
    ///
    /// Evidence (`-e` only): a flagged stop with `ev == E` adds E to both
    /// priors; with any other `ev` the sequence is ignored and the sums
    /// become `cd + ev`.
    fn fwd_stop(&mut self, c: usize) {
        let k = self.k;
        let idx = self.idx;
        let stop_pos = self.s[c].pos;
        let mut start_site = self.s[c].link.expect("stop without start");
        let first_start_pos = self.s[start_site].pos;
        let ln_cd1p = self.m.cd1p;
        let ln_cd2p = self.m.cd2p;
        let mut d = (stop_pos - first_start_pos) + 6;
        let u = stop_pos - k;
        let ev_flag = self.s[c].evf;
        let ev0 = self.par.ev_score;
        let use_ev = (self.s[c].flags & 0x10 != 0) && ev_flag != 0 && ev0 == self.s[c].ev;
        // note the operand order of the B prior (E + cB), as in the binary
        let (mut a, b0) = if use_ev { (ln_cd1p + ev0, ev0 + ln_cd2p) } else { (ln_cd1p, ln_cd2p) };
        let i = idx.at(u);
        a += self.t(3)[i];
        let mut b = self.t(0xb)[i] + b0;
        let k_mod3 = k % 3;
        // e = first base of the next codon to the right of the part left to sum
        let mut e;
        if k_mod3 == 0 {
            let i1 = idx.at(u - 1);
            a += self.t(1)[i1];
            b += self.t(9)[i1];
            let i2 = idx.at(u - 2);
            a += self.t(0)[i2];
            b += self.t(8)[i2];
            e = u - 2;
        } else if k_mod3 == 1 {
            e = u - 1;
            let i1 = idx.at(u - 1);
            a += self.t(0)[i1];
            b += self.t(8)[i1];
        } else {
            e = u;
        }
        // whole codons down to the nearest start
        while e > first_start_pos {
            let (i1, i2, i3) = (idx.at(e - 1), idx.at(e - 2), idx.at(e - 3));
            a += self.t(2)[i1];
            b += self.t(0xa)[i1];
            a = (a + self.t(1)[i2]) + self.t(0)[i3];
            b = (b + self.t(9)[i2]) + self.t(8)[i3];
            e -= 3;
        }
        if !(ev_flag == 0 || self.s[c].ev == ev0) {
            b = ln_cd2p + self.s[c].ev;
            a = ln_cd1p + self.s[c].ev;
        }
        let (na, nb) = self.fwd_eval(c, start_site, a, b, d);
        a = na;
        b = nb;
        // further starts to the left: extend the same sums codon by codon
        while let Some(nx) = self.s[start_site].link {
            let mut e = self.s[start_site].pos;
            start_site = nx;
            d = (stop_pos - self.s[nx].pos) + 6;
            while self.s[nx].pos < e {
                let (i1, i2, i3) = (idx.at(e - 1), idx.at(e - 2), idx.at(e - 3));
                a += self.t(2)[i1];
                b += self.t(0xa)[i1];
                a = (a + self.t(1)[i2]) + self.t(0)[i3];
                b = (b + self.t(9)[i2]) + self.t(8)[i3];
                e -= 3;
            }
            let (na, nb) = self.fwd_eval(c, start_site, a, b, d);
            a = na;
            b = nb;
        }
    }

    /// Reverse gene ending (in forward coordinates) at reverse start `c`
    /// (0x417760, viterbi.md 6.6). The gene runs from the ORF's reverse
    /// stop ("root", left, GB) to `c` (right, GE).
    ///
    /// Reverse starts of one ORF are visited left to right, so the coding
    /// sums are kept per frame in [`FrameState`]:
    ///
    /// * `c.link` is a start of the same ORF (incremental branch): add the
    ///   codons between the previous start and `c` to the frame's sums.
    ///   The loop starts at `c.pos - k`, so the three tables are used in an
    ///   order rotated by k mod 3: (6,4,5), (4,5,6), (5,6,4) for k%3 = 1, 2,
    ///   0, which equals the fresh sum's (4,5,6) pattern shifted.
    /// * `c.link` is the reverse stop (first start of the ORF, fresh
    ///   branch): reset the frame state, collect the overlap candidates of
    ///   this ORF (first forward stop of each frame, with penalties, and the
    ///   last reverse start of each other frame), then sum the coding region
    ///   from `c.pos` down to the stop.
    ///
    /// Both sums cover the conditional entries [stop.pos, c.pos - k - 1];
    /// the joint table X[7] at `c.pos - k` (right end) is added in the
    /// common tail. Then class, gene score, length and start terms as in
    /// the forward case, but with `t = (len + st) + s` (different
    /// association) and without carrying the hmin penalty over. With
    /// `-p 1` overlap candidates follow; unlike the forward case `t` is not
    /// accumulated and `cls` is set.
    fn rev_start(&mut self, c: usize) {
        let k = self.k;
        let k_mod3 = k % 3;
        let idx = self.idx;
        // the chain predecessor of c; after both branches below it is the
        // ORF's reverse stop (== rev_root)
        let mut left_site = self.s[c].link.expect("rev start without link");
        self.other_site = Some(left_site);
        let rev_root;
        if self.s[left_site].flags & 4 == 0 {
            // ---- incremental: left_site is the previous start of this ORF
            let mut e = self.s[c].pos - k;
            let lim = self.s[left_site].pos - k;
            let (mut a, mut b) = (0.0f64, 0.0f64);
            let tabs: (usize, usize, usize) = match k_mod3 {
                1 => (6, 4, 5),
                2 => (4, 5, 6),
                _ => (5, 6, 4),
            };
            while lim < e {
                let (i1, i2, i3) = (idx.at(e - 1), idx.at(e - 2), idx.at(e - 3));
                a += self.t(tabs.0)[i1];
                b += self.t(tabs.0 + 8)[i1];
                a = (a + self.t(tabs.1)[i2]) + self.t(tabs.2)[i3];
                b = (b + self.t(tabs.1 + 8)[i2]) + self.t(tabs.2 + 8)[i3];
                e -= 3;
            }
            let f = ((self.s[left_site].flags >> 6) - 1) as usize;
            let fs = &mut self.fs[f];
            // the C sum gets NEG added (no -u); operand order as in the binary
            fs.sum_c = NEG + fs.sum_c;
            fs.sum_a = a + fs.sum_a;
            fs.sum_b = b + fs.sum_b;
            // find the ORF root (reverse stop) through the links
            match self.s[left_site].link {
                None => rev_root = self.other_site.unwrap(),
                Some(mut r) => {
                    while let Some(l) = self.s[r].link {
                        r = l;
                    }
                    self.other_site = Some(r);
                    rev_root = r;
                }
            }
            left_site = rev_root;
        } else {
            // ---- fresh: left_site is the reverse stop, c the first start
            let f = ((self.s[left_site].flags >> 6) - 1) as usize;
            // the binary keeps the old penalties; they are only read
            // together with a fwd_stops entry, which is reset here
            self.fs[f] = FrameState::default();
            let cnt = self.s.len();
            // (ov) forward stops in the ORF region: start just left of the
            // sites sharing the stop's pos, scan right until all three
            // frames are found, the array ends, or the next reverse stop of
            // this frame (flags equal to the stop's flags)
            let mut scan = left_site;
            while scan > 0 {
                scan -= 1;
                if self.s[scan].pos != self.s[left_site].pos {
                    break;
                }
            }
            // always true here (scan <= left_site); kept from the binary
            if scan < cnt && (scan <= left_site || self.s[scan].flags != self.s[left_site].flags) {
                while !(self.fs[f].fwd_stops[0].is_some() && self.fs[f].fwd_stops[1].is_some() && self.fs[f].fwd_stops[2].is_some()) {
                    let scan_flags = self.s[scan].flags;
                    if scan_flags & 0x10 != 0 && self.s[left_site].pos <= self.s[scan].pos {
                        // exact flag compare: forward stop of frame tag 2, 3, 1
                        for (kk, fl) in [(1usize, 0x90u8), (2, 0xd0), (0, 0x50)] {
                            if scan_flags == fl && self.fs[f].fwd_stops[kk].is_none() {
                                self.fs[f].fwd_stops[kk] = Some(scan);
                                // computed with c = the FIRST start of the ORF
                                let pen = self.overlap_pen(scan, left_site, c, self.s[scan].back);
                                self.fs[f].fwd_stop_pen[kk] = pen;
                            }
                        }
                    }
                    scan += 1;
                    if scan >= cnt || (scan > left_site && self.s[scan].flags == self.s[left_site].flags) {
                        break;
                    }
                }
            }
            // (rs) reverse starts of other frames right of the stop: the
            // first reverse stop met fixes stop_frame; starts of that frame
            // (and of ours) are ignored after it, a reverse stop of a third
            // frame ends the scan; the last start per frame wins
            let mut stop_frame: i64 = -1;
            let mut rs_scan = left_site;
            loop {
                rs_scan += 1;
                if rs_scan >= cnt {
                    break;
                }
                let rs_flags = self.s[rs_scan].flags;
                if rs_flags == self.s[left_site].flags {
                    break;
                }
                // u8 arithmetic of the binary: frame tag 0 gives 255
                let rs_frame = (((rs_flags >> 6) as i64 - 1) & 0xff) as i64;
                let next_stop_frame;
                if rs_flags & 4 != 0 {
                    if !(stop_frame == -1 || rs_frame == stop_frame) {
                        break;
                    }
                    next_stop_frame = rs_frame;
                } else {
                    next_stop_frame = stop_frame;
                }
                stop_frame = next_stop_frame;
                if rs_flags & 8 != 0 && f as i64 != rs_frame && rs_frame != stop_frame && rs_frame < 3 {
                    self.fs[f].rev_starts[rs_frame as usize] = Some(rs_scan);
                }
            }
            // fresh coding sum from c.pos down to the stop (cf. fwd_stop)
            let ln_cd2p = self.m.cd2p;
            let ln_cd1p = self.m.cd1p;
            let e0 = self.s[c].pos - k;
            let cc = ln_cd2p - 1e100;
            let (mut a, mut b);
            let mut e;
            if k_mod3 == 0 {
                let i1 = idx.at(e0 - 1);
                a = self.t(5)[i1] + ln_cd1p;
                b = self.t(13)[i1] + ln_cd2p;
                let i2 = idx.at(e0 - 2);
                e = e0 - 2;
                b += self.t(14)[i2];
                a += self.t(6)[i2];
            } else if k_mod3 == 1 {
                e = e0 - 1;
                let i1 = idx.at(e0 - 1);
                a = self.t(6)[i1] + ln_cd1p;
                b = self.t(14)[i1] + ln_cd2p;
            } else {
                e = e0;
                b = ln_cd2p;
                a = ln_cd1p;
            }
            while self.s[left_site].pos < e {
                let (i1, i2, i3) = (idx.at(e - 1), idx.at(e - 2), idx.at(e - 3));
                a += self.t(4)[i1];
                b += self.t(12)[i1];
                a = (a + self.t(5)[i2]) + self.t(6)[i3];
                b = (b + self.t(13)[i2]) + self.t(14)[i3];
                e -= 3;
            }
            // evidence on the reverse stop (-e only). The replacement uses
            // c's ev (the start's, = E without evidence), not the stop's;
            // read from the binary [R], kept literally.
            if self.s[left_site].flags & 4 != 0 && self.s[left_site].evf != 0 {
                let ev = self.s[left_site].ev;
                if ev == self.par.ev_score {
                    a += self.par.ev_score;
                    b += self.par.ev_score;
                } else {
                    b = ln_cd2p + self.s[c].ev;
                    a = ln_cd1p + self.s[c].ev;
                }
            }
            self.fs[f].sum_c = cc;
            self.fs[f].sum_a = a;
            self.fs[f].sum_b = b;
            rev_root = left_site;
        }
        // ---- common tail: joint term, class, gene score, length, update
        let f = ((self.s[left_site].flags >> 6) - 1) as usize;
        let ki = idx.at(self.s[c].pos - k);
        let a = (self.fs[f].sum_a + 0.0) + self.t(7)[ki];
        let mut b = (self.fs[f].sum_b + 0.0) + self.t(15)[ki];
        let d = (self.s[c].pos - self.s[left_site].pos) + 6;
        if d < self.par.hmin {
            b += self.par.ev_score * -0.5;
        }
        let s;
        // `ucomisd` + `jb` (0x417fdd): a NaN a or b gives class 2
        if !(a >= b) {
            self.cls = 0x32;
            s = b;
        } else {
            self.cls = 0x31;
            s = a;
        }
        self.s[c].gs = s;
        let len = if self.p0 < self.s[rev_root].pos && self.s[c].pos < self.p1 { self.m.dur_c[d as usize] } else { self.m.dur_p[d as usize] };
        let t = (len + self.s[c].st) + s;
        let tot = self.s[rev_root].score + t;
        if better(self.s[c].score, tot) {
            self.s[c].score = tot;
            self.s[c].back = Some(rev_root);
            self.s[c].cls = self.cls;
        }
        if self.par.allow_overlap {
            let root_opt = self.other_site;
            // forward genes ending inside this ORF: only frame tags 1 and 2
            // (fwd_stops[0], [1]); the binary checks [2] only with -u
            for kk in 0..2 {
                if let Some(fstop) = self.fs[f].fwd_stops[kk] {
                    if fstop < c && self.s[fstop].back.is_some() && ptr_lt(self.s[fstop].back, root_opt) {
                        self.other_site = Some(fstop);
                        if self.s[fstop].pos < self.p1 {
                            let s2 = (self.s[fstop].score + t) - self.fs[f].fwd_stop_pen[kk];
                            if better(self.s[c].score, s2) {
                                self.s[c].score = s2;
                                self.s[c].back = Some(fstop);
                                self.s[c].cls = self.cls;
                            }
                        }
                    }
                }
            }
            // reverse genes of other frames: from the recorded start walk its
            // link chain leftwards, skip sites at or right of c, stop left of
            // our root
            for kk in 0..3 {
                let mut u = self.fs[f].rev_starts[kk];
                while let Some(uu) = u {
                    if ptr_lt(Some(uu), root_opt) {
                        break;
                    }
                    if uu >= c {
                        u = self.s[uu].link;
                        continue;
                    }
                    self.other_site = Some(uu);
                    let pen = self.overlap_pen(uu, root_opt.unwrap(), c, self.s[uu].back);
                    let s2 = (self.s[uu].score + t) - pen;
                    if !better(self.s[c].score, s2) {
                        u = self.s[uu].link;
                        continue;
                    }
                    self.s[c].score = s2;
                    self.s[c].back = Some(uu);
                    self.s[c].cls = self.cls;
                    u = self.s[uu].link;
                }
            }
        }
    }

    /// Termination (0x415f40, viterbi.md 6.7) into the end sentinel `cur`
    /// (last site, GB, pos n+3).
    ///
    /// 1. The sites at the right end with pos >= n-3 (the right boundary
    ///    pseudo stops and starts, genes running off the end) get `ln(1-p)`
    ///    added, after
    ///    `-ln 16` / `-ln 4` for the 2 / 1 bases left unexplained at pos
    ///    n-3 / n-2. These sites are never linked to the sentinel; they
    ///    compete with it directly in [`best_end`].
    /// 2. Noncoding candidates into the sentinel ([`Dp::nc_candidates`]
    ///    with `term`) from the GE sites left of them, except that a real
    ///    gene end at pos n-4 (stop codon ending at the last base) is
    ///    skipped there and offered afterwards as `ln(1-p) + score`.
    /// 3. The sentinel's noncoding score gets `ln p` added before that last
    ///    comparison.
    fn terminate(&mut self) {
        let last = self.cur;
        let p = self.par.p;
        // `ucomisd` + `jb` at 0x415f66 / 0x415f84: NaN takes the log
        // branch (see build_temp_sites), so `-i nan` gives NaN for both
        let lp = if p > 0.0 || p.is_nan() { p.ln() } else { NEG };
        let l1p = if (1.0 - p) > 0.0 || (1.0 - p).is_nan() { (1.0 - p).ln() } else { NEG };
        let p1 = self.p1;
        let mut k = last as i64 - 1;
        let mut end_site = k;
        if p1 <= self.s[k as usize].pos {
            loop {
                let ku = k as usize;
                if self.s[ku].pos == p1 {
                    self.s[ku].score -= LN16;
                } else if self.s[ku].pos == p1 + 1 {
                    self.s[ku].score -= LN4;
                }
                self.s[ku].score += l1p;
                end_site = k - 1;
                k -= 1;
                if !(p1 <= self.s[k as usize].pos) {
                    break;
                }
            }
        }
        // scan_from = first site from the right with pos < n-3
        let mut scan_from = last as i64 - 1;
        if p1 <= self.s[scan_from as usize].pos {
            loop {
                scan_from -= 1;
                if !(p1 <= self.s[scan_from as usize].pos) {
                    break;
                }
            }
        }
        if self.s[scan_from as usize].flags & 0x18 != 0 && self.s[scan_from as usize].pos == p1 - 1 {
            // complete gene ending at the last base: not a noncoding candidate
            end_site = scan_from;
            scan_from -= 1;
        }
        self.term = true;
        if scan_from >= 0 {
            self.nc_candidates(last, scan_from);
        }
        self.s[last].score = lp + self.s[last].score;
        let end_u = end_site as usize;
        if self.s[end_u].flags & 0x18 != 0 && self.s[end_u].pos == p1 - 1 {
            let v = l1p + self.s[end_u].score;
            if better(self.s[last].score, v) {
                self.s[last].score = v;
                self.s[last].back = Some(end_u);
                self.s[last].cls = self.cls;
            }
        }
    }

    /// 0x417660: init, noncoding segment sums, one left-to-right pass over
    /// the sites (each site's score depends only on sites to its left), and
    /// termination. The end sentinel is not visited by the main loop.
    fn run(&mut self) {
        self.init();
        self.nc_sums();
        let last = self.s.len() - 1;
        let mut c = self.cur;
        while c < last {
            self.cur = c;
            let f = self.s[c].flags;
            if f & 0x10 != 0 {
                self.fwd_stop(c);
            } else if f & 0x08 != 0 {
                self.rev_start(c);
            } else {
                // GB site: forward start or reverse stop
                self.nc_candidates(c, c as i64 - 1);
            }
            c += 1;
        }
        self.cur = c;
        self.terminate();
    }
}

/// One predicted gene (GeneRec, output.md 2.4), produced by [`traceback`].
#[derive(Debug, Clone)]
pub struct Gene {
    /// Running gene number over the whole input file (0x412640), set in main.
    pub id: i32,
    /// 1-based inclusive left end (start codon of '+', stop codon of '-').
    pub left: i64,
    /// 1-based inclusive right end.
    pub right: i64,
    /// b'+' or b'-'.
    pub strand: u8,
    /// Class byte b'1' / b'2' (printed in the LST Class column).
    pub cls: u8,
    /// Left / right partial marks: b' ', or b'<' / b'>' (gene runs off the
    /// sequence end or ends on N).
    pub lm: u8,
    pub rm: u8,
    /// RBS spacer of the start site (-1 none); replaced by the RBS position
    /// with `-r -K` (`output::k_positions`).
    pub spacer: i8,
    /// RBS score of the start: `(st - ln P(start type)) + (spacer + W) / NDEC`.
    pub rbs: f64,
    /// Evidence-supported end (`-e` only); GFF `stop_enforced`.
    pub evf: u8,
    /// Gene score (start site `gs`); GFF column 6 and `gene_score`.
    pub gs: f64,
    /// Flags of the start site; GFF `start_codon` = flags & 3.
    pub sflags: u8,
    /// Log-odds score ([`logodd`]).
    pub logodd: f64,
    /// Maximum logodd of the record (for `-b`), set by the output.
    pub maxlo: f64,
}

/// Run the DP over the site map of one record (0x4164a0 + main stores +
/// 0x417660). `sites` is updated in place; read the result with
/// [`traceback`]. `idx` = Markov image, `nuc` = codes before R/Y
/// resolution, `m` = the record's model with length tables already sized
/// for this record.
pub fn run_dp(sites: &mut Vec<Site>, idx: &Image, nuc: &[u8], m: &Model, par: &DpParams) {
    let n = idx.v.len() as i64;
    let mut dp = Dp {
        s: sites,
        idx,
        nuc,
        m,
        k: m.ordm as i64,
        fs: [FrameState::default(); 3],
        term: false,
        cls: 0,
        other_site: None,
        cur: 0,
        p0: 2,
        p1: n - 3,
        logq: par.q.ln(),
        par,
        ge_sites: Vec::new(),
        next_fstop: Vec::new(),
        rstarts: Vec::new(),
    };
    // speed-up lists (see the module docs)
    let len = dp.s.len();
    dp.ge_sites = (0..len).filter(|&i| dp.s[i].flags & 0x18 != 0).collect();
    dp.rstarts = (0..len).filter(|&i| dp.s[i].flags & 8 != 0).collect();
    let mut nf = vec![len; len + 1];
    for i in (0..len).rev() {
        nf[i] = if dp.s[i].flags & 0x10 != 0 { i } else { nf[i + 1] };
    }
    dp.next_fstop = nf;
    dp.run();
}

/// Best end of the parse (0x412310): the end sentinel, unless a site with
/// pos > n-4 (right boundary pseudo sites, genes running off the end) has
/// a strictly higher score. Scanned right to left, so ties keep the
/// sentinel or the rightmost site. The test is `ucomisd` + `jbe` (skip if
/// `score <= best` or unordered, 0x412377), so a NaN score is never taken
/// and a NaN best is never replaced (Ghidra shows it as `score <= best`;
/// verified with `-i nan`).
fn best_end(s: &[Site], n: i64) -> usize {
    let last = s.len() - 1;
    let mut best = last;
    let mut best_score = s[last].score;
    let mut k = last;
    while k > 0 && s[k - 1].pos > n - 4 {
        k -= 1;
        if s[k].score > best_score {
            best = k;
            best_score = s[k].score;
        }
    }
    best
}

/// Traceback (0x4123b0, viterbi.md 7): follow back pointers from the best
/// end and emit one [`Gene`] per GE site met; genes come out right to left
/// and are returned in ascending order.
///
/// Per gene: the right end is the GE site (`pos + 4`, 1-based last codon
/// base; for pos > n-4 a pseudo stop: `pos + 1` and mark '>'), the left end
/// its predecessor (`pos - 2`; for pos < 3 a pseudo start: `pos + 1` and
/// mark '<'). RBS spacer/score, gene score and start flags come from the
/// start site (right end for '-', left end for '+').
///
/// Reverse overlap entry: if a reverse start's back pointer is another
/// gene's right end (forward stop or reverse start), the gene's left end is
/// its ORF root (reverse stop), and that root's back pointer is rewritten
/// so that the walk continues into the overlapping gene. Forward overlaps
/// were already re-routed in the DP (start.back).
///
/// `m` supplies `log_start`, RBSS (W) and NDEC for the RBS score.
pub fn traceback(s: &mut [Site], n: i64, m: &Model) -> Vec<Gene> {
    let rbs = |x: &Site| (x.st - m.log_start[(x.flags & 3) as usize]) + (x.spacer as i32 + m.rbss) as f64 / m.ndec;
    let mut genes = Vec::new();
    let mut cur = best_end(s, n);
    loop {
        // walk back until a GE site with a predecessor; none -> done
        let mut next_back = s[cur].back;
        let mut gene_end = cur;
        let gene_begin;
        loop {
            match next_back {
                None => {
                    genes.reverse();
                    return genes;
                }
                Some(x) => {
                    if s[gene_end].flags & 0x18 != 0 {
                        gene_begin = x;
                        break;
                    }
                    next_back = s[x].back;
                    gene_end = x;
                }
            }
        }
        // the binary sets this here but recomputes it below
        let _ = gene_begin;
        let mut g = Gene {
            id: 0,
            left: 0,
            right: 0,
            strand: b'+',
            cls: 0,
            lm: b' ',
            rm: b' ',
            spacer: -1,
            rbs: 0.0,
            evf: 0,
            gs: 0.0,
            sflags: 0,
            logodd: 0.0,
            maxlo: 0.0,
        };
        let end_rec = s[gene_end].clone();
        if end_rec.flags & 0x10 != 0 && end_rec.evf != 0 {
            g.evf = 1;
        }
        if end_rec.pos > n - 4 {
            g.rm = b'>';
            g.right = end_rec.pos + 1;
        } else {
            g.right = end_rec.pos + 4;
        }
        g.strand = if end_rec.flags & 0x30 != 0 { b'+' } else { b'-' };
        g.cls = end_rec.cls;
        if end_rec.flags & 8 != 0 {
            // reverse start = start codon of a '-' gene
            g.spacer = end_rec.spacer;
            g.rbs = rbs(&end_rec);
            g.gs = end_rec.gs;
            g.sflags = end_rec.flags;
        }
        let gene_begin = if g.strand == b'-' {
            let b = end_rec.back.unwrap();
            if s[b].flags & 0x10 != 0 || s[b].flags & 8 != 0 {
                // overlap entry: left end = ORF root, continue from b
                let mut root = gene_end;
                while let Some(l) = s[root].link {
                    root = l;
                }
                s[root].back = Some(b);
                root
            } else {
                b
            }
        } else {
            end_rec.back.unwrap()
        };
        let begin_rec = s[gene_begin].clone();
        if begin_rec.pos < 3 {
            g.lm = b'<';
            g.left = begin_rec.pos + 1;
        } else {
            g.left = begin_rec.pos - 2;
        }
        if begin_rec.flags & 4 != 0 && begin_rec.evf != 0 {
            g.evf = 1;
        }
        if begin_rec.flags & 0x20 != 0 {
            // forward start = start codon of a '+' gene
            g.spacer = begin_rec.spacer;
            g.rbs = rbs(&begin_rec);
            g.gs = begin_rec.gs;
            g.sflags = begin_rec.flags;
        }
        genes.push(g);
        cur = gene_begin;
    }
}

/// Partial marks for genes ending on N (0x4120f0). `orig` = the codes
/// before gap splitting, with 5/6 at resolved R/Y (so those do not count
/// as N). If the left end base is N: mark '<' and move the left end right
/// by whole codons to the first non-N codon start (unchanged if none up to
/// the right end). Then the same for the right end, moving left, bounded
/// by the (possibly moved) left end. Coordinates are 1-based.
pub fn n_ends(genes: &mut [Gene], orig: &[u8]) {
    for g in genes.iter_mut() {
        let right = g.right;
        if orig[(g.left - 1) as usize] == 4 {
            g.lm = b'<';
            let mut i = g.left;
            loop {
                i += 3;
                if right < i {
                    break;
                }
                if orig[(i - 1) as usize] != 4 {
                    g.left = i;
                    break;
                }
            }
        }
        if orig[(right - 1) as usize] == 4 {
            g.rm = b'>';
            let mut j = right;
            loop {
                j -= 3;
                if j < g.left {
                    break;
                }
                if orig[(j - 1) as usize] != 4 {
                    g.right = j;
                    break;
                }
            }
        }
    }
}

/// Log-odds score of a gene (0x429d70 '+', 0x429e10 '-'; output.md 5):
/// coding minus noncoding log likelihood of the coding region (codons
/// excluded): per base `(x + X[p][idx]) - NC[idx]`, then the joint terms
/// `(x + X_joint[idx(end)]) - nl[k][idx(end)]`.
///
/// Always uses table set B (COD2): the binary picks set A only if the
/// gene's class byte equals the integer 1, but it holds ASCII '1'/'2'.
/// '-' genes use the reverse tables in the order 6, 5, 4.
///
/// The do-while (as in the binary) only stops when `a` meets `end` exactly,
/// i.e. it needs a gene length L >= k + 8 (after [`n_ends`]). A shorter
/// gene makes the binary's loop run past the image until it terminates; it
/// needs MINC < k + 8 or a gene that [`n_ends`] shrinks to a few bases.
/// Intended difference: the port then skips the per-base loop (logodd =
/// joint term only) and returns `short = true` so that main warns.
pub fn logodd(g: &Gene, idx: &Image, m: &Model) -> (f64, bool) {
    let k = &m.markov;
    let order = m.ordm as i64;
    let (p, last): ([&Vec<f64>; 3], &Vec<f64>) = if g.strand == b'+' { ([&k.c2[0], &k.c2[1], &k.c2[2]], &k.c2[3]) } else { ([&k.c2[6], &k.c2[5], &k.c2[4]], &k.c2[7]) };
    let nc = &k.nc;
    let nck = &k.nl[order as usize];
    // a = first base after the left codon (0-based), end = start of the
    // (k+1)-mer ending at the last base before the right codon
    let mut a = g.left - 1 + 3;
    let end = g.right - 1 - 3 - order;
    let mut x = 0.0f64;
    let ix = |i: i64| idx.at(i);
    if end - a < 1 {
        return ((x + last[ix(end)]) - nck[ix(end)], true);
    }
    loop {
        x = (x + p[0][ix(a)]) - nc[ix(a)];
        if a + 1 == end {
            break;
        }
        x = (x + p[1][ix(a + 1)]) - nc[ix(a + 1)];
        if a + 2 == end {
            break;
        }
        x = (x + p[2][ix(a + 2)]) - nc[ix(a + 2)];
        a += 3;
        if a == end {
            break;
        }
    }
    ((x + last[ix(end)]) - nck[ix(end)], false)
}
