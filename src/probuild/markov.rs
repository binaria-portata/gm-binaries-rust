//! Markov chain training (models.md §1) and the whole-file sequence stream
//! reader (models.md §0.2) of probuild 2.52.
//!
//! This is the numerical core of every model-building mode: `--mkmod_euk`,
//! `--mkmod_non`, `--mkmod` (COD1 / COD2 / NONC tables) and `--mkchain`.
//!
//! # Encoding
//!
//! Bases are coded A=0, C=1, G=2, T/U=3, anything else (N, IUPAC, letters
//! of later FASTA deflines) = 4.  A word of `order + 1` letters is indexed in
//! base 4 with the FIRST letter most significant, so the table of an order-k
//! chain has `4^(k+1)` entries and entry `16*b0 + 4*b1 + b2` is the word
//! `b0 b1 b2` for k = 2.  The last letter of the word is the "predicted"
//! base, the first k letters are its context.
//!
//! # What a trained chain contains
//!
//! The chain keeps one table per phase (`period` = 3 for coding sequence,
//! 1 for non-coding).  Counting adds 1.0 for every window of `order_count +
//! 1` letters without an N; phase `j` collects the windows whose last letter
//! is at offset `order_count + j (mod period)` from the start of the counted
//! range, so with period 3 the phases are the three codon positions (relative
//! to the start of each training sequence).
//!
//! Normalisation turns each phase table into JOINT word probabilities
//! `P(w) = (count(w) + pc) / (sum of counts + 4^(order+1) * pc)` with an
//! additive pseudocount `pc` per word (`--MKCHAIN_PSEUDOCOUNTS`, default 1).
//! Conditional probabilities P(last letter | context) are not formed here;
//! the consumer (gmhmmp / gmhmme3) derives them from the joint table.
//!
//! When `order_count < order` (only possible for the non-coding chain via
//! `--order_non`), words are counted at the shorter length and then extended
//! letter by letter to `order + 1` letters (see [`Chain::normalize`]).
//!
//! # Output formats
//!
//! [`Chain::word_values`] gives one line per word (GeneMark-ES `--mkmod_euk`
//! files), [`Chain::table`] the column table of GeneMarkS-T models and
//! `--mkchain` files.  Values are printed `%#.Pf` with tiny non-zero values
//! raised to `10^-P` so that no non-zero probability prints as zero.

use super::io;
use super::pm::{fixed, powi, P};

/// Base code of a letter (table at stream object +0x30, built by
/// 0x43bd94): A/a 0, C/c 1, G/g 2, T/t/U/u 3, any other letter 4.
pub fn code(c: u8) -> u8 {
    match c {
        b'a' | b'A' => 0,
        b'c' | b'C' => 1,
        b'g' | b'G' => 2,
        b't' | b'T' | b'u' | b'U' => 3,
        _ => 4,
    }
}

/// C locale `isspace`.
fn c_isspace(c: u8) -> bool {
    matches!(c, b' ' | b'\t' | b'\n' | 0x0b | 0x0c | b'\r')
}

/// Turns a loaded file into ONE coded sequence (0x43c254 + 0x43c33e).
///
/// Leading white space is skipped; if the next char is '>' the rest of that
/// line is skipped too.  After that every ASCII letter is kept as its
/// [`code`], everything else (digits, white space, '*', '-', '>') is
/// dropped.  Only the first defline is recognised: the letters of later
/// '>' lines become part of the sequence (models.md §0.2, behaviour that the
/// GeneMark-ES training files never trigger because they have no deflines).
pub fn decode_stream(d: &[u8]) -> Vec<u8> {
    let mut i = 0;
    while i < d.len() && c_isspace(d[i]) {
        i += 1;
    }
    if i < d.len() && d[i] == b'>' {
        // stop AT the '\n'; it is dropped below like any non-letter
        while i < d.len() && d[i] != b'\n' {
            i += 1;
        }
    }
    d[i..].iter().filter(|c| c.is_ascii_alphabetic()).map(|&c| code(c)).collect()
}

/// Loads a file with the core loader 0x40a306 (which prints its own error
/// messages) and decodes it as one sequence (0x43c44a).  `None` on a load
/// failure.
pub fn read_stream(name: &[u8]) -> Option<Vec<u8>> {
    io::load_file(name).map(|d| decode_stream(&d))
}

/// The word with table index `i` as `order + 1` letters "ACGT", most
/// significant (first) letter first (0x4406bc).
pub fn word(mut i: usize, order: i32) -> String {
    let n = (order + 1) as usize;
    let mut s = vec![b'?'; n];
    for l in 0..n {
        s[n - 1 - l] = b"ACGT"[i & 3];
        i >>= 2;
    }
    String::from_utf8(s).unwrap()
}

/// Index of the reverse complement of word `j`, where `j` has `k + 1`
/// letters (0x41c27a).  The last letter of `j` becomes the first of the
/// result, and `!j & 3` = `3 - (j & 3)` complements it (A<->T, C<->G).
pub fn revcomp_index(mut j: usize, k: i32) -> usize {
    let mut c = 0;
    for _ in 0..=k {
        c = c * 4 + (!j & 3);
        j >>= 2;
    }
    c
}

/// Strand symmetrisation of a table of `k + 1`-letter words (0x41c2bc):
/// every word and its reverse complement both get their mean.
///
/// Used for non-coding chains with `--revcomp_non`, where both strands
/// should look alike.  The loop visits every pair twice (at `j` and at its
/// partner); the second visit averages two equal values and changes
/// nothing.  The sum is formed as `v[r] + v[j]` in this order, as in the
/// original.
fn symmetrise(table: &mut [f64], k: i32) {
    for j in 0..(4usize << (2 * k)) {
        let r = revcomp_index(j, k);
        let d = (table[r] + table[j]) / 2.0;
        table[j] = d;
        table[r] = d;
    }
}

/// Markov chain object (0x41b982 ...; models.md §1, layout in the comments
/// as offsets of the original object).
pub struct Chain {
    /// +0: bases cut from the left end of every training sequence before
    /// counting (`--MKCHAIN_L_MARGING`, default 3; only used by set
    /// training, [`Chain::train_set`]).
    pub l_margin: i32,
    /// +4: bases cut from the right end (`--MKCHAIN_R_MARGING`, default 3).
    pub r_margin: i32,
    /// +8: also count the reverse complement strand (`--revcomp_non`).  Only
    /// honoured for period 1 (non-coding chains).
    pub revcomp: bool,
    /// +0xc: Markov order k; words have k + 1 letters.
    pub order: i32,
    /// +0x10: number of phase tables, 3 for coding (codon positions) or 1.
    pub period: usize,
    /// +0x18: one table of `4^(order+1)` doubles per phase.  Holds counts
    /// during training and joint word probabilities after
    /// [`Chain::normalize`].
    pub tables: Vec<Vec<f64>>,
    /// +0x30: word length - 1 actually counted.  Equals `order` except for
    /// a period-1 chain with `--order_non` smaller than the order.
    pub order_count: i32,
    /// +0x34: additive pseudocount per word (`--MKCHAIN_PSEUDOCOUNTS`,
    /// default 1).
    pub pseudocounts: i32,
    /// +0x38: decimals of the printed values (`--MKCHAIN_PRECISION`,
    /// default 5; must be 0..10 when printing).
    pub precision: i32,
}

/// Prints `m` to stderr and terminates with `exit(1)` (the original calls
/// `exit` directly from the chain code, bypassing the handler's return).
fn die(m: String) -> ! {
    io::err(m.as_bytes());
    io::c_exit(1)
}

impl Chain {
    /// Builds an empty chain of `order` and `period` from the parameters
    /// (0x41b982, checks 0x41b404, getters 0x41b66e).
    ///
    /// The options are re-read for every chain built, so an invalid value
    /// prints its error once per chain; getter failures only print and keep
    /// the default.  Range errors terminate the program with `exit(1)`.
    pub fn from_params(p: &mut P, order: i32, period: usize) -> Chain {
        // --order_non: -1 (the default, also when given explicitly) means
        // "same as order"; it is ignored for coding chains (period 3).
        let (_, order_non) = p.int_opt("--order_non", -1, true);
        let mut order_count = if order_non == -1 { order } else { order_non };
        if period != 1 {
            order_count = order;
        }
        if !(0..=14).contains(&order) {
            die(format!("error: order {order}\n"));
        }
        if !(0..=14).contains(&order_count) {
            die(format!("error: order {order_count}\n"));
        }
        if order < order_count {
            die(format!("error: order_count {order_count} is higher than order {order}\n"));
        }
        // 4^(order+1) words
        let size = 4usize << (2 * order);
        // Orders 13 and 14 need gigabytes per table. Where the memory is not
        // available the original's `new` throws std::bad_alloc ("Out of
        // memory!", exit 1); Rust would abort instead, so allocate fallibly.
        let mut tables = Vec::with_capacity(period);
        for _ in 0..period {
            let mut t: Vec<f64> = Vec::new();
            if t.try_reserve_exact(size).is_err() {
                io::throw("std::bad_alloc", "Markov chain table");
            }
            t.resize(size, 0.0);
            tables.push(t);
        }
        let mut c = Chain {
            l_margin: 3,
            r_margin: 3,
            revcomp: false,
            order,
            period,
            tables,
            order_count,
            pseudocounts: 1,
            precision: 5,
        };
        c.l_margin = p.int_opt("--MKCHAIN_L_MARGING", c.l_margin, false).1;
        c.r_margin = p.int_opt("--MKCHAIN_R_MARGING", c.r_margin, false).1;
        c.pseudocounts = p.int_opt("--MKCHAIN_PSEUDOCOUNTS", c.pseudocounts, false).1;
        c.precision = p.int_opt("--MKCHAIN_PRECISION", c.precision, false).1;
        c.revcomp = p.bool_opt("--revcomp_non", false).1;
        c
    }

    /// Number of words per phase table, `4^(order+1)`.
    pub fn size(&self) -> usize {
        self.tables[0].len()
    }

    /// Adds the words of `seq` to the counts (0x41bab2), or of its reverse
    /// complement if `rev` (0x41bce2, which walks the range backwards and
    /// complements each code as `3 - c`; N stays N).
    ///
    /// Only words of `order_count + 1` letters are counted, so with
    /// `order_count < order` only the first `4^(order_count+1)` entries of
    /// the table are used.  The first `order_count` letters prime the word
    /// index; from then on every position closes one window.  A window that
    /// contains an N is skipped but still advances the phase.  The phase
    /// starts at 0 for every call, i.e. at the first complete window.
    fn count(&mut self, seq: &[u8], rev: bool) {
        let oc = self.order_count as usize;
        if seq.len() < oc {
            return;
        }
        let n = seq.len();
        let at = |i: usize| if rev { seq[n - 1 - i] } else { seq[i] };
        // word_idx: base-4 index of the current window; n_mask: 3 in every
        // 2-bit slot that holds an N; keep_mask: 4^oc - 1, the slots of the
        // last oc letters kept when the window slides.
        let (mut word_idx, mut n_mask, mut keep_mask) = (0u32, 0u32, 0u32);
        for i in 0..oc {
            let c = at(i);
            word_idx <<= 2;
            n_mask <<= 2;
            keep_mask = keep_mask.wrapping_mul(4).wrapping_add(3);
            if c == 4 {
                n_mask += 3;
            } else {
                word_idx += if rev { 3 - c as u32 } else { c as u32 };
            }
        }
        let mut phase = 0;
        for i in oc..n {
            let c = at(i);
            word_idx = (word_idx & keep_mask) << 2;
            n_mask = (n_mask & keep_mask) << 2;
            if c == 4 {
                n_mask += 3;
            } else {
                word_idx += if rev { 3 - c as u32 } else { c as u32 };
            }
            if n_mask == 0 {
                self.tables[phase][word_idx as usize] += 1.0;
            }
            phase += 1;
            if phase == self.period {
                phase = 0;
            }
        }
    }

    /// Trains on one whole sequence stream (0x41d18e; `--mkmod_euk`,
    /// `--mkmod_non`, `--mkmod --cod/--non`).  Margins are ignored, the
    /// reverse strand is counted if `revcomp` and period 1, and the counts
    /// are NOT symmetrised (unlike [`Chain::train_set`]).  Then normalises.
    pub fn train_stream(&mut self, seq: &[u8]) {
        self.count(seq, false);
        if self.revcomp && self.period == 1 {
            self.count(seq, true);
        }
        self.normalize();
    }

    /// Trains on a set of sequences (0x41cfd2 -> 0x41bf06, then 0x41c3a2;
    /// `--mkmod` / `--mkchain` with `--geneset`).
    ///
    /// Each sequence longer than `lm + rm` is counted without its first `lm`
    /// and last `rm` bases (the phase restarts at 0 for every sequence, so
    /// for genes phase 0 is codon position `order + lm` mod 3).  With
    /// `revcomp` and period 1 the reverse complement of `s[rm .. len-lm)` is
    /// counted too (the window is mirrored, which only matters for lm != rm;
    /// the non-coding chain is always trained with 0, 0) and after all
    /// sequences the counts at level `order_count` are strand-symmetrised.
    /// Negative margins terminate with `exit(1)`.
    pub fn train_set(&mut self, seqs: &[Vec<u8>], lm: i32, rm: i32) {
        if lm < 0 || rm < 0 {
            die(format!("error: Lmarg {lm} Rmarg {rm}\n"));
        }
        let both_strands = self.revcomp && self.period == 1;
        let (l, r) = (lm as usize, rm as usize);
        for s in seqs {
            if l + r < s.len() {
                self.count(&s[l..s.len() - r], false);
                if both_strands {
                    self.count(&s[r..s.len() - l], true);
                }
            }
        }
        if both_strands {
            let k = self.order_count;
            symmetrise(&mut self.tables[0], k);
        }
        self.normalize();
    }

    /// Fills every table with a zero-order "random sequence" model of GC
    /// content `gc` (0x41d00c, `--RANDOM_NONCOD`): the joint probability of
    /// a word is the product of its letter probabilities, A = T = (1-gc)/2,
    /// C = G = gc/2.
    ///
    /// The product is formed from the LAST letter to the first, starting
    /// from 1.0, as `p * x` (operand order kept).  No range check: a GC
    /// percentage instead of a fraction (the --cod/--non path of --mkmod)
    /// gives values above 1 and negative ones.
    pub fn random_gc(&mut self, gc: f64) {
        let p = [(1.0 - gc) / 2.0, gc / 2.0, gc / 2.0, (1.0 - gc) / 2.0];
        let order = self.order;
        for table in self.tables.iter_mut() {
            for (word_idx, e) in table.iter_mut().enumerate() {
                let mut x = 1.0;
                let mut rest = word_idx;
                for _ in 0..=order {
                    x = p[rest & 3] * x;
                    rest >>= 2;
                }
                *e = x;
            }
        }
    }

    /// Converts the counts into joint word probabilities (0x41c3a2).
    ///
    /// `order_count == order` (0x41cdb8): per phase table, `t` = sum of all
    /// counts (ascending index), then `P(w) = (pc + n(w)) / (N*pc + t)` with
    /// `N = 4^(order+1)`.  `N*pc` is an int product converted to double
    /// afterwards (wrapping like the original's 32-bit multiply).
    ///
    /// `order_count < order` (0x41c3ea, then 0x41cb9c): only one table is
    /// allowed (period 1).  The `4^(order_count+1)` counted entries are
    /// normalised as above; then the word length is extended one letter at a
    /// time up to `order + 1` (see the comment in the loop) and the table is
    /// strand-symmetrised after every extension step.  (The original also
    /// checks `size <= n0` here with `error in vector size`; that cannot
    /// happen because order_count < order, so it is not reproduced.)
    ///
    /// A negative pseudocount terminates with `exit(1)`.
    pub fn normalize(&mut self) {
        if self.order_count == self.order {
            // 0x41cdb8
            if self.pseudocounts < 0 {
                die(format!("error: pseudocount {}\n", self.pseudocounts));
            }
            let pc = self.pseudocounts;
            let pc_total = (self.size() as i32).wrapping_mul(pc) as f64;
            for table in self.tables.iter_mut() {
                let mut t = 0.0;
                for &x in table.iter() {
                    t += x;
                }
                for x in table.iter_mut() {
                    *x = (pc as f64 + *x) / (pc_total + t);
                }
            }
        } else {
            // 0x41c3ea
            if self.period != 1 && self.tables.len() != 1 {
                die(format!("error: order_count is supported only for period == 1 {}\n", self.order_count));
            }
            let pc = self.pseudocounts;
            if pc < 0 {
                die(format!("error: pseudocount {pc}\n"));
            }
            // 0x41cb9c: normalise the counted (order_count+1)-letter words
            let n0 = 4usize << (2 * self.order_count);
            let table = &mut self.tables[0];
            let mut t = 0.0;
            for &x in table[..n0].iter() {
                t += x;
            }
            let den = (n0 as i32).wrapping_mul(pc) as f64 + t;
            for x in table[..n0].iter_mut() {
                *x = (pc as f64 + *x) / den;
            }
            // Extend from k+1 to k+2 letters.  `prev` holds P over k+1
            // letters; the four entries at block_start.. are the words
            // `w x` for one k-letter prefix w and x = A,C,G,T, so
            // q[i] = P(w i) / P(w .) = P(i | w).  The new entry for the word
            // `w a i` (index 4*(block_start + a) + i, written in ascending
            // order) is P(w a) * P(i | w).  As in the original: the context of the
            // new letter i is the prefix w, not the k letters directly in
            // front of i (models.md §16).  A prefix with total 0 gives 0.
            for k in self.order_count..self.order {
                let prev = table.clone();
                let mut out = 0;
                let mut block_start = 0;
                while block_start < (4usize << (2 * k)) {
                    let b = block_start;
                    // ((p0 + p1) + p2) + p3, left to right as in the original
                    let s = prev[b] + prev[b + 1] + prev[b + 2] + prev[b + 3];
                    let q: [f64; 4] = if s == 0.0 {
                        [0.0; 4]
                    } else {
                        [prev[b] / s, prev[b + 1] / s, prev[b + 2] / s, prev[b + 3] / s]
                    };
                    for a in 0..4 {
                        for qi in q.iter() {
                            table[out] = prev[b + a] * qi;
                            out += 1;
                        }
                    }
                    block_start += 4;
                }
                symmetrise(table, k + 1);
            }
        }
    }

    /// Checks the print precision and returns the smallest printable value
    /// `1 / 10^P` (with libgcc `__powidf2`, as the original).  A precision
    /// outside 0..10 prints `error: precision P` and exits with 1 (files
    /// already written stay as they are).
    fn check_precision(&self) -> f64 {
        let p = self.precision;
        if !(0..=10).contains(&p) {
            io::err(format!("error: precision {p}\n").as_bytes());
            io::c_exit(1);
        }
        1.0 / powi(10.0, p)
    }

    /// One value as `%#.Pf`; a non-zero value below `eps` is printed as
    /// `eps` so that it does not appear as 0.
    fn format_value(&self, eps: f64, x: f64) -> String {
        let x = if x < eps && x != 0.0 { eps } else { x };
        fixed(self.precision, x, true)
    }

    /// The values of word `idx` on one line (0x41d708): period 1 `v`,
    /// period 3 `v0 v1 v2` (single spaces).  Used for the one-line-per-word
    /// files of `--mkmod_euk` / `--mkmod_non`.
    pub fn word_values(&self, idx: usize) -> String {
        let eps = self.check_precision();
        match self.tables.len() {
            1 => self.format_value(eps, self.tables[0][idx]),
            3 => format!(
                "{} {} {}",
                self.format_value(eps, self.tables[0][idx]),
                self.format_value(eps, self.tables[1][idx]),
                self.format_value(eps, self.tables[2][idx])
            ),
            // unreachable: the constructor only allows period 1 or 3
            _ => String::new(),
        }
    }

    /// The whole chain as model-file table (0x41d6d8 -> 0x41d240): `"\n"`,
    /// then one line per word in index order, `v` (period 1) or
    /// `v0 v1 v2` (period 3).  Used for COD1/COD2/NONC of `--mkmod` and the
    /// `--mkchain` files.
    pub fn table(&self) -> String {
        let eps = self.check_precision();
        let mut o = String::from("\n");
        if self.tables.len() == 1 {
            for &x in &self.tables[0] {
                o += &self.format_value(eps, x);
                o.push('\n');
            }
        }
        if self.tables.len() == 3 {
            for i in 0..self.size() {
                o += &format!(
                    "{} {} {}\n",
                    self.format_value(eps, self.tables[0][i]),
                    self.format_value(eps, self.tables[1][i]),
                    self.format_value(eps, self.tables[2][i])
                );
            }
        }
        o
    }
}
