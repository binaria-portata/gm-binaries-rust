//! Training-set extraction from a genome and a gene list: the sequence
//! object (0x412f7e / 0x43c64c) and the gene parser object (ctor 0x405c0e).
//! models.md §0.2, §0.5.
//!
//! Given the `--seq` sequence (one coded stream, see
//! [`super::markov::decode_stream`]) and a `--geneset` list, the gene
//! parser cuts out the material the model modes train on:
//!
//! - coding sequences of the valid genes (set 1, and set 2 for class-2
//!   genes with `--SEPARATE_COD`), oriented 5'->3', stop codon included,
//!   [`GeneParser::collect_coding`];
//! - non-coding segments: the stream with every listed region masked,
//!   split into runs of at least `--MIN_NON` bases,
//!   [`GeneParser::collect_noncoding`];
//! - start codon frequencies of the coding set (ATG/GTG/TTG/CTG entries of
//!   a GeneMarkS-T model), [`GeneParser::start_stats`];
//! - pre-start sequences around the start codons and their positional
//!   nucleotide matrix (the RBS / upstream motif of `--mkmod --fixmotif`
//!   and the `--pre_start` file fed to Gibbs3),
//!   [`GeneParser::collect_pre_start`], [`GeneParser::pre_start_matrix`].
//!
//! Users: `--mkmod`, `--parse`, `--mkchain`, `--intergenic`, `--pre_start`,
//! `--STAT_NON_VS_ALL`, `--ORF_PER_SEQ`.
//!
//! Codon index of bases b0 b1 b2 (A0 C1 G2 T3) = 16*b0 + 4*b1 + b2, e.g.
//! ATG 14, GTG 46, TTG 62, CTG 30, TAA 48, TAG 50, TGA 56.

use super::geneset::Gene;
use super::io;
use super::markov::decode_stream;
use super::pm::{fixed, g, P};

/// Decoding table of the gene parser (0x40c52c): code 0..5 -> letter.  Code
/// 4 is N; 5 (the masking value of [`GeneParser::collect_noncoding`]) is
/// never printed.
pub const LETTERS: &[u8] = b"ACGTNX";

/// `std::cout << s << std::endl`: the text, a newline, and a flush.
pub fn out_endl(s: &str) {
    io::out(s.as_bytes());
    io::out(b"\n");
    io::flush_stdout();
}

/// Reverse complement of a coded sequence (0x43c756): reversed, 0<->3,
/// 1<->2, N (4) stays 4.
pub fn revcomp(x: &[u8]) -> Vec<u8> {
    x.iter().rev().map(|&c| if c < 4 { 3 - c } else { 4 }).collect()
}

/// Whole-file sequence object (0x412f7e): the `--seq` stream plus the stop
/// codon table used to reject coding sequences with in-frame stops.
pub struct SeqObj {
    /// The coded stream (codes 0..4), see [`decode_stream`].
    pub seq: Vec<u8>,
    /// +0x48: stop codons by codon index, from `--TAA_` (48), `--TGA_` (56),
    /// `--TAG_` (50); all false by default.
    pub stop_codon: [bool; 64],
}

impl SeqObj {
    pub fn new() -> SeqObj {
        SeqObj { seq: Vec::new(), stop_codon: [false; 64] }
    }

    /// Reads the stop codon options (0x43c006; failures only print), then
    /// the required `--seq` and the file (0x43c44a).  `false` if `--seq` is
    /// missing/invalid or the file cannot be loaded (0x43c64c).
    pub fn load(&mut self, p: &mut P) -> bool {
        for (k, c) in [("--TAA_", 48), ("--TGA_", 56), ("--TAG_", 50)] {
            self.stop_codon[c] = p.bool_opt(k, false).1;
        }
        let name = match p.str_req("--seq") {
            Some(n) => n,
            None => return false,
        };
        match io::load_file(&name) {
            Some(d) => {
                self.seq = decode_stream(&d);
                true
            }
            None => false,
        }
    }

    /// Scans `x` in frame 0 for the first stop codon (0x43c814).  Returns
    /// the 1-based index of the LAST base of that codon, or `len(x)` if
    /// there is none (also for `len(x) < 3`).  Codons containing N are not
    /// stops.
    pub fn stop_index(&self, x: &[u8]) -> usize {
        let n = x.len();
        if n < 3 {
            return n;
        }
        // pos_in_codon counts 0,1,2 within the codon; scanned is the
        // 1-based position reached when the loop ends
        let mut pos_in_codon = 0;
        let mut scanned = 0;
        for i in 0..n {
            scanned += 1;
            if pos_in_codon == 2 {
                pos_in_codon = 0;
                let (a, b, c) = (x[i - 2], x[i - 1], x[i]);
                if a < 4 && b < 4 && c < 4 && self.stop_codon[((a << 4) | (b << 2) | c) as usize] {
                    break;
                }
            } else {
                pos_in_codon += 1;
            }
        }
        scanned
    }
}

/// Gene parser object (ctor 0x405c0e; field offsets in the docs).  Holds
/// the options of the training-set extraction and the extracted sets.
pub struct GeneParser {
    /// +0x98 `--MIN_COD`: minimum coding length (after margins) for a gene
    /// to enter the coding set, and minimum gene length for a pre-start
    /// sequence.  Default 0 (GeneMarkS-T's par file: 450).
    pub min_cod: i32,
    /// +0x9c `--MIN_NON`: minimum length of a non-coding segment.
    pub min_non: i32,
    /// +0xa0 `--EXCLUDE_PARTIAL`: drop genes with a '<' or '>' end from the
    /// coding set (default true).
    pub exclude_partial: bool,
    /// +0xe8 `--START_PRECISION`: decimals of the start probabilities
    /// (default 3).
    pub start_precision: i32,
    /// +0xec `--START_PSEUDOCOUNTS`: pseudocount added to each enabled start
    /// codon (default 0).
    pub start_pseudocounts: i32,
    /// +0x78 `--SEPARATE_COD`: put class-1 and class-2 genes into separate
    /// coding sets (COD1 / COD2).  Default false: all valid genes in set 1.
    pub separate_cod: bool,
    /// +0xc0: start codons that are counted, by codon index (`--ATG_` 14,
    /// `--GTG_` 46, `--TTG_` 62, `--CTG_` 30).
    pub start_enabled: [bool; 64],
    /// +0: coding set 1 (codes, 5'->3' on the gene's strand).
    pub cod1: Vec<Vec<u8>>,
    /// +0x18: coding set 2 (class-2 genes with `--SEPARATE_COD`).
    pub cod2: Vec<Vec<u8>>,
    /// +0x30: non-coding segments in sequence order (codes 0..4).
    pub non: Vec<Vec<u8>>,
    /// +0x48: pre-start sequences (upstream + start region, 5'->3').
    pub pre_start: Vec<Vec<u8>>,
    /// +0x80: positional matrix of the pre-start sequences, one row per
    /// position, columns A C G T (N in column 4 is counted but ignored by
    /// the normalisation).  Counts, then probabilities.
    pub rbs_matrix: Vec<[f64; 5]>,
    /// +0x110: start codon probabilities by codon index (0.0 for disabled
    /// codons).
    pub start_prob: [f64; 64],
    /// The port's warning about genes outside the sequence was printed
    /// (not in the original).
    pub warned_outside: bool,
}

impl GeneParser {
    /// Defaults 0x40c5f2, then the getters of 0x40c666 in the original
    /// order (failures only print and keep the default).  Called before
    /// anything else in the handlers, so its messages come first.
    pub fn new(p: &mut P) -> GeneParser {
        let min_cod = p.int_opt("--MIN_COD", 0, false).1;
        let min_non = p.int_opt("--MIN_NON", 0, false).1;
        let exclude_partial = p.bool_opt("--EXCLUDE_PARTIAL", true).1;
        let start_precision = p.int_opt("--START_PRECISION", 3, false).1;
        let start_pseudocounts = p.int_opt("--START_PSEUDOCOUNTS", 0, false).1;
        let separate_cod = p.bool_opt("--SEPARATE_COD", false).1;
        let mut start_enabled = [false; 64];
        for (k, c) in [("--ATG_", 14), ("--GTG_", 46), ("--TTG_", 62), ("--CTG_", 30)] {
            start_enabled[c] = p.bool_opt(k, false).1;
        }
        GeneParser {
            min_cod,
            min_non,
            exclude_partial,
            start_precision,
            start_pseudocounts,
            separate_cod,
            start_enabled,
            cod1: Vec::new(),
            cod2: Vec::new(),
            non: Vec::new(),
            pre_start: Vec::new(),
            rbs_matrix: Vec::new(),
            start_prob: [0.0; 64],
            warned_outside: false,
        }
    }

    /// Builds the coding sets (0x40ce72).  `lm` / `rm` are bases removed
    /// at the 5' / 3' end of the gene (0, 0 for --mkmod and --mkchain;
    /// `--COD_L_MARGIN` / `--COD_R_MARGIN` for --parse).
    ///
    /// For every list entry in order, the gene is skipped if it is not
    /// valid, shorter than MIN_COD after the margins, partial (with
    /// EXCLUDE_PARTIAL), or shorter than the margins.  Otherwise its
    /// sequence (reverse complement for '-') is taken.  If an in-frame stop
    /// codon ends more than one base before the end (`k < len - 1`, k from
    /// [`SeqObj::stop_index`]; the regular stop codon at the end does not
    /// count), `warning, STOP codon found inside protein: a b k` is printed
    /// to stdout (a, b: 0-based cut ends; no flush) and the gene is
    /// dropped.  The final check `len == n` also drops any sequence that
    /// the out-of-range clamping shortened.
    pub fn collect_coding(&mut self, so: &SeqObj, genes: &[Gene], lm: i32, rm: i32) {
        for g in genes {
            if !g.valid {
                continue;
            }
            // coding length after the margins
            let n = g.right.wrapping_sub(g.left).wrapping_sub(rm).wrapping_sub(lm).wrapping_add(1);
            if n < self.min_cod {
                continue;
            }
            if self.exclude_partial && !(g.left_complete && g.right_complete) {
                continue;
            }
            // int arithmetic of the original (wraps for extreme margins)
            if g.right - g.left < rm.wrapping_add(lm) {
                continue;
            }
            // lm is the 5' margin: on '-' genes the 5' end is `right`
            let (a, b) = if g.strand == 1 {
                (g.left.wrapping_add(lm), g.right.wrapping_sub(rm))
            } else {
                (g.left.wrapping_add(rm), g.right.wrapping_sub(lm))
            };
            // 0x43ca9e copies the vector range [begin + a, begin + b + 1)
            // (64-bit pointer arithmetic). A negative length (only reachable
            // when the margins wrap) makes the allocation throw
            // std::bad_alloc: "Out of memory!", exit code 1.
            if (b as i64 + 1) < a as i64 {
                io::throw("std::bad_alloc", "vector");
            }
            // A range reaching outside the sequence (negative COD_R_MARGIN,
            // or a gene ending one base past the sequence, which the gene
            // list check lets through) makes the original read the heap
            // around its sequence buffer (undefined memory; far outside it
            // terminates). The port leaves such a gene out of the training set
            // (the length test below drops it) and warns once.
            let (ua, ub) = (a.max(0) as usize, ((b as i64 + 1).max(0) as usize).min(so.seq.len()));
            if a < 0 || b as i64 >= so.seq.len() as i64 {
                if !self.warned_outside {
                    self.warned_outside = true;
                    io::err(
                        b"probuild-rs: warning: a gene of --geneset reaches outside --seq after the margins \
                          (--COD_L_MARGIN / --COD_R_MARGIN) are applied; such genes are left out of the coding \
                          set (not handled by the original)\n",
                    );
                }
            }
            let mut x = if ua <= ub { so.seq[ua..ub].to_vec() } else { Vec::new() };
            if g.strand == 2 {
                x = revcomp(&x);
            }
            let k = so.stop_index(&x);
            // unsigned compare as in the original: for an empty x,
            // len - 1 wraps and the warning is printed
            if (k as u64) < (x.len() as u64).wrapping_sub(1) {
                io::out(format!("warning, STOP codon found inside protein: {a} {b} {k}\n").as_bytes());
                x.clear();
            }
            if x.len() as i64 == n as i64 {
                if !self.separate_cod || g.class == 1 {
                    self.cod1.push(x);
                } else if g.class == 2 {
                    self.cod2.push(x);
                }
                // with SEPARATE_COD another class is dropped silently
                // (unreachable: valid implies class 1 or 2)
            }
        }
    }

    /// Builds the non-coding segments (0x40d8c0 masking + 0x40d758 split).
    ///
    /// Every entry with both positions set (valid or not, so also genes
    /// rejected from the coding set and partial genes) masks
    /// `[left, right]` with code 5; an entry outside the sequence prints
    /// `error, if this is printed - error in code` to stderr instead.  Every
    /// maximal unmasked run of length >= MIN_NON becomes a segment.
    pub fn collect_noncoding(&mut self, so: &SeqObj, genes: &[Gene]) {
        let mut masked = so.seq.clone();
        let l = masked.len() as i64;
        for g in genes {
            if g.left == -1 || g.right == -1 {
                continue;
            }
            if g.left < 0 || g.right as i64 + 1 > l {
                io::err(b"error, if this is printed - error in code\n");
                continue;
            }
            if g.left <= g.right {
                for c in &mut masked[g.left as usize..=g.right as usize] {
                    *c = 5;
                }
            }
        }
        let n = masked.len();
        let mut i = 0;
        while i < n && masked[i] == 5 {
            i += 1;
        }
        while i < n {
            let mut j = i;
            while j < n && masked[j] != 5 {
                j += 1;
            }
            if self.min_non as i64 <= (j - i) as i64 {
                self.non.push(masked[i..j].to_vec());
            }
            i = j;
            while i < n && masked[i] == 5 {
                i += 1;
            }
        }
    }

    /// Collects pre-start sequences (0x40dd56, extraction 0x43cc40):
    /// `upstream` bases before the start codon plus `downstream` bases from
    /// the first base of the start codon on, oriented 5'->3'.  Appends them
    /// to `pre_start` and returns their names `" L R s\n"` (1-based gene
    /// ends, strand char).
    ///
    /// A gene is used if it is valid, at least MIN_COD long, has a complete
    /// 5' end and leaves room for the upstream part:
    /// '+': `left >= upstream`; '-': `right <= maxpos - upstream`
    /// (`maxpos` = sequence length, so this allows one base more than fits:
    /// such a sequence comes out one base short).
    ///
    /// With `intergenic`, the upstream gap to the neighbouring list ENTRY
    /// (previous entry for '+', next for '-'; whatever it is, a header
    /// line has right = -1 / left = -1) must be at least `upstream`; the
    /// first '+' entry and the last '-' entry of the list are always
    /// skipped.  Sequences containing N are skipped.
    pub fn collect_pre_start(
        &mut self,
        so: &SeqObj,
        genes: &[Gene],
        maxpos: i32,
        upstream: i32,
        downstream: i32,
        intergenic: bool,
    ) -> Vec<String> {
        let mut names = Vec::new();
        let s = &so.seq;
        let l = s.len() as i64;
        let n = genes.len();
        for (gi, g) in genes.iter().enumerate() {
            if !g.valid || g.right - g.left + 1 < self.min_cod {
                continue;
            }
            if g.strand == 1 && !(g.left_complete && upstream <= g.left) {
                continue;
            }
            if g.strand == 2 && !(g.right_complete && g.right <= maxpos.wrapping_sub(upstream)) {
                continue;
            }
            if intergenic {
                // distance to the neighbour; the gap is distance - 1 bases
                let dist = if g.strand == 1 {
                    if gi == 0 {
                        continue;
                    }
                    g.left - genes[gi - 1].right
                } else {
                    if gi + 1 >= n {
                        continue;
                    }
                    genes[gi + 1].left - g.right
                };
                if upstream > dist - 1 {
                    continue;
                }
            }
            // half-open window [a, b) of the stream, clamped to [0, l]
            let (frm, to) = (g.left as i64, g.right as i64 + 1);
            let (up, down) = (upstream as i64, downstream as i64);
            let (a, b) = if g.strand == 1 {
                (if frm - up > 0 { frm - up } else { 0 }, if frm + down < l { frm + down } else { l })
            } else {
                (if to - down > 0 { to - down } else { 0 }, if to + up < l { to + up } else { l })
            };
            let mut x = if a <= b { s[a as usize..b as usize].to_vec() } else { Vec::new() };
            if g.strand == 2 {
                x = revcomp(&x);
            }
            if x.contains(&4) {
                continue;
            }
            names.push(format!(" {} {} {}\n", g.left + 1, g.right + 1, b".+-"[g.strand as usize] as char));
            self.pre_start.push(x);
        }
        names
    }

    /// Start codon frequencies of the coding sets (0x40efc8).
    ///
    /// Counts the first codon of every set-1 and set-2 sequence (no N, at
    /// least 3 bases).  With `tot` = sum of the counts of the enabled
    /// codons and `k` = number of enabled codons:
    /// `P(c) = (pc + count(c)) / (tot + pc * k)` for enabled codons (int
    /// arithmetic, then one double division), disabled ones stay 0.0.  If
    /// no enabled start was seen, `warning, no start sites in training\n`
    /// plus endl (so an empty line follows) goes to stdout; if the
    /// denominator is 0 all probabilities stay 0.0.
    pub fn start_stats(&mut self) {
        let mut cnt = [0i32; 64];
        for x in self.cod1.iter().chain(self.cod2.iter()) {
            if x.len() > 2 && x[0] < 4 && x[1] < 4 && x[2] < 4 {
                cnt[((x[0] << 4) | (x[1] << 2) | x[2]) as usize] += 1;
            }
        }
        let mut tot: i32 = 0;
        let mut k: i32 = 0;
        for c in 0..64 {
            if self.start_enabled[c] {
                tot = tot.wrapping_add(cnt[c]);
                k += 1;
            }
        }
        if tot == 0 {
            out_endl("warning, no start sites in training\n");
        }
        tot = tot.wrapping_add(self.start_pseudocounts.wrapping_mul(k));
        if tot != 0 {
            for c in 0..64 {
                if self.start_enabled[c] {
                    self.start_prob[c] = self.start_pseudocounts.wrapping_add(cnt[c]) as f64 / tot as f64;
                }
            }
        }
    }

    /// Probability of start codon `c` as `%#.Pf`, P = START_PRECISION
    /// (0x40f3ee).
    pub fn start_prob_str(&self, c: usize) -> String {
        fixed(self.start_precision, self.start_prob[c], true)
    }

    /// GC fraction of the non-coding segments (0x411506), the default GC of
    /// `--RANDOM_NONCOD` on the gene-set path of --mkmod.  Double counts;
    /// `(C + G) / (T + A + C + G)` with the sums in exactly this order; 0.0
    /// without bases.  Prints `GC <%g>` + endl to stdout.
    pub fn gc_non(&self) -> f64 {
        let mut cnt = [0.0f64; 6];
        for x in &self.non {
            for &c in x {
                cnt[c as usize] += 1.0;
            }
        }
        let mut gc = 0.0;
        if cnt[3] + cnt[0] + cnt[1] + cnt[2] != 0.0 {
            gc = (cnt[1] + cnt[2]) / (cnt[3] + cnt[0] + cnt[1] + cnt[2]);
        }
        out_endl(&format!("GC {}", g(gc)));
        gc
    }

    /// Positional frequency matrix of the pre-start sequences (0x41184c),
    /// the RBS model of `--mkmod --fixmotif`.  `false` if a position has a
    /// zero total (only possible with pc <= 0; the caller fails silently).
    ///
    /// The row count is the length of the FIRST sequence; with no
    /// sequences the original dereferences NULL and terminates, the port stops
    /// with an error message.  Per
    /// position: `t = 4*pc + G + C + A + T` (int `pc << 2` first, then the
    /// counts in this order, N excluded), `P(b) = (pc + count(b)) / t`.
    pub fn pre_start_matrix(&mut self, pc: i32) -> bool {
        if self.pre_start.is_empty() {
            // The original takes the size of *begin() of the empty vector
            // (NULL dereference, SIGSEGV).
            io::fatal(
                "no pre-start (upstream) sequences to build the start motif model from; \
                 the gene set is empty or no gene has enough sequence upstream of its start \
                 (check --geneset and --PRE_START_WIDTH, or run without --fixmotif)",
            );
        }
        let n = self.pre_start[0].len();
        self.rbs_matrix = vec![[0.0; 5]; n];
        for x in &self.pre_start {
            for i in 0..n {
                // A later sequence shorter than the first is read past its
                // end in the original (undefined read); here the missing bases count
                // as N (4), which the normalisation ignores.
                let c = x.get(i).copied().unwrap_or(4);
                self.rbs_matrix[i][c as usize] += 1.0;
            }
        }
        for row in self.rbs_matrix.iter_mut() {
            let t = (pc.wrapping_shl(2)) as f64 + row[2] + row[1] + row[0] + row[3];
            if t == 0.0 {
                return false;
            }
            for j in 0..4 {
                row[j] = (pc as f64 + row[j]) / t;
            }
        }
        true
    }

    /// The matrix as text (0x411d2c), values `%#.Pf` joined by spaces.
    ///
    /// Mode 0 (the --log file): header `Pos   A   C   G   T`, then rows
    /// `i+1 A C G T` for the upstream positions 0..w-1, an empty line, the
    /// three start-codon rows w..w+2, an empty line, the rest.
    ///
    /// Mode 1 (RBSM of the model): rows only, from the last position down
    /// to w+3, then from w-1 down to 0, skipping the three start-codon
    /// rows (models.md §4.2 step 10).
    pub fn matrix_text(&self, mode: i32, w: i32, prec: i32) -> String {
        let row_text = |row: &[f64; 5]| -> String {
            (0..4).map(|j| fixed(prec, row[j], true)).collect::<Vec<_>>().join(" ")
        };
        let n = self.rbs_matrix.len() as i64;
        let w = w as i64;
        let mut o = String::new();
        if mode == 0 {
            o += "Pos   A   C   G   T\n";
            let mut i: i64 = 0;
            while i < n && i != w {
                o += &format!("{} {}\n", i + 1, row_text(&self.rbs_matrix[i as usize]));
                i += 1;
            }
            o += "\n";
            while i < n && i != w + 3 {
                o += &format!("{} {}\n", i + 1, row_text(&self.rbs_matrix[i as usize]));
                i += 1;
            }
            o += "\n";
            while i < n {
                o += &format!("{} {}\n", i + 1, row_text(&self.rbs_matrix[i as usize]));
                i += 1;
            }
        } else {
            let mut i = n - 1;
            while i != w + 2 && i >= 0 {
                o += &row_text(&self.rbs_matrix[i as usize]);
                o += "\n";
                i -= 1;
            }
            // skip the start codon rows w+2, w+1, w (if the first loop ran
            // out instead, i is already negative)
            i -= 3;
            while i >= 0 {
                o += &row_text(&self.rbs_matrix[i as usize]);
                o += "\n";
                i -= 1;
            }
        }
        o
    }
}

/// A coded sequence as letters (table [`LETTERS`]) plus '\n' (0x40d31a).
pub fn seqstr(x: &[u8]) -> String {
    let mut s: String = x.iter().map(|&c| LETTERS[c as usize] as char).collect();
    s.push('\n');
    s
}
