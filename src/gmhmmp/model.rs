//! gmhmmp model files and every transformation applied before scoring
//! (single model reader 0x425e50, meta reader 0x426230 and helpers;
//! specification `re-notes/gmhmmp/model.md`). All values are verified bit
//! for bit against gdb dumps of the loaded objects of the original.
//!
//! Two file kinds, chosen by the first byte (`_` = meta):
//! * single model: `PHMM ...` key/value text (optionally obfuscated by
//!   subtracting '0' from every byte);
//! * meta model (MetaGeneMark): sections `__GC30` .. `__GC70`, each a
//!   complete model body closed by `end`; the section is chosen per record
//!   by the record's GC percent.
//!
//! After parsing, [`Model::validate`] checks the parameters and derives the
//! codon flag table and the logs; the length tables and the [`Markov`]
//! tables are built from the raw COD1, COD2 and NONC matrices.
//!
//! Numerics: all f64; `ln`/`exp` are glibc libm (Rust's f64 methods call
//! the same functions on this platform). A probability of exactly 0.0
//! becomes [`NEG`] instead of -inf; NaN passes through to `ln`. Sums of
//! four are `((x0 + x1) + x2) + x3`.

use super::num::{c_f, strtod, strtol_i32};

/// "log of zero" (bits 0xd4b249ad2594c37d)
pub const NEG: f64 = -1e100;
/// Markov fill values for windows with non-ACGT (stored constants of the
/// binary, Model+0xf8 / +0x100; hard-coded bits, not computed at run time)
pub const FILL_C: f64 = f64::from_bits(0xbfe77cacec67bada); // ln 0.48
pub const FILL_N: f64 = f64::from_bits(0xbfe4ecf732f00389); // ln 0.52

/// Load failure. `stdout`/`stderr` hold the messages; `direct_exit` means the
/// binary exits without printing "error on load parameters from file".
#[derive(Debug, Default)]
pub struct LoadError {
    pub stdout: String,
    pub stderr: String,
    pub direct_exit: bool,
}

/// Load error whose message goes to stdout.
fn out_err(msg: String) -> LoadError {
    LoadError { stdout: msg, ..Default::default() }
}

/// Marker for errors of the port itself (inputs on which the original
/// terminates or reads uninitialised memory): such a message is printed to stderr
/// as `gmhmmp-rs: error: model file <path>: <message>` and the program
/// exits with rc 1 (see [`port_or_out_err`]).
const PORT_ERR: &str = "\u{1}port:";

/// Port error text for [`Model::parse_keys`] / [`Model::from_single_text`].
fn port_err(msg: String) -> String {
    format!("{PORT_ERR}{msg}")
}

/// Map a parser/validation error to a [`LoadError`]: messages of the
/// original go to stdout (the caller adds "error on load parameters from
/// file"), port errors (see [`PORT_ERR`]) go to stderr with exit 1.
/// `section` names the meta model section, if any.
fn port_or_out_err(msg: String, path: &str, section: Option<&str>) -> LoadError {
    match msg.strip_prefix(PORT_ERR) {
        Some(m) => {
            let sec = section.map(|x| format!(", section {x}")).unwrap_or_default();
            LoadError { stderr: format!("gmhmmp-rs: error: model file {path}{sec}: {m}\n"), direct_exit: true, ..Default::default() }
        }
        None => out_err(msg),
    }
}

/// Matrix row count of the binary for Markov order `ordm` (0x411940,
/// 0x4103d0, 0x40edc0): `4 << ((2 * ordm) & 31)` in 32-bit int arithmetic,
/// i.e. 4^(ordm+1) for ordm <= 13, wrapping for larger orders (14 gives
/// 2^30, 15 gives 0, 16 gives 4, ...).
pub fn rows_for_order(ordm: i32) -> i32 {
    4i32 << ((ordm.wrapping_mul(2)) & 31)
}

/// `x == 0.0 ? NEG : ln(x)` (0x40fb70, 0x40fa10).
#[inline]
fn lg(x: f64) -> f64 {
    if x == 0.0 {
        NEG
    } else {
        x.ln()
    }
}

/// Markov object (0x170 bytes in the binary), built by 0x40fe00 from the
/// raw COD1, COD2 and NONC matrices (model.md 1.3, 6.5).
///
/// With k = ORDM and N = 4^(k+1), table index u is a (k+1)-mer with the
/// FIRST base as the most significant base-4 digit (the Markov image of
/// `predict.rs` produces exactly these indices). Every table has N+1
/// entries; entry N is the fill value used for windows containing N.
#[derive(Debug, Clone)]
pub struct Markov {
    pub order: i32,
    /// 4^(k+1)
    pub n: usize,
    /// Coding set from COD1 (set A, K+0x00..0x38): eight ln tables, see
    /// `coding_set`; each has n+1 entries, [n] = FILL_C.
    pub c1: Vec<Vec<f64>>,
    /// Coding set from COD2 (set B, K+0x40..0x78), same layout.
    pub c2: Vec<Vec<f64>>,
    /// ln NONC conditional `P(first base | next k bases)` (normalised over
    /// the most significant digit), [n] = FILL_N.
    pub nc: Vec<f64>,
    /// ln lower-order JOINT marginals of NONC: `nl[j]` holds (j+1)-mers
    /// (4^(j+1) entries + fill FILL_N); `nl[k]` = ln of the raw NONC.
    /// Used for noncoding gaps shorter than k+1 and for the joint term at
    /// the right end of a long gap.
    pub nl: Vec<Vec<f64>>,
    /// NONC order-0 marginals (A, C, G, T): sums of the raw NONC, linear,
    /// not normalised. RBS background.
    pub nonc0: [f64; 4],
    /// COD1 per-phase order-0 marginals (A, C, G, T per phase column),
    /// linear. RBS context-mode background after the start codon.
    pub cod1mono: [[f64; 4]; 3],
}

/// One gene model (0x110-byte Model object, model.md 1.2). Field names are
/// the model file keys.
#[derive(Debug, Clone)]
pub struct Model {
    /// RBSS: RBS motif width W; -1 = no RBS (also set by main for every
    /// record when `-r` is absent).
    pub rbss: i32,
    /// RBSM: W x 4 motif matrix, row-major (row r, column A,C,G,T). Rows
    /// are stored right to left: row 0 is the motif's last base.
    pub rbsm: Vec<f64>,
    /// RBSD: spacer distribution (RBSL values, 1 if RBSL < 0);
    /// renormalised in place on every record with `-r` (spacer mode).
    pub rbsd: Vec<f64>,
    /// RBSL: number of spacer lengths (>= 0, prokaryotic spacer mode), or
    /// minus the number of bases scored after the start codon (< 0, context
    /// mode, e.g. -3 in GeneMarkS-T trained models).
    pub rbsl: i32,
    /// ORDM: Markov order k.
    pub ordm: i32,
    /// ln P(ATG,GTG,TTG,CTG) or NEG
    pub log_start: [f64; 4],
    /// ATG_, GTG_, TTG_, CTG_ as read (not normalised).
    pub start: [f64; 4],
    /// TAA, TAG, TGA (raw); only `> 0` is ever tested (codon flags, R/Y
    /// resolution, TGA translation).
    pub stop: [f64; 3],
    /// codon flags, index 16*b0+4*b1+b2, [64] = 0
    pub flags: [u8; 65],
    /// MINC: minimum gene length incl. codons; shorter records are skipped.
    pub minc: i32,
    /// MAXN: maximum noncoding length (limit of the noncoding scan).
    pub maxn: i32,
    /// MAXC: maximum coding length (minimum size of the coding tables).
    pub maxc: i32,
    /// NDEC: noncoding length decay (geometric distribution).
    pub ndec: f64,
    /// CDEC: coding length decay (gamma-like distribution).
    pub cdec: f64,
    /// CDCD: only validated (>= 0), never used.
    pub cdcd: f64,
    /// ln CD1P / ln CD2P (if > 0): prior of coding class 1 / class 2,
    /// added to the coding log likelihoods.
    pub cd1p: f64,
    pub cd2p: f64,
    /// Raw matrices in file order: COD1/COD2 3 * 4^(k+1) values (row u =
    /// (k+1)-mer, 3 phase columns), NONC 4^(k+1).
    pub cod1: Vec<f64>,
    pub cod2: Vec<f64>,
    pub nonc: Vec<f64>,
    /// NAME (one token); printed as "(null)" if absent.
    pub name: Option<String>,
    /// Warnings of the port about this model (printed by [`load`]).
    pub warnings: Vec<String>,
    /// An ORDM key was read (the binary's matrix pointers are set).
    pub ordm_seen: bool,
    /// noncoding duration table (0xa8), len = Nn + 9 (last entry unused)
    pub dur_n: Vec<f64>,
    /// coding duration tables (0xb0 complete, 0xb8 partial), table[L-1] for length L
    pub dur_c: Vec<f64>,
    pub dur_p: Vec<f64>,
    pub markov: Markov,
}

impl Model {
    /// Defaults of the Model init 0x4110f0.
    fn empty() -> Model {
        Model {
            rbss: -1,
            rbsm: Vec::new(),
            rbsd: Vec::new(),
            rbsl: -1,
            ordm: -1,
            log_start: [NEG; 4],
            start: [0.0; 4],
            stop: [0.0; 3],
            flags: [0; 65],
            minc: 0,
            maxn: 0,
            maxc: 0,
            ndec: 0.0,
            cdec: 0.0,
            cdcd: 0.0,
            cd1p: 0.0,
            cd2p: 0.0,
            cod1: Vec::new(),
            cod2: Vec::new(),
            nonc: Vec::new(),
            name: None,
            warnings: Vec::new(),
            ordm_seen: false,
            dur_n: Vec::new(),
            dur_c: Vec::new(),
            dur_p: Vec::new(),
            markov: Markov {
                order: 0,
                n: 0,
                c1: Vec::new(),
                c2: Vec::new(),
                nc: Vec::new(),
                nl: Vec::new(),
                nonc0: [0.0; 4],
                cod1mono: [[0.0; 4]; 3],
            },
        }
    }

    /// 0x4103d0 key parser over a text body (model.md 2.3): `strtok` on
    /// " \t\n\r", exact key match, the last occurrence of a scalar wins.
    /// ORDM must precede the matrices, RBSS/RBSL must precede RBSM/RBSD
    /// (with RBSS <= 0 the RBS values are not consumed and then fail as
    /// unknown keys). Tokens `end` and anything containing `__GC` are
    /// ignored so that swallowed meta sections still parse.
    fn parse_keys(&mut self, text: &str) -> Result<(), String> {
        let mut toks = text.split(|c: char| c == ' ' || c == '\t' || c == '\n' || c == '\r').filter(|s| !s.is_empty());
        let mut have_ordm = false;
        // COD1, COD2, NONC read since the last ORDM (ORDM re-allocates them)
        let mut have_mat = [false; 3];
        while let Some(t) = toks.next() {
            let mut val = || -> Result<&str, String> {
                // Intended difference: a key without a value at the end of
                // the text makes the binary pass NULL to strtod/strtol/strcpy
                // and terminate (SIGSEGV, rc 139 [V]); the port stops with an error.
                toks.next().ok_or_else(|| port_err(format!("key {t} has no value (end of the model text); add the value or remove the key")))
            };
            match t {
                "PHMM" => {
                    // the version token is never used: a missing one is no
                    // abnormal termination in the binary (`PHMM` alone: "error, model order -1")
                    toks.next();
                }
                "NAME" => self.name = Some(val()?.to_string()),
                "ATG_" => self.start[0] = strtod(val()?),
                "GTG_" => self.start[1] = strtod(val()?),
                "TTG_" => self.start[2] = strtod(val()?),
                "CTG_" => self.start[3] = strtod(val()?),
                "TAA_" => self.stop[0] = strtod(val()?),
                "TAG_" => self.stop[1] = strtod(val()?),
                "TGA_" => self.stop[2] = strtod(val()?),
                "ORDM" => {
                    self.ordm = strtol_i32(val()?);
                    if self.ordm < 0 {
                        return Err(format!("error, model order {}\n", self.ordm));
                    }
                    // 0x411940 allocates 3N, 3N and N doubles with the int N
                    // of rows_for_order. Intended difference: for N = 2^30
                    // (ORDM 14, 30, ...) 3N overflows to a negative size and
                    // the binary aborts with std::bad_alloc (rc 134 [V]).
                    if rows_for_order(self.ordm).wrapping_mul(3) < 0 {
                        return Err(port_err(format!("ORDM {} is not a usable Markov order (gmhmmp supports orders 0 to 8)", self.ordm)));
                    }
                    // the matrices are filled as their tokens are read (the
                    // binary's arrays are uninitialised until then)
                    self.cod1.clear();
                    self.cod2.clear();
                    self.nonc.clear();
                    have_ordm = true;
                    have_mat = [false; 3];
                }
                "MINC" => self.minc = strtol_i32(val()?),
                "MAXC" => self.maxc = strtol_i32(val()?),
                "MAXN" => self.maxn = strtol_i32(val()?),
                "NDEC" => self.ndec = strtod(val()?),
                "CDEC" => self.cdec = strtod(val()?),
                "CDCD" => self.cdcd = strtod(val()?),
                "CD1P" => self.cd1p = strtod(val()?),
                "CD2P" => self.cd2p = strtod(val()?),
                "COD1" | "COD2" | "NONC" => {
                    if !have_ordm {
                        return Err("error, missing model order before matrix\n".into());
                    }
                    let rows = rows_for_order(self.ordm);
                    let (arr, what, cnt, mi) = match t {
                        "COD1" => (&mut self.cod1, "coding1", rows.wrapping_mul(3), 0),
                        "COD2" => (&mut self.cod2, "coding2", rows.wrapping_mul(3), 1),
                        _ => (&mut self.nonc, "noncoding", rows, 2),
                    };
                    arr.clear();
                    for _ in 0..cnt.max(0) {
                        match toks.next() {
                            Some(x) => arr.push(strtod(x)),
                            None => return Err(format!("error, out of matrix element for {what}\n")),
                        }
                    }
                    have_mat[mi] = true;
                }
                "COD3" => {
                    // third coding matrix: tokens skipped (its enable flag is never set)
                    if !have_ordm {
                        return Err("error, missing model order before matrix\n".into());
                    }
                    let n = rows_for_order(self.ordm).wrapping_mul(3).max(0);
                    for _ in 0..n {
                        if toks.next().is_none() {
                            return Err("error, out of matrix element for coding3\n".into());
                        }
                    }
                }
                "RBSS" => self.rbss = strtol_i32(val()?),
                "RBSL" => self.rbsl = strtol_i32(val()?),
                "RBSM" => {
                    if self.rbss > 0 {
                        // int count 4*RBSS as in the binary; when it wraps
                        // negative `new double[]` throws std::bad_alloc and the
                        // binary aborts. Intended difference: the port stops.
                        let cnt = self.rbss.wrapping_mul(4);
                        if cnt < 0 {
                            return Err(port_err(format!("RBSS {} is too large for an RBS motif matrix", self.rbss)));
                        }
                        let mut v = Vec::new();
                        for _ in 0..cnt {
                            match toks.next() {
                                Some(x) => v.push(strtod(x)),
                                None => return Err("error, out of matrix element for rbs\n".into()),
                            }
                        }
                        self.rbsm = v;
                    }
                }
                "RBSD" => {
                    if self.rbss > 0 {
                        let n = if self.rbsl >= 0 { self.rbsl } else { 1 };
                        let mut v = Vec::new();
                        for _ in 0..n {
                            match toks.next() {
                                Some(x) => v.push(strtod(x)),
                                None => return Err("error, out of matrix element for rbs\n".into()),
                            }
                        }
                        self.rbsd = v;
                    }
                }
                "end" => {}
                x if x.contains("__GC") => {}
                x => return Err(format!("error in model file, unknown: {x}\n")),
            }
        }
        if have_ordm {
            // Intended difference: after ORDM a missing COD1/COD2/NONC (or
            // one given only before a repeated ORDM) leaves the binary's
            // matrix uninitialised (uninitialised heap data, in practice often zeros).
            // The port uses zeros and warns.
            let rows = rows_for_order(self.ordm).max(0) as usize;
            for (mi, key) in ["COD1", "COD2", "NONC"].iter().enumerate() {
                if !have_mat[mi] {
                    let n = if mi < 2 { 3 * rows } else { rows };
                    let arr = match mi {
                        0 => &mut self.cod1,
                        1 => &mut self.cod2,
                        _ => &mut self.nonc,
                    };
                    *arr = vec![0.0; n];
                    self.warnings.push(format!("no {key} matrix after the last ORDM; using probability 0 for all its entries"));
                }
            }
        }
        self.ordm_seen = have_ordm;
        Ok(())
    }

    /// 0x4118f0 -> 0x411540 validation (model.md 3; checks and messages in
    /// the binary's order, including its spelling), then the codon flags
    /// (0x411360) and the logs (0x411250). NaN passes every `<`/`<=` test.
    fn validate(&mut self) -> Result<(), String> {
        if self.ordm < 0 {
            return Err(format!("error, model order {}\n", self.ordm));
        }
        for (i, nm) in ["ATG", "GTG", "TTG", "CTG"].iter().enumerate() {
            if self.start[i] < 0.0 {
                return Err(format!("error, {nm} probability {}\n", c_f(self.start[i])));
            }
        }
        if ((self.start[0] + self.start[1]) + self.start[2]) + self.start[3] == 0.0 {
            return Err("error, missing START codons in model\n".into());
        }
        for (i, nm) in ["TAA", "TAG", "TGA"].iter().enumerate() {
            if self.stop[i] < 0.0 {
                return Err(format!("error, {nm} probability {}\n", c_f(self.stop[i])));
            }
        }
        // TAG + TAA + TGA, in this order
        if (self.stop[1] + self.stop[0]) + self.stop[2] == 0.0 {
            return Err("error, missing STOP codons in model\n".into());
        }
        if self.minc < 1 {
            return Err(format!("error, min gene length in model {}\n", self.minc));
        }
        if self.maxc < 1 {
            return Err(format!("errror, max Coding length in model {}\n", self.maxc));
        }
        if self.cdec <= 0.0 {
            return Err(format!("error, Coding length distr parameter in model {}\n", c_f(self.cdec)));
        }
        if self.maxn < 1 {
            return Err(format!("errror, max Non Coding length in model {}\n", self.maxn));
        }
        if self.ndec <= 0.0 {
            return Err(format!("error, Non Coding length distr parameter in model {}\n", c_f(self.ndec)));
        }
        if self.cdcd < 0.0 {
            return Err(format!("error, coding to coding parameter in model {}\n", c_f(self.cdcd)));
        }
        if self.cd1p <= 0.0 {
            return Err(format!("error, model parameter for Coding 1 {}\n", c_f(self.cd1p)));
        }
        if self.cd2p <= 0.0 {
            return Err(format!("error, model parameter for Coding 2 {}\n", c_f(self.cd2p)));
        }
        if self.minc < self.ordm {
            return Err(format!("error, min gene length {} less than model order {}\n", self.minc, self.ordm));
        }
        // the binary tests the three matrix pointers, set by ORDM
        if !self.ordm_seen {
            return Err("error, no matrix loaded\n".into());
        }
        if self.rbss >= 1 && (self.rbsd.is_empty() || self.rbsm.is_empty()) {
            return Err("error in rbs model\n".into());
        }
        // 0x411360 flags: forward codon and its reverse complement, only for
        // probabilities > 0 (bits: see predict.rs module docs)
        let f = &mut self.flags;
        *f = [0; 65];
        if self.start[0] > 0.0 {
            f[0x0e] = 0x20; // ATG
            f[0x13] = 0x08; // CAT
        }
        if self.start[1] > 0.0 {
            f[0x2e] = 0x21; // GTG
            f[0x11] = 0x09; // CAC
        }
        if self.start[2] > 0.0 {
            f[0x3e] = 0x22; // TTG
            f[0x10] = 0x0a; // CAA
        }
        if self.start[3] > 0.0 {
            f[0x1e] = 0x23; // CTG
            f[0x12] = 0x0b; // CAG
        }
        if self.stop[0] > 0.0 {
            f[0x30] = 0x10; // TAA
            f[0x3c] = 0x04; // TTA
        }
        if self.stop[2] > 0.0 {
            f[0x38] = 0x10; // TGA
            f[0x34] = 0x04; // TCA
        }
        if self.stop[1] > 0.0 {
            f[0x32] = 0x10; // TAG
            f[0x1c] = 0x04; // CTA
        }
        // 0x411250 logs (no normalisation of the start probabilities)
        for i in 0..4 {
            if self.start[i] > 0.0 {
                self.log_start[i] = self.start[i].ln();
            }
        }
        if self.cd1p > 0.0 {
            self.cd1p = self.cd1p.ln();
        }
        if self.cd2p > 0.0 {
            self.cd2p = self.cd2p.ln();
        }
        Ok(())
    }

    /// 0x411b70(m, P_req): coding length tables (model.md 6.3), indexed
    /// `table[L-1]` for gene length L (codons included), ln values:
    ///
    /// * `dur_c` (complete genes): `L^2 exp(-L/C) / (2 C^3)` on multiples
    ///   of 3, normalised over L = MINC..P;
    /// * `dur_p` (partial genes): `(L^2 + 2LC + 2C^2) exp(-L/C) / (2 C^3)`
    ///   on every L (tail form), normalised the same way.
    ///
    /// C = CDEC. Stateful: rebuilt only if `p_req >= MAXC` or the tables do
    /// not exist yet, then with P = max(p_req, MAXC). Since the
    /// normalisation range depends on P, a record with a long ORF changes
    /// all values for the following records (per model, also per GC model
    /// in meta mode) [V]. Entry [P] stays 0.0 (never written by the binary).
    pub fn coding_lengths(&mut self, p_req: i32) {
        if p_req < self.maxc && !self.dur_c.is_empty() {
            return;
        }
        let p = p_req.max(self.maxc) as usize;
        let c = self.cdec;
        let minc = self.minc.max(1) as usize;
        let mut f = vec![0.0f64; p + 1];
        let mut sum = 0.0;
        for l in minc..=p {
            if l % 3 == 0 {
                // recomputed per length, as in the binary
                let k = 1.0 / (((c + c) * c) * c);
                let v = ((k * l as f64) * l as f64) * (-(l as f64) / c).exp();
                f[l - 1] = v;
                sum = v + sum;
            } else {
                f[l - 1] = 0.0;
            }
        }
        for v in f.iter_mut().take(p) {
            *v = if *v == 0.0 { NEG } else { (*v / sum).ln() };
        }
        let mut g = vec![0.0f64; p + 1];
        let mut sum = 0.0;
        for l in minc..=p {
            let d = (c + c) * c;
            let k2 = 1.0 / (c * d);
            // L*L is an i32 multiply in the binary (wraps for L > 46340)
            let sq = (l as i32).wrapping_mul(l as i32) as f64;
            let v = (k2 * ((sq + (2 * l) as f64 * c) + d)) * (-(l as f64) / c).exp();
            g[l - 1] = v;
            sum = v + sum;
        }
        for v in g.iter_mut().take(p) {
            *v = if *v == 0.0 { NEG } else { (*v / sum).ln() };
        }
        self.dur_c = f;
        self.dur_p = g;
    }

    /// 0x411f60(m, P2_req): noncoding length table (model.md 6.4), indexed
    /// by `d = noncoding length + 7`: `[8 + i] = i ln(e) + ln(1-e)` with
    /// `e = exp(-1/NDEC)` (geometric, not normalised), `[0] = NEG`,
    /// `[1..8] = [8]` (overlapping or adjacent codons). Same caching rule
    /// as [`Model::coding_lengths`] with MAXN; a rebuild only extends the
    /// table. Size `max(P2_req, MAXN) + 9`, last entry unused.
    pub fn noncoding_lengths(&mut self, p2_req: i32) {
        if p2_req < self.maxn && !self.dur_n.is_empty() {
            return;
        }
        let nn = p2_req.max(self.maxn) as usize;
        let p = (-1.0 / self.ndec).exp();
        let a = (1.0 - p).ln();
        let b = p.ln();
        let mut t = vec![0.0f64; nn + 9];
        for i in 0..nn {
            t[8 + i] = i as f64 * b + a;
        }
        t[0] = NEG;
        for j in 1..8 {
            t[j] = t[8];
        }
        self.dur_n = t;
    }

    fn build_markov(&mut self) {
        self.markov = build_markov(self.ordm, &self.cod1, &self.cod2, &self.nonc);
    }

    /// Single model text (after magic check / de-obfuscation): parse,
    /// validate, initial length tables (P = 0, i.e. MAXC/MAXN), Markov.
    fn from_single_text(text: &str) -> Result<Model, String> {
        let mut m = Model::empty();
        m.parse_keys(text)?;
        m.validate()?;
        // Intended difference: 0x411b70 zeroes F[0 .. MINC-2] of a table
        // with MAXC+1 entries; for MINC > MAXC + 2 that writes past the heap
        // block and the binary aborts in malloc or terminates with SIGSEGV (rc 134/139
        // [V], MAXC 15000: MINC 15002 runs, 15003 aborts). The port stops.
        if m.minc as i64 > m.maxc as i64 + 2 {
            return Err(port_err(format!("MINC {} (minimum gene length) is larger than MAXC {} (maximum coding length) + 2; use MINC <= MAXC", m.minc, m.maxc)));
        }
        m.coding_lengths(0);
        m.noncoding_lengths(0);
        // Orders 9 to 13 pass the loader but every record then stops with
        // ": fatal error, [order] > [short int]" (0x40c170) before any
        // Markov table is read, so the tables are not built. Intended
        // difference: for ORDM >= 15 the binary's table size has wrapped
        // (rows_for_order) and its Markov builder 0x40fe00 reads far past
        // the matrices (about 30 s and GiBs for ORDM 16); the port stops.
        if m.ordm > 13 {
            return Err(port_err(format!("ORDM {} is not a usable Markov order (gmhmmp supports orders 0 to 8)", m.ordm)));
        }
        if m.ordm <= 8 {
            m.build_markov();
        }
        Ok(m)
    }
}

// ---------------------------------------------------------------- Markov
//
// Primitives of 0x40fe00 (model.md 6.5). Index u in 0..N has base-4 digits
// d_0 (least significant) .. d_k (most significant = first base, "MSD").

/// Reverse complement of a (k+1)-mer index (0x40ee60).
fn rc_index(mut u: usize, order: i32) -> usize {
    let mut v = 0;
    for _ in 0..=order {
        v = (3 - (u & 3)) + 4 * v;
        u >>= 2;
    }
    v
}

/// 0x40ef10: normalise each group {j, j+N/4, j+2N/4, j+3N/4} (varying first
/// base) into a new array: P(first base | following k bases). A zero sum
/// gives 0.0 for all four.
fn cond_msd(src: &[f64]) -> Vec<f64> {
    let n = src.len();
    let q = n >> 2;
    let mut out = vec![0.0; n];
    for j in 0..q {
        let s = ((src[j] + src[j + q]) + src[j + 2 * q]) + src[j + 3 * q];
        for d in 0..4 {
            out[j + d * q] = if s == 0.0 { 0.0 } else { src[j + d * q] / s };
        }
    }
    out
}

/// 0x40f000: normalise each group {4j..4j+3} (varying last base) in place:
/// P(last base | preceding k bases).
fn norm_lsd(a: &mut [f64]) {
    for g in a.chunks_mut(4) {
        let s = ((g[0] + g[1]) + g[2]) + g[3];
        for x in g.iter_mut() {
            *x = if s == 0.0 { 0.0 } else { *x / s };
        }
    }
}

/// `out[u] = src[rc(u)]` (0x40ee60).
fn rcmap(src: &[f64], order: i32) -> Vec<f64> {
    (0..src.len()).map(|u| src[rc_index(u, order)]).collect()
}

/// 0x40f0c0: marginals of lower order: `lv[k]` = src, `lv[j-1]` sums out the
/// last base of `lv[j]`; `lv[0]` = the 4 single-base sums.
fn levels(src: &[f64], order: i32) -> Vec<Vec<f64>> {
    let mut lv = vec![Vec::new(); order as usize + 1];
    lv[order as usize] = src.to_vec();
    for k in (1..=order as usize).rev() {
        let hi = &lv[k];
        let lo: Vec<f64> = (0..hi.len() / 4).map(|j| ((hi[4 * j] + hi[4 * j + 1]) + hi[4 * j + 2]) + hi[4 * j + 3]).collect();
        lv[k - 1] = lo;
    }
    lv
}

/// Take `lg` of every entry and append the fill slot.
fn finish(mut v: Vec<f64>, fill: f64) -> Vec<f64> {
    for x in v.iter_mut() {
        *x = lg(*x);
    }
    v.push(fill);
    v
}

/// One coding set (0x40f2b0) from COD1 or COD2, eight ln tables:
///
/// * s[0], s[1], s[2]: forward conditional `P(base | next k bases)` for
///   (k+1)-mers whose first base is at codon position 1, 2, 3 (COD
///   columns 0, 1, 2);
/// * s[3]: raw joint (k+1)-mer of the column whose first base has the
///   codon position of `stop.pos - k` (the joint term at the right end of
///   a forward gene): column 1, 0, 2 for k mod 3 = 1, 2, 0;
/// * s[4], s[5], s[6]: reverse-strand conditional tables: the columns
///   normalised over their LAST base, then reverse complemented (so again
///   conditional on the bases to the right in forward coordinates), and
///   rotated by k mod 3 so that `predict.rs` can use them in the fixed
///   order 4, 5, 6;
/// * s[7]: raw joint of column 0 reverse complemented (joint term at the
///   right end of a reverse gene, next to its start codon).
fn coding_set(cod: &[f64], order: i32) -> Vec<Vec<f64>> {
    let n = cod.len() / 3;
    let mut c: Vec<Vec<f64>> = (0..3).map(|ph| (0..n).map(|i| cod[3 * i + ph]).collect()).collect();
    let r = order % 3;
    let s3 = match r {
        1 => c[1].clone(),
        2 => c[0].clone(),
        _ => c[2].clone(),
    };
    let s0 = cond_msd(&c[0]);
    let s1 = cond_msd(&c[1]);
    let s2 = cond_msd(&c[2]);
    // before norm_lsd: raw joint values
    let s7 = rcmap(&c[0], order);
    for col in c.iter_mut() {
        norm_lsd(col);
    }
    let mut s4 = rcmap(&c[0], order);
    let mut s5 = rcmap(&c[1], order);
    let mut s6 = rcmap(&c[2], order);
    if r == 1 {
        let (a, b, cc) = (s6, s4, s5);
        s4 = a;
        s5 = b;
        s6 = cc;
    } else if r == 2 {
        let (a, b, cc) = (s5, s6, s4);
        s4 = a;
        s5 = b;
        s6 = cc;
    }
    vec![s0, s1, s2, s3, s4, s5, s6, s7].into_iter().map(|v| finish(v, FILL_C)).collect()
}

/// 0x40fe00: the whole Markov object. Build order follows the binary
/// (COD1 phase compositions, set A, set B, NONC tables).
fn build_markov(order: i32, cod1: &[f64], cod2: &[f64], nonc: &[f64]) -> Markov {
    let n = 4usize << (2 * order);
    let mut cod1mono = [[0.0; 4]; 3];
    for (ph, mono) in cod1mono.iter_mut().enumerate() {
        let col: Vec<f64> = (0..n).map(|i| cod1[3 * i + ph]).collect();
        let lv = levels(&col, order);
        mono.copy_from_slice(&lv[0][..4]);
    }
    let c1 = coding_set(cod1, order);
    let c2 = coding_set(cod2, order);
    let lv = levels(nonc, order);
    let mut nonc0 = [0.0; 4];
    nonc0.copy_from_slice(&lv[0][..4]);
    let nc = finish(cond_msd(nonc), FILL_N);
    let nl = lv.into_iter().map(|v| finish(v, FILL_N)).collect();
    Markov { order, n, c1, c2, nc, nl, nonc0, cod1mono }
}

// ---------------------------------------------------------------- loading

/// The loaded `-m` file: one model, or 41 GC models (meta).
pub enum Models {
    Single(Box<Model>),
    /// sections for gc 30..=70, index gc-30
    Meta(Vec<Model>),
}

impl Models {
    /// Section for a GC percent 0..100 (0x425f70 / 0x4260d0 vectors): with
    /// all 41 keys present the 101-entry lookup is `clamp(gc, 30, 70)`.
    pub fn index_for_gc(gc: i32) -> usize {
        (gc.clamp(30, 70) - 30) as usize
    }
}

/// 0x4058b0: first byte decides the reader ('_' = meta).
pub fn is_meta(bytes: &[u8]) -> bool {
    bytes.first() == Some(&b'_')
}

/// 0x425e50 (single) or 0x426230 (meta). Errors carry the binary's
/// messages; the caller adds "error on load parameters from file" unless
/// `direct_exit`.
pub fn load(path: &str) -> Result<Models, LoadError> {
    let data = std::fs::read(path);
    let meta = data.as_ref().map(|d| is_meta(d)).unwrap_or(false);
    if meta {
        return load_meta(path, &data.unwrap()).map(Models::Meta);
    }
    // single: 0x40d110
    let data = match data {
        Ok(d) => d,
        Err(_) => {
            let msg = if std::path::Path::new(path).exists() {
                format!("error: reading file <{path}>\n")
            } else {
                format!("error: file <{path}> is not found\n")
            };
            return Err(out_err(msg));
        }
    };
    // behaviour of the binary: "size 0" fires when exactly 1 byte was read (an
    // empty file instead fails the magic check below)
    if data.len() == 1 {
        return Err(out_err(format!("error: file <{path}> has size 0\n")));
    }
    let mut buf = data;
    // C string
    if let Some(z) = buf.iter().position(|&b| b == 0) {
        buf.truncate(z);
    }
    // 0x410fc0 magic: "PHMM", or an obfuscated file (every byte - '0'),
    // which is de-obfuscated as a whole
    if !buf.starts_with(b"PHMM") {
        let ob = buf.len() >= 4 && buf[..4].iter().zip(b"PHMM").all(|(&a, &b)| a.wrapping_add(0x30) == b);
        if ob {
            for b in buf.iter_mut() {
                *b = b.wrapping_add(0x30);
            }
            if let Some(z) = buf.iter().position(|&b| b == 0) {
                buf.truncate(z);
            }
        } else {
            return Err(out_err("error, not valid model file format\n".into()));
        }
    }
    let text = String::from_utf8_lossy(&buf);
    let m = Model::from_single_text(&text).map_err(|e| port_or_out_err(e, path, None))?;
    for w in &m.warnings {
        eprintln!("gmhmmp-rs: warning: model file {path}: {w}");
    }
    Ok(Models::Single(Box::new(m)))
}

/// Meta model (0x426230, model.md 2.2): `#` comments blanked, sections split
/// by [`split_sections`], then `__GC30`..`__GC70` each parsed like a single
/// model body (no PHMM magic check). For a duplicate label the last one
/// wins.
fn load_meta(path: &str, raw: &[u8]) -> Result<Vec<Model>, LoadError> {
    let mut buf = raw.to_vec();
    if let Some(z) = buf.iter().position(|&b| b == 0) {
        buf.truncate(z);
    }
    // 0x419460: comments
    let mut i = 0;
    while i < buf.len() {
        if buf[i] == b'#' {
            while i < buf.len() && buf[i] != b'\n' && buf[i] != b'\r' {
                buf[i] = b' ';
                i += 1;
            }
        } else {
            i += 1;
        }
    }
    let text = String::from_utf8_lossy(&buf).into_owned();
    let sections = split_sections(&text).map_err(|e| LoadError {
        stderr: format!("{e}error in format of file {path}\n"),
        direct_exit: true,
        ..Default::default()
    })?;
    let mut models = Vec::with_capacity(41);
    for gc in 30..=70 {
        let key = format!("__GC{gc}");
        let body = match sections.iter().rev().find(|(k, _)| *k == key) {
            Some((_, b)) => b.clone(),
            None => {
                return Err(LoadError { stderr: format!("error: missing parameter {key}\n"), ..Default::default() });
            }
        };
        let m = Model::from_single_text_meta(&body).map_err(|e| port_or_out_err(e, path, Some(&key)))?;
        for w in &m.warnings {
            eprintln!("gmhmmp-rs: warning: model file {path}, section {key}: {w}");
        }
        models.push(m);
    }
    Ok(models)
}

impl Model {
    /// meta sections: no magic check
    fn from_single_text_meta(text: &str) -> Result<Model, String> {
        Model::from_single_text(text)
    }
}

/// 0x41a430: label -> body pairs, in file order (model.md 2.2). Kept as in
/// the original: the body ends at the first SUBSTRING "end" after the
/// label (must follow whitespace, else "error3"), the next label is
/// searched from the end of the LABEL (so a section without `end`
/// swallows the following ones, which are still parsed on their own too),
/// a body in double quotes loses them. The error strings are written to
/// stderr without a newline, followed by "error in format of file".
fn split_sections(text: &str) -> Result<Vec<(String, String)>, String> {
    let b = text.as_bytes();
    let ws = |c: u8| c == b' ' || c == b'\t' || c == b'\n' || c == b'\r';
    let mut out = Vec::new();
    let mut pos = match b.iter().position(|&c| !ws(c)) {
        Some(p) => p,
        None => return Err("error1".into()),
    };
    if !text[pos..].starts_with("__") {
        return Err("error2".into());
    }
    loop {
        // label
        let mut le = pos;
        while le < b.len() && !ws(b[le]) {
            le += 1;
        }
        let label = text[pos..le].to_string();
        // body ends at the first "end" after the label
        let body_end = match text[le..].find("end") {
            Some(off) => {
                let e = le + off;
                if e == 0 || !ws(b[e - 1]) {
                    return Err("error3".into());
                }
                e
            }
            None => b.len(),
        };
        let mut body = text[le.min(body_end)..body_end].trim_matches(|c: char| c == ' ' || c == '\t' || c == '\n' || c == '\r').to_string();
        if body.len() > 1 && body.starts_with('"') && body.ends_with('"') {
            body = body[1..body.len() - 1].to_string();
        }
        out.push((label, body));
        // next label searched from the end of the LABEL
        match text[le..].find("__") {
            Some(off) => {
                let p = le + off;
                if p > 0 && !ws(b[p - 1]) {
                    return Err("error3".into());
                }
                pos = p;
            }
            None => break,
        }
    }
    Ok(out)
}
