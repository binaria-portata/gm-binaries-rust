//! Model file reader and parameter transformations.
//!
//! Port of gmhmme3 0x412fe0 (load) / 0x413390 (read and tokenise) / 0x4123c0
//! (all keys) and their subtree. See `re-notes/model.md` for the
//! binary-level specification (object offsets, order of operations). All
//! arithmetic is f64 with the same summation order as the binary; this was
//! verified bit-exact against gdb dumps of the loaded model.
//!
//! # The model file
//!
//! A GeneMark `.mod` file is a list of `$KEY value...` entries (`#` starts a
//! comment). Values are numbers, quoted strings, 1-based index lists
//! (`idx value` pairs) or k-mer tables (`KMER v0 v1 ...` rows). A file
//! starting with `?` is obfuscated (bytes shifted by 0x30).
//!
//! # What the reader produces ([`Model`])
//!
//! Every probability is converted to a natural log here, with [`NEG`]
//! (-1e120) as "log of zero". The later stages only add these values.
//!
//! * Length (duration) distributions ([`LenDist`]) for intergenic, intron
//!   and the four exon types, normalised, plus a "tail" table (normalised
//!   survival function) for introns truncated by the sequence end.
//! * Markov chains of order K = `$ORDER` from the `$MARKOV` table, whose
//!   columns are joint (K+1)-mer distributions: columns 0..2 = coding, one
//!   per codon position ([`Coding`]); column 3 = intron, column 4 =
//!   intergenic ([`NonCoding`]). From each column the reader derives the
//!   conditional tables the score tracks need (see [`norm_first`]),
//!   reverse-complement tables and lower-order tables.
//! * Signal site matrices ([`Site`]): positional weight matrices with an
//!   internal Markov order; with `$NORM_SITES` they become log-odds against
//!   the Markov background of the flanking regions ([`site_norm`]).
//! * Exon phase tables, gene-structure transition probabilities
//!   (`ToSingleGene` ...), initial state probabilities, stop codon switches,
//!   GC donor parameters, and optionally the branch point block
//!   ([`BranchPoint`], `$BP`).
//!
//! # Errors
//!
//! The reader stops at the first error. Most errors print a specific message
//! and then main adds `error in model file <F>`; a few call `exit(1)`
//! directly in the original (no extra line), see [`ModelError::fatal_exit`].
//! All messages go to stdout.

use std::fmt;
use std::path::Path;

/// "log of zero" used throughout the reader and the DP (bit pattern
/// 0xd8d8c8dac6a0342a). Finite on purpose: sums of a few NEG values stay
/// finite, and `x - 1e120` with ordinary `x` rounds to exactly NEG.
pub const NEG: f64 = -1e120;
/// ln 3: phase tables of initial/terminal exons are stored relative to the
/// uniform distribution over 3 phases (0x598050).
const LN3: f64 = 1.0986122886681098;
/// ln 9: the same for the 9 phase pairs of internal exons.
const LN9: f64 = 2.1972245773362196;

/// Errors carry the exact text gmhmme3 prints to stdout before
/// `error in model file <F>`.
#[derive(Debug)]
pub struct ModelError {
    /// lines to print, in order (each with its own newline)
    pub messages: Vec<String>,
    /// true if the binary calls exit(1) directly (no "error in model file").
    pub fatal_exit: bool,
    /// Error of the port itself (printed to stderr as `gmhmm-rs: error:
    /// ...`, exit code 1) for model values on which the original terminates
    /// (abort or SIGSEGV) or computes with undefined memory.
    /// None for the original's own messages.
    pub port: Option<String>,
}

impl fmt::Display for ModelError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        for m in &self.messages {
            f.write_str(m)?;
        }
        Ok(())
    }
}

/// Error after which main prints `error in model file <F>`.
fn fail(msg: String) -> ModelError {
    ModelError { messages: vec![msg], fatal_exit: false, port: None }
}

/// Error on which the original calls exit(1) itself.
fn fatal(msg: String) -> ModelError {
    ModelError { messages: vec![msg], fatal_exit: true, port: None }
}

/// Port error (see [`ModelError::port`]): `msg` names the model file, the
/// key and what is wrong with its value.
fn port_err(msg: String) -> ModelError {
    ModelError { messages: Vec::new(), fatal_exit: true, port: Some(msg) }
}

/// Size of a k-mer table as the original computes it, `4 << (2*order &
/// 31)` in 32-bit int arithmetic (0x411c50, 0x411a90, 0x41bce0, 0x41c1f0):
/// `4^(order+1)` for orders 0..14, but the shift wraps for other values, so
/// for example order 15 or -1 gives 0 and order 16 gives 4.
pub fn c_table_size(order: i32) -> i32 {
    4i32.wrapping_shl((order.wrapping_mul(2) & 31) as u32)
}

/// `vec![0.0; n]` that reports an allocation failure instead of aborting.
/// The original allocates its tables with `operator new`, which throws
/// `std::bad_alloc` (abort) when the memory is not available.
fn zeros(n: usize, file: &str, what: &str) -> Result<Vec<f64>, ModelError> {
    let mut v: Vec<f64> = Vec::new();
    if v.try_reserve_exact(n).is_err() {
        return Err(port_err(format!(
            "model file {file}: not enough memory for {what} ({n} values); check the value of the key that sets its size"
        )));
    }
    v.resize(n, 0.0);
    Ok(v)
}

/// The getters' common failure: missing key, duplicate key, wrong number of
/// tokens, bad syntax.
fn param_err(key: &str) -> ModelError {
    fail(format!("error reading parameter {key}\n"))
}

// ---------------------------------------------------------------------------
// number parsing with C prefix semantics (strtol / strtod)

/// C `isspace` in the "C" locale (includes \v and \f).
fn is_c_space(b: u8) -> bool {
    matches!(b, b' ' | b'\t' | b'\n' | 0x0b | 0x0c | b'\r')
}

/// `strtol(s, NULL, 10)`: leading whitespace, optional sign, the longest
/// run of digits; anything after it is ignored, no digits give 0. The
/// result saturates like a 64-bit `long`; callers truncate it to C `int`
/// (`as i32`) where the binary stores 4 bytes.
pub fn strtol(s: &str) -> i64 {
    let b = s.as_bytes();
    let mut i = 0;
    while i < b.len() && is_c_space(b[i]) {
        i += 1;
    }
    let mut neg = false;
    if i < b.len() && (b[i] == b'+' || b[i] == b'-') {
        neg = b[i] == b'-';
        i += 1;
    }
    let mut v: i64 = 0;
    while i < b.len() && b[i].is_ascii_digit() {
        v = v.saturating_mul(10).saturating_add((b[i] - b'0') as i64);
        i += 1;
    }
    if neg {
        -v
    } else {
        v
    }
}

/// Exactly representable powers of ten for Clinger's fast path.
const POW10: [f64; 23] = [1e0, 1e1, 1e2, 1e3, 1e4, 1e5, 1e6, 1e7, 1e8, 1e9, 1e10, 1e11, 1e12, 1e13, 1e14, 1e15, 1e16, 1e17, 1e18, 1e19, 1e20, 1e21, 1e22];

/// Clinger's fast path for a complete plain decimal token ([+-]digits[.digits]):
/// mantissa < 2^53 and |power of ten| <= 22 give an exact operand pair, so one
/// IEEE multiplication/division is correctly rounded, i.e. identical to strtod.
///
/// This is a speed trick only (model files hold hundreds of thousands of
/// such numbers). Anything else (exponent, too many digits, trailing
/// uninitialised data) returns None and goes through the general path of [`strtod`].
#[inline]
fn fast_decimal(b: &[u8]) -> Option<f64> {
    let mut i = 0;
    let neg = match b.first() {
        Some(b'-') => {
            i = 1;
            true
        }
        Some(b'+') => {
            i = 1;
            false
        }
        _ => false,
    };
    let mut mantissa: u64 = 0;
    let mut sig_digits = 0;
    let mut frac_digits = 0i32;
    let mut seen_dot = false;
    let mut any_digit = false;
    while i < b.len() {
        let c = b[i];
        if c.is_ascii_digit() {
            any_digit = true;
            if mantissa == 0 && c == b'0' {
                // leading zeros do not count as significant digits
            } else {
                sig_digits += 1;
                // 19 digits always fit a u64; more would overflow
                if sig_digits > 19 {
                    return None;
                }
                mantissa = mantissa * 10 + (c - b'0') as u64;
            }
            // every digit after the dot (zeros included) scales by 1/10
            if seen_dot {
                frac_digits += 1;
            }
        } else if c == b'.' && !seen_dot {
            seen_dot = true;
        } else {
            return None;
        }
        i += 1;
    }
    if !any_digit || mantissa > (1u64 << 53) || frac_digits > 22 {
        return None;
    }
    let v = mantissa as f64 / POW10[frac_digits as usize];
    Some(if neg { -v } else { v })
}

/// `strtod(s, NULL)`: value of the longest valid decimal floating-point
/// prefix of `s` (0.0 if there is none), correctly rounded like glibc.
///
/// Three tiers, fastest first: [`fast_decimal`]; Rust's correctly rounded
/// parser for plain numeric tokens; and a manual scan for the longest
/// numeric prefix (leading whitespace, sign, digits, '.', exponent) or
/// inf/infinity/nan, which emulates C's prefix semantics for odd input.
pub fn strtod(s: &str) -> f64 {
    let bs = s.as_bytes();
    if let Some(v) = fast_decimal(bs) {
        return v;
    }
    // other plain tokens: Rust's parser is correctly rounded like glibc
    if !bs.is_empty() && bs.iter().all(|&c| c.is_ascii_digit() || matches!(c, b'.' | b'-' | b'+' | b'e' | b'E')) {
        if let Ok(v) = s.parse::<f64>() {
            return v;
        }
    }
    let t = s.trim_start_matches(|c: char| c.is_ascii() && is_c_space(c as u8));
    if let Ok(v) = t.parse::<f64>() {
        return v;
    }
    let b = t.as_bytes();
    let mut i = 0;
    if i < b.len() && (b[i] == b'+' || b[i] == b'-') {
        i += 1;
    }
    let lower = t[i..].to_ascii_lowercase();
    for (word, val) in [("infinity", f64::INFINITY), ("inf", f64::INFINITY), ("nan", f64::NAN)] {
        if lower.starts_with(word) {
            return if b.first() == Some(&b'-') { -val } else { val };
        }
    }
    // mantissa: digits [. digits]; a lone "." is not a number
    let digits_start = i;
    while i < b.len() && b[i].is_ascii_digit() {
        i += 1;
    }
    if i < b.len() && b[i] == b'.' {
        i += 1;
        while i < b.len() && b[i].is_ascii_digit() {
            i += 1;
        }
    }
    if i == digits_start || (i == digits_start + 1 && b[digits_start] == b'.') {
        return 0.0;
    }
    // the exponent is only part of the prefix if it has at least one digit
    // ("2e" parses as 2)
    let mut end = i;
    if i < b.len() && (b[i] == b'e' || b[i] == b'E') {
        let mut j = i + 1;
        if j < b.len() && (b[j] == b'+' || b[j] == b'-') {
            j += 1;
        }
        let exp_start = j;
        while j < b.len() && b[j].is_ascii_digit() {
            j += 1;
        }
        if j > exp_start {
            end = j;
        }
    }
    t[..end].parse::<f64>().unwrap_or(0.0)
}

/// Table logarithm (0x410e10, 0x41bce0, inline in 0x4111f0): 0 -> NEG,
/// otherwise ln. No sign check: a negative entry gives NaN, as in the
/// original.
fn log0(x: f64) -> f64 {
    if x == 0.0 {
        NEG
    } else {
        x.ln()
    }
}

/// LOG (0x411040), used for scalar probabilities: x > 0 -> ln x, x == 0 ->
/// NEG, negative or NaN -> fatal error (exit 1).
fn checked_log(x: f64) -> Result<f64, ModelError> {
    if x > 0.0 {
        Ok(x.ln())
    } else if x == 0.0 {
        Ok(NEG)
    } else {
        Err(fatal(format!("error, negative argument in LOG function : {}\n", c_fmt_f(x))))
    }
}

/// printf("%f"): six decimals, glibc spelling of NaN and infinities.
pub fn c_fmt_f(x: f64) -> String {
    if x.is_nan() {
        if x.is_sign_negative() { "-nan".into() } else { "nan".into() }
    } else if x.is_infinite() {
        if x < 0.0 { "-inf".into() } else { "inf".into() }
    } else {
        format!("{x:.6}")
    }
}

/// Whitespace tokenizer over a value (C isspace), without allocation.
struct TokenIter<'a> {
    s: &'a str,
    i: usize,
}

impl<'a> TokenIter<'a> {
    fn new(s: &'a str) -> Self {
        TokenIter { s, i: 0 }
    }
}

impl<'a> Iterator for TokenIter<'a> {
    type Item = &'a str;
    #[inline]
    fn next(&mut self) -> Option<&'a str> {
        let b = self.s.as_bytes();
        let mut i = self.i;
        while i < b.len() && is_c_space(b[i]) {
            i += 1;
        }
        if i >= b.len() {
            self.i = i;
            return None;
        }
        let s0 = i;
        while i < b.len() && !is_c_space(b[i]) {
            i += 1;
        }
        self.i = i;
        Some(&self.s[s0..i])
    }
}

// ---------------------------------------------------------------------------
// key list (0x413390, 0x413210) and typed getters (0x413a40..0x414520)

/// A key/value range pair into [`KeyList::text`].
type Entry = (std::ops::Range<usize>, std::ops::Range<usize>);

/// The model file split into `$KEY value` entries, in file order.
struct KeyList {
    /// the decoded file text (comments blanked); entries are ranges into it
    text: String,
    /// (key range, value range) per entry; the value runs up to the next `$`
    entries: Vec<Entry>,
    /// model file name, for the port's own error messages
    file: String,
}

impl KeyList {
    /// 0x413390 + 0x413210. None means `error in file format` (the original
    /// rejects the whole file).
    fn parse(raw: &[u8]) -> Option<KeyList> {
        let mut buf = raw.to_vec();
        // obfuscated format: "?" + one more byte, then every byte shifted
        if buf.first() == Some(&b'?') {
            buf[0] = b' ';
            if buf.len() > 1 {
                buf[1] = b' ';
            }
            for b in buf.iter_mut().skip(2) {
                *b = b.wrapping_add(0x30);
            }
        }
        // comments: '#' up to (not including) the next '\n' or '\r' become
        // spaces, so a '$' inside a comment is harmless
        let mut i = 0;
        while let Some(off) = memchr::memchr(b'#', &buf[i..]) {
            let j = i + off;
            let mut k = j;
            while k < buf.len() && buf[k] != b'\n' && buf[k] != b'\r' {
                buf[k] = b' ';
                k += 1;
            }
            i = k;
            if i >= buf.len() {
                break;
            }
        }
        // a NUL in the file ends the C string
        if let Some(z) = memchr::memchr(0, &buf) {
            buf.truncate(z);
        }
        // Bytes that are not valid UTF-8 are kept in escaped form (they only
        // matter in values printed back, such as a `_TAG` name in a
        // warning). Every byte the parser looks at below is ASCII, and an
        // escaped byte becomes non-ASCII bytes only, so the scan is the same.
        let text = crate::sequence::decode(&buf);
        let buf = text.into_bytes();
        // the first non-space character must be '$'
        let mut p = 0;
        while p < buf.len() && is_c_space(buf[p]) {
            p += 1;
        }
        if p >= buf.len() || buf[p] != b'$' {
            return None;
        }
        // every '$' (anywhere, even inside quotes) starts a new entry and
        // must follow a whitespace character (which the original overwrites
        // with NUL to end the previous value)
        let mut starts = vec![p];
        for d in memchr::memchr_iter(b'$', &buf[p + 1..]) {
            let d = d + p + 1;
            if !is_c_space(buf[d - 1]) {
                return None;
            }
            starts.push(d);
        }
        // valid UTF-8 by construction (decoded above)
        let text = String::from_utf8(buf).unwrap_or_default();
        let buf = text.as_bytes();
        let mut entries = Vec::with_capacity(starts.len());
        for (n, &s) in starts.iter().enumerate() {
            // segment ends at the whitespace char before the next '$' (made NUL)
            let end = if n + 1 < starts.len() { starts[n + 1] - 1 } else { buf.len() };
            let seg0 = s + 1;
            let seg_end = end.max(seg0);
            let seg = &buf[seg0..seg_end];
            // key = up to the first whitespace; a key without following
            // whitespace, an empty key ("$ ") or an empty value is an error
            let klen = seg.iter().position(|&b| is_c_space(b))?;
            if klen == 0 {
                return None;
            }
            let mut v = klen;
            while v < seg.len() && is_c_space(seg[v]) {
                v += 1;
            }
            if v >= seg.len() {
                return None;
            }
            entries.push((seg0..seg0 + klen, seg0 + v..seg_end));
        }
        Some(KeyList { text, entries, file: String::new() })
    }

    #[inline]
    fn key(&self, e: &Entry) -> &str {
        &self.text[e.0.clone()]
    }

    /// 0x4136c0: does the key exist (no message)? Used for optional keys.
    fn has(&self, key: &str) -> bool {
        self.entries.iter().any(|e| self.key(e) == key)
    }

    /// Raw value text of a key. Missing and duplicate keys are both errors
    /// (the getters scan the rest of the list after the first match).
    fn raw(&self, key: &str) -> Result<&str, ModelError> {
        let mut it = self.entries.iter().filter(|e| self.key(e) == key);
        let first = it.next().ok_or_else(|| param_err(key))?;
        if it.next().is_some() {
            return Err(param_err(key));
        }
        Ok(&self.text[first.1.clone()])
    }

    /// The value must be exactly one token.
    fn single_token(&self, key: &str) -> Result<&str, ModelError> {
        let v = self.raw(key)?;
        let mut t = v.split(|c: char| c.is_ascii() && is_c_space(c as u8)).filter(|s| !s.is_empty());
        let tok = t.next().ok_or_else(|| param_err(key))?;
        if t.next().is_some() {
            return Err(param_err(key));
        }
        Ok(tok)
    }

    /// 0x413a40: int, `(int)strtol` of one token.
    fn int(&self, key: &str) -> Result<i32, ModelError> {
        Ok(strtol(self.single_token(key)?) as i32)
    }

    /// 0x413b30: double, `strtod` of one token.
    fn dbl(&self, key: &str) -> Result<f64, ModelError> {
        Ok(strtod(self.single_token(key)?))
    }

    /// 0x413c10: bool, `strtod(v) != 0.0`.
    fn boolean(&self, key: &str) -> Result<bool, ModelError> {
        Ok(self.dbl(key)? != 0.0)
    }

    /// Optional bool with a default for a missing key.
    fn opt_bool(&self, key: &str, default: bool) -> Result<bool, ModelError> {
        if self.has(key) {
            self.boolean(key)
        } else {
            Ok(default)
        }
    }

    /// 0x413d20: string, `"..."` with nothing after the closing quote.
    fn string(&self, key: &str) -> Result<String, ModelError> {
        let v = self.raw(key)?;
        let v = v.trim_end_matches(|c: char| c.is_ascii() && is_c_space(c as u8));
        let rest = v.strip_prefix('"').ok_or_else(|| param_err(key))?;
        let close = rest.find('"').ok_or_else(|| param_err(key))?;
        if close + 1 != rest.len() {
            return Err(param_err(key));
        }
        Ok(rest[..close].to_string())
    }

    /// All whitespace-separated tokens of a value.
    fn tokens<'a>(&'a self, key: &str) -> Result<Vec<&'a str>, ModelError> {
        // C isspace: ' ' \t \n \v \f \r (split_ascii_whitespace lacks \v)
        let v = self.raw(key)?;
        let b = v.as_bytes();
        let mut out = Vec::with_capacity(b.len() / 8);
        let mut i = 0;
        while i < b.len() {
            while i < b.len() && is_c_space(b[i]) {
                i += 1;
            }
            let s0 = i;
            while i < b.len() && !is_c_space(b[i]) {
                i += 1;
            }
            if i > s0 {
                out.push(&v[s0..i]);
            }
        }
        Ok(out)
    }

    /// 0x413fa0 / 0x413e30: 1-based index list (`idx value` pairs) into an
    /// array of n doubles, zero-initialised. Index < 1 is an error, index > n
    /// only a warning (printed at once, entry ignored); a repeated index
    /// overwrites. `warn` only counts warnings (empty strings), it is
    /// never printed.
    fn index_list(&self, key: &str, n: usize, warn: &mut Vec<String>) -> Result<Vec<f64>, ModelError> {
        // allocated before the key is looked up, as in the original
        let mut arr = zeros(n, &self.file, &format!("the table {key}"))?;
        let v = self.raw(key)?;
        let mut it = TokenIter::new(v);
        while let Some(idx_tok) = it.next() {
            let val_tok = match it.next() {
                Some(x) => x,
                None => return Err(param_err(key)),
            };
            let idx = strtol(idx_tok);
            let val = strtod(val_tok);
            if idx < 1 {
                return Err(ModelError {
                    messages: vec![
                        format!("gmhmme error: index {idx} is out of range (less then 1) for key {key}\n"),
                        format!("error reading parameter {key}\n"),
                    ],
                    fatal_exit: false,
                    port: None,
                });
            }
            if idx as usize > n {
                // printed immediately, like the binary (before any later error)
                crate::out(&format!("gmhmme warning: index {idx} is out of range for max {n} for key {key}\n"));
                warn.push(String::new());
            } else {
                arr[idx as usize - 1] = val;
            }
        }
        Ok(arr)
    }

    /// 0x414520 / 0x414350: k-mer rows `KMER v1 .. v_ncol` (the `$MARKOV`
    /// layout). Returns `ncol` raw (not logged) tables of n+1 doubles,
    /// n = 4^(order+1), indexed by the k-mer index (first base most
    /// significant); the extra last slot (N slot) stays 0 here. Missing rows
    /// stay 0, a repeated row overwrites, a row with too few values is an
    /// error.
    ///
    /// `order_key` names the key that holds `order` (for messages). The
    /// table size is the original's [`c_table_size`]: for orders outside
    /// 0..14 it wraps. Size 0 (e.g. order 15 or -1) rejects every row, which
    /// gives the original's `error reading parameter KEY`. Any other order
    /// outside 0..14 makes the original terminate after reading the rows (it
    /// allocates and fills tables of 4^(order+1) entries, which fails or
    /// overruns); the port checks the rows the same way and then stops with
    /// an error naming the order key.
    fn kmer_rows(&self, key: &str, order_key: &str, order: i32, ncol: usize) -> Result<Vec<Vec<f64>>, ModelError> {
        let n = c_table_size(order);
        let in_range = (0..=14).contains(&order);
        let mut cols = Vec::new();
        if in_range {
            for c in 0..ncol {
                cols.push(zeros(n as usize + 1, &self.file, &format!("column {} of {key}", c + 1))?);
            }
        }
        let t = self.tokens(key)?;
        let mut i = 0;
        while i < t.len() {
            let idx = kmer_index(t[i]).ok_or_else(|| param_err(key))?;
            // the second test is "fewer than ncol values left after the
            // k-mer" (i + ncol >= t.len(); t.len() >= 1 inside the loop)
            if idx < 0 || idx >= n || i + ncol >= t.len() {
                return Err(param_err(key));
            }
            for (c, col) in cols.iter_mut().enumerate() {
                col[idx as usize] = strtod(t[i + 1 + c]);
            }
            i += 1 + ncol;
        }
        if !in_range {
            return Err(port_err(format!(
                "model file {}: ${order_key} {order} is out of range; the Markov order must be between 0 and 14 \
                 (not handled by the original)",
                self.file
            )));
        }
        Ok(cols)
    }
}

/// Index of a k-mer key in a model table (0x4140d0 / 0x414350): base-4
/// number, A=0 C=1 G=2 T=3, first character most significant. Only upper
/// case ACGT is accepted. No length check: a shorter k-mer maps to a small
/// index, as in the original. The index is a 32-bit int that wraps for
/// k-mers of 16 or more letters (e.g. `C` followed by 16 `A` gives 0), and
/// the callers reject `idx < 0 || idx >= size`, as in the original.
fn kmer_index(s: &str) -> Option<i32> {
    let mut v = 0i32;
    for b in s.bytes() {
        let d = match b {
            b'A' => 0,
            b'C' => 1,
            b'G' => 2,
            b'T' => 3,
            _ => return None,
        };
        v = v.wrapping_mul(4).wrapping_add(d);
    }
    Some(v)
}

// ---------------------------------------------------------------------------
// length distributions (0x4111f0)

/// Explicit duration distribution of one segment type (LenDist, 0x40 bytes).
///
/// `tab[i]` is log P(length index i+1) for the file's 1-based indices
/// 1..=MAX, renormalised over MIN..=MAX (entries below MIN are zeroed, i.e.
/// NEG). `tail[i]` is the log of the normalised survival function
/// `sum_{j >= i} tab[j]` over i < MAX-1, used for introns that run into a
/// sequence end.
///
/// The Viterbi looks entries up by the map distance between two records,
/// which differs from the true segment length by a type-specific constant.
/// The original absorbs that constant into shifted base pointers
/// (`tab_base = tab - p4`, `tail_base = tail - p5`) with exclusive bounds;
/// [`LenDist::base`], [`LenDist::lo`], [`LenDist::hi`] etc. reproduce them.
/// The resulting length conventions are listed in re-notes/viterbi.md
/// section 2.
#[derive(Debug, Clone)]
pub struct LenDist {
    /// `<NAME>_MAX` (+0x00)
    pub max: i32,
    /// `<NAME>_MIN` (+0x04)
    pub min: i32,
    /// log P, tab[i] = file index i+1 (MAX entries)
    pub tab: Vec<f64>,
    /// log of the normalised survival function, MAX-1 entries
    pub tail: Vec<f64>,
    /// base pointer offset of `tab`: tab_base[x] = tab[x - p4]
    pub p4: i32,
    /// base pointer offset of `tail`: tail_base[x] = tail[x - p5]
    pub p5: i32,
}

impl LenDist {
    /// exclusive lower bound of tab_base (struct +0x20)
    /// (32-bit int arithmetic of the original: wraps for MIN/MAX near the
    /// int limits)
    pub fn lo(&self) -> i32 {
        self.min.wrapping_add(self.p4).wrapping_sub(2)
    }
    /// exclusive upper bound of tab_base (struct +0x24)
    pub fn hi(&self) -> i32 {
        self.max.wrapping_add(self.p4)
    }
    /// exclusive lower bound of tail_base (struct +0x28)
    pub fn tail_lo(&self) -> i32 {
        self.p5 - 1
    }
    /// exclusive upper bound of tail_base (struct +0x2c)
    pub fn tail_hi(&self) -> i32 {
        self.max.wrapping_sub(1).wrapping_add(self.p5)
    }
    /// `tab_base[x]`; callers check `lo() < x < hi()` first. With MIN <= 0
    /// that range starts below the table (lengths below 1): the original
    /// reads the heap in front of it there (only reachable through the BP
    /// distributions, which warn at load time). The port returns NEG
    /// (impossible length).
    #[inline]
    pub fn base(&self, x: i32) -> f64 {
        let i = x - self.p4;
        if i < 0 {
            return NEG;
        }
        self.tab[i as usize]
    }
    /// `tail_base[x]`; callers check `tail_lo() < x < tail_hi()` first.
    #[inline]
    pub fn tail_base(&self, x: i32) -> f64 {
        self.tail[(x - self.p5) as usize]
    }
}

/// 0x4111f0(M, off, NAME, p4, p5): read and transform one length
/// distribution. Keys: `NAME_MAX`, `NAME_MIN`, `NAME_TAG` (name of the key
/// holding the index-list table; required even for parametric types),
/// `NAME_TYPE` (CONSTANT, EXP, GAMMA, or anything else = use the table).
fn read_len_dist(kl: &KeyList, name: &str, p4: i32, p5: i32, warn: &mut Vec<String>) -> Result<LenDist, ModelError> {
    let max = kl.int(&format!("{name}_MAX"))?;
    let min = kl.int(&format!("{name}_MIN"))?;
    if max < min {
        // the binary prints MIN and the address of the struct via cout; the
        // address is not reproducible, the port prints 0x0
        return Err(fail(format!("{min}\n0x0\n")));
    }
    let tag = kl.string(&format!("{name}_TAG"))?;
    // MAX < 0: the original's table allocation throws
    // std::bad_array_new_length (abort); MAX == 0 passes here but the tail
    // table of MAX-1 entries throws the same a few lines further down.
    let max_err = || {
        port_err(format!(
            "model file {}: {name}_MAX is {max}; a length distribution needs a maximum of at least 1 \
             (not handled by the original)",
            kl.file
        ))
    };
    if max < 0 {
        return Err(max_err());
    }
    let mut tab = kl.index_list(&tag, max as usize, warn)?;
    let typ = kl.string(&format!("{name}_TYPE"))?;
    // first 0-based index inside [MIN, MAX]; i below is the 0-based index
    // (file index i+1), used as is in the formulas, as in the original
    // (the original's fill loops run while `MIN-1 < MAX` in int arithmetic,
    // so MIN == INT_MIN fills nothing)
    let fills = min.wrapping_sub(1) < max;
    let first = if fills { (min as i64 - 1).max(0) as usize } else { tab.len() };
    // entries zeroed below (memset of MIN-1 entries, only for MIN > 1)
    let below_min = (min as i64 - 1).max(0) as usize;
    // The parametric types fill indices MIN-1 .. MAX-1; with MIN <= 0 the
    // original writes 1 - MIN values in front of its table (heap
    // out-of-bounds write, it aborts when the table is freed at the end, after the
    // output has been written). The port fills only the table itself, which
    // gives the same output, and says so.
    if min <= 0 && fills && max > 0 && matches!(typ.as_str(), "CONSTANT" | "EXP" | "GAMMA") {
        crate::warn(&format!(
            "model file {}: {name}_MIN is {min} with {name}_TYPE {typ}; lengths below 1 do not exist and are \
             left out (not handled by the original)",
            kl.file
        ));
    }
    match typ.as_str() {
        "CONSTANT" => {
            for v in tab.iter_mut().skip(first) {
                *v = 1.0;
            }
        }
        "EXP" => {
            let d = kl.dbl(&format!("{name}_EXP_D"))?;
            for (i, v) in tab.iter_mut().enumerate().skip(first) {
                *v = (-(i as f64) / d).exp();
            }
        }
        "GAMMA" => {
            let a = kl.dbl(&format!("{name}_GAMMA_A"))?;
            let b = kl.dbl(&format!("{name}_GAMMA_B"))?;
            for (i, v) in tab.iter_mut().enumerate().skip(first) {
                *v = (-(i as f64) / b).exp() * (i as f64).powf(a);
            }
        }
        _ => {}
    }
    if max == 0 {
        return Err(max_err());
    }
    // lengths below MIN are impossible (overwrites file values too)
    for v in tab.iter_mut().take(below_min) {
        *v = 0.0;
    }
    // normalise; sequential left-to-right sum, as in the original
    let mut sum = 0.0;
    for &v in &tab {
        sum += v;
    }
    if sum == 0.0 {
        return Err(fatal("error in length distribution\n".into()));
    }
    for v in tab.iter_mut() {
        *v /= sum;
    }
    // survival function accumulated from the top:
    // tail[i] = tab[MAX-1] + tab[MAX-2] + ... + tab[i]
    let n = (max - 1).max(0) as usize;
    let mut tail = vec![0.0; n];
    let mut survival = tab[max as usize - 1];
    for i in (0..n).rev() {
        survival += tab[i];
        tail[i] = survival;
    }
    let mut sum = 0.0;
    for &v in &tail {
        sum += v;
    }
    if sum == 0.0 {
        return Err(fatal("error in length distribution\n".into()));
    }
    for v in tail.iter_mut() {
        *v /= sum;
    }
    for v in tab.iter_mut() {
        *v = log0(*v);
    }
    for v in tail.iter_mut() {
        *v = log0(*v);
    }
    Ok(LenDist { max, min, tab, tail, p4, p5 })
}

// ---------------------------------------------------------------------------
// Markov tables (0x40f960, 0x4100f0, 0x410460)
//
// All helpers work on raw joint (K+1)-mer probabilities (index = base-4
// number, first base most significant) of n = 4^(K+1) entries plus an
// N slot at index n; logs are taken last, in `finish`.

/// "nf" table: P(first nt | following K nts), 0x410d70 + 0x410b20.
///
/// The k-mers sharing their last K bases are `i + j*q` (j = first base,
/// q = 4^K). Each such group of four is divided by its sum, summed in the
/// original's order `((t[i+q] + t[i]) + t[i+2q]) + t[i+3q]`; an all-zero
/// group stays 0. The score tracks read sequences right to left with these
/// tables (see re-notes/viterbi.md 3.1).
fn norm_first(t: &[f64], k: i32) -> Vec<f64> {
    let q = 4usize.pow(k as u32);
    let mut out = vec![0.0; 4 * q + 1];
    for i in 0..q {
        let s = ((t[i + q] + t[i]) + t[i + 2 * q]) + t[i + 3 * q];
        if s == 0.0 {
            continue;
        }
        for j in 0..4 {
            out[i + j * q] = t[i + j * q] / s;
        }
    }
    out[4 * q] = t[4 * q];
    out
}

/// Reverse complement of a k-mer index: reverse the `digits` base-4 digits
/// and complement each (`3 - d`).
#[inline]
pub fn rc_index(mut v: usize, digits: i32) -> usize {
    let mut r = 0;
    for _ in 0..digits {
        r = (3 - (v & 3)) + r * 4;
        v >>= 2;
    }
    r
}

/// 0x410c80: table for the reverse strand, `out[u] = t[rc(u)]` (N slot copied).
fn revcomp(t: &[f64], k: i32) -> Vec<f64> {
    let n = 4usize.pow((k + 1) as u32);
    let mut out = vec![0.0; n + 1];
    for (u, o) in out.iter_mut().enumerate().take(n) {
        *o = t[rc_index(u, k + 1)];
    }
    out[n] = t[n];
    out
}

/// 0x410e80 + 0x4107e0: lower-order tables. `low[K]` = the joint table
/// itself (not normalised); `low[j]` for j < K = P(last nt | j previous
/// nts), 4^(j+1)+1 entries. All marginals are computed first (summing out
/// the last base, `((a0+a1)+a2)+a3`), then each lower table is normalised
/// in groups of four (no zero check, as in the original). Used as
/// background for the site normalisation ([`site_norm`]).
fn lower_orders(t: &[f64], k: i32) -> Vec<Vec<f64>> {
    let mut low: Vec<Vec<f64>> = vec![Vec::new(); (k + 1) as usize];
    low[k as usize] = t.to_vec();
    for j in (1..=k as usize).rev() {
        let m = 4usize.pow(j as u32);
        let src = &low[j];
        let mut out = vec![0.0; m + 1];
        for i in 0..m {
            out[i] = ((src[4 * i] + src[4 * i + 1]) + src[4 * i + 2]) + src[4 * i + 3];
        }
        out[m] = src[4 * m];
        low[j - 1] = out;
    }
    for j in (0..k as usize).rev() {
        let n = 4usize.pow(j as u32 + 1);
        let p = &mut low[j];
        for g in (0..n).step_by(4) {
            let s = ((p[g] + p[g + 1]) + p[g + 2]) + p[g + 3];
            for x in &mut p[g..g + 4] {
                *x /= s;
            }
        }
    }
    low
}

/// Final step for every Markov table: set the N slot (last entry) to the
/// model's N probability (`N_PROB_COD` or `N_PROB_NON`), then take logs
/// element-wise (0 -> NEG). Windows containing an N have k-mer index n and
/// read this slot.
fn finish(mut t: Vec<f64>, nprob: f64) -> Vec<f64> {
    let last = t.len() - 1;
    t[last] = nprob;
    for v in t.iter_mut() {
        *v = log0(*v);
    }
    t
}

/// Non-coding Markov chain (MARKOV column 3 at 0x948 with reverse complement,
/// or column 4 at 0x930 forward only). Tables are logs, length 4^(K+1)+1, the
/// last entry is the N slot.
#[derive(Debug, Clone)]
pub struct NonCoding {
    /// log joint P((K+1)-mer)
    pub fwd: Vec<f64>,
    /// log P(first nt | following K nts)
    pub fwd_nf: Vec<f64>,
    /// `fwd` of the reverse complement (empty for column 4)
    pub rc: Vec<f64>,
    /// `fwd_nf` of the reverse complement (empty for column 4)
    pub rc_nf: Vec<f64>,
    /// logs of [`lower_orders`], index 0..=K
    pub low: Vec<Vec<f64>>,
}

/// 0x4100f0 (with reverse complement) / 0x410460 (without).
fn build_noncoding(raw: &[f64], k: i32, nprob: f64, with_rc: bool) -> NonCoding {
    let fwd_nf = finish(norm_first(raw, k), nprob);
    let (rc, rc_nf) = if with_rc {
        let r = revcomp(raw, k);
        (finish(r.clone(), nprob), finish(norm_first(&r, k), nprob))
    } else {
        (Vec::new(), Vec::new())
    };
    let low = lower_orders(raw, k).into_iter().map(|t| finish(t, nprob)).collect();
    NonCoding { fwd: finish(raw.to_vec(), nprob), fwd_nf, rc, rc_nf, low }
}

/// Coding Markov block at M+0x7e8: three codon-position columns (MARKOV
/// columns 0..2), each with the same derived tables as [`NonCoding`]. The
/// score tracks walk these columns as a ring, advancing one column per
/// base (`tracks::ring`).
#[derive(Debug, Clone)]
pub struct Coding {
    pub fwd: [Vec<f64>; 3],
    pub fwd_nf: [Vec<f64>; 3],
    pub rc: [Vec<f64>; 3],
    pub rc_nf: [Vec<f64>; 3],
    pub low: [Vec<Vec<f64>>; 3],
    /// Markov order K
    pub order: i32,
}

impl Coding {
    /// Reverse ring column order (records Ra, Rb, Rc of 0x40f960):
    /// a = (-K) mod 3, b = (1-K) mod 3, c = (2-K) mod 3 (for K = 5: 1, 2,
    /// 0). Values as found in the original's pointer layout (verified for
    /// K = 5 with gdb).
    pub fn rc_ring(&self) -> [usize; 3] {
        let k = self.order;
        [(-k).rem_euclid(3) as usize, (1 - k).rem_euclid(3) as usize, (2 - k).rem_euclid(3) as usize]
    }
}

/// 0x40f960: build the coding block from MARKOV columns 0..2.
fn build_coding(raw: &[Vec<f64>], k: i32, nprob: f64) -> Coding {
    let build_column = |c: usize| -> (Vec<f64>, Vec<f64>, Vec<f64>, Vec<f64>, Vec<Vec<f64>>) {
        let r = revcomp(&raw[c], k);
        (
            finish(raw[c].clone(), nprob),
            finish(norm_first(&raw[c], k), nprob),
            finish(r.clone(), nprob),
            finish(norm_first(&r, k), nprob),
            lower_orders(&raw[c], k).into_iter().map(|t| finish(t, nprob)).collect(),
        )
    };
    let (a0, b0, c0, d0, e0) = build_column(0);
    let (a1, b1, c1, d1, e1) = build_column(1);
    let (a2, b2, c2, d2, e2) = build_column(2);
    Coding {
        fwd: [a0, a1, a2],
        fwd_nf: [b0, b1, b2],
        rc: [c0, c1, c2],
        rc_nf: [d0, d1, d2],
        low: [e0, e1, e2],
        order: k,
    }
}

// ---------------------------------------------------------------------------
// site matrices (0x411a90, 0x41bce0, 0x41c3e0)

/// Signal site model (0x50 bytes): a position-specific matrix of WIDTH
/// positions with an internal Markov order `order`.
///
/// Column p scores the base at window position p given the `order`
/// preceding bases: `mat[p*n + idx]` with idx = k-mer index of (context,
/// base) and the base least significant. The window is placed relative to
/// the site position by MARGIN (see `tracks::site_fwd` for the exact
/// placement). Values are logs, or log-odds after NORM_SITES.
#[derive(Debug, Clone)]
pub struct Site {
    /// key prefix `<S>` (e.g. "INI"), for the port's messages
    pub name: String,
    /// `<S>_MARGIN`: window offset (bases left of the site's anchor)
    pub margin: i32,
    /// `<S>_WIDTH`: window length
    pub width: i32,
    /// `<S>_ORDER`: Markov order inside the window
    pub order: i32,
    /// 4^(order+1); for an order outside 0..14 the original's wrapped size
    /// ([`c_table_size`]). The scoring mask is `n - 1`.
    pub n: usize,
    /// mat[p*n + idx], idx = k-mer of (order context nts, nt at p), nt at p
    /// least significant; plus one trailing raw 0.0 element.
    pub mat: Vec<f64>,
    /// reverse complement of the final `mat` (window reversed, k-mers
    /// complemented); read by nothing in the gene finder, kept for the
    /// model dump
    pub rc: Vec<f64>,
    /// `<S>_MARGIN_P` (default MARGIN); only used for the tr output windows
    pub margin_p: i32,
    /// `<S>_WIDTH_P` (default WIDTH); only used for the tr output windows
    pub width_p: i32,
}

/// Reverse-complement matrix: `rc[(W-1-p)*n + rc(u)] = mat[p*n + u]`.
fn site_rc(mat: &[f64], n: usize, w: i32, o: i32) -> Vec<f64> {
    let w = w as usize;
    let mut rc = vec![0.0; n * w];
    for p in 0..w {
        for u in 0..n {
            rc[(w - 1 - p) * n + rc_index(u, o + 1)] = mat[p * n + u];
        }
    }
    rc
}

/// 0x411a90 + 0x41bce0: read one site model `name` (keys `name_WIDTH`,
/// `_MARGIN`, `_ORDER`, `_MAT`, optional `_WIDTH_P`, `_MARGIN_P`).
fn read_site(kl: &KeyList, name: &str) -> Result<Site, ModelError> {
    let width = kl.int(&format!("{name}_WIDTH"))?;
    let margin = kl.int(&format!("{name}_MARGIN"))?;
    let order = kl.int(&format!("{name}_ORDER"))?;
    let key = format!("{name}_MAT");
    let file = &kl.file;
    // 0x411a90 -> 0x414290: the matrix has `n * WIDTH + 1` doubles with the
    // original's wrapped table size n (see c_table_size), computed in 32-bit
    // int and allocated before the key is looked up. A negative count
    // throws std::bad_array_new_length (abort); a product beyond the int
    // range wraps to a too small allocation that the rows then overrun.
    let nc = c_table_size(order);
    let cells = nc.wrapping_mul(width).wrapping_add(1);
    if width < 0 {
        return Err(port_err(format!(
            "model file {file}: {name}_WIDTH is {width}; a site width cannot be negative (not handled by the original)"
        )));
    }
    if cells < 0 || nc as i64 * width as i64 + 1 != cells as i64 {
        return Err(port_err(format!(
            "model file {file}: {name}_WIDTH {width} with {name}_ORDER {order} gives a matrix of more than 2^31 values \
             (not handled by the original)"
        )));
    }
    let n = nc as usize;
    let w = width as usize;
    let mut mat = zeros(cells as usize, file, &format!("the matrix {key}"))?;
    let t = kl.tokens(&key)?;
    // whether the value ends with white space (see the width 0 case below)
    let trailing_space = kl.raw(&key)?.as_bytes().last().is_some_and(|&c| is_c_space(c));
    let mut i = 0;
    while i < t.len() {
        let idx = kmer_index(t[i]).ok_or_else(|| param_err(&key))?;
        // A row needs its k-mer plus W values, i.e. tokens i ..= i + w. The
        // original (0x4140d0) reads the W values one by one and fails at the
        // first missing token, which prints "error reading parameter
        // <S>_MAT" and stops; same rule as in `kmer_rows`. It also fails
        // when a k-mer is the very last token of the value (its end-of-token
        // search 0x408670 returns NULL at the terminating NUL); for W >= 1
        // that is covered by the count test, for W == 0 it is tested here.
        if idx < 0 || idx as usize >= n || i + w >= t.len() || (w == 0 && i + 1 == t.len() && !trailing_space) {
            return Err(param_err(&key));
        }
        for j in 0..w {
            mat[idx as usize + j * n] = strtod(t[i + 1 + j]);
        }
        i += 1 + w;
    }
    // Orders outside 0..14 with a non-zero wrapped size (n == 0 has already
    // failed above, every row is rejected). 0x41bce0 builds the reverse
    // matrix with order+1 digits: for order >= 16 that is more digits than
    // the table has and the last row is written in front of the buffer (heap
    // out-of-bounds write, the original aborts). For order <= -2 no digit is used and
    // everything else works on the wrapped size, i.e. the site behaves like
    // one of order `order & 15`; that case is reproduced (NORM_SITES is
    // checked in site_norm).
    if order > 14 && width > 0 {
        return Err(port_err(format!(
            "model file {file}: {name}_ORDER is {order}; a site order must be between 0 and 14 \
             (not handled by the original)"
        )));
    }
    // log of the W*n matrix entries; the extra element stays raw 0.0
    for v in mat.iter_mut().take(n * w) {
        *v = log0(*v);
    }
    let rc = site_rc(&mat, n, width, order);
    let mut s = Site { name: name.to_string(), margin, width, order, n, mat, rc, margin_p: margin, width_p: width };
    let kw = format!("{name}_WIDTH_P");
    if kl.has(&kw) {
        s.width_p = kl.int(&kw)?;
    }
    let km = format!("{name}_MARGIN_P");
    if kl.has(&km) {
        s.margin_p = kl.int(&km)?;
    }
    Ok(s)
}

/// 0x41c3e0: subtract the Markov background (NORM_SITES), turning the site
/// matrix into log-odds against the content model of the flanking regions.
///
/// The window is split into three parts: MARGIN positions left of the
/// consensus, `clen` consensus positions (GT/AG: 2, codons: 3, branch point:
/// 0) and `R = W - MARGIN - clen` positions right of it. One side is
/// non-coding (background `bg` = lower-order tables of MARKOV column 3 or 4),
/// the other coding (`cod` = lower-order tables of the three codon columns)
/// when `left_coding` says so; consensus rows get 0 (nothing subtracted).
/// Each row p gets the background table of order `min(p, cap)` tiled over
/// the n matrix columns (`cap` = largest usable order, limited by the site
/// order, the stretch length and K).
///
/// As in the original, for exact output: in the coding part the codon column is NOT
/// rotated for the first `last+1` positions (the ones using growing orders);
/// only the later positions advance it by one per base. Both "never rotate"
/// and "always rotate" disagree with the binary (re-notes/model.md 3.4).
///
/// NEG entries (log 0) are left unchanged. The `rc` matrix is rebuilt from
/// the result.
#[allow(clippy::too_many_arguments)]
fn site_norm(
    s: &mut Site,
    bg: &[Vec<f64>],
    cod: &[Vec<Vec<f64>>; 3],
    phase: usize,
    left_coding: bool,
    clen: i32,
    k: i32,
    file: &str,
) -> Result<(), ModelError> {
    let n = s.n;
    let w = s.width;
    let mg = s.margin;
    let site_order = s.order;
    let r = w - mg - clen;
    // Rows written below: 0 .. MARGIN-1 (if MARGIN > 0) and MARGIN+clen ..
    // WIDTH-1 (if R > 0). The original's background matrix has WIDTH rows
    // and no bounds check: MARGIN > WIDTH, or MARGIN + clen < 0 with R > 0,
    // writes outside it (out-of-bounds write: the original aborts or terminates
    // with SIGSEGV). A negative site order (accepted for an unnormalised
    // site, see read_site) makes it index its table of background pointers
    // with a negative order (undefined memory).
    let name = &s.name;
    if (mg > 0 && mg > w) || (r > 0 && mg + clen < 0) {
        return Err(port_err(format!(
            "model file {file}: {name}_MARGIN {mg} does not fit {name}_WIDTH {w}; with $NORM_SITES the margin must be \
             between {} and the width (not handled by the original)",
            -clen
        )));
    }
    if site_order < 0 && (mg > 0 || r > 0) {
        return Err(port_err(format!(
            "model file {file}: {name}_ORDER is {site_order}; with $NORM_SITES a site order must be between 0 and 14 \
             (not handled by the original)"
        )));
    }
    let mut b = vec![0.0; n * w as usize];
    // largest usable order for a stretch of l positions
    let cap = |l: i32| -> i32 {
        let mut o = site_order;
        if l <= site_order {
            o = l - 1;
        }
        if k < o {
            o = k;
        }
        o
    };
    // tile a lower-order table (4^(o+1) entries) over the n columns of a row
    let mut fill = |row: i32, table: &[f64], o: i32| {
        let m = 4usize.pow((o + 1) as u32);
        for j in 0..n {
            b[row as usize * n + j] = table[j % m];
        }
    };
    // codon column c -> c+1 (pointer chain cod0 -> cod1 -> cod2 -> cod0)
    let nxt = |c: usize| (c + 1) % 3;
    if !left_coding {
        // [non-coding | consensus | coding]: acceptors, start codon, BP
        if mg > 0 {
            let om = cap(mg);
            for p in 0..mg {
                let o = p.min(om);
                fill(p, &bg[o as usize], o);
            }
        }
        if r > 0 {
            let om = cap(r);
            let mut c = phase;
            let base = mg + clen;
            let last = om.min(r - 1);
            // no rotation of c here (as in the original, see above)
            for q in 0..=last {
                fill(base + q, &cod[c][q as usize], q);
            }
            for q in last + 1..r {
                c = nxt(c);
                fill(base + q, &cod[c][om as usize], om);
            }
        }
    } else {
        // [coding | consensus | non-coding]: donors, stop codons
        if mg > 0 {
            let om = cap(mg);
            // codon column of window position 0, from the phase at the
            // consensus counted back over MARGIN bases
            let mut c = ((phase as i32 + 3 + ((mg / 3) * 3 - mg)) % 3) as usize;
            let last = om.min(mg - 1);
            // no rotation of c here (as in the original, see above)
            for p in 0..=last {
                fill(p, &cod[c][p as usize], p);
            }
            for p in last + 1..mg {
                c = nxt(c);
                fill(p, &cod[c][om as usize], om);
            }
        }
        if r > 0 {
            let om = cap(r);
            for q in 0..r {
                let o = q.min(om);
                fill(mg + clen + q, &bg[o as usize], o);
            }
        }
    }
    for u in 0..n * w as usize {
        if s.mat[u] != NEG {
            s.mat[u] -= b[u];
        }
    }
    s.rc = site_rc(&s.mat, n, w, site_order);
    Ok(())
}

// ---------------------------------------------------------------------------

/// Branch point section (key BP), M+0x5d0..0x7e8. Verified bit-exactly
/// against a gdb dump of fungus.mod (re-notes/bp.md). `limit_duration`
/// (LIMIT_BP_DURATION) and the spacer `fwd`/`rc`/`low` tables are read but
/// never used by the binary.
///
/// The BP intron model splits an intron into donor..branch point (DON_BP
/// length distribution), the branch point site (BRANCH matrix) and a
/// branch point..acceptor spacer (BP_ACC length distribution, own Markov
/// chain MARKOV_BP_SPACER), ending in a BP-specific acceptor (ACC_BP_0..2).
#[allow(dead_code)]
#[derive(Debug, Clone)]
pub struct BranchPoint {
    /// BRANCH: branch point site matrix
    pub branch: Site,
    /// DON_BP: distance donor .. branch point (p4 = p5 = 0)
    pub don_bp: LenDist,
    /// BP_ACC: spacer length branch point .. acceptor (p4 = p5 = 0)
    pub bp_acc: LenDist,
    /// LOG(TO_BP): added to the BP acceptor score
    pub log_to_bp: f64,
    /// LOG(AROUND_BP): added to every normal acceptor site score
    pub log_around_bp: f64,
    /// LIMIT_BP_DURATION (unused)
    pub limit_duration: i32,
    /// BP_SPACER_ORDER: Markov order of the spacer chain
    pub spacer_order: i32,
    /// MARKOV_BP_SPACER (one column, N slot = N_PROB_NON); only `fwd_nf`
    /// and `rc_nf` enter the scores
    pub spacer: NonCoding,
    /// ACC_BP_0..2: acceptor following a branch point, per intron phase
    pub acc_bp: [Site; 3],
}

/// The loaded and transformed model (model object `M`, 0x9f0 bytes).
/// All probabilities are natural logs.
#[allow(dead_code)]
#[derive(Debug, Clone)]
pub struct Model {
    pub file: String,
    /// `$ORGANISM` (unused by the gene finder)
    pub organism: String,
    /// `$ORDER`: Markov order K of the content models (M+0x38)
    pub order: i32,
    /// INTERGENIC_* (p4 = 2, p5 = 1)
    pub intergenic: LenDist,
    /// INTRON_* (p4 = p5 = 0)
    pub intron: LenDist,
    /// SINGLE_* (p4 = p5 = 0)
    pub single: LenDist,
    /// INITIAL_* (p4 = 1, p5 = 0)
    pub initial: LenDist,
    /// EXON_* (internal exons, p4 = 2, p5 = 1); its MAX is also the scan
    /// cut-off for every exon type in the Viterbi
    pub internal: LenDist,
    /// TERMINAL_* (p4 = 1, p5 = 0)
    pub terminal: LenDist,
    /// DONOR_0..2: GT donor per intron phase
    pub donor: [Site; 3],
    /// ACCEPTOR_0..2: AG acceptor per intron phase
    pub acceptor: [Site; 3],
    /// INI: start codon
    pub start: Site,
    /// TAA, TAG, TGA
    pub stop: [Site; 3],
    /// DONOR_GC_0..2, only with GC_DONOR_ON
    pub donor_gc: Option<[Site; 3]>,
    /// branch point block, only with `$BP` set
    pub bp: Option<BranchPoint>,
    /// MARKOV columns 0..2
    pub coding: Coding,
    /// MARKOV column 4 (M+0x930): forward only; intergenic content model
    pub nc4: NonCoding,
    /// MARKOV column 3 (M+0x948): with reverse complement; intron content
    /// model
    pub nc3: NonCoding,
    /// INITIAL_EXON_PHASE: log(P(phase) / (1/3))
    pub phase_initial: [f64; 3],
    /// TERMINAL_EXON_PHASE: log(P(phase) / (1/3))
    pub phase_terminal: [f64; 3],
    /// INTERNAL_EXON_PHASE: log(P(phase pair) / (1/9)), index right + 3*left
    pub phase_internal: [f64; 9],
    /// log ToSingleGene (M+0x990)
    pub log_to_single: f64,
    /// log ToMultiGene (M+0x998)
    pub log_to_multi: f64,
    /// log ToInternalExon (M+0x9a0)
    pub log_to_internal: f64,
    /// log ToTerminalExon (M+0x9a8)
    pub log_to_terminal: f64,
    /// INITIAL_DIRECT_INTRON, INITIAL_REVERSE_INTRON, INITIAL_INTERGENIC,
    /// INITIAL_DIRECT_CODING, INITIAL_REVERSE_CODING (logs)
    pub log_initial: [f64; 5],
    /// `$TAA_ON` (default on)
    pub taa_on: bool,
    /// `$TAG_ON` (default on)
    pub tag_on: bool,
    /// `$TGA_ON` (default on)
    pub tga_on: bool,
    /// `$NORM_SITES`: site matrices were turned into log-odds
    pub norm_sites: bool,
    /// `$GC_DONOR_ON`
    pub gc_on: bool,
    /// log GT_DONOR_PROB (0.0 = log 1 when GC donors are off)
    pub log_gt_donor: f64,
    /// raw 0.0 when GC donors are off
    pub log_gc_donor: f64,
    /// one empty entry per index-list warning (the warnings themselves are
    /// printed immediately); not read anywhere
    pub warnings: Vec<String>,
}

/// 0x411d80 (with 0x411020): normalise a phase table and store it as
/// `log(a[i] / sum) + c`, i.e. log of the probability relative to the
/// uniform distribution (c = ln 3 or ln 9). Zero entries become NEG; an
/// all-zero table is fatal. Sum from 0.0, left to right.
fn phase_table<const N: usize>(a: Vec<f64>, c: f64) -> Result<[f64; N], ModelError> {
    let mut s = 0.0;
    for &x in &a {
        s += x;
    }
    if s == 0.0 {
        return Err(fatal("error in normlog\n".into()));
    }
    let mut out = [0.0; N];
    for i in 0..N {
        out[i] = if a[i] == 0.0 { NEG } else { (a[i] / s).ln() + c };
    }
    Ok(out)
}

impl Model {
    /// 0x412fe0. On error, `ModelError.messages` holds everything the binary
    /// prints; the caller appends `error in model file <F>` unless fatal.
    /// `name` is the file name as given (escaped bytes, see
    /// `sequence::decode`).
    pub fn load(name: &str) -> Result<Model, ModelError> {
        let pb = crate::sequence::os_path(name);
        let path: &Path = &pb;
        let name = name.to_string();
        let raw = match std::fs::read(path) {
            Ok(r) => r,
            Err(_) => {
                let msg = if path.exists() {
                    format!("gmhmme3: error, reading file <{name}>\n")
                } else {
                    format!("gmhmme3: error, file <{name}> is not found\n")
                };
                return Err(fatal(msg));
            }
        };
        let mut kl = KeyList::parse(&raw).ok_or_else(|| fatal(format!("error in file format, {name}\n")))?;
        kl.file = name.clone();
        Self::from_keys(&kl, name)
    }

    /// 0x4123c0: read every key in the original's order. The order matters
    /// because the first failing key decides the error message.
    fn from_keys(kl: &KeyList, file: String) -> Result<Model, ModelError> {
        let mut warn = Vec::new();
        let organism = kl.string("ORGANISM")?;
        let order = kl.int("ORDER")?;
        // length distributions, with the p4/p5 base-pointer offsets
        let single = read_len_dist(kl, "SINGLE", 0, 0, &mut warn)?;
        let internal = read_len_dist(kl, "EXON", 2, 1, &mut warn)?;
        let initial = read_len_dist(kl, "INITIAL", 1, 0, &mut warn)?;
        let terminal = read_len_dist(kl, "TERMINAL", 1, 0, &mut warn)?;
        let intron = read_len_dist(kl, "INTRON", 0, 0, &mut warn)?;
        let intergenic = read_len_dist(kl, "INTERGENIC", 2, 1, &mut warn)?;
        let norm_sites = kl.opt_bool("NORM_SITES", false)?;

        // 0x411c50: MARKOV columns 0..2 coding, 3 intron, 4 intergenic
        let cols = kl.kmer_rows("MARKOV", "ORDER", order, 5)?;
        let n_non = kl.dbl("N_PROB_NON")?;
        let n_cod = kl.dbl("N_PROB_COD")?;
        let coding = build_coding(&cols[0..3], order, n_cod);
        let nc3 = build_noncoding(&cols[3], order, n_non, true);
        let nc4 = build_noncoding(&cols[4], order, n_non, false);

        let start = read_site(kl, "INI")?;
        let taa = read_site(kl, "TERM_TAA")?;
        let tag = read_site(kl, "TERM_TAG")?;
        let tga = read_site(kl, "TERM_TGA")?;
        let d0 = read_site(kl, "DONOR_0")?;
        let d1 = read_site(kl, "DONOR_1")?;
        let d2 = read_site(kl, "DONOR_2")?;
        let a0 = read_site(kl, "ACCEPTOR_0")?;
        let a1 = read_site(kl, "ACCEPTOR_1")?;
        let a2 = read_site(kl, "ACCEPTOR_2")?;

        let mut bp = None;
        if kl.opt_bool("BP", false)? {
            let branch = read_site(kl, "BRANCH")?;
            let bp_acc = read_len_dist(kl, "BP_ACC", 0, 0, &mut warn)?;
            let don_bp = read_len_dist(kl, "DON_BP", 0, 0, &mut warn)?;
            let to_bp = kl.dbl("TO_BP")?;
            let around = kl.dbl("AROUND_BP")?;
            let limit_duration = kl.int("LIMIT_BP_DURATION")?;
            let spacer_order = kl.int("BP_SPACER_ORDER")?;
            let sp = kl.kmer_rows("MARKOV_BP_SPACER", "BP_SPACER_ORDER", spacer_order, 1)?;
            // N_PROB_NON is read a second time here, as in the original
            let n_non2 = kl.dbl("N_PROB_NON")?;
            let spacer = build_noncoding(&sp[0], spacer_order, n_non2, true);
            let acc_bp = [read_site(kl, "ACC_BP_0")?, read_site(kl, "ACC_BP_1")?, read_site(kl, "ACC_BP_2")?];
            // The BP candidates look these distributions up from length
            // MIN-1 on (0x41fad0 / 0x41fec0); with MIN <= 0 that is in front
            // of the table and the original reads heap memory there. The
            // port treats lengths below 1 as impossible (LenDist::base).
            // (MIN near INT_MIN wraps the lookup bound, then nothing is read)
            for (d, key) in [(&bp_acc, "BP_ACC"), (&don_bp, "DON_BP")] {
                if d.min <= 0 && d.lo() < 0 {
                    crate::warn(&format!(
                        "model file {}: {key}_MIN is {}; lengths below 1 are treated as impossible \
                         (not handled by the original)",
                        kl.file, d.min
                    ));
                }
            }
            bp = Some(BranchPoint {
                branch,
                don_bp,
                bp_acc,
                log_to_bp: checked_log(to_bp)?,
                log_around_bp: checked_log(around)?,
                limit_duration,
                spacer_order,
                spacer,
                acc_bp,
            });
        }

        // 0x411d80
        let phase_initial = phase_table::<3>(kl.index_list("INITIAL_EXON_PHASE", 3, &mut warn)?, LN3)?;
        let phase_terminal = phase_table::<3>(kl.index_list("TERMINAL_EXON_PHASE", 3, &mut warn)?, LN3)?;
        let phase_internal = phase_table::<9>(kl.index_list("INTERNAL_EXON_PHASE", 9, &mut warn)?, LN9)?;

        // gene-structure transitions; all four are read before any LOG, and
        // they are not renormalised against each other
        let ts = kl.dbl("ToSingleGene")?;
        let tm = kl.dbl("ToMultiGene")?;
        let ti = kl.dbl("ToInternalExon")?;
        let tt = kl.dbl("ToTerminalExon")?;
        let log_to_single = checked_log(ts)?;
        let log_to_multi = checked_log(tm)?;
        let log_to_internal = checked_log(ti)?;
        let log_to_terminal = checked_log(tt)?;

        // initial state probabilities: main's defaults, unless the file has
        // INITIAL_INTERGENIC, which makes all five keys required
        let mut init = [0.33, 0.33, 0.33, 0.005, 0.005];
        if kl.has("INITIAL_INTERGENIC") {
            init = [
                kl.dbl("INITIAL_DIRECT_INTRON")?,
                kl.dbl("INITIAL_REVERSE_INTRON")?,
                kl.dbl("INITIAL_INTERGENIC")?,
                kl.dbl("INITIAL_DIRECT_CODING")?,
                kl.dbl("INITIAL_REVERSE_CODING")?,
            ];
        }
        let mut log_initial = [0.0; 5];
        for i in 0..5 {
            log_initial[i] = checked_log(init[i])?;
        }

        let taa_on = kl.opt_bool("TAA_ON", true)?;
        let tag_on = kl.opt_bool("TAG_ON", true)?;
        let tga_on = kl.opt_bool("TGA_ON", true)?;
        let gc_on = kl.opt_bool("GC_DONOR_ON", false)?;
        // without GC donors: log GT prob = LOG(1.0) = 0, and the GC slot
        // keeps its raw default 0.0 (not NEG)
        let (donor_gc, log_gt_donor, log_gc_donor) = if gc_on {
            let g = [read_site(kl, "DONOR_GC_0")?, read_site(kl, "DONOR_GC_1")?, read_site(kl, "DONOR_GC_2")?];
            let gt = kl.dbl("GT_DONOR_PROB")?;
            let gc = kl.dbl("GC_DONOR_PROB")?;
            (Some(g), checked_log(gt)?, checked_log(gc)?)
        } else {
            (None, checked_log(1.0)?, 0.0)
        };

        let mut m = Model {
            file,
            organism,
            order,
            intergenic,
            intron,
            single,
            initial,
            internal,
            terminal,
            donor: [d0, d1, d2],
            acceptor: [a0, a1, a2],
            start,
            stop: [taa, tag, tga],
            donor_gc,
            bp,
            coding,
            nc4,
            nc3,
            phase_initial,
            phase_terminal,
            phase_internal,
            log_to_single,
            log_to_multi,
            log_to_internal,
            log_to_terminal,
            log_initial,
            taa_on,
            tag_on,
            tga_on,
            norm_sites,
            gc_on,
            log_gt_donor,
            log_gc_donor,
            warnings: warn,
        };
        if norm_sites {
            m.normalise_sites()?;
        }
        Ok(m)
    }

    /// NORM_SITES: background subtraction for every site, in the original's
    /// order (acceptors, donors, start, stops, GC donors, BP sites). Each
    /// call only touches its own site, so the order does not change values,
    /// but it is kept for traceability.
    fn normalise_sites(&mut self) -> Result<(), ModelError> {
        let k = self.order;
        let cod = self.coding.low.clone();
        let bg3 = self.nc3.low.clone();
        let bg4 = self.nc4.low.clone();
        let f = self.file.clone();
        for ph in 0..3 {
            site_norm(&mut self.acceptor[ph], &bg3, &cod, ph, false, 2, k, &f)?;
        }
        for ph in 0..3 {
            site_norm(&mut self.donor[ph], &bg3, &cod, ph, true, 2, k, &f)?;
        }
        site_norm(&mut self.start, &bg4, &cod, 0, false, 3, k, &f)?;
        for s in self.stop.iter_mut() {
            site_norm(s, &bg4, &cod, 0, true, 3, k, &f)?;
        }
        if let Some(g) = self.donor_gc.as_mut() {
            for (ph, s) in g.iter_mut().enumerate() {
                site_norm(s, &bg3, &cod, ph, true, 2, k, &f)?;
            }
        }
        if let Some(bp) = self.bp.as_mut() {
            site_norm(&mut bp.branch, &bg3, &cod, 0, false, 0, k, &f)?;
            for (ph, s) in bp.acc_bp.iter_mut().enumerate() {
                site_norm(s, &bg3, &cod, ph, false, 2, k, &f)?;
            }
        }
        Ok(())
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn fast_decimal_matches_parse() {
        for t in ["0.0000071028", "1", "-0.5", "0", "-0", "12345.6789", "0.1", "000.000120", "9007199254740992", "1.", ".5", "+3.25", "0.0001204770"] {
            let a = strtod(t);
            let b: f64 = t.parse().unwrap();
            assert_eq!(a.to_bits(), b.to_bits(), "{t}");
        }
        // pseudo-random decimals with up to 15 digits
        let mut x: u64 = 88172645463325252;
        for _ in 0..200000 {
            x ^= x << 13;
            x ^= x >> 7;
            x ^= x << 17;
            let m = x % 1_000_000_000_000_000;
            let f = (x >> 50) % 16;
            let s = format!("{}", m);
            let t = if (f as usize) < s.len() { format!("{}.{}", &s[..s.len() - f as usize], &s[s.len() - f as usize..]) } else { format!("0.{:0>w$}", s, w = f as usize) };
            let a = strtod(&t);
            let b: f64 = t.parse().unwrap();
            assert_eq!(a.to_bits(), b.to_bits(), "{t}");
        }
    }

    #[test]
    fn number_prefix_parsing() {
        assert_eq!(strtol("12abc"), 12);
        assert_eq!(strtol("-7"), -7);
        assert_eq!(strtod("0.5x"), 0.5);
        assert_eq!(strtod("1e-3"), 0.001);
        assert_eq!(strtod("2e"), 2.0);
        assert_eq!(strtod("abc"), 0.0);
    }

    #[test]
    fn keylist_rules() {
        let kl = KeyList::parse(b"  $A 1\n$B \"x\" # c $ in comment\n$C 1 2\n").unwrap();
        assert_eq!(kl.int("A").unwrap(), 1);
        assert_eq!(kl.string("B").unwrap(), "x");
        assert!(kl.int("C").is_err());
        assert!(KeyList::parse(b"$A 1$B 2").is_none());
        assert!(KeyList::parse(b"x $A 1").is_none());
        let dup = KeyList::parse(b"$A 1\n$A 2\n").unwrap();
        assert!(dup.int("A").is_err());
    }

    #[test]
    fn rc_index_roundtrip() {
        // AC (0*4+1) -> revcomp GT (2*4+3)
        assert_eq!(rc_index(1, 2), 11);
    }
}
