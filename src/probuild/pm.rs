//! Parameter access and number formatting for the model modes
//! (re-notes/probuild/models.md section 0.1; getter rules as in core.md
//! section 3).
//!
//! The model handlers work on a mutable copy of the parameter map (`P`)
//! because the original reads values through `std::map::operator[]`
//! (0x413238), which inserts missing keys: after e.g. a failed string check
//! of a missing key, that key is "present" for later `has` tests.  `P::val`
//! reproduces this insertion.  Unlike the `Params` getters, these return the
//! value instead of writing through an out parameter, and some checks are
//! exposed separately because the model code calls them on their own.  All
//! messages go to stderr.
//!
//! Number formatting: the original prints doubles through libstdc++
//! `ostream`, which formats via the C library's printf.  `cfmt` calls the
//! host libc's `snprintf` for the same format string, so `%f`, `%g` and
//! their rounding are produced by glibc itself; `fixed` and `g` are the two
//! `ostream` modes used.  `powi` reproduces the compiler's integer power
//! helper with its exact multiplication order.

use super::io;
use super::params::{atof, atoi, Bytes, Params};
use std::collections::BTreeMap;

/// White space rejected by the string check 0x41fd58.
const WS: &[u8] = b" \t\n\r";

/// Parameter map of a model handler (copy of `Params::map`, mutated by
/// `val` exactly where the original's `operator[]` inserts).
pub struct P {
    /// Option key (with `--`) to raw value.
    pub m: BTreeMap<Bytes, Bytes>,
}

impl P {
    /// Copy of the command line parameters.
    pub fn new(p: &Params) -> P {
        P { m: p.map.clone() }
    }
    /// Wrap an existing map (e.g. a model file parsed with
    /// `Params::parse_par_text`).
    pub fn from_map(m: BTreeMap<Bytes, Bytes>) -> P {
        P { m }
    }
    /// Presence test 0x405988 (no insertion).
    pub fn has(&self, k: &str) -> bool {
        self.m.contains_key(k.as_bytes())
    }
    /// `map[k]`: the value, inserting an empty value if `k` is missing.
    pub fn val(&mut self, k: &str) -> Bytes {
        self.m.entry(k.as_bytes().to_vec()).or_default().clone()
    }
    /// 0x41fa98: `error in value of parameter KEY: VALUE` (inserts `k`).
    pub fn bad(&mut self, k: &str) {
        let mut m = b"error in value of parameter ".to_vec();
        m.extend_from_slice(k.as_bytes());
        m.extend_from_slice(b": ");
        m.extend_from_slice(&self.val(k));
        m.push(b'\n');
        io::err(&m);
    }
    /// 0x41fb36: `error: missing parameter KEY`.
    pub fn missing(&self, k: &str) {
        io::err(format!("error: missing parameter {k}\n").as_bytes());
    }
    /// Every byte of the value is in `set` (an empty value passes); reports
    /// a bad value once.
    fn only(&mut self, k: &str, set: &[u8]) -> bool {
        let ok = self.val(k).iter().all(|c| set.contains(c));
        if !ok {
            self.bad(k);
        }
        ok
    }
    /// 0x41fbd6: digits only (empty passes).
    pub fn chk_int(&mut self, k: &str) -> bool {
        self.only(k, b"0123456789")
    }
    /// 0x41fd58: non-empty, no white space.  On a missing key the inserted
    /// empty value fails, giving `error in value of parameter KEY: `.
    pub fn chk_str(&mut self, k: &str) -> bool {
        let v = self.val(k);
        let ok = !v.is_empty() && !v.iter().any(|c| WS.contains(c));
        if !ok {
            self.bad(k);
        }
        ok
    }
    /// 0x41fc38: float characters `+-eE.0123456789` (one message).
    pub fn chk_float(&mut self, k: &str) -> bool {
        self.only(k, b"+-eE.0123456789")
    }
    /// 0x41fd10: float characters; on failure the message is printed twice
    /// (once by 0x41fc38, once by 0x41fd10 itself).
    pub fn chk_float2(&mut self, k: &str) -> bool {
        if !self.chk_float(k) {
            self.bad(k);
            return false;
        }
        true
    }
    /// 0x41fc9a: float characters and value >= 0.  Wrong characters print
    /// the message twice, a negative value once.
    pub fn chk_float_nonneg(&mut self, k: &str) -> bool {
        let ok = self.chk_float(k) && atof(&self.val(k)) >= 0.0;
        if !ok {
            self.bad(k);
        }
        ok
    }
    /// 0x4204ba: required int (digits only), converted with `atoi`.
    pub fn int_req(&mut self, k: &str) -> Option<i32> {
        if !self.has(k) {
            self.missing(k);
            return None;
        }
        if !self.chk_int(k) {
            return None;
        }
        Some(atoi(&self.val(k)))
    }
    /// 0x41ff86 (`signed == false`, digits) / 0x41ff1a (`signed == true`,
    /// `-` allowed): optional int.  Returns (ok, value); the value is `cur`
    /// when the key is absent or bad.
    pub fn int_opt(&mut self, k: &str, cur: i32, signed: bool) -> (bool, i32) {
        if !self.has(k) {
            return (true, cur);
        }
        let set: &[u8] = if signed { b"-0123456789" } else { b"0123456789" };
        if !self.only(k, set) {
            return (false, cur);
        }
        (true, atoi(&self.val(k)))
    }
    /// 0x4203f4: required word (non-empty, no white space).
    pub fn str_req(&mut self, k: &str) -> Option<Bytes> {
        if !self.has(k) {
            self.missing(k);
            return None;
        }
        if !self.chk_str(k) {
            return None;
        }
        Some(self.val(k))
    }
    /// 0x41fe70: optional word.  Some(value, or `cur` if absent), None on a
    /// bad value.
    pub fn str_opt(&mut self, k: &str, cur: Bytes) -> Option<Bytes> {
        if !self.has(k) {
            return Some(cur);
        }
        if !self.chk_str(k) {
            return None;
        }
        Some(self.val(k))
    }
    /// 0x4200de: optional bool.  Empty value (plain flag) = true; else float
    /// characters (message twice on failure) and `atof != 0`.  Returns
    /// (ok, value).
    pub fn bool_opt(&mut self, k: &str, cur: bool) -> (bool, bool) {
        if !self.has(k) {
            return (true, cur);
        }
        if self.val(k).is_empty() {
            return (true, true);
        }
        if !self.chk_float2(k) {
            return (false, cur);
        }
        (true, atof(&self.val(k)) != 0.0)
    }
    /// 0x41fff2: optional double (float characters, message once).  Returns
    /// (ok, value); a bad value keeps `cur`.
    pub fn f64_opt(&mut self, k: &str, cur: f64) -> (bool, f64) {
        if !self.has(k) {
            return (true, cur);
        }
        if !self.chk_float(k) {
            return (false, cur);
        }
        (true, atof(&self.val(k)))
    }
    /// 0x4206b2: `strtod(value) != 0`, without any character check (inserts
    /// a missing key, which then reads as false).
    pub fn bool_val(&mut self, k: &str) -> bool {
        atof(&self.val(k)) != 0.0
    }
    /// `strtod(value)` without any check (inserts a missing key).
    pub fn f64_val(&mut self, k: &str) -> f64 {
        atof(&self.val(k))
    }
}

/// Format one double with the C format `fmt` (e.g. `"%.4f"`) through the
/// host libc's `snprintf`.  glibc's printf rounding is exact (it rounds the
/// exact binary value), so the result is identical to the original's static
/// glibc 2.17 for `%f` / `%g`.  The 512-byte first attempt covers all
/// normal cases; longer results (huge `%f` values) are formatted again into
/// a buffer of the reported size.
pub fn cfmt(fmt: &str, x: f64) -> String {
    let f = std::ffi::CString::new(fmt).unwrap();
    let mut buf = vec![0u8; 512];
    let n = unsafe { libc::snprintf(buf.as_mut_ptr() as *mut libc::c_char, buf.len(), f.as_ptr(), x) };
    if n < 0 {
        // The result would be longer than INT_MAX characters (a precision
        // option near 2^31). The original formats into a stack buffer of that
        // size and ends with SIGSEGV.
        io::fatal(&format!("cannot format a number with \"{fmt}\": the precision is too large"));
    }
    let n = n as usize;
    if n >= buf.len() {
        buf = vec![0u8; n + 1];
        unsafe { libc::snprintf(buf.as_mut_ptr() as *mut libc::c_char, buf.len(), f.as_ptr(), x) };
    }
    buf.truncate(n);
    String::from_utf8(buf).unwrap()
}

/// `ostream << fixed << setprecision(p) << x`, printf `%.pf`; with
/// `showpoint` it is `%#.pf` (differs only for p = 0: `1.` instead of `1`).
/// A negative precision is treated as 6, as libstdc++ `_M_insert_float`
/// does.
pub fn fixed(p: i32, x: f64, showpoint: bool) -> String {
    let p = if p < 0 { 6 } else { p };
    cfmt(&format!("%{}.{}f", if showpoint { "#" } else { "" }, p), x)
}

/// `ostream << x` with default flags and precision: printf `%g`.
pub fn g(x: f64) -> String {
    cfmt("%g", x)
}

/// libgcc `__powidf2` (the helper behind GCC's `__builtin_powi`):
/// square-and-multiply over the bits of `|n|`, then one reciprocal for
/// negative `n`.  The order of the multiplications decides the last bits of
/// the result, so it is written out here instead of relying on how Rust
/// lowers `f64::powi`.
pub fn powi(mut x: f64, n: i32) -> f64 {
    let mut m = n.unsigned_abs();
    let mut y = if m % 2 == 1 { x } else { 1.0 };
    loop {
        m >>= 1;
        if m == 0 {
            break;
        }
        x *= x;
        if m % 2 == 1 {
            y *= x;
        }
    }
    if n < 0 {
        1.0 / y
    } else {
        y
    }
}
