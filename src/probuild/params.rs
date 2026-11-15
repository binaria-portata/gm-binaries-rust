//! Parameter object of probuild 2.52: the `std::map<string,string>` of all
//! command line options, its parsing from argv and from a `--par` file, and
//! the typed getters the mode handlers use (re-notes/probuild/core.md
//! sections 2 and 3).
//!
//! Every option is a key starting with `--`, optionally followed by one
//! value.  Keys are stored with their dashes (`"--seq"`), values as raw bytes;
//! a key given without a value maps to the empty string.  Keys compare
//! bytewise, so a `BTreeMap` iterates in the same order as the C++ map.
//!
//! The getters (`req_str`, `opt_int`, `opt_bool`, ...) validate a value by
//! its character set, print the original's messages to stderr on failure and
//! convert with C `atoi` / `atof` semantics (reproduced here as `atoi`,
//! `strtol`, `atof`).  They return `false` on failure; handlers chain them
//! with `&&`, so the first failing getter stops the handler and later options
//! are never checked.  The model modes use the copy-on-read variant in
//! `pm.rs` instead.

use super::io::{err, load_file};
use std::collections::BTreeMap;

/// Byte string (C++ `std::string`): file names, keys and values are not
/// required to be UTF-8.
pub type Bytes = Vec<u8>;

/// Delimiters of the key test on argv (0x5655ae): a key may not be followed
/// by a blank, so `"-- x"` given as one argument is not a key.
const ARG_DELIMS: &[u8] = b" \t";
/// Delimiters of the `--par` file format and of the "value contains white
/// space" check of string getters (0x5655f3).
const PAR_DELIMS: &[u8] = b" \t\n\r";

/// Characters allowed by the signed integer getter (0x41ff1a).
const SIGNED_INT_CHARS: &[u8] = b"-0123456789";
/// Characters allowed by the unsigned integer getters (0x41fbd6).
const UNSIGNED_INT_CHARS: &[u8] = b"0123456789";
/// Characters allowed by the float and bool getters (0x41fc38).
const FLOAT_CHARS: &[u8] = b"+-eE.0123456789";

/// The parameter object: option key (with leading `--`) to raw value.
#[derive(Default, Clone)]
pub struct Params {
    /// All options; argv entries override `--par` file entries.
    pub map: BTreeMap<Bytes, Bytes>,
}

/// `s.find(pat, from)`: first occurrence of `pat` at or after `from`.
fn find_from(s: &[u8], pat: &[u8], from: usize) -> Option<usize> {
    if from > s.len() {
        return None;
    }
    if pat.is_empty() {
        return Some(from);
    }
    s[from..].windows(pat.len()).position(|w| w == pat).map(|p| p + from)
}
/// `s.find_first_of(set, from)`.
fn first_of(s: &[u8], set: &[u8], from: usize) -> Option<usize> {
    (from..s.len()).find(|&i| set.contains(&s[i]))
}
/// `s.find_first_not_of(set, from)`.
fn first_not_of(s: &[u8], set: &[u8], from: usize) -> Option<usize> {
    (from..s.len()).find(|&i| !set.contains(&s[i]))
}
/// `find_last_not_of(set, pos)`: search backwards from `min(pos, len-1)`.
fn last_not_of(s: &[u8], set: &[u8], pos: usize) -> Option<usize> {
    if s.is_empty() {
        return None;
    }
    let mut i = pos.min(s.len() - 1);
    loop {
        if !set.contains(&s[i]) {
            return Some(i);
        }
        if i == 0 {
            return None;
        }
        i -= 1;
    }
}

/// Key test 0x41e83c: does a key (the prefix `--`) start at `pos` of `s`?
///
/// Besides the prefix itself, the byte before it (if any) must be a
/// delimiter, and the byte after it (if any) must NOT be a delimiter.  So
/// `--` alone is a key (the C++ code accepts it), `-- x` is not, and a `--`
/// inside a word (`a--b`) does not start a key.
fn is_key(s: &[u8], pos: usize, delims: &[u8]) -> bool {
    if pos > s.len() || s.len() - pos < 2 || &s[pos..pos + 2] != b"--" {
        return false;
    }
    if pos != 0 && !delims.contains(&s[pos - 1]) {
        return false;
    }
    if pos + 2 < s.len() && delims.contains(&s[pos + 2]) {
        return false;
    }
    true
}

impl Params {
    /// Presence test 0x405988 (`map.find(key) != end`).  Flags such as
    /// `--allow_x` are often tested this way, so their value is irrelevant.
    pub fn has(&self, k: &str) -> bool {
        self.map.contains_key(k.as_bytes())
    }

    /// Value of `k`, or empty if absent.  Read-only: the original's
    /// `map[key]` inserts a missing key, which only matters for later
    /// presence tests; the model modes, where it does matter, use `pm::P`.
    pub fn val(&self, k: &str) -> Bytes {
        self.map.get(k.as_bytes()).cloned().unwrap_or_default()
    }

    /// 0x41f6c4: parse argv (skipping the program name), then, if `--par` is
    /// present, the parameter file named by it.  The file is parsed into an
    /// empty map and the argv entries are then written over it, so the
    /// command line wins over the file.  `None` = parse failure (the message
    /// is already printed; main exits with code 1).
    pub fn parse(args: &[Bytes]) -> Option<Params> {
        let mut p = Params::default();
        if !p.parse_argv(args) {
            return None;
        }
        if p.has("--par") {
            let mut par_file = Bytes::new();
            if !p.opt_str(&mut par_file, "--par") {
                return None;
            }
            let argv_entries = p.map.clone();
            p.map.clear();
            if !p.parse_par_file(&par_file) {
                return None;
            }
            for (k, v) in argv_entries {
                p.map.insert(k, v);
            }
        }
        Some(p)
    }

    /// argv parser 0x41f000.  Each argument is either a key (see `is_key`) or
    /// the value of the preceding key:
    ///
    /// - no pending key: the argument must be a key, else
    ///   `error in command line` / `error, key expected, this found: ARG`;
    /// - pending key and the argument is a key: the pending key stays with
    ///   an empty value and the new key becomes pending;
    /// - pending key and the argument is not a key: it is the value.
    ///
    /// Consequences: values may start with `-` (`--min_contig -5`), a
    /// repeated key keeps the last value, and `--a v w` fails at `w`.
    fn parse_argv(&mut self, args: &[Bytes]) -> bool {
        // Key waiting for its value; empty = none pending.  (An argument that
        // is itself empty can never be a key, so emptiness is a safe marker.)
        let mut pending_key: Bytes = Vec::new();
        for a in args.iter().skip(1) {
            if pending_key.is_empty() {
                pending_key = a.clone();
                if is_key(&pending_key, 0, ARG_DELIMS) {
                    self.map.insert(pending_key.clone(), Vec::new());
                } else {
                    err(b"error in command line\n");
                    let mut m = b"error, key expected, this found: ".to_vec();
                    m.extend_from_slice(&pending_key);
                    m.push(b'\n');
                    err(&m);
                    return false;
                }
            } else if is_key(a, 0, ARG_DELIMS) {
                pending_key = a.clone();
                self.map.insert(pending_key.clone(), Vec::new());
            } else {
                self.map.insert(pending_key.clone(), a.clone());
                pending_key.clear();
            }
        }
        true
    }

    /// `--par` file reader 0x41f428: load the whole file (0x40a306), blank
    /// out comments (0x41ef00), then parse the text (0x41e978).  A format
    /// error prints `error in format of file NAME`.
    fn parse_par_file(&mut self, name: &[u8]) -> bool {
        let mut text = match load_file(name) {
            Some(t) => t,
            None => return false,
        };
        // Comments: from '#' up to (not including) the next '\n' or '\r'
        // every byte becomes a blank.  There is no quoting, so a '#' inside a
        // quoted value also starts a comment.
        let mut search_from = 0;
        while let Some(hash) = (search_from..text.len()).find(|&j| text[j] == b'#') {
            let comment_end = first_of(&text, b"\n\r", hash).unwrap_or(text.len());
            for c in &mut text[hash..comment_end] {
                *c = b' ';
            }
            search_from = comment_end;
        }
        if !self.parse_par_text(&text) {
            let mut m = b"error in format of file ".to_vec();
            m.extend_from_slice(name);
            m.push(b'\n');
            err(&m);
            return false;
        }
        true
    }

    /// Key/value text parser 0x41e978 (also used by `--combine` to read model
    /// files, see `mkmod.rs`).  The text is a sequence of `--KEY [VALUE]` items
    /// separated by white space (`" \t\n\r"`); line breaks have no special
    /// meaning.  For each key:
    ///
    /// - the key runs up to the next white space;
    /// - the value is everything between the key and the next key (the next
    ///   `--` that passes `is_key`), with surrounding white space trimmed;
    /// - a value longer than one byte enclosed in `"..."` or `'...'` loses
    ///   the quotes and may contain white space; an unquoted value containing
    ///   white space (i.e. two words) is a format error.
    ///
    /// Returns false at the first item not starting with a key.  Entries
    /// parsed before a failure stay in the map.
    pub fn parse_par_text(&mut self, s: &[u8]) -> bool {
        let delims = PAR_DELIMS;
        let mut item_start = first_not_of(s, delims, 0);
        while let Some(key_start) = item_start {
            if !is_key(s, key_start, delims) {
                return false;
            }
            let key_end = first_of(s, delims, key_start).unwrap_or(s.len());
            let key = s[key_start..key_end].to_vec();
            // Find the next key: the next "--" that passes the key test.  A
            // rejected hit continues the search one byte further.
            let mut search_from = key_end;
            let mut next_key = None;
            while let Some(hit) = find_from(s, b"--", search_from) {
                if is_key(s, hit, delims) {
                    next_key = Some(hit);
                    break;
                }
                search_from = hit + 1;
            }
            // The value ends at the last non-blank before the next key (or
            // before the end of the text).  `next_key - 1` is the delimiter in
            // front of it, guaranteed by the key test.  When there is no
            // non-blank at all, the original computes `npos + 1 == 0`.
            let value_limit = match next_key {
                Some(f) => f - 1,
                None => s.len(),
            };
            let mut value_end = match last_not_of(s, delims, value_limit) {
                Some(x) => x + 1,
                None => 0,
            };
            let value_start = first_not_of(s, delims, key_end).unwrap_or(s.len());
            // Key without value: the backward search found the key itself.
            if value_end < value_start {
                value_end = value_start;
            }
            let mut v = s[value_start..value_end].to_vec();
            if v.len() > 1 {
                let (first, last) = (v[0], v[v.len() - 1]);
                if (first == b'"' && last == b'"') || (first == b'\'' && last == b'\'') {
                    v = v[1..v.len() - 1].to_vec();
                } else if first_of(&v, delims, 0).is_some() {
                    return false;
                }
            }
            self.map.insert(key, v);
            item_start = next_key;
        }
        true
    }

    // ---- error messages ------------------------------------------------

    /// 0x41fa98: `error in value of parameter KEY: VALUE` on stderr.
    fn bad_value(&self, k: &str) {
        let mut m = b"error in value of parameter ".to_vec();
        m.extend_from_slice(k.as_bytes());
        m.extend_from_slice(b": ");
        m.extend_from_slice(&self.val(k));
        m.push(b'\n');
        err(&m);
    }

    /// 0x41fb36: `error: missing parameter KEY` on stderr.
    pub fn missing(&self, k: &str) {
        err(format!("error: missing parameter {k}\n").as_bytes());
    }

    /// Character check: every byte of the value is in `set` (an empty value
    /// passes).  Reports a bad value.
    fn only(&self, k: &str, set: &[u8]) -> bool {
        let ok = first_not_of(&self.val(k), set, 0).is_none();
        if !ok {
            self.bad_value(k);
        }
        ok
    }
    /// 0x41fd58: value is a "word": non-empty and without white space.
    fn check_word(&self, k: &str) -> bool {
        let v = self.val(k);
        let ok = !v.is_empty() && first_of(&v, PAR_DELIMS, 0).is_none();
        if !ok {
            self.bad_value(k);
        }
        ok
    }
    /// 0x41fd10: calls the float character check 0x41fc38 (which reports the
    /// error) and then reports the same error a second time.  The duplicated
    /// stderr line is observable and kept.
    fn check_float(&self, k: &str) -> bool {
        let ok = self.only(k, FLOAT_CHARS);
        if !ok {
            self.bad_value(k);
        }
        ok
    }

    // ---- getters ---------------------------------------------------------
    //
    // Common contract: `out` is only written on success with a present key,
    // so it keeps the caller's default otherwise.  The return value is false
    // on a missing required key or a bad value (message already printed).

    /// 0x41fde6: optional string, any value (may be empty or contain
    /// blanks).
    pub fn opt_raw(&self, out: &mut Bytes, k: &str) -> bool {
        if let Some(v) = self.map.get(k.as_bytes()) {
            *out = v.clone();
        }
        true
    }
    /// 0x41fe70: optional word (a present value must be non-empty and free
    /// of white space).
    pub fn opt_str(&self, out: &mut Bytes, k: &str) -> bool {
        if self.has(k) {
            if !self.check_word(k) {
                return false;
            }
            *out = self.val(k);
        }
        true
    }
    /// 0x42034c: required string, any value.
    pub fn req_raw(&self, out: &mut Bytes, k: &str) -> bool {
        if self.has(k) {
            *out = self.val(k);
            true
        } else {
            self.missing(k);
            false
        }
    }
    /// 0x4203f4: required word.
    pub fn req_str(&self, out: &mut Bytes, k: &str) -> bool {
        if self.has(k) {
            if !self.check_word(k) {
                return false;
            }
            *out = self.val(k);
            true
        } else {
            self.missing(k);
            false
        }
    }
    /// 0x41ff1a: optional signed int: only `-0123456789`, then `atoi`.  The
    /// check does not validate the position of `-` (`5-3` gives 5) and an
    /// empty value passes (`atoi("") == 0`).
    pub fn opt_int(&self, out: &mut i32, k: &str) -> bool {
        if self.has(k) {
            if !self.only(k, SIGNED_INT_CHARS) {
                return false;
            }
            *out = atoi(&self.val(k));
        }
        true
    }
    /// 0x41ff86: optional "unsigned" int: only digits, then `atoi`.  The
    /// result is still a signed 32-bit int: a long digit string wraps (it can
    /// become negative), exactly like the original.
    pub fn opt_uint(&self, out: &mut i32, k: &str) -> bool {
        if self.has(k) {
            if !self.only(k, UNSIGNED_INT_CHARS) {
                return false;
            }
            *out = atoi(&self.val(k));
        }
        true
    }
    /// 0x4204ba: required "unsigned" int (as `opt_uint`).
    pub fn req_uint(&self, out: &mut i32, k: &str) -> bool {
        if self.has(k) {
            if !self.only(k, UNSIGNED_INT_CHARS) {
                return false;
            }
            *out = atoi(&self.val(k));
            true
        } else {
            self.missing(k);
            false
        }
    }
    #[allow(dead_code)]
    /// 0x41fff2: optional double: only `+-eE.0123456789`, then `atof`.
    pub fn opt_f64(&self, out: &mut f64, k: &str) -> bool {
        if self.has(k) {
            if !self.only(k, FLOAT_CHARS) {
                return false;
            }
            *out = atof(&self.val(k));
        }
        true
    }
    /// 0x4200de: optional bool.  A present key with an empty value (a plain
    /// flag) is true; otherwise the value must pass the float character
    /// check (reported twice on failure, see `check_float`) and is true if
    /// `atof(value) != 0` (0x4206b2), so `--soft 0` is false.
    pub fn opt_bool(&self, out: &mut bool, k: &str) -> bool {
        if self.has(k) {
            if self.val(k).is_empty() {
                *out = true;
                return true;
            }
            if !self.check_float(k) {
                return false;
            }
            *out = atof(&self.val(k)) != 0.0;
        }
        true
    }
    #[allow(dead_code)]
    /// 0x420540: required bool.  Unlike `opt_bool` there is no empty-value
    /// shortcut: an empty value passes the character check and gives
    /// `atof("") == 0`, i.e. false.
    pub fn req_bool(&self, out: &mut bool, k: &str) -> bool {
        if self.has(k) {
            if !self.check_float(k) {
                return false;
            }
            *out = atof(&self.val(k)) != 0.0;
            true
        } else {
            self.missing(k);
            false
        }
    }
    /// 0x4201ae: required non-empty list; items are separated by runs of any
    /// byte of `delims` (empty items vanish) and appended to `out`.
    /// `--reformat_fasta --in` uses `"\t"` only, so file names may contain
    /// blanks.  A value consisting of delimiters only gives an empty list
    /// and still succeeds.
    pub fn req_list(&self, out: &mut Vec<Bytes>, delims: &[u8], k: &str) -> bool {
        if !self.has(k) {
            self.missing(k);
            return false;
        }
        let v = self.val(k);
        if v.is_empty() {
            self.bad_value(k);
            return false;
        }
        let mut item_start = first_not_of(&v, delims, 0);
        while let Some(b) = item_start {
            let e = first_of(&v, delims, b).unwrap_or(v.len());
            out.push(v[b..e].to_vec());
            item_start = first_not_of(&v, delims, e);
        }
        true
    }
}

/// glibc `atoi` on LP64 = `(int) strtol(s, NULL, 10)`: the saturated 64-bit
/// result is truncated to its low 32 bits, so `"4294967295"` gives -1 and
/// `"99999999999"` gives 1215752191.
pub fn atoi(s: &[u8]) -> i32 {
    strtol(s) as i32
}

/// glibc `strtol(s, NULL, 10)`: optional leading C white space, optional
/// sign, then decimal digits up to the first non-digit (no digits: 0).
/// Out-of-range values saturate at `i64::MIN` / `i64::MAX`.
pub fn strtol(b: &[u8]) -> i64 {
    let mut i = 0;
    while i < b.len() && matches!(b[i], b' ' | b'\t' | b'\n' | 0x0b | 0x0c | b'\r') {
        i += 1;
    }
    let mut neg = false;
    if i < b.len() && (b[i] == b'+' || b[i] == b'-') {
        neg = b[i] == b'-';
        i += 1;
    }
    // Accumulate in i128.  Once the magnitude exceeds |i64::MIN| it is pinned
    // to the sentinel i64::MAX + 2, which saturates in either direction
    // below and keeps the i128 from overflowing on very long digit strings.
    let mut v: i128 = 0;
    while i < b.len() && b[i].is_ascii_digit() {
        v = v * 10 + (b[i] - b'0') as i128;
        if v > (i64::MAX as i128) + 1 {
            v = (i64::MAX as i128) + 2;
        }
        i += 1;
    }
    if neg {
        let n = -v;
        if n < i64::MIN as i128 {
            i64::MIN
        } else {
            n as i64
        }
    } else if v > i64::MAX as i128 {
        i64::MAX
    } else {
        v as i64
    }
}

/// glibc `atof` = `strtod(s, NULL)`, restricted to plain decimal syntax:
/// optional leading white space and sign, digits with an optional '.', and
/// an optional exponent that is only consumed if it has at least one digit
/// (`"1e"` is 1).  No digits at all gives 0.  The longest such prefix is
/// converted with Rust's correctly rounded parser, which agrees with glibc's
/// correctly rounded `strtod`.
///
/// Not covered: the `inf`, `nan` and hexadecimal forms `strtod` also accepts.
/// The getters' character checks exclude them; the unchecked callers in
/// `pm.rs` would give 0 for such text where glibc would not.
pub fn atof(b: &[u8]) -> f64 {
    let s = String::from_utf8_lossy(b);
    let t = s.trim_start_matches([' ', '\t', '\n', '\x0b', '\x0c', '\r']);
    let tb = t.as_bytes();
    let mut i = 0;
    if i < tb.len() && (tb[i] == b'+' || tb[i] == b'-') {
        i += 1;
    }
    let digits_start = i;
    while i < tb.len() && tb[i].is_ascii_digit() {
        i += 1;
    }
    // Number of mantissa digits before and after the '.'.
    let mut n_digits = i - digits_start;
    if i < tb.len() && tb[i] == b'.' {
        i += 1;
        let frac_start = i;
        while i < tb.len() && tb[i].is_ascii_digit() {
            i += 1;
        }
        n_digits += i - frac_start;
    }
    if n_digits == 0 {
        return 0.0;
    }
    let mut end = i;
    if i < tb.len() && (tb[i] == b'e' || tb[i] == b'E') {
        let mut j = i + 1;
        if j < tb.len() && (tb[j] == b'+' || tb[j] == b'-') {
            j += 1;
        }
        let exp_digits_start = j;
        while j < tb.len() && tb[j].is_ascii_digit() {
            j += 1;
        }
        if j > exp_digits_start {
            end = j;
        }
    }
    // Rust's f64 parser rejects a mantissa that starts or ends with '.', so
    // "1." becomes "1.0" and ".5" / "-.5" become "0.5" / "-0.5" (same value).
    let mut num = t[..end].to_string();
    if num.ends_with('.') {
        num.push('0');
    }
    if num.starts_with('.') || num.starts_with("+.") || num.starts_with("-.") {
        num = num.replacen('.', "0.", 1);
    }
    num.parse::<f64>().unwrap_or(0.0)
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn keys() {
        assert!(is_key(b"--a", 0, ARG_DELIMS));
        assert!(is_key(b"--", 0, ARG_DELIMS));
        assert!(!is_key(b"-- a", 0, ARG_DELIMS));
        assert!(!is_key(b"-a", 0, ARG_DELIMS));
    }
    #[test]
    fn nums() {
        assert_eq!(atoi(b"12"), 12);
        assert_eq!(atoi(b"99999999999"), 1215752191);
        assert_eq!(atoi(b"99999999999999999999"), -1);
        assert_eq!(atof(b"1e3"), 1000.0);
        assert_eq!(atof(b".5"), 0.5);
        assert_eq!(atof(b"-"), 0.0);
        assert_eq!(atof(b"1.e"), 1.0);
    }
    #[test]
    fn par() {
        let mut p = Params::default();
        assert!(p.parse_par_text(b"--A  1\n--B \"x y\"\n--C\n--D 'q'"));
        assert_eq!(p.val("--A"), b"1");
        assert_eq!(p.val("--B"), b"x y");
        assert_eq!(p.val("--C"), b"");
        assert_eq!(p.val("--D"), b"q");
        let mut p = Params::default();
        assert!(!p.parse_par_text(b"--A 1 2"));
    }
}
