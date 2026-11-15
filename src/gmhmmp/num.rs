//! C number parsing (glibc strtod / strtol prefix semantics) and printf %f,
//! as used by the binary for model files, options, defline and GFF values:
//! the longest valid prefix is converted, text after it is ignored, and no
//! valid prefix gives 0.

/// C `isspace` in the C locale.
#[inline]
pub fn is_c_space(b: u8) -> bool {
    matches!(b, b' ' | b'\t' | b'\n' | 0x0b | 0x0c | b'\r')
}

/// `strtol(s, NULL, 10)` of glibc: the magnitude is accumulated without
/// sign; if it exceeds LONG_MAX (or 2^63 for a '-' sign) the result
/// saturates to LONG_MAX / LONG_MIN (all digits are still consumed).
pub fn strtol_i64(s: &str) -> i64 {
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
    let mut v: u64 = 0;
    let mut over = false;
    while i < b.len() && b[i].is_ascii_digit() {
        if !over {
            match v.checked_mul(10).and_then(|x| x.checked_add((b[i] - b'0') as u64)) {
                Some(x) => v = x,
                None => over = true,
            }
        }
        i += 1;
    }
    if neg {
        if over || v > 1u64 << 63 {
            i64::MIN
        } else {
            (v as i64).wrapping_neg()
        }
    } else if over || v > i64::MAX as u64 {
        i64::MAX
    } else {
        v as i64
    }
}

/// `(int) strtol(s, NULL, 10)`: glibc strtol (saturating at LONG_MIN /
/// LONG_MAX, see [`strtol_i64`]), then the conversion to int keeps the low
/// 32 bits, as the binary does (`movl` of the long result). Verified on
/// MINC values: 500000000000000000 gives -743309312, 2^62 gives 0,
/// LONG_MAX and above give -1, LONG_MIN and below give 0.
pub fn strtol_i32(s: &str) -> i32 {
    strtol_i64(s) as i32
}

/// `strtod(s, NULL)`: longest valid prefix (decimal, hex, inf, nan).
pub fn strtod(s: &str) -> f64 {
    let t = s.trim_start_matches(|c: char| c.is_ascii() && is_c_space(c as u8));
    // fast path: the whole token parses as a Rust f64. Tokens of 9 or
    // more characters that start with a sign, 'i' or 'n' still go through
    // the prefix scanner below.
    if let Ok(v) = t.parse::<f64>() {
        if !t.starts_with(['i', 'I', 'n', 'N', '+', '-']) || t.len() < 9 {
            return v;
        }
    }
    let b = t.as_bytes();
    let mut i = 0;
    let neg = !b.is_empty() && b[0] == b'-';
    if i < b.len() && (b[i] == b'+' || b[i] == b'-') {
        i += 1;
    }
    let rest = &t[i..];
    let lower = rest.to_ascii_lowercase();
    if lower.starts_with("inf") {
        return if neg { f64::NEG_INFINITY } else { f64::INFINITY };
    }
    if lower.starts_with("nan") {
        return if neg { -f64::NAN } else { f64::NAN };
    }
    if lower.starts_with("0x") {
        // hexadecimal float
        let h = &lower[2..];
        let mut mant: f64 = 0.0;
        let mut exp: i32 = 0;
        let mut seen = false;
        let mut frac = false;
        let mut j = 0;
        let hb = h.as_bytes();
        while j < hb.len() {
            let c = hb[j];
            if c == b'.' && !frac {
                frac = true;
            } else if let Some(d) = (c as char).to_digit(16) {
                mant = mant * 16.0 + d as f64;
                if frac {
                    exp -= 4;
                }
                seen = true;
            } else {
                break;
            }
            j += 1;
        }
        if !seen {
            return 0.0;
        }
        if j < hb.len() && hb[j] == b'p' {
            let e = strtol_i32(&h[j + 1..]);
            exp += e;
        }
        let v = mant * (2f64).powi(exp);
        return if neg { -v } else { v };
    }
    let ds = i;
    while i < b.len() && b[i].is_ascii_digit() {
        i += 1;
    }
    if i < b.len() && b[i] == b'.' {
        i += 1;
        while i < b.len() && b[i].is_ascii_digit() {
            i += 1;
        }
    }
    if i == ds || (i == ds + 1 && b[ds] == b'.') {
        return 0.0;
    }
    let mut end = i;
    if i < b.len() && (b[i] == b'e' || b[i] == b'E') {
        let mut j = i + 1;
        if j < b.len() && (b[j] == b'+' || b[j] == b'-') {
            j += 1;
        }
        let es = j;
        while j < b.len() && b[j].is_ascii_digit() {
            j += 1;
        }
        if j > es {
            end = j;
        }
    }
    t[..end].parse::<f64>().unwrap_or(0.0)
}

/// printf("%f"): 6 decimals, "nan"/"-nan" and "inf"/"-inf" spelled as glibc.
pub fn c_f(x: f64) -> String {
    if x.is_nan() {
        if x.is_sign_negative() { "-nan".into() } else { "nan".into() }
    } else if x.is_infinite() {
        if x < 0.0 { "-inf".into() } else { "inf".into() }
    } else {
        format!("{x:.6}")
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn nums() {
        assert_eq!(strtol_i32("120.9"), 120);
        assert_eq!(strtol_i32("0x78"), 0);
        assert_eq!(strtol_i32("500000000000000000"), -743309312);
        assert_eq!(strtol_i32("-500000000000000000"), 743309312);
        assert_eq!(strtol_i32("4611686018427387904"), 0);
        assert_eq!(strtol_i32("9223372036854775807"), -1);
        assert_eq!(strtol_i32("99999999999999999999"), -1);
        assert_eq!(strtol_i32("-9223372036854775808"), 0);
        assert_eq!(strtol_i32("-99999999999999999999"), 0);
        assert_eq!(strtod("5e-1"), 0.5);
        assert_eq!(strtod("0x1p-1"), 0.5);
        assert_eq!(strtod("1.5abc"), 1.5);
        assert!(strtod("inf").is_infinite());
    }
}
