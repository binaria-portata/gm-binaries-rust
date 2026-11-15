//! Sequence input: FASTA reader (0x4095f0 with 0x4087a0, 0x4088c0, 0x408b20,
//! 0x408e20) and the k-mer index used by all Markov chains (0x4099a0 ->
//! 0x409830). See re-notes/input.md sections 2-3.
//!
//! Second stage of the pipeline. The sequence is stored as one byte per base
//! ("code"): A=0, C=1, G=2, T/U=3, anything else (N, IUPAC, other letters)
//! = 4. Case is discarded, so soft-masking is not visible here (it only
//! reaches the scores through `-d` `mask` lines and `-k`).
//!
//! The reader is deliberately primitive, like the original: only the first
//! defline is recognised, and every alphabetic byte after it becomes
//! sequence, including the letters of any later `>name` line. Coordinates
//! everywhere in the program are 0-based positions into [`Sequence::codes`].

use std::path::Path;

/// Base code of A (also the digit value in k-mer indices).
pub const A: u8 = 0;
/// Base code of C.
pub const C: u8 = 1;
/// Base code of G.
pub const G: u8 = 2;
/// Base code of T and U.
pub const T: u8 = 3;
/// Base code of every other letter (N, IUPAC codes, X, ...).
pub const N: u8 = 4;

/// Sequence object (in main, 0xb8 bytes in the original).
#[allow(dead_code)]
#[derive(Debug, Clone)]
pub struct Sequence {
    /// file name as given on the command line
    pub file: String,
    /// per-base codes 0..4 (seq+0x10); its length is the sequence length L
    pub codes: Vec<u8>,
    /// k-mer index per position for the model order (seq+0x18), filled by
    /// [`Sequence::build_kmer`]; see [`kmer_index`]
    pub kmer: Vec<i32>,
    /// evidence marks per position (seq+0x28), only with -d; bit layout in
    /// the `evidence` module
    pub marks: Option<Vec<u32>>,
    /// FASTA definition line including '>', excluding the line end
    /// (seq+0x30); None without one (printed as "(null)")
    pub defline: Option<String>,
    /// sequence id (seq+0x38): the `-s` tag, else defline token 1, else "seq"
    pub seqid: String,
    /// defline token 2 (seq+0x58, used by -z as the sequence id)
    pub token2: String,
    /// defline token 3 as int (seq+0x78, used by -z as coordinate offset + 1)
    pub token3: i32,
    /// number of bases per code 0..4
    pub counts: [i64; 5],
    /// G+C percent over A/C/G/T only (seq+0x98); code-4 bases are excluded
    /// from numerator and denominator
    pub gc: f64,
}

/// C `isspace` in the "C" locale.
#[inline]
fn is_c_space(b: u8) -> bool {
    matches!(b, b' ' | b'\t' | b'\n' | 0x0b | 0x0c | b'\r')
}

/// Letter to base code (table built by 0x4088c0).
fn code_of(b: u8) -> u8 {
    match b {
        b'A' | b'a' => A,
        b'C' | b'c' => C,
        b'G' | b'g' => G,
        b'T' | b't' | b'U' | b'u' => T,
        _ => N,
    }
}

/// Error message for stdout; every input error is fatal (exit 1).
pub type InputError = String;

/// First char of the escape range for bytes that are not valid UTF-8
/// (plane 16 private use area; byte b is U+10FF00 + b).
const ESC_BASE: u32 = 0x10FF00;

/// Bytes to `String`, keeping every byte: valid UTF-8 stays as is, each byte
/// of an invalid sequence becomes the char U+10FF00 + byte. The original
/// works on raw bytes (file names from argv, FASTA definition lines, -d
/// lines), and [`encode`] turns such strings back into the same bytes when
/// they are printed or used as file names. (A literal U+10FF80..U+10FFFF in
/// the input would be written back as a single byte; such chars do not
/// occur in GeneMark data.)
pub fn decode(b: &[u8]) -> String {
    let mut s = String::with_capacity(b.len());
    for chunk in b.utf8_chunks() {
        s.push_str(chunk.valid());
        for &x in chunk.invalid() {
            s.push(char::from_u32(ESC_BASE + x as u32).unwrap_or('?'));
        }
    }
    s
}

/// Inverse of [`decode`].
pub fn encode(s: &str) -> Vec<u8> {
    if s.is_ascii() {
        return s.as_bytes().to_vec();
    }
    let mut v = Vec::with_capacity(s.len());
    for c in s.chars() {
        let u = c as u32;
        if (ESC_BASE + 0x80..=ESC_BASE + 0xff).contains(&u) {
            v.push((u - ESC_BASE) as u8);
        } else {
            let mut buf = [0u8; 4];
            v.extend_from_slice(c.encode_utf8(&mut buf).as_bytes());
        }
    }
    v
}

/// File name (possibly with escaped bytes, see [`decode`]) as a path.
pub fn os_path(s: &str) -> std::path::PathBuf {
    use std::os::unix::ffi::OsStringExt;
    std::path::PathBuf::from(std::ffi::OsString::from_vec(encode(s)))
}

impl Sequence {
    /// 0x4095f0(seq, file, tag): read and parse the FASTA file.
    pub fn read(path: &str, tag: &str) -> Result<Sequence, InputError> {
        let pb = os_path(path);
        let p: &Path = &pb;
        // the original calls stat() first: a missing file gives "not found",
        // a file that exists but cannot be read gives "reading file"
        let meta = std::fs::metadata(p).map_err(|_| format!("gmhmme3: error, file <{path}> is not found\n"))?;
        let _ = meta;
        let data = std::fs::read(p).map_err(|_| format!("gmhmme3: error, reading file <{path}>\n"))?;
        // the buffer is a C string: stop at an embedded NUL
        let data = match data.iter().position(|&b| b == 0) {
            Some(z) => &data[..z],
            None => &data[..],
        };
        Ok(Self::parse(path, data, tag)?)
    }

    /// 0x408b20 + 0x408e20: parse the file contents. `tag` is the `-s`
    /// value (empty if not given).
    pub fn parse(file: &str, data: &[u8], tag: &str) -> Result<Sequence, InputError> {
        // skip leading whitespace; a defline is recognised only if the first
        // non-space byte is '>'
        let mut i = 0;
        while i < data.len() && is_c_space(data[i]) {
            i += 1;
        }
        let mut defline = None;
        if i < data.len() && data[i] == b'>' {
            let mut j = i;
            while j < data.len() && data[j] != b'\n' && data[j] != b'\r' {
                j += 1;
            }
            defline = Some(decode(&data[i..j]));
            i = j;
        }
        // Every alphabetic byte from here to EOF is a base; everything else
        // (digits, '-', '*', whitespace, '>', bytes >= 0x80) is dropped.
        // A second record's header letters therefore become sequence.
        let mut codes = Vec::with_capacity(data.len() - i);
        let mut counts = [0i64; 5];
        for &b in &data[i..] {
            // isalpha in the C locale, called with a signed char: bytes
            // >= 0x80 are negative and never alphabetic
            if b.is_ascii_alphabetic() {
                let c = code_of(b);
                counts[c as usize] += 1;
                codes.push(c);
            }
        }
        if codes.is_empty() {
            return Err("gmhmme3 : error, sequence has size 0\n".into());
        }
        let total = counts[0] + counts[1] + counts[2] + counts[3];
        let gc = if total != 0 { (counts[1] + counts[2]) as f64 * 100.0 / total as f64 } else { 0.0 };

        // 0x408e20: sequence id. With -s the defline is not tokenised at all
        // (token2 empty, token3 0).
        let mut seqid = String::new();
        let mut token2 = String::new();
        let mut token3 = 0;
        if !tag.is_empty() {
            seqid = tag.to_string();
        } else if let Some(d) = &defline {
            let toks = defline_tokens(d);
            if let Some(t1) = toks.first() {
                seqid = t1.clone();
            }
            if let Some(t2) = toks.get(1) {
                token2 = t2.clone();
            }
            if let Some(t3) = toks.get(2) {
                token3 = crate::model::strtol(t3) as i32;
                // cleared only when a third token exists and parses to 0
                if token3 == 0 {
                    token2.clear();
                }
            }
        }
        if seqid.is_empty() {
            seqid = "seq".into();
        }
        Ok(Sequence {
            file: file.to_string(),
            codes,
            kmer: Vec::new(),
            marks: None,
            defline,
            seqid,
            token2,
            token3,
            counts,
            gc,
        })
    }

    /// Sequence length L in bases.
    pub fn len(&self) -> usize {
        self.codes.len()
    }

    /// 0x4099a0 + 0x409830: check the length against the Markov order and
    /// build [`Sequence::kmer`].
    pub fn build_kmer(&mut self, order: i32) -> Result<(), InputError> {
        let l = self.codes.len() as i64;
        if l <= order as i64 {
            return Err(format!("gmhmme3 : error, size of sequence {l}, order {order}\n"));
        }
        // 4^(order+1) must fit an int (index of the invalid / N slot)
        if 2 * order + 2 >= 32 {
            return Err(format!("gmhmme3 : error, can't convert address for order {order} to 'int'\n"));
        }
        self.kmer = kmer_index(&self.codes, order);
        Ok(())
    }
}

/// 0x409830(seq, order, invalid = 4^(order+1), tail = 0): per-position
/// index of the (order+1)-mer starting at that position.
///
/// `index[i]` is the base-4 number of `codes[i..=i+order]` with the first
/// base as the most significant digit. A window containing an N gets
/// `4^(order+1)`, which is the "N slot" at the end of every Markov table.
/// The last `order` positions have no complete window: they get the N slot
/// if any base from there to the end is N, else 0 (index of A...A; as in
/// the original, reproduced because such entries can be read by the score
/// tracks near the sequence end).
pub fn kmer_index(codes: &[u8], order: i32) -> Vec<i32> {
    let len = codes.len();
    let invalid = 4i32.pow((order + 1) as u32);
    let ord = order as usize;
    let mut index = vec![0i32; len];
    for i in 0..len {
        if i + ord >= len {
            index[i] = if codes[i..].contains(&N) { invalid } else { 0 };
        } else {
            let window = &codes[i..=i + ord];
            if window.contains(&N) {
                index[i] = invalid;
            } else {
                let mut v = 0i32;
                for &x in window {
                    v = (v << 2) | x as i32;
                }
                index[i] = v;
            }
        }
    }
    index
}

/// Tokenisation of 0x408e20: up to three tokens separated by runs of ' ' or
/// '\t', starting at position 1 (skipping '>', so `>  X1` gives "X1").
fn defline_tokens(d: &str) -> Vec<String> {
    let b = d.as_bytes();
    let mut out = Vec::new();
    let mut pos = 1;
    while out.len() < 3 {
        while pos < b.len() && (b[pos] == b' ' || b[pos] == b'\t') {
            pos += 1;
        }
        if pos >= b.len() {
            break;
        }
        let s = pos;
        while pos < b.len() && b[pos] != b' ' && b[pos] != b'\t' {
            pos += 1;
        }
        // (split at ASCII bytes only, so escaped chars stay whole)
        out.push(d[s..pos].to_string());
    }
    out
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn parse_basic() {
        let s = Sequence::parse("x", b"  >id two 3\r\nACgtN\nU-*1\n>s2\n", "").unwrap();
        assert_eq!(s.defline.as_deref(), Some(">id two 3"));
        assert_eq!(s.seqid, "id");
        assert_eq!(s.token2, "two");
        assert_eq!(s.token3, 3);
        // second header letters become sequence: s2 -> s (N)
        assert_eq!(s.codes, vec![0, 1, 2, 3, 4, 3, 4]);
    }

    #[test]
    fn kmer_tail() {
        let c = [0u8, 1, 2, 3, 4, 0];
        let k = kmer_index(&c, 1);
        assert_eq!(k, vec![1, 6, 11, 16, 16, 0]);
    }
}
