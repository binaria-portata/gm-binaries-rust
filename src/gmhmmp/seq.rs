//! gmhmmp FASTA reading and per-sequence treatment (input.md 2-3):
//! 0x408270 (records), 0x407ad0 (encode, GC), 0x407650 (id), 0x407420 (gap
//! split), 0x407260 (R/Y stops), 0x40c200 (codon image).
//!
//! Base codes: A=0, C=1, G=2, T/U=3, any other letter 4 (N, IUPAC, X;
//! case ignored). Three views of a record's codes are used downstream:
//!
//! * `Record::orig`: before gap splitting, with 5 / 6 written at R / Y
//!   positions resolved to a stop (used for N ends and the -A/-D output);
//! * the working `Record::codes` after gap splitting, copied by main
//!   before R/Y resolution (Markov image and RBS scores see R/Y as 4);
//! * `Record::codes` after R/Y resolution, the input of the codon image.

/// One FASTA record as returned by a call of 0x408270.
#[derive(Debug, Clone)]
pub struct Record {
    /// defline including '>' (">empty-fasta-def-line" for text before the
    /// first '>')
    pub defline: String,
    /// sequence id (0x407650), the key for `-e` evidence
    pub id: String,
    /// working codes 0..4 (gap split applied unless -n)
    pub codes: Vec<u8>,
    /// copy before gap split; 5/6 at resolved R/Y (the -A/-D source).
    /// Filled by main (empty after `encode`).
    pub orig: Vec<u8>,
    /// 0-based letter positions of R/r and Y/y (encoded as 4)
    pub r_pos: Vec<usize>,
    pub y_pos: Vec<usize>,
    /// GC percent of the A/C/G/T letters, rounded half up (0 if none);
    /// selects the meta model section, may be overridden by `[gc=]` (-c)
    pub gc: i32,
}

/// Letter code table 0x4077f0.
#[inline]
fn code_of(c: u8) -> u8 {
    match c {
        b'a' | b'A' => 0,
        b'c' | b'C' => 1,
        b'g' | b'G' => 2,
        b't' | b'T' | b'u' | b'U' => 3,
        _ => 4,
    }
}

/// Line-based record reader with the binary's getline semantics: a line
/// is a header only if it starts with '>' in column 0; all other lines are
/// appended without separator (so `;comment` lines count as sequence).
/// Text before the first header forms a record of its own.
pub struct Reader {
    lines: Vec<Vec<u8>>,
    /// index of the next unread line
    next: usize,
    /// the header line that ended the previous record, plus "\n"
    pending: Vec<u8>,
}

impl Reader {
    pub fn new(data: Vec<u8>) -> Reader {
        let mut lines: Vec<Vec<u8>> = data.split(|&b| b == b'\n').map(|l| l.to_vec()).collect();
        if data.ends_with(b"\n") {
            lines.pop();
        }
        if data.is_empty() {
            lines.clear();
        }
        Reader { lines, next: 0, pending: Vec::new() }
    }

    /// 0x408270: next raw record text, `">defline\n"` followed by all its
    /// sequence lines concatenated. Empty when nothing is left; a record
    /// may also have 0 letters.
    pub fn next_raw(&mut self) -> Vec<u8> {
        let mut raw = std::mem::take(&mut self.pending);
        while self.next < self.lines.len() {
            let ln = &self.lines[self.next];
            self.next += 1;
            if ln.first() == Some(&b'>') {
                let mut p = ln.clone();
                p.push(b'\n');
                self.pending = p;
                break;
            }
            raw.extend_from_slice(ln);
        }
        raw
    }

    /// A header is waiting: more records follow (main's eof test).
    pub fn has_pending(&self) -> bool {
        !self.pending.is_empty()
    }
}

/// 0x407650: id = from the first non-blank after '>' up to the first
/// SPACE (a TAB does not end it). Edge cases as in input.md 3.4 (`>`
/// alone gives `>`, an all-blank defline gives "" or "\t").
fn seqid(s: &[u8]) -> Vec<u8> {
    let n = s.len();
    let start = if n >= 2 && (s[1] == b' ' || s[1] == b'\t') {
        let mut i = 1;
        loop {
            i += 1;
            if i >= n {
                break n - 1;
            }
            if s[i] != b' ' && s[i] != b'\t' {
                break i;
            }
        }
    } else if n >= 2 {
        1
    } else {
        n.saturating_sub(1)
    };
    if n == 0 {
        return Vec::new();
    }
    if s[start] == b' ' {
        return Vec::new();
    }
    let mut j = start;
    loop {
        j += 1;
        if j >= n {
            return s[start..n].to_vec();
        }
        if s[j] == b' ' {
            return s[start..j].to_vec();
        }
    }
}

/// 0x407ad0: encode a raw record (None if the raw text is empty).
///
/// The defline ends at the first '\r' or '\n'; letters are collected from
/// there on, so text after a '\r' inside the defline line counts as
/// sequence. Only `isalpha` bytes (C locale) are kept. GC:
/// `(int)((C+G) * 100.0 / (A+C+G+T) + 0.5)`.
pub fn encode(raw: &[u8]) -> Option<Record> {
    if raw.is_empty() {
        return None;
    }
    let (defl, p): (Vec<u8>, usize) = if raw[0] == b'>' {
        let mut j = 1;
        while j < raw.len() && raw[j] != b'\r' && raw[j] != b'\n' {
            j += 1;
        }
        (raw[..j].to_vec(), j)
    } else {
        (b">empty-fasta-def-line".to_vec(), 0)
    };
    let mut codes = Vec::with_capacity(raw.len());
    let mut r_pos = Vec::new();
    let mut y_pos = Vec::new();
    let mut cnt = [0i64; 5];
    for &c in &raw[p..] {
        if c.is_ascii_alphabetic() {
            if c == b'r' || c == b'R' {
                r_pos.push(codes.len());
            }
            if c == b'y' || c == b'Y' {
                y_pos.push(codes.len());
            }
            let v = code_of(c);
            cnt[v as usize] += 1;
            codes.push(v);
        }
    }
    let tot = cnt[0] + cnt[1] + cnt[2] + cnt[3];
    let gc = if tot > 0 { ((cnt[1] + cnt[2]) as f64 * 100.0 / tot as f64 + 0.5) as i32 } else { 0 };
    let id = seqid(&defl);
    Some(Record {
        defline: String::from_utf8_lossy(&defl).into_owned(),
        id: String::from_utf8_lossy(&id).into_owned(),
        orig: Vec::new(),
        codes,
        r_pos,
        y_pos,
        gc,
    })
}

/// First 24 bases of a long N run: CATC x3 then GTAA x3. Gives a reverse
/// start (CAT) and a forward stop (TAA) in every frame, so genes left of
/// the gap end there.
const GAP_HEAD: [u8; 24] = [1, 0, 3, 1, 1, 0, 3, 1, 1, 0, 3, 1, 2, 3, 0, 0, 2, 3, 0, 0, 2, 3, 0, 0];
/// Last 24 bases of a long N run: TTAC x3 then GATG x3. Gives a reverse
/// stop (TTA) and a forward start (ATG) in every frame, so genes right of
/// the gap can begin there.
const GAP_TAIL: [u8; 24] = [3, 3, 0, 1, 3, 3, 0, 1, 3, 3, 0, 1, 2, 0, 3, 2, 2, 0, 3, 2, 2, 0, 3, 2];

/// 0x407420 (default, off with -n): every maximal run of code 4 of length
/// >= 49 gets the stop/start patterns at both ends, the middle stays 4, so
/// incomplete genes end and start at the gap. The settings value 48 and
/// the usage text ("longer then 20") do not describe the real threshold.
/// After a rewritten run the scan resumes at j + 1 (j = first base after
/// the run, which is not 4 anyway).
pub fn split_gaps(seq: &mut [u8]) {
    let l = seq.len();
    let mut i = 0;
    while i < l {
        if seq[i] != 4 {
            i += 1;
            continue;
        }
        let mut j = i + 1;
        while j < l && seq[j] == 4 {
            j += 1;
        }
        if j - i >= 49 {
            seq[i..i + 24].copy_from_slice(&GAP_HEAD);
            seq[j - 24..j].copy_from_slice(&GAP_TAIL);
            i = j + 1;
        } else {
            i = j;
        }
    }
}

/// 0x407260: resolve R/Y that make a certain stop codon. Every rule needs
/// TAA > 0; R in TAR (needs TAG > 0) or TRA (needs TGA > 0) becomes A,
/// Y in YTA (reverse TAR) or TYA (reverse TRA) becomes T. `orig` gets 5
/// (R) / 6 (Y) at resolved positions. R rules run first; neighbours are
/// read from the working sequence, so earlier rewrites are seen.
/// Arguments: the model's TAA_, TAG_, TGA_ values (only `> 0` matters).
pub fn resolve_ry(rec: &mut Record, taa: f64, tag: f64, tga: f64) {
    let l = rec.codes.len();
    let w = &mut rec.codes;
    for &p in &rec.r_pos {
        if taa > 0.0 {
            if tag > 0.0 && p > 1 && w[p - 2] == 3 && w[p - 1] == 0 {
                w[p] = 0;
                rec.orig[p] = 5;
            }
            if tga > 0.0 && p > 0 && p + 1 < l && w[p - 1] == 3 && w[p + 1] == 0 {
                w[p] = 0;
                rec.orig[p] = 5;
            }
        }
    }
    for &p in &rec.y_pos {
        if taa > 0.0 {
            if tag > 0.0 && p + 2 < l && w[p + 1] == 3 && w[p + 2] == 0 {
                w[p] = 3;
                rec.orig[p] = 6;
            }
            if tga > 0.0 && p > 0 && p + 1 < l && w[p - 1] == 3 && w[p + 1] == 0 {
                w[p] = 3;
                rec.orig[p] = 6;
            }
        }
    }
}

/// 0x40c200: codon image. Entry i = flags of the codon ending at i
/// (`table[16*a + 4*b + c]`, `table[64]` = 0 if any base is 4), with the
/// frame tag `((i % 3) + 1) << 6` ORed in when nonzero. Entries 0 and 1 are
/// `table[64]`. Returns the image and the count of nonzero entries (the
/// binary overwrites the encoded sequence in place).
pub fn codon_image(seq: &[u8], table: &[u8; 65]) -> (Vec<u8>, usize) {
    let l = seq.len();
    let mut out = vec![0u8; l];
    let mut cnt = 0;
    if l > 0 {
        out[0] = table[64];
    }
    if l > 1 {
        out[1] = table[64];
    }
    for i in 2..l {
        let (a, b, c) = (seq[i - 2], seq[i - 1], seq[i]);
        let mut v = if a == 4 || b == 4 || c == 4 { table[64] } else { table[16 * a as usize + 4 * b as usize + c as usize] };
        if v != 0 {
            cnt += 1;
            v |= match i % 3 {
                2 => 0xc0,
                0 => 0x40,
                _ => 0x80,
            };
        }
        out[i] = v;
    }
    (out, cnt)
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn id_rules() {
        assert_eq!(seqid(b">  rec1\tsecond field"), b"rec1\tsecond".to_vec());
        assert_eq!(seqid(b">a b"), b"a".to_vec());
    }
    #[test]
    fn gaps() {
        let mut s = vec![4u8; 60];
        s.insert(0, 0);
        split_gaps(&mut s);
        assert_eq!(&s[1..25], &GAP_HEAD);
        assert_eq!(&s[37..61], &GAP_TAIL);
    }
}
