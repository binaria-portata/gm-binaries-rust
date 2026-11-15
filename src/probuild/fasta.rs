//! FASTA reader classes shared by the utility modes (re-notes/probuild/core.md
//! section 4), the report-section writer (section 5) and small sequence
//! helpers (line writer, run finder, definition-line word).
//!
//! probuild has two independent FASTA classes with different parsing rules;
//! which one a mode uses is fixed by main, and the differences are visible
//! in the output, so both are reproduced:
//!
//! - `SeqUtil` (object built by 0x405aa0; records at +0x48, letter
//!   translation table at +0x60, reader 0x423a8a).  Used by the modes in
//!   `seqmodes.rs`.  The handler chooses the table (default upper case,
//!   soft-masking, allow-X, details, protein) before reading.  A definition
//!   line ends at '\n' or '\r'; bytes that translate to 0 are dropped.
//! - `FastaX` (objects 0x405dce = FastaX(0) and 0x405e8a = FastaX(allow_x);
//!   table 0x44b2b8, reader 0x44b624).  Used by the modes in `statmodes.rs`.
//!   Every letter is kept (unknown letters become `N`), a definition line
//!   ends at '\n' only (a CR stays in it), and every chunk, even an empty
//!   one, becomes a record.
//!
//! Both objects carry a `Report` (map at +0x18), the original's generic
//! "named output sections" container used to write results.
//!
//! A record's sequence is stored after translation, so it only contains the
//! table's output letters.  Coordinates elsewhere are 0-based half-open
//! `[start, end)` pairs of `u64` unless stated otherwise.

use super::io::{self, err, load_file, OFile};
use std::collections::BTreeMap;

/// Byte string (C++ `std::string`).
pub type Bytes = Vec<u8>;

/// One FASTA record.
#[derive(Clone, Default, Debug)]
pub struct Record {
    /// Definition line including the leading '>' and without the line end
    /// (empty if the record has none, e.g. text before the first '>').
    pub def: Bytes,
    /// Sequence after translation through the reader's table.
    pub seq: Bytes,
}

/// C `isspace` in the C locale: blank, \t, \n, \v, \f, \r.
#[inline]
pub fn c_isspace(c: u8) -> bool {
    matches!(c, b' ' | b'\t' | b'\n' | 0x0b | 0x0c | b'\r')
}

// ---------------------------------------------------------------------------
// SeqUtil translation tables.  Each is a 256-entry byte map; 0 means "drop
// this byte".  Handlers build the default table and then patch it with the
// soft / allow-X functions, or replace it entirely.

/// 0x422ad4: default nucleotide table.  `acgt`/`ACGT` to upper case, `u`/`U`
/// to `T`, IUPAC ambiguity codes and `x`/`X` (both cases) to `N`; every other
/// byte (digits, '*', '-', ...) is dropped.
pub fn table_default() -> [u8; 256] {
    let mut t = [0u8; 256];
    for (a, b) in [(b'a', b'A'), (b'c', b'C'), (b'g', b'G'), (b't', b'T'), (b'u', b'T')] {
        t[a as usize] = b;
        t[a.to_ascii_uppercase() as usize] = b;
    }
    for c in b"nxrykmswbdhv" {
        t[*c as usize] = b'N';
        t[c.to_ascii_uppercase() as usize] = b'N';
    }
    t
}

/// 0x422e2c: "soft" patch, keeps lower case (soft-masked) letters lower
/// case: `acgt` stay, `u` -> `t`, `n` and the lower case ambiguity codes
/// (including `x`) -> `n`.  Upper case entries are not touched.
pub fn table_soft(t: &mut [u8; 256]) {
    for (a, b) in [(b'a', b'a'), (b'c', b'c'), (b'g', b'g'), (b't', b't'), (b'u', b't'), (b'n', b'n')] {
        t[a as usize] = b;
    }
    for c in b"xrykmswbdhv" {
        t[*c as usize] = b'n';
    }
}

/// 0x422fd2: allow the letter X (masked sequence) instead of mapping it to
/// N.  `X` -> `X`; `x` -> `x` only when `soft` (keep lower case), else `X`.
pub fn table_allow_x(t: &mut [u8; 256], soft: bool) {
    t[b'x' as usize] = if soft { b'x' } else { b'X' };
    t[b'X' as usize] = b'X';
}

/// 0x42304e: identity for all nucleotide and IUPAC letters in both cases
/// (`--stat --details`, which counts every letter separately).  This is a
/// fresh table: it replaces any soft / allow-X patches.
pub fn table_details() -> [u8; 256] {
    let mut t = [0u8; 256];
    for c in b"acgtunxrykmswbdhv" {
        t[*c as usize] = *c;
        t[c.to_ascii_uppercase() as usize] = c.to_ascii_uppercase();
    }
    t
}

/// 0x4233a6: protein alphabet (`--split_fasta`, `--select_fasta`): the 20
/// amino acids plus `u x b o z j`, both cases mapped to upper case.  As in the original:
/// the upper case `Z` entry is missing (the original sets `z` twice), so a
/// `Z` is dropped while a `z` becomes `Z`.
pub fn table_protein() -> [u8; 256] {
    let mut t = [0u8; 256];
    for c in b"arndceqghilkmfpswtyvuxbozj" {
        t[*c as usize] = c.to_ascii_uppercase();
        if *c != b'z' {
            t[c.to_ascii_uppercase() as usize] = c.to_ascii_uppercase();
        }
    }
    t
}

/// SeqUtil table lookup.  The original indexes with `(long)(signed char)c`,
/// so bytes >= 0x80 read the 128 bytes BEFORE the table (heap memory, not
/// deterministic).  The port treats them as unknown (0, dropped); such input
/// is outside the supported range (core.md section 10).
#[inline]
fn translate_byte(t: &[u8; 256], c: u8) -> u8 {
    if c >= 0x80 {
        0
    } else {
        t[c as usize]
    }
}

/// The SeqUtil object: translation table, records and report sections.
pub struct SeqUtil {
    /// Byte translation table used by `load` (0 = drop the byte).
    pub table: [u8; 256],
    /// Records read so far.  `load` appends; handlers that read several
    /// files clear it first.
    pub recs: Vec<Record>,
    /// Report sections written by `--stat` / `--gc`.
    pub report: Report,
}

impl SeqUtil {
    /// Fresh object with the default table.
    pub fn new() -> SeqUtil {
        SeqUtil { table: table_default(), recs: Vec::new(), report: Report::default() }
    }

    /// 0x423a8a: read FASTA file `name` and append its records to `recs`.
    /// False if the loader fails (message already printed).
    pub fn load(&mut self, name: &[u8]) -> bool {
        let buf = match load_file(name) {
            Some(b) => b,
            None => return false,
        };
        self.parse(&buf);
        true
    }

    /// Parse a whole file image.  Each round reads one record:
    ///
    /// 1. skip white space; if the next byte is '>', the definition runs up
    ///    to (not including) the next '\n' or '\r';
    /// 2. skip white space; the sequence runs up to the next '>' or the end;
    /// 3. if the definition or the sequence is non-empty, push the record
    ///    with the sequence translated through `table` (0 results dropped).
    ///
    /// Consequences: text before the first '>' forms a record with an empty
    /// definition; a '>' anywhere (even inside a sequence line) starts a new
    /// record; with CRLF line ends the CR ends the definition.
    pub fn parse(&mut self, b: &[u8]) {
        let end = b.len();
        // The original tests the character before the bound; past the end it
        // reads the string's NUL terminator.  `at` returns 0 there so that
        // the same operand order is safe.
        let at = |i: usize| if i < end { b[i] } else { 0 };
        let mut pos = 0;
        while pos < end {
            // 0x4238ae: definition line
            let mut cursor = pos;
            while c_isspace(at(cursor)) && cursor < end {
                cursor += 1;
            }
            let def_start = cursor;
            if at(cursor) == b'>' && cursor < end {
                while cursor < end && b[cursor] != b'\n' && b[cursor] != b'\r' {
                    cursor += 1;
                }
            }
            let def_end = cursor;
            // 0x423952: sequence
            while c_isspace(at(cursor)) && cursor < end {
                cursor += 1;
            }
            let seq_start = cursor;
            while at(cursor) != b'>' && cursor < end {
                cursor += 1;
            }
            let seq_end = cursor;
            pos = seq_end;
            if def_start != def_end || seq_start != seq_end {
                let mut seq = Vec::with_capacity(seq_end - seq_start);
                for &c in &b[seq_start..seq_end] {
                    let code = translate_byte(&self.table, c);
                    if code != 0 {
                        seq.push(code);
                    }
                }
                self.recs.push(Record { def: b[def_start..def_end].to_vec(), seq });
            }
        }
    }
}

/// 0x44b2b8: FastaX table.  Starts with `N` in all 256 entries, then
/// `acgtACGT` map to themselves, `u` -> `t`, `U` -> `T` and, with
/// `allow_x`, `x`/`X` -> `X`.  Only letters are looked up (the reader drops
/// non-letters first), so in effect every other letter becomes `N`.  Note
/// that lower case is preserved here, unlike the SeqUtil default table.
pub fn table_x(allow_x: bool) -> [u8; 256] {
    let mut t = [b'N'; 256];
    for c in b"acgtACGT" {
        t[*c as usize] = *c;
    }
    t[b'u' as usize] = b't';
    t[b'U' as usize] = b'T';
    if allow_x {
        t[b'x' as usize] = b'X';
        t[b'X' as usize] = b'X';
    }
    t
}

/// The FastaX object: translation table, records and report sections.
pub struct FastaX {
    /// Letter translation table (`table_x`).
    pub table: [u8; 256],
    /// Records read so far (`load` appends).
    pub recs: Vec<Record>,
    /// Report sections written by `--stat_fasta` / `--substring_n_distr`.
    pub report: Report,
}

impl FastaX {
    /// FastaX(0) (`allow_x == false`) or FastaX(1).
    pub fn new(allow_x: bool) -> FastaX {
        FastaX { table: table_x(allow_x), recs: Vec::new(), report: Report::default() }
    }

    /// 0x44b624: read FASTA file `name` and append its records.  Each round:
    ///
    /// 1. skip white space; if the next byte is '>', the definition runs up
    ///    to the next '\n' only (a '\r' stays in the definition);
    /// 2. the sequence runs from there to the next '>' or the end;
    /// 3. push the record (always, even if both parts are empty) with every
    ///    ASCII letter (`isalpha`) translated through `table`; all other
    ///    bytes, including white space, digits and bytes >= 0x80, are
    ///    dropped.
    ///
    /// A file of white space only therefore gives one empty record.
    pub fn load(&mut self, name: &[u8]) -> bool {
        let buf = match load_file(name) {
            Some(b) => b,
            None => return false,
        };
        let end = buf.len();
        let mut pos = 0;
        while pos < end {
            // 0x44b40c: definition line
            let mut cursor = pos;
            while cursor < end && c_isspace(buf[cursor]) {
                cursor += 1;
            }
            let def_start = cursor;
            if cursor < end && buf[cursor] == b'>' {
                while cursor < end && buf[cursor] != b'\n' {
                    cursor += 1;
                }
            }
            let def_end = cursor;
            // 0x44b4ac: sequence up to the next '>' (no white space skip:
            // white space is dropped below anyway)
            let seq_start = def_end;
            let mut seq_end = seq_start;
            while seq_end < end && buf[seq_end] != b'>' {
                seq_end += 1;
            }
            pos = seq_end;
            // 0x44b516: letters only (isalpha), through the table
            let mut seq = Vec::with_capacity(seq_end - seq_start);
            for &c in &buf[seq_start..seq_end] {
                if c.is_ascii_alphabetic() {
                    seq.push(self.table[c as usize]);
                }
            }
            self.recs.push(Record { def: buf[def_start..def_end].to_vec(), seq });
        }
        true
    }
}

/// Report sections: the original's `map<string, pair<string,string>>` at
/// object +0x18 (key -> (head, body)) plus a request list at +0 (0x44822c
/// appends), written by 0x421806.  Modes put their result text into a
/// section and then write the report to a file or stdout.
#[derive(Default)]
pub struct Report {
    /// Section key -> (head, body).  Iterated in key order when writing all.
    pub map: BTreeMap<Bytes, (Bytes, Bytes)>,
    /// Keys to write, in this order; empty = write every section.
    pub requests: Vec<Bytes>,
}

impl Report {
    /// 0x420f68: set (or replace) section `key`.
    pub fn set(&mut self, key: &str, head: &[u8], body: &[u8]) {
        self.map.insert(key.as_bytes().to_vec(), (head.to_vec(), body.to_vec()));
    }

    /// 0x44822c: add `key` to the request list.
    pub fn request(&mut self, key: &str) {
        self.requests.push(key.as_bytes().to_vec());
    }

    /// 0x421806: write the report to file `name` (truncated; an open failure
    /// is silently ignored) or to stdout if `name` is empty.  Each section is
    /// written as `head + body + "\n"`.
    ///
    /// - Empty request list (0x421290): all sections in key order.
    /// - Otherwise every request is validated first (0x4210ee) and, if one
    ///   fails, an error goes to stderr, nothing is written and the result
    ///   is false: empty key
    ///   `error: no value in request`, unknown key
    ///   `error: no entry to report, KEY`, empty body
    ///   `error: no value to report, KEY`.  Then the requested sections are
    ///   written in request order (0x4214e0).
    pub fn write(&mut self, name: &[u8]) -> bool {
        let mut text = Vec::new();
        if self.requests.is_empty() {
            // 0x421290: all entries in key order
            for (head, body) in self.map.values() {
                text.extend_from_slice(head);
                text.extend_from_slice(body);
                text.push(b'\n');
            }
        } else {
            // 0x4210ee: validation
            for r in &self.requests {
                if r.is_empty() {
                    err(b"error: no value in request\n");
                    return false;
                }
                match self.map.get(r) {
                    None => {
                        let mut m = b"error: no entry to report, ".to_vec();
                        m.extend_from_slice(r);
                        m.push(b'\n');
                        err(&m);
                        return false;
                    }
                    Some((_, body)) if body.is_empty() => {
                        let mut m = b"error: no value to report, ".to_vec();
                        m.extend_from_slice(r);
                        m.push(b'\n');
                        err(&m);
                        return false;
                    }
                    _ => {}
                }
            }
            // 0x4214e0
            for r in &self.requests {
                let (head, body) = &self.map[r];
                text.extend_from_slice(head);
                text.extend_from_slice(body);
                text.push(b'\n');
            }
        }
        if name.is_empty() {
            io::out(&text);
        } else {
            let mut f = OFile::open(name, false);
            f.w(&text);
            f.close();
        }
        true
    }
}

/// Line writer 0x42c30a: append `s` to `o` in lines of `lpl` letters, each
/// followed by '\n'.  `lpl == 0` means one line of the whole sequence.  The
/// arithmetic is unsigned 64-bit as in the original: a negative `lpl`
/// becomes a huge width, so it also gives a single line.  An empty sequence
/// writes nothing (not even a newline).
pub fn write_lines(s: &[u8], lpl: i32, o: &mut Vec<u8>) {
    let len = s.len() as u64;
    // Line width.  For lpl == 0 the original uses the length converted to
    // int (only differs for sequences of 2 GiB or more).
    let width: i64 = if lpl == 0 { s.len() as i32 as i64 } else { lpl as i64 };
    let mut line_start: u64 = 0;
    while line_start < len {
        // `width + line_start < len` with wrapping u64 arithmetic, as in the
        // original; a wrapped (negative) width ends the line at `len`.
        let line_end =
            if (width as u64).wrapping_add(line_start) < len { line_start.wrapping_add(width as u64) } else { len };
        let (a, b) = (line_start.min(len) as usize, line_end.min(len) as usize);
        if a <= b {
            o.extend_from_slice(&s[a..b]);
        }
        o.push(b'\n');
        line_start = line_end;
    }
}

/// Maximal runs of byte `c` in `s` as `[start, end)` (the run loops of
/// 0x42ca2c / 0x44dc16).
pub fn runs_of(s: &[u8], c: u8) -> Vec<(u64, u64)> {
    let mut v = Vec::new();
    let mut i = 0;
    while i < s.len() {
        if s[i] == c {
            let b = i;
            while i < s.len() && s[i] == c {
                i += 1;
            }
            v.push((b as u64, i as u64));
        } else {
            i += 1;
        }
    }
    v
}

/// Maximal runs of bytes from `set` as `[start, end)` (0x42cb2e:
/// `find_first_of` / `find_first_not_of`).  Bytes of the set are
/// interchangeable within a run (e.g. `acgtnx` for soft-masked regions).
pub fn runs_of_set(s: &[u8], set: &[u8]) -> Vec<(u64, u64)> {
    let mut v = Vec::new();
    let mut i = 0;
    while i < s.len() {
        if set.contains(&s[i]) {
            let b = i;
            while i < s.len() && set.contains(&s[i]) {
                i += 1;
            }
            v.push((b as u64, i as u64));
        } else {
            i += 1;
        }
    }
    v
}

/// Sequence id taken from a definition line by 0x430e8a / 0x42d304: skip the
/// leading '>' characters, then blanks (space, tab), and return the word up
/// to the next blank.  `None` means there is no such word (e.g. the line is
/// `>` only, or `>` followed by blanks): there the original's `substr`
/// throws `std::out_of_range`, and callers reproduce that.
pub fn first_word_after_gt(d: &[u8]) -> Option<Bytes> {
    let after_gt = d.iter().position(|&c| c != b'>');
    let word_start = match after_gt {
        Some(p) => d[p..].iter().position(|&c| c != b' ' && c != b'\t').map(|x| x + p),
        None => None,
    };
    let word_start = word_start?;
    let word_end =
        d[word_start..].iter().position(|&c| c == b' ' || c == b'\t').map(|x| x + word_start).unwrap_or(d.len());
    Some(d[word_start..word_end].to_vec())
}
