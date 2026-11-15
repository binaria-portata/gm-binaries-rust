//! Gene list reader of probuild 2.52 (`--geneset`; 0x409afe -> 0x409a66 ->
//! 0x4098dc, line parsers 0x4094a2 / 0x4092ea).  models.md §0.3.
//!
//! The gene lists are GeneMarkS-T / gmhmmp `.lst` predictions run without
//! RBS output, one gene per line with exactly six tokens:
//!
//! ```text
//! <gene number>  <strand +|->  <left>  <right>  <length>  <class>
//! ```
//!
//! with 1-based inclusive coordinates; a leading `<` on left / `>` on right
//! marks a partial gene.  Used by `--mkmod`, `--parse`, `--mkchain`,
//! `--pre_start`, `--STAT_NON_VS_ALL`, `--GENE_SET_STAT`, `--ORF_PER_SEQ`,
//! `--compare`; `--intergenic` uses the two-token "pair" format.
//!
//! IMPORTANT: every line of the file becomes one [`Gene`] entry, valid or
//! not.  Header lines, blank lines and column titles are entries with
//! strand 0 and left = right = -1.  Several modes depend on this (the
//! neighbour test of the pre-start sequences, the entry counts of
//! `--compare`, the pairs of `--GENE_SET_STAT`), so invalid lines must not
//! be filtered out.

use super::io;
use super::pm::P;

/// One gene list entry (16 bytes in the original; offsets in the field
/// docs).  Coordinates are 0-based and inclusive after parsing.
#[derive(Clone, Debug)]
pub struct Gene {
    /// +0: 0 = none (not a gene line), 1 = '+', 2 = '-'.
    pub strand: i32,
    /// +4: 0-based first base, -1 if the line is not a usable gene.
    pub left: i32,
    /// +8: 0-based last base (inclusive), -1 if not usable.
    pub right: i32,
    /// +0xc: left end complete (false if the coordinate had a '<').
    pub left_complete: bool,
    /// +0xd: right end complete (false if the coordinate had a '>').
    pub right_complete: bool,
    /// +0xe: a complete, well-formed gene line of class 1 or 2.
    pub valid: bool,
    /// +0xf: gene class from the 6th column (1 or 2 for valid genes; used
    /// by `--SEPARATE_COD` to split two coding models).  Stored as a signed
    /// char, so the parsed int is truncated to 8 bits.
    pub class: i8,
}

impl Default for Gene {
    /// The entry of a line that is not a gene: no strand, no positions,
    /// both ends "complete", not valid, class 0.
    fn default() -> Gene {
        Gene { strand: 0, left: -1, right: -1, left_complete: true, right_complete: true, valid: false, class: 0 }
    }
}

/// Splits a line into maximal runs of chars not in " \t\n\r" (0x40912c).
fn tokens(s: &[u8]) -> Vec<&[u8]> {
    s.split(|c| b" \t\n\r".contains(c)).filter(|t| !t.is_empty()).collect()
}

/// Token to int (0x409222): -1 if the token has a non-digit, otherwise the
/// `istringstream >> int` value (empty -> 0, overflow saturates at
/// INT_MAX).
fn parse_int(s: &[u8]) -> i32 {
    if !s.iter().all(|c| c.is_ascii_digit()) {
        return -1;
    }
    let mut v: i64 = 0;
    for &c in s {
        v = v * 10 + (c - b'0') as i64;
        if v > i32::MAX as i64 {
            return i32::MAX;
        }
    }
    v as i32
}

/// Parses one `.lst` line (pair_format 0, 0x4094a2).
///
/// The checks run in this order and each early return keeps whatever was
/// set so far (that is why e.g. a gene with a length not divisible by 3
/// keeps its strand and positions but is not valid):
///
/// 1. not exactly 6 tokens -> default entry;
/// 2. strand token not "+"/"-" -> default entry;
/// 3. left - 1 < 0 (also for a non-numeric token) -> left = right = -1;
/// 4. right - 1 < 0 or right < left -> left = right = -1;
/// 5. length not a multiple of 3 -> positions kept, not valid;
/// 6. `maxpos >= 1 && right > maxpos` -> positions kept, not valid.
///    `right` is 0-based, so right == maxpos (one past the sequence end)
///    still passes (inherited off-by-one);
/// 7. class = int(token 6) as signed char; valid iff class is 1 or 2.
///
/// The gene number (token 1) and length (token 5) are ignored.
fn parse_lst_line(line: &[u8], maxpos: i32) -> Gene {
    let mut g = Gene::default();
    let t = tokens(line);
    if t.len() != 6 {
        return g;
    }
    g.strand = match t[1] {
        b"+" => 1,
        b"-" => 2,
        _ => return g,
    };
    let mut left_tok = t[2];
    if left_tok[0] == b'<' {
        g.left_complete = false;
        left_tok = &left_tok[1..];
    }
    g.left = parse_int(left_tok).wrapping_sub(1);
    if g.left < 0 {
        g.left = -1;
        g.right = -1;
        return g;
    }
    let mut right_tok = t[3];
    if right_tok[0] == b'>' {
        g.right_complete = false;
        right_tok = &right_tok[1..];
    }
    g.right = parse_int(right_tok).wrapping_sub(1);
    if g.right < 0 || g.right < g.left {
        g.left = -1;
        g.right = -1;
        return g;
    }
    if (g.right - g.left + 1) % 3 != 0 {
        return g;
    }
    if maxpos >= 1 && g.right > maxpos {
        return g;
    }
    g.class = parse_int(t[5]) as i8;
    g.valid = (1..=2).contains(&g.class);
    g
}

/// Parses one line of the pair format (0x4092ea, `--intergenic` only):
/// exactly two tokens `a b` (1-based inclusive region); left =
/// offset + a - 1, right = offset + b - 1 in wrapping int arithmetic;
/// valid iff right >= left; strand stays 0.  A non-numeric token counts as
/// -1 (see [`parse_int`]).
fn parse_pair_line(line: &[u8], offset: i32) -> Gene {
    let mut g = Gene::default();
    let t = tokens(line);
    if t.len() != 2 {
        return g;
    }
    g.left = offset.wrapping_add(parse_int(t[0])).wrapping_sub(1);
    g.right = offset.wrapping_add(parse_int(t[1])).wrapping_sub(1);
    g.valid = g.right >= g.left;
    g
}

/// Splits a loaded gene list into lines and parses every line into one
/// entry (0x4098dc).
///
/// Line 1 is `[0, first '\n')`; every later line STARTS with its
/// preceding '\n' (harmless, '\n' is a token delimiter).  A final '\n'
/// does not start another line, so "" gives 0 entries and "\n" gives 1.
/// `pair` selects the `--intergenic` pair format with `offset`
/// (`--shift_zero`); `maxpos` is the sequence length or -1 for "no check".
pub fn parse(d: &[u8], maxpos: i32, pair: bool, offset: i32) -> Vec<Gene> {
    let mut genes = Vec::new();
    let n = d.len();
    let mut start = 0;
    let mut i = 0;
    while i < n {
        while i < n && d[i] != b'\n' {
            i += 1;
        }
        let line = &d[start..i];
        genes.push(if pair { parse_pair_line(line, offset) } else { parse_lst_line(line, maxpos) });
        // the next line begins at this '\n'
        start = i;
        i += 1;
    }
    genes
}

/// Loads and parses a `.lst` gene list (0x409a66).  `None` if the core
/// loader fails (it prints its own message).
pub fn load(name: &[u8], maxpos: i32) -> Option<Vec<Gene>> {
    io::load_file(name).map(|d| parse(&d, maxpos, false, 0))
}

/// Reads `--geneset` and loads it (0x409afe).  The getter's result is
/// ignored by the original: a missing or bad key prints its getter message,
/// the name stays empty and the loader then prints `error: file name is
/// missing` and fails.
pub fn load_param(p: &mut P, maxpos: i32) -> Option<Vec<Gene>> {
    let name = p.str_req("--geneset").unwrap_or_default();
    load(&name, maxpos)
}
