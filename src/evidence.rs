//! External evidence (`-d`, "GeneMark.hmm plus"): 0x40cc40 (file loop),
//! 0x40c340 (line parser), 0x409d10 (marks). See re-notes/input.md section 5.
//!
//! Third stage of the pipeline (optional). Hints from a GFF-like file
//! (introns, exons, start/stop codons, splice sites, masked regions) are
//! turned into one u32 bitmask per base. The marks are copied into every
//! site-map record (`Rec::mark`) and used in two places:
//! * the site track adds +100 to a signal with matching point evidence and
//!   -100 to a non-stop signal inside an intron hint ([`crate::tracks`]);
//! * the coding tracks add the `-k` penalty per base carrying [`MASK`].
//!
//! Bit layout (names from the dead debug dumper 0x40b120):
//! point marks use the site type bits of [`crate::sitemap`]
//! (0x1 DirStart, 0x2 DirEnd, 0x4 DirDn, 0x8 DirAc, 0x10 RevStart,
//! 0x20 RevEnd, 0x40 RevDn, 0x80 RevAc); region marks are 0x1000 CleanMark,
//! 0x2000 CodMark, 0x4000 NonMark, 0x8000 IntronMark, 0x10000 InterMark,
//! 0x20000 MaskMark, combined into the region constants below.
//!
//! Messages go to stdout; a bad line is reported and ignored, never fatal.
//! Only an unopenable file is fatal.

/// Coding region: CleanMark | CodMark.
pub const COD: u32 = 0x3000;
/// Intron region: CleanMark | NonMark | IntronMark.
pub const INTRON: u32 = 0xd000;
/// Non-coding region (UTR, non-coding RNA): CleanMark | NonMark.
pub const NONCOD: u32 = 0x5000;
/// Intergenic region: CleanMark | NonMark | InterMark.
pub const INTERGENIC: u32 = 0x15000;
/// Soft-masked region (`mask` lines, used with `-k`): MaskMark.
pub const MASK: u32 = 0x20000;

/// The per-base mark array being filled (seq+0x28, u32[L]).
struct Marks<'a>(&'a mut [u32]);

impl Marks<'_> {
    /// OR `v` into positions `a..=b` (0-based, inclusive; empty if a > b).
    /// The range is clipped to the sequence (the caller has already checked
    /// that the feature lies inside it).
    fn range(&mut self, a: i64, b: i64, v: u32) {
        let l = self.0.len() as i64;
        let mut i = a.max(0);
        while i <= b && i < l {
            self.0[i as usize] |= v;
            i += 1;
        }
    }
    /// OR `v` into one position. The original has no bounds check here:
    /// `Internal+` starting at base 1 writes marks[-1], `Initial+` ending at
    /// L writes marks[L]. marks[-1] is the upper half of the heap chunk's
    /// size field: the original aborts ("double free or corruption") when it
    /// frees the array at the end, after the output was written (for arrays
    /// of 128 KiB or more, i.e. L >= 32768, it survives). marks[L] lands in
    /// the allocation's slack (harmless) or, for L = 2 (mod 4), on the next
    /// chunk's size field (SIGSEGV before any output). The port drops such
    /// out-of-range point marks, which is what the surviving runs produce,
    /// and warns once where the original would terminate.
    fn point(&mut self, p: i64, v: u32) {
        if p >= 0 && (p as usize) < self.0.len() {
            self.0[p as usize] |= v;
        } else if original_not_handled(self.0.len() as u64, p) {
            // only where the original terminates; elsewhere its write is
            // harmless and the output identical, so nothing is said
            crate::warn(
                "-d: a hint touches the first or last base of the sequence so that one of its splice site or \
                 codon marks falls just outside the sequence; that mark is ignored, the rest of the hint is \
                 used (not handled by the original)",
            );
        }
    }
}

/// Does the original terminate after writing marks[p] (p = -1 or L) of its
/// `malloc(4 * L)` array (glibc 2.17 malloc, 64-bit)? See [`Marks::point`].
fn original_not_handled(l: u64, p: i64) -> bool {
    let req = 4 * l;
    // request2size: 8 bytes overhead, 16-byte alignment, minimum 32
    let chunk = ((req + 8 + 15) & !15).max(32);
    if chunk >= 128 * 1024 {
        // mmap()ed chunk: whole pages, 16-byte header; marks[-1] only
        // changes high bits of the size field (survives); marks[L] faults
        // if it lies past the mapping
        let mapped = (chunk + 8 + 4095) & !4095;
        return p >= 0 && req + 4 > mapped - 16;
    }
    // heap chunk: usable size = chunk - 8; marks[-1] is the high half of
    // the size field (abort at free), marks[L] past the usable size hits
    // the next chunk's size field
    p < 0 || req + 4 > chunk - 8
}

/// 0x409d10: apply one feature with 0-based inclusive coordinates `s..=e`
/// and strand byte `strand` ('+', '-' or '.'). Returns false for a
/// feature/strand combination the original does not know.
///
/// For exons the region mark covers the coding bases between the signals
/// and the point marks sit where the corresponding map record would be
/// (for example `Initial +`: start codon A at `s`, donor G at `e+1`).
/// The checks fall through in the original's order: strand-specific
/// structural features first, then strand-independent regions, then single
/// sites, then `mask`.
fn apply_feature(m: &mut Marks, feature: &str, strand: u8, s: i64, e: i64) -> bool {
    match (strand, feature) {
        (b'+', "Single") => { m.range(s + 1, e - 1, COD); m.point(s, 1); m.point(e, 2); return true; }
        (b'+', "Initial") => { m.range(s + 1, e, COD); m.point(s, 1); m.point(e + 1, 4); return true; }
        (b'+', "Internal") => { m.range(s, e, COD); m.point(s - 1, 8); m.point(e + 1, 4); return true; }
        (b'+', "Terminal") => { m.range(s, e - 1, COD); m.point(s - 1, 8); m.point(e, 2); return true; }
        (b'+', "intron") => { m.range(s + 1, e - 1, INTRON); m.point(s, 4); m.point(e, 8); return true; }
        (b'+', "Part_don") => { m.range(s + 1, e, COD); m.point(e + 1, 4); return true; }
        (b'+', "Part_acc") => { m.range(s, e - 1, COD); m.point(s - 1, 8); return true; }
        (b'+', "Part_ini") => { m.range(s + 1, e - 1, COD); m.point(s, 1); return true; }
        (b'+', "Part_ter") => { m.range(s + 1, e - 1, COD); m.point(e, 2); return true; }
        (b'+', "Part_exon") => { m.range(s + 1, e - 1, COD); return true; }
        (b'-', "Single") => { m.range(s + 1, e - 1, COD); m.point(s, 0x20); m.point(e, 0x10); return true; }
        (b'-', "Initial") => { m.range(s, e - 1, COD); m.point(s - 1, 0x40); m.point(e, 0x10); return true; }
        (b'-', "Internal") => { m.range(s, e, COD); m.point(s - 1, 0x40); m.point(e + 1, 0x80); return true; }
        (b'-', "Terminal") => { m.range(s + 1, e, COD); m.point(s, 0x20); m.point(e + 1, 0x80); return true; }
        (b'-', "intron") => { m.range(s + 1, e - 1, INTRON); m.point(s, 0x80); m.point(e, 0x40); return true; }
        (b'-', "Part_don") => { m.range(s, e - 1, COD); m.point(s - 1, 0x40); return true; }
        (b'-', "Part_acc") => { m.range(s + 1, e, COD); m.point(e + 1, 0x80); return true; }
        (b'-', "Part_ini") => { m.range(s + 1, e - 1, COD); m.point(e, 0x10); return true; }
        (b'-', "Part_ter") => { m.range(s + 1, e - 1, COD); m.point(s, 0x20); return true; }
        (b'-', "Part_exon") => { m.range(s + 1, e - 1, COD); return true; }
        // 0x33 = DirStart | DirEnd | RevStart | RevEnd: all four codon-site
        // bits on the two bases flanking the intergenic hint
        (b'.', "Intergenic") => { m.range(s, e, INTERGENIC); m.point(s - 1, 0x33); m.point(e + 1, 0x33); return true; }
        _ => {}
    }
    // strand-independent regions
    match feature {
        "UTR_Intron" | "UTR_Exon" | "Non_Coding" => { m.range(s, e, NONCOD); return true; }
        "Part_intron" => { m.range(s, e, INTRON); return true; }
        "Part_intergenic" => { m.range(s, e, INTERGENIC); return true; }
        _ => {}
    }
    // single sites (only with strand + or -)
    let pt = match (strand, feature) {
        (b'+', "Donor") => Some((s, 4)),
        (b'+', "Acceptor") => Some((e, 8)),
        (b'+', "start_codon") => Some((s, 1)),
        (b'+', "stop_codon") => Some((e, 2)),
        (b'-', "Donor") => Some((e, 0x40)),
        (b'-', "Acceptor") => Some((s, 0x80)),
        (b'-', "start_codon") => Some((e, 0x10)),
        (b'-', "stop_codon") => Some((s, 0x20)),
        _ => None,
    };
    if let Some((p, v)) = pt {
        m.point(p, v);
        return true;
    }
    if feature == "mask" {
        m.range(s, e, MASK);
        return true;
    }
    false
}

/// 0x40c340: parse one line. `Ok(())` also for skipped lines (blank,
/// comment, other sequence id); `Err(message)` for a rejected line (the
/// caller adds the "line ignored" warning).
fn parse_line(line: &str, seqid: &str, m: &mut Marks) -> Result<(), String> {
    let b = line.as_bytes();
    // blank or only spaces/tabs: skipped. A line that is just "\r" is not
    // blank here and fails the column count, as in the original.
    let mut pos = match b.iter().position(|&c| c != b' ' && c != b'\t') {
        Some(p) => p,
        None => return Ok(()),
    };
    if b[pos] == b'#' {
        return Ok(());
    }
    // Tab-separated tokens; runs of tabs collapse (empty fields vanish),
    // spaces are not separators, at most 9 tokens are taken (a 10th column
    // is ignored, the 9th token ends at the next tab).
    let mut toks: Vec<&str> = Vec::with_capacity(9);
    while toks.len() < 9 {
        let end = b[pos..].iter().position(|&c| c == b'\t').map(|x| x + pos).unwrap_or(b.len());
        toks.push(&line[pos..end]);
        match b[end..].iter().position(|&c| c != b'\t') {
            Some(q) => pos = end + q,
            None => break,
        }
    }
    if toks.len() < 9 {
        return Err(format!("Data format error, 9 columns is expected but {} was found in: {line}\n", toks.len()));
    }
    // lines for other sequences are silently ignored (exact, case-sensitive)
    if toks[0] != seqid {
        return Ok(());
    }
    // 1-based inclusive GFF coordinates, parsed like (int)strtol
    let s = crate::model::strtol(toks[3]) as i32 as i64;
    let e = crate::model::strtol(toks[4]) as i32 as i64;
    let l = m.0.len() as i64;
    if e < s || s < 1 || e > l {
        return Err(format!("Data format error, check start and/or end coordinates: {line}\n"));
    }
    let strand = toks[6].as_bytes().first().copied().unwrap_or(0);
    if strand != b'+' && strand != b'-' && strand != b'.' {
        return Err(format!("Data format error, unexpected symbol found in strand field: {line}\n"));
    }
    // the only feature-name normalisations of the original
    let feature = match toks[2] {
        "Intron" => "intron",
        "Start_codon" => "start_codon",
        "Stop_codon" => "stop_codon",
        x => x,
    };
    if !apply_feature(m, feature, strand, s - 1, e - 1) {
        return Err(format!("Data format error, unexpected feature {feature} in : {line}\n"));
    }
    Ok(())
}

/// 0x40cc40: read the evidence file for a sequence with id `seqid` and
/// length `len`. Returns the per-base marks; all messages are appended to
/// `out` (the caller prints them). `Err(())` only if the file cannot be
/// opened (fatal in main).
pub fn read(path: &str, seqid: &str, len: usize, out: &mut String) -> Result<Vec<u32>, ()> {
    // A directory: std::ifstream opens it, getline then fails without
    // setting eof, and the original's `while (!eof)` loop never ends (it
    // hangs). The port reports it like a file that cannot be opened.
    let data = match std::fs::read(crate::sequence::os_path(path)) {
        Ok(d) => d,
        Err(_) => {
            out.push_str(&format!("gmhmme3 : error, opening file {path}\n"));
            return Err(());
        }
    };
    // raw bytes, invalid UTF-8 kept in escaped form for the messages
    let text = crate::sequence::decode(&data);
    let mut marks = vec![0u32; len];
    let mut m = Marks(&mut marks);
    // getline(.., '\n'): a trailing '\r' stays in the line
    let mut lines: Vec<&str> = text.split('\n').collect();
    // getline: a final newline does not produce an extra empty line
    if text.ends_with('\n') {
        lines.pop();
    }
    for line in lines {
        if let Err(msg) = parse_line(line, seqid, &mut m) {
            out.push_str(&msg);
            out.push_str(&format!("gmhmme3 : warning, file {path} line ignored : {line}\n"));
        }
    }
    Ok(marks)
}
