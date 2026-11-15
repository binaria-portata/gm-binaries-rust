//! gmhmmp external evidence (-e): GFF reader 0x41f5d0/0x41ec60/0x41eb00,
//! grouping 0x4236c0, per record 0x421690 (range), 0x421280 (M mask),
//! 0x421570 (N), 0x421450 (R), 0x423d60 (F), 0x4240e0 (C) with 0x40ce80,
//! 0x40cd10, 0x40c8a0. Port of re-notes/gmhmmp/input.md section 5.
//!
//! The evidence GFF steers predictions per sequence (GFF column 1 = FASTA
//! id). Types used: M (mask: bases become N), N (noncoding: remove starts
//! whose ORF reaches into the region), R (remove starts of ORFs whose stop
//! lies in the region), F and C (coding segments: validated against the
//! codon image, alternative starts removed, the stop recorded in
//! [`Maps`] for the DP and the output). H, T, E are read but unused.
//! Not used by GeneMark-ETP; the DP side of it is read from the binary
//! only (viterbi.md [R]).
//!
//! Coordinates: GFF start/end are 1-based; codon image index k is the
//! codon ENDING at 0-based base k, so the codon starting at 1-based s ends
//! at index s+1.

use super::num::{strtod, strtol_i32};
use std::collections::BTreeMap;

/// One grouped GFF record (0x20 bytes in the binary).
#[derive(Debug, Clone)]
pub struct EvRec {
    /// type code M=1 N=2 C=3 F=4 R=5 H=6 T=7 E=8; set to 0 when the record
    /// lies outside its sequence (permanently)
    pub code: u8,
    /// b'+', b'-' or b'.'
    pub strand: u8,
    /// 1-based inclusive
    pub start: i32,
    pub end: i32,
    /// `ID=` attribute value ("" if none)
    pub id: String,
    /// GFF score column ("." = 0.0); becomes the site's `ev` in the DP
    pub score: f64,
}

/// Per-sequence maps (+0x98/+0xc8 ids, +0xf8/+0x128 records), keyed by
/// the codon image index of the recorded stop codon (forward stop for '+',
/// reverse stop for '-').
#[derive(Default)]
pub struct Maps {
    pub ids_plus: BTreeMap<i64, Vec<String>>,
    pub ids_minus: BTreeMap<i64, Vec<String>>,
    pub rec_plus: BTreeMap<i64, EvRec>,
    pub rec_minus: BTreeMap<i64, EvRec>,
}

/// All selected GFF records grouped by sequence id, in file order.
pub struct Evidence {
    pub by_seq: BTreeMap<String, Vec<EvRec>>,
}

/// Type letter (already upper-cased) to code (0x4224e0).
fn type_code(t: &str) -> Option<u8> {
    match t {
        "M" => Some(1),
        "N" => Some(2),
        "C" => Some(3),
        "F" => Some(4),
        "R" => Some(5),
        "H" => Some(6),
        "T" => Some(7),
        "E" => Some(8),
        _ => None,
    }
}

/// 0x41eb00: value of the first "ID=" anywhere in the attributes (so
/// "geneID=" matches), from the first non-space up to ';', trailing spaces
/// removed.
fn attr_id(a: &str) -> String {
    let p = match a.find("ID=") {
        Some(p) => p,
        None => return String::new(),
    };
    let b = a.as_bytes();
    let mut st = p + 3;
    while st < b.len() && b[st] == b' ' {
        st += 1;
    }
    if st >= b.len() {
        return String::new();
    }
    let v = match a[st..].find(';') {
        Some(e) => &a[st..st + e],
        None => &a[st..],
    };
    v.trim_end_matches(' ').to_string()
}

/// One parsed GFF line (fields used later).
struct Line {
    seqid: String,
    typ: String,
    start: i32,
    end: i32,
    score: f64,
    strand: u8,
    attr: String,
}

/// 0x41ec60: Ok(None) = silently skipped (blank or '#'); Err(msg) = logged
/// and ignored. Fields are tokenised and checked one by one, in order;
/// runs of TABs collapse (empty fields vanish), spaces are not separators,
/// a 10th field or a trailing TAB is an error.
fn parse_line(line: &str) -> Result<Option<Line>, &'static str> {
    let b = line.as_bytes();
    if !b.iter().any(|&c| c != b' ' && c != b'\t') {
        return Ok(None);
    }
    if b[0] == b'#' {
        return Ok(None);
    }
    if b[0] == b'\t' {
        return Err("ignored GFF line; line starts with TAB");
    }
    let mut pos = 0usize;
    let mut en: Option<usize> = None;
    let mut l = Line { seqid: String::new(), typ: String::new(), start: 0, end: 0, score: 0.0, strand: 0, attr: String::new() };
    for f in 1..=9 {
        let st = if pos <= b.len() { (pos..b.len()).find(|&i| b[i] != b'\t') } else { None };
        let st = match st {
            Some(s) => s,
            None => return Err("ignored GFF line; 9 TAB delimited fields are expected"),
        };
        en = line[st..].find('\t').map(|x| x + st);
        let tok = match en {
            Some(e) => &line[st..e],
            None => &line[st..],
        };
        pos = en.unwrap_or(b.len() + 1);
        match f {
            1 => l.seqid = tok.to_string(),
            3 => l.typ = tok.to_string(),
            4 => {
                l.start = strtol_i32(tok);
                if l.start < 1 {
                    return Err("ignored GFF line; start position less then 1");
                }
            }
            5 => {
                l.end = strtol_i32(tok);
                if l.end < 1 {
                    return Err("ignored GFF line; end position less then 1");
                }
                if l.end < l.start {
                    return Err("ignored GFF line; end position is less than start position");
                }
            }
            6 => l.score = if tok == "." { 0.0 } else { strtod(tok) },
            7 => {
                l.strand = b[st];
                if !(tok.len() == 1 && "+-.".contains(tok)) {
                    return Err("ignored GFF line; unexpected value for strand");
                }
            }
            8 => {
                if !(tok.len() == 1 && "012.".contains(tok)) {
                    return Err("ignored GFF line; unexpected value for phase");
                }
            }
            9 => l.attr = tok.to_string(),
            _ => {}
        }
    }
    if en.is_some() {
        return Err("ignored GFF line; 9 TAB delimited fields are expected");
    }
    Ok(Some(l))
}

/// 0x41f5d0 + 0x4236c0: read and group the GFF. Log lines are appended to
/// `log` (written by main to -L or, with -v, stdout). C, F, R, E records
/// need a strand.
pub fn load(path: &str, log: &mut Vec<String>) -> Result<Evidence, String> {
    let data = std::fs::read(path).map_err(|_| format!("error on open GFF file: {path}\n"))?;
    let text = String::from_utf8_lossy(&data);
    let mut lines: Vec<&str> = text.split('\n').collect();
    if lines.last() == Some(&"") {
        lines.pop();
    }
    let mut sel: Vec<Line> = Vec::new();
    for (n, line) in lines.iter().enumerate() {
        match parse_line(line) {
            Ok(None) => {}
            Err(m) => {
                log.push(m.to_string());
                log.push(line.to_string());
            }
            Ok(Some(mut l)) => {
                l.typ = l.typ.to_ascii_uppercase();
                if type_code(&l.typ).is_none() {
                    log.push(format!("ignored GFF line; unsupported type found on line: {}", n + 1));
                    log.push(line.to_string());
                    continue;
                }
                sel.push(l);
            }
        }
    }
    log.push(format!("GFF lines in file: {}", lines.len()));
    log.push(format!("GFF lines selected: {}", sel.len()));
    let mut by_seq: BTreeMap<String, Vec<EvRec>> = BTreeMap::new();
    for l in sel {
        let code = type_code(&l.typ).unwrap();
        if matches!(code, 3 | 4 | 5 | 8) && l.strand == b'.' {
            log.push("GFF record ignored, strand must be defined for this feature".into());
            continue;
        }
        by_seq.entry(l.seqid.clone()).or_default().push(EvRec {
            code,
            strand: l.strand,
            start: l.start,
            end: l.end,
            id: attr_id(&l.attr),
            score: l.score,
        });
    }
    Ok(Evidence { by_seq })
}

/// 0x421690 + 0x421280 inside the reader (before copy and gap split):
/// records ending past the sequence get code 0 (logged); M records set
/// their bases to code 4. `codes` = the record's freshly encoded codes.
pub fn reader_step(recs: &mut [EvRec], codes: &mut [u8], log: &mut Vec<String>) {
    let l = codes.len() as i32;
    for r in recs.iter_mut() {
        if l < r.end {
            r.code = 0;
            log.push("GFF record ignored, end position is out of sequence".into());
        }
    }
    for r in recs.iter() {
        if r.code == 1 {
            if l < r.end {
                log.push("GFF record ignored, end position is out of sequence".into());
            } else {
                for i in (r.start - 1)..r.end {
                    codes[i as usize] = 4;
                }
            }
        }
    }
}

/// 0x40cd10 (R records, also the first step of N): `s`, `e` = GFF 1-based
/// start/end. Forward: for every forward stop ending at k in [s+1,
/// min(L-1, e+1)], clear the forward starts upstream in its frame down to
/// the previous stop (not below s+1). Reverse: for every reverse stop in
/// [s-1, e-1], clear the reverse starts downstream in its frame up to the
/// next reverse stop (not past e-1). Invalid ranges print
/// "error in external data".
fn clear_starts_of_region_stops(img: &mut [u8], s: i64, e: i64, out: &mut String) {
    let l = img.len() as i64;
    if s < 1 || e < s || l < e {
        out.push_str(&format!("error in external data: {s} {e}\n"));
        return;
    }
    let lo = s + 1;
    let hi = (l - 1).min(e + 1);
    for k in lo..=hi {
        if img[k as usize] & 0x10 != 0 && k - 3 >= lo {
            let mut j = k - 3;
            let mut b = img[j as usize];
            while b & 0x10 == 0 {
                if b & 0x20 != 0 {
                    img[j as usize] ^= 0x20;
                }
                j -= 3;
                if j < lo {
                    break;
                }
                b = img[j as usize];
            }
        }
    }
    for k in (s - 1)..e {
        if img[k as usize] & 4 != 0 && k + 3 <= e - 1 {
            let mut j = k + 3;
            let mut b = img[j as usize];
            while b & 4 == 0 {
                if b & 8 != 0 {
                    img[j as usize] ^= 8;
                }
                j += 3;
                if j > e - 1 {
                    break;
                }
                b = img[j as usize];
            }
        }
    }
}

/// 0x40ce80 (N records): [`clear_starts_of_region_stops`], then from the
/// three codon ends s+1, s, s-1 walk left in steps of 3 clearing forward
/// starts until a forward stop, and from e-1, e, e+1 walk right clearing
/// reverse starts until a reverse stop: no ORF reaching into the region
/// keeps a start.
fn clear_starts_noncoding(img: &mut [u8], s: i64, e: i64, out: &mut String) {
    let l = img.len() as i64;
    clear_starts_of_region_stops(img, s, e, out);
    for p in [s + 1, s, s - 1] {
        let mut j = p;
        while j >= 0 {
            if j < l {
                if img[j as usize] & 0x20 != 0 {
                    img[j as usize] ^= 0x20;
                }
                if img[j as usize] & 0x10 != 0 {
                    break;
                }
            }
            j -= 3;
        }
    }
    for p in [e - 1, e, e + 1] {
        let mut j = p;
        while j < l {
            if j >= 0 {
                if img[j as usize] & 8 != 0 {
                    img[j as usize] ^= 8;
                }
                if img[j as usize] & 4 != 0 {
                    break;
                }
            }
            j += 3;
        }
    }
}

/// 0x40c8a0: check a CDS segment against the codon image, remove the
/// alternative starts of its ORF and return the codon image index of its
/// stop, or -1 with a warning (printf to stdout). As in the original: "stop not
/// found" on '-' returns 0, which callers accept as a stop index.
///
/// * `s0`, `e0`: 0-based first and last base of the GFF segment.
/// * `start_only`: check only the first codon; `stop_only`: only the last.
/// * '+': no in-frame stop before the last codon; the stop is searched from
///   the last codon onwards; the start is the nearest in-frame start at or
///   before the first codon (not across a stop), else the first one inside
///   the segment; starts strictly between start and stop are cleared.
/// * '-': mirror image (reverse stop on the left, reverse starts to the
///   right); its internal-stop message says "for" and prints a different
///   position (input.md 5.4).
#[allow(clippy::too_many_arguments)]
fn check_cds_segment(img: &mut [u8], s0: i64, e0: i64, strand: u8, sid: &str, start_only: bool, stop_only: bool, out: &mut String) -> i64 {
    let l = img.len() as i64;
    // 1-based start/end for the messages
    let start1 = s0 + 1;
    let end1 = e0 + 1;
    // one past the checked part (0-based); also reused as the position in
    // the '-' internal-stop message
    let mut seg_end = end1;
    let (mut seg_lo, mut seg_hi) = (s0, e0);
    if start_only {
        seg_hi = s0 + 2;
        seg_end = s0 + 3;
    }
    let last_codon = seg_hi - 2;
    if stop_only {
        seg_lo = last_codon;
    }
    let d = seg_end - seg_lo;
    let sc = strand as char;
    let at = |img: &[u8], i: i64| -> u8 { if i >= 0 && i < l { img[i as usize] } else { 0 } };
    if d - (d / 3) * 3 != 0 {
        out.push_str(&format!("Warning, segment was ignored. Length of coding sequence is not divisible by three: {sid}\t{start1}\t{end1}\t{sc}\n"));
        return -1;
    }
    if strand == b'+' {
        // from here seg_lo is the codon END index of the first codon
        seg_lo += 2;
        let mut j = seg_lo;
        while j < last_codon {
            if at(img, j) & 0x10 != 0 {
                out.push_str(&format!("Warning, segment was ignored. Internal stop codon was found at position {} : {sid}\t{start1}\t{end1}\t+\n", j + 1));
                return -1;
            }
            j += 3;
        }
        if seg_hi < l {
            let mut k = seg_hi;
            let mut found = true;
            while at(img, k) & 0x10 == 0 {
                k += 3;
                if k >= l {
                    found = false;
                    break;
                }
            }
            if found {
                let mut st: Option<i64> = None;
                if seg_lo >= 0 && at(img, seg_lo) & 0x10 == 0 {
                    let mut j = seg_lo;
                    let mut ok = true;
                    while at(img, j) & 0x20 == 0 {
                        j -= 3;
                        if j < 0 || at(img, j) & 0x10 != 0 {
                            ok = false;
                            break;
                        }
                    }
                    st = if ok { Some(j) } else { None };
                }
                if st.is_none() && seg_lo <= seg_hi {
                    let mut j = seg_lo;
                    let mut ok = true;
                    while at(img, j) & 0x20 == 0 {
                        j += 3;
                        if j > seg_hi {
                            ok = false;
                            break;
                        }
                    }
                    st = if ok { Some(j) } else { None };
                }
                let st = match st {
                    Some(x) => x,
                    None => {
                        out.push_str(&format!("Warning, segment was ignored. Start codon was not found: {sid}\t{start1}\t{end1}\t+\n"));
                        return -1;
                    }
                };
                let mut m = st + 3;
                while m < k {
                    if at(img, m) & 0x20 != 0 {
                        img[m as usize] ^= 0x20;
                    }
                    m += 3;
                }
                return k;
            }
        }
        out.push_str(&format!("Warning, segment was ignored. Stop codon was not found: {sid}\t{start1}\t{end1}\t+\n"));
        -1
    } else if strand == b'-' {
        if seg_lo + 4 < seg_hi {
            let mut k = seg_hi;
            let mut internal = None;
            if at(img, seg_hi) & 4 == 0 {
                loop {
                    let prev = k;
                    k -= 3;
                    if k <= seg_lo + 4 {
                        break;
                    }
                    if at(img, k) & 4 != 0 {
                        seg_end = prev - 2;
                        internal = Some(k);
                        break;
                    }
                }
            } else {
                internal = Some(k);
            }
            if internal.is_some() {
                out.push_str(&format!("Warning, segment was ignored. Internal stop codon was found at position {seg_end} for {sid}\t{start1}\t{end1}\t-\n"));
                return -1;
            }
        }
        let mut k = seg_lo + 2;
        if k < 2 {
            out.push_str(&format!("Warning, segment was ignored. Stop codon was not found: {sid}\t{start1}\t{end1}\t-\n"));
            return 0;
        }
        while at(img, k) & 4 == 0 {
            k -= 3;
            if k < 2 {
                out.push_str(&format!("Warning, segment was ignored. Stop codon was not found: {sid}\t{start1}\t{end1}\t-\n"));
                return 0;
            }
        }
        let mut st: Option<i64> = None;
        if seg_hi < l && at(img, seg_hi) & 4 == 0 {
            let mut j = seg_hi;
            let mut ok = true;
            while at(img, j) & 8 == 0 {
                j += 3;
                if j >= l || at(img, j) & 4 != 0 {
                    ok = false;
                    break;
                }
            }
            st = if ok { Some(j) } else { None };
        }
        if st.is_none() && seg_lo + 1 < seg_hi {
            let mut j = seg_hi;
            let mut ok = true;
            while at(img, j) & 8 == 0 {
                j -= 3;
                if j <= seg_lo + 1 {
                    ok = false;
                    break;
                }
            }
            st = if ok { Some(j) } else { None };
        }
        let st = match st {
            Some(x) => x,
            None => {
                out.push_str(&format!("Warning, segment was ignored. Start codon was not found: {sid}\t{start1}\t{end1}\t-\n"));
                return -1;
            }
        };
        let mut m = st - 3;
        while m > k {
            if at(img, m) & 8 != 0 {
                img[m as usize] ^= 8;
            }
            m -= 3;
        }
        k
    } else {
        out.push_str(&format!("Warning, segment was ignored. Strand was not specified: {sid}\t{start1}\t{end1}\t{sc}\n"));
        -1
    }
}

/// After the codon image (main, only with -e): N, R, (clear maps), F, C
/// records in file order, editing the codon image `img` in place.
/// F: start and stop codon checked separately, both indices recorded (only
/// if `start < end - 1`). C: whole segment; an existing record is replaced
/// only by one reaching further (+: smaller start, -: larger end) and never
/// if it is an F record. `out` = printf text (stdout).
pub fn image_step(recs: &[EvRec], img: &mut [u8], sid: &str, out: &mut String) -> Maps {
    for r in recs {
        if r.code == 2 {
            clear_starts_noncoding(img, r.start as i64, r.end as i64, out);
        }
    }
    for r in recs {
        if r.code == 5 {
            clear_starts_of_region_stops(img, r.start as i64, r.end as i64, out);
        }
    }
    let mut maps = Maps::default();
    for r in recs {
        if r.code != 4 || !(r.strand == b'+' || r.strand == b'-') || !(r.start < r.end - 1) {
            continue;
        }
        let (s0, e0) = (r.start as i64 - 1, r.end as i64 - 1);
        let (a, b) = if r.strand == b'+' {
            (check_cds_segment(img, s0, e0, b'+', sid, true, false, out), check_cds_segment(img, s0, e0, b'+', sid, false, true, out))
        } else {
            (check_cds_segment(img, s0, e0, b'-', sid, false, true, out), check_cds_segment(img, s0, e0, b'-', sid, true, false, out))
        };
        if a < 0 || b < 0 {
            continue;
        }
        let (ids, rec) = if r.strand == b'+' { (&mut maps.ids_plus, &mut maps.rec_plus) } else { (&mut maps.ids_minus, &mut maps.rec_minus) };
        ids.entry(a).or_default().push(r.id.clone());
        if a != b {
            ids.entry(b).or_default().push(r.id.clone());
            rec.insert(a, r.clone());
            rec.insert(b, r.clone());
        } else {
            rec.insert(a, r.clone());
        }
    }
    for r in recs {
        if r.code != 3 {
            continue;
        }
        let k = check_cds_segment(img, r.start as i64 - 1, r.end as i64 - 1, r.strand, sid, false, false, out);
        if k < 0 {
            continue;
        }
        if r.strand == b'+' {
            maps.ids_plus.entry(k).or_default().push(r.id.clone());
            let replace = match maps.rec_plus.get(&k) {
                None => true,
                Some(ex) => r.start < ex.start && ex.code != 4,
            };
            if replace {
                maps.rec_plus.insert(k, r.clone());
            }
        } else if r.strand == b'-' {
            maps.ids_minus.entry(k).or_default().push(r.id.clone());
            let replace = match maps.rec_minus.get(&k) {
                None => true,
                Some(ex) => ex.end < r.end && ex.code != 4,
            };
            if replace {
                maps.rec_minus.insert(k, r.clone());
            }
        }
    }
    maps
}
