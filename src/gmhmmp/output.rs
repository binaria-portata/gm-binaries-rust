//! gmhmmp writers: headers (0x4251e0, 0x425030), per-record block
//! (0x413a70, 0x4133f0, 0x412fc0, 0x413510), -A/-D (0x412cf0, 0x412a60,
//! 0x412b40), -K adjustment (0x412210). See re-notes/gmhmmp/output.md.
//!
//! Texts are built per record and written by main as one block each
//! (main file, -A, -D), every block followed by one extra newline. Gene
//! coordinates are 1-based inclusive. Two gene numberings exist: the
//! running id over the whole file (`Gene::id`) and, with `-b`, a separate
//! counter of the kept genes for the LST/GFF lines.

use super::cli::Settings;
use super::evidence::{EvRec, Maps};
use super::num::c_f;
use super::predict::Gene;

pub const VERSION: &str = "3.20";

/// File header (LST 0x4251e0 or GFF 0x425030). `date` is ctime output
/// (ends in '\n'); `single_name` is the NAME of a single model (meta models
/// print it per record instead, see [`model_info`]).
pub fn header(s: &Settings, date: &str, single_name: Option<&str>) -> String {
    let rbs = if s.rbs { "true" } else { "false" };
    let mut o = String::new();
    if s.format == b'G' {
        o.push_str("##gff-version 2\n");
        o.push_str(&format!("##source-version GeneMark.hmm_PROKARYOTIC {VERSION}\n"));
        o.push_str(&format!("##date: {date}"));
        o.push_str(&format!("# Sequence file name: {}\n", s.seq));
        o.push_str(&format!("# Model file name: {}\n", s.model));
        o.push_str(&format!("# RBS: {rbs}\n"));
        if let Some(n) = single_name {
            o.push_str(&format!("# Model information: {n}\n"));
        }
    } else {
        o.push_str(&format!("GeneMark.hmm PROKARYOTIC (Version {VERSION})\n"));
        o.push_str(&format!("Date: {date}"));
        o.push_str(&format!("Sequence file name: {}\n", s.seq));
        o.push_str(&format!("Model file name: {}\n", s.model));
        o.push_str(&format!("RBS: {rbs}\n"));
        if let Some(n) = single_name {
            o.push_str(&format!("Model information: {n}\n"));
        }
    }
    o
}

/// Per-record model line of a meta model (0x424ef0).
pub fn model_info(s: &Settings, name: &str) -> String {
    if s.format == b'G' {
        format!("# Model information: {name}\n")
    } else {
        format!("Model information: {name}\n")
    }
}

/// 0x412210 (`-r -k -K`): RBS position instead of spacer, computed on
/// bytes as in the binary, so the printed value is the true position
/// modulo 256 as a signed byte (e.g. 266 -> 10). `w` = RBSS.
pub fn k_positions(genes: &mut [Gene], w: i32) {
    for g in genes.iter_mut() {
        let v = if g.strand == b'+' {
            (g.left as u8).wrapping_sub(g.spacer as u8).wrapping_sub(w as u8)
        } else {
            (g.right as u8).wrapping_add(g.spacer as u8).wrapping_add(w as u8)
        };
        g.spacer = v as i8;
    }
}

/// Genetic code, index 16*b0 + 4*b1 + b2 (A=0 C=1 G=2 T=3); entry 56
/// (TGA) becomes 'W' if the model's TGA_ is 0.
const CODE: &[u8; 64] = b"KNKNTTTTRSRSIIMIQHQHPPPPRRRRLLLLEDEDAAAAGGGGVVVV*Y*YSSSS*CWCLFLF";

/// Sequence body of -A/-D: a newline after every 60th character and one at
/// the end (so a length divisible by 60 ends in a blank line, as in the
/// binary).
fn wrap(o: &mut String, body: &[u8]) {
    for (i, &c) in body.iter().enumerate() {
        o.push(c as char);
        if i != 0 && (i + 1) % 60 == 0 {
            o.push('\n');
        }
    }
    o.push('\n');
}

/// Output texts of one record.
pub struct Block {
    /// main file (-o) text
    pub main: String,
    /// -A protein text
    pub aa: String,
    /// -D nucleotide text
    pub nt: String,
    /// messages to stdout ("Internal error!", evidence join warnings)
    pub stdout: String,
}

/// 0x413a70 for one record. `best_counter` is the global -b counter.
///
/// * `genes` get `maxlo` = max logodd of the record (start -1000, strict
///   `>`).
/// * `-b`: in LST/GFF a gene with `logodd < max - 0.1` leaves an EMPTY line,
///   kept genes are numbered by `best_counter`; in -A/-D the threshold is
///   `max - 0.01` and the headers use `Gene::id` (both thresholds as in the
///   binary).
/// * A record without genes still gets the 4 LST header lines.
/// * `orig` = codes before gap split with 5/6 at resolved R/Y; reverse
///   genes map every code >= 4 to 5, so N prints as 'R' in -D.
/// * `tga_prob` = the model's TGA_ (translation table switch).
/// * `ev`: evidence maps, only with -e (extra columns, map entries may be
///   inserted as by C++ `operator[]`).
#[allow(clippy::too_many_arguments)]
pub fn record_block(s: &Settings, defline: &str, genes: &mut [Gene], orig: &[u8], tga_prob: f64, best_counter: &mut i32, mut ev: Option<&mut Maps>) -> Block {
    let mut b = Block { main: String::new(), aa: String::new(), nt: String::new(), stdout: String::new() };
    let lst = s.format != b'G';
    let defl_nogt = defline.get(1..).unwrap_or("");
    if lst {
        b.main.push_str(&format!("FASTA definition line: {defl_nogt}\n"));
        b.main.push_str("Predicted genes\n");
        if s.rbs && s.motif_info {
            if s.motif_pos {
                b.main.push_str("   Gene    Strand    LeftEnd    RightEnd       Gene     Class     RBS      RBS\n");
                b.main.push_str("    #                                         Length            position  score\n");
            } else {
                b.main.push_str("   Gene    Strand    LeftEnd    RightEnd       Gene     Class   Spacer     RBS\n");
                b.main.push_str("    #                                         Length            Length    score\n");
            }
        } else {
            b.main.push_str("   Gene    Strand    LeftEnd    RightEnd       Gene     Class\n");
            b.main.push_str("    #                                         Length\n");
        }
    }
    // logodd max (computed by the caller into g.logodd); M starts at -1000
    let mut mx = -1000.0f64;
    for g in genes.iter() {
        if g.logodd > mx {
            mx = g.logodd;
        }
    }
    for g in genes.iter_mut() {
        g.maxlo = mx;
    }
    for g in genes.iter() {
        let len = g.right + 1 - g.left;
        let line = if s.best && g.logodd < mx - 0.1 {
            String::new()
        } else {
            let number = if s.best {
                *best_counter += 1;
                *best_counter
            } else {
                g.id
            };
            if lst {
                let l = format!("{}{}", g.lm as char, g.left);
                let r = format!("{}{}", g.rm as char, g.right);
                let mut x = format!("{:5} {:>8} {:>11} {:>11} {:12} {:>8}", number, g.strand as char, l, r, len, g.cls as char);
                if s.rbs && s.motif_info {
                    x.push_str(&format!(" {:8} {:8.4}", g.spacer as i32, g.rbs));
                }
                if !s.evidence.is_empty() {
                    // label of 0x412260
                    let lab = if g.evf == 0 { "__" } else if g.strand == b'+' { "_S" } else { "S_" };
                    x.push_str("   ");
                    x.push_str(lab);
                }
                x
            } else {
                let attr = format!(
                    "gene_id={}, length={}, gene_score={}, rbs_score={}, rbs_spacer={}, stop_enforced={}, start_codon={}, logodd={}",
                    number,
                    len,
                    c_f(g.gs),
                    c_f(g.rbs),
                    g.spacer as i32,
                    if g.evf != 0 { 'Y' } else { 'N' },
                    ((g.sflags & 3) + b'0') as char,
                    c_f(g.logodd)
                );
                format!("{}\t{}\t{}\t{}\t{}\t{}\t{}\t0\t{}", defl_nogt, "GeneMark.hmm", "CDS", g.left, g.right, c_f(g.gs), g.strand as char, attr)
            }
        };
        let mut line = line;
        if let Some(m) = ev.as_deref_mut() {
            if g.evf != 0 {
                evidence_columns(&mut line, g, m, &mut b.stdout);
            }
        }
        b.main.push_str(&line);
        b.main.push('\n');
    }
    if s.print_aa || s.print_nt {
        let mut table = *CODE;
        if tga_prob == 0.0 {
            table[56] = b'W';
        }
        for g in genes.iter() {
            if s.best && g.logodd < g.maxlo - 0.01 {
                continue;
            }
            let mut nt: Vec<u8> = orig[(g.left - 1) as usize..g.right as usize].to_vec();
            if g.strand == b'-' {
                nt.reverse();
                for c in nt.iter_mut() {
                    *c = if *c < 4 { 3 - *c } else { 5 };
                }
            }
            if s.print_aa {
                let mut aa = Vec::new();
                let mut i = 0;
                while i + 2 < nt.len() {
                    let c = &nt[i..i + 3];
                    if c.iter().any(|&x| x >= 4) {
                        aa.push(b'x');
                    } else {
                        aa.push(table[c[0] as usize * 16 + c[1] as usize * 4 + c[2] as usize]);
                    }
                    i += 3;
                }
                if aa.last() == Some(&b'*') {
                    aa.pop();
                }
                for &c in &aa {
                    if c == b'*' {
                        b.stdout.push_str("Internal error! Please report to code developers!");
                    }
                }
                b.aa.push_str(&format!("{}\tgene_{}|GeneMark.hmm|{}_aa|{}|{}|{}\t\n", defline, g.id, aa.len(), g.strand as char, g.left, g.right));
                wrap(&mut b.aa, &aa);
            }
            if s.print_nt {
                let body: Vec<u8> = nt.iter().map(|&c| b"ACGTNRY"[c as usize]).collect();
                b.nt.push_str(&format!("{}\tgene_{}|GeneMark.hmm|{}_nt|{}|{}|{}\n", defline, g.id, nt.len(), g.strand as char, g.left, g.right));
                wrap(&mut b.nt, &body);
            }
        }
    }
    b
}

/// Default-constructed map value (C++ `operator[]` on a missing key).
fn zero_rec() -> EvRec {
    EvRec { code: 0, strand: 0, start: 0, end: 0, id: String::new(), score: 0.0 }
}

/// join of 0x412e90 (operator[] of the id map inserts an empty vector)
fn join_ids(ids: &mut std::collections::BTreeMap<i64, Vec<String>>, key: i64, stdout: &mut String) -> String {
    let v = ids.entry(key).or_default();
    if v.is_empty() {
        stdout.push_str("vector is empty, check code\n");
        return String::new();
    }
    // assign while the result is empty, else append " " + id
    let mut r = String::new();
    for id in v.iter() {
        if r.is_empty() {
            r = id.clone();
        } else {
            r.push(' ');
            r.push_str(id);
        }
    }
    r
}

/// Evidence columns appended in 0x413a70 (only with -e and gene+0x20 != 0)
fn evidence_columns(line: &mut String, g: &Gene, m: &mut Maps, stdout: &mut String) {
    if g.strand == b'+' {
        let key = g.right - 1;
        let e = m.rec_plus.entry(key).or_insert_with(zero_rec).clone();
        line.push_str(&format!(" {:11} {:11}", e.start, e.end));
        if e.code == 4 {
            let j = join_ids(&mut m.ids_plus, key, stdout);
            line.push_str(&format!(" F {j}"));
            return;
        }
        let mut c_col = e.code == 3;
        if !c_col {
            let o = m.rec_minus.entry(g.left + 1).or_insert_with(zero_rec);
            c_col = o.code == b'H';
        }
        if c_col {
            let j = join_ids(&mut m.ids_plus, key, stdout);
            line.push_str(&format!(" C {j}"));
        }
    } else if g.strand == b'-' {
        let key = g.left + 1;
        let e = m.rec_minus.entry(key).or_insert_with(zero_rec).clone();
        line.push_str(&format!(" {:11} {:11}", e.start, e.end));
        if e.code == 4 {
            let j = join_ids(&mut m.ids_minus, key, stdout);
            line.push_str(&format!(" F {j}"));
            return;
        }
        if e.code == 3 || e.code == b'T' {
            let j = join_ids(&mut m.ids_minus, key, stdout);
            line.push_str(&format!(" C {j}"));
        }
    }
}
