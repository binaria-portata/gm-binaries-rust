//! `--gibbs` (0x455cc4): turns the output of an external Gibbs motif
//! sampler into the RBS block of a GeneMarkS-T model.  models.md §14.
//!
//! gmst.pl uses this route only with `--fixmotif 0`: `probuild --mkmod
//! --pre_start startseq` writes the pre-start sequences (all of the same
//! length L), `Gibbs3` finds one motif site of width W in each, and this
//! mode reads the sites and their positions and appends to the model:
//!
//! - `RBSS W`, `RBSM` the W x 4 nucleotide frequency matrix of the sites
//!   (rows in reverse order);
//! - `RBSL NPOS` (NPOS = L - W + 1 possible positions), `RBSD` the
//!   distribution of the motif's distance to the start codon, smoothed
//!   with a 0.2 / 0.6 / 0.2 kernel.
//!
//! Three Gibbs output formats are recognised (version 3 = Gibbs3, 1 = NCBI
//! gibbs, 10 = a simple `site position` list).  The parsers keep the
//! original's tokenising rules exactly, including which odd lines produce
//! warnings on stdout.  No randomness is involved in probuild itself.

use super::io::{self, OFile};
use super::params::{atoi, Params};
use super::pm::{fixed, g, P};

/// Smoothing kernel of the position distribution (0x454c5a; doubles
/// 0x3fc999999999999a, 0x3fe3333333333333).
const KERNEL: [f64; 3] = [0.2, 0.6, 0.2];

/// stdout write followed by a flush (`cout << ... << endl`).
fn out(s: &[u8]) {
    io::out(s);
    io::flush_stdout();
}

/// The lines of a file as `while (!eof) getline()` delivers them: the data
/// split at '\n', so a trailing "\n" gives a final empty line and '\r' is
/// kept.  `None` if the file cannot be read.
fn getlines(path: &[u8]) -> Option<Vec<Vec<u8>>> {
    let d = std::fs::read(io::path(path)).ok()?;
    Some(d.split(|&c| c == b'\n').map(|s| s.to_vec()).collect())
}

/// `std::string::find_first_not_of(set, pos)`; `None` = npos (also for a
/// `None` start position, which lets the tokeniser chain the calls).
fn find_first_not_of(s: &[u8], set: &[u8], pos: Option<usize>) -> Option<usize> {
    let pos = pos?;
    (pos..s.len()).find(|&i| !set.contains(&s[i]))
}

/// `std::string::find_first_of(set, pos)`; `None` = npos.
fn find_first_of(s: &[u8], set: &[u8], pos: Option<usize>) -> Option<usize> {
    let pos = pos?;
    (pos..s.len()).find(|&i| set.contains(&s[i]))
}

/// Tokenising loop of the parsers 0x4531a6 / 0x453bb8 / 0x452e00: at most
/// `limit` tokens separated by runs of `delims`.  Returns the tokens and a
/// count `n`.  With `stop_at_f`, the loop ends right after pushing a token
/// equal to `F` without counting it, so `n` is the index of that `F` (or
/// the token count if there is none); `toks` may then hold one more token
/// than `n`.
fn tokens(line: &[u8], delims: &[u8], limit: usize, stop_at_f: bool) -> (Vec<Vec<u8>>, usize) {
    let mut toks = Vec::new();
    let mut a = find_first_not_of(line, delims, Some(0));
    let mut b = find_first_of(line, delims, a);
    let mut n = 0;
    while n < limit {
        let Some(av) = a else { break };
        toks.push(match b {
            None => line[av..].to_vec(),
            Some(bv) => line[av..bv].to_vec(),
        });
        a = find_first_not_of(line, delims, b);
        b = find_first_of(line, delims, a);
        if stop_at_f && toks[n] == b"F" {
            break;
        }
        n += 1;
    }
    (toks, n)
}

/// All chars are digits (empty passes; 0x457592).
fn isint(t: &[u8]) -> bool {
    t.iter().all(|c| c.is_ascii_digit())
}

/// All chars are nucleotide letters "ATCGNUatcgnu" (empty passes;
/// 0x4575d2).
fn isnuc(t: &[u8]) -> bool {
    t.iter().all(|c| b"ATCGNUatcgnu".contains(c))
}

/// Site rule shared by the version 3 and version 1 parsers, for a line of
/// `k` = 6 or 7 fields.  Requires t0, t1, t2 int and t3, t4 nucleotides.
///
/// - k == 7, t5 nuc, t6 int: the Gibbs3 line with both flanks
///   `n, s  pos  left SITE right  end ...`: site = t4, position = t2.
/// - k == 6, t5 int: one flank is missing; at position 1 there is no left
///   flank, so the site is t3, otherwise t4.
/// - anything else: the line is dropped silently.
fn take_site(k: usize, t: &[Vec<u8>], sites: &mut Vec<Vec<u8>>, pos: &mut Vec<i32>) {
    if !(isint(&t[0]) && isint(&t[1]) && isint(&t[2]) && isnuc(&t[3]) && isnuc(&t[4])) {
        return;
    }
    if k == 7 && isnuc(&t[5]) && isint(&t[6]) {
        sites.push(t[4].clone());
        pos.push(atoi(&t[2]));
        return;
    }
    if k == 6 && isint(&t[5]) {
        if atoi(&t[2]) == 1 {
            sites.push(t[3].clone());
        } else {
            sites.push(t[4].clone());
        }
        pos.push(atoi(&t[2]));
    }
}

/// `prefix + line + "\n"` as bytes.
fn with_line(prefix: &str, line: &[u8]) -> Vec<u8> {
    let mut m = prefix.as_bytes().to_vec();
    m.extend_from_slice(line);
    m.push(b'\n');
    m
}

/// Format detection (0x452c0a): reads up to 11 lines, stopping at the
/// first that starts with `NAME` or `Gibbs`, and classifies the LAST line
/// read: `Gibbs` -> 3, `NAME` -> 10, else 1.  0 = file cannot be opened
/// (`error: opening file <f>` printed).
fn detect(path: &[u8]) -> i32 {
    let Some(lines) = getlines(path) else {
        io::err(&with_line("error: opening file ", path));
        return 0;
    };
    let mut cur: &[u8] = b"";
    for (cnt, ln) in lines.iter().enumerate() {
        cur = ln;
        if ln.starts_with(b"NAME") || ln.starts_with(b"Gibbs") || cnt + 1 > 10 {
            break;
        }
    }
    if cur.starts_with(b"Gibbs") {
        3
    } else if cur.starts_with(b"NAME") {
        10
    } else {
        1
    }
}

/// Gibbs3 format (0x4531a6).  Sites start after the line beginning with
/// `Num Motifs`; blank lines (" \t" only) are skipped, a line starting with
/// `*` ends the list.  Each line is split on ", \t" into at most 10
/// tokens up to the `F` column; k = (index of `F`) - 1 fields.  k other
/// than 6 or 7 prints `warning 1: unexpected line in gibbs paring: ...`
/// (sic); a site with no flank on either side gives k = 5 and this
/// warning.
fn parse_v3(path: &[u8], sites: &mut Vec<Vec<u8>>, pos: &mut Vec<i32>) {
    let mut started = false;
    for line in getlines(path).unwrap_or_default() {
        if !started {
            if line.starts_with(b"Num Motifs") {
                started = true;
            }
            continue;
        }
        let Some(a) = find_first_not_of(&line, b" \t", Some(0)) else { continue };
        if line[a] == b'*' {
            break;
        }
        let (t, n) = tokens(&line, b", \t", 10, true);
        let k = n as i64 - 1;
        if k == 6 || k == 7 {
            take_site(k as usize, &t, sites, pos);
        } else {
            out(&with_line("warning 1: unexpected line in gibbs paring: ", &line));
        }
    }
}

/// `pat` occurs in `s` at or after `from` (`std::string::find != npos`).
fn find(s: &[u8], pat: &[u8], from: usize) -> bool {
    from <= s.len() && s[from..].windows(pat.len()).any(|w| w == pat)
}

/// NCBI gibbs format (0x453bb8), verified on hand-made files only.  Every
/// line is examined: blank lines and `MOTIF A` lines are skipped, a line
/// with `sites:`, `(` or `)` ends the list.  Split on "- \t" into at most
/// 8 tokens; 6 or 7 tokens are a site line (rule of [`take_site`]),
/// followed by `warning: unexpected format in gibbs line ...` if its first
/// number is not the count of sites stored so far; any other count prints
/// `warning: unexpected line in gibbs paring: ...`.
fn parse_v1(path: &[u8], sites: &mut Vec<Vec<u8>>, pos: &mut Vec<i32>) {
    for line in getlines(path).unwrap_or_default() {
        let Some(a) = find_first_not_of(&line, b" \t", Some(0)) else { continue };
        if find(&line, b"MOTIF A", a) {
            continue;
        }
        if find(&line, b"sites:", 0) || find_first_of(&line, b"()", Some(0)).is_some() {
            break;
        }
        let (t, k) = tokens(&line, b"- \t", 8, false);
        if k == 6 || k == 7 {
            take_site(k, &t, sites, pos);
            if atoi(&t[0]) as i64 != sites.len() as i64 {
                out(&with_line("warning: unexpected format in gibbs line ", &line));
            }
        } else {
            out(&with_line("warning: unexpected line in gibbs paring: ", &line));
        }
    }
}

/// Simple list format (0x452e00): after the line starting with `VERSION`,
/// every non-blank line is `site position` (split on " \t", at most 2
/// tokens, further ones ignored).  Sites are not validated.  A line with
/// one token prints `warning 10: unexpected line in gibbs paring: ...`.
fn parse_v10(path: &[u8], sites: &mut Vec<Vec<u8>>, pos: &mut Vec<i32>) {
    let mut started = false;
    for line in getlines(path).unwrap_or_default() {
        if !started {
            if line.starts_with(b"VERSION") {
                started = true;
            }
            continue;
        }
        if find_first_not_of(&line, b" \t", Some(0)).is_none() {
            continue;
        }
        let (t, n) = tokens(&line, b" \t", 2, false);
        if n == 2 {
            sites.push(t[0].clone());
            pos.push(atoi(&t[1]));
        } else {
            out(&with_line("warning 10: unexpected line in gibbs paring: ", &line));
        }
    }
}

/// Record count and letters per record of the `--seq` FASTA (0x43c50c).
/// A record starts after optional white space with an optional `>` line
/// and runs to the next `>`; only ASCII letters are counted.  `None` if
/// the core loader fails (it prints the message).
fn record_lengths(path: &[u8]) -> Option<Vec<usize>> {
    let d = io::load_file(path)?;
    let n = d.len();
    let is_space = |c: u8| matches!(c, b' ' | b'\t' | b'\n' | 0x0b | 0x0c | b'\r');
    let mut i = 0;
    let mut lens = Vec::new();
    while i < n {
        let mut a = i;
        while a < n && is_space(d[a]) {
            a += 1;
        }
        let mut b = a;
        if a < n && d[a] == b'>' {
            while b < n && d[b] != b'\n' {
                b += 1;
            }
        }
        let mut e = b;
        while e < n && d[e] != b'>' {
            e += 1;
        }
        lens.push(d[b..e].iter().filter(|c| c.is_ascii_alphabetic()).count());
        i = e;
    }
    Some(lens)
}

/// Writes report entries `label + value + "\n"` in the given order,
/// appending to `fname` (0x4214e0 with append), or to stdout if the name
/// is empty.  An open failure is silent.
fn write_report(fname: &[u8], items: &[(&[u8], Vec<u8>)]) {
    let mut s = Vec::new();
    for (l, v) in items {
        s.extend_from_slice(l);
        s.extend_from_slice(v);
        s.push(b'\n');
    }
    if fname.is_empty() {
        out(&s);
    } else {
        let mut f = OFile::open(fname, true);
        f.w(&s);
        f.close();
    }
}

/// `--gibbs` handler (0x455cc4).  See the module docs and models.md §14
/// for the options, checks and output layout.
pub fn run(pp: &Params) -> bool {
    let mut p = P::new(pp);
    let p = &mut p;
    // decimals and pseudocounts of the matrix (_mat) and of the position
    // distribution (_dur)
    let (mut prec_mat, mut prec_dur, mut pc_mat, mut pc_dur) = (3i32, 3i32, 1i32, 1i32);
    let mut model_file = Vec::new();
    let mut log_file = Vec::new();
    let mut chi_file = Vec::new();
    if p.has("--mod") {
        match p.str_req("--mod") {
            Some(m) => model_file = m,
            None => return false,
        }
    }
    if p.has("--log") {
        match p.str_opt("--log", Vec::new()) {
            Some(l) => log_file = l,
            None => return false,
        }
    }
    if p.has("--chi") {
        // result ignored: an invalid value prints and continues without chi
        if let Some(c) = p.str_req("--chi") {
            chi_file = c;
        }
    }
    let Some(gibbs_file) = p.str_req("--gibbs") else { return false };
    let Some(seq_file) = p.str_req("--seq") else { return false };
    for (k, dst) in [
        ("--GIBBS_PRECISION_MAT", &mut prec_mat),
        ("--GIBBS_PRECISION_DUR", &mut prec_dur),
        ("--GIBBS_PSEUDOCOUNTS_MAT", &mut pc_mat),
        ("--GIBBS_PSEUDOCOUNTS_DUR", &mut pc_dur),
    ] {
        let (ok, v) = p.int_opt(k, *dst, false);
        if !ok {
            return false;
        }
        *dst = v;
    }
    // sites (motif instances) and their 1-based start positions within the
    // pre-start sequences, as reported by the sampler
    let mut sites: Vec<Vec<u8>> = Vec::new();
    let mut pos: Vec<i32> = Vec::new();
    match detect(&gibbs_file) {
        3 => parse_v3(&gibbs_file, &mut sites, &mut pos),
        10 => parse_v10(&gibbs_file, &mut sites, &mut pos),
        1 => parse_v1(&gibbs_file, &mut sites, &mut pos),
        _ => return false,
    }
    // 1. no sites: note in the log (or on stdout), success, no model
    if sites.is_empty() {
        write_report(&log_file, &[(b"Gibbs output is empty: ", gibbs_file.clone())]);
        return true;
    }
    // 2. all sites must have the width of the first (0x45462c)
    let w = sites[0].len();
    if sites.iter().any(|s| s.len() != w) {
        io::err(&with_line("parse gibbs, error in file format ", &gibbs_file));
        return false;
    }
    // 3./4. one equally long sequence per site (0x43c50c, 0x454716).  The
    // swapped labels of the count message are the original's.
    let Some(lens) = record_lengths(&seq_file) else { return false };
    if lens.len() != sites.len() || lens.len() != pos.len() {
        io::err(format!("gibbs lines: {}\nseq lines: {}\nerror: in gibbs parse input data\n", lens.len(), sites.len()).as_bytes());
        return false;
    }
    if lens.is_empty() {
        return false;
    }
    let l = lens[0];
    if lens.iter().any(|&x| x != l) {
        return false;
    }
    // Number of possible motif positions.  Each position is converted to
    // the distance index NPOS - pos: 0 = the motif ends at the end of the
    // pre-start sequence, next to the start codon.
    let npos = l as i64 - w as i64 + 1;
    for x in pos.iter_mut() {
        if npos < *x as i64 {
            io::err(b"error: contradiction in sequence size and motif position in gibbs\n");
            return false;
        }
        if *x < 1 {
            io::err(b"error: unexpected motif position in gibbs\n");
            return false;
        }
        *x = (npos - *x as i64) as i32;
    }
    let npos = npos as usize;
    // 5. position distribution with pseudocount (0x454a16); the total is
    // accumulated as `(pc + x) + tot` in index order
    let mut dur = vec![0.0f64; npos];
    for &x in &pos {
        dur[x as usize] += 1.0;
    }
    let mut tot = 0.0;
    for &x in &dur {
        tot = pc_dur as f64 + x + tot;
    }
    if tot == 0.0 {
        return false;
    }
    let dur: Vec<f64> = dur.iter().map(|&x| (pc_dur as f64 + x) / tot).collect();
    // 6. optional chi-square statistic of the distribution against its
    // mean (0x455b66), sums in the original's `x + acc` order
    if !chi_file.is_empty() {
        let mut mean = 0.0;
        for &x in &dur {
            mean = x + mean;
        }
        mean /= dur.len() as f64;
        let mut chi = 0.0;
        for &x in &dur {
            chi = ((x - mean) * (x - mean)) / mean + chi;
        }
        let mut f = OFile::open(&chi_file, false);
        if !f.is_open() {
            io::err(&with_line("error: opening file ", &chi_file));
            return false;
        }
        f.w(format!("{}\n", g(chi)).as_bytes());
        f.close();
    }
    // 7. smoothing (0x454c5a): truncated at the edges, not renormalised;
    // `new = new + dur * k` for j = -1, 0, 1
    let n = dur.len() as i64;
    let mut smoothed = vec![0.0f64; npos];
    for i in 0..n {
        for j in -1i64..=1 {
            if i + j > -1 && i + j < n {
                smoothed[i as usize] = smoothed[i as usize] + dur[(i + j) as usize] * KERNEL[(j + 1) as usize];
            }
        }
    }
    let dur = smoothed;
    // 8. site matrix (0x454ede): only UPPERCASE A C G T count, everything
    // else (lowercase, N, U) goes to column 4 and is ignored by the total
    // `4*pc + G + C + A + T` (int shift first, then this order)
    let mut cnt = vec![[0.0f64; 5]; w];
    for s in &sites {
        for j in 0..w {
            let c = match s[j] {
                b'A' => 0,
                b'C' => 1,
                b'G' => 2,
                b'T' => 3,
                _ => 4,
            };
            cnt[j][c] += 1.0;
        }
    }
    let mut mat = Vec::new();
    for c in &cnt {
        let t = (pc_mat.wrapping_shl(2)) as f64 + c[2] + c[1] + c[0] + c[3];
        if t == 0.0 {
            return false;
        }
        mat.push([
            (pc_mat as f64 + c[0]) / t,
            (pc_mat as f64 + c[1]) / t,
            (pc_mat as f64 + c[2]) / t,
            (pc_mat as f64 + c[3]) / t,
        ]);
    }
    let row = |r: &[f64; 4]| r.iter().map(|&v| fixed(prec_mat, v, true)).collect::<Vec<_>>().join(" ");
    // model block, appended: matrix rows last position first (0x455712
    // mode 1), distribution in index order (0x45548a mode 1)
    if !model_file.is_empty() {
        let rbsm: String = mat.iter().rev().map(|r| row(r) + "\n").collect();
        let rbsd: String = dur.iter().map(|&v| fixed(prec_dur, v, true) + "\n").collect();
        write_report(
            &model_file,
            &[
                (b"RBSS ", w.to_string().into_bytes()),
                (b"RBSM \n", rbsm.into_bytes()),
                (b"RBSL ", npos.to_string().into_bytes()),
                (b"RBSD \n", rbsd.into_bytes()),
            ],
        );
    }
    // log block (only if the key is present; an empty name means stdout):
    // distribution as `-index\tvalue` from the farthest position down to 0,
    // then the matrix forward with a header
    if p.has("--log") {
        let mut s = String::new();
        for k in (0..npos).rev() {
            s += &format!("{}\t{}\n", -(k as i64), fixed(prec_dur, dur[k], true));
        }
        s += "Pos   A   C   G   T\n";
        for (j, r) in mat.iter().enumerate() {
            s += &format!("{} {}\n", j + 1, row(r));
        }
        write_report(&log_file, &[(b"", s.into_bytes())]);
    }
    true
}
