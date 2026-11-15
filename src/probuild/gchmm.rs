//! `--gc_hmm` (0x457e34): two-state GC segmentation by Viterbi decoding.
//! models.md §15.  Not used by the GeneMark pipelines.
//!
//! Each FASTA record is reduced to a two-letter alphabet (A/T -> 'A',
//! C/G -> 'C', other letters 'N') and decoded with a two-state HMM: a
//! high-GC state 'f' (emits C with probability `--gchigh`, default 0.52)
//! and a low-GC state 's' (`--gclow`, 0.26), with self-transition
//! probabilities `--trhigh` / `--trlow` (0.999 each).  N emits with
//! log-probability 0 in both states.  The low-GC ('s') segments of the
//! most probable path are appended to the `--gc_hmm` file as
//! `name\tgc_hmm\tlow\tstart\tend\t.\t.\t.\t.` (1-based, inclusive).
//!
//! All arithmetic is in natural-log space with glibc `log` and the exact
//! operand grouping of the original; ties go to 's'.  The original's debug
//! printing (object flag +0x48) is always off and not reproduced.

use super::io::{self, OFile};
use super::params::Params;
use super::pm::P;

/// ln(1/2), the uniform initial state probability (double
/// 0xbfe62e42fefa39ef).
const LN_HALF: f64 = -0.6931471805599453;

/// glibc `log` (ln 0 = -inf, ln(<0) = nan), called through libc so that
/// the results match the original's static glibc.
fn ln(x: f64) -> f64 {
    extern "C" {
        fn log(x: f64) -> f64;
    }
    unsafe { log(x) }
}

/// FASTA reader 0x44b624 with the 'A'/'C'/'N' table of 0x459ac6: a list of
/// (defline, sequence) in file order.  A record: optional white space, an
/// optional '>' line (kept with '>' and any '\r', without '\n'), then
/// everything up to the next '>'; only ASCII letters are kept and mapped
/// a/A/t/T -> 'A', c/C/g/G -> 'C', others -> 'N'.
fn read_fasta_pairs(d: &[u8]) -> Vec<(Vec<u8>, Vec<u8>)> {
    let is_space = |c: u8| matches!(c, b' ' | b'\t' | b'\n' | 0x0b | 0x0c | b'\r');
    let n = d.len();
    let mut recs = Vec::new();
    let mut pos = 0;
    while pos < n {
        let mut a = pos;
        while a < n && is_space(d[a]) {
            a += 1;
        }
        let mut b = a;
        if a < n && d[a] == b'>' {
            while b < n && d[b] != b'\n' {
                b += 1;
            }
        }
        let defl = d[a..b].to_vec();
        let mut e = b;
        while e < n && d[e] != b'>' {
            e += 1;
        }
        let seq: Vec<u8> = d[b..e]
            .iter()
            .filter(|c| c.is_ascii_alphabetic())
            .map(|&c| match c {
                b'a' | b'A' | b't' | b'T' => b'A',
                b'c' | b'C' | b'g' | b'G' => b'C',
                _ => b'N',
            })
            .collect();
        recs.push((defl, seq));
        pos = e;
    }
    recs
}

/// `--gc_hmm` handler (0x457e34).
///
/// Order of events (models.md §15): read `--gc_hmm`, create/truncate the
/// output file at once (failure: `error, on open file: <f>`, exit(1)),
/// read `--seq` and load the FASTA, only then the four double options (so
/// a bad option still leaves an empty output file), then decode every
/// record and append its segments.
pub fn run(pp: &Params) -> bool {
    let mut p = P::new(pp);
    let p = &mut p;
    let Some(out) = p.str_req("--gc_hmm") else { return false };
    {
        let f = OFile::open(&out, false);
        if !f.is_open() {
            let mut m = b"error, on open file: ".to_vec();
            m.extend_from_slice(&out);
            m.push(b'\n');
            io::err(&m);
            io::c_exit(1);
        }
    }
    let Some(seq_file) = p.str_req("--seq") else { return false };
    let Some(d) = io::load_file(&seq_file) else { return false };
    let recs = read_fasta_pairs(&d);
    // double getter 0x41fff2 for each: float characters, then strtod
    let mut vals = [0.52, 0.26, 0.999, 0.999];
    for (i, k) in ["--gchigh", "--gclow", "--trhigh", "--trlow"].iter().enumerate() {
        if p.has(k) {
            if !p.chk_float(k) {
                return false;
            }
            vals[i] = p.f64_val(k);
        }
    }
    let [gc_high, gc_low, tr_high, tr_low] = vals;
    // Log transition probabilities.  As in the original: a row-stochastic HMM
    // would use ln(1 - trhigh) for high -> low and ln(1 - trlow) for
    // low -> high; the original uses them the other way round (hl for
    // low -> high, ll for high -> low, see the scores below).  Identical
    // with the default trhigh == trlow.
    let hl = ln(1.0 - tr_high);
    let ll = ln(1.0 - tr_low);
    let hh = ln(tr_high);
    let lolo = ln(tr_low);
    // log emissions: at_* for 'A' (A/T), gc_* for 'C' (C/G)
    let at_hi = ln(1.0 - gc_high);
    let at_lo = ln(1.0 - gc_low);
    let gc_hi = ln(gc_high);
    let gc_lo = ln(gc_low);
    for (defl, s) in &recs {
        let l = s.len();
        if l == 0 {
            // not handled by the original (index 0 of empty vectors)
            io::err(&[b"probuild-rs: error: --gc_hmm: record ".as_slice(), defl, b" has no sequence letters\n"].concat());
            return false;
        }
        // Viterbi scores of the paths ending in the high / low state, and
        // the back pointers (state at i-1 on the best path into state x at
        // i)
        let mut hi = vec![0.0f64; l];
        let mut lo = vec![0.0f64; l];
        let mut back_hi = vec![b'f'; l];
        let mut back_lo = vec![b'f'; l];
        hi[0] = LN_HALF;
        lo[0] = LN_HALF;
        // operand order per letter as in the original
        match s[0] {
            b'C' => {
                hi[0] = gc_hi + hi[0];
                lo[0] = gc_lo + lo[0];
            }
            b'A' => {
                hi[0] = hi[0] + at_hi;
                lo[0] = lo[0] + at_lo;
            }
            _ => {
                hi[0] = hi[0] + 0.0;
                lo[0] = lo[0] + 0.0;
            }
        }
        for i in 1..l {
            // a: hi -> hi, b: lo -> hi, c: hi -> lo, d: lo -> lo; always
            // grouped as prev + (emission + transition) with the operand
            // order of the original
            let (a, b, c, d) = match s[i] {
                b'C' => (
                    hi[i - 1] + (gc_hi + hh),
                    lo[i - 1] + (gc_hi + hl),
                    hi[i - 1] + (gc_lo + ll),
                    lo[i - 1] + (gc_lo + lolo),
                ),
                b'A' => (
                    hi[i - 1] + (hh + at_hi),
                    lo[i - 1] + (at_hi + hl),
                    hi[i - 1] + (at_lo + ll),
                    lo[i - 1] + (lolo + at_lo),
                ),
                _ => (hi[i - 1] + (hh + 0.0), lo[i - 1] + (hl + 0.0), hi[i - 1] + (ll + 0.0), lo[i - 1] + (lolo + 0.0)),
            };
            // strict '>': a tie or NaN picks the low state 's'
            // (ucomisd / jbe in the original)
            if a > b {
                hi[i] = a;
                back_hi[i] = b'f';
            } else {
                hi[i] = b;
                back_hi[i] = b's';
            }
            if c > d {
                lo[i] = c;
                back_lo[i] = b'f';
            } else {
                lo[i] = d;
                back_lo[i] = b's';
            }
        }
        // traceback from the better final state
        let mut path = vec![b'f'; l];
        let mut k = l - 1;
        path[k] = if hi[k] > lo[k] { b'f' } else { b's' };
        while k >= 1 {
            path[k - 1] = if path[k] == b'f' { back_hi[k] } else { back_lo[k] };
            k -= 1;
        }
        // name = defline without its first char ('>')
        let name: &[u8] = if defl.is_empty() { b"" } else { &defl[1..] };
        let mut buf = Vec::new();
        // opens a low-GC line at 1-based position `pos`; the end is written
        // when the segment closes
        let open_low = |buf: &mut Vec<u8>, pos: usize| {
            buf.extend_from_slice(name);
            buf.extend_from_slice(format!("\tgc_hmm\tlow\t{pos}\t").as_bytes());
        };
        if path[0] == b's' {
            open_low(&mut buf, 1);
        }
        for i in 1..l {
            if path[i - 1] == b's' && path[i] == b'f' {
                // last low base is i-1 (0-based) = i (1-based)
                buf.extend_from_slice(format!("{i}\t.\t.\t.\t.\n").as_bytes());
            }
            if path[i - 1] == b'f' && path[i] == b's' {
                open_low(&mut buf, i + 1);
            }
        }
        if path[l - 1] == b's' {
            buf.extend_from_slice(format!("{l}\t.\t.\t.\t.\n").as_bytes());
        }
        let mut f = OFile::open(&out, true);
        f.w(&buf);
        f.close();
    }
    true
}
