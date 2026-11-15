//! gmhmmp-rs: re-implementation of GeneMark.hmm PROKARYOTIC 3.20 (`gmhmmp`,
//! the gene finder of GeneMarkS-T inside GeneMark-ETP), reconstructed from
//! the binary (BuildID 5107402e52b7b0a67d36869b4c5daf557f678ff1). Output
//! files, stdout, stderr and exit codes are byte-identical to the original;
//! every floating-point value is bit-identical. Specifications:
//! `re-notes/gmhmmp/` (input.md, model.md, viterbi.md, output.md).
//!
//! # Pipeline (main 0x403380, implemented in [`run`])
//!
//! 1. Command line (`cli.rs`): settings, usage, `-V` dumps, validation.
//! 2. Optional log (`-L` file or stdout with `-v`), model file (`model.rs`:
//!    one model, or 41 GC-specific models for a meta file), optional
//!    evidence GFF (`evidence.rs`, `-e`).
//! 3. Output file header (`output.rs`).
//! 4. Per FASTA record (`seq.rs`): encode, skip records with < 20 letters,
//!    gap split, `-c` defline parameters, choice of the GC model (meta),
//!    skip records shorter than MINC, then `predict.rs`: RBS tables,
//!    Markov image, R/Y stop resolution, codon flags, evidence edits, site
//!    map, start/RBS scores, length tables, DP, traceback, N ends, gene
//!    numbers, logodd; finally the record's output block (`output.rs`).
//!
//! Several model fields are per-record state, as in the binary: `-r`
//! absent sets RBSS to -1, the RBS spacer distribution is renormalised in
//! place, and the length tables grow when a record needs longer ones
//! (model.md 6.3/6.4).

pub mod cli;
pub mod evidence;
pub mod model;
pub mod num;
pub mod output;
pub mod predict;
pub mod seq;

use std::io::Write;

/// Debug hook (`gmhmmp-rs --dump-model MODEL DIR`, not in the original):
/// dump the loaded model(s) as raw little-endian f64 arrays in the layout
/// of the gdb dumper (re-notes/gmhmmp/model.md section 9), for bit-exact
/// comparison of the model transformations.
fn dump_models(path: &str, dir: &str) {
    let ms = match model::load(path) {
        Ok(m) => m,
        Err(e) => {
            print!("{}", e.stdout);
            eprint!("{}", e.stderr);
            std::process::exit(1);
        }
    };
    let list: Vec<(String, &model::Model)> = match &ms {
        model::Models::Single(m) => vec![("single".into(), m.as_ref())],
        model::Models::Meta(v) => v.iter().enumerate().map(|(i, m)| (format!("gc{}", i + 30), m)).collect(),
    };
    std::fs::create_dir_all(dir).unwrap();
    let w = |name: &str, v: &[f64]| {
        let b: Vec<u8> = v.iter().flat_map(|x| x.to_le_bytes()).collect();
        std::fs::write(format!("{dir}/{name}.bin"), b).unwrap();
    };
    for (tag, m) in list {
        let k = &m.markov;
        for j in 0..8 {
            w(&format!("{tag}_c1_{j}"), &k.c1[j]);
            w(&format!("{tag}_c2_{j}"), &k.c2[j]);
        }
        w(&format!("{tag}_nc"), &k.nc);
        for (j, t) in k.nl.iter().enumerate() {
            w(&format!("{tag}_nl_{j}"), t);
        }
        w(&format!("{tag}_cod1"), &m.cod1);
        w(&format!("{tag}_cod2"), &m.cod2);
        w(&format!("{tag}_nonc"), &m.nonc);
        if m.rbss > 0 {
            w(&format!("{tag}_rbsm"), &m.rbsm);
            w(&format!("{tag}_rbsd"), &m.rbsd);
        }
        w(&format!("{tag}_durC"), &m.dur_c[..m.maxc as usize]);
        w(&format!("{tag}_durP"), &m.dur_p[..m.maxc as usize]);
        w(&format!("{tag}_durN"), &m.dur_n[..m.maxn as usize + 8]);
        std::fs::write(format!("{dir}/{tag}_flags.bin"), m.flags).unwrap();
        let mut sc = vec![m.cd1p, m.cd2p, m.ndec, m.cdec, m.cdcd];
        sc.extend_from_slice(&m.log_start);
        sc.extend_from_slice(&k.nonc0);
        for p in &k.cod1mono {
            sc.extend_from_slice(p);
        }
        w(&format!("{tag}_scalars"), &sc);
        std::fs::write(format!("{dir}/{tag}_ints.txt"), format!("{} {} {} {} {} {} {:?}\n", m.rbss, m.rbsl, m.ordm, m.minc, m.maxn, m.maxc, m.name)).unwrap();
    }
}

/// `ctime(time(NULL))` of libc (local time zone, ends in '\n'), as in the
/// header's Date line; the only output line that differs between runs.
fn ctime_now() -> String {
    unsafe {
        let mut t: libc::time_t = 0;
        libc::time(&mut t);
        let mut buf = [0 as libc::c_char; 64];
        let p = libc::ctime_r(&t, buf.as_mut_ptr());
        if p.is_null() {
            return "\n".into();
        }
        std::ffi::CStr::from_ptr(p).to_string_lossy().into_owned()
    }
}

/// Print `out` to stdout and `err` to stderr, flush, exit with code 1.
fn die(out: &str, err: &str) -> ! {
    print!("{out}");
    eprint!("{err}");
    let _ = std::io::stdout().flush();
    std::process::exit(1);
}

/// Output file for -o / -A / -D. The binary opens the file once with trunc
/// (so it exists even if empty), then appends each pending text block
/// followed by `std::endl` (FLUSH, output.md 3.1): every block ends with
/// one extra '\n'.
struct OutFile {
    f: std::fs::File,
}

impl OutFile {
    fn create(name: &str) -> OutFile {
        match std::fs::File::create(name) {
            Ok(f) => OutFile { f },
            Err(_) => die("", &format!("Error on open output file {name}\n")),
        }
    }
    /// FLUSH: write `text`, then the endl newline.
    fn flush_block(&mut self, text: &str) {
        let _ = self.f.write_all(text.as_bytes());
        let _ = self.f.write_all(b"\n");
    }
}

/// Accumulated stage times in seconds (developer switch GMHMMP_RS_TIMING):
/// image+map, sites, lengths, dp, traceback, output.
static mut TIMES: [f64; 6] = [0.0; 6];
/// Add the time since `*t` to stage `k` and restart the clock.
fn tstamp(t: &mut std::time::Instant, k: usize) {
    let now = std::time::Instant::now();
    unsafe {
        TIMES[k] += (now - *t).as_secs_f64();
    }
    *t = now;
}

/// The whole program (main 0x403380); returns the exit code. Fatal errors
/// exit directly through [`die`].
fn run(argv: &[String]) -> i32 {
    use model::Models;
    let timing = std::env::var_os("GMHMMP_RS_TIMING").is_some();
    let (s, pre) = match cli::settings(argv) {
        Ok(x) => x,
        Err(m) => die(&m, ""),
    };
    // -V settings dumps
    print!("{pre}");
    // logging (main 0x403380)
    let mut log: Option<std::fs::File> = None;
    if !s.log_file.is_empty() {
        match std::fs::File::create(&s.log_file) {
            Ok(mut f) => {
                let _ = f.write_all(b"Start of log session\n");
                log = Some(f);
            }
            Err(_) => die(&format!("error on open logfile: {}\n", s.log_file), ""),
        }
    }
    if s.debug {
        // "loging" (sic) as in the binary
        print!("debug mode\nCLASS LOGINFO\n");
        if log.is_some() {
            print!("loging into file: {}\n", s.log_file);
        } else {
            print!("loging into stdout\n");
        }
    }
    let mut models = match model::load(&s.model) {
        Ok(m) => m,
        Err(e) => {
            if e.direct_exit {
                die(&e.stdout, &e.stderr);
            }
            die(&format!("{}error on load parameters from file: {}\n", e.stdout, s.model), &e.stderr);
        }
    };
    if s.three_coding {
        // the original reads through NULL table pointers and terminates with SIGSEGV
        die("", "gmhmmp-rs: error: option -u (three coding models) is not supported (not handled by the original); run without -u\n");
    }
    let data = match std::fs::read(&s.seq) {
        Ok(d) => d,
        Err(_) => {
            eprint!("Can't open input file {}\n", s.seq);
            return 1;
        }
    };
    // log stream: -L file, else stdout with -v, else discarded
    let mut log_lines: Vec<String> = Vec::new();
    let flush_log = |lines: &mut Vec<String>, log: &mut Option<std::fs::File>| {
        for l in lines.drain(..) {
            if let Some(f) = log.as_mut() {
                let _ = writeln!(f, "{l}");
            } else if s.verbose {
                println!("{l}");
            }
        }
    };
    let mut evidence = if s.evidence.is_empty() {
        None
    } else {
        match evidence::load(&s.evidence, &mut log_lines) {
            Ok(e) => {
                flush_log(&mut log_lines, &mut log);
                Some(e)
            }
            Err(m) => die("", &m),
        }
    };
    let mut out = OutFile::create(&s.out);
    let mut aa_f = if s.print_aa { Some(OutFile::create(&s.aa_file)) } else { None };
    let mut nt_f = if s.print_nt { Some(OutFile::create(&s.nt_file)) } else { None };
    let date = ctime_now();
    let single_name = match &models {
        Models::Single(m) => Some(m.name.clone().unwrap_or_else(|| "(null)".into())),
        Models::Meta(_) => None,
    };
    out.flush_block(&output::header(&s, &date, single_name.as_deref()));
    let mut reader = seq::Reader::new(data);
    // gene number over the whole file (0x412640, a counter in main)
    let mut gene_counter = 0i32;
    // -b numbering of kept genes (a global counter), never reset
    let mut best_counter = 0i32;
    // the RBS model warning of the port is printed once per run
    let mut warned_rbs_model = false;
    loop {
        let raw = reader.next_raw();
        let mut rec = seq::encode(&raw);
        // -e: out-of-sequence check and M masking inside the reader, before
        // the copy and the gap split (input.md 5.3)
        if let (Some(r), Some(ev)) = (rec.as_mut(), evidence.as_mut()) {
            if let Some(list) = ev.by_seq.get_mut(&r.id) {
                evidence::reader_step(list, &mut r.codes, &mut log_lines);
                flush_log(&mut log_lines, &mut log);
            }
        }
        let letters = rec.as_ref().map(|r| r.codes.len()).unwrap_or(0);
        // eof only after a read with no letters and no pending defline, so
        // the last real record is always processed
        let eof = letters == 0 && !reader.has_pending();
        if letters < 20 {
            if eof {
                break;
            }
            continue;
        }
        let mut rec = rec.unwrap();
        // copy before the gap split (+0x1a0): N ends, -A/-D, gets 5/6 at R/Y
        rec.orig = rec.codes.clone();
        if s.split_gaps != 0 {
            seq::split_gaps(&mut rec.codes);
        }
        // -c defline parameters (0x407f70, 0x408120)
        let mut p = s.ini_noncoding;
        if s.parse_defline {
            if let Some(i) = rec.defline.find("[shape=") {
                if let Some(j) = rec.defline[i + 7..].find(']') {
                    match &rec.defline[i + 7..i + 7 + j] {
                        "circular" | "partial" => p = 0.5,
                        "linear" => p = 1.0,
                        _ => {}
                    }
                }
            }
            if let Some(i) = rec.defline.find("[gc=") {
                if let Some(j) = rec.defline[i + 4..].find(']') {
                    let v = &rec.defline[i + 4..i + 4 + j];
                    if !v.is_empty() {
                        let d = num::strtod(v) as i32;
                        // `(unsigned)(int)d > 100`: negative values fail too
                        if d as u32 > 100 {
                            die("error in format of GC parameter inside defline\n", "");
                        }
                        rec.gc = d;
                    }
                }
            }
        }
        // meta: GC model choice and its "Model information" block, BEFORE
        // the MINC length test
        let m: &mut model::Model = match &mut models {
            Models::Single(m) => m.as_mut(),
            Models::Meta(v) => {
                let i = Models::index_for_gc(rec.gc);
                let name = v[i].name.clone().unwrap_or_else(|| "(null)".into());
                out.flush_block(&output::model_info(&s, &name));
                &mut v[i]
            }
        };
        let n = rec.codes.len();
        if (n as i64) < m.minc as i64 {
            continue;
        }
        if !s.rbs {
            m.rbss = -1;
        }
        let mut stage_clock = std::time::Instant::now();
        // codes after the gap split, before R/Y resolution: what the Markov
        // image and the RBS scoring see (image byte 0)
        let pre_ry = rec.codes.clone();
        let rbs_tables = if m.rbss > 0 && m.rbsl >= 0 && m.rbss <= predict::RBS_MAX_SPACER_W { Some(predict::rbs_prepare(m, &pre_ry)) } else { None };
        // 0x40c170: `16 < (unsigned)(2 * ordm)`
        if m.ordm > 8 {
            die(": fatal error, [order] > [short int]\n", "");
        }
        let markov_idx = predict::markov_image(&pre_ry, m.ordm);
        seq::resolve_ry(&mut rec, m.stop[0], m.stop[1], m.stop[2]);
        let (mut codon, _) = seq::codon_image(&rec.codes, &m.flags);
        let mut ev_maps = None;
        if let Some(ev) = evidence.as_ref() {
            let mut out_txt = String::new();
            let empty = Vec::new();
            let list = ev.by_seq.get(&rec.id).unwrap_or(&empty);
            ev_maps = Some(evidence::image_step(list, &mut codon, &rec.id, &mut out_txt));
            print!("{out_txt}");
        }
        let mut map = predict::build_map(&codon, p, m.minc as i64, s.strand, s.evidence_score, ev_maps.as_ref());
        tstamp(&mut stage_clock, 0);
        let issues = predict::site_scores(&mut map.sites, &pre_ry, m, rbs_tables.as_ref());
        // fallbacks of the port where the binary reads memory it does not
        // own (predict::site_scores); model-level causes are reported once
        if issues.rows_missing && !warned_rbs_model {
            warned_rbs_model = true;
            eprintln!(
                "gmhmmp-rs: warning: -r requested, but the RBS context model in {} has RBSS {} < -RBSL {} (motif rows missing); predicting without RBS scores (not handled by the original)",
                s.model, m.rbss, -m.rbsl
            );
        }
        if issues.too_wide && !warned_rbs_model {
            warned_rbs_model = true;
            eprintln!(
                "gmhmmp-rs: warning: -r requested, but the RBS motif width RBSS {} in {} is above {}; predicting without RBS scores (not handled by the original)",
                m.rbss, s.model, predict::RBS_MAX_SPACER_W
            );
        }
        if issues.window_outside {
            eprintln!(
                "gmhmmp-rs: warning: record {}: the RBS context window of model {} (RBSS {}, RBSL {}) leaves the sequence for start codons next to its ends; those starts get no RBS score (not handled by the original)",
                rec.id, s.model, m.rbss, m.rbsl
            );
        }
        if !s.motif_file.is_empty() && m.rbss > 0 {
            // rewritten per record: the file keeps only the last record
            if std::fs::write(&s.motif_file, predict::k_file(&map.sites)).is_err() {
                die(&format!("Error on open file {}\n", s.motif_file), "");
            }
        }
        tstamp(&mut stage_clock, 1);
        // length tables (stateful per model, model.md 6.3/6.4)
        m.coding_lengths(map.max_gene_len as i32);
        m.noncoding_lengths(map.max_nc_span as i32 + 8);
        // 0x40dcc0: GC of ORFs; only its error exit is observable
        if predict::gc_error(&map.sites) {
            die("", "error in code GCofInterval\n");
        }
        if map.sites.iter().any(|x| x.link_stale) {
            // Intended difference: with MINC < 7 the binary's DP follows
            // pointers into the freed temporary site array (predict::compact)
            die(
                "",
                &format!(
                    "gmhmmp-rs: error: model file {}: MINC {} allows genes shorter than 7 bp (start and stop codons overlapping) in record {} (not handled by the original); use a model with MINC >= 7\n",
                    s.model, m.minc, rec.id
                ),
            );
        }
        let dp_params = predict::DpParams { hmin: s.min_heu_len as i64, ev_score: s.evidence_score, allow_overlap: s.overlap, q: s.overlap_penalty, p };
        tstamp(&mut stage_clock, 2);
        predict::run_dp(&mut map.sites, &markov_idx, &pre_ry, m, &dp_params);
        tstamp(&mut stage_clock, 3);
        let mut genes = predict::traceback(&mut map.sites, n as i64, m);
        predict::n_ends(&mut genes, &rec.orig);
        for g in genes.iter_mut() {
            gene_counter += 1;
            g.id = gene_counter;
        }
        let mut short_lo = false;
        for g in genes.iter_mut() {
            let (lo, short) = predict::logodd(g, &markov_idx, m);
            g.logodd = lo;
            short_lo |= short;
        }
        if markov_idx.oob.get() {
            eprintln!(
                "gmhmmp-rs: warning: record {}: genes shorter than ORDM + 9 bp (model {} has MINC {}, ORDM {}) need Markov windows outside the sequence; scoring those windows as unknown bases (N) (not handled by the original)",
                rec.id, s.model, m.minc, m.ordm
            );
        }
        if short_lo {
            eprintln!(
                "gmhmmp-rs: warning: record {}: a predicted gene is shorter than ORDM + 8 bp (after trimming N ends), too short for the logodd score; using only its last Markov window (not handled by the original)",
                rec.id
            );
        }
        if s.rbs && s.motif_info && s.motif_pos {
            output::k_positions(&mut genes, m.rbss);
        }
        tstamp(&mut stage_clock, 4);
        let block = output::record_block(&s, &rec.defline, &mut genes, &rec.orig, m.stop[2], &mut best_counter, ev_maps.as_mut());
        print!("{}", block.stdout);
        out.flush_block(&block.main);
        if let Some(f) = aa_f.as_mut() {
            f.flush_block(&block.aa);
        }
        if let Some(f) = nt_f.as_mut() {
            f.flush_block(&block.nt);
        }
        tstamp(&mut stage_clock, 5);
        if eof {
            break;
        }
    }
    if timing {
        let t = unsafe { TIMES };
        eprintln!("[timing] image+map {:.3} sites {:.3} lengths {:.3} dp {:.3} traceback {:.3} output {:.3}", t[0], t[1], t[2], t[3], t[4], t[5]);
    }
    0
}

/// Entry point of `gmhmmp-rs`.
pub fn main() {
    let args: Vec<String> = std::env::args().collect();
    if args.len() == 4 && args[1] == "--dump-model" {
        dump_models(&args[2], &args[3]);
        return;
    }
    let rc = run(&args);
    let _ = std::io::stdout().flush();
    std::process::exit(rc);
}
