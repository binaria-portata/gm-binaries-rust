//! gmhmm-rs: Rust re-implementation of GeneMark.hmm eukaryotic (gmhmme3 3.68),
//! reconstructed from the GeneMark-ETP binary for behavioural parity: same
//! stdout, same output files, same exit codes and bit-identical floating-point
//! results.
//!
//! # What the program computes
//!
//! gmhmme3 is the gene finder of GeneMark-ES/ET/EP/ETP. It predicts
//! protein-coding genes on both strands of one DNA sequence with a generalized
//! (explicit-duration, "semi-Markov") hidden Markov model:
//!
//! * hidden states: intergenic, intron (+/-, three phases), and exons
//!   (initial, internal, terminal, single; +/-);
//! * state content: Markov chains of order `$ORDER` (three-periodic for
//!   coding, homogeneous for introns and intergenic regions);
//! * state boundaries: signal sites (start codon, stop codons, donor and
//!   acceptor splice sites) scored with position-specific matrices;
//! * state durations: explicit length distributions from the model file.
//!
//! The Viterbi algorithm finds the single most probable parse of the
//! sequence; the parse is written as lst, GTF, GFF3 or "tr" (training) text.
//!
//! # Pipeline (follows `main()` at 0x404ce0 in the original binary)
//!
//! 1. [`cli`]: parse the command line (0x4075e0, 0x407c20).
//! 2. [`sequence`]: read the FASTA file into base codes 0..4 (0x4095f0).
//! 3. [`evidence`]: with `-d`, read external hints (GFF) into a per-base
//!    bitmask (0x40cc40). Done before the model is read, as in the binary.
//! 4. [`model`]: read and transform the model file (0x412fe0).
//! 5. [`sitemap`]: list every candidate signal (ATG, stops, GT/GC, AG and
//!    their reverse complements) sorted by position (0x40ed90). This list of
//!    "map records" is the time axis of the Viterbi.
//! 6. [`sequence::Sequence::build_kmer`]: per-base k-mer index for the Markov
//!    chains (0x4099a0).
//! 7. [`bp`]: only for branch-point models (`$BP 1`, GeneMark-ES `--fungus`).
//! 8. [`tracks`]: precompute per-record content and site scores (0x41aae0,
//!    0x41aa30, 0x419f60, 0x41abe0).
//! 9. [`viterbi`]: segment Viterbi over the map records (0x41ef70).
//! 10. [`output`]: traceback, gene assembly and the writers.
//!
//! All messages, including errors, go to stdout, as in the original. Every
//! fatal error exits with code 1.
//!
//! Developer hooks (not part of the original): `--dump-model-json MODEL` and
//! `--dump-input ...` (see [`dump`]); environment variables `GMHMM_RS_TIMING`
//! (stage timings on stderr), `GMHMM_RS_NO_PRUNE` (disable the BP pruning in
//! the Viterbi) and `GMHMM_RS_BPDUMP=<dir>` (raw dump of the BP data).

mod bp;
mod cli;
mod dump;
mod evidence;
mod model;
mod output;
mod sequence;
mod sitemap;
mod tracks;
mod viterbi;

use std::io::Write;

/// ctime(3) of the current local time, including the trailing newline.
/// The original calls `time` + `ctime` in every writer; the date line is the
/// only part of the output that differs between runs.
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

/// Write `s` to stdout (C stdio semantics: write errors are ignored, unlike
/// Rust's `print!`, which panics on them). Bytes that were not valid UTF-8
/// in the input are written back as the original bytes
/// ([`sequence::encode`]).
pub fn out(s: &str) {
    let _ = std::io::stdout().write_all(&sequence::encode(s));
}

/// Print `msg` to stdout (no newline added), flush, and exit with code 1.
/// This is how every fatal error of the original ends.
fn fail(msg: &str) -> ! {
    out(msg);
    let _ = std::io::stdout().flush();
    std::process::exit(1);
}

/// The port's own error (not in the original, which terminates or computes
/// with undefined memory at this point): `gmhmm-rs: error: MSG` on stderr,
/// exit code 1. stdout is flushed first.
pub fn port_fail(msg: &str) -> ! {
    let _ = std::io::stdout().flush();
    let mut e = std::io::stderr();
    let _ = e.write_all(&sequence::encode(&format!("gmhmm-rs: error: {msg}\n")));
    std::process::exit(1);
}

/// The port's own warning on stderr, `gmhmm-rs: warning: MSG`, printed
/// where the port takes a defined fallback on input for which the original
/// terminates or reads undefined memory. Each distinct message is printed once.
pub fn warn(msg: &str) {
    static SEEN: std::sync::Mutex<Vec<String>> = std::sync::Mutex::new(Vec::new());
    let mut seen = SEEN.lock().unwrap_or_else(|e| e.into_inner());
    if seen.iter().any(|m| m == msg) {
        return;
    }
    seen.push(msg.to_string());
    let mut e = std::io::stderr();
    let _ = e.write_all(&sequence::encode(&format!("gmhmm-rs: warning: {msg}\n")));
}

/// Make stdout behave like C stdio for unusual destinations: SIGPIPE gets
/// its default action again (Rust ignores it, so writing to a closed pipe
/// would fail with EPIPE instead of ending the process as it ends the
/// original), and a closed file descriptor 1 is replaced by /dev/null
/// (printf to a closed stdout fails silently in C).
fn setup_stdout() {
    unsafe {
        libc::signal(libc::SIGPIPE, libc::SIG_DFL);
        if libc::fcntl(1, libc::F_GETFD) == -1 {
            let fd = libc::open(c"/dev/null".as_ptr(), libc::O_WRONLY);
            if fd >= 0 && fd != 1 {
                libc::dup2(fd, 1);
                libc::close(fd);
            }
        }
    }
}

/// Stage timer: prints to stderr only when GMHMM_RS_TIMING is set.
struct Timer(Option<std::time::Instant>);
impl Timer {
    fn new() -> Timer {
        Timer(std::env::var_os("GMHMM_RS_TIMING").map(|_| std::time::Instant::now()))
    }
    fn lap(&mut self, what: &str) {
        if let Some(t) = self.0 {
            eprintln!("[timing] {what}: {:.3} s", t.elapsed().as_secs_f64());
            self.0 = Some(std::time::Instant::now());
        }
    }
}

/// The whole gene-finding run for one command line. Returns the exit code
/// (fatal errors exit directly through [`fail`]).
fn run(argv: &[String]) -> i32 {
    let mut timer = Timer::new();
    let opts = match cli::settings(argv) {
        Ok(s) => s,
        Err(m) => fail(&m),
    };
    let mut seq = match sequence::Sequence::read(&opts.seq, &opts.tag) {
        Ok(q) => q,
        Err(m) => fail(&m),
    };
    // -d evidence: read after the sequence (it needs the sequence id and
    // length) and before the model, so its messages come first.
    if !opts.evidence.is_empty() {
        let mut msgs = String::new();
        let result = evidence::read(&opts.evidence, &seq.seqid, seq.len(), &mut msgs);
        out(&msgs);
        match result {
            Ok(marks) => seq.marks = Some(marks),
            Err(()) => fail(""),
        }
    }
    timer.lap("sequence");
    let model = match model::Model::load(&opts.model) {
        Ok(m) => m,
        Err(e) => {
            out(&e.to_string());
            if let Some(m) = &e.port {
                port_fail(m);
            }
            // errors that the original handles with exit(1) inside the
            // reader do not get the "error in model file" line
            if !e.fatal_exit {
                out(&format!("error in model file {}\n", opts.model));
            }
            fail("");
        }
    };
    let flags = sitemap::Flags { taa: model.taa_on, tag: model.tag_on, tga: model.tga_on, gc_donor: model.gc_on };
    timer.lap("model");
    // Same order as the binary: the site map is built before the
    // sequence-length/order check of the k-mer index (0x4099a0), so a too
    // short sequence fails only here.
    let map = sitemap::build(&seq.codes, seq.marks.as_deref(), &flags);
    if let Err(m) = seq.build_kmer(model.order) {
        fail(&m);
    }
    // 0x41f6d0: branch point data (only for $BP models)
    let mut bp_data = match &model.bp {
        Some(b) => match bp::build(&model, b, &seq.codes, &seq.kmer) {
            Ok(d) => Some(d),
            Err(m) => fail(&m),
        },
        None => None,
    };
    timer.lap("map+kmer+bp");
    let tracks = tracks::build(&map, &seq.codes, &seq.kmer, seq.marks.as_deref(), &model, opts.mask_penalty);
    // The BP acceptor term is relative to the finished site track
    // (0x41f9a0 / 0x41b230 / 0x41fa40 run after 0x41abe0).
    if let (Some(d), Some(b)) = (bp_data.as_mut(), model.bp.as_ref()) {
        bp::acceptor_terms(d, b, &map, &seq.codes, &tracks.site);
    }
    tracks::warn_outside(&opts.model);
    timer.lap("tracks");
    let cells = viterbi::run(&map, &tracks, &model, &seq.codes, bp_data.as_ref(), opts.verbose);
    if let (Ok(dir), Some(d)) = (std::env::var("GMHMM_RS_BPDUMP"), bp_data.as_ref()) {
        dump::bp_bins(&dir, d, &tracks, &cells.bp, &cells.score);
    }
    timer.lap("viterbi");
    let records = output::traceback(&map, &cells.bp, &cells.score);
    // No reachable end state (e.g. a sequence of ORDER+1 or ORDER+2 bases,
    // or a model whose initial probabilities are all 0): the parse is empty.
    // The gtf, gff3 and tr writers handle that; the original's lst writer
    // reads the state record list without checking and terminates (SIGSEGV).
    // The port writes the lst header without genes, as for any parse
    // without genes.
    if records.is_empty() && opts.format != "gff3" && opts.format != "gtf" && opts.format != "tr" {
        warn(&format!(
            "no valid parse of sequence {} with model {} (sequence too short for the model, or no allowed \
             state at the sequence ends); writing a prediction without genes (not handled by the original)",
            opts.seq, opts.model
        ));
    }
    let genes = output::genes(&records, &seq.codes);
    let ctx = output::Ctx {
        seqfile: &opts.seq,
        modfile: &opts.model,
        seq: &seq,
        model: &model,
        min_gene: opts.min_gene,
        trace_z: opts.trace_z,
        protein: opts.protein,
        nucleotide: opts.nucleotide,
        date: ctime_now(),
    };
    // Dispatch as in main: exact string compare, anything else (including
    // accepted oddities such as `-f 'gff3|gtf'`) falls through to lst.
    let text = match opts.format.as_str() {
        "gff3" => output::gff3(&ctx, &records, &genes),
        "tr" => match output::tr(&ctx, &records, &genes) {
            Ok(t) => t,
            Err(output::TrError::Frame) => fail(""),
            Err(output::TrError::Port(m)) => port_fail(&m),
        },
        "gtf" => output::gtf(&ctx, &records, &genes),
        _ => output::lst(&ctx, &records, &genes),
    };
    // An unwritable output file is reported but is not fatal (exit code 0).
    if std::fs::write(sequence::os_path(&opts.out), sequence::encode(&text)).is_err() {
        if opts.format == "tr" {
            out(&format!("error on open file: {}\n", opts.out));
        } else {
            out(&format!("can't open a file {}\n", opts.out));
        }
    }
    // 0x418900: -b intron statistics, only for BP models
    if let Some(d) = bp_data.as_ref() {
        if !opts.bp_out.is_empty() {
            let is_tr = opts.format == "tr";
            let text = bp::intron_stats(d, &model, &map, &tracks.intron, &records, is_tr);
            // the original opens with "a" for tr (and then writes nothing),
            // with "w" otherwise
            let mut open_opts = std::fs::OpenOptions::new();
            if is_tr {
                open_opts.append(true).create(true);
            } else {
                open_opts.write(true).create(true).truncate(true);
            }
            match open_opts.open(sequence::os_path(&opts.bp_out)) {
                Ok(mut f) => {
                    let _ = f.write_all(text.as_bytes());
                }
                Err(_) => out(&format!("can't open a file {}\n", opts.bp_out)),
            }
        }
    }
    // 0x4153c0: -r, printed after the output file has been written
    if opts.best_path {
        out(&format!("# best path probability: {}\n", model::c_fmt_f(viterbi::best_path(&cells.score))));
    }
    0
}

fn main() {
    setup_stdout();
    // Arguments as raw bytes, like C's argv; bytes that are not valid UTF-8
    // are carried through the program in escaped form (`sequence::decode`).
    // (std::env::args() would panic on them.)
    let args: Vec<String> = {
        use std::os::unix::ffi::OsStringExt;
        std::env::args_os().map(|a| sequence::decode(&a.into_vec())).collect()
    };
    // developer hooks
    if args.len() == 3 && args[1] == "--dump-model-json" {
        match model::Model::load(&args[2]) {
            Ok(m) => println!("{}", dump::model_json(&m)),
            Err(e) => {
                out(&e.to_string());
                if let Some(m) = &e.port {
                    port_fail(m);
                }
                std::process::exit(1);
            }
        }
        return;
    }
    if args.len() >= 5 && args[1] == "--dump-input" {
        dump::input_bins(&args[2..]);
        return;
    }
    let rc = run(&args);
    let _ = std::io::stdout().flush();
    std::process::exit(rc);
}
