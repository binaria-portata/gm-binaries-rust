//! probuild-rs: re-implementation of the GeneMark utility probuild 2.52.
//!
//! probuild is a multi-purpose command line tool of the GeneMark pipelines
//! (GeneMark-ES/ET/EP/ETP, GeneMarkS-T).  One invocation performs exactly one
//! "mode", selected by a mode key such as `--reformat_fasta`, `--split`,
//! `--stat` or `--mkmod_euk`.  All options are `--key [value]` pairs and are
//! collected into one parameter map (`params::Params`, the original's
//! `std::map<string,string>`), optionally merged with a `--par` file.
//!
//! This module is the entry point (`main` 0x402ef5 of the original, a static
//! non-PIE x86-64 binary, so the addresses quoted throughout the port are
//! stable provenance).  Module map:
//!
//! - `params`: parameter map, argv and `--par` parsing, typed getters.
//! - `io`: emulation of the C/C++ runtime (stdout buffering, `ofstream`,
//!   exceptions, abnormal terminations, whole-file loader, `getline` loops).
//! - `fasta`: the two FASTA reader classes (`SeqUtil`, `FastaX`), their
//!   letter translation tables and the report-section writer.
//! - `seqmodes`: modes working on `SeqUtil` (`--reformat_fasta`, `--split`,
//!   `--cut`, `--stat`, `--gc`, `--clean_join`, ...).
//! - `statmodes`: modes working on `FastaX` (`--clean_fasta`, `--stat_fasta`,
//!   `--CLEAN_CODING`, ...).
//! - `pm`: parameter access and number formatting for the model modes.
//! - `model` and the remaining modules: the model-building and gene-set
//!   modes (`--mkmod*`, `--parse`, `--gibbs`, `--compare`, ...).
//!
//! Specifications: `re-notes/probuild/core.md` (parameter object, dispatch,
//! FASTA classes, utility modes) and `re-notes/probuild/models.md`.

pub mod compare;
pub mod fasta;
pub mod gchmm;
pub mod genemodes;
pub mod geneset;
pub mod gibbs;
pub mod gp;
pub mod markov;
pub mod mkmod;
pub mod orf;
pub mod pm;
pub mod io;
pub mod model;
pub mod params;
pub mod seqmodes;
pub mod statmodes;

use params::{Bytes, Params};

/// Usage text printed to stdout when probuild is called without any argument
/// (0x40130e).  Stored verbatim in `usage.txt`, which is output text and must
/// not be edited.
const USAGE: &str = include_str!("usage.txt");

/// Mode keys in the order main (0x402ef5) tests them.  The first key present
/// in the parameter map wins, so e.g. `--stat --gc` runs `--gc`.  The order is
/// part of the observable behaviour and must not be changed.
pub const MODE_KEYS: &[&str] = &[
    "--clean_fasta", "--stat_fasta", "--clean_join", "--gc", "--stat", "--clean", "--split",
    "--parse", "--intergenic", "--mkchain", "--mkmod", "--mkmod_euk", "--mkmod_non",
    "--CLEAN_CODING", "--compare", "--pre_start", "--gibbs", "--mask", "--combine",
    "--STAT_NON_VS_ALL", "--ORF_PER_SEQ", "--orf", "--GENE_SET_STAT", "--count_substring_n",
    "--substring_n_distr", "--count_substring_atcg", "--split_fasta", "--reformat_fasta",
    "--gc_hmm", "--select_fasta", "--cut",
];

/// Program body of main (0x402ef5).  `args` is the full argv including the
/// program name.  Returns the process exit code: 1 for the usage text, else
/// the result of `dispatch` as mapped by the exception handler.
pub fn run(args: &[Bytes]) -> i32 {
    // argc == 1: usage and exit code 1, outside the try block.
    if args.len() == 1 {
        io::out(USAGE.as_bytes());
        return 1;
    }
    let args = args.to_vec();
    // Everything else runs inside main's try/catch: a C++ exception thrown by
    // a handler ends the program with "Out of memory!" or "error" and code 1.
    io::catch_main(move || dispatch(&args))
}

/// Parse the command line and run the handler of the first mode key present.
/// The handler's `true` becomes exit code 0, `false` exit code 1.
fn dispatch(args: &[Bytes]) -> i32 {
    let params = match Params::parse(args) {
        Some(p) => p,
        None => return 1,
    };
    let mode_key = MODE_KEYS.iter().find(|k| params.has(k));
    let ok = match mode_key.copied() {
        None => {
            io::err(b"error in command line - action key is not found\n");
            false
        }
        Some("--clean_fasta") => {
            // Only this mode has its own help text; it is printed instead of
            // running the mode and still gives exit code 1.
            if params.has("--help") {
                io::out(
                    b"Usage: --clean_fasta [filename] required input sequence FASTA\n  \
--clean_fasta_out [filename] required output sequence FASTA\n  \
--clean_fasta_index [filename] required outpt table with new vs old FASTA deflines\n  \
--clean_fasta_label [label] required label for new FASTA defline\n",
                );
                false
            } else {
                statmodes::clean_fasta(&params)
            }
        }
        // main builds FastaX(1) instead of FastaX(0) when --allow_x is
        // present; the flag decides the FastaX table for these modes.
        Some("--stat_fasta") => statmodes::stat_fasta(&params, params.has("--allow_x")),
        Some("--clean_join") => seqmodes::clean_join(&params),
        Some("--gc") => seqmodes::gc(&params),
        Some("--stat") => seqmodes::stat(&params),
        Some("--clean") => seqmodes::clean(&params),
        Some("--split") => seqmodes::split(&params),
        Some("--parse") => model::run_parse(&params),
        Some("--intergenic") => model::run_intergenic(&params),
        Some("--mkchain") => model::run_mkchain(&params),
        Some("--mkmod") => model::run_mkmod(&params),
        Some("--mkmod_euk") => model::run_mkmod_euk(&params),
        Some("--mkmod_non") => model::run_mkmod_non(&params),
        Some("--CLEAN_CODING") => statmodes::clean_coding(&params),
        Some("--compare") => model::run_compare(&params),
        Some("--pre_start") => model::run_pre_start(&params),
        Some("--gibbs") => model::run_gibbs(&params),
        Some("--mask") => seqmodes::mask(&params),
        Some("--combine") => model::run_combine(&params),
        Some("--STAT_NON_VS_ALL") => model::run_stat_non_vs_all(&params),
        Some("--ORF_PER_SEQ") => model::run_orf_per_seq(&params),
        Some("--orf") => model::run_orf(&params),
        Some("--GENE_SET_STAT") => model::run_gene_set_stat(&params),
        Some("--count_substring_n") => statmodes::count_substring_n(&params, params.has("--allow_x")),
        Some("--substring_n_distr") => statmodes::substring_n_distr(&params, params.has("--allow_x")),
        Some("--count_substring_atcg") => statmodes::count_substring_atcg(&params),
        Some("--split_fasta") => seqmodes::split_fasta(&params),
        Some("--reformat_fasta") => seqmodes::reformat_fasta(&params),
        Some("--gc_hmm") => model::run_gc_hmm(&params),
        Some("--select_fasta") => seqmodes::select_fasta(&params),
        Some("--cut") => seqmodes::cut(&params),
        Some(_) => unreachable!(),
    };
    if ok {
        0
    } else {
        1
    }
}

/// Process entry point used by `src/bin/probuild-rs.rs`.  Arguments are taken
/// as raw bytes (no UTF-8 requirement, like C's `argv`).  The emulated stdout
/// buffer is written out before exiting, as glibc does at normal termination.
pub fn main() {
    use std::os::unix::ffi::OsStringExt;
    let args: Vec<Bytes> = std::env::args_os().map(|a| a.into_vec()).collect();
    let rc = run(&args);
    io::flush_stdout();
    std::process::exit(rc);
}
