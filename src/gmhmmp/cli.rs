//! gmhmmp command line: 0x406ef0 (flow), 0x405950 (defaults), 0x405ac0
//! (parser), 0x406d50 (validation), 0x4052e0 (usage), 0x405520 (dump).
//! See re-notes/gmhmmp/input.md section 1. All messages go to stdout.
//!
//! Behaviour of the original kept here: the documented `-a`/`-d` are rejected
//! (only `-A file`/`-D file` exist), `-s` validates the CURRENT value
//! instead of the new one, option arguments may not start with '-',
//! `-K` implies `-k`, no `-o` forces LST output to `<seq>.lst`, and `-V`
//! prints the settings dump twice (before and after validation).

use super::num::{c_f, strtod, strtol_i32};

/// Settings object (in main, offsets in input.md 1.2). Several
/// fields are only parsed and printed in the `-V`/`-H` dump.
#[derive(Debug, Clone)]
pub struct Settings {
    /// positional: FASTA file
    pub seq: String,
    /// -m model file (single or meta)
    pub model: String,
    pub gc_model: String,  // -M (only checked for uniqueness)
    pub rbs_model: String, // -R (same)
    /// -o output file (default `<seq>.lst`)
    pub out: String,
    pub format: u8, // b'L' / b'G'
    pub print_aa: bool,
    pub print_nt: bool,
    pub motif_info: bool, // -k (and -K)
    pub motif_pos: bool,  // -K
    pub best: bool,       // -b
    /// -A protein file, -D nucleotide file, -K RBS score file
    pub aa_file: String,
    pub nt_file: String,
    pub motif_file: String,
    /// -v: log to stdout (only the evidence code logs)
    pub verbose: bool,
    /// -V: settings dumps and log info; forces -v
    pub debug: bool,
    /// -H: settings dump, exit 1
    pub help: bool,
    /// -L log file
    pub log_file: String,
    /// -e evidence GFF
    pub evidence: String,
    /// -E evidence score (default 100.0)
    pub evidence_score: f64,
    /// -T, -S: stored and printed only
    pub key_path: String,
    pub shape: String,
    /// -s first char of the argument: b'd', b'r', anything else = both
    pub strand: u8,
    /// -c: `[shape=...]` and `[gc=...]` from the defline
    pub parse_defline: bool,
    /// -i probability p of starting/ending in noncoding state (default 0.5)
    pub ini_noncoding: f64,
    /// -r: use the RBS model
    pub rbs: bool,
    /// -p: overlaps allowed (default true; `-p` with a value parsing to
    /// 0.0 switches it off, any other value leaves it)
    pub overlap: bool,
    pub partial_dur: bool, // -P (inverted in the dump, unused)
    /// -N: never used
    pub mask_lowcase: bool,
    /// 48, 0 with -n; only `!= 0` matters (the gap threshold is 49)
    pub split_gaps: i32,
    /// -x overlap parameter q (default 0.26)
    pub overlap_penalty: f64,
    /// never set, never used
    pub gc_window: i32,
    pub model_type: i32, // 1 single, 2 meta
    /// -h minimum "heuristic" gene length (default 0)
    pub min_heu_len: i32,
    /// -u three coding models (terminates the original; refused here)
    pub three_coding: bool,
    /// argv joined with TAB, for the dump
    pub cmd_line: String,
}

impl Settings {
    /// 0x405950
    fn defaults() -> Settings {
        Settings {
            seq: String::new(),
            model: String::new(),
            gc_model: String::new(),
            rbs_model: String::new(),
            out: String::new(),
            format: b'L',
            print_aa: false,
            print_nt: false,
            motif_info: false,
            motif_pos: false,
            best: false,
            aa_file: String::new(),
            nt_file: String::new(),
            motif_file: String::new(),
            verbose: false,
            debug: false,
            help: false,
            log_file: String::new(),
            evidence: String::new(),
            evidence_score: 100.0,
            key_path: String::new(),
            shape: "partial".into(),
            strand: b'.',
            parse_defline: false,
            ini_noncoding: 0.5,
            rbs: false,
            overlap: true,
            partial_dur: false,
            mask_lowcase: false,
            split_gaps: 48,
            overlap_penalty: 0.26,
            gc_window: 0,
            model_type: 1,
            min_heu_len: 0,
            three_coding: false,
            cmd_line: String::new(),
        }
    }
}

pub const USAGE: &str = "
GeneMark.hmm version 3.20

Usage:
\tgmhmmp [options] [sequence file]

\t\t[sequence file] - file with sequence/s in FASTA format

Required options:
\t-m\t[model file] file with gene finding parameters
\t\t               metagenome style model file is allowed

Optional parameters:

*** Output options

\t-o\t[output file] output file name
\t\tdefault = [sequence file].lst
\t-f\toutput format: [L] for LST, [G] for GFF
\t\tdefault = L
\t-a\toutput protein sequences of predicted genes
\t-d\toutput nucleotide sequences of predicted genes
\t-k\toutput RBS score and spacer
\t-K\t[file] output all RBS scores and location to file
\t-b\toutput one best prediction per sequence

*** Prediction options

\t-r\tUse RBS for gene start prediction
\t-s\t[d|r|.] Strand to predict genes on
\t\tdefault = '.' predict genes on both strands
\t-p\t[1|0] Gene overlap are allowed [1] or prohibited [0]
\t\tdefault = 1
\t\t[0] is used for gene finding in intron-less eukaryotic genomes
\t-e\t[filename] use this file with external infromation to profile predictions
\t-c\tparse parameters from FASTA defline

*** Developer options

\t-M\t[model file] file with gene finding parameters for metagenomes
\t-R\t[model file] file with parameters of RBS sub model
\t-L\t[log file] output information into loginfo file
\t-i\t[number] Probability of initiation and termination in non-coding state
\t\tdefault = 0.5
\t-n\tTurn OFF prediction of the incomplete genes at sequence gaps
\t\t(gap : string of unknown letters 'N' longer then 20 nt)
\t-E\t[value] evidence score for coding regions
\t-v\tverbose
\t-V\tdebug mode

Examples:

\tgmhmmp -r -m bsub.mod -o sequence.lst sequence.fasta
\tgmhmmp -r -p 0 -m yeast.mod -o sequence.lst sequence.fasta

";

/// 0x405ac0: false = "error in command line". Options are exactly two
/// characters; string options may be given only once (except -S); a second
/// positional argument fails.
fn parse(s: &mut Settings, argv: &[String]) -> bool {
    let argc = argv.len();
    let mut i = 1;
    while i < argc {
        let a = &argv[i];
        let b = a.as_bytes();
        if b.first() == Some(&b'-') {
            if b.len() != 2 {
                return false;
            }
            let c = b[1];
            let takes_arg = matches!(c, b'm' | b'o' | b'M' | b'R' | b'A' | b'D' | b'K' | b'L' | b'T' | b'e' | b'S' | b'f' | b's' | b'E' | b'i' | b'x' | b'h' | b'p');
            let mut arg = "";
            if takes_arg {
                i += 1;
                if i >= argc {
                    return false;
                }
                arg = &argv[i];
                if arg.is_empty() || arg.starts_with('-') {
                    return false;
                }
            }
            let set_once = |f: &mut String, v: &str| -> bool {
                if !f.is_empty() {
                    return false;
                }
                *f = v.to_string();
                true
            };
            let ok = match c {
                b'm' => set_once(&mut s.model, arg),
                b'o' => set_once(&mut s.out, arg),
                b'M' => set_once(&mut s.gc_model, arg),
                b'R' => set_once(&mut s.rbs_model, arg),
                b'A' => {
                    s.print_aa = true;
                    set_once(&mut s.aa_file, arg)
                }
                b'D' => {
                    s.print_nt = true;
                    set_once(&mut s.nt_file, arg)
                }
                b'K' => {
                    s.motif_pos = true;
                    s.motif_info = true;
                    set_once(&mut s.motif_file, arg)
                }
                b'L' => set_once(&mut s.log_file, arg),
                b'T' => set_once(&mut s.key_path, arg),
                b'e' => set_once(&mut s.evidence, arg),
                b'S' => {
                    s.shape = arg.to_string();
                    true
                }
                b'f' => {
                    let f = arg.as_bytes()[0];
                    if f == b'L' || f == b'G' {
                        s.format = f;
                        true
                    } else {
                        false
                    }
                }
                b's' => {
                    // as in the original: checks the current value, stores any first char
                    if matches!(s.strand, b'd' | b'r' | b'.') {
                        s.strand = arg.as_bytes()[0];
                        true
                    } else {
                        false
                    }
                }
                b'E' => {
                    s.evidence_score = strtod(arg);
                    true
                }
                b'i' => {
                    s.ini_noncoding = strtod(arg);
                    true
                }
                b'x' => {
                    s.overlap_penalty = strtod(arg);
                    true
                }
                b'h' => {
                    s.min_heu_len = strtol_i32(arg);
                    true
                }
                b'p' => {
                    if strtod(arg) == 0.0 {
                        s.overlap = false;
                    }
                    true
                }
                b'b' => {
                    s.best = true;
                    true
                }
                b'c' => {
                    s.parse_defline = true;
                    true
                }
                b'k' => {
                    s.motif_info = true;
                    true
                }
                b'n' => {
                    s.split_gaps = 0;
                    true
                }
                b'r' => {
                    s.rbs = true;
                    true
                }
                b'u' => {
                    s.three_coding = true;
                    true
                }
                b'v' => {
                    s.verbose = true;
                    true
                }
                b'V' => {
                    s.debug = true;
                    true
                }
                b'H' => {
                    s.help = true;
                    true
                }
                b'N' => {
                    s.mask_lowcase = true;
                    true
                }
                b'P' => {
                    s.partial_dur = true;
                    true
                }
                _ => false,
            };
            if !ok {
                return false;
            }
        } else {
            if a.is_empty() || !s.seq.is_empty() {
                return false;
            }
            s.seq = a.clone();
        }
        i += 1;
    }
    true
}

fn tf(b: bool) -> &'static str {
    if b {
        "true"
    } else {
        "false"
    }
}

/// 0x405520 settings dump (-V, -H)
pub fn dump(s: &Settings) -> String {
    let mut o = String::new();
    o.push_str("debug mode\nCLASS SETTINGS\n");
    o.push_str(&format!("cmd line          : {}\n\n", s.cmd_line));
    o.push_str("required parameters\n");
    o.push_str(&format!("     sequence file  :  {}\n", s.seq));
    o.push_str(&format!("  m  model file     :  {}\n", s.model));
    o.push_str("special modes\n");
    o.push_str(&format!("  M  gc model file  :  {}\n", s.gc_model));
    o.push_str(&format!("  R  rbs model file :  {}\n", s.rbs_model));
    o.push_str("output parameters\n");
    o.push_str(&format!("  o  output file    :  {}\n", s.out));
    o.push_str(&format!("  f  output file format     :  {}\n", s.format as char));
    o.push_str(&format!("  a  print aa to output file:  {}\n", tf(s.print_aa)));
    o.push_str(&format!("  d  print nt to output file:  {}\n", tf(s.print_nt)));
    o.push_str(&format!("  b  print best per sequence:  {}\n", tf(s.best)));
    o.push_str(&format!("  k  print motif info to output file:  {}\n", tf(s.motif_info)));
    o.push_str(&format!("  -  print motif position instead of spacer length:  {}\n", tf(s.motif_pos)));
    o.push_str(&format!("  A  aa out file   :  {}\n", s.aa_file));
    o.push_str(&format!("  D  nt out file   :  {}\n", s.nt_file));
    o.push_str(&format!("  K  motif out file:  {}\n", s.motif_file));
    o.push_str(&format!("  v  verbose status:  {}\n", tf(s.verbose)));
    o.push_str(&format!("  V  debug status  :  {}\n", tf(s.debug)));
    o.push_str(&format!("  H  help status   :  {}\n", tf(s.help)));
    o.push_str(&format!("  L  logfile       :  {}\n", s.log_file));
    o.push_str("PLUS interface\n");
    o.push_str(&format!("  e  evidence file :  {}\n", s.evidence));
    o.push_str(&format!("  E  evidence score:  {}\n", c_f(s.evidence_score)));
    o.push_str("miscellaneous\n");
    o.push_str(&format!("  T  key path      :  {}\n", s.key_path));
    o.push_str("algorithm parameters\n");
    o.push_str("sequence treatment\n");
    o.push_str(&format!("  S  sequence shape: {}\n", s.shape));
    o.push_str(&format!("  s  strand to look: {}\n", s.strand as char));
    o.push_str(&format!("  c  parse from fasta defline :  {}\n", tf(s.parse_defline)));
    o.push_str(&format!("  i  ini noncoding probability:  {}\n", c_f(s.ini_noncoding)));
    o.push_str("algorithm switches\n");
    o.push_str(&format!("  r  motif status  :  {}\n", tf(s.rbs)));
    o.push_str(&format!("  p  overlap status:  {}\n", tf(s.overlap)));
    o.push_str(&format!("  P  partial duration status:  {}\n", tf(!s.partial_dur)));
    o.push_str("sequence modification\n");
    o.push_str(&format!("  N  mask lowcase  :  {}\n", tf(s.mask_lowcase)));
    o.push_str(&format!("  n  split at gaps :  {}\n", s.split_gaps));
    o.push_str("better not to change parameters\n");
    o.push_str(&format!("  x  overlap penalty:  {}\n", c_f(s.overlap_penalty)));
    o.push_str(&format!("  -  gc window:  {}\n", s.gc_window));
    o.push_str("input indicators\n");
    o.push_str(&format!("  -  model type    :  {}\n", s.model_type));
    o.push_str("parameters in test mode\n");
    o.push_str(&format!("  h  minimum heu gene length: {}\n", s.min_heu_len));
    o.push_str(&format!("  u  use 3 coding state model: {}\n", tf(s.three_coding)));
    o.push('\n');
    o
}

/// Result of 0x406ef0: Ok(settings, text to print first (the `-V`
/// dumps)), or Err(stdout text) for exit code 1. Validation 0x406d50 in
/// the binary's order: sequence file, model file, default output name,
/// uniqueness of all file names (plain string compare), model type.
pub fn settings(argv: &[String]) -> Result<(Settings, String), String> {
    if argv.len() == 1 {
        return Err(USAGE.to_string());
    }
    let mut s = Settings::defaults();
    s.cmd_line = argv.join("\t");
    if !parse(&mut s, argv) {
        return Err("gmhmmp: error in command line\n".into());
    }
    if s.help {
        return Err(dump(&s));
    }
    let mut pre = String::new();
    if s.debug {
        pre.push_str(&dump(&s));
    }
    // 0x406d50
    if s.debug {
        s.verbose = true;
    }
    if s.seq.is_empty() {
        return Err(pre + "gmhmmp: error, sequence file name is missing\n");
    }
    if std::fs::File::open(&s.seq).is_err() {
        return Err(pre + &format!("gmhmmp: error, on open file: <{}>\n", s.seq));
    }
    if s.model.is_empty() {
        return Err(pre + "gmhmmp: error, model file name is missing\n");
    }
    let mfile = std::fs::read(&s.model);
    if mfile.is_err() && std::fs::File::open(&s.model).is_err() {
        return Err(pre + &format!("gmhmmp: error, on open file <{}>\n", s.model));
    }
    if s.out.is_empty() {
        s.out = format!("{}.lst", s.seq);
        s.format = b'L';
    }
    // 0x406280 uniqueness
    let mut names = vec![s.seq.clone(), s.out.clone(), s.model.clone()];
    for x in [&s.gc_model, &s.rbs_model, &s.evidence, &s.nt_file, &s.aa_file, &s.log_file, &s.motif_file] {
        if !x.is_empty() {
            names.push(x.clone());
        }
    }
    let mut set = names.clone();
    set.sort();
    set.dedup();
    if set.len() != names.len() {
        return Err(pre + "gmhmmp: error, file names on command line should be unique\n");
    }
    // 0x4058b0 model type
    let first = mfile.ok().and_then(|d| d.first().copied());
    s.model_type = if first == Some(b'_') { 2 } else { 1 };
    if s.debug {
        pre.push_str(&dump(&s));
    }
    Ok((s, pre))
}
