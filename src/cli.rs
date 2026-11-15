//! Command line: parser (0x4075e0), post-processing and validation
//! (0x407c20) and the usage text (0x4079e0). First stage of the pipeline;
//! the resulting [`Settings`] drive every later stage.
//!
//! The parser is a hand-written loop over `argv`, not getopt, and has several
//! details that are reproduced on purpose (re-notes/input.md section 1):
//! options must be exactly two characters, `-r` checks the `-n` flag for
//! duplicates, `-b` as last argument ends parsing with success, `-k`/`-w`
//! accept values starting with '-' and may be repeated, and the format check
//! is a substring test. All messages go to stdout; every error exits with 1.

/// Parsed command line (settings object in main, 0x128 bytes).
/// Empty strings mean "option not given".
#[derive(Debug, Default, Clone)]
pub struct Settings {
    /// positional argument: FASTA sequence file (+0x00)
    pub seq: String,
    /// `-m`: model file (+0x20)
    pub model: String,
    /// `-o`: output file; defaults to `<seq>.lst` for every format (+0x40)
    pub out: String,
    /// `-d`: external evidence GFF ("GeneMark.hmm plus" hints) (+0x60)
    pub evidence: String,
    /// `-b`: intron / branch point statistics file; only written for `$BP`
    /// models (+0x80)
    pub bp_out: String,
    /// `-p`: protein translations (only the lst writer uses it) (+0xa0)
    pub protein: bool,
    /// `-n`: nucleotide sequences of genes (only the lst writer) (+0xa1)
    pub nucleotide: bool,
    /// `-f`: output format, "lst" (default), "gff3", "gtf" or "tr" (+0xa8)
    pub format: String,
    /// `-s`: sequence tag replacing the FASTA-derived sequence id (+0xc8)
    pub tag: String,
    /// `-T`: undocumented; parsed and stored but never read (+0xe8)
    pub t_opt: String,
    /// `-r`: print the best path log-probability to stdout (+0x108)
    pub best_path: bool,
    /// `-k`: soft-mask penalty, added (as given, usually negative) per
    /// masked base of coding segments; needs `-d` with `mask` lines (+0x110)
    pub mask_penalty: f64,
    /// `-z`: "trace back": use defline token 2 as sequence id and shift all
    /// coordinates by token 3 - 1 (+0x118)
    pub trace_z: bool,
    /// `-w`: minimum gene length for output (genes with evidence are always
    /// kept); never enters the Viterbi (+0x11c)
    pub min_gene: i32,
    /// `-v`: progress messages from the Viterbi (+0x120)
    pub verbose: bool,
}

/// Usage text printed when the program is called without arguments
/// (0x4079e0, exit code 1). Byte-exact, including the 5-space indent of the
/// developer options.
pub const USAGE: &str = "GeneMark.hmm eukaryotic, version 3.68

Usage: gmhmme3 [options] <sequence file>

  required parameters:
    -m <model file>

  optional parameters:
    -o <output file>
    -p write protein translation
    -n write nucleotide sequence
    -b <output file> output statistics of predicted introns
    -d <file name> provide input for GeneMark.hmm plus
    -s <string> sequence tag in GFF output format
    -f <format> output prediction in [lst|gff3|gtf] format; default [lst]
    -k <number> value for soft-mask penalty

developer options:
     -z trace back seqid and position
     -w <number> minimum gene length
     -r report best path probability
     -v verbose
";

/// Option loop of 0x4075e0. Returns false for any syntax error (the caller
/// then prints "error in command line").
fn parse(s: &mut Settings, argv: &[String]) -> bool {
    let argc = argv.len();
    let mut i = 1;
    while i < argc {
        let a = argv[i].as_bytes();
        // anything not starting with '-' is the (single) sequence file
        if a.first() != Some(&b'-') {
            if !s.seq.is_empty() {
                return false;
            }
            s.seq = argv[i].clone();
            i += 1;
            continue;
        }
        // options are exactly "-x"; "-mx" or "--m" are errors
        if a.len() != 2 {
            return false;
        }
        let next = argv.get(i + 1);
        // String options: error on a duplicate, on a missing value and on a
        // value starting with '-'. An empty value ("") is stored as empty,
        // so the option then counts as not given.
        let string_opt = |field: &mut String, next: Option<&String>| -> bool {
            if !field.is_empty() {
                return false;
            }
            match next {
                Some(n) if !n.starts_with('-') => {
                    *field = n.clone();
                    true
                }
                _ => false,
            }
        };
        match a[1] {
            b'm' => {
                if !string_opt(&mut s.model, next) {
                    return false;
                }
                i += 2;
            }
            b'o' => {
                if !string_opt(&mut s.out, next) {
                    return false;
                }
                i += 2;
            }
            b'd' => {
                if !string_opt(&mut s.evidence, next) {
                    return false;
                }
                i += 2;
            }
            b'f' => {
                if !string_opt(&mut s.format, next) {
                    return false;
                }
                i += 2;
            }
            b's' => {
                if !string_opt(&mut s.tag, next) {
                    return false;
                }
                i += 2;
            }
            b'T' => {
                if !string_opt(&mut s.t_opt, next) {
                    return false;
                }
                i += 2;
            }
            b'b' => {
                if !s.bp_out.is_empty() {
                    return false;
                }
                match next {
                    // as in the original: `-b` as the last argument is accepted without
                    // a value (parsing ends with success)
                    None => return true,
                    // as in the original: a following option is not consumed, `-b`
                    // then takes no value
                    Some(n) if n.starts_with('-') => i += 1,
                    Some(n) => {
                        s.bp_out = n.clone();
                        i += 2;
                    }
                }
            }
            b'p' => {
                if s.protein {
                    return false;
                }
                s.protein = true;
                i += 1;
            }
            b'n' => {
                if s.nucleotide {
                    return false;
                }
                s.nucleotide = true;
                i += 1;
            }
            b'r' => {
                // the binary tests the -n byte here, not its own: `-n -r`
                // fails, `-r -n` and `-r -r` are accepted
                if s.nucleotide {
                    return false;
                }
                s.best_path = true;
                i += 1;
            }
            b'v' => {
                if s.verbose {
                    return false;
                }
                s.verbose = true;
                i += 1;
            }
            b'z' => {
                if s.trace_z {
                    return false;
                }
                s.trace_z = true;
                i += 1;
            }
            // -k and -w: no '-' check (gmes passes `-k -0.03`), no duplicate
            // check (last one wins), C prefix parsing (non-numeric text gives 0)
            b'k' => match next {
                None => return false,
                Some(n) => {
                    s.mask_penalty = crate::model::strtod(n);
                    i += 2;
                }
            },
            b'w' => match next {
                None => return false,
                Some(n) => {
                    s.min_gene = crate::model::strtol(n) as i32;
                    i += 2;
                }
            },
            _ => return false,
        }
    }
    true
}

/// 0x407c20: parse and validate the command line, fill in defaults.
/// `Err(message)` means: print the message to stdout and exit 1. The checks
/// run in the order of the original, so the first failing one is reported.
pub fn settings(argv: &[String]) -> Result<Settings, String> {
    if argv.len() == 1 {
        return Err(USAGE.to_string());
    }
    let mut s = Settings::default();
    if !parse(&mut s, argv) {
        return Err("gmhmme3 : error in command line\n".into());
    }
    if s.seq.is_empty() {
        return Err("gmhmme3 : error, sequence file's name is not provided\n".into());
    }
    // default output name, regardless of -f (so -f gtf writes GTF to .lst)
    if s.out.is_empty() {
        s.out = format!("{}.lst", s.seq);
    }
    if s.model.is_empty() {
        return Err("gmhmme3 : error, model file's name is not provided\n".into());
    }
    // exact string equality only; different spellings of one path pass
    if s.out == s.seq || s.out == s.model || s.out == s.evidence {
        return Err("gmhmme3 : error, check file names\n".into());
    }
    if s.format.is_empty() {
        s.format = "lst".into();
    }
    // substring test as in the original: "gff3|gtf" is accepted here and
    // later falls through to the lst writer
    if !"|lst|gff3|gtf|tr|".contains(&format!("|{}|", s.format)) {
        return Err("gmhmme3 : error, unknown output format found\n".into());
    }
    Ok(s)
}
