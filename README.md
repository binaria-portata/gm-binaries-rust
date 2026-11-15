# gm-binaries-rust

> **Unofficial software.** This is an independent re-implementation of three
> programs that ship as binaries with GeneMark-ETP. It is not developed,
> reviewed, endorsed or supported by the GeneMark authors or by the Georgia
> Institute of Technology. Please do not send questions or bug reports about
> these programs to the GeneMark team; use this repository's issue tracker.
> The name GeneMark is used here only to say which programs these replace.
>
> **Non-commercial use only.** The original programs are distributed under the
> Creative Commons Attribution-NonCommercial-ShareAlike 4.0 International
> license (CC BY-NC-SA 4.0). This re-implementation is an adaptation of them
> and is distributed under the same license. Read
> [License and terms of use](#license-and-terms-of-use) before using,
> copying or sharing it.

Rust re-implementations of three programs that ship as binaries with
GeneMark-ETP. They are meant as drop-in replacements: same command line, same
predictions, same output files, same messages and exit codes on the tested
inputs. On a few inputs that the original programs do not handle, these
programs print a warning or an error instead (see
[Known differences](#known-differences-from-the-originals)).

All credit for the gene-finding methods, the models and the algorithms goes
to the GeneMark authors. This repository only re-implements the existing
programs in another language. If you use them, please cite the original
publications listed under [References](#references).

| Binary | Original | Status |
| --- | --- | --- |
| `gmhmm-rs` | `gmhmme3` 3.68, GeneMark.hmm eukaryotic (GeneMark-ES/ET/EP/ETP gene finder) | complete |
| `gmhmmp-rs` | `gmhmmp` 3.20, GeneMark.hmm prokaryotic (GeneMarkS-T gene finder) | complete |
| `probuild-rs` | `probuild` 2.52 (sequence/model utility) | complete |

Build: `cargo build --release` (binaries in `target/release/`).

## gmhmm-rs (gmhmme3)

Covers the command line (all options), the model file reader, FASTA input,
`-d` evidence (all features, marks, `mask` + `-k`), score tracks, Viterbi
and traceback, branch-point models (`$BP`, GeneMark-ES `--fungus`) and all
output formats (lst with `-p`/`-n`, gtf, gff3, tr, `-r`, `-w`, `-s`, `-z`,
`-b`). The Viterbi follows the original order of operations, so every
floating-point result is the same.

Developer switches: `GMHMM_RS_TIMING=1` prints stage timings to stderr;
`GMHMM_RS_NO_PRUNE=1` disables the branch-point candidate pruning (for
checking that pruning does not change any result).

## gmhmmp-rs (gmhmmp)

Covers single and meta (`__GC30..70`) models, RBS in spacer and context mode
(`-r`), overlap (`-p`), strand (`-s`), gap split (`-n`), defline parameters
(`-c`), `-b`, `-k`/`-K`, `-A`/`-D`, LST and GFF, logging (`-v`, `-V`, `-L`,
`-H`) and external evidence (`-e`/`-E`). The documented `-a`/`-d` options are
rejected by the original and by the port alike; `-u` is not supported.

## probuild-rs (probuild)

Covers all modes of the usage text: the sequence utilities used by the
GeneMark pipelines (`--reformat_fasta`, `--split`, `--stat`, `--stat_fasta`,
`--substring_n_distr`, `--clean_join`, `--cut`, ...) and the model modes
(`--mkmod_euk`, `--mkmod_non`, `--mkmod`, `--parse`, `--mkchain`,
`--intergenic`, `--pre_start`, `--compare`, `--orf`, `--ORF_PER_SEQ`,
`--gibbs`, `--gc_hmm`, `--combine`, ...), including the parameter file
(`--par`), messages and output buffering.

## Verification

Each test runs a program on the same input as the original and compares
stdout, stderr, exit code and every output file byte for byte.

- `tests/compare.sh` (gmhmme3, `tests/cases.tsv`, 96 cases) and
  `SUITE=gmhmmp tests/compare.sh` (gmhmmp, `tests/gmhmmp/cases.tsv`, 60
  cases): against stored outputs in `tests/oracle/` and
  `tests/gmhmmp/oracle/` (date lines masked). All pass. For the cases listed
  under Known differences, the stored output is the port's.
- `tests/fuzz.py`, `tests/fuzz_gmhmmp.py`: randomised differential testing
  against the original (random genome and transcript segments, N runs, IUPAC
  codes, case, CRLF, heuristic models, random evidence and options).
- `tests/cmp_model_dump.py`: compares every transformed gmhmme3 model value
  with a debugger dump of the original's model object.
- `tests/probuild/run_cases.py --ref PROBUILD --rs target/release/probuild-rs
  --cases FILE`: `tests/probuild/cases.tsv` (253 calls) and
  `tests/probuild/model_cases.tsv` (415 calls), each run in a fresh directory
  next to the original. All identical except the 5 cases listed under Known
  differences.
- `tests/probuild/replay.py --ref ORIGINAL --rs PORT CAPDIR...`: replays
  calls recorded during real pipeline runs (argv, working directory, input
  files) through both programs and compares all outputs. Works for all three
  binaries. Recorded runs of GeneMark-ETP, GeneMark-ES, ES `--fungus` and
  GeneMarkS-T replay without a difference.

Regenerate an oracle with `GMHMME3_REF=/path/to/gmhmme3 tests/make_oracle.sh`
(gmhmmp: `SUITE=gmhmmp GMHMME3_REF=/path/to/gmhmmp tests/make_oracle.sh`).

## Known differences from the originals

This section lists the modifications from the original material, as the
license requires (Section 3(a)(1)(B)). It describes only the three binaries
named under [Source of the reconstruction](#source-of-the-reconstruction),
as observed in the tests of this repository; other GeneMark versions may
behave differently.

On all tested inputs that the original programs handle, the port gives the
same output. For a few inputs that the original programs do not handle
(mostly unusual model values, very short sequences, evidence at the very
ends of a sequence, empty gene sets), the port:

- continues where a sensible result exists, and prints one
  `<program>-rs: warning: ... (not handled by the original)` line per cause
  to stderr, saying what was requested, what is missing in which input, and
  what it does instead;
- otherwise stops with `<program>-rs: error: ... (not handled by the
  original)` on stderr and exit code 1.

gmhmm-rs (gmhmme3):

- Evidence (`-d`) marks that would fall just outside the sequence (hints
  starting at the first or ending at the last base) are dropped, with a
  warning; the rest of the hint is used. A directory given to `-d` gives
  "error, opening file".
- Sequences of `ORDER+1` or `ORDER+2` bases: the header without genes, with
  a warning.
- Site windows (also `tr` windows) reaching outside the sequence: bases
  there count as N, with a warning; windows 2^30 or more bases outside:
  error.
- Model files with values outside their valid range (negative or
  overflowing widths and orders, `_MAX` <= 0, margins outside the window):
  error naming the key. `_MIN` <= 0 in a length distribution and
  `BP_ACC_MIN`/`DON_BP_MIN` <= 0: same output, with a warning.
- Bases past the sequence end read as `A`. On `MAX < MIN` in a length
  distribution the port prints `0x0` where the original prints an address.
  The `-v` progress line format is inferred.

gmhmmp-rs (gmhmmp):

- An N inside a noncoding gap of 1 to ORDM bases between a stop and a start
  uses the Markov table's N value.
- Error for: a model key without a value at the end of the model text;
  MINC > MAXC + 2; ORDM >= 14; RBSM with 4 * RBSS above the int range;
  MINC < 7 when a gene shorter than 7 bp occurs; option `-u`.
- Warning, with a defined result, for: a model without COD1, COD2 or NONC
  after its last ORDM (all-zero matrix); MINC < ORDM + 9 (Markov windows
  outside the sequence score as N); `-r` with an RBS context window leaving
  the sequence, RBSS < -RBSL, or a spacer-mode RBSS above 13 (no RBS score
  for the affected starts).
- Very large MAXC/MAXN values allocate as much memory as in the original.

probuild-rs (probuild):

- Error for: `--mkmod --fixmotif` without any pre-start sequence (for
  example GeneMarkS-T on very few transcripts); a `--mask` interval with its
  end before its start; a `--reformat_fasta` `--mask_margin` out of range;
  `--orf` margins moving an ORF 2^30 or more bases outside the sequence;
  precisions around 2^31; `--gc_hmm` on a record without sequence letters.
- Defined results, with a warning where it applies: `--orf` on an empty
  sequence reports no ORFs; `--reformat_fasta` without `--in` lists the
  options with an empty `--in`; positions outside the sequence (from
  margins above the int range) count as N or are left out; an invalid
  `--RANDOM_NONCOD_GC` gives GC 0.
- Bytes >= 0x80 in sequences and output files of runs that end with an
  error may differ. `--count_substring_atcg` on lower case input reports
  `Out of memory!` at once.

## Source of the reconstruction

- Binaries from GeneMark-ETP as packaged in the `teambraker/braker3`
  container (v3.0.10): `gmhmme3` BuildID
  77f1ce8f11da86a48e2adf5c7b88871217402b03, `gmhmmp` BuildID
  5107402e52b7b0a67d36869b4c5daf557f678ff1, `probuild` BuildID
  e606888692eb361e968ed00705d8ad38834edab1. These binaries are distributed
  under CC BY-NC-SA 4.0.
- Method: comparison of inputs and outputs, and analysis of the binaries.
  Comments in the source refer to addresses in these binaries and to
  analysis notes (`re-notes/`) that are not part of this repository.

## License and terms of use

This section summarises the license. It is not legal advice. Where it and
`LICENSE` differ, `LICENSE` applies.

### The original programs

`gmhmme3`, `gmhmmp` and `probuild` are part of GeneMark-ETP
(<https://github.com/gatech-genemark/GeneMark-ETP>), developed by the
Borodovsky lab at the Georgia Institute of Technology, Atlanta, Georgia, USA.
GeneMark-ETP is distributed under CC BY-NC-SA 4.0. The GeneMark-ES suite
README that ships with GeneMark-ETP names the Georgia Institute of Technology
as copyright holder. The publications describing the programs are listed
under [References](#references).

### This re-implementation

- It is an adaptation of the original programs in the sense of CC BY-NC-SA
  4.0 and is distributed under the same license, as Section 3(b)
  (ShareAlike) requires. `LICENSE` is the license text exactly as it ships
  with GeneMark-ETP.
- The Rust code was written for this repository from the observed behaviour
  of the binaries and from analysis of the binaries
  ([Source of the reconstruction](#source-of-the-reconstruction)). The
  GeneMark source code was not available and was not used. This repository
  contains no GeneMark binaries.
- Changes from the original material (Section 3(a)(1)(B)): the programs are
  rewritten in another language. On a few inputs that the originals do not
  handle, these programs print a warning or an error instead ([Known differences](#known-differences-from-the-originals)).
  Everywhere else they aim to give identical output. That is a test result
  for the inputs in `tests/`, not a guarantee.
- Some files in `tests/` (model files and the expected outputs in the oracle
  directories) come from GeneMark-ETP or were produced by the original
  programs. They are under the same license and attribution.

### What the license allows

You may use, copy, modify and share these programs and adaptations of them
for non-commercial purposes, if you:

- give credit to the original authors and the Georgia Institute of
  Technology, keep the copyright and license notices, and link to the
  original GeneMark-ETP repository;
- state that this is a modified re-implementation, and state your own
  changes if you make any;
- share adaptations under CC BY-NC-SA 4.0.

### What the license does not allow or grant

- **Commercial use.** The license defines NonCommercial as "not primarily
  intended for or directed towards commercial advantage or monetary
  compensation" (Section 1(k)). This repository cannot grant any other
  license. For commercial use of GeneMark, contact the copyright holder of
  the original programs.
- **Endorsement.** You may not state or imply that your use, or results from
  these programs, are connected with, sponsored or endorsed by the GeneMark
  authors or the Georgia Institute of Technology (Section 2(a)(6)). When you
  publish results, name this re-implementation and its version or commit,
  and cite the original GeneMark publications for the methods.
- **Patent and trademark rights.** These are not licensed (Section 2(b)(2)).
- **Warranty.** There is none. The material is provided as-is, and liability
  is limited as far as the law allows (Section 5).

## References

If you use these programs, please cite the original GeneMark publications.

GeneMark.hmm eukaryotic (`gmhmme3`) and the GeneMark-E* pipelines:

- GeneMark-ES: Lomsadze A, Ter-Hovhannisyan V, Chernoff YO, Borodovsky M.
  Gene identification in novel eukaryotic genomes by self-training
  algorithm. Nucleic Acids Res. 2005;33(20):6494-6506.
  doi:10.1093/nar/gki937
- GeneMark-ES `--fungus` (branch point model): Ter-Hovhannisyan V,
  Lomsadze A, Chernoff YO, Borodovsky M. Gene prediction in novel fungal
  genomes using an ab initio algorithm with unsupervised training. Genome
  Res. 2008;18(12):1979-1990. doi:10.1101/gr.081612.108
- GeneMark-ET: Lomsadze A, Burns PD, Borodovsky M. Integration of mapped
  RNA-Seq reads into automatic training of eukaryotic gene finding
  algorithm. Nucleic Acids Res. 2014;42(15):e119. doi:10.1093/nar/gku557
- GeneMark-EP+: Brůna T, Lomsadze A, Borodovsky M. GeneMark-EP+: eukaryotic
  gene prediction with self-training in the space of genes and proteins.
  NAR Genom Bioinform. 2020;2(2):lqaa026. doi:10.1093/nargab/lqaa026
- GeneMark-ETP: Brůna T, Lomsadze A, Borodovsky M. GeneMark-ETP
  significantly improves the accuracy of automatic annotation of large
  eukaryotic genomes. Genome Res. 2024;34(5):757-768.
  doi:10.1101/gr.278373.123

GeneMark.hmm prokaryotic (`gmhmmp`) and GeneMarkS-T:

- GeneMark.hmm: Lukashin AV, Borodovsky M. GeneMark.hmm: new solutions for
  gene finding. Nucleic Acids Res. 1998;26(4):1107-1115.
  doi:10.1093/nar/26.4.1107
- Heuristic models: Besemer J, Borodovsky M. Heuristic approach to deriving
  models for gene finding. Nucleic Acids Res. 1999;27(19):3911-3920.
  doi:10.1093/nar/27.19.3911
- GeneMarkS (RBS models): Besemer J, Lomsadze A, Borodovsky M. GeneMarkS: a
  self-training method for prediction of gene starts in microbial genomes.
  Implications for finding sequence motifs in regulatory regions. Nucleic
  Acids Res. 2001;29(12):2607-2618. doi:10.1093/nar/29.12.2607
- MetaGeneMark (GC-binned meta models): Zhu W, Lomsadze A, Borodovsky M. Ab
  initio gene identification in metagenomic sequences. Nucleic Acids Res.
  2010;38(12):e132. doi:10.1093/nar/gkq275
- GeneMarkS-T: Tang S, Lomsadze A, Borodovsky M. Identification of protein
  coding regions in RNA transcripts. Nucleic Acids Res. 2015;43(12):e78.
  doi:10.1093/nar/gkv227
