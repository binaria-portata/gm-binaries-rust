//! Entry points of the model modes of probuild 2.52 (spec:
//! re-notes/probuild/models.md), called by the dispatcher in `mod.rs` once
//! main has found the mode key.  Each returns the handler result: `true`
//! gives exit code 0, `false` exit code 1 (models.md §0).  Some paths end
//! the process directly (`exit(1)` in the chain code, `io::fatal` where the original terminates).
//!
//! Where the modes live:
//!
//! | mode | handler | module |
//! |---|---|---|
//! | --mkmod_euk, --mkmod_non, --mkmod, --combine | 0x440730, 0x442010, 0x442432, 0x442c58 | `mkmod` (chains in `markov`) |
//! | --parse, --mkchain, --intergenic, --pre_start, --STAT_NON_VS_ALL, --GENE_SET_STAT | 0x4100ec, 0x4109c2, 0x40fb32, 0x40eda4, 0x40f654, 0x409d16 | `genemodes` (training sets in `gp`, gene lists in `geneset`) |
//! | --compare | 0x452264 | `compare` |
//! | --orf, --ORF_PER_SEQ | 0x44a182, 0x44a734 | `orf` |
//! | --gibbs | 0x455cc4 | `gibbs` |
//! | --gc_hmm | 0x457e34 | `gchmm` |
//!
//! (--CLEAN_CODING is in statmodes.rs.)

use super::params::Params;
use super::{compare, gchmm, genemodes, gibbs, mkmod, orf};

/// `--parse`: coding and non-coding sequences of a genome + gene list.
pub fn run_parse(p: &Params) -> bool {
    genemodes::parse(p)
}
/// `--intergenic`: sequence outside a list of regions.
pub fn run_intergenic(p: &Params) -> bool {
    genemodes::intergenic(p)
}
/// `--mkchain`: coding / non-coding chain tables from a genome + gene list.
pub fn run_mkchain(p: &Params) -> bool {
    genemodes::mkchain(p)
}
/// `--mkmod`: GeneMarkS-T model.
pub fn run_mkmod(p: &Params) -> bool {
    mkmod::mkmod(p)
}
/// `--mkmod_euk`: GeneMark-ES coding + non-coding Markov table.
pub fn run_mkmod_euk(p: &Params) -> bool {
    mkmod::mkmod_euk(p)
}
/// `--mkmod_non`: GeneMark-ES non-coding (spacer) Markov table.
pub fn run_mkmod_non(p: &Params) -> bool {
    mkmod::mkmod_non(p)
}
/// `--compare`: agreement of two gene lists.
pub fn run_compare(p: &Params) -> bool {
    compare::run(p)
}
/// `--pre_start`: sequences upstream of the start codons (Gibbs3 input).
pub fn run_pre_start(p: &Params) -> bool {
    genemodes::pre_start(p)
}
/// `--gibbs`: RBS block of a model from Gibbs sampler output.
pub fn run_gibbs(p: &Params) -> bool {
    gibbs::run(p)
}
/// `--combine` (0x442c58 on the model object 0x4062a2): reads two model
/// files and writes a combined one through the model writer 0x445560.
pub fn run_combine(p: &Params) -> bool {
    mkmod::combine(p)
}
/// `--STAT_NON_VS_ALL`: non-coding fraction and GC content.
pub fn run_stat_non_vs_all(p: &Params) -> bool {
    genemodes::stat_non_vs_all(p)
}
/// `--ORF_PER_SEQ`: ORF density of the intergenic sequence.
pub fn run_orf_per_seq(p: &Params) -> bool {
    orf::orf_per_seq(p)
}
/// `--orf`: ORF list of a sequence.
pub fn run_orf(p: &Params) -> bool {
    orf::orf(p)
}
/// `--GENE_SET_STAT`: neighbouring gene pair statistics.
pub fn run_gene_set_stat(p: &Params) -> bool {
    genemodes::gene_set_stat(p)
}
/// `--gc_hmm`: low-GC segments by two-state Viterbi decoding.
pub fn run_gc_hmm(p: &Params) -> bool {
    gchmm::run(p)
}
