//! gmhmmp-rs: Rust re-implementation of GeneMark.hmm PROKARYOTIC 3.20.
//! Drop-in replacement for `gmhmmp` (same command line, output, messages
//! and exit codes). The program lives in `src/gmhmmp/` (see its `mod.rs`
//! for the pipeline); this file only provides the binary entry point.

#[path = "../gmhmmp/mod.rs"]
mod gmhmmp;

fn main() {
    gmhmmp::main();
}
