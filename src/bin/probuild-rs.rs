//! probuild-rs: Rust re-implementation of the GeneMark utility probuild 2.52.
//!
//! Thin binary wrapper: the whole program lives in `src/probuild/`
//! (start reading at `src/probuild/mod.rs`).  The crate has no library
//! target, so each binary includes its module tree by path (as `gmhmmp-rs`
//! does with `src/gmhmmp/`).

#[path = "../probuild/mod.rs"]
mod probuild;

fn main() {
    probuild::main();
}
