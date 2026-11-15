//! Developer dumps, not part of the original program. They write internal
//! data in the layouts of the gdb dumps made from the original binary, so
//! that the port can be compared value by value (bit patterns, not printed
//! decimals):
//! * [`model_json`] (`gmhmm-rs --dump-model-json MODEL`): the transformed
//!   model, JSON layout of the model dumper (re-notes/model.md section 6),
//!   compared by `tests/cmp_model_dump.py`;
//! * [`input_bins`] (`--dump-input ...`): codes, k-mer index, site map and
//!   evidence marks (re-notes/input.md);
//! * [`bp_bins`] (`GMHMM_RS_BPDUMP=<dir>`): branch point data, intron and
//!   site tracks and the Viterbi cells (re-notes/bp.md).
//!
//! Doubles are written as 16-digit hex bit patterns (JSON) or raw
//! little-endian bytes.

use crate::model::{LenDist, Model, Site};
use std::fmt::Write;

/// f64 as a quoted 16-digit hex bit pattern.
fn hx(x: f64) -> String {
    format!("\"{:016x}\"", x.to_bits())
}

/// JSON array of bit patterns.
fn arr(v: &[f64]) -> String {
    let parts: Vec<String> = v.iter().map(|&x| hx(x)).collect();
    format!("[{}]", parts.join(","))
}

/// JSON array of arrays of bit patterns (lower-order Markov tables).
fn arr2(v: &[Vec<f64>]) -> String {
    let parts: Vec<String> = v.iter().map(|x| arr(x)).collect();
    format!("[{}]", parts.join(","))
}

/// One LenDist in the dumper's layout: the struct ints (+0x00..+0x0c,
/// +0x20..+0x2c), the base pointer offsets and both tables.
fn len_dist_json(d: &LenDist) -> String {
    format!(
        "{{\"ints_0_4_8_c\":[{},{},1,{}],\"ints_20_24_28_2c\":[{},{},{},{}],\"tab_base_minus_tab\":{},\"cum_base_minus_cum\":{},\"tab_vals\":{},\"cum_vals\":{}}}",
        d.max, d.min, d.max - 1, d.lo(), d.hi(), d.tail_lo(), d.tail_hi(), -d.p4, -d.p5, arr(&d.tab), arr(&d.tail)
    )
}

/// One site model in the dumper's layout.
fn site_json(s: &Site) -> String {
    format!(
        "{{\"margin\":{},\"width\":{},\"order\":{},\"margin_p\":{},\"width_p\":{},\"n\":{},\"mat_vals\":{},\"rc_vals\":{}}}",
        s.margin, s.width, s.order, s.margin_p, s.width_p, s.n, arr(&s.mat), arr(&s.rc)
    )
}

/// The whole transformed model as one JSON object. Markov tables are keyed
/// by their offset in the original model object (e.g. "0x948").
pub fn model_json(m: &Model) -> String {
    let mut o = String::new();
    write!(o, "{{\"organism\":\"{}\",\"order\":{},\"len\":{{", m.organism, m.order).unwrap();
    let l = [
        ("SINGLE", &m.single),
        ("EXON", &m.internal),
        ("INITIAL", &m.initial),
        ("TERMINAL", &m.terminal),
        ("INTRON", &m.intron),
        ("INTERGENIC", &m.intergenic),
    ];
    let mut l = l.to_vec();
    if let Some(b) = &m.bp {
        l.push(("DON_BP", &b.don_bp));
        l.push(("BP_ACC", &b.bp_acc));
    }
    let v: Vec<String> = l.iter().map(|(n, d)| format!("\"{n}\":{}", len_dist_json(d))).collect();
    o += &v.join(",");
    o += "},\"scalars\":{";
    let sc = [
        (0x990, m.log_to_single),
        (0x998, m.log_to_multi),
        (0x9a0, m.log_to_internal),
        (0x9a8, m.log_to_terminal),
        (0x9b0, m.log_initial[0]),
        (0x9b8, m.log_initial[1]),
        (0x9c0, m.log_initial[2]),
        (0x9c8, m.log_initial[3]),
        (0x9d0, m.log_initial[4]),
        (0x9e0, m.log_gt_donor),
        (0x9e8, m.log_gc_donor),
    ];
    let v: Vec<String> = sc.iter().map(|(k, x)| format!("\"0x{k:x}\":{}", hx(*x))).collect();
    o += &v.join(",");
    write!(
        o,
        "}},\"phase\":{{\"INITIAL\":{},\"TERMINAL\":{},\"INTERNAL\":{}}},\"markov\":{{",
        arr(&m.phase_initial),
        arr(&m.phase_terminal),
        arr(&m.phase_internal)
    )
    .unwrap();
    let c = &m.coding;
    let mut mk = vec![
        ("0x948".to_string(), arr(&m.nc3.fwd)),
        ("0x950".into(), arr(&m.nc3.fwd_nf)),
        ("0x958".into(), arr(&m.nc3.rc)),
        ("0x960".into(), arr(&m.nc3.rc_nf)),
        ("0x968".into(), arr2(&m.nc3.low)),
        ("0x930".into(), arr(&m.nc4.fwd)),
        ("0x938".into(), arr(&m.nc4.fwd_nf)),
        ("0x940".into(), arr2(&m.nc4.low)),
    ];
    for k in 0..3 {
        let b = 0x7e8 + 0x10 * k;
        mk.push((format!("0x{:x}", b), arr(&c.fwd[k])));
        mk.push((format!("0x{:x}", b + 8), arr(&c.fwd_nf[k])));
        mk.push((format!("0x{:x}", 0x818 + 0x10 * k), arr(&c.rc[k])));
        mk.push((format!("0x{:x}", 0x820 + 0x10 * k), arr(&c.rc_nf[k])));
        mk.push((format!("0x{:x}", 0x848 + 8 * k), arr2(&c.low[k])));
    }
    let v: Vec<String> = mk.iter().map(|(k, x)| format!("\"{k}\":{x}")).collect();
    o += &v.join(",");
    o += "},\"sites\":{";
    let mut s = vec![
        ("INI", &m.start),
        ("TERM_TAA", &m.stop[0]),
        ("TERM_TAG", &m.stop[1]),
        ("TERM_TGA", &m.stop[2]),
        ("DONOR_0", &m.donor[0]),
        ("DONOR_1", &m.donor[1]),
        ("DONOR_2", &m.donor[2]),
        ("ACCEPTOR_0", &m.acceptor[0]),
        ("ACCEPTOR_1", &m.acceptor[1]),
        ("ACCEPTOR_2", &m.acceptor[2]),
    ];
    if let Some(g) = &m.donor_gc {
        s.push(("DONOR_GC_0", &g[0]));
        s.push(("DONOR_GC_1", &g[1]));
        s.push(("DONOR_GC_2", &g[2]));
    }
    if let Some(b) = &m.bp {
        s.push(("BRANCH", &b.branch));
        s.push(("ACC_BP_0", &b.acc_bp[0]));
        s.push(("ACC_BP_1", &b.acc_bp[1]));
        s.push(("ACC_BP_2", &b.acc_bp[2]));
    }
    let v: Vec<String> = s.iter().map(|(n, x)| format!("\"{n}\":{}", site_json(x))).collect();
    o += &v.join(",");
    o += "}";
    if let Some(b) = &m.bp {
        write!(
            o,
            ",\"bp_scalars\":{{\"0x6a8\":{},\"0x6b0\":{}}},\"bp\":{{\"limit_6b8\":{},\"spacer_order_6f0\":{},\"0x6c0\":{},\"0x6c8\":{},\"0x6d0\":{},\"0x6d8\":{},\"0x6e0\":{}}}",
            hx(b.log_to_bp),
            hx(b.log_around_bp),
            b.limit_duration,
            b.spacer_order,
            arr(&b.spacer.fwd),
            arr(&b.spacer.fwd_nf),
            arr(&b.spacer.rc),
            arr(&b.spacer.rc_nf),
            arr2(&b.spacer.low)
        )
        .unwrap();
    }
    o += "}";
    o
}

/// Write codes, k-mer index, final site map (and evidence marks) in the
/// binary layouts of the debugger dumps of the original (re-notes/input.md).
/// args: seq.fa order outprefix taa,tag,tga,gc [evidence.gff]
/// (the flags as 0/1, e.g. `1,1,1,0`). Map records are written as 16-byte
/// records with zeroed padding.
pub fn input_bins(args: &[String]) {
    use crate::{evidence, sequence::Sequence, sitemap};
    let mut s = Sequence::read(&args[0], "").expect("read");
    let order: i32 = args[1].parse().unwrap();
    let pre = &args[2];
    let fl: Vec<bool> = args.get(3).map(|x| x.split(',').map(|v| v == "1").collect()).unwrap_or(vec![true; 4]);
    s.build_kmer(order).expect("kmer");
    if let Some(g) = args.get(4) {
        let mut out = String::new();
        s.marks = Some(evidence::read(g, &s.seqid, s.len(), &mut out).expect("ev"));
        print!("{out}");
        let b: Vec<u8> = s.marks.as_ref().unwrap().iter().flat_map(|v| v.to_le_bytes()).collect();
        std::fs::write(format!("{pre}_marks.bin"), b).unwrap();
    }
    std::fs::write(format!("{pre}_codes.bin"), &s.codes).unwrap();
    let b: Vec<u8> = s.kmer.iter().flat_map(|v| v.to_le_bytes()).collect();
    std::fs::write(format!("{pre}_kmer.bin"), b).unwrap();
    let f = sitemap::Flags { taa: fl[0], tag: fl[1], tga: fl[2], gc_donor: fl[3] };
    let m = sitemap::build(&s.codes, s.marks.as_deref(), &f);
    let mut b = Vec::with_capacity(m.len() * 16);
    for r in &m {
        b.extend(r.typ.to_le_bytes());
        b.extend(r.pos.to_le_bytes());
        b.push(r.frame);
        b.push(r.sub);
        b.extend([0u8, 0]);
        b.extend(r.mark.to_le_bytes());
    }
    std::fs::write(format!("{pre}_C.bin"), b).unwrap();
}

/// Write the BP object, the intron and site tracks and the Viterbi cells in
/// the raw layouts of the debugger dump of the original (re-notes/bp.md), into `dir`.
/// Enabled with the environment variable GMHMM_RS_BPDUMP=<dir>.
pub fn bp_bins(dir: &str, d: &crate::bp::BpData, t: &crate::tracks::Tracks, bp: &[i32], score: &[f64]) {
    let f64s = |v: &[f64]| -> Vec<u8> { v.iter().flat_map(|x| x.to_le_bytes()).collect() };
    let _ = std::fs::create_dir_all(dir);
    let w = |n: &str, b: Vec<u8>| std::fs::write(format!("{dir}/{n}"), b).unwrap();
    w("bp_arr.bin", f64s(&d.arr));
    w("bp_fwd.bin", f64s(&d.fwd));
    w("bp_rev.bin", f64s(&d.rev));
    w("bp_csum_f.bin", f64s(&d.csum_f));
    w("bp_csum_r.bin", f64s(&d.csum_r));
    w("trk_in.bin", f64s(&t.intron));
    w("trk_site.bin", f64s(&t.site));
    let mut c = Vec::with_capacity(score.len() * 24);
    for (b, s) in bp.iter().zip(score) {
        c.extend(b.to_le_bytes());
        c.extend([0u8; 4]);
        c.extend(s.to_le_bytes());
        c.extend(crate::model::NEG.to_le_bytes());
    }
    w("cells.bin", c);
    std::fs::write(format!("{dir}/ints.txt"), format!("{} {} {} {}\n", d.arr.len(), d.l, d.acc_max, d.acc_min1)).unwrap();
}
