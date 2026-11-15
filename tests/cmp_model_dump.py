#!/usr/bin/env python3
"""Compare `gmhmm-rs --dump-model-json` with a gdb dump of gmhmme3's model.

usage: cmp_model_dump.py rust.json gdb.json
The gdb dump comes from a script in the analysis notes (not distributed).
Prints mismatch counts per table; exit 1 on any mismatch.
"""
import json
import sys

r = json.load(open(sys.argv[1]))
g = json.load(open(sys.argv[2]))
bad = 0
tot = 0


def cmp(label, a, b):
    global bad, tot
    if isinstance(a, list) and a and isinstance(a[0], list):
        if len(a) != len(b):
            print("LEN", label, len(a), len(b)); bad += 1; return
        for i, (x, y) in enumerate(zip(a, b)):
            cmp(f"{label}[{i}]", x, y)
        return
    tot += len(a)
    if len(a) != len(b):
        print("LEN", label, len(a), len(b)); bad += 1; return
    n = sum(1 for x, y in zip(a, b) if x != y)
    if n:
        i = next(i for i, (x, y) in enumerate(zip(a, b)) if x != y)
        print(f"MISMATCH {label}: {n}/{len(a)} first {i}: {a[i]} vs {b[i]}")
        bad += n


for name, d in r["len"].items():
    e = g["len"][name]
    if d["ints_0_4_8_c"] != e["ints_0_4_8_c"] or d["ints_20_24_28_2c"] != e["ints_20_24_28_2c"]:
        print("INTS", name, d["ints_0_4_8_c"], e["ints_0_4_8_c"], d["ints_20_24_28_2c"], e["ints_20_24_28_2c"]); bad += 1
    cmp(f"len {name} tab", d["tab_vals"], e["tab_vals"])
    cmp(f"len {name} cum", d["cum_vals"], e["cum_vals"])
for k, v in r["scalars"].items():
    cmp(f"scalar {k}", [v], [g["scalars"][k]])
for k, v in r["phase"].items():
    cmp(f"phase {k}", v, g["phase"][k])
for k, v in r["markov"].items():
    cmp(f"markov {k}", v, g["markov"][k])
for k, v in r["sites"].items():
    e = g["sites"][k]
    gi = e["ints"]  # margin, width, order, ..., +0x20 margin_p, +0x24 width_p, +0x4c n
    want = [gi[0], gi[1], gi[2], gi[8], gi[9], gi[0x13]]
    have = [v["margin"], v["width"], v["order"], v["margin_p"], v["width_p"], v["n"]]
    if want != have:
        print("SITE INTS", k, have, want); bad += 1
    cmp(f"site {k} mat", v["mat_vals"], e["mat_vals"])
    cmp(f"site {k} rc", v["rc_vals"], e["rc_vals"])
if "bp" in r or "bp" in g:
    if "bp" not in r or "bp" not in g or g.get("bp_flag_5d0") != 1:
        print("BP section present in only one dump"); bad += 1
    else:
        for k in ("0x6a8", "0x6b0"):
            cmp(f"bp scalar {k}", [r["bp_scalars"][k]], [g["bp_scalars"][k]])
        for k in ("limit_6b8", "spacer_order_6f0"):
            if r["bp"][k] != g["bp"][k]:
                print("BP INT", k, r["bp"][k], g["bp"][k]); bad += 1
        for k in ("0x6c0", "0x6c8", "0x6d0", "0x6d8", "0x6e0"):
            cmp(f"bp spacer {k}", r["bp"][k], g["bp"][k])
print(f"compared {tot} values, mismatches {bad}")
sys.exit(1 if bad else 0)
