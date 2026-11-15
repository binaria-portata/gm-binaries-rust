#!/usr/bin/env python3
"""Replay captured calls (probuild, gmhmme3, gmhmmp) through the original and the port.

usage: replay.py --ref PROBUILD --rs PROBUILD_RS CAPDIR [CAPDIR ...]
                 [--filter REGEX] [--keep DIR] [--jobs N]

A capture directory holds call_*/ subdirectories written by
probuild_wrapper.sh: `args0` (NUL-separated argv), `cwd`, `in/<i>` (copy of
argument i if it was a file before the run), `out/<i>`, stdout, stderr, rc.

For each call: two fresh directories get the input files at the argument
paths (absolute paths are rewritten to relative ones), both programs run
there, and stdout, stderr, exit code and every file in the directory tree
afterwards are compared byte for byte (ctime date lines masked). Differences
are kept in --keep.
"""
import argparse
import os
import re
import shutil
import subprocess
import sys
import tempfile
from concurrent.futures import ThreadPoolExecutor


def load_call(d):
    args = open(os.path.join(d, "args0"), "rb").read().split(b"\0")[:-1]
    args = [a.decode() for a in args]
    cwd = open(os.path.join(d, "cwd")).read().strip()
    ins = {}
    for f in os.listdir(os.path.join(d, "in")):
        ins[int(f)] = os.path.join(d, "in", f)
    return args, cwd, ins


def rel(path, cwd):
    if os.path.isabs(path):
        r = os.path.relpath(path, cwd)
        if r.startswith(".."):
            r = path.lstrip("/")
        return r
    return path


def setup(base, args, cwd, ins):
    new = []
    for i, a in enumerate(args):
        if i in ins:
            r = rel(a, cwd)
            p = os.path.join(base, r)
            os.makedirs(os.path.dirname(p) or base, exist_ok=True)
            shutil.copy(ins[i], p)
            new.append(r)
        elif os.path.isabs(a):
            r = rel(a, cwd)
            os.makedirs(os.path.join(base, os.path.dirname(r)) or base, exist_ok=True)
            new.append(r)
        else:
            d = os.path.dirname(a)
            if d and not a.startswith("-"):
                os.makedirs(os.path.join(base, d), exist_ok=True)
            new.append(a)
    return new


DATE = re.compile(rb"^(#+ ?|Date: |##date: )?(Mon|Tue|Wed|Thu|Fri|Sat|Sun) [A-Z][a-z]{2} [ 0-9][0-9] [0-9:]{8} [0-9]{4}$", re.M)


def snapshot(base):
    out = {}
    for root, _, files in os.walk(base):
        for f in files:
            p = os.path.join(root, f)
            out[os.path.relpath(p, base)] = DATE.sub(b"<DATE>", open(p, "rb").read())
    return out


def run(prog, base, args):
    import time
    t0 = time.perf_counter()
    p = subprocess.run([prog] + args, cwd=base, capture_output=True, timeout=3600)
    dt = time.perf_counter() - t0
    return p.returncode, p.stdout, p.stderr, snapshot(base), dt


def one(call, o):
    args, cwd, ins = load_call(call)
    with tempfile.TemporaryDirectory() as t:
        a_dir = os.path.join(t, "ref"); b_dir = os.path.join(t, "rs")
        os.makedirs(a_dir); os.makedirs(b_dir)
        new = setup(a_dir, args, cwd, ins)
        setup(b_dir, args, cwd, ins)
        ra = run(o.ref, a_dir, new)
        rb = run(o.rs, b_dir, new)
        diffs = []
        for k, name in ((0, "rc"), (1, "stdout"), (2, "stderr")):
            if ra[k] != rb[k]:
                diffs.append(name)
        for f in sorted(set(ra[3]) | set(rb[3])):
            if ra[3].get(f) != rb[3].get(f):
                diffs.append(f)
        if diffs and o.keep:
            k = os.path.join(o.keep, os.path.basename(call))
            shutil.rmtree(k, ignore_errors=True)
            shutil.copytree(t, k)
            with open(os.path.join(k, "args"), "w") as f:
                f.write(" ".join(new) + "\n")
        return call, " ".join(new), diffs, ra[4], rb[4]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--ref", required=True)
    ap.add_argument("--rs", required=True)
    ap.add_argument("--filter", default=".")
    ap.add_argument("--keep")
    ap.add_argument("--jobs", type=int, default=4)
    ap.add_argument("caps", nargs="+")
    o = ap.parse_args()
    o.ref = os.path.abspath(o.ref); o.rs = os.path.abspath(o.rs)
    calls = []
    for c in o.caps:
        for d in sorted(os.listdir(c)):
            if d.startswith("call_"):
                p = os.path.join(c, d)
                if re.search(o.filter, open(os.path.join(p, "args")).read()):
                    calls.append(p)
    nbad = 0
    tref = trs = 0.0
    with ThreadPoolExecutor(o.jobs) as ex:
        for call, args, diffs, ta, tb in ex.map(lambda c: one(c, o), calls):
            tref += ta
            trs += tb
            if diffs:
                nbad += 1
            print(("ok   " if not diffs else "FAIL ") + os.path.basename(os.path.dirname(call)) + "/" + os.path.basename(call)
                  + " " + args[:150] + ("" if not diffs else "  diff: " + ",".join(diffs[:6])), flush=True)
    print(f"calls {len(calls)}, failed {nbad}; total runtime original {tref:.1f} s, port {trs:.1f} s")
    sys.exit(1 if nbad else 0)


if __name__ == "__main__":
    main()
