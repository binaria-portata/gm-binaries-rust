#!/usr/bin/env python3
"""Run hand-written probuild calls through the original and the port.

usage: run_cases.py --ref PROBUILD --rs PROBUILD_RS [--cases FILE]
                    [--filter REGEX] [--keep DIR] [--jobs N]

cases.tsv: one case per line, `name<TAB>arguments` (shell syntax, parsed
with shlex, the text <TAB> becomes a tab character); '#' lines and blank lines are ignored.  Each case runs in two
fresh directories holding a copy of tests/probuild/data/; stdout, stderr,
exit code and every file in the directory afterwards are compared byte for
byte (no stored expectations: the original is the oracle).  Both programs
run with a 1 GB address-space limit.
"""
import argparse
import os
import re
import resource
import shlex
import shutil
import subprocess
import sys
import tempfile
from concurrent.futures import ThreadPoolExecutor

HERE = os.path.dirname(os.path.abspath(__file__))


def snapshot(base):
    out = {}
    for root, _, files in os.walk(base):
        for f in files:
            p = os.path.join(root, f)
            if os.path.islink(p):
                continue
            out[os.path.relpath(p, base)] = open(p, "rb").read()
    return out


def limit_memory():
    # --count_substring_atcg loops forever on lower-case letters while its
    # segment vector grows; with 1 GB of address space the original ends
    # with std::bad_alloc after a few seconds instead of eating the machine.
    resource.setrlimit(resource.RLIMIT_AS, (1 << 30, 1 << 30))


def run(prog, base, args):
    try:
        p = subprocess.run([prog] + args, cwd=base, capture_output=True, timeout=600,
                           preexec_fn=limit_memory)
        return p.returncode, p.stdout, p.stderr, snapshot(base)
    except subprocess.TimeoutExpired:
        return "timeout", b"", b"", snapshot(base)


def one(case, o):
    name, args = case
    with tempfile.TemporaryDirectory() as t:
        res = []
        for tag, prog in (("ref", o.ref), ("rs", o.rs)):
            d = os.path.join(t, tag)
            shutil.copytree(os.path.join(HERE, "data"), d)
            res.append(run(prog, d, args))
        ra, rb = res
        diffs = []
        for k, what in ((0, "rc"), (1, "stdout"), (2, "stderr")):
            if ra[k] != rb[k]:
                diffs.append(what)
        for f in sorted(set(ra[3]) | set(rb[3])):
            if ra[3].get(f) != rb[3].get(f):
                diffs.append(f)
        if diffs and o.keep:
            k = os.path.join(o.keep, name)
            shutil.rmtree(k, ignore_errors=True)
            shutil.copytree(t, k)
            for tag, r in (("ref", ra), ("rs", rb)):
                for ext, v in (("rc", str(r[0]).encode()), ("stdout", r[1]), ("stderr", r[2])):
                    open(os.path.join(k, f"{tag}.{ext}"), "wb").write(v)
        return name, " ".join(shlex.quote(a) for a in args), diffs


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--ref", required=True)
    ap.add_argument("--rs", required=True)
    ap.add_argument("--cases", default=os.path.join(HERE, "cases.tsv"))
    ap.add_argument("--filter", default=".")
    ap.add_argument("--keep")
    ap.add_argument("--jobs", type=int, default=8)
    o = ap.parse_args()
    o.ref = os.path.abspath(o.ref)
    o.rs = os.path.abspath(o.rs)
    cases = []
    for line in open(o.cases):
        line = line.rstrip("\n")
        if not line.strip() or line.startswith("#"):
            continue
        name, _, args = line.partition("\t")
        if re.search(o.filter, name + " " + args):
            cases.append((name, [a.replace("<TAB>", "\t") for a in shlex.split(args)]))
    names = [c[0] for c in cases]
    dup = {n for n in names if names.count(n) > 1}
    if dup:
        sys.exit(f"duplicate case names: {sorted(dup)}")
    nbad = 0
    with ThreadPoolExecutor(o.jobs) as ex:
        for name, args, diffs in ex.map(lambda c: one(c, o), cases):
            nbad += bool(diffs)
            print(("ok   " if not diffs else "FAIL ") + name + " " + args[:140]
                  + ("" if not diffs else "  diff: " + ",".join(diffs[:6])), flush=True)
    print(f"cases {len(cases)}, failed {nbad}")
    sys.exit(1 if nbad else 0)


if __name__ == "__main__":
    main()
