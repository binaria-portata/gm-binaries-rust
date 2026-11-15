#!/usr/bin/env python3
"""Randomised differential test: gmhmm-rs vs the original gmhmme3.

usage: fuzz.py --ref GMHMME3 --rs GMHMM_RS --genome A.fa [--genome B.fa]
               --models DIR_OR_FILES... --n 200 --jobs 4 --work DIR [--seed 1]

Each case cuts a random segment (200 bp .. 400 kb) from a random genome,
optionally injects N runs, IUPAC codes, case changes, CRLF line ends,
then runs both programs with a random model, format and options (-d with
random evidence/mask lines, -k, -w, -s, -z, -p/-n) and compares output
file, stdout and exit code (date line masked). Failing cases are kept in
WORK/fail_<i>/ with a run.sh to reproduce.
"""
import argparse
import os
import random
import re
import shutil
import subprocess
import sys
from concurrent.futures import ThreadPoolExecutor

DATE = re.compile(rb"^(#+ ?)?(Mon|Tue|Wed|Thu|Fri|Sat|Sun) [A-Z][a-z]{2} [ 0-9][0-9] [0-9:]{8} [0-9]{4}$", re.M)


def read_fasta(path):
    seqs, name, buf = {}, None, []
    with open(path) as f:
        for line in f:
            if line.startswith(">"):
                if name:
                    seqs[name] = "".join(buf)
                name, buf = line[1:].split()[0], []
            else:
                buf.append(line.strip())
    if name:
        seqs[name] = "".join(buf)
    return seqs


def make_case(rng, genomes, models, i, work):
    d = os.path.join(work, f"case_{i}")
    os.makedirs(d, exist_ok=True)
    g = rng.choice(genomes)
    name = rng.choice(list(g))
    s = g[name]
    ln = min(len(s), int(10 ** rng.uniform(2.3, 5.6)))
    a = rng.randrange(0, len(s) - ln + 1)
    seq = list(s[a:a + ln])
    if rng.random() < 0.3:
        for _ in range(rng.randint(1, 5)):
            p = rng.randrange(len(seq)); q = min(len(seq), p + rng.randint(1, 2000))
            for k in range(p, q):
                seq[k] = "N"
    if rng.random() < 0.2:
        for _ in range(rng.randint(1, 20)):
            seq[rng.randrange(len(seq))] = rng.choice("RYKMSWBDHVNXn")
    if rng.random() < 0.3:
        for _ in range(rng.randint(1, 10)):
            p = rng.randrange(len(seq)); q = min(len(seq), p + rng.randint(10, 5000))
            for k in range(p, q):
                seq[k] = seq[k].lower()
    seq = "".join(seq)
    sid = rng.choice(["ctg", "X1", "chr_2"])
    defl = rng.choice([f">{sid}", f">{sid} desc", f">{sid}\tcut {a + 1} more", f"> {sid}"])
    eol = "\r\n" if rng.random() < 0.05 else "\n"
    width = rng.choice([60, 80, 100000])
    with open(os.path.join(d, "s.fa"), "w", newline="") as f:
        f.write(defl + eol)
        for k in range(0, len(seq), width):
            f.write(seq[k:k + width] + eol)
    model = rng.choice(models)
    shutil.copy(model, os.path.join(d, "m.mod"))
    args = ["-m", "m.mod"]
    fmt = rng.choice(["lst", "lst", "gtf", "gtf", "gff3", "tr"])
    if fmt != "lst":
        args += ["-f", fmt]
    if rng.random() < 0.3:
        L = len(seq)
        lines = []
        feats_pm = ["Single", "Initial", "Internal", "Terminal", "Intron", "intron", "Part_exon",
                    "Part_don", "Part_acc", "Part_ini", "Part_ter", "Donor", "Acceptor",
                    "Start_codon", "Stop_codon", "UTR_Exon", "Non_Coding", "CDS"]
        for _ in range(rng.randint(1, 40)):
            st = rng.randrange(1, L + 1); en = min(L, st + rng.randint(0, 3000))
            if rng.random() < 0.4:
                lines.append(f"{sid}\tsrc\tmask\t{st}\t{en}\t0\t.\t.\tsoftmask")
            else:
                f_ = rng.choice(feats_pm + ["Intergenic", "Part_intergenic", "Part_intron"])
                strand = "." if f_ in ("Intergenic", "Part_intergenic") else rng.choice("+-")
                lines.append(f"{sid}\tsrc\t{f_}\t{st}\t{en}\t.\t{strand}\t.\tx")
        with open(os.path.join(d, "e.gff"), "w") as f:
            f.write("\n".join(lines) + "\n")
        args += ["-d", "e.gff"]
        if rng.random() < 0.7:
            args += ["-k", rng.choice(["-0.03", "-0.5", "0.1", "-2"])]
    elif rng.random() < 0.1:
        args += ["-k", "-0.03"]
    if rng.random() < 0.15:
        args += ["-w", str(rng.choice([100, 300, 600, 1500]))]
    if rng.random() < 0.1:
        args += ["-s", "TAG1"]
    if rng.random() < 0.1:
        args += ["-z"]
    if fmt == "lst" and rng.random() < 0.3:
        args += ["-p", "-n"] if rng.random() < 0.5 else ["-p"]
    if rng.random() < 0.5:
        args = ["-r"] + args
    args += ["s.fa"]
    return d, args


def run(prog, d, args, tag):
    out = os.path.join(d, f"{tag}.out")
    if os.path.exists(out):
        os.remove(out)
    p = subprocess.run([prog] + args + ["-o", f"{tag}.out"], cwd=d, capture_output=True, timeout=3600)
    data = open(out, "rb").read() if os.path.exists(out) else None
    if data is not None:
        data = DATE.sub(b"<DATE>", data)
    return p.returncode, p.stdout, data


def one(i, rng_seed, opts, genomes, models):
    rng = random.Random(rng_seed)
    d, args = make_case(rng, genomes, models, i, opts.work)
    a = run(opts.ref, d, args, "ref")
    b = run(opts.rs, d, args, "rs")
    ok = a == b
    if a[0] < 0:
        # the original did not finish (e.g. evidence marks written past the
        # sequence end); nothing to compare
        shutil.rmtree(d)
        return i, None, " ".join(args), a[0]
    if ok:
        shutil.rmtree(d)
    else:
        with open(os.path.join(d, "run.sh"), "w") as f:
            f.write(" ".join(["$PROG"] + args) + "\n")
        os.rename(d, os.path.join(opts.work, f"fail_{i}"))
    return i, ok, " ".join(args), a[0]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--ref", required=True)
    ap.add_argument("--rs", required=True)
    ap.add_argument("--genome", action="append", required=True)
    ap.add_argument("--models", nargs="+", required=True)
    ap.add_argument("--n", type=int, default=100)
    ap.add_argument("--jobs", type=int, default=4)
    ap.add_argument("--work", required=True)
    ap.add_argument("--seed", type=int, default=1)
    o = ap.parse_args()
    o.ref = os.path.abspath(o.ref); o.rs = os.path.abspath(o.rs)
    os.makedirs(o.work, exist_ok=True)
    genomes = [read_fasta(g) for g in o.genome]
    models = []
    for m in o.models:
        if os.path.isdir(m):
            models += [os.path.join(m, x) for x in sorted(os.listdir(m)) if x.endswith(".mod")]
        else:
            models.append(m)
    master = random.Random(o.seed)
    seeds = [master.randrange(1 << 30) for _ in range(o.n)]
    nfail = n_unfinished = 0
    with ThreadPoolExecutor(o.jobs) as ex:
        for i, ok, args, rc in ex.map(lambda t: one(t[0], t[1], o, genomes, models), enumerate(seeds)):
            if ok is None:
                n_unfinished += 1
                tag = "REFFAIL"
            elif not ok:
                nfail += 1
                tag = "FAIL"
            else:
                tag = "ok"
            print(f"{tag:8} {i} rc={rc} {args}", flush=True)
    print(f"cases {o.n}, failed {nfail}, original not finished {n_unfinished}")
    sys.exit(1 if nfail else 0)


if __name__ == "__main__":
    main()
