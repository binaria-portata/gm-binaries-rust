#!/usr/bin/env python3
"""Randomised differential test: gmhmmp-rs vs the original gmhmmp.

usage: fuzz_gmhmmp.py --ref GMHMMP --rs GMHMMP_RS --genome A.fa [--genome B.fa]
                      --models M1.mod M2.mod ... --n 200 --jobs 4 --work DIR [--seed 1]

Each case writes a multi-record FASTA of random genome segments (10 bp ..
20 kb, some reverse-complemented, some with N runs, IUPAC codes, lower case,
CRLF, [gc=..]/[shape=..] defline tags), picks a model and random options
(-p, -r, -f, -k, -K, -b, -s, -i, -x, -h, -n, -c, -A, -D) and compares the
output file, side files, stdout and exit code (date lines masked). Cases where
the original does not finish are counted separately. Failures are kept in
WORK/fail_<i>/ with run.sh.
"""
import argparse
import os
import random
import re
import shutil
import subprocess
import sys
from concurrent.futures import ThreadPoolExecutor

DATE = re.compile(rb"^(#+ ?|Date: |##date: )?(Mon|Tue|Wed|Thu|Fri|Sat|Sun) [A-Z][a-z]{2} [ 0-9][0-9] [0-9:]{8} [0-9]{4}$", re.M)
COMP = str.maketrans("ACGTacgtNn", "TGCAtgcaNn")


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
    names = list(g)
    nrec = rng.randint(1, 12)
    eol = "\r\n" if rng.random() < 0.05 else "\n"
    parse_c = rng.random() < 0.15
    with open(os.path.join(d, "s.fa"), "w", newline="") as f:
        if rng.random() < 0.05:
            f.write("ACGTACGTACGTACGTACGTACGTACGT" + eol)
        for r in range(nrec):
            s = g[rng.choice(names)]
            ln = min(len(s), int(10 ** rng.uniform(1.0, 4.3)))
            a = rng.randrange(0, len(s) - ln + 1)
            seq = s[a:a + ln]
            if rng.random() < 0.4:
                seq = seq.translate(COMP)[::-1]
            seq = list(seq)
            if rng.random() < 0.25 and seq:
                for _ in range(rng.randint(1, 3)):
                    p = rng.randrange(len(seq)); q = min(len(seq), p + rng.randint(1, 120))
                    for k in range(p, q):
                        seq[k] = "N"
            if rng.random() < 0.15 and seq:
                for _ in range(rng.randint(1, 8)):
                    seq[rng.randrange(len(seq))] = rng.choice("RYKMSWBDHVNXnry")
            if rng.random() < 0.2:
                seq = [c.lower() for c in seq]
            seq = "".join(seq)
            defl = f">r{r} seg {a + 1}"  # id = r<k> (evidence key)
            if parse_c and rng.random() < 0.7:
                defl += rng.choice([" [gc=%d]" % rng.randint(0, 100), " [shape=linear]", " [shape=circular]", " [shape=xx]", " [gc=12.7][shape=partial]"])
            width = rng.choice([60, 70, 100000])
            f.write(defl + eol)
            for k in range(0, len(seq), width):
                f.write(seq[k:k + width] + eol)
    if rng.random() < 0.15:
        # random evidence GFF for the records written above
        lines = []
        for _ in range(rng.randint(1, 40)):
            r = rng.randrange(nrec)
            st = rng.randint(1, 3000); en = st + rng.randint(0, 900)
            t = rng.choice("MNCFRHTEmncfx")
            strand = rng.choice("+-.")
            ph = rng.choice("0.")
            sc = rng.choice([".", "5", "-3.5", "0"])
            lines.append(f"r{r}\tsrc\t{t}\t{st}\t{en}\t{sc}\t{strand}\t{ph}\tID=e{rng.randint(0, 99)};")
        with open(os.path.join(d, "e.gff"), "w") as f:
            f.write("\n".join(lines) + "\n")
        ev_args = ["-e", "e.gff"]
        if rng.random() < 0.3:
            ev_args += ["-E", rng.choice(["0", "5", "250"])]
        if rng.random() < 0.3:
            ev_args += ["-v"]
    else:
        ev_args = []
    model = rng.choice(models)
    shutil.copy(model, os.path.join(d, "m.mod"))
    args = ["-m", "m.mod"]
    if rng.random() < 0.5:
        args += ["-p", "0"]
    if rng.random() < 0.5:
        args += ["-r"]
    if rng.random() < 0.5:
        args += ["-f", "G"]
    side = None
    roll = rng.random()
    if roll < 0.15:
        args += ["-K", "side.b"]; side = "side.b"
    elif roll < 0.3:
        args += ["-A", "side.b"]; side = "side.b"
    elif roll < 0.45:
        args += ["-D", "side.b"]; side = "side.b"
    if rng.random() < 0.3:
        args += ["-k"]
    if rng.random() < 0.15:
        args += ["-b"]
    if rng.random() < 0.15:
        args += ["-s", rng.choice(["d", "r"])]
    if rng.random() < 0.1:
        args += ["-i", rng.choice(["0.2", "0.9", "0.5"])]
    if rng.random() < 0.1:
        args += ["-x", rng.choice(["0.1", "0.5"])]
    if rng.random() < 0.1:
        args += ["-h", rng.choice(["300", "600"])]
    if rng.random() < 0.15:
        args += ["-n"]
    if parse_c:
        args += ["-c"]
    args += ev_args
    args += ["s.fa"]
    return d, args, side


def run(prog, d, args, tag, side):
    for f in (f"{tag}.out", "side.b"):
        p = os.path.join(d, f)
        if os.path.exists(p):
            os.remove(p)
    p = subprocess.run([prog] + args + ["-o", f"{tag}.out"], cwd=d, capture_output=True, timeout=3600)
    res = [p.returncode, p.stdout, p.stderr]
    for f in (f"{tag}.out", "side.b"):
        q = os.path.join(d, f)
        res.append(DATE.sub(b"<DATE>", open(q, "rb").read()) if os.path.exists(q) else None)
    return res


def one(i, seed, o, genomes, models):
    rng = random.Random(seed)
    d, args, side = make_case(rng, genomes, models, i, o.work)
    a = run(o.ref, d, args, "ref", side)
    b = run(o.rs, d, args, "rs", side)
    if a[0] < 0:
        shutil.rmtree(d)
        return i, None, " ".join(args), a[0]
    ok = a == b
    if ok:
        shutil.rmtree(d)
    else:
        with open(os.path.join(d, "run.sh"), "w") as f:
            f.write(" ".join(["$PROG"] + args) + "\n")
        os.rename(d, os.path.join(o.work, f"fail_{i}"))
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
    models = [os.path.abspath(m) for m in o.models]
    master = random.Random(o.seed)
    seeds = [master.randrange(1 << 30) for _ in range(o.n)]
    nfail = n_unfinished = 0
    with ThreadPoolExecutor(o.jobs) as ex:
        for i, ok, args, rc in ex.map(lambda t: one(t[0], t[1], o, genomes, models), enumerate(seeds)):
            tag = "REFFAIL" if ok is None else ("ok" if ok else "FAIL")
            if ok is None:
                n_unfinished += 1
            elif not ok:
                nfail += 1
            print(f"{tag:8} {i} rc={rc} {args}", flush=True)
    print(f"cases {o.n}, failed {nfail}, original not finished {n_unfinished}")
    sys.exit(1 if nfail else 0)


if __name__ == "__main__":
    main()
