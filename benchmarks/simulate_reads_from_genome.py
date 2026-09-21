#!/usr/bin/env python3
"""
simulate_reads_from_genome.py — sample long reads from a mock genome.

This is a revision of simulate_fp_read.py, the module-1 simulator.This version samples reads
from the mock genomes built by build_mock_genome.py, for the genome-scale module.

Four changes from simulate_fp_read.py
-------------------------------------
1. The original samples only the forward strand, which is not a realistic library. 

2. Depth is specified as -c/--coverage and the read count is derived: 
n = round(coverage * genome_len / read_len). 

3. HERALD accepts either. Quality is a constant Phred consistent with --error-rate,
   since the error model is uniform. Use --fasta for FASTA instead.

4. Every read records the 0-based start coordinate itwas drawn from in the mock genome, plus its strand:
       @pair1_rep1_read000042 pos=1893004 strand=- len=20000

Error model is substitutions only, matching module 1: stand-in for HiFi but not for ONT.

Examples
--------
  # Long-read arm, 50x, matching module-1 read conditions
  python simulate_reads_from_genome.py -i mock/pair1_rep1.fasta \\
      -o reads/pair1_rep1.long.fq.gz -c 50 -l 20000 -e 0.005 -s 101

"""

import argparse
import gzip
import math
import os
import random
import sys

BASES = "ACGT"
COMP = str.maketrans("ACGTNacgtn", "TGCANtgcan")


def open_out(path):
    if path.endswith(".gz"):
        return gzip.open(path, "wt")
    return open(path, "w")


def read_single_fasta(path):
    recs, name, chunks = [], None, []
    with open(path) as fh:
        for line in fh:
            if line.startswith(">"):
                if name is not None:
                    recs.append((name, "".join(chunks)))
                name = line[1:].split()[0]
                chunks = []
            else:
                chunks.append(line.strip())
    if name is not None:
        recs.append((name, "".join(chunks)))
    if not recs:
        sys.exit(f"No FASTA records in {path}")
    if len(recs) > 1:
        recs.sort(key=lambda r: -len(r[1]))
        sys.stderr.write(f"WARNING: {path} has {len(recs)} records; "
                         f"using the longest ({recs[0][0]}).\n")
    return recs[0][0], recs[0][1].upper()


def revcomp(seq):
    return seq.translate(COMP)[::-1]


def mutate(seq, rate, rng):
    """Substitute exactly round(len*rate) positions. Exact count, not a draw,
    so a given seed always yields the same library."""
    k = int(round(len(seq) * rate))
    if k <= 0:
        return seq
    s = list(seq)
    for i in rng.sample(range(len(s)), min(k, len(s))):
        b = s[i]
        if b in BASES:
            s[i] = rng.choice([x for x in BASES if x != b])
    return "".join(s)


def phred_char(error_rate):
    if error_rate <= 0:
        return "I"
    q = int(round(-10 * math.log10(error_rate)))
    return chr(min(max(q, 2), 40) + 33)


def main():
    p = argparse.ArgumentParser(
        description="Sample long reads from a mock genome, both strands, "
                    "coverage-based depth, provenance in headers")
    p.add_argument("-i", "--input", required=True, help="Mock genome FASTA")
    p.add_argument("-o", "--output", required=True,
                   help="Output .fq / .fq.gz / .fa / .fa.gz")
    p.add_argument("-c", "--coverage", type=float, default=50.0)
    p.add_argument("-l", "--read_len", type=int, default=20000)
    p.add_argument("-e", "--error_rate", type=float, default=0.005,
                   help="Substitution rate (0.005 = 0.5%%)")
    p.add_argument("-s", "--seed", type=int, required=True)
    p.add_argument("--prefix", default=None,
                   help="Read name prefix (default: output basename)")
    p.add_argument("--fasta", action="store_true", help="Write FASTA not FASTQ")
    p.add_argument("-n", "--n_reads", type=int, default=None,
                   help="Override coverage with an explicit read count")
    args = p.parse_args()

    rng = random.Random(args.seed)
    _, genome = read_single_fasta(args.input)
    L = len(genome)

    if L < args.read_len:
        sys.exit(f"Genome is {L:,} bp, shorter than read length {args.read_len:,}")

    if args.n_reads is not None:
        n = args.n_reads
    else:
        n = int(round(args.coverage * L / args.read_len))

    prefix = args.prefix or os.path.basename(args.output).split(".")[0]
    qual = phred_char(args.error_rate) * args.read_len
    width = max(6, len(str(n)))

    outdir = os.path.dirname(args.output)
    if outdir:
        os.makedirs(outdir, exist_ok=True)

    n_fwd = 0
    with open_out(args.output) as out:
        for i in range(n):
            start = rng.randint(0, L - args.read_len)
            read = genome[start:start + args.read_len]
            strand = "+" if rng.random() < 0.5 else "-"
            if strand == "-":
                read = revcomp(read)
                # NOTE: pos stays in forward-genome coordinates regardless of
                # strand, so it can be compared directly to the truth table.
            else:
                n_fwd += 1
            read = mutate(read, args.error_rate, rng)
            name = (f"{prefix}_read{i:0{width}d} pos={start} "
                    f"strand={strand} len={len(read)}")
            if args.fasta:
                out.write(f">{name}\n{read}\n")
            else:
                out.write(f"@{name}\n{read}\n+\n{qual}\n")

    achieved = n * args.read_len / L
    sys.stderr.write(
        f"{prefix}: {n:,} reads x {args.read_len:,} bp from {L:,} bp genome\n"
        f"  coverage {achieved:.2f}x | {n_fwd:,} forward / {n - n_fwd:,} reverse "
        f"| error {args.error_rate:.3%} | seed {args.seed}\n"
        f"  -> {args.output}\n")


if __name__ == "__main__":
    main()