#!/usr/bin/env python3
"""
read_census.py — build read-level truth for the genome-scale module.

Why this is needed
------------------
all_truth.tsv records GENOME-level events: where each insert sits in the chromosome. 
But recall is measured over READS, and the number of reads carrying a given insert depends on coverage, read length, and insert size. 

Without this census you cannot tell three things apart:
  - HERALD missed an event it had reads for      (algorithmic miss)
  - no read carried enough donor to be callable  (structural, not a miss)
  - no read covered the junction at all          (coverage)

The census is derived from the pos= field written into every read header by simulate_reads_from_genome.py, 
cross-referenced against mock_start/mock_end in the truth table. Both are in MOCK genome coordinates, 
so they compare directly. (ref_pos / ref_pos_vcf are acceptor-reference coordinates and are not
interchangeable with these.)

Denominator columns, weakest to strongest
----------------------------------------
  reads_any        read overlaps the insert by >=1 bp. Includes reads that
                   barely clip an edge. 
  reads_ge200      overlap >=200 bp. A reasonable "carries donor" definition.
  reads_ge_frag    overlap >= one fragment length (--frag-len). Below one
                   fragment, host sequence dominates every fragment the insert
                   touches and the read cannot be called.
  reads_flanked    insert fully inside the read with >=--flank bp of host on
                   both sides. The strictest, closest to an ideal mosaic read.

Usage
-----
  python read_census.py --truth mock/all_truth.tsv --reads reads/ \\
      --frag-len 2000 --out-events census_events.tsv --out-samples census_samples.tsv

  # exclude the calibration sample from the reported set
  python read_census.py ... --exclude pair1_phylum_rep1
"""

import argparse
import glob
import gzip
import os
import sys
from collections import defaultdict


def open_maybe_gz(path):
    return gzip.open(path, "rt") if path.endswith(".gz") else open(path)


def read_truth(path):
    """Return {sample: [event dicts]} keyed by sample name."""
    events = defaultdict(list)
    with open(path) as fh:
        hdr = fh.readline().rstrip("\n").split("\t")
        need = {"sample", "event_id", "insert_len", "mock_start", "mock_end"}
        missing = need - set(hdr)
        if missing:
            sys.exit(f"{path} is missing column(s): {', '.join(sorted(missing))}")
        for line in fh:
            if not line.strip():
                continue
            d = dict(zip(hdr, line.rstrip("\n").split("\t")))
            for k in ("insert_len", "mock_start", "mock_end"):
                d[k] = int(d[k])
            events[d["sample"]].append(d)
    for s in events:
        events[s].sort(key=lambda e: e["mock_start"])
    return events


def read_headers(path):
    """Yield (pos, length) for every read. Requires pos= in the header."""
    out = []
    warned = False
    with open_maybe_gz(path) as fh:
        first = fh.readline()
        if not first:
            return out
        fastq = first.startswith("@")
        fh.seek(0)
        if fastq:
            while True:
                h = fh.readline()
                if not h:
                    break
                seq = fh.readline().rstrip("\n")
                fh.readline()
                fh.readline()
                kv = dict(x.split("=", 1) for x in h.strip().split()[1:] if "=" in x)
                if "pos" not in kv:
                    if not warned:
                        sys.stderr.write(
                            f"WARNING: {os.path.basename(path)} headers lack pos=; "
                            f"cannot build a census for this file.\n")
                        warned = True
                    return []
                out.append((int(kv["pos"]), len(seq)))
        else:
            name, ln = None, 0
            for line in fh:
                if line.startswith(">"):
                    if name is not None:
                        out.append((name, ln))
                    kv = dict(x.split("=", 1) for x in line.strip().split()[1:]
                              if "=" in x)
                    name = int(kv["pos"]) if "pos" in kv else None
                    ln = 0
                else:
                    ln += len(line.strip())
            if name is not None:
                out.append((name, ln))
            if any(n is None for n, _ in out):
                sys.stderr.write(f"WARNING: {os.path.basename(path)} lacks pos=\n")
                return []
    return out


def census(reads, events, frag_len, flank):
    """Per-event counts plus a per-sample host/donor split."""
    rows = []
    donor_reads = set()

    for e in events:
        ms, me = e["mock_start"], e["mock_end"]
        any_ = ge200 = gefrag = flanked = 0
        for i, (p, L) in enumerate(reads):
            ov = min(p + L, me) - max(p, ms)
            if ov <= 0:
                continue
            any_ += 1
            donor_reads.add(i)
            if ov >= 200:
                ge200 += 1
            if ov >= frag_len:
                gefrag += 1
            if p <= ms - flank and p + L >= me + flank:
                flanked += 1
        rows.append({
            "event_id": e["event_id"],
            "insert_len": e["insert_len"],
            "mock_start": ms,
            "mock_end": me,
            "reads_any": any_,
            "reads_ge200": ge200,
            "reads_ge_frag": gefrag,
            "reads_flanked": flanked,
        })
    return rows, len(donor_reads)


def main():
    p = argparse.ArgumentParser(
        description="Read-level truth census for HGT detection benchmarking")
    p.add_argument("--truth", required=True, help="all_truth.tsv")
    p.add_argument("--reads", required=True,
                   help="Directory of <sample>.long.fq.gz files")
    p.add_argument("--pattern", default="*.long.fq.gz",
                   help="Glob within --reads (default *.long.fq.gz)")
    p.add_argument("--frag-len", type=int, default=2000,
                   help="Fragment length in bp = read_len / --size (default 2000)")
    p.add_argument("--flank", type=int, default=1000,
                   help="Host bp required each side for reads_flanked (1000)")
    p.add_argument("--exclude", action="append", default=[],
                   help="Sample name to omit; repeatable")
    p.add_argument("--out-events", default="census_events.tsv")
    p.add_argument("--out-samples", default="census_samples.tsv")
    args = p.parse_args()

    truth = read_truth(args.truth)
    files = sorted(glob.glob(os.path.join(args.reads, args.pattern)))
    if not files:
        sys.exit(f"No files matching {args.pattern} in {args.reads}")

    ev_rows, s_rows = [], []
    skipped = []

    for f in files:
        sample = os.path.basename(f)
        for suf in (".long.fq.gz", ".long.fq", ".long.fa.gz", ".long.fa"):
            if sample.endswith(suf):
                sample = sample[: -len(suf)]
                break
        if sample in args.exclude:
            skipped.append(sample)
            continue

        reads = read_headers(f)
        if not reads:
            continue

        events = truth.get(sample, [])
        rows, n_donor = census(reads, events, args.frag_len, args.flank)
        for r in rows:
            r["sample"] = sample
            ev_rows.append(r)

        n = len(reads)
        total_bp = sum(L for _, L in reads)
        gl = max(p + L for p, L in reads) if reads else 0
        s_rows.append({
            "sample": sample,
            "is_negative": "yes" if not events else "no",
            "n_reads": n,
            "mean_read_len": round(total_bp / n, 1) if n else 0,
            "approx_genome_len": gl,
            "approx_coverage": round(total_bp / gl, 2) if gl else 0,
            "n_events": len(events),
            "reads_with_donor": n_donor,
            "reads_pure_host": n - n_donor,
            "pct_with_donor": round(100 * n_donor / n, 3) if n else 0,
        })

    ev_cols = ["sample", "event_id", "insert_len", "mock_start", "mock_end",
               "reads_any", "reads_ge200", "reads_ge_frag", "reads_flanked"]
    s_cols = ["sample", "is_negative", "n_reads", "mean_read_len",
              "approx_genome_len", "approx_coverage", "n_events",
              "reads_pure_host", "reads_with_donor", "pct_with_donor"]

    for path, rows, cols in ((args.out_events, ev_rows, ev_cols),
                             (args.out_samples, s_rows, s_cols)):
        with open(path, "w") as fh:
            fh.write("\t".join(cols) + "\n")
            for r in rows:
                fh.write("\t".join(str(r.get(c, "")) for c in cols) + "\n")

    # ── console summary ───────────────────────────────────────────────────────
    npos = sum(1 for r in s_rows if r["is_negative"] == "no")
    nneg = len(s_rows) - npos
    print(f"{len(s_rows)} samples ({npos} with events, {nneg} negative), "
          f"{len(ev_rows)} events")
    if skipped:
        print(f"excluded: {', '.join(skipped)}")

    by_size = defaultdict(lambda: [0, 0, 0, 0, 0])
    for r in ev_rows:
        a = by_size[r["insert_len"]]
        a[0] += 1
        for i, k in enumerate(("reads_any", "reads_ge200",
                               "reads_ge_frag", "reads_flanked"), start=1):
            a[i] += r[k]

    print(f"\nMean usable reads per event, by insert size "
          f"(fragment length {args.frag_len:,} bp):\n")
    print(f"{'insert':>9}{'n':>5}{'any':>9}{'>=200bp':>10}"
          f"{'>=1 frag':>10}{'flanked':>9}")
    for size in sorted(by_size):
        n, a, b, c, d = by_size[size]
        print(f"{size:>9,}{n:>5}{a/n:>9.1f}{b/n:>10.1f}{c/n:>10.1f}{d/n:>9.1f}")

    print("\nUse >=1 frag as the denominator for algorithmic recall: an insert "
          "\nsmaller than one fragment cannot be resolved regardless of depth.")
    print(f"\nWrote {args.out_events} and {args.out_samples}")


if __name__ == "__main__":
    main()