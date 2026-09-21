#!/usr/bin/env python3
"""
length.py — read-length bookkeeping for HERALD runs.

Three subcommands:

  table   Build a read_id -> length table from FASTA/FASTQ (plain or .gz). Run this before fragmenting.

  plan    given the length table, propose length bins and the -s value for each so that fragments 
  land at or above a target length. Reports the median (not the mean) per bin, since the long tail drags the
  mean upward and that is what produced shorter-than-intended fragments previously.

  join    Add a read_len column to candidates_annotated.tsv (or any TSV with a 'read' column) and report the candidate rate as a function of
  read length -- i.e. whether false positives come from short reads.

Examples
--------
  python3 length.py table reads.fastq.gz -o read_lengths.tsv
  python3 length.py plan read_lengths.tsv --bins 6 --min-frag 2000 --min-len 3000
  python3 length.py join candidates_annotated.tsv read_lengths.tsv -o candidates_with_len.tsv
"""

import argparse
import gzip
import os
import sys
from collections import Counter


def opener(path):
    return gzip.open(path, "rt") if path.endswith(".gz") else open(path, "r")


# ---------------------------------------------------------------- table

def cmd_table(args):
    n = 0
    with open(args.output, "w") as out:
        out.write("read\tlength\n")
        for path in args.reads:
            with opener(path) as fh:
                first = fh.readline()
                if not first:
                    continue
                fh.seek(0)
                if first.startswith("@"):        # FASTQ
                    for i, line in enumerate(fh):
                        if i % 4 == 0:
                            rid = line[1:].split()[0]
                        elif i % 4 == 1:
                            out.write(f"{rid}\t{len(line.strip())}\n")
                            n += 1
                else:                             # FASTA
                    rid, ln = None, 0
                    for line in fh:
                        if line.startswith(">"):
                            if rid is not None:
                                out.write(f"{rid}\t{ln}\n"); n += 1
                            rid, ln = line[1:].split()[0], 0
                        else:
                            ln += len(line.strip())
                    if rid is not None:
                        out.write(f"{rid}\t{ln}\n"); n += 1
    print(f"Wrote {args.output}: {n:,} reads")
    return 0


def load_lengths(path):
    d = {}
    with open(path) as fh:
        header = fh.readline()
        if not header.lower().startswith("read"):
            fh.seek(0)
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if len(f) >= 2:
                try:
                    d[f[0]] = int(f[1])
                except ValueError:
                    pass
    return d


def quantile(sorted_vals, q):
    if not sorted_vals:
        return 0
    i = q * (len(sorted_vals) - 1)
    lo, hi = int(i), min(int(i) + 1, len(sorted_vals) - 1)
    return sorted_vals[lo] + (sorted_vals[hi] - sorted_vals[lo]) * (i - lo)


# ---------------------------------------------------------------- plan

def cmd_plan(args):
    lengths = load_lengths(args.table)
    vals = sorted(lengths.values())
    total = len(vals)
    kept = [v for v in vals if v >= args.min_len]
    dropped = total - len(kept)

    print(f"Reads: {total:,}")
    print(f"  min {vals[0]:,}   Q1 {quantile(vals,.25):,.0f}   "
          f"median {quantile(vals,.5):,.0f}   Q3 {quantile(vals,.75):,.0f}   "
          f"max {vals[-1]:,}")
    print(f"  mean {sum(vals)/total:,.0f}  (median is lower -- the tail pulls "
          f"the mean up)")
    print(f"\nLength floor {args.min_len:,} bp drops {dropped:,} reads "
          f"({dropped/total:.2%}); {len(kept):,} remain.\n")
    if not kept:
        print("Nothing left after the floor.", file=sys.stderr)
        return 1

    # equal-count bins over the retained reads
    edges = [kept[0]]
    for i in range(1, args.bins):
        edges.append(quantile(kept, i / args.bins))
    edges.append(kept[-1])

    print(f"{'bin':>4s} {'range (bp)':>21s} {'reads':>10s} {'median':>9s} "
          f"{'-s':>4s} {'frag @median':>13s} {'frag @min':>10s}")
    print("-" * 78)
    plan = []
    for b in range(args.bins):
        lo = edges[b] if b == 0 else edges[b] + 1
        hi = edges[b + 1]
        sub = [v for v in kept if lo <= v <= hi]
        if not sub:
            continue
        med = quantile(sub, .5)
        # largest -s such that fragments at the bin MEDIAN stay >= min_frag,
        # and never below 2 (HERALD requires -s > 1)
        s = max(2, int(med // args.min_frag))
        frag_med = med / s
        frag_min = min(sub) / s
        warn = "  <-- short" if frag_min < args.min_frag * 0.5 else ""
        print(f"{b+1:4d} {int(lo):>9,d}-{int(hi):<11,d} {len(sub):10,d} "
              f"{med:9,.0f} {s:4d} {frag_med:13,.0f} {frag_min:10,.0f}{warn}")
        plan.append((b + 1, int(lo), int(hi), len(sub), s))

    print(f"\nFragments at each bin's shortest read are shown so you can see "
          f"the worst case,\nnot just the typical one. Narrow a bin or raise "
          f"the floor if that column is small.")

    if args.write_plan:
        with open(args.write_plan, "w") as out:
            out.write("bin\tmin_len\tmax_len\treads\ts_value\n")
            for row in plan:
                out.write("\t".join(str(x) for x in row) + "\n")
        print(f"\nWrote {args.write_plan}")
    return 0


# ---------------------------------------------------------------- join

def cmd_join(args):
    lengths = load_lengths(args.table)
    rows, missing = [], 0
    with open(args.candidates) as fh:
        header = fh.readline().rstrip("\n").split("\t")
        if "read" not in header:
            print("No 'read' column found.", file=sys.stderr)
            return 1
        ridx = header.index("read")
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if len(f) < len(header):
                continue
            L = lengths.get(f[ridx])
            if L is None:
                missing += 1
            rows.append((f, L))

    with open(args.output, "w") as out:
        out.write("\t".join(header + ["read_len"]) + "\n")
        for f, L in rows:
            out.write("\t".join(f + [str(L) if L is not None else "NA"]) + "\n")
    print(f"Wrote {args.output}: {len(rows):,} rows"
          + (f"  ({missing:,} with no length match)" if missing else ""))

    cand = sorted(L for _f, L in rows if L)
    allr = sorted(lengths.values())
    if not cand:
        return 0

    print(f"\nRead length: candidates vs all reads")
    print(f"  median, all reads  {quantile(allr,.5):9,.0f} bp")
    print(f"  median, candidates {quantile(cand,.5):9,.0f} bp"
          f"   ({quantile(cand,.5)/max(quantile(allr,.5),1):.2f}x)")

    # candidate rate per length decile of the full read set
    print(f"\nCandidate rate by read-length decile (of all reads)")
    print(f"{'decile':>7s} {'range (bp)':>21s} {'reads':>10s} {'cands':>7s} "
          f"{'rate':>9s}")
    print("-" * 60)
    edges = [quantile(allr, i / 10) for i in range(11)]
    stats = []
    for i in range(10):
        lo, hi = edges[i], edges[i + 1]
        nr = sum(1 for v in allr if lo <= v <= hi)
        nc = sum(1 for v in cand if lo <= v <= hi)
        stats.append((lo, hi, nr, nc, nc / nr if nr else 0))
    peak = max(r for *_x, r in stats) or 1
    for i, (lo, hi, nr, nc, rate) in enumerate(stats):
        bar = "#" * int(round(30 * rate / peak))
        print(f"{i+1:7d} {lo:>9,.0f}-{hi:<11,.0f} {nr:10,d} {nc:7,d} "
              f"{rate:8.4%} {bar}")
    print("\nA rate that climbs toward the short end means false positives are "
          "\ndriven by short reads; a flat profile means length is not the "
          "driver.")
    return 0


def main():
    ap = argparse.ArgumentParser(description=__doc__,
            formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest="cmd", required=True)

    t = sub.add_parser("table", help="build read_id -> length table")
    t.add_argument("reads", nargs="+")
    t.add_argument("-o", "--output", default="read_lengths.tsv")
    t.set_defaults(func=cmd_table)

    p = sub.add_parser("plan", help="propose length bins and -s values")
    p.add_argument("table")
    p.add_argument("--bins", type=int, default=6)
    p.add_argument("--min-frag", type=int, default=2000,
                   help="target minimum fragment length [2000]")
    p.add_argument("--min-len", type=int, default=3000,
                   help="discard reads shorter than this [3000]")
    p.add_argument("--write-plan", default=None)
    p.set_defaults(func=cmd_plan)

    j = sub.add_parser("join", help="add read_len to a candidate table")
    j.add_argument("candidates")
    j.add_argument("table")
    j.add_argument("-o", "--output", default="candidates_with_len.tsv")
    j.set_defaults(func=cmd_join)

    args = ap.parse_args()
    return args.func(args)


if __name__ == "__main__":
    try:
        sys.exit(main())
    except BrokenPipeError:
        os.dup2(os.open(os.devnull, os.O_WRONLY), sys.stdout.fileno())
        sys.exit(0)