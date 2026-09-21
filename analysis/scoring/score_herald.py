#!/usr/bin/env python3
"""
score_herald.py — score HERALD fragment results against ground truth.

What it does
------------
Reads 1+ HERALD post-processing result files, matches every flagged read back to the truth 
table using the pos= field in the read header, and reports recall, precision, misattribution, and false-positive rates.

Classification of each flagged read
-----------------------------------
  TRUE POSITIVE   overlaps a true event by >=--min-overlap bp AND carries at
                  least one fragment assigned to that event's donor genome
  MISATTRIBUTED   overlaps a true event, but no fragment matches the true
                  donor - the read is real mosaic, credited to the wrong genome.
                  Counted separately because it is neither a clean hit nor a
                  spurious call, and it is the failure mode expected to appear
                  when the reference panel grows.
  FALSE POSITIVE  overlaps no true event

Negative-control samples have no rows in the truth table, so every flagged read
there is a false positive by construction.

Event-level detection
---------------------
An event counts as detected when at least --min-reads reads support it. 


Usage
-----
  python score_herald.py --truth mock/all_truth.tsv --reads reads/ \\
      --results herald_out/panel/ --arm panel \\
      --min-overlap 200 --min-reads 2 --exclude pair1_phylum_rep1

  # a parameter sweep: each result file is scored as its own run
  python score_herald.py --truth mock/all_truth.tsv --reads reads/ \\
      --results "sweep/*.txt" --arm sweep
"""

import argparse
import glob
import gzip
import os
import re
import sys
from collections import defaultdict

READ_RE = re.compile(r"^(?P<sample>.+?)_read(?P<num>\d+)_frag_(?P<frag>\d+)$")


# ── inputs ────────────────────────────────────────────────────────────────────

def open_maybe_gz(p):
    return gzip.open(p, "rt") if p.endswith(".gz") else open(p)


def read_truth(path):
    events = defaultdict(list)
    with open(path) as fh:
        hdr = fh.readline().rstrip("\n").split("\t")
        need = {"sample", "event_id", "insert_len", "mock_start", "mock_end",
                "acceptor_acc", "donor_acc"}
        miss = need - set(hdr)
        if miss:
            sys.exit(f"{path} missing column(s): {', '.join(sorted(miss))}")
        for line in fh:
            if not line.strip():
                continue
            d = dict(zip(hdr, line.rstrip("\n").split("\t")))
            for k in ("insert_len", "mock_start", "mock_end"):
                d[k] = int(d[k])
            events[d["sample"]].append(d)
    return events


def read_positions(reads_dir, sample, pattern="{sample}.long.fq.gz"):
    """Return {read_name: (pos, len)} for one sample, or None if unavailable."""
    path = os.path.join(reads_dir, pattern.format(sample=sample))
    if not os.path.exists(path):
        for alt in (".long.fq", ".long.fa.gz", ".long.fa"):
            cand = os.path.join(reads_dir, sample + alt)
            if os.path.exists(cand):
                path = cand
                break
        else:
            return None
    out = {}
    with open_maybe_gz(path) as fh:
        fastq = path.endswith((".fq", ".fq.gz"))
        while True:
            h = fh.readline()
            if not h:
                break
            if fastq:
                seq = fh.readline().rstrip("\n")
                fh.readline()
                fh.readline()
            else:
                if not h.startswith(">"):
                    continue
                seq = fh.readline().rstrip("\n")
            parts = h[1:].strip().split()
            kv = dict(x.split("=", 1) for x in parts[1:] if "=" in x)
            if "pos" not in kv:
                return None
            out[parts[0]] = (int(kv["pos"]), len(seq))
    return out


def parse_results(path):
    """Return {sample: {read_name: [genome, ...]}}. Blank lines are separators."""
    per = defaultdict(lambda: defaultdict(list))
    bad = 0
    with open(path) as fh:
        for line in fh:
            line = line.rstrip("\n")
            if not line.strip():
                continue
            f = line.split("\t")
            if len(f) < 2:
                bad += 1
                continue
            m = READ_RE.match(f[0])
            if not m:
                bad += 1
                continue
            read = f[0][: f[0].rindex("_frag_")]
            per[m.group("sample")][read].append(f[1])
    return per, bad


# ── scoring ───────────────────────────────────────────────────────────────────

def score_sample(sample, calls, positions, events, min_overlap):
    """calls: {read: [genomes]}. Returns (per-event counts, read tallies)."""
    ev_hits = {e["event_id"]: 0 for e in events}
    tp = mis = fp = unknown = 0
    fp_detail = defaultdict(int)

    for read, genomes in calls.items():
        if positions is None or read not in positions:
            unknown += 1
            continue
        p, L = positions[read]
        overlapping = [e for e in events
                       if min(p + L, e["mock_end"]) - max(p, e["mock_start"])
                       >= min_overlap]
        if not overlapping:
            fp += 1
            fp_detail["+".join(sorted(set(genomes)))] += 1
            continue
        gset = set(genomes)
        credited = [e for e in overlapping if e["donor_acc"] in gset]
        if credited:
            tp += 1
            for e in credited:
                ev_hits[e["event_id"]] += 1
        else:
            mis += 1

    return ev_hits, {"tp": tp, "misattributed": mis, "fp": fp,
                     "unresolved": unknown, "fp_detail": dict(fp_detail)}


def denominator(positions, e, min_overlap):
    if positions is None:
        return 0
    return sum(1 for p, L in positions.values()
               if min(p + L, e["mock_end"]) - max(p, e["mock_start"]) >= min_overlap)


# ── output ────────────────────────────────────────────────────────────────────

def write_tsv(path, rows, cols):
    with open(path, "w") as fh:
        fh.write("\t".join(cols) + "\n")
        for r in rows:
            fh.write("\t".join(str(r.get(c, "")) for c in cols) + "\n")


def main():
    p = argparse.ArgumentParser(description="Score HERALD results vs truth")
    p.add_argument("--truth", required=True)
    p.add_argument("--reads", required=True, help="Directory of long-read files")
    p.add_argument("--results", required=True,
                   help="Directory of result .txt files, or a quoted glob")
    p.add_argument("--arm", default="run", help="Label for this analysis arm")
    p.add_argument("--min-overlap", type=int, default=200,
                   help="Donor bp in a read to count as carrying the event (200)")
    p.add_argument("--min-reads", type=int, default=2,
                   help="Supporting reads for event-level detection (2)")
    p.add_argument("--exclude", action="append", default=[])
    p.add_argument("--out-prefix", default="panel",
                   help="Output filename prefix; writes <prefix>_events.tsv and\n<prefix>_samples.tsv (default: panel)")
    args = p.parse_args()

    # --results accepts a single file, a directory, or a glob pattern.
    if any(c in args.results for c in "*?["):
        files = sorted(glob.glob(args.results))
    elif os.path.isfile(args.results):
        files = [args.results]
    elif os.path.isdir(args.results):
        files = sorted(glob.glob(os.path.join(args.results, "*.txt")))
        if not files:
            sys.exit(f"{args.results} is a directory but contains no *.txt files")
    else:
        sys.exit(f"--results path does not exist: {args.results}")
    if not files:
        sys.exit(f"No result files found for {args.results}")

    truth = read_truth(args.truth)
    pos_cache = {}

    ev_rows, s_rows = [], []
    total_bad = 0

    for rf in files:
        run = os.path.basename(rf)
        for suf in (".txt", ".tsv"):
            if run.endswith(suf):
                run = run[: -len(suf)]
        per_sample, bad = parse_results(rf)
        total_bad += bad

        for sample, calls in per_sample.items():
            if sample in args.exclude:
                continue
            if sample not in pos_cache:
                pos_cache[sample] = read_positions(args.reads, sample)
            positions = pos_cache[sample]
            if positions is None:
                sys.stderr.write(
                    f"WARNING: no readable long-read file with pos= for "
                    f"{sample}; its reads cannot be scored.\n")
            events = truth.get(sample, [])
            ev_hits, tally = score_sample(sample, calls, positions,
                                          events, args.min_overlap)

            detected = 0
            for e in events:
                n = ev_hits[e["event_id"]]
                den = denominator(positions, e, args.min_overlap)
                ok = n >= args.min_reads
                detected += int(ok)
                ev_rows.append({
                    "arm": args.arm, "run": run, "sample": sample,
                    "event_id": e["event_id"], "insert_len": e["insert_len"],
                    "reads_available": den, "reads_detected": n,
                    "read_recall_pct": round(100 * n / den, 2) if den else "",
                    "detected": "yes" if ok else "no",
                })

            nreads = len(positions) if positions else 0
            flagged = len(calls)
            tp, fp = tally["tp"], tally["fp"]
            s_rows.append({
                "arm": args.arm, "run": run, "sample": sample,
                "is_negative": "yes" if not events else "no",
                "total_reads": nreads, "reads_flagged": flagged,
                "true_positive": tp, "misattributed": tally["misattributed"],
                "false_positive": fp, "unresolved": tally["unresolved"],
                "precision_pct": round(100 * tp / (tp + fp), 2) if tp + fp else "",
                "fp_per_100k_reads": round(1e5 * fp / nreads, 2) if nreads else "",
                "events_total": len(events), "events_detected": detected,
                "fp_genomes": "; ".join(f"{k}:{v}" for k, v in
                                        sorted(tally["fp_detail"].items(),
                                               key=lambda x: -x[1])[:3]),
            })

    ev_cols = ["arm", "run", "sample", "event_id", "insert_len",
               "reads_available", "reads_detected", "read_recall_pct", "detected"]
    s_cols = ["arm", "run", "sample", "is_negative", "total_reads",
              "reads_flagged", "true_positive", "misattributed",
              "false_positive", "unresolved", "precision_pct",
              "fp_per_100k_reads", "events_total", "events_detected",
              "fp_genomes"]

    events_path = f"{args.out_prefix}_events.tsv"
    samples_path = f"{args.out_prefix}_samples.tsv"
    write_tsv(events_path, ev_rows, ev_cols)
    write_tsv(samples_path, s_rows, s_cols)

    # ── console summary ───────────────────────────────────────────────────────
    if total_bad:
        print(f"note: {total_bad} unparsed lines across {len(files)} file(s) "
              f"(blank separators are ignored, not counted here)\n")

    for run in sorted({r["run"] for r in s_rows}):
        srun = [r for r in s_rows if r["run"] == run]
        erun = [r for r in ev_rows if r["run"] == run]
        pos_s = [r for r in srun if r["is_negative"] == "no"]
        neg_s = [r for r in srun if r["is_negative"] == "yes"]

        print(f"=== {run} ===")
        print(f"{len(srun)} sample(s): {len(pos_s)} with events, "
              f"{len(neg_s)} negative")

        by = defaultdict(lambda: [0, 0, 0, 0])
        for r in erun:
            a = by[r["insert_len"]]
            a[0] += 1
            a[1] += r["reads_available"]
            a[2] += r["reads_detected"]
            a[3] += int(r["detected"] == "yes")
        if by:
            print(f"\n{'insert':>9}{'events':>8}{'avail':>8}{'found':>8}"
                  f"{'read recall':>13}{'events found':>14}")
            for size in sorted(by):
                n, av, fd, ok = by[size]
                print(f"{size:>9,}{n:>8}{av:>8}{fd:>8}"
                      f"{(100*fd/av if av else 0):>12.1f}%"
                      f"{f'{ok}/{n}':>14}")

        tp = sum(r["true_positive"] for r in srun)
        mis = sum(r["misattributed"] for r in srun)
        fp = sum(r["false_positive"] for r in srun)
        nr = sum(r["total_reads"] for r in srun)
        ed = sum(r["events_detected"] for r in srun)
        et = sum(r["events_total"] for r in srun)
        print(f"\nreads: {tp} true positive, {mis} misattributed, "
              f"{fp} false positive")
        if tp + fp:
            print(f"precision {100*tp/(tp+fp):.2f}%")
        if nr:
            print(f"FP rate  {1e5*fp/nr:.2f} per 100,000 reads "
                  f"({fp}/{nr:,})")
        if et:
            print(f"events   {ed}/{et} detected "
                  f"(>={args.min_reads} supporting reads)")

        if neg_s:
            bad_negs = [r for r in neg_s if r["false_positive"] > 0]
            nrn = sum(r["total_reads"] for r in neg_s)
            fpn = sum(r["false_positive"] for r in neg_s)
            print(f"\nnegative controls: {fpn} FP across {len(neg_s)} sample(s)"
                  + (f", {1e5*fpn/nrn:.2f} per 100,000 reads" if nrn else ""))
            for r in bad_negs:
                print(f"  {r['sample']}: {r['false_positive']} "
                      f"({r['fp_genomes']})")
            if not bad_negs:
                print("  all clean")
        print()

    print(f"Wrote {events_path} and {samples_path}")


if __name__ == "__main__":
    main()