#!/usr/bin/env python3
"""
donor_content.py — stratify reads by how much donor material they carry.

For each read, compute the number of base pairs overlapping a true insert interval, then 
cross-tabulate against HERALD's flagged candidates to give recall as a function of donor content.

INPUTS
------
--inserts   Either a truth TSV (auto-detected by its mock_start/mock_end
            header) or a BED of insert locations in the MOCK genome:
                mock_contig <tab> start <tab> end   [<tab> label]
            0-based, half-open. One line per simulated insert.

--reads     Where read coordinates come from. Either:
              a SAM/BAM of reads aligned back to the MOCK genome
              (recommended, and exact), or
              a FASTA whose headers encode position (--coords header).

--flagged   HERALD's candidate read IDs, one per line (or a SAM/TSV whose
            first whitespace-delimited field is the read ID). Optional; if
            omitted the script reports the donor-content distribution only.

USAGE
-----
  # Recommended: coordinates from an alignment to the mock genome
  python donor_content.py --inserts pair1_rep1.inserts.bed \
                          --reads pair1_rep1.tomock.bam \
                          --flagged pair1_rep1.herald_candidates.txt \
                          --fragment-length 2000 \
                          --sample pair1_rep1 --out-prefix results/pair1_rep1

  # Batch across all 28 samples
  python donor_content.py --manifest samples.tsv --fragment-length 2000 \
                          --out-prefix results/all

  # FASTA headers instead of an alignment
  python donor_content.py --inserts x.bed --reads x.fasta --coords header \
      --header-regex '(?P<contig>\\S+)_(?P<start>\\d+)_(?P<end>\\d+)'

MANIFEST FORMAT (--manifest, TSV with header)
  sample <tab> inserts <tab> reads <tab> flagged
"""

import argparse, csv, os, re, sys
from collections import defaultdict

# ---------------------------------------------------------------- intervals

def _merge(by_contig):
    for c, iv in by_contig.items():
        iv.sort()
        merged = [list(iv[0])]
        for s, e in iv[1:]:
            if s <= merged[-1][1]:
                merged[-1][1] = max(merged[-1][1], e)
            else:
                merged.append([s, e])
        by_contig[c] = [tuple(m) for m in merged]
    return dict(by_contig)


def load_truth(path, contig_col="acceptor_acc",
               start_col="mock_start", end_col="mock_end"):
    """Truth TSV -> {contig: [(start, end), ...]}.

    Uses mock_start/mock_end, which are coordinates in the MODIFIED (mock) genome — the same coordinate space the reads were sampled 
    from. Don't use ref_pos/donor_start: ref_pos is the insertion point in the unmodified acceptor, 
    and donor_* are coordinates in the donor genome. Only mock_* tells you where the insert actually sits in the sequence reads came from.

    The contig name must match the reference name in your alignment. If your mock FASTA headers are not the bare acceptor 
    accession, pass --contig-col or --contig-name.
    """
    by_contig = defaultdict(list)
    n = 0
    with open(path) as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            if not row.get(start_col):
                continue
            by_contig[row[contig_col]].append((
                int(row[start_col]), int(row[end_col]),
                row.get("event_id", f"event{n + 1}"),
                int(row["insert_len"]) if row.get("insert_len") else None))
            n += 1
    if not n:
        raise SystemExit(f"no intervals parsed from {path}")
    for c in by_contig:
        by_contig[c].sort()
        iv = by_contig[c]
        for a, b in zip(iv, iv[1:]):
            if b[0] < a[1]:
                print(f"  warning: events {a[2]} and {b[2]} overlap; donor bp "
                      "will be attributed to the first", file=sys.stderr)
    return dict(by_contig)


def load_inserts(path):
    """BED -> {contig: [(start, end), ...]} with overlapping intervals merged.

    If two simulated inserts abut or overlap, counting them separately would double-count shared bases and inflate donor content.
    """
    by_contig = defaultdict(list)
    with open(path) as fh:
        for line in fh:
            if not line.strip() or line.startswith(("#", "track", "browser")):
                continue
            f = line.split()
            by_contig[f[0]].append((int(f[1]), int(f[2])))
    merged = _merge(by_contig)
    return {c: [(s_, e_, f"{c}:{s_}-{e_}", e_ - s_) for s_, e_ in iv]
            for c, iv in merged.items()}


def _looks_like_truth(path):
    """True if the file has the truth-TSV header (mock_start/mock_end)."""
    with open(path) as fh:
        head = fh.readline()
    return "mock_start" in head and "mock_end" in head


def overlap_bp(read_start, read_end, intervals):
    """Donor bp in [read_start, read_end), plus the event contributing most.

    Returns (total_bp, event_id, insert_len, bp_from_that_event). Attribution goes to the single largest 
    contributor, which matters only for reads spanning two events -- rare here, since inserts are >=100 kb apart.
    """
    total = 0
    best = (0, None, None)
    for s, e, eid, ilen in intervals:
        if s >= read_end:
            break            # intervals are sorted; nothing further can overlap
        if e <= read_start:
            continue
        ov = min(e, read_end) - max(s, read_start)
        total += ov
        if ov > best[0]:
            best = (ov, eid, ilen)
    return total, best[1], best[2], best[0]

# ---------------------------------------------------------------- read coords

# Reference-consuming CIGAR ops. Soft clips (S) and insertions (I) advance the read but not the reference, 
# so they must not extend the reference span.
_CIG = re.compile(r"(\d+)([MIDNSHP=X])")

def _ref_span(cigar):
    return sum(int(n) for n, op in _CIG.findall(cigar) if op in "MDN=X")


def reads_from_alignment(path):
    """Yield (read_id, contig, start, end) for primary, mapped alignments.

    Secondary (0x100) and supplementary (0x800) records are skipped so each read contributes exactly one interval. 
    Requires `samtools` for BAM/CRAM.
    """
    if path.endswith((".bam", ".cram")):
        import subprocess
        proc = subprocess.Popen(["samtools", "view", path],
                                stdout=subprocess.PIPE, text=True)
        stream, close = proc.stdout, proc
    else:
        stream, close = open(path), None

    seen = set()
    try:
        for line in stream:
            if line.startswith("@"):
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 6:
                continue
            flag = int(f[1])
            if flag & 0x4 or flag & 0x100 or flag & 0x800:
                continue
            rid = f[0]
            if rid in seen:                 # defensive: duplicate primaries
                continue
            seen.add(rid)
            start = int(f[3]) - 1           # SAM is 1-based
            yield rid, f[2], start, start + _ref_span(f[5])
    finally:
        if close:
            close.stdout.close(); close.wait()
        else:
            stream.close()


def reads_from_fasta(path, regex):
    """Parse coordinates out of FASTA headers using a named-group regex.

    Only use this if your simulator writes true positions into the header. 
    """
    rx = re.compile(regex)
    with open(path) as fh:
        for line in fh:
            if not line.startswith(">"):
                continue
            head = line[1:].rstrip("\n")
            m = rx.search(head)
            if not m:
                print(f"  warning: header did not match regex: {head[:70]}",
                      file=sys.stderr)
                continue
            g = m.groupdict()
            yield head.split()[0], g["contig"], int(g["start"]), int(g["end"])


_FRAG_SUFFIX = re.compile(r"_frag_\d+$")


def load_flagged(path):
    """Read IDs HERALD flagged.

    Accepts a plain ID list, a SAM, or HERALD's fragment-level results file:

        pair4_..._read000035_frag_1 <tab> CP024307.1 <tab> AS:i:3940
        pair4_..._read000035_frag_2 <tab> CP024307.1 <tab> AS:i:3910
        <blank line between reads>

    A trailing _frag_N is stripped so every fragment collapses to its parent read, and the set deduplicates. 
    The parent ID must match the read name in the alignment. Make sure the SAM/BAM was built from the same FASTA.
    """
    ids = set()
    if not path:
        return ids
    with open(path) as fh:
        for line in fh:
            if line.startswith("@") or not line.strip():
                continue
            ids.add(_FRAG_SUFFIX.sub("", line.split()[0]))
    return ids

# ---------------------------------------------------------------- binning

def make_bins(frag_len):
    """Bin edges in bp. The final boundary is one fragment length."""
    edges = [1, 200, 500, 1000, frag_len, float("inf")]
    labels = []
    for i in range(len(edges) - 1):
        lo, hi = edges[i], edges[i + 1]
        labels.append(f">={int(lo)}" if hi == float("inf")
                      else f"{int(lo)}-{int(hi) - 1}")
    return edges, labels


def assign_bin(bp, edges):
    for i in range(len(edges) - 1):
        if edges[i] <= bp < edges[i + 1]:
            return i
    return None

# ---------------------------------------------------------------- Wilson CI

def wilson(k, n, z=1.96):
    """Wilson score interval — behaves sensibly at 0% and 100%, unlike the normal approximation, 
    which matters for the empty low-donor bins."""
    if n == 0:
        return (float("nan"), float("nan"))
    p = k / n
    d = 1 + z * z / n
    c = (p + z * z / (2 * n)) / d
    h = z * ((p * (1 - p) / n + z * z / (4 * n * n)) ** 0.5) / d
    return (max(0.0, c - h), min(1.0, c + h))

# ---------------------------------------------------------------- core

def run_sample(sample, inserts_path, reads_path, flagged_path,
               frag_len, coords, header_regex, per_read_writer,
               contig_col="acceptor_acc", contig_name=None):
    if inserts_path.lower().endswith((".tsv", ".txt")) and _looks_like_truth(inserts_path):
        inserts = load_truth(inserts_path, contig_col=contig_col)
    else:
        inserts = load_inserts(inserts_path)
    if contig_name:
        pooled = []
        for iv in inserts.values():
            pooled.extend(iv)
        pooled.sort()
        inserts = {contig_name: pooled}
    flagged = load_flagged(flagged_path)
    edges, labels = make_bins(frag_len)

    reader = (reads_from_alignment(reads_path) if coords == "align"
              else reads_from_fasta(reads_path, header_regex))

    seen_contigs = set()

    n_reads = n_carrying = 0
    tot = [0] * (len(edges) - 1)
    det = [0] * (len(edges) - 1)
    fp_no_donor = 0
    ev_reads = defaultdict(int)
    ev_det = defaultdict(int)
    ev_len = {}

    for rid, contig, start, end in reader:
        n_reads += 1
        seen_contigs.add(contig)
        bp, eid, ilen, ev_bp = overlap_bp(start, end, inserts.get(contig, []))
        hit = rid in flagged

        if bp == 0:
            # Flagged with zero donor material = false positive.
            if hit:
                fp_no_donor += 1
            continue

        n_carrying += 1
        b = assign_bin(bp, edges)
        tot[b] += 1
        if hit:
            det[b] += 1
        ev_reads[eid] += 1
        ev_det[eid] += int(hit)
        ev_len[eid] = ilen

        if per_read_writer:
            per_read_writer.writerow([sample, rid, contig, start, end,
                                      end - start, bp,
                                      round(bp / max(1, end - start), 4),
                                      labels[b], eid, ilen, ev_bp, int(hit)])

    # Silent zero-overlap is almost always a contig-naming mismatch, not a real absence of donor material. 
    if n_carrying == 0 and n_reads:
        refs = sorted(seen_contigs)
        hint = refs[0] if len(refs) == 1 else "<name from Alignment refs>"
        raise SystemExit(
            f"\n[{sample}] ERROR: no read overlapped any insert interval.\n"
            f"  Insert contigs : {sorted(inserts)[:4]}\n"
            f"  Alignment refs : {refs[:4]}\n"
            "  Your truth file and your alignment use different names for the\n"
            "  same sequence. Rerun with the name from 'Alignment refs':\n"
            f"      --contig-name {hint}\n"
            "  (That is the contig name in your mock FASTA header, not the\n"
            "  SAM filename and not the acceptor accession.)")

    return {"sample": sample, "n_reads": n_reads, "n_carrying": n_carrying,
            "labels": labels, "total": tot, "detected": det,
            "fp_no_donor": fp_no_donor,
            "events": {e: (ev_reads[e], ev_det[e], ev_len[e]) for e in ev_reads}}

# ---------------------------------------------------------------- reporting

def report(results, frag_len, out_prefix):
    labels = results[0]["labels"]
    nb = len(labels)
    tot = [sum(r["total"][i] for r in results) for i in range(nb)]
    det = [sum(r["detected"][i] for r in results) for i in range(nb)]
    fp = sum(r["fp_no_donor"] for r in results)
    n_reads = sum(r["n_reads"] for r in results)
    n_carry = sum(r["n_carrying"] for r in results)

    rows = []
    for i in range(nb):
        lo, hi = wilson(det[i], tot[i])
        rows.append([labels[i], tot[i], det[i],
                     (det[i] / tot[i] if tot[i] else float("nan")), lo, hi])

    w = max(12, max(len(l) for l in labels) + 2)
    print(f"\n{'Donor bp in read':<{w}} {'Reads':>9} {'Detected':>9} "
          f"{'Recall':>8}  95% CI")
    print("-" * (w + 46))
    for lab, n, k, r, lo, hi in rows:
        ci = "         —" if n == 0 else f"  [{lo:.3f}, {hi:.3f}]"
        rr = "     —" if n == 0 else f"{r:7.1%}"
        print(f"{lab:<{w}} {n:>9,} {k:>9,} {rr}{ci}")
    print("-" * (w + 46))

    uncond_k, uncond_n = sum(det), sum(tot)
    hi_bin = nb - 1
    print(f"{'TOTAL (>=1 bp)':<{w}} {uncond_n:>9,} {uncond_k:>9,} "
          f"{uncond_k / uncond_n if uncond_n else float('nan'):7.1%}")
    if tot[hi_bin]:
        print(f"{'>=1 frag len':<{w}} {tot[hi_bin]:>9,} {det[hi_bin]:>9,} "
              f"{det[hi_bin] / tot[hi_bin]:7.1%}")

    ev_tot = ev_hit = 0
    by_len = defaultdict(lambda: [0, 0])
    for r in results:
        for eid, (nr, nd, ilen) in r["events"].items():
            ev_tot += 1
            ev_hit += int(nd > 0)
            by_len[ilen][0] += 1
            by_len[ilen][1] += int(nd > 0)
    if by_len:
        print(f"\n{'Insert length':<{w}} {'Events':>9} {'Detected':>9} {'Recall':>8}")
        print("-" * (w + 30))
        for ilen in sorted(by_len, key=lambda v: (v is None, v)):
            n, k = by_len[ilen]
            lab = f"{ilen:,} bp" if ilen is not None else "unknown"
            print(f"{lab:<{w}} {n:>9,} {k:>9,} {k / n:7.1%}")
        print("-" * (w + 30))
        print(f"{'ALL EVENTS':<{w}} {ev_tot:>9,} {ev_hit:>9,} "
              f"{ev_hit / ev_tot if ev_tot else float('nan'):7.1%}")

    print(f"\nReads examined:                  {n_reads:,}")
    print(f"Reads overlapping an insert:     {n_carry:,} "
          f"({n_carry / n_reads:.2%})" if n_reads else "")
    print(f"Flagged with zero donor content: {fp:,}   <- false positives")
    print(f"Fragment length assumed:         {frag_len:,} bp")

    if uncond_n and tot[hi_bin]:
        print(f"\nUnconditional recall (>=1 bp of donor): "
              f"{uncond_k / uncond_n:.1%}")
        print(f"Recall among reads carrying >=1 fragment length: "
              f"{det[hi_bin] / tot[hi_bin]:.1%}")
        print("Report both. The first is end-to-end sensitivity; the second\n"
              "describes what the method can resolve. Neither replaces the other.")

    if out_prefix:
        os.makedirs(os.path.dirname(out_prefix) or ".", exist_ok=True)
        with open(f"{out_prefix}.strata.tsv", "w", newline="") as fh:
            wtr = csv.writer(fh, delimiter="\t")
            wtr.writerow(["donor_bp_bin", "reads", "detected", "recall",
                          "ci_low", "ci_high"])
            wtr.writerows(rows)
        with open(f"{out_prefix}.by_sample.tsv", "w", newline="") as fh:
            wtr = csv.writer(fh, delimiter="\t")
            wtr.writerow(["sample", "donor_bp_bin", "reads", "detected",
                          "recall"])
            for r in results:
                for i, lab in enumerate(labels):
                    n, k = r["total"][i], r["detected"][i]
                    wtr.writerow([r["sample"], lab, n, k,
                                  f"{k / n:.4f}" if n else ""])
        with open(f"{out_prefix}.events.tsv", "w", newline="") as fh:
            wtr = csv.writer(fh, delimiter="\t")
            wtr.writerow(["sample", "event_id", "insert_len", "reads_overlapping",
                          "reads_detected", "event_detected"])
            for r in results:
                for eid, (nr, nd, ilen) in sorted(r["events"].items(),
                                                  key=lambda kv: str(kv[0])):
                    wtr.writerow([r["sample"], eid, ilen if ilen is not None else "",
                                  nr, nd, int(nd > 0)])
        print(f"\nWrote {out_prefix}.strata.tsv, {out_prefix}.by_sample.tsv, "
              f"and {out_prefix}.events.tsv")


def main():
    ap = argparse.ArgumentParser(
        description="Stratify reads by donor content and report recall per stratum.")
    ap.add_argument("--manifest", help="TSV: sample, inserts, reads, flagged")
    ap.add_argument("--inserts"); ap.add_argument("--reads")
    ap.add_argument("--flagged")
    ap.add_argument("--sample", default="sample")
    ap.add_argument("--fragment-length", type=int, default=2000,
                    help="read length / n_fragments (default 2000 = 20kb/10)")
    ap.add_argument("--coords", choices=["align", "header"], default="align")
    ap.add_argument("--contig-col", default="acceptor_acc",
                    help="truth-TSV column holding the mock contig name")
    ap.add_argument("--contig-name", default=None,
                    help="force all inserts onto this contig name (use when "
                         "the mock FASTA header is not the accession)")
    ap.add_argument("--header-regex",
                    default=r"(?P<contig>\S+)[_:](?P<start>\d+)[-_](?P<end>\d+)")
    ap.add_argument("--out-prefix", default="")
    a = ap.parse_args()

    jobs = []
    if a.manifest:
        with open(a.manifest) as fh:
            for row in csv.DictReader(fh, delimiter="\t"):
                jobs.append((row["sample"], row["inserts"], row["reads"],
                             row.get("flagged")))
    elif a.inserts and a.reads:
        jobs.append((a.sample, a.inserts, a.reads, a.flagged))
    else:
        ap.error("provide --manifest, or both --inserts and --reads")

    per_read_fh = per_read_writer = None
    if a.out_prefix:
        os.makedirs(os.path.dirname(a.out_prefix) or ".", exist_ok=True)
        per_read_fh = open(f"{a.out_prefix}.per_read.tsv", "w", newline="")
        per_read_writer = csv.writer(per_read_fh, delimiter="\t")
        per_read_writer.writerow(["sample", "read_id", "contig", "start", "end",
                                  "read_len", "donor_bp", "donor_frac",
                                  "bin", "event_id", "insert_len",
                                  "event_donor_bp", "flagged"])

    results = []
    for sample, ins, rds, flg in jobs:
        print(f"[{sample}] {os.path.basename(rds)}", file=sys.stderr)
        results.append(run_sample(sample, ins, rds, flg, a.fragment_length,
                                  a.coords, a.header_regex, per_read_writer,
                                  a.contig_col, a.contig_name))

    if per_read_fh:
        per_read_fh.close()
    report(results, a.fragment_length, a.out_prefix)


if __name__ == "__main__":
    main()