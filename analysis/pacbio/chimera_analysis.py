#!/usr/bin/env python3
"""
chimera_analysis.py — measure the chimeric-read rate directly from a first-round minimap2 SAM, 
independently of HERALD.

Secondary alignments (FLAG 256) are ignored: they are alternative placements of the same piece of the read,
not evidence of a junction.

Usage
-----
  python3 chimera_analysis.py first_round_*.sam
  python3 chimera_analysis.py aln.sam.gz --min-seg-frac 0.15 --min-mapq 20
  python3 chimera_analysis.py *.sam --out chimeras.tsv

Notes
-----
* Assumes all alignment lines for a given read are adjacent.  Pass --unsorted to buffer the whole file by read
  name instead (uses more memory).
* Handles both soft (S) and hard (H) clipping, so supplementary alignments are placed correctly on the original read.
* Reads may be gzipped.
"""

import argparse
import gzip
import os
import re
import sys
from collections import Counter, defaultdict

CIGAR_RE = re.compile(r'(\d+)([MIDNSHP=X])')

# Operations that consume bases on the query (the read)
QUERY_OPS = set("MIS=X")
# Operations that consume bases on the reference
REF_OPS = set("MDN=X")


def opener(path):
    if path.endswith(".gz"):
        return gzip.open(path, "rt")
    return open(path, "r")


def parse_cigar(cigar):
    """Return (query_start, query_end, ref_span, full_read_len).

    query_start/end are 0-based coordinates on the ORIGINAL read, correct for both soft- and hard-clipped records.
    """
    ops = CIGAR_RE.findall(cigar)
    if not ops:
        return None
    lead = 0
    for n, op in ops:
        if op in "SH":
            lead += int(n)
        else:
            break
    trail = 0
    for n, op in reversed(ops):
        if op in "SH":
            trail += int(n)
        else:
            break
    aligned_q = sum(int(n) for n, op in ops if op in QUERY_OPS and op != "S")
    ref_span = sum(int(n) for n, op in ops if op in REF_OPS)
    full_len = lead + aligned_q + trail
    return lead, lead + aligned_q, ref_span, full_len


def genome_of(rname, contig_map):
    """Map a reference sequence name to a genome name."""
    return contig_map.get(rname, rname)


def load_contig_map(path):
    """Optional two-column TSV: reference_name <TAB> genome_name."""
    m = {}
    if not path:
        return m
    with open(path) as fh:
        for line in fh:
            if not line.strip() or line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) >= 2:
                m[f[0]] = f[1]
    return m


def classify_read(alns, min_seg_frac, min_mapq, min_total_cov):
    """
    alns: list of dicts for one read (primary + supplementary only).

    Returns None if the read is not a usable multi-segment read, else a dict describing it.
    """
    if len(alns) < 2:
        return None

    read_len = max(a["full_len"] for a in alns)
    if read_len <= 0:
        return None

    segs = [a for a in alns
            if (a["qend"] - a["qstart"]) >= min_seg_frac * read_len
            and a["mapq"] >= min_mapq]
    if len(segs) < 2:
        return None

    segs.sort(key=lambda a: a["qstart"])

    # Total fraction of the read explained by the retained segments
    covered, cursor = 0, -1
    for s in segs:
        start = max(s["qstart"], cursor)
        if s["qend"] > start:
            covered += s["qend"] - start
            cursor = s["qend"]
    if covered < min_total_cov * read_len:
        return None

    genomes = [s["genome"] for s in segs]
    distinct = sorted(set(genomes))

    # Junction = midpoint of the gap/overlap between consecutive segments
    junctions = []
    for a, b in zip(segs, segs[1:]):
        junctions.append((a["qend"] + b["qstart"]) / 2.0 / read_len)

    if len(distinct) > 1:
        kind = "inter_genome"
    else:
        # Same genome: is it a real rearrangement/repeat jump, or just a
        # contiguous alignment split into two records?
        gaps = []
        for a, b in zip(segs, segs[1:]):
            if a["strand"] != b["strand"]:
                gaps.append(float("inf"))
            else:
                gaps.append(abs(b["rpos"] - a["rend"]))
        kind = "intra_genome" if max(gaps) > 10000 else "contiguous_split"

    return dict(
        read=segs[0]["qname"],
        read_len=read_len,
        n_segments=len(segs),
        kind=kind,
        genomes=",".join(distinct),
        pair=",".join(distinct[:2]) if len(distinct) > 1 else "",
        junction_frac=round(junctions[0], 4) if junctions else None,
        junction_bp=int(junctions[0] * read_len) if junctions else None,
        covered_frac=round(covered / read_len, 4),
        strands="".join(s["strand"] for s in segs),
    )


def iter_reads(paths, unsorted_mode=False):
    """Yield (qname, [alignment dicts]) groups, skipping secondary records."""
    for path in paths:
        buf, cur = [], None
        pool = defaultdict(list) if unsorted_mode else None
        with opener(path) as fh:
            for line in fh:
                if line[0] == "@":
                    continue
                f = line.rstrip("\n").split("\t")
                if len(f) < 11:
                    continue
                qname, flag, rname, pos, mapq, cigar = (
                    f[0], int(f[1]), f[2], int(f[3]), int(f[4]), f[5])

                if flag & 4 or rname == "*" or cigar == "*":
                    continue
                if flag & 256:            # secondary: not a junction
                    continue

                parsed = parse_cigar(cigar)
                if not parsed:
                    continue
                qstart, qend, ref_span, full_len = parsed

                rec = dict(
                    qname=qname, rname=rname, mapq=mapq,
                    qstart=qstart, qend=qend, full_len=full_len,
                    rpos=pos, rend=pos + ref_span,
                    strand="-" if flag & 16 else "+",
                    supplementary=bool(flag & 2048),
                )

                if unsorted_mode:
                    pool[qname].append(rec)
                    continue

                if cur is not None and qname != cur:
                    yield cur, buf
                    buf = []
                cur = qname
                buf.append(rec)

        if unsorted_mode:
            for k, v in pool.items():
                yield k, v
        elif buf:
            yield cur, buf


def main():
    ap = argparse.ArgumentParser(
        description="Measure chimeric-read rate from a first-round SAM.")
    ap.add_argument("sam", nargs="+", help="SAM file(s), optionally .gz")
    ap.add_argument("--min-seg-frac", type=float, default=0.15,
                    help="minimum fraction of the read a segment must cover "
                         "to count as a real piece [0.15]")
    ap.add_argument("--min-total-cov", type=float, default=0.70,
                    help="segments together must explain at least this much "
                         "of the read [0.70]")
    ap.add_argument("--min-mapq", type=int, default=1,
                    help="minimum MAPQ per segment [1]")
    ap.add_argument("--contig-map", default=None,
                    help="optional TSV mapping reference name -> genome name, "
                         "needed if one genome spans several FASTA records")
    ap.add_argument("--unsorted", action="store_true",
                    help="buffer whole file by read name (use if alignment "
                         "records for a read are not adjacent)")
    ap.add_argument("--out", default="chimeric_reads.tsv",
                    help="per-read output table [chimeric_reads.tsv]")
    args = ap.parse_args()

    contig_map = load_contig_map(args.contig_map)

    total_reads = 0
    multi_record = 0
    counts = Counter()
    pairs = Counter()
    junction_hist = Counter()
    junction_vals = []
    lengths = {"all": [], "inter": []}

    rows = []
    for qname, alns in iter_reads(args.sam, args.unsorted):
        total_reads += 1
        for a in alns:
            a["genome"] = genome_of(a["rname"], contig_map)
        if alns:
            lengths["all"].append(max(a["full_len"] for a in alns))
        if len(alns) > 1:
            multi_record += 1
        res = classify_read(alns, args.min_seg_frac,
                            args.min_mapq, args.min_total_cov)
        if not res:
            counts["single_segment"] += 1
            continue
        counts[res["kind"]] += 1
        if res["kind"] == "inter_genome":
            pairs[res["pair"]] += 1
            lengths["inter"].append(res["read_len"])
            if res["junction_frac"] is not None:
                junction_vals.append(res["junction_frac"])
                junction_hist[round(res["junction_frac"] * 10) / 10] += 1
            rows.append(res)

    if total_reads == 0:
        print("No alignment records found. Are these SAM files?", file=sys.stderr)
        return 1

    inter = counts["inter_genome"]
    intra = counts["intra_genome"]

    print("=" * 72)
    print("CHIMERA COUNT FROM FIRST-ROUND ALIGNMENT")
    print("=" * 72)
    print(f"  reads with >=1 usable alignment : {total_reads:,}")
    print(f"  reads with >1 alignment record  : {multi_record:,}")
    print()
    print(f"  inter-genome chimeras           : {inter:,}"
          f"   ({inter/total_reads:.4%})   <-- compare to HERALD candidates")
    print(f"  intra-genome (same genome, far) : {intra:,}"
          f"   ({intra/total_reads:.4%})")
    print(f"  contiguous split (not chimeric) : {counts['contiguous_split']:,}")
    print(f"  single-segment reads            : {counts['single_segment']:,}")

    if inter:
        print()
        print("-" * 72)
        print("JUNCTION POSITION along the read (0.5 = midpoint)")
        print("-" * 72)
        for k in sorted(junction_hist):
            n = junction_hist[k]
            print(f"  {k:0.1f}  {n:6d}  {n/inter:6.1%}  "
                  + "#" * int(50 * n / inter))
        mid = sum(1 for v in junction_vals if 0.4 <= v <= 0.6)
        print(f"\n  within 0.4-0.6 of read length: {mid}/{inter} "
              f"({mid/inter:.1%})")

        print()
        print("-" * 72)
        print("READ LENGTH: chimeric vs all")
        print("-" * 72)
        def med(v):
            v = sorted(v)
            return v[len(v)//2] if v else 0
        print(f"  median length, all reads      : {med(lengths['all']):,} bp")
        print(f"  median length, inter-genome   : {med(lengths['inter']):,} bp")
        if med(lengths["all"]):
            print(f"  ratio                         : "
                  f"{med(lengths['inter'])/med(lengths['all']):.2f}x  "
                  f"(expect ~2x if two molecules were ligated)")

        print()
        print("-" * 72)
        print("TOP GENOME PAIRS")
        print("-" * 72)
        for p, c in pairs.most_common(15):
            a, b = p.split(",")
            print(f"  {c:6d}  {a[:34]:34s} / {b[:34]}")

        cols = ["read", "read_len", "n_segments", "kind", "genomes",
                "junction_frac", "junction_bp", "covered_frac", "strands"]
        with open(args.out, "w") as out:
            out.write("\t".join(cols) + "\n")
            for r in rows:
                out.write("\t".join(str(r[c]) for c in cols) + "\n")
        print(f"\nWrote {len(rows)} chimeric reads to {args.out}")
        print("\nCompare the inter-genome count above with HERALD "
              "candidates.\nSimilar numbers mean two independent methods "
              "agree the reads are chimeric.")
    return 0


if __name__ == "__main__":
    sys.exit(main())