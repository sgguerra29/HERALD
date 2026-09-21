#!/usr/bin/env python3
"""
genome_overlap.py — measure pairwise sequence sharing across a genome panel, select donor/acceptor 
pairs for HGT-detection benchmarking, and emit exclusion masks marking regions that are not 
safe to use as insertion sites or as insert source material.


Usage
-----
  # One shot: run minimap2 and analyse (needs minimap2 on PATH)
  python genome_overlap.py --fasta panel_100.fasta --outdir overlap/ \
      --run-minimap2 --threads 16

  # Or reuse an existing PAF
  minimap2 -c -x asm20 -t 16 --secondary=no panel_100.fasta panel_100.fasta > panel.paf
  python genome_overlap.py --fasta panel_100.fasta --paf panel.paf --outdir overlap/

Multi-contig genomes
--------------------
By default each FASTA record is treated as one genome (for a panel of complete single-chromosome 
assemblies). If your genomes have multiple contigs, supply --map, a two-column TSV:
 <contig_id>\t<genome_id>

Outputs (in --outdir)
---------------------
  pairs.tsv            every genome pair with any detected sharing
  genomes.tsv          per-genome summary incl. a "promiscuity" score
  suggested_pairs.tsv  candidate acceptor/donor pairs binned by sharing stratum
  masks/<genome>.bed   regions to avoid (shared with any other genome, or
                       self-repeated within this genome)
"""

import argparse
import os
import subprocess
import sys
from collections import defaultdict


# ── FASTA ─────────────────────────────────────────────────────────────────────

def read_fasta_lengths(path):
    """Return {record_id: length}. Streams; never holds a sequence in memory."""
    lengths = {}
    cur, n = None, 0
    with open(path) as fh:
        for line in fh:
            if line.startswith(">"):
                if cur is not None:
                    lengths[cur] = n
                cur = line[1:].split()[0]
                n = 0
            else:
                n += len(line.strip())
    if cur is not None:
        lengths[cur] = n
    return lengths


def read_contig_map(path):
    m = {}
    with open(path) as fh:
        for line in fh:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            parts = line.split("\t")
            if len(parts) < 2:
                sys.exit(f"--map line is not two columns: {line!r}")
            m[parts[0]] = parts[1]
    return m


# ── Intervals ─────────────────────────────────────────────────────────────────

def merge_intervals(intervals):
    """Merge a list of [start, end) tuples. Returns sorted, disjoint list."""
    if not intervals:
        return []
    intervals = sorted(intervals)
    merged = [list(intervals[0])]
    for s, e in intervals[1:]:
        if s <= merged[-1][1]:
            if e > merged[-1][1]:
                merged[-1][1] = e
        else:
            merged.append([s, e])
    return [tuple(x) for x in merged]


def total_bp(intervals):
    return sum(e - s for s, e in intervals)


# ── minimap2 ──────────────────────────────────────────────────────────────────

def run_minimap2(fasta, paf_out, preset, threads, extra=None):
    """
    -X puts minimap2 in all-vs-all mode: it drops the trivial self-diagonal but retains secondary alignments. 
    Do not add --secondary=no. In a self alignment the full-length diagonal is the primary hit and every repeat copy
    is a secondary, so suppressing secondaries silently deletes intra-genome repeats and most cross-genome homology.
    """
    cmd = ["minimap2", "-c", "-x", preset, "-X", "-t", str(threads)]
    if extra:
        cmd += extra.split()
    cmd += [fasta, fasta]
    sys.stderr.write("Running: " + " ".join(cmd) + "\n")
    with open(paf_out, "w") as fh:
        rc = subprocess.call(cmd, stdout=fh)
    if rc != 0:
        sys.exit(f"minimap2 exited with status {rc}")
    return paf_out


# ── PAF parsing ───────────────────────────────────────────────────────────────

def is_self_alignment(qname, tname, strand, qs, qe, ts, te):
    """
    A record is the trivial diagonal (a contig aligned to itself at the same  coordinates) if names match, strand 
    is forward, and the query and target spans overlap by most of their length. Genuine internal repeats also 
    haveqname == tname but sit OFF the diagonal, and we want to keep those.
    """
    if qname != tname or strand != "+":
        return False
    ov = min(qe, te) - max(qs, ts)
    if ov <= 0:
        return False
    shorter = min(qe - qs, te - ts)
    return shorter > 0 and ov / shorter > 0.9


def parse_paf(paf_path, contig2genome, min_identity, min_block):
    """
    Returns
      inter[(gA, gB)] -> {'a': [intervals on gA], 'b': [intervals on gB],
                          'nmatch': int, 'alen': int}
        with gA < gB lexicographically
      selfrep[g] -> [intervals] regions of g that are repeated elsewhere in g
      seen_contigs -> set
    """
    inter = defaultdict(lambda: {"a": [], "b": [], "nmatch": 0, "alen": 0})
    selfrep = defaultdict(list)
    seen = set()
    kept = dropped_id = dropped_len = diagonal = 0

    with open(paf_path) as fh:
        for line in fh:
            if not line.strip():
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 12:
                continue
            qname, tname, strand = f[0], f[5], f[4]
            qs, qe = int(f[2]), int(f[3])
            ts, te = int(f[7]), int(f[8])
            nmatch, alen = int(f[9]), int(f[10])

            seen.add(qname)
            seen.add(tname)

            if alen < min_block:
                dropped_len += 1
                continue
            if alen == 0 or nmatch / alen < min_identity:
                dropped_id += 1
                continue
            if is_self_alignment(qname, tname, strand, qs, qe, ts, te):
                diagonal += 1
                continue

            gq = contig2genome.get(qname, qname)
            gt = contig2genome.get(tname, tname)

            if gq == gt:
                # internal repeat within one genome
                selfrep[gq].append((qs, qe))
                selfrep[gq].append((ts, te))
            else:
                key = (gq, gt) if gq < gt else (gt, gq)
                rec = inter[key]
                if key[0] == gq:
                    rec["a"].append((qs, qe))
                    rec["b"].append((ts, te))
                else:
                    rec["a"].append((ts, te))
                    rec["b"].append((qs, qe))
                rec["nmatch"] += nmatch
                rec["alen"] += alen
            kept += 1

    sys.stderr.write(
        f"PAF: kept {kept:,} alignments; dropped {dropped_len:,} short, "
        f"{dropped_id:,} low-identity, {diagonal:,} self-diagonal\n"
    )
    return inter, selfrep, seen


# ── Analysis ──────────────────────────────────────────────────────────────────

def block_stats(intervals):
    lens = [e - s for s, e in intervals]
    return {
        "n_blocks": len(lens),
        "max_block": max(lens) if lens else 0,
        "n_ge500": sum(1 for x in lens if x >= 500),
        "n_ge1000": sum(1 for x in lens if x >= 1000),
        "n_ge5000": sum(1 for x in lens if x >= 5000),
    }


def analyse(inter, selfrep, glen):
    """Collapse raw intervals into per-pair rows and per-genome masks."""
    pair_rows = []
    mask = defaultdict(list)          # genome -> intervals to avoid
    shared_partners = defaultdict(set)
    shared_bp = defaultdict(int)

    for (ga, gb), rec in inter.items():
        ia = merge_intervals(rec["a"])
        ib = merge_intervals(rec["b"])
        bp_a, bp_b = total_bp(ia), total_bp(ib)
        st = block_stats(ia)
        ident = rec["nmatch"] / rec["alen"] if rec["alen"] else 0.0

        la, lb = glen.get(ga, 0), glen.get(gb, 0)
        pair_rows.append({
            "genome_a": ga,
            "genome_b": gb,
            "len_a": la,
            "len_b": lb,
            "shared_bp_a": bp_a,
            "shared_bp_b": bp_b,
            "pct_a": round(100.0 * bp_a / la, 4) if la else 0.0,
            "pct_b": round(100.0 * bp_b / lb, 4) if lb else 0.0,
            "mean_identity": round(ident, 4),
            **st,
        })

        mask[ga].extend(ia)
        mask[gb].extend(ib)
        shared_partners[ga].add(gb)
        shared_partners[gb].add(ga)
        shared_bp[ga] += bp_a
        shared_bp[gb] += bp_b

    # internal repeats also go into the mask
    for g, iv in selfrep.items():
        mask[g].extend(iv)

    mask = {g: merge_intervals(iv) for g, iv in mask.items()}

    genome_rows = []
    for g, L in sorted(glen.items()):
        m = mask.get(g, [])
        rep = merge_intervals(selfrep.get(g, []))
        genome_rows.append({
            "genome": g,
            "length": L,
            "n_partners": len(shared_partners.get(g, ())),
            "shared_bp_total": shared_bp.get(g, 0),
            "self_repeat_bp": total_bp(rep),
            "masked_bp": total_bp(m),
            "pct_masked": round(100.0 * total_bp(m) / L, 4) if L else 0.0,
        })

    pair_rows.sort(key=lambda r: -r["shared_bp_a"])
    return pair_rows, genome_rows, mask


def suggest_pairs(pair_rows, genome_rows, n_per_stratum, max_pct_masked):
    """
    Pick candidate acceptor/donor pairs spanning a range of sharing, excluding genomes that are broadly similar to
    the rest of the panel (make  donor-attribution problem ambiguous regardless of partner).
    """
    clean = {r["genome"] for r in genome_rows
             if r["pct_masked"] <= max_pct_masked}

    strata = [
        ("low",  0,      1_000),
        ("mid",  1_000,  10_000),
        ("high", 10_000, float("inf")),
    ]
    out = []
    for name, lo, hi in strata:
        picks = [r for r in pair_rows
                 if lo <= r["shared_bp_a"] < hi
                 and r["genome_a"] in clean and r["genome_b"] in clean]
        picks.sort(key=lambda r: -r["shared_bp_a"])
        # spread across the stratum rather than taking the top n
        if len(picks) > n_per_stratum:
            step = len(picks) / n_per_stratum
            picks = [picks[int(i * step)] for i in range(n_per_stratum)]
        for r in picks:
            row = dict(r)
            row["stratum"] = name
            out.append(row)
    return out


# ── Output ────────────────────────────────────────────────────────────────────

def write_tsv(path, rows, cols):
    with open(path, "w") as fh:
        fh.write("\t".join(cols) + "\n")
        for r in rows:
            fh.write("\t".join(str(r.get(c, "")) for c in cols) + "\n")


def write_masks(outdir, mask):
    d = os.path.join(outdir, "masks")
    os.makedirs(d, exist_ok=True)
    for g, iv in mask.items():
        safe = g.replace("/", "_").replace("|", "_")
        with open(os.path.join(d, f"{safe}.bed"), "w") as fh:
            for s, e in iv:
                fh.write(f"{g}\t{s}\t{e}\n")
    return d


# ── Main ──────────────────────────────────────────────────────────────────────

def main():
    p = argparse.ArgumentParser(
        description="Pairwise genome sharing analysis for HGT benchmark design"
    )
    p.add_argument("--fasta", required=True,
                   help="Multi-FASTA panel (one record per genome unless --map)")
    p.add_argument("--paf", default=None,
                   help="Precomputed all-vs-all PAF (skips minimap2)")
    p.add_argument("--run-minimap2", action="store_true",
                   help="Run minimap2 all-vs-all on --fasta")
    p.add_argument("--preset", default="asm20",
                   help="minimap2 preset (default asm20, ~5%% divergence)")
    p.add_argument("--mm2-extra", default=None,
                   help="Extra minimap2 flags, e.g. '-w 5' for more sensitivity "
                        "to short shared blocks. Never pass --secondary=no.")
    p.add_argument("--threads", type=int, default=8)
    p.add_argument("--map", default=None,
                   help="TSV contig_id<TAB>genome_id for multi-contig genomes")
    p.add_argument("--min-identity", type=float, default=0.95,
                   help="Minimum alignment identity to count as shared (0.95)")
    p.add_argument("--min-block", type=int, default=200,
                   help="Minimum alignment length in bp to count (200)")
    p.add_argument("--max-pct-masked", type=float, default=5.0,
                   help="Exclude genomes with more than this %% masked (5.0)")
    p.add_argument("--n-per-stratum", type=int, default=4,
                   help="Suggested pairs per sharing stratum (4)")
    p.add_argument("--outdir", required=True)
    args = p.parse_args()

    os.makedirs(args.outdir, exist_ok=True)

    contig_len = read_fasta_lengths(args.fasta)
    if not contig_len:
        sys.exit(f"No FASTA records found in {args.fasta}")
    contig2genome = read_contig_map(args.map) if args.map else {}

    glen = defaultdict(int)
    for c, L in contig_len.items():
        glen[contig2genome.get(c, c)] += L
    glen = dict(glen)
    sys.stderr.write(f"Panel: {len(contig_len):,} records → {len(glen):,} genomes, "
                     f"{sum(glen.values()):,} bp total\n")

    if args.paf:
        paf = args.paf
    elif args.run_minimap2:
        paf = run_minimap2(args.fasta,
                           os.path.join(args.outdir, "allvall.paf"),
                           args.preset, args.threads, args.mm2_extra)
    else:
        sys.exit("Supply either --paf or --run-minimap2")

    inter, selfrep, _ = parse_paf(paf, contig2genome,
                                  args.min_identity, args.min_block)
    pair_rows, genome_rows, mask = analyse(inter, selfrep, glen)

    pair_cols = ["genome_a", "genome_b", "len_a", "len_b",
                 "shared_bp_a", "shared_bp_b", "pct_a", "pct_b",
                 "mean_identity", "n_blocks", "max_block",
                 "n_ge500", "n_ge1000", "n_ge5000"]
    genome_cols = ["genome", "length", "n_partners", "shared_bp_total",
                   "self_repeat_bp", "masked_bp", "pct_masked"]

    write_tsv(os.path.join(args.outdir, "pairs.tsv"), pair_rows, pair_cols)
    write_tsv(os.path.join(args.outdir, "genomes.tsv"), genome_rows, genome_cols)

    sugg = suggest_pairs(pair_rows, genome_rows,
                         args.n_per_stratum, args.max_pct_masked)
    write_tsv(os.path.join(args.outdir, "suggested_pairs.tsv"),
              sugg, ["stratum"] + pair_cols)

    maskdir = write_masks(args.outdir, mask)

    n_sharing = len(pair_rows)
    n_possible = len(glen) * (len(glen) - 1) // 2
    sys.stderr.write(
        f"\n{n_sharing:,} of {n_possible:,} possible pairs share "
        f"\u2265{args.min_block} bp at \u2265{args.min_identity:.0%} identity\n"
        f"Wrote pairs.tsv, genomes.tsv, suggested_pairs.tsv and "
        f"{len(mask):,} mask files to {maskdir}\n"
    )

    n_with_rep = sum(1 for r in genome_rows if r["self_repeat_bp"] > 0)
    frac = n_with_rep / len(genome_rows) if genome_rows else 0
    if frac < 0.5:
        sys.stderr.write(
            f"\n*** ACCEPTANCE TEST FAILED ***\n"
            f"Only {n_with_rep}/{len(genome_rows)} genomes show any internal repeat.\n"
            f"Bacterial genomes carry 2-7 near-identical rRNA operons, so this\n"
            f"should be near 100%. The alignment step is losing secondary\n"
            f"alignments. Check that --secondary=no is NOT in the minimap2 call\n"
            f"and that -X is present.\n")
    else:
        med = sorted(r["self_repeat_bp"] for r in genome_rows)[len(genome_rows)//2]
        sys.stderr.write(
            f"\nAcceptance test passed: {n_with_rep}/{len(genome_rows)} genomes "
            f"show internal repeats (median {med:,} bp).\n")

    if sugg:
        sys.stderr.write("\nSuggested pairs:\n")
        sys.stderr.write(f"{'stratum':8} {'acceptor':22} {'donor':22} "
                         f"{'shared_bp':>10} {'max_blk':>8} {'ident':>6}\n")
        for r in sugg:
            sys.stderr.write(
                f"{r['stratum']:8} {r['genome_a'][:22]:22} {r['genome_b'][:22]:22} "
                f"{r['shared_bp_a']:>10,} {r['max_block']:>8,} "
                f"{r['mean_identity']:>6.3f}\n"
            )


if __name__ == "__main__":
    main()