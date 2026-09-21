#!/usr/bin/env python3
"""
build_mock_genome.py — construct a mock HGT organism by inserting donor segments into an acceptor 
genome, and record ground truth in both coordinatesystems needed for scoring.

-----------------------------------------------------
An insertion point can be named in two different coordinate spaces: 
- a position in the *acceptor reference*: the unmodified genome
- a position in the *mock genome* that was built from it

The truth table therefore carries:
  ref_pos       insertion point in the ORIGINAL acceptor, 0-based
  ref_pos_vcf   the same point, 1-based, matching VCF convention
  mock_start / mock_end   the insert's span in the MOCK genome, 0-based
                          half-open, for read-level checks

Masking
-------
Pass the BED files produced by genome_overlap.py. Insertion sites are rejected if anything 
within --site-clearance bp is masked, because a breakpoint with ambiguous flanking sequence 
cannot be scored fairly. Donor source segments are rejected if they overlap the donor's mask,
because sequence the donor shares with other panel members can be attributed elsewhere.

Divergence
----------
Acceptor and donor lineages diverge independently, so divergence is applied to each separately: 
--snp-acceptor to host sequence, --snp-donor to insert sequence. SNP counts are exact (round(length * rate)) 
rather than drawn, so runs are reproducible from the seed alone.

Negative controls
-----------------
--insert-sizes "" produces a divergence-only genome with no insertions. Ground truth is 
written as an empty table with the same header so downstream scoring needs no special case.

Examples
--------
  # One HGT mock genome: five inserts, one of each size
  python build_mock_genome.py \
      --acceptor panel/CP000302.1.fa --donor panel/AP014956.1.fa \
      --acceptor-mask overlap/masks/CP000302.1.bed \
      --donor-mask    overlap/masks/AP014956.1.bed \
      --insert-sizes 500,1000,2000,4000,8000 \
      --sample pair1_rep1 --seed 101 --outdir mock/

  # Matching negative control
  python build_mock_genome.py \
      --acceptor panel/CP000302.1.fa --donor panel/AP014956.1.fa \
      --insert-sizes "" --sample pair1_neg --seed 201 --outdir mock/
"""

import argparse
import os
import random
import sys

TRUTH_COLS = [
    "sample", "event_id", "acceptor_acc", "donor_acc",
    "ref_pos", "ref_pos_vcf", "insert_len",
    "donor_start", "donor_end", "donor_strand",
    "mock_start", "mock_end",
    "snp_acceptor", "snp_donor", "seed",
]

COMP = str.maketrans("ACGTNacgtn", "TGCANtgcan")


# ── IO ────────────────────────────────────────────────────────────────────────

def read_single_fasta(path):
    """Return (record_id, sequence). If multi-record, take the longest and warn."""
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
        sys.stderr.write(
            f"WARNING: {path} has {len(recs)} records; using the longest "
            f"({recs[0][0]}, {len(recs[0][1]):,} bp) and ignoring the rest.\n")
    return recs[0][0], recs[0][1].upper()


def read_bed(path):
    """Return sorted merged [start, end) intervals, or [] if path is None."""
    if not path:
        return []
    iv = []
    with open(path) as fh:
        for line in fh:
            if not line.strip() or line.startswith(("#", "track")):
                continue
            f = line.split("\t")
            iv.append((int(f[1]), int(f[2])))
    if not iv:
        return []
    iv.sort()
    merged = [list(iv[0])]
    for s, e in iv[1:]:
        if s <= merged[-1][1]:
            merged[-1][1] = max(merged[-1][1], e)
        else:
            merged.append([s, e])
    return [tuple(x) for x in merged]


def write_fasta(path, name, seq, width=80):
    with open(path, "w") as fh:
        fh.write(f">{name}\n")
        for i in range(0, len(seq), width):
            fh.write(seq[i:i + width] + "\n")


# ── Interval helpers ──────────────────────────────────────────────────────────

def overlaps(intervals, start, end):
    """True if [start, end) intersects any interval. Linear; masks are small."""
    for s, e in intervals:
        if s >= end:
            return False
        if e > start:
            return True
    return False


def revcomp(seq):
    return seq.translate(COMP)[::-1]


# ── Selection ─────────────────────────────────────────────────────────────────

def pick_donor_segments(donor, mask, sizes, rng, max_tries=10000):
    """Choose one non-overlapping, unmasked source segment per requested size."""
    chosen, used = [], []
    L = len(donor)
    for size in sizes:
        if L < size + 1:
            sys.exit(f"Donor is {L:,} bp, too short for a {size:,} bp insert.")
        for _ in range(max_tries):
            s = rng.randint(0, L - size)
            e = s + size
            if overlaps(mask, s, e) or overlaps(used, s, e):
                continue
            if donor[s:e].count("N") > 0:
                continue
            used.append((s, e))
            used.sort()
            strand = rng.choice(["+", "-"])
            chosen.append((s, e, strand))
            break
        else:
            sys.exit(
                f"Could not place a {size:,} bp donor segment after {max_tries} "
                f"tries. The donor mask may cover too much of the genome; check "
                f"pct_masked in genomes.tsv.")
    return chosen


def pick_sites(acceptor_len, mask, n, rng, min_spacing, edge_buffer,
               clearance, acceptor, max_tries=100000):
    """Choose n insertion points with clear, unique, well-separated flanks."""
    sites = []
    lo, hi = edge_buffer, acceptor_len - edge_buffer
    if hi - lo < n * min_spacing:
        sys.exit(
            f"Acceptor is {acceptor_len:,} bp; cannot fit {n} sites at "
            f"{min_spacing:,} bp spacing inside a {edge_buffer:,} bp edge "
            f"buffer. Use a larger acceptor or fewer inserts per genome.")
    for _ in range(max_tries):
        if len(sites) == n:
            break
        p = rng.randint(lo, hi)
        if overlaps(mask, p - clearance, p + clearance):
            continue
        if any(abs(p - q) < min_spacing for q in sites):
            continue
        if acceptor[p - clearance:p + clearance].count("N") > 0:
            continue
        sites.append(p)
    if len(sites) != n:
        sys.exit(
            f"Only placed {len(sites)}/{n} insertion sites. Loosen "
            f"--min-spacing or --site-clearance, or pick a larger acceptor.")
    return sorted(sites)


# ── Divergence ────────────────────────────────────────────────────────────────

def apply_snps(seq_list, offset, length, rate, rng):
    """
    Substitute exactly round(length * rate) positions inside
    seq_list[offset : offset+length]. Mutates in place, returns count.
    """
    k = int(round(length * rate))
    if k <= 0 or length <= 0:
        return 0
    positions = rng.sample(range(offset, offset + length), min(k, length))
    bases = "ACGT"
    n = 0
    for i in positions:
        b = seq_list[i]
        if b not in bases:
            continue
        seq_list[i] = rng.choice([x for x in bases if x != b])
        n += 1
    return n


# ── Build ─────────────────────────────────────────────────────────────────────

def build(acceptor, donor, sites, segments, snp_a, snp_d, rng):
    """
    Assemble the mock genome left to right, tracking the running offset so mock
    coordinates stay correct as earlier insertions shift later sequence.
    Returns (mock_sequence, events, snp_counts).
    """
    order = sorted(zip(sites, segments), key=lambda x: x[0])
    pieces, events = [], []
    prev, shift = 0, 0
    host_spans, insert_spans = [], []
    cursor = 0

    for site, (ds, de, strand) in order:
        host_chunk = acceptor[prev:site]
        pieces.append(host_chunk)
        host_spans.append((cursor, len(host_chunk)))
        cursor += len(host_chunk)

        seg = donor[ds:de]
        if strand == "-":
            seg = revcomp(seg)
        pieces.append(seg)
        insert_spans.append((cursor, len(seg)))

        events.append({
            "ref_pos": site,
            "ref_pos_vcf": site + 1,
            "insert_len": de - ds,
            "donor_start": ds,
            "donor_end": de,
            "donor_strand": strand,
            "mock_start": cursor,
            "mock_end": cursor + len(seg),
        })
        cursor += len(seg)
        prev = site
        shift += len(seg)

    tail = acceptor[prev:]
    pieces.append(tail)
    host_spans.append((cursor, len(tail)))

    mock = list("".join(pieces))

    # verify placement before divergence is applied
    for ev in events:
        seg = donor[ev["donor_start"]:ev["donor_end"]]
        if ev["donor_strand"] == "-":
            seg = revcomp(seg)
        placed = "".join(mock[ev["mock_start"]:ev["mock_end"]])
        if placed != seg:
            sys.exit("Internal error: insert does not match its source segment.")

    n_a = sum(apply_snps(mock, off, ln, snp_a, rng) for off, ln in host_spans)
    n_d = sum(apply_snps(mock, off, ln, snp_d, rng) for off, ln in insert_spans)

    return "".join(mock), events, (n_a, n_d)


# ── Main ──────────────────────────────────────────────────────────────────────

def parse_sizes(s):
    s = s.strip()
    if not s:
        return []
    return [int(x) for x in s.split(",") if x.strip()]


def main():
    p = argparse.ArgumentParser(
        description="Build a mock HGT genome with ground truth in reference "
                    "and mock coordinates")
    p.add_argument("--acceptor", required=True)
    p.add_argument("--donor", required=True)
    p.add_argument("--acceptor-mask", default=None,
                   help="BED from genome_overlap.py masks/")
    p.add_argument("--donor-mask", default=None)
    p.add_argument("--insert-sizes", default="500,1000,2000,4000,8000",
                   help="Comma-separated bp. Empty string = negative control.")
    p.add_argument("--reps", type=int, default=1,
                   help="Copies of each size in THIS genome (default 1)")
    p.add_argument("--min-spacing", type=int, default=100000,
                   help="Minimum bp between insertion sites (default 100000)")
    p.add_argument("--edge-buffer", type=int, default=50000,
                   help="Keep sites this far from genome ends (default 50000)")
    p.add_argument("--site-clearance", type=int, default=2000,
                   help="Required unmasked window either side of a site (2000)")
    p.add_argument("--snp-acceptor", type=float, default=0.001,
                   help="Substitution rate applied to host sequence (0.001)")
    p.add_argument("--snp-donor", type=float, default=0.001,
                   help="Substitution rate applied to insert sequence (0.001)")
    p.add_argument("--sample", required=True, help="Sample name / output prefix")
    p.add_argument("--seed", type=int, required=True)
    p.add_argument("--outdir", required=True)
    args = p.parse_args()

    os.makedirs(args.outdir, exist_ok=True)
    rng = random.Random(args.seed)

    acc_id, acceptor = read_single_fasta(args.acceptor)
    don_id, donor = read_single_fasta(args.donor)
    acc_mask = read_bed(args.acceptor_mask)
    don_mask = read_bed(args.donor_mask)

    sizes = parse_sizes(args.insert_sizes) * args.reps

    if sizes:
        segments = pick_donor_segments(donor, don_mask, sizes, rng)
        sites = pick_sites(len(acceptor), acc_mask, len(sizes), rng,
                           args.min_spacing, args.edge_buffer,
                           args.site_clearance, acceptor)
        mock, events, (n_a, n_d) = build(
            acceptor, donor, sites, segments,
            args.snp_acceptor, args.snp_donor, rng)
    else:
        seq = list(acceptor)
        n_a = apply_snps(seq, 0, len(seq), args.snp_acceptor, rng)
        n_d = 0
        mock, events = "".join(seq), []

    expected = len(acceptor) + sum(e["insert_len"] for e in events)
    if len(mock) != expected:
        sys.exit(f"Internal error: mock is {len(mock):,} bp, expected {expected:,}")

    fasta_path = os.path.join(args.outdir, f"{args.sample}.fasta")
    write_fasta(fasta_path, f"{args.sample} acceptor={acc_id} donor={don_id}", mock)

    truth_path = os.path.join(args.outdir, f"{args.sample}_truth.tsv")
    with open(truth_path, "w") as fh:
        fh.write("\t".join(TRUTH_COLS) + "\n")
        for i, ev in enumerate(events, 1):
            row = {
                "sample": args.sample,
                "event_id": f"{args.sample}_e{i}",
                "acceptor_acc": acc_id,
                "donor_acc": don_id,
                "snp_acceptor": args.snp_acceptor,
                "snp_donor": args.snp_donor,
                "seed": args.seed,
                **ev,
            }
            fh.write("\t".join(str(row[c]) for c in TRUTH_COLS) + "\n")

    print(f"{args.sample}: acceptor {acc_id} ({len(acceptor):,} bp) + "
          f"{len(events)} insert(s) from {don_id} = {len(mock):,} bp")
    if events:
        print("  sizes: " + ", ".join(f"{e['insert_len']:,}" for e in events))
        print("  ref positions: " + ", ".join(f"{e['ref_pos']:,}" for e in events))
    print(f"  SNPs: {n_a:,} host, {n_d:,} insert")
    print(f"  masks: acceptor {len(acc_mask)} interval(s), donor {len(don_mask)}")
    print(f"  -> {fasta_path}\n  -> {truth_path}")


if __name__ == "__main__":
    main()