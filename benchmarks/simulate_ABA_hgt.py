"""
simulate_ABA_hgt.py

Generates a mixed FASTA file of:
  - Normal (non-HGT) reads sampled across all genomes in a multi-genome FASTA
  - A specified number of synthetic ABA chimeric HGT reads, where a donor
    insert is flanked by host sequence on both sides (contiguous split)

Excluded genomes (high false-positive offenders from false-positive testing):
  AP019314.1   Microcystis aeruginosa NIES-102
  CP012153.2  Francisella noatunensis subsp. orientalis FNO01

All synthetic HGT events are written to a ground-truth TSV so detection
results can be compared back to known truth.

Usage example:
  python simulate_ABA_hgt.py \\
      --genomes genomes_500.fasta \\
      --read_len 20000 \\
      --insert_pct 0.35 \\
      --n_hgt 3000 \\
      --reads_per_genome 20 \\
      --error_rate 0.005 \\
      -o spiked_aba \\
      --seed 42

Outputs:
  spiked_aba.fasta        All reads (normal + ABA HGT), shuffled
  spiked_aba_truth.tsv    Ground-truth record of every HGT read
"""

import argparse
import random
import sys
from Bio import SeqIO

# ── Constants ─────────────────────────────────────────────────────────────────

BASES = ["A", "C", "G", "T"]
EXCLUDED_ACCESSIONS = {"AP019314.1", "CP012153.2"}


# ── Sequence utilities ────────────────────────────────────────────────────────

def mutate_base(b):
    return random.choice([x for x in BASES if x != b])


def introduce_sub(seq, error_rate):
    if error_rate == 0:
        return seq
    seq = list(seq)
    for i, b in enumerate(seq):
        if random.random() < error_rate:
            seq[i] = mutate_base(b)
    return "".join(seq)


def sample_segment(genome, length):
    """Draw a random contiguous segment of `length` bp from `genome`."""
    if len(genome) < length:
        return None
    start = random.randint(0, len(genome) - length)
    return genome[start:start + length]


# ── Read simulators ───────────────────────────────────────────────────────────

def make_normal_read(genome, read_len, error_rate):
    """A single non-chimeric read from one genome."""
    seg = sample_segment(genome, read_len)
    if seg is None:
        return None
    return introduce_sub(seg, error_rate)


def make_hgt_read_ABA(host_genome, donor_genome, read_len, insert_len, error_rate):
    """
    ABA chimera: [host_seg1][donor_insert][host_seg2]

    Host segments are drawn as a single contiguous region from the host genome
    and split at the midpoint — simulating a real insertion site where the
    donor sequence is inserted between two adjacent host regions.

    Returns (read_seq, host1_len, insert_len, host2_len) or None if too short.
    """
    host_total = read_len - insert_len
    host1_len = host_total // 2
    host2_len = host_total - host1_len  # absorbs odd remainder

    # draw one contiguous host region and split it
    host_seg_full = sample_segment(host_genome, host_total)
    donor_seg = sample_segment(donor_genome, insert_len)

    if host_seg_full is None or donor_seg is None:
        return None

    host_seg1 = host_seg_full[:host1_len]
    host_seg2 = host_seg_full[host1_len:]

    read = introduce_sub(host_seg1 + donor_seg + host_seg2, error_rate)
    return read, host1_len, insert_len, host2_len


# ── Genome loading ────────────────────────────────────────────────────────────

def load_genomes(fasta_path, min_length, excluded=EXCLUDED_ACCESSIONS):
    """
    Parse a multi-record FASTA. Returns a dict {accession: sequence_str}.
    Skips excluded accessions and genomes shorter than min_length.
    """
    genomes = {}
    skipped_excluded = 0
    skipped_short = 0

    for record in SeqIO.parse(fasta_path, "fasta"):
        acc = record.id
        if acc in excluded:
            skipped_excluded += 1
            continue
        seq = str(record.seq).upper()
        if len(seq) < min_length:
            skipped_short += 1
            continue
        genomes[acc] = seq

    print(f"Loaded {len(genomes):,} genomes from {fasta_path}")
    if skipped_excluded:
        print(f"  Excluded {skipped_excluded} high-FP genome(s): "
              f"{', '.join(sorted(excluded))}")
    if skipped_short:
        print(f"  Skipped {skipped_short} genome(s) shorter than {min_length:,} bp")
    return genomes


# ── Main ──────────────────────────────────────────────────────────────────────

def main():
    parser = argparse.ArgumentParser(
        description="Simulate a spiked long-read dataset with ABA HGT events"
    )
    parser.add_argument(
        "--genomes", required=True,
        help="Multi-genome FASTA (e.g. your 500-genome file)"
    )
    parser.add_argument(
        "-l", "--read_len", type=int, default=20000,
        help="Total read length in bp (default: 20000)"
    )
    parser.add_argument(
        "--insert_pct", type=float, default=0.35,
        help="Donor insert size as fraction of read length (default: 0.35 = 35%%)"
    )
    parser.add_argument(
        "--n_hgt", type=int, default=3000,
        help="Number of synthetic ABA HGT reads to spike in (default: 3000)"
    )
    parser.add_argument(
        "--reads_per_genome", type=int, default=20,
        help="Number of normal reads per genome (default: 20)"
    )
    parser.add_argument(
        "-e", "--error_rate", type=float, default=0.005,
        help="Substitution error rate, e.g. 0.005 = 0.5%% (default: 0.005)"
    )
    parser.add_argument(
        "-o", "--output_prefix", required=True,
        help="Output prefix; produces <prefix>.fasta and <prefix>_truth.tsv"
    )
    parser.add_argument(
        "-s", "--seed", type=int, default=42,
        help="Random seed for reproducibility (default: 42)"
    )
    parser.add_argument(
        "--exclude", nargs="*", default=[],
        help="Additional accessions to exclude (space-separated)"
    )

    args = parser.parse_args()
    random.seed(args.seed)

    insert_len = int(args.read_len * args.insert_pct)
    if insert_len < 1:
        sys.exit("Error: insert_pct too small — insert_len rounds to 0 bp.")
    if insert_len >= args.read_len:
        sys.exit("Error: insert_pct must be < 1.0.")

    excluded = EXCLUDED_ACCESSIONS | set(args.exclude)

    # ── Load genomes ──────────────────────────────────────────────────────────
    genomes = load_genomes(args.genomes, min_length=args.read_len, excluded=excluded)
    if len(genomes) < 2:
        sys.exit("Error: need at least 2 genomes after filtering.")

    accessions = list(genomes.keys())
    n_genomes = len(accessions)

    print(f"\nParameters:")
    print(f"  Read length  : {args.read_len:,} bp")
    print(f"  Insert       : {insert_len:,} bp ({args.insert_pct*100:.1f}% of read)")
    print(f"  Host flanks  : {(args.read_len - insert_len)//2:,} bp each side")
    print(f"  Error rate   : {args.error_rate*100:.2f}%")
    print(f"  HGT reads    : {args.n_hgt:,}")
    print(f"  Normal reads : {args.reads_per_genome} per genome × "
          f"{n_genomes:,} genomes = ~{args.reads_per_genome * n_genomes:,}")
    print()

    all_reads = []
    truth_rows = []

    # ── Generate normal reads ─────────────────────────────────────────────────
    normal_count = 0
    skipped_normal = 0
    for acc in accessions:
        genome = genomes[acc]
        for i in range(args.reads_per_genome):
            seq = make_normal_read(genome, args.read_len, args.error_rate)
            if seq is None:
                skipped_normal += 1
                continue
            header = f"NORMAL_{acc}_read{i}"
            all_reads.append((header, seq))
            normal_count += 1

    print(f"Generated {normal_count:,} normal reads "
          f"({skipped_normal} skipped — genome too short)")

    # ── Generate ABA HGT reads ────────────────────────────────────────────────
    hgt_count = 0
    hgt_skipped = 0

    for idx in range(args.n_hgt):
        host_acc, donor_acc = random.sample(accessions, 2)
        host_genome = genomes[host_acc]
        donor_genome = genomes[donor_acc]

        result = make_hgt_read_ABA(
            host_genome, donor_genome, args.read_len, insert_len, args.error_rate
        )

        if result is None:
            hgt_skipped += 1
            continue

        seq, host1_len, insert_len_out, host2_len = result

        header = (
            f"HGT_ABA_host-{host_acc}_donor-{donor_acc}"
            f"_insert{insert_len}bp_read{idx}"
        )
        all_reads.append((header, seq))

        truth_rows.append({
            "read_id":      header,
            "orientation":  "ABA",
            "host_acc":     host_acc,
            "donor_acc":    donor_acc,
            "read_len":     args.read_len,
            "insert_len":   insert_len,
            "insert_pct":   round(args.insert_pct, 4),
            "donor_start":  host1_len,
            "donor_end":    host1_len + insert_len,
            "error_rate":   args.error_rate,
        })
        hgt_count += 1

    print(f"Generated {hgt_count:,} ABA HGT reads "
          f"({hgt_skipped} skipped — genome pair too short)")

    # ── Shuffle and write FASTA ───────────────────────────────────────────────
    random.shuffle(all_reads)

    fasta_path = f"{args.output_prefix}.fasta"
    with open(fasta_path, "w") as fh:
        for header, seq in all_reads:
            fh.write(f">{header}\n{seq}\n")

    total = len(all_reads)
    spike_pct = hgt_count / total * 100 if total else 0
    print(f"\nWrote {total:,} reads ({normal_count:,} normal + {hgt_count:,} HGT, "
          f"{spike_pct:.2f}% spike rate) → {fasta_path}")

    # ── Write ground-truth TSV ────────────────────────────────────────────────
    truth_path = f"{args.output_prefix}_truth.tsv"
    tsv_cols = [
        "read_id", "orientation", "host_acc", "donor_acc",
        "read_len", "insert_len", "insert_pct",
        "donor_start", "donor_end", "error_rate",
    ]
    with open(truth_path, "w") as fh:
        fh.write("\t".join(tsv_cols) + "\n")
        for row in truth_rows:
            fh.write("\t".join(str(row[c]) for c in tsv_cols) + "\n")

    print(f"Wrote {len(truth_rows):,} ground-truth HGT records → {truth_path}")


if __name__ == "__main__":
    main()