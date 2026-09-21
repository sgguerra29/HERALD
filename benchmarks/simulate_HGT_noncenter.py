"""
simulate_HGT_noncenter.py

Generates a mixed FASTA file of:
  - Normal (non-HGT) reads sampled across all genomes in a multi-genome FASTA
  - A specified number of synthetic chimeric HGT reads in AB or ABA layout

Excluded genomes (high false-positive offenders from false-positive testing):
  AP019314.1  Microcystis aeruginosa NIES-102
  CP012153.2  Francisella noatunensis subsp. orientalis FNO01

All synthetic HGT events are written to a ground-truth TSV so detection
results can be compared back to known truth.

------------------
By default the donor insert is centered in the read (equal host flanks on
each side, ABA layout).  Use --flank1_pct to slide the insert left or right:

  --flank1_pct 0.0    → [donor][host]              AB  layout (insert at left end)
  --flank1_pct 0.20   → [host 20%][donor][host 45%]  ABA layout (insert near left)
  --flank1_pct 0.325  → [host 32.5%][donor][host 32.5%]  ABA layout (centered, default)
  --flank1_pct 0.50   → [host 50%][donor][host 15%]  ABA layout (insert near right)
  --flank1_pct 0.65   → [host][donor]              BA  layout (insert at right end)
                        (where 0.65 = 1.0 - insert_pct for insert_pct=0.35)

Setting flank1_pct=0 produces an AB read (donor flush against the left end, host fills the right).  
Setting flank1_pct = 1 - insert_pct produces a BA read (host fills the left, donor flush against the right end).

The right-flank fraction is derived automatically:
  flank2_pct = 1.0 - flank1_pct - insert_pct

A validation error is raised if insert_pct + flank1_pct > 1.0 (no room for the right flank) or if flank1_pct < 0.

You can also sweep a range of positions across your HGT reads with --flank1_pct_sweep.  
When  flag is set, each HGT read is assigned a flank1_pct drawn uniformly at random from [0, 1 - insert_pct], 
covering the full range from AB through ABA to BA.

Usage examples:
  # Fixed centered insert / ABA (original behaviour)
  python simulate_HGT_noncenter.py --genomes genomes_500.fasta --n_hgt 3000 -o spiked

  # AB layout — donor flush against the left end
  python simulate_HGT_noncenter.py --genomes genomes_500.fasta --n_hgt 3000 \\
      --flank1_pct 0.0 -o spiked_ab

  # Insert starts 20% into the read (ABA, near-left)
  python simulate_HGT_noncenter.py --genomes genomes_500.fasta --n_hgt 3000 \\
      --flank1_pct 0.20 -o spiked_left

  # Sweep insert position uniformly from AB through ABA to BA
  python simulate_HGT_noncenter.py --genomes genomes_500.fasta --n_hgt 3000 \\
      --flank1_pct_sweep -o spiked_sweep

Outputs:
  <prefix>.fasta        All reads (normal + chimeric HGT), shuffled
  <prefix>_truth.tsv    Ground-truth record of every HGT read
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


def make_hgt_read(host_genome, donor_genome,
                  read_len, insert_len, flank1_len,
                  error_rate):
    """
    Chimeric HGT read: [host_flank1][donor_insert][host_flank2]

    Supports AB, ABA, and BA layouts depending on flank lengths:
      flank1_len = 0                      → AB  (donor at left end)
      0 < flank1_len < read_len-insert_len → ABA (donor in middle)
      flank1_len = read_len - insert_len  → BA  (donor at right end)

    When both flanks are non-zero, they are drawn as a single contiguous region from the host 
    genome and split at flank1_len.  When one flank is zero the host contributes only the remaining single segment.

    Parameters
    ----------
    host_genome : str
    donor_genome : str
    read_len    : int   total read length in bp
    insert_len  : int   donor insert length in bp
    flank1_len  : int   left host flank length in bp (0 = AB layout)
                        (flank2 fills the remainder; may also be 0 = BA layout)
    error_rate  : float per-base substitution rate

    Returns
    -------
    (read_seq, flank1_len, insert_len, flank2_len)  or  None if too short
    """
    flank2_len = read_len - insert_len - flank1_len
    if flank1_len < 0 or flank2_len < 0 or insert_len < 1:
        return None

    donor_seg = sample_segment(donor_genome, insert_len)
    if donor_seg is None:
        return None

    host_total = flank1_len + flank2_len
    if host_total > 0:
        host_seg_full = sample_segment(host_genome, host_total)
        if host_seg_full is None:
            return None
        host_seg1 = host_seg_full[:flank1_len]
        host_seg2 = host_seg_full[flank1_len:]
    else:
        host_seg1 = ""
        host_seg2 = ""

    read = introduce_sub(host_seg1 + donor_seg + host_seg2, error_rate)
    return read, flank1_len, insert_len, flank2_len


# ── Flank length resolver ─────────────────────────────────────────────────────

def resolve_flank1_len(read_len, insert_len, flank1_pct, sweep):
    """
    Return the left-flank length (int) for one HGT read.

    sweep=True  → sample flank1_pct uniformly from [0, 1 - insert_pct],
                  covering AB (0) through ABA to BA (1 - insert_pct).
    sweep=False → use the fixed flank1_pct (may be 0 for AB layout).
    """
    if sweep:
        insert_pct = insert_len / read_len
        flank1_pct = random.uniform(0.0, 1.0 - insert_pct)

    flank1_len = int(read_len * flank1_pct)
    # guard: flank1 cannot exceed what's left after the insert
    flank1_len = max(0, min(flank1_len, read_len - insert_len))
    return flank1_len


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
        description="Simulate a spiked long-read dataset with chimeric HGT events (AB, ABA, or BA)"
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

    # ── Insert-position arguments ─────────────────────────────────────────────
    pos_group = parser.add_mutually_exclusive_group()
    pos_group.add_argument(
        "--flank1_pct", type=float, default=None,
        help=(
            "Left host-flank size as a fraction of read length. "
            "0.0 = AB layout (donor at left end); "
            "1.0 - insert_pct = BA layout (donor at right end); "
            "values in between give ABA layout. "
            "The right flank fills the remainder: "
            "flank2 = 1 - flank1_pct - insert_pct. "
            "Default: centered ABA insert, i.e. (1 - insert_pct) / 2."
        )
    )
    pos_group.add_argument(
        "--flank1_pct_sweep", action="store_true",
        help=(
            "Randomly sweep the insert position for each HGT read, drawing "
            "flank1_pct uniformly from the full valid range. "
            "Mutually exclusive with --flank1_pct."
        )
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

    # Resolve fixed flank1_pct (used when not sweeping)
    if args.flank1_pct is not None:
        flank1_pct_fixed = args.flank1_pct
    else:
        flank1_pct_fixed = (1.0 - args.insert_pct) / 2  # centered default

    # Validate fixed position (skip when sweeping — every value is valid)
    if not args.flank1_pct_sweep:
        flank2_pct = 1.0 - flank1_pct_fixed - args.insert_pct
        if flank1_pct_fixed < 0 or flank2_pct < 0:
            sys.exit(
                f"Error: flank1_pct={flank1_pct_fixed:.4f} is out of range. "
                f"Must be between 0.0 and {1.0 - args.insert_pct:.4f} "
                f"(= 1 - insert_pct) inclusive."
            )

    excluded = EXCLUDED_ACCESSIONS | set(args.exclude)

    # ── Load genomes ──────────────────────────────────────────────────────────
    genomes = load_genomes(args.genomes, min_length=args.read_len, excluded=excluded)
    if len(genomes) < 2:
        sys.exit("Error: need at least 2 genomes after filtering.")

    accessions = list(genomes.keys())
    n_genomes  = len(accessions)

    flank1_len_fixed = resolve_flank1_len(
        args.read_len, insert_len, flank1_pct_fixed, sweep=False
    )
    flank2_len_fixed = args.read_len - insert_len - flank1_len_fixed

    print(f"\nParameters:")
    print(f"  Read length  : {args.read_len:,} bp")
    print(f"  Insert       : {insert_len:,} bp ({args.insert_pct*100:.1f}% of read)")
    if args.flank1_pct_sweep:
        print(f"  Host flanks  : SWEEP — flank1 drawn uniformly from "
              f"[0 bp (AB), {args.read_len - insert_len:,} bp (BA)] per read")
    else:
        print(f"  Host flanks  : flank1={flank1_len_fixed:,} bp "
              f"({flank1_pct_fixed*100:.1f}%)  |  "
              f"flank2={flank2_len_fixed:,} bp "
              f"({(1-flank1_pct_fixed-args.insert_pct)*100:.1f}%)")
    print(f"  Error rate   : {args.error_rate*100:.2f}%")
    print(f"  HGT reads    : {args.n_hgt:,}")
    print(f"  Normal reads : {args.reads_per_genome} per genome × "
          f"{n_genomes:,} genomes = ~{args.reads_per_genome * n_genomes:,}")
    print()

    all_reads  = []
    truth_rows = []

    # ── Generate normal reads ─────────────────────────────────────────────────
    normal_count   = 0
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

    # ── Generate chimeric HGT reads ───────────────────────────────────────────
    hgt_count   = 0
    hgt_skipped = 0

    for idx in range(args.n_hgt):
        host_acc, donor_acc = random.sample(accessions, 2)
        host_genome  = genomes[host_acc]
        donor_genome = genomes[donor_acc]

        # Determine this read's left-flank length
        flank1_len = resolve_flank1_len(
            args.read_len, insert_len,
            flank1_pct_fixed,
            sweep=args.flank1_pct_sweep,
        )

        result = make_hgt_read(
            host_genome, donor_genome,
            args.read_len, insert_len, flank1_len,
            args.error_rate,
        )

        if result is None:
            hgt_skipped += 1
            continue

        seq, h1, ins, h2 = result

        # Derive orientation label from flank lengths
        if h1 == 0:
            orientation = "AB"
        elif h2 == 0:
            orientation = "BA"
        else:
            orientation = "ABA"

        header = (
            f"HGT_{orientation}_host-{host_acc}_donor-{donor_acc}"
            f"_insert{insert_len}bp_flank1-{h1}bp_read{idx}"
        )
        all_reads.append((header, seq))

        truth_rows.append({
            "read_id":      header,
            "orientation":  orientation,
            "host_acc":     host_acc,
            "donor_acc":    donor_acc,
            "read_len":     args.read_len,
            "insert_len":   ins,
            "insert_pct":   round(args.insert_pct, 4),
            "flank1_len":   h1,
            "flank2_len":   h2,
            "flank1_pct":   round(h1 / args.read_len, 4),
            "donor_start":  h1,
            "donor_end":    h1 + ins,
            "error_rate":   args.error_rate,
        })
        hgt_count += 1

    print(f"Generated {hgt_count:,} chimeric HGT reads "
          f"({hgt_skipped} skipped — genome pair too short)")

    # ── Shuffle and write FASTA ───────────────────────────────────────────────
    random.shuffle(all_reads)

    fasta_path = f"{args.output_prefix}.fasta"
    with open(fasta_path, "w") as fh:
        for header, seq in all_reads:
            fh.write(f">{header}\n{seq}\n")

    total     = len(all_reads)
    spike_pct = hgt_count / total * 100 if total else 0
    print(f"\nWrote {total:,} reads ({normal_count:,} normal + {hgt_count:,} HGT, "
          f"{spike_pct:.2f}% spike rate) → {fasta_path}")

    # ── Write ground-truth TSV ────────────────────────────────────────────────
    truth_path = f"{args.output_prefix}_truth.tsv"
    tsv_cols = [
        "read_id", "orientation", "host_acc", "donor_acc",
        "read_len", "insert_len", "insert_pct",
        "flank1_len", "flank2_len", "flank1_pct",
        "donor_start", "donor_end", "error_rate",
    ]
    with open(truth_path, "w") as fh:
        fh.write("\t".join(tsv_cols) + "\n")
        for row in truth_rows:
            fh.write("\t".join(str(row[c]) for c in tsv_cols) + "\n")

    print(f"Wrote {len(truth_rows):,} ground-truth HGT records → {truth_path}")


if __name__ == "__main__":
    main()