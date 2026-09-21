#!/usr/bin/env python3
"""
Analyze read/fragment alignment files to determine how many reads have
0, 1, or >1 fragments that aligned to a "foreign" genome (i.e. a genome
different from the one the read was simulated from), broken down by
read length (10000bp, 15000bp, 20000bp, 25000bp, etc.).

Expected line format (tab or whitespace separated):
    <read_id>_frag_<N>    <aligned_genome>    AS:i:<score>

Where <read_id> looks like:
    <TRUE_GENOME>_read<N>_<LENGTH>bp
e.g. CP061472.1_read8_10000bp

A fragment counts as "foreign" if <aligned_genome> != <TRUE_GENOME>.

Usage:
    python results_single_frag.py file1.txt file2.txt ...
    python results_single_frag.py *.txt
"""

import sys
import re
import glob
import csv
from collections import defaultdict

FRAG_RE = re.compile(r'^(.*)_frag_\d+$')
READ_RE = re.compile(r'^(.+?)_read\d+_(\d+)bp$')


def parse_file(path, reads):
    """
    Populate `reads`: dict mapping read_id -> {
        "true_genome": str,
        "length": int,
        "foreign_count": int,
        "total_frags": int,
    }
    """
    with open(path) as f:
        for line in f:
            line = line.strip()
            if not line:
                continue

            parts = line.split('\t')
            if len(parts) < 2:
                parts = line.split()
            if len(parts) < 2:
                continue

            frag_name, aligned_genome = parts[0], parts[1]

            m = FRAG_RE.match(frag_name)
            if not m:
                continue
            read_id = m.group(1)  # e.g. CP061472.1_read8_10000bp

            m2 = READ_RE.match(read_id)
            if not m2:
                continue
            true_genome, length = m2.group(1), int(m2.group(2))

            if read_id not in reads:
                reads[read_id] = {
                    "true_genome": true_genome,
                    "length": length,
                    "foreign_count": 0,
                    "total_frags": 0,
                }
            reads[read_id]["total_frags"] += 1
            if aligned_genome != true_genome:
                reads[read_id]["foreign_count"] += 1


def bucket_for(foreign_count):
    if foreign_count == 0:
        return "0"
    elif foreign_count == 1:
        return "1"
    else:
        return ">1"


def main():
    if len(sys.argv) < 2:
        print("Usage: python results_single_frag.py <file1.txt> [file2.txt ...]")
        print("       (wildcards like *.txt are fine too)")
        sys.exit(1)

    # Expand any wildcards the shell didn't already expand
    files = []
    for arg in sys.argv[1:]:
        matched = glob.glob(arg)
        files.extend(matched if matched else [arg])

    reads = {}
    for path in files:
        parse_file(path, reads)

    if not reads:
        print("No reads parsed. Check that your file(s) match the expected "
              "naming pattern <GENOME>_read<N>_<LENGTH>bp_frag_<N>.")
        sys.exit(1)

    overall = {"0": 0, "1": 0, ">1": 0}
    per_length = defaultdict(lambda: {"0": 0, "1": 0, ">1": 0})

    for read_id, info in reads.items():
        b = bucket_for(info["foreign_count"])
        overall[b] += 1
        per_length[info["length"]][b] += 1

    print(f"Files parsed: {len(files)}")
    print(f"Total reads analyzed: {len(reads)}\n")

    print("=== Overall ===")
    print(f"  0 foreign fragments   : {overall['0']}")
    print(f"  1 foreign fragment    : {overall['1']}")
    print(f"  >1 foreign fragments  : {overall['>1']}")
    print()

    print("=== Per read length ===")
    for length in sorted(per_length):
        b = per_length[length]
        total = sum(b.values())
        print(f"  {length}bp (n={total}):")
        print(f"      0 foreign  : {b['0']}")
        print(f"      1 foreign  : {b['1']}")
        print(f"      >1 foreign : {b['>1']}")
    print()

    out_csv = "fragment_alignment_summary.csv"
    with open(out_csv, "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(["read_id", "true_genome", "length_bp",
                          "total_frags", "foreign_frags", "bucket"])
        for read_id, info in sorted(reads.items()):
            b = bucket_for(info["foreign_count"])
            writer.writerow([read_id, info["true_genome"], info["length"],
                              info["total_frags"], info["foreign_count"], b])
    print(f"Per-read details written to {out_csv}")


if __name__ == "__main__":
    main()