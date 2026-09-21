# Benchmark and simulation scripts

Scripts used to generate the test data for the HERALD manuscript. They are not
part of the HERALD tool and nothing in `src/` depends on them. In particular,
HERALD itself does not require Biopython; three scripts here do.

```sh
pip3 install -r requirements.txt
```

External tools: minimap2 (2.30-r1287 was used), and seqkit for the genome
subsetting described in Methods.

## Read-level modules (modules 1 and 2)

Reads are constructed directly rather than sampled from a genome, so every read
is a known quantity.

| Script | Produces |
|---|---|
| `simulate_fp_read.py` | Non-chimeric reads from each genome in a multi-FASTA. Baseline false-positive testing. |
| `simulate_ABA_hgt.py` | Chimeric reads with a centered donor insert (ABA), plus a ground-truth TSV. |
| `simulate_HGT_noncenter.py` | Chimeric reads with insert position parameterized by `--flank1_pct`, covering AB, ABA and BA layouts, or swept across the full range. Supersedes the above for position experiments. |

## Genome-scale module (module 4)

Run in this order:

1. `genome_overlap.py` — all-vs-all minimap2 over the genome panel. Emits
   per-genome exclusion masks (BED) and the pair table, marking regions unsafe
   to use as insertion sites or insert sources.
2. `build_mock_genome.py` — inserts donor segments into an acceptor genome at
   unmasked sites, applies divergence, and writes ground truth in both reference
   and mock coordinates. `--insert-sizes ""` produces the matched negative
   control.
3. `simulate_reads_from_genome.py` — samples reads from the mock genome at a
   specified coverage, both strands, recording each read's origin coordinate in
   the header.

## PacBio module (module 3)

`build_msa1003_reference.py` — assembles the ATCC MSA-1003 reference from
genomes downloaded from NCBI, one record per organism. Run
`--list-accessions` first to get the accession list.

## Panel selection

`phylogeny_check.py` scans the genome collection for same-species and same-genus
pairs above 3 Mb, which is how the acceptor/donor pairs spanning the
shared-sequence gradient were chosen. The input filename is hardcoded near the
top; edit it before running.

## Reproducibility

All read generation in the manuscript used seed 42. Note that the two error
models differ: `simulate_fp_read.py`, `simulate_ABA_hgt.py` and
`simulate_HGT_noncenter.py` draw each position independently against the rate,
so substitution counts vary between reads. `simulate_reads_from_genome.py` and
`build_mock_genome.py` substitute exactly `round(length * rate)` positions. Both
are reproducible from the seed; only the first matches the description in
Methods 2.5.

## Scoring

The code that compares HERALD output against the `*_truth.tsv` tables is in
`analysis/scoring/`.
