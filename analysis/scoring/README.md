# Scoring (module 4)

These produce the scored TSVs in `analysis/`, which the R figure scripts read.

| `score_herald.py` | `panel_events.tsv`, `panel_samples.tsv`, and by re-running with a different `--arm`, `limitedref_events.tsv` and `sweep_events.tsv` | HERALD result files, `all_truth.tsv` from `build_mock_genome.py`, and the simulated read files |
| `donor_content.py` | `all_strata.tsv` | Truth intervals, read coordinates, and HERALD's flagged read IDs |
| `read_census.py` | `census_events.tsv`, `census_samples.tsv` | `all_truth.tsv` and the read files |

## How reads are classified

`score_herald.py` matches every flagged read back to the truth table using the
`pos=` field written into read headers by `simulate_reads_from_genome.py`, then
labels it:

* **true positive** — overlaps a true event by at least `--min-overlap` bp and  carries at least one fragment assigned to that event's donor genome
* **misattributed** — overlaps a true event but no fragment matches the true donor: a real mosaic read credited to the wrong genome
* **false positive** — overlaps no true event

An event counts as detected at `--min-reads` supporting reads. Negative-control
samples have no truth rows, so every flagged read there is a false positive by
construction.

`read_census.py` exists because recall is measured over reads but the truth
table records genome-level events. It reports four denominators of increasing
strictness (`reads_any`, `reads_ge200`, `reads_ge_frag`, `reads_flanked`) so a
detection failure can be told apart from an event no read could have carried.

Output filenames follow `--out-prefix`, which defaults to `panel`, so each arm
writes its own files without overwriting the others:

```sh
python3 score_herald.py --arm panel  --out-prefix panel      ...
python3 score_herald.py --arm panel  --out-prefix limitedref ...
python3 score_herald.py --arm sweep  --out-prefix sweep      --results "sweep/*.txt"
```

## Note on `all_strata.tsv`

The deposited file has `run`, `ci_low` and `ci_high` columns that
`donor_content.py`'s `by_sample` output does not write; the script emits
`sample, donor_bp_bin, reads, detected, recall`. The confidence intervals are
recomputed in R. 
