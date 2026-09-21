# PacBio mock-community analysis (module 3)

Supporting analysis for the ATCC MSA-1003 module. These produce the counts that
were assembled into `HERALD_pacbio_figure_data.xlsx`.

| Script | Does |
|---|---|
| `analyze_candidates.py` | Classifies each HERALD candidate read: host genome, donor, architecture (AB / BA / ABA / complex), insert size in fragments. Reports candidate rates normalized by each organism's nominal abundance in the mix. Aims to separate reference-quality effects from coverage effects. Writes `candidates_annotated.tsv`. |
| `length.py` | Three subcommands. `table` builds a read-length table from FASTA/FASTQ; `plan` proposes length bins and `-s` values so fragments land above a target length; `join` adds read length to the candidate table and reports candidate rate by read-length decile. |
| `chimera_analysis.py` | Measures the chimeric-read rate directly from a first-round SAM using supplementary alignments, independently of HERALD. |

`analyze_candidates.py` carries the MSA-1003 nominal composition as a hardcoded table, including a flag for whether the reference genome used is the strain actually in the mix. That table is specific to this community and would need replacing for any other dataset.
