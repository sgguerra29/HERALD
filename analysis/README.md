# Analysis and figure scripts

These are the R scripts used to produce the figures and statistics. They are provided as-run, for transparency. They are not packaged for reuse and are not part of the HERALD tool: nothing in `src/`
depends on anything here.

## Requirements

R with: tidyverse, readxl, patchwork, gridExtra, grid, showtext, ggplot2.

`showtext` downloads the Lato typeface from Google Fonts at run time, so the
first run needs a network connection.

## Running

The data files are in this directory, so the scripts run as-is from here:

```sh
cd analysis
Rscript HERALD_fig1.R
Rscript HERALD_fig2.R
Rscript HERALD_pacbio.R
Rscript HERALD_module4.R
```

Figures are written to the same directory. To read and write elsewhere, set
`HERALD_DIR`; `HERALD_theme.R` and the data files must be there too.

For the module 4 figure as published, source the panel D replacement in the
same session rather than running `HERALD_module4.R` alone:

```r
source("HERALD_module4.R")
source("HERALD_paneldonor.R")
```

## Included data

| File | Used by |
|---|---|
| `HERALD_plottingdata.xlsx` | `HERALD_fig1.R`, `HERALD_fig2.R` |
| `HERALD_pacbio_figure_data.xlsx` | `HERALD_pacbio.R` |
| `panel_events.tsv`, `limitedref_events.tsv`, `sweep_events.tsv`, `panel_samples.tsv`, `overlap/pairs.tsv` | `HERALD_module4.R` |
| `all_strata.tsv` | `HERALD_paneldonor.R` |

These are scored results, not raw HERALD output. The code that produced them
from HERALD output and the simulator ground-truth tables is in `scoring/`.

`panel_samples.tsv` contains the 20 insert-bearing samples only. The eight
negative controls produced no flagged reads and therefore no output rows;
`HERALD_module4.R` prints a note about this. Their read counts are stated in the
manuscript text.

## Scripts and their inputs

| Script | Produces | Reads |
|---|---|---|
| `HERALD_theme.R` | Shared styling; sourced by the others, not run directly | — |
| `HERALD_fig1.R` | Module 1 figure, four panels, plus binomial GLMs | `HERALD.plottingdata.xlsx` sheets: false positive, spike rate, var insert length, var insert placement |
| `HERALD_fig2.R` | Module 2 parameter figure, four panels plus one supplementary, plus GLMs | `HERALD.plottingdata.xlsx` sheets: mmf, frag number |
| `HERALD_pacbio.R` | PacBio mock-community figure and Table 1 | `HERALD_pacbio_figure_data.xlsx` sheets: filter_ladder, genome_pairs, candidates_with_length, shared_blocks_staph, length_deciles, length_multiples |
| `HERALD_module4.R` | Module 4 figure, plus logistic regression, Fisher, chi-square, Clopper-Pearson intervals, and supplementary CSVs | `panel_events.tsv`, `limitedref_events.tsv`, `sweep_events.tsv`, `panel_samples.tsv`, `overlap/pairs.tsv` |
| `HERALD_paneldonor.R` | Replacement panel D for the module 4 figure | `all_strata.tsv`, plus the environment left by `HERALD_module4.R` |

`HERALD_paneldonor.R` is a patch, not a standalone script. Source it in the same
R session immediately after `HERALD_module4.R`; it checks for the objects it
needs and stops with a clear message if run on its own.

## Subdirectories

| Directory | Contents |
|---|---|
| `scoring/` | Produces the scored TSVs above from HERALD output and ground truth (module 4) |
| `pacbio/` | Candidate classification, read-length analysis, and an independent chimera cross-check (module 3) |
| `readlevel/` | The single-fragment correction described in Methods 2.7 (modules 1 and 2) |

Each has its own README.
