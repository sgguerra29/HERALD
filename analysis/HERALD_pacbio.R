#!/usr/bin/env Rscript
# =============================================================================
# HERALD — PacBio mock community (ATCC MSA-1003, SRR11606871)
# =============================================================================

library(tidyverse)
library(readxl)
library(patchwork)
library(gridExtra)
library(grid)

# Directory holding the input data files and receiving the output figures.
# Defaults to the current working directory; override without editing this file:
#     HERALD_DIR=/path/to/data Rscript HERALD_pacbio.R
outdir <- Sys.getenv("HERALD_DIR", unset = getwd())
source(file.path(outdir, "HERALD_theme.R"))

# fail early if stale theme file is on disk
if (!"mult" %in% names(formals(herald_table_theme))) {
  stop("HERALD_theme.R on disk is out of date: herald_table_theme() has no ",
       "`mult` argument. Replace the file, then re-source it.", call. = FALSE)
}

W <- 13; H <- 6.5
SC <- herald_sc(W)

path <- file.path(outdir, "HERALD_pacbio_figure_data.xlsx")

TOTAL_READS <- 5646150
MEDIAN_READ <- 10061

BOTH_LAB  <- "All candidates"
STAPH_LAB <- "Staphylococcal pair"
OTHER_LAB <- "Other pairs"

COL_BOTH   <- HCOL$data      # all candidates
COL_STAPH  <- HCOL$accent    # the staphylococcal pair, amber
COL_OTHER  <- HCOL$foreign   # all other pairs, green
FILL_STAPH <- "#FBF3D5"      # pale amber, marks the outlier column in table B

ladder  <- read_excel(path, sheet = "filter_ladder")
pairs   <- read_excel(path, sheet = "genome_pairs")
cand    <- read_excel(path, sheet = "candidates_with_length")
blocks  <- read_excel(path, sheet = "shared_blocks_staph")
deciles <- read_excel(path, sheet = "length_deciles")

st <- cand %>% filter(is_staph_pair)
ot <- cand %>% filter(!is_staph_pair)
sp <- pairs %>% filter(is_staph_pair)

# =============================================================================
# TABLE A — filter ladder
# =============================================================================
tabA <- ladder %>%
  rowwise() %>%
  mutate(ci_low  = 1e5 * binom.test(n, TOTAL_READS)$conf.int[1],
         ci_high = 1e5 * binom.test(n, TOTAL_READS)$conf.int[2]) %>%
  ungroup() %>%
  transmute(` ` = c("All candidates",
                    "Flanked",
                    "Flanked, insert \u2265 2",
                    "Flanked, insert \u2265 2,\nexcluding staphylococcal pair"),
            Total = scales::comma(n),
            `Per 100k reads` = sprintf("%.1f", 1e5 * n / TOTAL_READS),
            `95% CI` = sprintf("%.1f\u2013%.1f", ci_low, ci_high))

tabA

fillA <- matrix("white", nrow = nrow(tabA), ncol = ncol(tabA))
fillA[, 1] <- HCOL$fill_head                  # row-label column
fillA[nrow(tabA), 2] <- HCOL$fill_hilite      # the final count

gA <- tableGrob(tabA, rows = NULL,
                theme = herald_table_theme(W, fillA, mult = 0.72))

# =============================================================================
# TABLE B — the two classes
# =============================================================================
tabB <- tibble(
  ` ` = c("Candidate reads (n)",
          "Percent of all candidates",
          "Enrichment vs abundance",
          "Identical sequence shared (bp)",
          "Median alignment score",
          "Median read length, relative\nto all reads",
          "Flanked architecture (n)",
          "Terminal architecture (n)"),
  `S. aureus /\nS. epidermidis` = c(
    scales::comma(nrow(st)),
    sprintf("%.1f", 100 * nrow(st) / nrow(cand)),
    sprintf("%.1fx", sp$obs_over_exp),
    scales::comma(sp$shared_bp_est),
    scales::comma(median(st$mean_AS)),
    sprintf("%.2fx", median(st$read_len) / MEDIAN_READ),
    as.character(sum(st$is_flanked)),
    as.character(sum(!st$is_flanked))),
  `All other pairs` = c(
    scales::comma(nrow(ot)),
    sprintf("%.1f", 100 * nrow(ot) / nrow(cand)),
    "~1x",
    "< 5,000",
    scales::comma(median(ot$mean_AS)),
    sprintf("%.2fx", median(ot$read_len) / MEDIAN_READ),
    as.character(sum(ot$is_flanked)),
    as.character(sum(!ot$is_flanked))))

tabB

fillB <- matrix("white", nrow = nrow(tabB), ncol = ncol(tabB))
fillB[, 1] <- HCOL$fill_head
fillB[, 2] <- FILL_STAPH

gB <- tableGrob(tabB, rows = NULL,
                theme = herald_table_theme(W, fillB, mult = 0.72))

# =============================================================================
# PANEL C — the U curve
# =============================================================================
dec_long <- deciles %>%
  select(decile, min_bp, reads, cand_staph, cand_other, cand_total) %>%
  pivot_longer(c(cand_staph, cand_other, cand_total),
               names_to = "class", values_to = "count") %>%
  mutate(class = case_when(class == "cand_staph" ~ STAPH_LAB,
                           class == "cand_other" ~ OTHER_LAB,
                           TRUE                  ~ BOTH_LAB),
         class = factor(class, levels = c(BOTH_LAB, STAPH_LAB, OTHER_LAB))) %>%
  rowwise() %>%
  mutate(rate    = 1e5 * count / reads,
         ci_low  = 1e5 * binom.test(count, reads)$conf.int[1],
         ci_high = 1e5 * binom.test(count, reads)$conf.int[2]) %>%
  ungroup()

dec_long

dec_labels <- sprintf("%.1f", deciles$min_bp / 1000)

p_u <- ggplot(dec_long, aes(x = decile, y = rate,
                            colour = class, shape = class, linetype = class)) +
  geom_errorbar(aes(ymin = ci_low, ymax = ci_high),
                width = 0.15, linewidth = 0.45 * SC, alpha = 0.75) +
  geom_line(linewidth = 1.0 * SC) +
  geom_point(size = 2.4 * SC) +
  scale_colour_manual(values = setNames(c(COL_BOTH, COL_STAPH, COL_OTHER),
                                        c(BOTH_LAB, STAPH_LAB, OTHER_LAB)),
                      name = NULL) +
  scale_shape_manual(values = setNames(c(15, 16, 17),
                                       c(BOTH_LAB, STAPH_LAB, OTHER_LAB)),
                     name = NULL) +
  scale_linetype_manual(values = setNames(c("solid", "dotted", "dashed"),
                                          c(BOTH_LAB, STAPH_LAB, OTHER_LAB)),
                        name = NULL) +
  scale_x_continuous(breaks = 1:10, labels = dec_labels,
                     sec.axis = sec_axis(~ ., breaks = 1:10, labels = 1:10,
                                         name = LAB$decile)) +
  scale_y_log10(breaks = c(0.01, 0.1, 1, 10, 100),
                labels = c("0.01", "0.1", "1", "10", "100")) +
  labs(x = LAB$decile_lo, y = LAB$cand_rate) +
  guides(colour   = guide_legend(nrow = 2, byrow = TRUE),
         shape    = guide_legend(nrow = 2, byrow = TRUE),
         linetype = guide_legend(nrow = 2, byrow = TRUE)) +
  herald_theme(W) +
  theme(legend.position = "top", legend.direction = "horizontal")

print(p_u)

# =============================================================================
# ASSEMBLE — tables stacked on the left, plot on the right
# =============================================================================
left_col <- wrap_elements(full = gA) / wrap_elements(full = gB) +
  plot_layout(heights = c(1, 1.6))

fig <- left_col | wrap_elements(full = ggplotGrob(p_u))

fig <- fig +
  plot_layout(widths = c(1.4, 1)) +
  plot_annotation(tag_levels = "A") &
  herald_tag(W)

print(fig)

herald_save(fig, "herald_pacbio_figure", W, H, outdir)

write_csv(tabA, file.path(outdir, "herald_table1A.csv"))
write_csv(tabB, file.path(outdir, "herald_table1B.csv"))

# =============================================================================
# STATISTICS
# =============================================================================
pairs %>% filter(candidates > 0) %>%
  summarise(pairs_with_candidates = n(), top_pair = max(candidates),
            top_share = max(candidates) / sum(candidates),
            second = sort(candidates, decreasing = TRUE)[2])

chisq.test(round(matrix(c(sp$candidates, sum(pairs$candidates) - sp$candidates,
                          sp$expected_candidates,
                          sum(pairs$expected_candidates) - sp$expected_candidates),
                        nrow = 2)))

blocks %>% summarise(n_blocks = n(), n_ge_500 = sum(block_length_bp >= 500),
                     n_over_fragment = sum(exceeds_fragment_length_1050bp),
                     longest = max(block_length_bp))

wilcox.test(mean_AS  ~ is_staph_pair, data = cand, conf.int = TRUE)
wilcox.test(read_len ~ is_staph_pair, data = cand, conf.int = TRUE)

deciles %>% select(decile, min_bp, max_bp, reads, cand_staph, cand_other,
                   cand_total, rate_staph_per_100k, rate_other_per_100k,
                   rate_total_per_100k, flanked, terminal) %>% as.data.frame()

mult <- read_excel(path, sheet = "length_multiples")
long_tail <- mult %>% filter(bin_lower_x_median >= 1.0)
prop.trend.test(long_tail$candidates, long_tail$reads)
prop.trend.test(long_tail$candidates[long_tail$reads >= 10000],
                long_tail$reads[long_tail$reads >= 10000])
long_tail %>% filter(bin_lower_x_median >= 1.5) %>%
  summarise(terminal = sum(terminal), flanked = sum(flanked))

cand %>% filter(is_flanked) %>%
  summarise(median_ratio = median(donor_host_ratio),
            n_ge_0.8 = sum(donor_host_ratio >= 0.8), n = n())

cand %>% filter(is_flanked, insert_frags >= 2, !is_staph_pair) %>%
  select(read, pattern, insert_frags, host, donor, read_len,
         host_AS_mean, donor_AS_mean) %>% as.data.frame()

