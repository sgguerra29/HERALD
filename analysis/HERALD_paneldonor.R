#!/usr/bin/env Rscript
# =============================================================================
# HERALD — Panel D replacement: read-level recall by donor content
#
# INPUT
#   all_strata.tsv from donor_content.py, pooled across the 20 samples:
#     run, sample, donor_bp_bin, reads, detected, recall, ci_low, ci_high
#   Counts are pooled before computing recall. CIs assume pooling
#
# CIs are Clopper-Pearson via binom.test
# =============================================================================

library(tidyverse)
library(patchwork)

FRAG_LEN <- READ_LEN / N_FRAG          

DONOR_LAB <- if (!is.null(LAB$donor_bp)) LAB$donor_bp else
  "Donor sequence contained in read (kb)"

.needed <- c("rd", "kb_lab", "herald_theme", "herald_save", "herald_tag",
             "W", "SC", "HCOL", "LAB", "READ_LEN", "N_FRAG")
.missing <- .needed[!vapply(.needed, exists, logical(1))]
if (length(.missing))
  stop("run the genome-scale script first; missing: ",
       paste(.missing, collapse = ", "), call. = FALSE)


BIN_LEVELS <- c("1-199", "200-499", "500-999", "1000-1999", ">=2000")

BIN_LABELS <- c("< 0.2", "0.2\u20130.5", "0.5\u20131", "1\u20132", "\u2265 2")
PARTIAL_BIN <- BIN_LABELS[4]

strata_raw <- rd("all_strata.tsv")

by_donor <- strata_raw %>%
  mutate(reads = as.integer(reads), detected = as.integer(detected)) %>%
  group_by(donor_bp_bin) %>%
  summarise(reads = sum(reads), detected = sum(detected), .groups = "drop") %>%
  mutate(bin = factor(donor_bp_bin, levels = BIN_LEVELS, labels = BIN_LABELS)) %>%
  arrange(bin) %>%
  rowwise() %>%
  mutate(rate    = detected / reads,
         ci_low  = binom.test(detected, reads)$conf.int[1],
         ci_high = binom.test(detected, reads)$conf.int[2]) %>%
  ungroup()

by_donor

stopifnot(nrow(by_donor) == length(BIN_LEVELS), !any(is.na(by_donor$bin)))


by_donor <- by_donor %>%
  mutate(above_floor = factor(ifelse(bin == tail(BIN_LABELS, 1),
                                     "at_or_above", "below"),
                              levels = c("below", "at_or_above")))


FILL_BELOW <- HCOL$cat2[[2]]
FILL_ABOVE <- HCOL$cat2[[1]]

p_donor <- ggplot(by_donor, aes(bin, rate, fill = above_floor)) +
  geom_col(width = 0.66, colour = "grey35", linewidth = 0.25 * SC) +
  geom_errorbar(aes(ymin = ci_low, ymax = ci_high), width = 0.14,
                linewidth = 0.45 * SC, colour = "grey25") +
  geom_text(aes(y = ci_high + 0.035,
                label = sprintf("%.1f%%\n%s/%s", 100 * rate,
                                scales::comma(detected),
                                scales::comma(reads))),
            vjust = 0, size = 2.9 * SC, family = "lato", lineheight = 0.95,
            colour = "grey15") +
  annotate("text", x = 4.42, y = 1.1,
           label = sprintf("2kb fragment length", kb_lab(FRAG_LEN)),
           hjust = 1, vjust = 0, size = 2.9 * SC, family = "lato",
           fontface = 3, colour = "grey35") +
  scale_fill_manual(values = c(below = FILL_BELOW, at_or_above = FILL_ABOVE),
                    guide = "none") +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1),
                     limits = c(0, 1.22), breaks = seq(0, 1, 0.25),
                     expand = c(0, 0)) +
  labs(x = DONOR_LAB, y = LAB$recall) +
  herald_theme(W) +
  theme(plot.margin = margin(6, 14, 6, 6))

print(p_donor)

if (all(vapply(c("p_size", "p_mmf", "p_pair"), exists, logical(1)))) {
  
  right_col <- p_mmf / p_pair + plot_layout(heights = c(1, 1))
  
  fig <- (p_size / p_donor + plot_layout(heights = c(1, 1))) | right_col
  
  fig <- fig +
    plot_layout(widths = c(1.2, 1)) +
    plot_annotation(tag_levels = "A") &
    herald_tag(W)
  
  print(fig)
  
  herald_save(fig, "herald_genomescale_figure", W, H * 1.15, outdir)
  
} else {
  message("p_size / p_mmf / p_pair not in scope: panel built, figure not ",
          "reassembled. Source the genome-scale script first to rebuild it.")
}
write_csv(by_donor %>%
            transmute(`Donor bp in read` = as.character(bin),
                      `Reads` = reads,
                      `Reads detected` = detected,
                      `Recall (%)` = sprintf("%.1f", 100 * rate),
                      `95% CI` = sprintf("%.1f\u2013%.1f",
                                         100 * ci_low, 100 * ci_high)),
          file.path(outdir, "herald_suppS6_recall_by_donor_content.csv"))

# =============================================================================
# STATISTICS
# =============================================================================
cat("\n--- recall by donor content, pooled across 20 samples ---\n")
by_donor %>%
  transmute(bin, reads, detected, recall = sprintf("%.1f%%", 100 * rate),
            ci = sprintf("%.3f\u2013%.3f", ci_low, ci_high)) %>%
  as.data.frame()

cat(sprintf("\nTOTAL  %s / %s = %.1f%%\n",
            scales::comma(sum(by_donor$detected)),
            scales::comma(sum(by_donor$reads)),
            100 * sum(by_donor$detected) / sum(by_donor$reads)))

fl <- by_donor %>% filter(above_floor == "at_or_above")
cat(sprintf("At or above one fragment length: %s / %s = %.1f%%\n",
            scales::comma(fl$detected), scales::comma(fl$reads),
            100 * fl$rate))

cat("\n--- below vs at/above one fragment length ---\n")
regime <- by_donor %>%
  mutate(regime = ifelse(above_floor == "at_or_above",
                         ">=1 frag len", "<1 frag len")) %>%
  group_by(regime) %>%
  summarise(detected = sum(detected),
            missed = sum(reads) - sum(detected), .groups = "drop")
regime_m <- as.matrix(regime[, -1]); rownames(regime_m) <- regime$regime
print(regime_m)
print(fisher.test(regime_m))

cat("\n--- the partial bin ---\n")
by_donor %>% filter(bin == PARTIAL_BIN) %>%
  transmute(reads, detected, recall = sprintf("%.1f%%", 100 * rate),
            ci = sprintf("%.1f\u2013%.1f%%", 100 * ci_low, 100 * ci_high)) %>%
  as.data.frame()