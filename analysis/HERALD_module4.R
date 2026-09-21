#!/usr/bin/env Rscript
# =============================================================================
# HERALD — genome-scale synthetic module
# =============================================================================
library(tidyverse)
library(patchwork)
library(gridExtra)
library(grid)
library(showtext)
# Directory holding the input data files and receiving the output figures.
# Defaults to the current working directory; override without editing this file:
#     HERALD_DIR=/path/to/data Rscript HERALD_module4.R
outdir <- Sys.getenv("HERALD_DIR", unset = getwd())
source(file.path(outdir, "HERALD_theme.R"))

W <- 12.5; H <- 8.6
SC <- herald_sc(W)

PANEL_LAB <- "100-genome reference"
LIMIT_LAB <- "2-genome reference"

MMF_OP   <- 0.95
N_FRAG   <- 10
READ_LEN <- 20000

PAIR_KEY <- tribble(
  ~pair,                ~acceptor,      ~donor,
  "pair1_phylum",       "AE017282.2",   "CP000976.1",
  "pair2_family",       "CP047394.1",   "CP033049.1",
  "pair3_genus_small",  "CP000891.1",   "CP000302.1",
  "pair4_genus_large",  "CP024307.1",   "CP013107.1")

PAIR_LAB <- c(pair1_phylum      = "Different phyla",
              pair2_family      = "Same family",
              pair3_genus_small = "Same genus (i)",
              pair4_genus_large = "Same genus (ii)")

PAIR_ORDER <- c("pair1_phylum", "pair2_family",
                "pair3_genus_small", "pair4_genus_large")

rd <- function(f, required = TRUE) {
  p <- file.path(outdir, f)
  if (!file.exists(p)) {
    if (required) stop("missing input: ", p, call. = FALSE)
    warning("optional input not found: ", p, call. = FALSE)
    return(NULL)
  }
  read.delim(p, stringsAsFactors = FALSE)
}

tidy_ev <- function(d, lab) {
  d %>% mutate(arm = lab,
               insert_len      = as.integer(insert_len),
               reads_available = as.integer(reads_available),
               reads_detected  = as.integer(reads_detected),
               found = detected == "yes",
               pair  = sub("_rep[0-9]+$", "", sample))
}

pev <- tidy_ev(rd("panel_events.tsv"),      PANEL_LAB)
lev <- tidy_ev(rd("limitedref_events.tsv"), LIMIT_LAB)
sev <- tidy_ev(rd("sweep_events.tsv"),      "sweep") %>%
  mutate(mmf = as.integer(sub(".*mmf\\.?([0-9]+).*", "\\1", run)) / 100)
psa <- rd("panel_samples.tsv")
ovl <- rd("overlap/pairs.tsv", required = FALSE)

stopifnot(!any(is.na(sev$mmf)))

size_levels <- kb_lab(sort(unique(pev$insert_len)))
as_size_lab <- function(x) factor(kb_lab(x), levels = size_levels)

# =============================================================================
# PANEL A — recall by insert length, both reference sets, shared samples only
# =============================================================================
shared_samples <- intersect(unique(pev$sample), unique(lev$sample))
cat("samples scored under both reference sets:", length(shared_samples), "\n")
stopifnot(length(shared_samples) == 8)

by_arm <- bind_rows(pev, lev) %>%
  filter(sample %in% shared_samples) %>%
  group_by(arm, insert_len) %>%
  summarise(events = n(), avail = sum(reads_available),
            found_reads = sum(reads_detected),
            events_found = sum(found), .groups = "drop") %>%
  rowwise() %>%
  mutate(rate    = found_reads / avail,
         ci_low  = binom.test(found_reads, avail)$conf.int[1],
         ci_high = binom.test(found_reads, avail)$conf.int[2]) %>%
  ungroup() %>%
  mutate(arm = factor(arm, levels = c(PANEL_LAB, LIMIT_LAB)),
         size_lab = as_size_lab(insert_len))

by_arm    # the two arms should agree exactly

lab_dat <- by_arm %>%
  group_by(size_lab) %>%
  summarise(rate = max(rate), top = max(ci_high), .groups = "drop")

p_size <- ggplot(by_arm, aes(size_lab, rate, fill = arm)) +
  geom_col(position = position_dodge(0.72), width = 0.66,
           colour = "grey35", linewidth = 0.25 * SC) +
  geom_errorbar(aes(ymin = ci_low, ymax = ci_high),
                position = position_dodge(0.72), width = 0.14,
                linewidth = 0.45 * SC, colour = "grey25") +
  geom_text(data = lab_dat, inherit.aes = FALSE,
            aes(size_lab, top + 0.05, label = sprintf("%.1f%%", 100 * rate)),
            size = 3.2 * SC, family = "lato", fontface = 2, colour = "grey15") +
  scale_fill_manual(values = setNames(HCOL$cat2,
                                      c(PANEL_LAB, LIMIT_LAB)), name = NULL) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1),
                     limits = c(0, 1.12), breaks = seq(0, 1, 0.25),
                     expand = c(0, 0)) +
  labs(x = LAB$insert, y = LAB$recall) +
  herald_theme(W) +
  theme(legend.position = "top", legend.direction = "horizontal",
        legend.justification = "left",
        legend.margin = margin(l = -10, r = 0),
        legend.text = element_text(size = 0.70 * herald_base(W)),
        plot.margin = margin(6, 14, 6, 6))

print(p_size)

# =============================================================================
# Shared tile aesthetics — insert length on x in BOTH panels
# =============================================================================
tile_base <- list(
  geom_tile(colour = "white", linewidth = 1.1 * SC),
  geom_text(aes(label = lab, colour = ifelse(rate > 0.5, "hi", "lo")),
            size = 2.9 * SC, family = "lato", lineheight = 0.95),
  scale_colour_manual(values = c(hi = "white", lo = "grey25"),
                      guide = "none"),
  scale_x_discrete(expand = c(0, 0)),
  scale_y_discrete(expand = c(0, 0)),
  herald_fill_scale()
)

tile_theme <- herald_theme(W) +
  theme(axis.line = element_blank(), axis.ticks = element_blank(),
        plot.margin = margin(6, 6, 4, 6))

# =============================================================================
# PANEL B — mmf on y, insert length on x
# =============================================================================
sw <- sev %>%
  group_by(mmf, insert_len) %>%
  summarise(avail = sum(reads_available), found_reads = sum(reads_detected),
            events = n(), events_found = sum(found), .groups = "drop") %>%
  mutate(rate = ifelse(avail > 0, found_reads / avail, 0),
         mmf_lab  = factor(sprintf("%.2f", mmf)),
         size_lab = as_size_lab(insert_len),
         lab = sprintf("%.0f%%\n%d/%d", 100 * rate, events_found, events))

sw

p_mmf <- ggplot(sw, aes(size_lab, mmf_lab, fill = rate)) +
  tile_base +
  labs(x = NULL, y = LAB$mmf) +
  tile_theme +
  theme(
    axis.text.x  = element_blank(),
    axis.title.y = element_text(margin = margin(r = -25 * SC)),
    legend.position = "right")

print(p_mmf)

# =============================================================================
# PANEL C — acceptor/donor pair on y, insert length on x
# =============================================================================
by_pair <- pev %>%
  group_by(pair, insert_len) %>%
  summarise(avail = sum(reads_available), found_reads = sum(reads_detected),
            events = n(), events_found = sum(found), .groups = "drop") %>%
  mutate(rate = ifelse(avail > 0, found_reads / avail, 0),
         size_lab = as_size_lab(insert_len),
         # reversed so the least-similar pair sits at the top
         pair_lab = factor(PAIR_LAB[pair], levels = rev(PAIR_LAB[PAIR_ORDER])),
         lab = sprintf("%.0f%%\n%d/%d", 100 * rate, events_found, events))

by_pair %>% select(pair, insert_len, avail, found_reads, rate, events_found)

p_pair <- ggplot(by_pair, aes(size_lab, pair_lab, fill = rate)) +
  tile_base +
  labs(x = LAB$insert, y = NULL) +
  tile_theme +
  theme(legend.position = "none")       # legend lives on panel B

print(p_pair)

# =============================================================================
# PANEL D — primary results, four columns plus a Total row
# =============================================================================
by_size <- pev %>%
  group_by(insert_len) %>%
  summarise(events = n(), avail = sum(reads_available),
            found_reads = sum(reads_detected),
            events_found = sum(found), .groups = "drop") %>%
  rowwise() %>%
  mutate(rate    = found_reads / avail,
         ci_low  = binom.test(found_reads, avail)$conf.int[1],
         ci_high = binom.test(found_reads, avail)$conf.int[2]) %>%
  ungroup()

by_size

tot_av <- sum(by_size$avail)
tot_fd <- sum(by_size$found_reads)
tot_ci <- binom.test(tot_fd, tot_av)$conf.int

tabD <- by_size %>%
  transmute(`Insert (kb)`    = kb_lab(insert_len),
            `Reads detected` = paste0(scales::comma(found_reads), " / ",
                                      scales::comma(avail)),
            `Recall (%)`     = sprintf("%.1f", 100 * rate),
            `Events`         = paste0(events_found, "/", events)) %>%
  bind_rows(tibble(
    `Insert (kb)`    = "Total",
    `Reads detected` = paste0(scales::comma(tot_fd), " / ",
                              scales::comma(tot_av)),
    `Recall (%)`     = sprintf("%.1f", 100 * tot_fd / tot_av),
    `Events`         = paste0(sum(by_size$events_found), "/",
                              sum(by_size$events))))

tabD

n_row <- nrow(tabD)
fillD <- matrix("white", nrow = n_row, ncol = ncol(tabD))
fillD[, 1] <- HCOL$fill_head
fillD[, 3] <- HCOL$fill_hilite    

faceD <- matrix(1, nrow = n_row, ncol = ncol(tabD))
faceD[n_row, ] <- 2

gD <- tableGrob(tabD, rows = NULL,
                theme = herald_table_theme(W, fillD, faceD, mult = 0.85))

# Full version, with the confidence intervals, for the supplement
tabD_full <- by_size %>%
  transmute(`Insert (kb)`      = kb_lab(insert_len),
            `Reads with donor` = scales::comma(avail),
            `Reads detected`   = scales::comma(found_reads),
            `Recall (%)`       = sprintf("%.1f", 100 * rate),
            `95% CI`           = sprintf("%.1f\u2013%.1f",
                                         100 * ci_low, 100 * ci_high),
            `Events detected`  = paste0(events_found, "/", events)) %>%
  bind_rows(tibble(
    `Insert (kb)`      = "Total",
    `Reads with donor` = scales::comma(tot_av),
    `Reads detected`   = scales::comma(tot_fd),
    `Recall (%)`       = sprintf("%.1f", 100 * tot_fd / tot_av),
    `95% CI`           = sprintf("%.1f\u2013%.1f",
                                 100 * tot_ci[1], 100 * tot_ci[2]),
    `Events detected`  = paste0(sum(by_size$events_found), "/",
                                sum(by_size$events))))

tabD_full

# =============================================================================
# SUPPLEMENTARY — the four acceptor/donor pairs, full detail
# =============================================================================
shared_for <- function(acc, don) {
  if (is.null(ovl)) return(c(NA_real_, NA_real_))
  hit <- ovl %>% filter((genome_a == acc & genome_b == don) |
                          (genome_a == don & genome_b == acc))
  if (nrow(hit) == 0) return(c(0, 0))
  h <- hit[1, ]
  # the ACCEPTOR's shared bp; the two totals differ because merged intervals
  # are not symmetric when a region is repeated in one of the genomes
  bp <- if (h$genome_a == acc) h$shared_bp_a else h$shared_bp_b
  c(as.numeric(bp), as.numeric(h$max_block))
}

pair_overlap <- PAIR_KEY %>%
  rowwise() %>%
  mutate(shared_bp = shared_for(acceptor, donor)[1],
         max_block = shared_for(acceptor, donor)[2]) %>%
  ungroup()

pair_overlap  

pair_stats <- pev %>%
  group_by(pair) %>%
  summarise(events = n(), events_found = sum(found),
            avail = sum(reads_available), found_reads = sum(reads_detected),
            .groups = "drop") %>%
  left_join(
    psa %>% mutate(pair = sub("_rep[0-9]+$", "", sample)) %>%
      group_by(pair) %>%
      summarise(samples = n(), reads = sum(total_reads),
                flagged = sum(reads_flagged), tp = sum(true_positive),
                mis = sum(misattributed), fp = sum(false_positive),
                .groups = "drop"),
    by = "pair")

tabPairs <- pair_overlap %>%
  left_join(pair_stats, by = "pair") %>%
  transmute(Acceptor = acceptor, Donor = donor,
            Relation = PAIR_LAB[pair],
            `Shared (bp)` = ifelse(is.na(shared_bp), "n/a",
                                   scales::comma(shared_bp)),
            `Max block (bp)` = ifelse(is.na(max_block) | max_block == 0,
                                      "\u2013", scales::comma(max_block)),
            `Reads with donor` = scales::comma(avail),
            `Reads detected` = scales::comma(found_reads),
            `Recall (%)` = sprintf("%.1f", 100 * found_reads / avail),
            `Reads flagged` = scales::comma(flagged),
            `False positives` = as.character(fp),
            `Misattributed` = as.character(mis),
            `Events detected` = paste0(events_found, "/", events))

tabPairs

# =============================================================================
# ASSEMBLE
#   left  column: A over D
#   right column: B over C, axes aligned via native patchwork "/"
# =============================================================================
right_col <- p_mmf / p_pair + plot_layout(heights = c(1, 1))

fig <- (p_size / wrap_elements(full = gD) + plot_layout(heights = c(1.35, 1))) |
  right_col

fig <- fig +
  plot_layout(widths = c(1.2, 1)) +
  plot_annotation(tag_levels = "A") &
  herald_tag(W)

print(fig)

herald_save(fig, "herald_genomescale_figure", W, H, outdir)

write_csv(tabD_full, file.path(outdir, "herald_suppS2_recall_by_size.csv"))
write_csv(tabPairs,  file.path(outdir, "herald_suppS3_pairs.csv"))
by_arm %>%
  transmute(arm, `Insert (kb)` = kb_lab(insert_len),
            `Reads with donor` = avail, `Reads detected` = found_reads,
            `Recall (%)` = sprintf("%.1f", 100 * rate),
            `Events detected` = paste0(events_found, "/", events)) %>%
  write_csv(file.path(outdir, "herald_suppS4_reference_arms.csv"))

# =============================================================================
# STATISTICS
# =============================================================================
cat("\n--- primary arm, overall (all 20 samples) ---\n")
pev %>% summarise(samples = n_distinct(sample), events = n(),
                  reads_with_donor = sum(reads_available),
                  reads_detected = sum(reads_detected),
                  recall = sum(reads_detected) / sum(reads_available),
                  events_detected = sum(found)) %>% as.data.frame()

cat("\n--- reads detected, >=2 kb vs <=1 kb inserts ---\n")
# Structural zeros at 0.5 kb make a five-level test degenerate.
ct <- pev %>%
  mutate(regime = ifelse(insert_len >= 2000, ">=2kb", "<=1kb")) %>%
  group_by(regime) %>%
  summarise(detected = sum(reads_detected),
            missed = sum(reads_available) - sum(reads_detected),
            .groups = "drop")
ct_m <- as.matrix(ct[, -1]); rownames(ct_m) <- ct$regime
print(ct_m); print(fisher.test(ct_m))

cat("\n--- reference-set size: paired on matched events ---\n")
matched <- inner_join(
  pev %>% select(sample, event_id, p_reads = reads_detected, p_found = found),
  lev %>% select(sample, event_id, l_reads = reads_detected, l_found = found),
  by = c("sample", "event_id"))
matched %>% summarise(matched_events = n(), samples = n_distinct(sample),
                      read_counts_differ = sum(p_reads != l_reads),
                      detection_discordant = sum(p_found != l_found)) %>%
  as.data.frame()
if (all(matched$p_reads == matched$l_reads)) {
  cat("Detection IDENTICAL read-for-read across reference sets.\n",
      "Report as an identity in the caption; no test is applicable.\n", sep = "")
} else {
  print(mcnemar.test(table(matched$p_found, matched$l_found)))
  print(wilcox.test(matched$p_reads, matched$l_reads, paired = TRUE))
}

cat("\n--- pair x insert size, recall (%) ---\n")
by_pair %>% mutate(v = sprintf("%.1f", 100 * rate)) %>%
  select(pair, insert_len, v) %>%
  pivot_wider(names_from = insert_len, values_from = v) %>% as.data.frame()

cat("\n--- pair x insert size, events detected ---\n")
by_pair %>% mutate(v = paste0(events_found, "/", events)) %>%
  select(pair, insert_len, v) %>%
  pivot_wider(names_from = insert_len, values_from = v) %>% as.data.frame()

cat("\n--- spread within each insert size, across pairs ---\n")
by_pair %>% group_by(insert_len) %>%
  summarise(min = 100 * min(rate), max = 100 * max(rate),
            spread_pts = 100 * (max(rate) - min(rate)), .groups = "drop") %>%
  as.data.frame()

cat("\n--- pair effect at 2 kb and above ---\n")
bp2 <- pev %>% filter(insert_len >= 2000) %>%
  group_by(pair) %>%
  summarise(detected = sum(reads_detected),
            missed = sum(reads_available) - sum(reads_detected),
            .groups = "drop")
bp2_m <- as.matrix(bp2[, -1]); rownames(bp2_m) <- bp2$pair
print(bp2_m); print(chisq.test(bp2_m))
cat("A non-significant result supports acceptor/donor similarity not affecting\n",
    "detection across the shared-sequence gradient.\n", sep = "")

cat("\n--- mmf sweep: recall (%) ---\n")
sw %>% mutate(v = sprintf("%.0f", 100 * rate)) %>%
  select(mmf, insert_len, v) %>%
  pivot_wider(names_from = insert_len, values_from = v) %>% as.data.frame()

cat("\n--- mmf sweep: events detected ---\n")
sw %>% mutate(v = paste0(events_found, "/", events)) %>%
  select(mmf, insert_len, v) %>%
  pivot_wider(names_from = insert_len, values_from = v) %>% as.data.frame()

sev %>% group_by(mmf) %>%
  summarise(pooled_recall = sum(reads_detected) / sum(reads_available),
            events_detected = sum(found), events = n(), .groups = "drop") %>%
  as.data.frame()

cat("\n--- trend across mmf (event-level logistic) ---\n")
# Event-level, not read-level: reads within an event overlap the same insert and are not independent.
m_mmf <- glm(found ~ mmf, data = sev, family = binomial)
print(summary(m_mmf)$coefficients)
cat(sprintf("odds ratio per 0.10 increase in mmf: %.2f\n",
            exp(coef(m_mmf)[["mmf"]] * 0.10)))

cat("\n--- false positives and misattribution, all samples ---\n")
psa %>% summarise(samples = n(), total_reads = sum(total_reads),
                  flagged = sum(reads_flagged), tp = sum(true_positive),
                  misattributed = sum(misattributed),
                  false_positive = sum(false_positive),
                  precision = 100 * sum(true_positive) /
                    (sum(true_positive) + sum(false_positive))) %>%
  as.data.frame()

if (!any(psa$is_negative == "yes")) {
  cat("\nNOTE: negative controls are absent from panel_samples.tsv because\n",
      "samples yielding zero flagged reads produce no output rows. Zero output\n",
      "IS zero false positives, but state the read counts for all eight\n",
      "negative controls explicitly in the manuscript text.\n", sep = "")
}