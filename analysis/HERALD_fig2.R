#!/usr/bin/env Rscript
# =============================================================================
# HERALD — parameter module (minimum matched fraction and fragment count)
# =============================================================================

library(tidyverse)
library(readxl)
library(patchwork)

# Directory holding the input data files and receiving the output figures.
# Defaults to the current working directory; override without editing this file:
#     HERALD_DIR=/path/to/data Rscript HERALD_fig2.R
outdir <- Sys.getenv("HERALD_DIR", unset = getwd())
source(file.path(outdir, "HERALD_theme.R"))

W <- 12.5; H <- 8.0
SC <- herald_sc(W)

path <- file.path(outdir, "HERALD_plottingdata.xlsx")
FRAG_COLS <- setNames(ramp_blue(3), c("5", "10", "15"))

# =============================================================================
# MINIMUM MATCHED FRACTION
# =============================================================================
mmf_data <- read_excel(path, sheet = "mmf", skip = 1)
glimpse(mmf_data)

# ── Panel A: false positives ────────────────────────────────────────────────
fp_mmf_clean <- mmf_data %>%
  filter(insert == "0", `read length` != "total") %>%
  mutate(read_length     = as.numeric(str_remove(`read length`, "kb")),
         mmf             = `min-fraction-aligned`,
         normal_reads    = `normal reads`,
         false_positives = `# false positives`,
         corrected_fp    = `#false corrected`) %>%
  select(read_length, mmf, normal_reads, false_positives, corrected_fp) %>%
  rowwise() %>%
  mutate(fp_rate = false_positives / normal_reads,
         ci_low  = binom.test(false_positives, normal_reads)$conf.int[1],
         ci_high = binom.test(false_positives, normal_reads)$conf.int[2]) %>%
  ungroup()

fp_mmf_clean

stopifnot(!any(is.na(fp_mmf_clean$corrected_fp)))

model_fp_mmf <- glm(cbind(false_positives, normal_reads - false_positives) ~
                      mmf + read_length,
                    family = binomial, data = fp_mmf_clean)
summary(model_fp_mmf)

fp_mmf_long <- fp_mmf_clean %>%
  select(read_length, mmf, false_positives, corrected_fp) %>%
  pivot_longer(c(false_positives, corrected_fp),
               names_to = "type", values_to = "count") %>%
  mutate(type = factor(type, levels = c("false_positives", "corrected_fp")))

p_fp_mmf <- ggplot(fp_mmf_long, aes(x = mmf, y = count,
                                    colour = factor(read_length),
                                    linetype = type)) +
  geom_line(linewidth = 0.9 * SC) +
  geom_point(size = 2.2 * SC) +
  scale_colour_manual(values = HCOL$readlen, name = LAB$readlen) +
  scale_linetype_manual(
    breaks = c("false_positives", "corrected_fp"),
    values = c("22", "solid"),
    labels = c("Raw", "Corrected"),
    name = NULL) +
  scale_x_continuous(breaks = c(0.15, 0.35, 0.55, 0.75)) +
  labs(x = LAB$mmf, y = LAB$fp_n) +
  guides(linetype = guide_legend(
    override.aes = list(linewidth = 0.5, colour = "grey20"))) +
  herald_theme(W) +
  theme(legend.position = "right",
        legend.key.width = unit(1.1, "cm"))

print(p_fp_mmf)

# ── Panel B: recall ─────────────────────────────────────────────────────────
detect_mmf_clean <- mmf_data %>%
  filter(insert != "0", `read length` != "total") %>%
  mutate(insert_length  = as.numeric(str_remove(insert, "kb")),
         mmf            = `min-fraction-aligned`,
         hgt_reads      = `hgt reads`,
         true_positives = `# true positives`) %>%
  select(insert_length, mmf, hgt_reads, true_positives) %>%
  rowwise() %>%
  mutate(recall  = true_positives / hgt_reads,
         ci_low  = binom.test(true_positives, hgt_reads)$conf.int[1],
         ci_high = binom.test(true_positives, hgt_reads)$conf.int[2]) %>%
  ungroup()

detect_mmf_clean

# 1.6 kb sits at 0% for every mmf level and causes quasi-separation, so the
# models are fitted without it and it is excluded from panel B. The floor is
# reported in the module-1 figure and in the text instead.
detect_mmf_no_floor <- detect_mmf_clean %>% filter(insert_length != 1.6)

model_full <- glm(cbind(true_positives, hgt_reads - true_positives) ~
                    factor(insert_length) * mmf,
                  family = quasibinomial, data = detect_mmf_no_floor)
model_reduced <- glm(cbind(true_positives, hgt_reads - true_positives) ~ mmf,
                     family = quasibinomial, data = detect_mmf_no_floor)

anova(model_reduced, model_full, test = "F")  # insert length can be dropped
summary(model_reduced)                        # the model we report

insert_lens <- sort(unique(detect_mmf_no_floor$insert_length))
INSERT_COLS <- setNames(ramp_n(length(insert_lens)),
                        as.character(insert_lens))

p_detect_mmf <- ggplot(detect_mmf_no_floor,
                       aes(x = mmf, y = recall,
                           colour = factor(insert_length))) +
  geom_errorbar(aes(ymin = ci_low, ymax = ci_high),
                width = 0.02, linewidth = 0.5 * SC) +
  geom_line(linewidth = 0.9 * SC) +
  geom_point(size = 2.2 * SC) +
  scale_colour_manual(values = INSERT_COLS, name = LAB$insert) +
  scale_y_continuous(labels = scales::percent, limits = c(0.4, 1.0)) +
  scale_x_continuous(breaks = c(0.15, 0.35, 0.55, 0.75)) +
  labs(x = LAB$mmf, y = LAB$recall) +
  herald_theme(W) +
  theme(legend.position = "right")

print(p_detect_mmf)

# =============================================================================
# FRAGMENT COUNT
# =============================================================================
frag_data <- read_excel(path, sheet = "frag number", skip = 1,
                        na = c("na", "NA")) %>%
  fill(`fragments made`, .direction = "down")
glimpse(frag_data)

# ── Panel C: placement, insert fixed at 7.6 kb ──────────────────────────────
frag_placement_clean <- frag_data %>%
  filter(insert == "7.6kb", !is.na(`left host flank`)) %>%
  mutate(left_flank     = as.numeric(str_remove(`left host flank`, "kb")),
         fragments      = `fragments made`,
         hgt_reads      = `hgt reads`,
         true_positives = `# true positives`) %>%
  select(fragments, left_flank, hgt_reads, true_positives) %>%
  rowwise() %>%
  mutate(recall  = true_positives / hgt_reads,
         ci_low  = binom.test(true_positives, hgt_reads)$conf.int[1],
         ci_high = binom.test(true_positives, hgt_reads)$conf.int[2]) %>%
  ungroup()

frag_placement_clean

# binomial first to check dispersion, then quasibinomial for the reported p
model_placement_frag <- glm(cbind(true_positives, hgt_reads - true_positives) ~
                              left_flank * fragments,
                            family = binomial, data = frag_placement_clean)
summary(model_placement_frag)

model_placement_frag_q <- glm(cbind(true_positives, hgt_reads - true_positives) ~
                                left_flank * fragments,
                              family = quasibinomial,
                              data = frag_placement_clean)
summary(model_placement_frag_q)   # ns, all p > 0.4

p_placement_frag <- ggplot(frag_placement_clean,
                           aes(x = left_flank, y = recall,
                               colour = factor(fragments))) +
  geom_errorbar(aes(ymin = ci_low, ymax = ci_high),
                position = position_dodge(width = 1), width = 0.6,
                linewidth = 0.5 * SC) +
  geom_point(position = position_dodge(width = 1), size = 2.2 * SC) +
  scale_colour_manual(values = FRAG_COLS, name = LAB$frags) +
  scale_y_continuous(labels = scales::percent) +
  scale_x_continuous(breaks = frag_placement_clean$left_flank) +
  labs(x = LAB$flank, y = LAB$recall) +
  herald_theme(W) +
  theme(legend.position = "right")

print(p_placement_frag)

# ── Panel D: insert length ──────────────────────────────────────────────────
frag_insertlen_clean <- frag_data %>%
  filter(insert != "7.6kb", !is.na(insert)) %>%
  mutate(insert_length  = as.numeric(str_remove(insert, "kb")),
         fragments      = `fragments made`,
         hgt_reads      = `hgt reads`,
         true_positives = `# true positives`) %>%
  select(fragments, insert_length, hgt_reads, true_positives) %>%
  rowwise() %>%
  mutate(recall  = true_positives / hgt_reads,
         ci_low  = binom.test(true_positives, hgt_reads)$conf.int[1],
         ci_high = binom.test(true_positives, hgt_reads)$conf.int[2]) %>%
  ungroup()

frag_insertlen_clean

# Two targeted tests instead of the full factorial, which quasi-separated:
#   Q1 does the 1.6 kb floor lift between low (5/10) and high (15) fragments?
#   Q2 at 15 fragments, is 1.6 kb still distinguishable from the plateau?
floor_only <- frag_insertlen_clean %>% filter(insert_length == 1.6)
floor_only
prop.test(x = c(0 + 0, 204), n = c(250 + 250, 250))   # Q1

frag15_only <- frag_insertlen_clean %>% filter(fragments == 15)
frag15_only
table_15 <- rbind(frag15_only$true_positives,
                  frag15_only$hgt_reads - frag15_only$true_positives)
colnames(table_15) <- frag15_only$insert_length
chisq.test(table_15)                                   # Q2
pairwise.prop.test(x = frag15_only$true_positives,
                   n = frag15_only$hgt_reads,
                   p.adjust.method = "holm")

frag_insertlen_no_floor <- frag_insertlen_clean %>% filter(insert_length != 1.6)
model_insertlen_frag_q <- glm(cbind(true_positives, hgt_reads - true_positives) ~
                                insert_length * fragments,
                              family = quasibinomial,
                              data = frag_insertlen_no_floor)
summary(model_insertlen_frag_q)   # ns, all p > 0.5

# The full data are plotted, including 1.6 kb: the 0% floor at 5 and 10
# fragments and its recovery at 15 is the point of the panel.
p_insertlen_frag <- ggplot(frag_insertlen_clean,
                           aes(x = insert_length, y = recall,
                               colour = factor(fragments))) +
  geom_errorbar(aes(ymin = ci_low, ymax = ci_high),
                position = position_dodge(width = 0.5), width = 0.3,
                linewidth = 0.5 * SC) +
  geom_point(position = position_dodge(width = 0.5), size = 2.2 * SC) +
  scale_colour_manual(values = FRAG_COLS, name = LAB$frags) +
  scale_y_continuous(labels = scales::percent, limits = c(0, 1)) +
  scale_x_continuous(breaks = frag_insertlen_clean$insert_length) +
  labs(x = LAB$insert, y = LAB$recall) +
  herald_theme(W) +
  theme(legend.position = "right")

print(p_insertlen_frag)

# =============================================================================
# SUPPLEMENTARY — false positives x fragment count
# =============================================================================
fp_frag_clean <- frag_data %>%
  filter(is.na(`left host flank`), is.na(insert)) %>%
  mutate(source_block    = case_when(`hgt reads` == 2500 ~ "placement",
                                     `hgt reads` == 1750 ~ "insert_length"),
         fragments       = `fragments made`,
         normal_reads    = `normal reads`,
         false_positives = `# false positives`,
         corrected_fp    = `false positives: corrected`) %>%
  select(source_block, fragments, normal_reads, false_positives, corrected_fp)

fp_frag_clean

model_fp_frag <- glm(cbind(false_positives, normal_reads - false_positives) ~
                       fragments + source_block,
                     family = binomial, data = fp_frag_clean)
summary(model_fp_frag)   # fragments p = 0.034, source_block p = 0.005; n = 6

# At 5 fragments the correction does not applry. 
fp_frag_long <- fp_frag_clean %>%
  mutate(corrected_display = if_else(fragments == 5,
                                     false_positives, corrected_fp)) %>%
  select(source_block, fragments, false_positives, corrected_display) %>%
  pivot_longer(c(false_positives, corrected_display),
               names_to = "type", values_to = "count")

p_fp_frag <- ggplot(fp_frag_long, aes(x = fragments, y = count,
                                      colour = source_block,
                                      linetype = type)) +
  geom_line(linewidth = 0.9 * SC) +
  geom_point(size = 2.5 * SC) +
  scale_colour_manual(
    values = setNames(HCOL$cat2, c("placement", "insert_length")),
    labels = c(placement = "Placement design",
               insert_length = "Insert-length design"), name = NULL) +
  scale_linetype_manual(
    values = c(false_positives = "dashed", corrected_display = "solid"),
    labels = c(false_positives = "Raw", corrected_display = "Corrected*"),
    name = NULL) +
  annotate("text", x = 10,
           y = max(fp_frag_long$count, na.rm = TRUE) * 0.95,
           label = "Fragment count: p = 0.034 (n = 6)",
           size = 3.0 * SC, fontface = "italic", colour = "grey40",
           family = "lato") +
  scale_x_continuous(breaks = c(5, 10, 15)) +
  labs(x = LAB$frags, y = LAB$fp_n) +
  herald_theme(W) +
  theme(legend.position = "right")

print(p_fp_frag)

# =============================================================================
# ASSEMBLE
# =============================================================================
fig <- (p_fp_mmf + p_detect_mmf) / (p_placement_frag + p_insertlen_frag) +
  plot_annotation(tag_levels = "A") &
  herald_tag(W)

print(fig)

herald_save(fig, "herald_module2_figure", W, H, outdir)

# supplementary, on its own
herald_save(p_fp_frag + herald_tag(W),
            "herald_suppS5_fp_by_fragments", 6.5, 4.5, outdir)