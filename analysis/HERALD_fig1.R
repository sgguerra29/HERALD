#!/usr/bin/env Rscript
# =============================================================================
# HERALD — synthetic read-level module (module 1)
# =============================================================================

library(tidyverse)
library(readxl)
library(patchwork)

# Directory holding the input data files and receiving the output figures.
# Defaults to the current working directory; override without editing this file:
#     HERALD_DIR=/path/to/data Rscript HERALD_fig1.R
outdir <- Sys.getenv("HERALD_DIR", unset = getwd())
source(file.path(outdir, "HERALD_theme.R"))

W <- 12.5; H <- 8.0
SC <- herald_sc(W)

path <- file.path(outdir, "HERALD_plottingdata.xlsx")

# =============================================================================
# FALSE POSITIVES x ERROR RATE
# =============================================================================
fp_data <- read_excel(path, sheet = "false positive", skip = 1)
# skip = 1: the sheet has a merged header row above the real column names
glimpse(fp_data)

fp_clean <- fp_data %>%
  filter(`read length` != "total") %>%          # drop the 4 summary rows
  mutate(read_length     = as.numeric(str_remove(`read length`, "kb")),
         error_rate      = `error rate`,
         normal_reads    = `normal reads`,
         false_positives = `# false positives`,
         single_frag_fp  = `# false + w single frag`,
         corrected_fp    = `# false corrected`) %>%
  select(read_length, error_rate, normal_reads,
         false_positives, single_frag_fp, corrected_fp)

glimpse(fp_clean)

stopifnot(!any(is.na(fp_clean$corrected_fp)))

fp_long <- fp_clean %>%
  select(read_length, error_rate, false_positives, corrected_fp) %>%
  pivot_longer(c(false_positives, corrected_fp),
               names_to = "type", values_to = "count") %>%
  mutate(type = factor(type, levels = c("false_positives", "corrected_fp")))

p_fp <- ggplot(fp_long, aes(x = error_rate, y = count,
                            colour = factor(read_length), linetype = type)) +
  geom_line(linewidth = 0.9 * SC) +
  geom_point(size = 2.0 * SC) +
  scale_x_log10(breaks = c(0.1, 0.5, 1, 5)) +
  scale_colour_manual(values = HCOL$readlen, name = LAB$readlen) +
  scale_linetype_manual(
    breaks = c("false_positives", "corrected_fp"),
    values = c("22", "solid"),
    labels = c("Raw", "Corrected"),
    name = NULL) +
  labs(x = LAB$error, y = LAB$fp_n) +
  guides(linetype = guide_legend(
    override.aes = list(linewidth = 0.5, colour = "grey20"))) +
  herald_theme(W) +
  theme(legend.position = "right",
        legend.key.width = unit(1.1, "cm"))

print(p_fp)

model_raw <- glm(cbind(false_positives, normal_reads - false_positives) ~
                   log10(error_rate) + read_length,
                 family = binomial, data = fp_clean)
summary(model_raw)

model_corrected <- glm(cbind(corrected_fp, normal_reads - corrected_fp) ~
                         log10(error_rate) + read_length,
                       family = binomial, data = fp_clean)
summary(model_corrected)

# =============================================================================
# RECALL x SPIKE RATE
# =============================================================================
spike_data <- read_excel(path, sheet = "spike rate")
glimpse(spike_data)

spike_clean <- spike_data %>%
  mutate(spike_rate     = `spike rate`,
         hgt_reads      = `hgt reads`,
         true_positives = `# true positives`) %>%
  select(spike_rate, hgt_reads, true_positives) %>%
  rowwise() %>%
  mutate(recall  = true_positives / hgt_reads,
         ci_low  = binom.test(true_positives, hgt_reads)$conf.int[1],
         ci_high = binom.test(true_positives, hgt_reads)$conf.int[2]) %>%
  ungroup()

spike_clean

model_spike <- glm(cbind(true_positives, hgt_reads - true_positives) ~
                     log10(spike_rate), family = binomial, data = spike_clean)
summary(model_spike)

p_spike <- ggplot(spike_clean, aes(x = spike_rate, y = recall)) +
  geom_line(colour = HCOL$data, linewidth = 0.9 * SC) +
  geom_errorbar(aes(ymin = ci_low, ymax = ci_high), width = 0.02,
                colour = HCOL$data, linewidth = 0.5 * SC) +
  geom_point(colour = HCOL$data, size = 2.2 * SC) +
  scale_x_log10(breaks = spike_clean$spike_rate) +
  scale_y_continuous(labels = scales::percent, limits = c(0, 1)) +
  labs(x = LAB$spike, y = LAB$recall) +
  herald_theme(W)

print(p_spike)

# =============================================================================
# RECALL x INSERT LENGTH
# =============================================================================
insert_len_data <- read_excel(path, sheet = "var insert length", skip = 1)
glimpse(insert_len_data)

insert_len_clean <- insert_len_data %>%
  filter(!is.na(`error rate`)) %>%   # the total row is the only NA error rate
  mutate(insert_length  = as.numeric(str_remove(insert, "kb")),
         hgt_reads      = `hgt reads`,
         true_positives = `# true positives`) %>%
  select(insert_length, hgt_reads, true_positives) %>%
  rowwise() %>%
  mutate(recall  = true_positives / hgt_reads,
         ci_low  = binom.test(true_positives, hgt_reads)$conf.int[1],
         ci_high = binom.test(true_positives, hgt_reads)$conf.int[2]) %>%
  ungroup() %>%
  mutate(label_y = ci_high + 0.05)

insert_len_clean

overall_table <- rbind(insert_len_clean$true_positives,
                       insert_len_clean$hgt_reads - insert_len_clean$true_positives)
colnames(overall_table) <- insert_len_clean$insert_length
chisq.test(overall_table)

pairwise.prop.test(x = insert_len_clean$true_positives,
                   n = insert_len_clean$hgt_reads,
                   p.adjust.method = "holm")

p_insert_len <- ggplot(insert_len_clean, aes(x = insert_length, y = recall)) +
  geom_pointrange(aes(ymin = ci_low, ymax = ci_high),
                  colour = HCOL$data, size = 0.7 * SC,
                  linewidth = 0.6 * SC) +
  geom_text(aes(y = label_y, label = paste0(true_positives, "/", hgt_reads)),
            size = 3.0 * SC, colour = "grey30", family = "lato") +
  # bracket 1: 1.6 vs 3
  annotate("segment", x = 1.6, xend = 3,   y = 1.05, yend = 1.05) +
  annotate("segment", x = 1.6, xend = 1.6, y = 1.05, yend = 1.02) +
  annotate("segment", x = 3,   xend = 3,   y = 1.05, yend = 1.02) +
  annotate("text", x = 2.3, y = 1.08, label = "***",
           size = 4.0 * SC, family = "lato") +
  # bracket 2: the plateau, 3 to 8.6
  annotate("segment", x = 3,   xend = 8.6, y = 1.05, yend = 1.05) +
  annotate("segment", x = 3,   xend = 3,   y = 1.05, yend = 1.02) +
  annotate("segment", x = 8.6, xend = 8.6, y = 1.05, yend = 1.02) +
  annotate("text", x = 5.8, y = 1.08, label = "ns",
           size = 3.5 * SC, fontface = "italic", family = "lato") +
  # bracket 3: 8.6 vs 10
  annotate("segment", x = 8.6, xend = 10,  y = 1.05, yend = 1.05) +
  annotate("segment", x = 8.6, xend = 8.6, y = 1.05, yend = 1.02) +
  annotate("segment", x = 10,  xend = 10,  y = 1.05, yend = 1.02) +
  annotate("text", x = 9.3, y = 1.08, label = "*",
           size = 4.0 * SC, family = "lato") +
  scale_y_continuous(labels = scales::percent, limits = c(0, 1.18)) +
  scale_x_continuous(breaks = insert_len_clean$insert_length) +
  labs(x = LAB$insert, y = LAB$recall) +
  herald_theme(W)

print(p_insert_len)

# =============================================================================
# RECALL x INSERT PLACEMENT
# =============================================================================
insert_place_data <- read_excel(path, sheet = "var insert placement", skip = 1)
glimpse(insert_place_data)

insert_place_clean <- insert_place_data %>%
  filter(`read length` != "total") %>%
  mutate(left_flank     = as.numeric(str_remove(`left host flank`, "kb")),
         hgt_reads      = `hgt reads`,
         true_positives = `# true positives`) %>%
  select(left_flank, hgt_reads, true_positives) %>%
  rowwise() %>%
  mutate(recall  = true_positives / hgt_reads,
         ci_low  = binom.test(true_positives, hgt_reads)$conf.int[1],
         ci_high = binom.test(true_positives, hgt_reads)$conf.int[2]) %>%
  ungroup()

insert_place_clean

overall_table_place <- rbind(
  insert_place_clean$true_positives,
  insert_place_clean$hgt_reads - insert_place_clean$true_positives)
colnames(overall_table_place) <- insert_place_clean$left_flank
chisq_place <- chisq.test(overall_table_place)
chisq_place

pairwise.prop.test(x = insert_place_clean$true_positives,
                   n = insert_place_clean$hgt_reads,
                   p.adjust.method = "holm")

pooled_rate <- sum(insert_place_clean$true_positives) /
  sum(insert_place_clean$hgt_reads)
pooled_ci <- binom.test(sum(insert_place_clean$true_positives),
                        sum(insert_place_clean$hgt_reads))$conf.int

p_insert_place <- ggplot(insert_place_clean, aes(x = left_flank, y = recall)) +
  annotate("rect", xmin = -1, xmax = 13,
           ymin = pooled_ci[1], ymax = pooled_ci[2],
           fill = HCOL$data, alpha = 0.08) +
  annotate("segment", x = -1, xend = 13, y = pooled_rate, yend = pooled_rate,
           linetype = "dashed", colour = HCOL$data, linewidth = 0.5 * SC) +
  geom_pointrange(aes(ymin = ci_low, ymax = ci_high),
                  colour = HCOL$data, size = 0.7 * SC,
                  linewidth = 0.6 * SC) +
  geom_text(aes(y = ci_high + 0.03,
                label = paste0(true_positives, "/", hgt_reads)),
            size = 3.0 * SC, colour = "grey30", family = "lato") +
  annotate("text", x = 6, y = 1.05,
           label = paste0("p = ", round(chisq_place$p.value, 2), " (ns)"),
           size = 3.5 * SC, fontface = "italic", colour = "grey40",
           family = "lato") +
  scale_y_continuous(labels = scales::percent, limits = c(0, 1.18)) +
  scale_x_continuous(breaks = insert_place_clean$left_flank,
                     limits = c(-1, 13)) +
  labs(x = LAB$flank, y = LAB$recall) +
  herald_theme(W)

print(p_insert_place)

# =============================================================================
# ASSEMBLE
# =============================================================================
fig <- (p_fp + p_spike) / (p_insert_len + p_insert_place) +
  plot_annotation(tag_levels = "A") &
  herald_tag(W)

print(fig)

herald_save(fig, "herald_module1_figure", W, H, outdir)