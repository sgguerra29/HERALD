#!/usr/bin/env Rscript
# =============================================================================
# HERALD_theme.R — shared styling for every HERALD figure.
#
# Source this at the top of each figure script:
#
#   source("HERALD_theme.R")
#   W <- 12.5; H <- 8.6          # this figure's output size in inches
#   SC <- herald_sc(W)           # font/size scale factor for that width
#
# Then:
#   - replace theme_classic(base_size = ..., base_family = "lato")
#     with herald_theme(W)
#   - multiply every hard-coded geom size by SC:
#       size = 3.2      ->  size = 3.2 * SC
#       linewidth = 0.9 ->  linewidth = 0.9 * SC
#   - use herald_table_theme(W, fill_matrix, face_matrix, mult) for tableGrob
#   - use LAB$<name> for axis titles instead of literal strings
# =============================================================================

suppressPackageStartupMessages({
  library(ggplot2)
  library(gridExtra)
  library(grid)
  library(showtext)
})

if (!"lato" %in% sysfonts::font_families()) {
  font_add_google(name = "Lato", family = "lato")
}
showtext_auto()

# ── Sizing ───────────────────────────────────────────────────────────────────
HERALD_MIN_PT    <- 6.5      # smallest acceptable printed point size
HERALD_DISPLAY_W <- 6.5    # inches the figure occupies on the printed page
HERALD_REF_BASE  <- 10.5   # the base_size the scripts were originally written at

herald_base <- function(w) {
  # base_size such that axis text (0.8 * base) prints at >= HERALD_MIN_PT
  HERALD_MIN_PT * (w / HERALD_DISPLAY_W) / 0.8
}

herald_sc <- function(w) herald_base(w) / HERALD_REF_BASE

herald_theme <- function(w, classic = TRUE) {
  b <- herald_base(w)
  base <- if (classic) theme_classic(base_size = b, base_family = "lato")
  else         theme_bw(base_size = b, base_family = "lato")
  base + theme(
    axis.text         = element_text(colour = "black"),
    legend.background = element_blank(),
    legend.key.height = unit(1.35 * b, "pt"),
    legend.key.width  = unit(0.90 * b, "pt"),
    legend.title      = element_text(size = 0.85 * b),
    legend.text       = element_text(size = 0.80 * b),
    plot.title        = element_text(hjust = 0.5),
    plot.subtitle     = element_text(size = 0.80 * b, colour = "grey35"),
    plot.caption      = element_text(hjust = 0, size = 0.75 * b,
                                     colour = "grey40")
  )
}

# Panel tags. Apply with:  & herald_tag(W)
herald_tag <- function(w) {
  theme(plot.tag = element_text(size = 1.30 * herald_base(w),
                                family = "lato", face = "bold"))
}
herald_table_theme <- function(w, fill_matrix, face_matrix = NULL,
                               mult = 0.95) {
  b <- mult * herald_base(w)
  fg <- list(hjust = 0.5, x = 0.5)
  if (!is.null(face_matrix)) fg$fontface <- as.vector(face_matrix)
  ttheme_minimal(
    base_size = b, base_family = "lato",
    core = list(fg_params = fg,
                bg_params = list(fill = as.vector(fill_matrix),
                                 col = "grey75", lwd = 0.5)),
    colhead = list(fg_params = list(fontface = 2, hjust = 0.5, x = 0.5),
                   bg_params = list(fill = HCOL$fill_head,
                                    col = "grey75", lwd = 0.5)))
}

# ── Palette ──────────────────────────────────────────────────────────────────
HCOL <- list(
  # schematics only
  host        = "#782360",   # host / acceptor sequence
  foreign     = "#8AC223",   # foreign / donor insert
  readlen     = c(`10` = "#501DA3", `15` = "#2a6f97",
                  `20` = "#8AC223", `25` = "#ffed7d"),
  ramp_lo     = "#FFFDF5",
  ramp_hi     = "#782360",
  data        = "#256bb3",
  cat2        = c("#2a6f97", "#7CBAF7"),
  accent      = "#FFC700",
  neutral     = "#BDBDBD",
  
  # table cell fills
  fill_head   = "#EDEDED",
  fill_hilite = "#EDE3EA",   # pale tint of `host`
  fill_total  = "#DCD0D8"    # slightly deeper, for a Total row
)

ramp_n <- function(n) {
  grDevices::colorRampPalette(c("#E7D9E2", HCOL$ramp_hi))(n)
}

ramp_blue <- function(n) {
  grDevices::colorRampPalette(c("#103459", "#BBDBFA"))(n)
}

herald_fill_scale <- function(name = LAB$recall_short) {
  scale_fill_gradient(low = HCOL$ramp_lo, high = HCOL$ramp_hi,
                      limits = c(0, 1),
                      labels = scales::percent_format(accuracy = 1),
                      name = name)
}

# ── Canonical labels ─────────────────────────────────────────────────────────
LAB <- list(
  recall        = "Recall",
  recall_short  = "Recall",
  insert        = "Insert length (kb)",
  mmf           = "Minimum matched fraction",
  readlen       = "Read length (kb)",
  frags         = "Fragments per read",
  flank         = "Left host flank (kb)",
  error         = "Error rate (%)",
  spike         = "Spike rate",
  fp_n          = "False positives (n)",
  cand_rate     = "Candidate reads per 100,000 reads",
  decile        = "Read-length decile",
  decile_lo     = "Decile lower bound (kb)",
  events        = "Events detected"
)

# ── Unit helper ──────────────────────────────────────────────────────────────
kb_lab <- function(bp) {
  v <- bp / 1000
  ifelse(v == round(v), sprintf("%.0f", v), sprintf("%.1f", v))
}

# Factor with levels in numeric order, so 10 does not sort before 2.
kb_factor <- function(bp) {
  u <- sort(unique(bp))
  factor(kb_lab(bp), levels = kb_lab(u))
}

# ── Save helper ──────────────────────────────────────────────────────────────
herald_save <- function(fig, path_noext, w, h, outdir = ".") {
  showtext::showtext_opts(dpi = 400)
  ggsave(file.path(outdir, paste0(path_noext, ".png")),
         fig, width = w, height = h, dpi = 400, bg = "white")
  showtext::showtext_opts(dpi = 300)
  ggsave(file.path(outdir, paste0(path_noext, ".pdf")),
         fig, width = w, height = h, bg = "white")
  showtext::showtext_opts(dpi = 96)
  message(sprintf("wrote %s.{png,pdf} at %.1f x %.1f in (base_size %.1f)",
                  path_noext, w, h, herald_base(w)))
}

# ── Report ───────────────────────────────────────────────────────────────────
message(sprintf(
  "HERALD_theme.R loaded. Design base_size by width (for %.0f pt minimum at %.1f in display):",
  HERALD_MIN_PT, HERALD_DISPLAY_W))
for (w in c(6.5, 10, 11, 12.5, 13)) {
  message(sprintf("  %5.1f in  ->  base_size %5.1f   scale x%.2f",
                  w, herald_base(w), herald_sc(w)))
}