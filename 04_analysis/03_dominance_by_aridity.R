#!/usr/bin/env Rscript

# ==============================================================================
# Script: 04_analysis/03_dominance_by_aridity.R
#
# Project: Aridity and timescale bound the separability of vegetation water stress
#
# Purpose (EXPLORATORY — not a paper figure yet):
#   Cross the FDR dominance classification (4 categories) with the aridity
#   gradient, area-weighted, pooled over the 12 months, to check whether the
#   significant controls concentrate in the same semi-arid band (AI ~ 0.28-0.34)
#   as the five indicators.
#
#   Analysis 1 — control composition per aridity class:
#     For each Zomer aridity class, the % of (cellSize-weighted) area in each of
#     the 4 categories (no-signal / soil / demand / joint). Reports the class
#     where the WITH-SIGNAL fraction (non-"no-signal") is maximal.
#
#   Analysis 2 — position of each control along the gradient:
#     For each with-signal category (soil, demand, joint), the area-weighted AI
#     distribution of its cells (median, p25-p75, plus p5/p95).
#
#   Plus a quick exploratory plot: category shares vs AI over fine bins (0.05).
#
#   Classification is IDENTICAL to scripts 42/43/44 (ALPHA = 0.05, classify_month)
#   so categories match the dominance figures. Aridity bins follow Zomer et al.
#   (2022), the same thresholds as 01_prepare_data/07_classify_aridity.R.
#
# Inputs (READ-ONLY; mtimes recorded and confirmed unchanged at the end):
#   outputs/varpart_global/varpart_signif_global_2blocks.nc   (subds: shared)
#   outputs/varpart_global/fdr_adjusted_pvalues.nc            (p_full/soil/demand_adj)
#   data/processed/aridity/ai_1982_2022_period.nc             (mean annual AI)
#
# Outputs (the script stops if they already exist):
#   outputs/tables/dominance_aridity_composition.csv
#   outputs/tables/dominance_aridity_ai_position.csv
#   outputs/exploratory/dominance_vs_aridity_finebins.png
#
# Run from repo root (AFTER 03_variance_partitioning/02_apply_fdr.R):
#   Rscript 04_analysis/03_dominance_by_aridity.R
#
# Dependencies: terra, ggplot2, dplyr
# ==============================================================================

suppressPackageStartupMessages({
  library(terra)
  library(ggplot2)
  library(dplyr)
})

# ------------------------------------------------------------------------------
# Paths
# ------------------------------------------------------------------------------
config_file <- file.path("R", "config.R")
if (file.exists(config_file)) source(config_file)

file_nc <- if (exists("paths") && !is.null(paths$varpart_signif)) {
  paths$varpart_signif
} else {
  file.path("outputs", "varpart_global", "varpart_signif_global_2blocks.nc")
}
file_adj <- file.path("outputs", "varpart_global", "fdr_adjusted_pvalues.nc")
file_ai <- if (exists("paths") && !is.null(paths$aridity_index)) {
  paths$aridity_index
} else {
  file.path("data", "processed", "aridity", "ai_1982_2022_period.nc")
}

tables_dir <- if (exists("paths") && !is.null(paths$tables)) {
  paths$tables
} else {
  file.path("outputs", "tables")
}
explor_dir <- file.path("outputs", "exploratory")
dir.create(tables_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(explor_dir, recursive = TRUE, showWarnings = FALSE)

out_comp <- file.path(tables_dir, "dominance_aridity_composition.csv")
out_posi <- file.path(tables_dir, "dominance_aridity_ai_position.csv")
out_plot <- file.path(explor_dir, "dominance_vs_aridity_finebins.png")

for (f in c(out_comp, out_posi, out_plot))
  if (file.exists(f)) stop("Output already exists (not overwriting): ", f, call. = FALSE)

for (f in c(file_nc, file_adj, file_ai))
  if (!file.exists(f)) stop("Required input not found: ", f, call. = FALSE)

mt_before <- file.info(c(file_nc, file_adj, file_ai))$mtime

# ------------------------------------------------------------------------------
# Settings — dominance (as in 42/43/44) and aridity (Zomer 2022, as in 06_*)
# ------------------------------------------------------------------------------
ALPHA <- 0.05

cat_levels <- c("sin_senal", "soil", "demand", "coupled")   # codes 1..4
cat_labels <- c("No significant signal", "Soil moisture", "AED", "Joint (non-separable)")

# Zomer et al. (2022) thresholds; left-closed [lower, upper). Same as 06_*.
ai_breaks <- c(-Inf, 0.03, 0.20, 0.50, 0.65, Inf)
ai_class_labels <- c("Hyper-arid", "Arid", "Semi-arid", "Dry sub-humid", "Humid")

fine_step <- 0.05                                           # fine AI bins for the plot only

# ------------------------------------------------------------------------------
# Helpers (dominance logic)
# ------------------------------------------------------------------------------
classify_set_na <- function(r) { r[r == -9999] <- NA; r }

classify_month <- function(pf, ps, pd, sh) {
  out   <- rep(NA_integer_, length(sh))
  valid <- !is.na(sh)
  full_sig <- valid & !is.na(pf) & pf < ALPHA
  soil_sig <- !is.na(ps) & ps < ALPHA
  dem_sig  <- !is.na(pd) & pd < ALPHA
  out[valid & !full_sig]                 <- 1L
  out[full_sig &  soil_sig & !dem_sig]   <- 2L
  out[full_sig & !soil_sig &  dem_sig]   <- 3L
  out[full_sig & ((!soil_sig & !dem_sig) | ( soil_sig &  dem_sig))] <- 4L
  out
}

# Area-weighted quantile (type-7-like).
wquantile <- function(x, w, probs) {
  ok <- !is.na(x) & !is.na(w) & w > 0
  x <- x[ok]; w <- w[ok]
  if (length(x) == 0) return(rep(NA_real_, length(probs)))
  o <- order(x); x <- x[o]; w <- w[o]
  cw <- (cumsum(w) - 0.5 * w) / sum(w)
  approx(cw, x, xout = probs, rule = 2, ties = "ordered")$y
}

# ------------------------------------------------------------------------------
# Read inputs and build monthly dominance classification (codes 1..4)
# ------------------------------------------------------------------------------
message("Reading dominance inputs and classifying (FDR, 4 categories, 12 months)...")
r_sh <- classify_set_na(rast(file_nc,  subds = "shared"))
r_pf <- classify_set_na(rast(file_adj, subds = "p_full_adj"))
r_ps <- classify_set_na(rast(file_adj, subds = "p_soil_adj"))
r_pd <- classify_set_na(rast(file_adj, subds = "p_demand_adj"))

cls_stack <- rast(r_sh); values(cls_stack) <- NA
for (m in 1:12)
  values(cls_stack[[m]]) <- classify_month(values(r_pf[[m]]), values(r_ps[[m]]),
                                            values(r_pd[[m]]), values(r_sh[[m]]))

if (is.na(crs(cls_stack)) || crs(cls_stack) == "") crs(cls_stack) <- "EPSG:4326"

# ------------------------------------------------------------------------------
# Aridity: resample AI (0.5 deg) onto the dominance grid (0.1 deg), then classify
# ------------------------------------------------------------------------------
message("Resampling AI to the dominance grid and classifying (Zomer 2022)...")
ai <- rast(file_ai)[[1]]
if (is.na(crs(ai)) || crs(ai) == "") crs(ai) <- "EPSG:4326"
if (xmax(ai) > 180 + 1e-6) ai <- rotate(ai)                # ensure -180/180 like GLEAM grid

ai_r  <- resample(ai, cls_stack[[1]], method = "bilinear") # continuous AI on 0.1 deg grid
ai_cl <- classify(ai_r, matrix(c(head(ai_breaks, -1), tail(ai_breaks, -1), 1:5),
                               ncol = 3), right = FALSE)    # class codes 1..5

# ------------------------------------------------------------------------------
# Vectors (area-weighted, pooled over 12 months)
# ------------------------------------------------------------------------------
area   <- cellSize(cls_stack[[1]], unit = "km", mask = FALSE)
area_v <- as.vector(values(area))
ai_v   <- as.vector(values(ai_r))
aicl_v <- as.vector(values(ai_cl))

code_all <- unlist(lapply(1:12, function(m) as.vector(values(cls_stack[[m]]))))
area_all <- rep(area_v,  12)
ai_all   <- rep(ai_v,   12)
aicl_all <- rep(aicl_v, 12)

keep <- !is.na(code_all) & !is.na(aicl_all)
code_all <- code_all[keep]; area_all <- area_all[keep]
ai_all   <- ai_all[keep];   aicl_all <- aicl_all[keep]

# ------------------------------------------------------------------------------
# Analysis 1 — control composition per aridity class (row-% of area)
# ------------------------------------------------------------------------------
message("\n=== Analysis 1: control composition per aridity class (area-weighted %) ===")
area_mat <- tapply(area_all, list(factor(aicl_all, levels = 1:5),
                                  factor(code_all, levels = 1:4)), sum)
area_mat[is.na(area_mat)] <- 0
row_tot <- rowSums(area_mat)
pct_mat <- sweep(area_mat, 1, ifelse(row_tot > 0, row_tot, NA), "/") * 100

with_signal <- rowSums(pct_mat[, 2:4, drop = FALSE])       # soil+demand+joint

comp <- data.frame(
  ai_class   = ai_class_labels,
  ai_range   = c("<0.03", "0.03-0.20", "0.20-0.50", "0.50-0.65", ">=0.65"),
  no_signal      = round(pct_mat[, 1], 1),
  soil_moisture  = round(pct_mat[, 2], 1),
  AED            = round(pct_mat[, 3], 1),
  joint          = round(pct_mat[, 4], 1),
  with_signal    = round(with_signal, 1),
  area_Mkm2      = round(row_tot / 12 / 1e6, 3),           # /12: back to single-month area
  row.names = NULL
)
print(comp, row.names = FALSE)

peak_i <- which.max(ifelse(is.na(with_signal), -Inf, with_signal))
message(sprintf("\n  -> WITH-SIGNAL fraction is MAXIMAL in: %s (%.1f%%)",
                ai_class_labels[peak_i], with_signal[peak_i]))
write.csv(comp, out_comp, row.names = FALSE)

# ------------------------------------------------------------------------------
# Analysis 2 — AI distribution of each with-signal control (area-weighted)
# ------------------------------------------------------------------------------
message("\n=== Analysis 2: AI position of each control (area-weighted) ===")
pos_rows <- lapply(2:4, function(k) {
  sel <- code_all == k
  q <- wquantile(ai_all[sel], area_all[sel], c(0.05, 0.25, 0.50, 0.75, 0.95))
  data.frame(category = cat_labels[k],
             n_cellmonths = sum(sel),
             AI_p5 = round(q[1], 3), AI_p25 = round(q[2], 3),
             AI_median = round(q[3], 3),
             AI_p75 = round(q[4], 3), AI_p95 = round(q[5], 3))
})
posi <- do.call(rbind, pos_rows)
print(posi, row.names = FALSE)
message("\n  (Semi-arid convergence band of the five indicators: AI ~ 0.28-0.34)")
write.csv(posi, out_posi, row.names = FALSE)

# ------------------------------------------------------------------------------
# Quick exploratory plot — category shares vs AI over fine bins (0.05)
# ------------------------------------------------------------------------------
message("\nBuilding exploratory plot (fine AI bins)...")
fbreaks <- seq(0, 1, by = fine_step)
fbin    <- cut(ai_all, breaks = fbreaks, right = FALSE, include.lowest = TRUE)
fmid    <- head(fbreaks, -1) + fine_step / 2

fine_area <- tapply(area_all, list(fbin, factor(code_all, levels = 1:4)), sum)
fine_area[is.na(fine_area)] <- 0
fine_tot  <- rowSums(fine_area)
fine_pct  <- sweep(fine_area, 1, ifelse(fine_tot > 0, fine_tot, NA), "/") * 100

plot_df <- do.call(rbind, lapply(1:4, function(k) data.frame(
  ai = fmid[seq_len(nrow(fine_pct))], category = cat_labels[k], pct = fine_pct[, k])))
plot_df$category <- factor(plot_df$category, levels = cat_labels)

p <- ggplot(plot_df, aes(ai, pct, colour = category)) +
  annotate("rect", xmin = 0.28, xmax = 0.34, ymin = -Inf, ymax = Inf,
           fill = "grey85", alpha = 0.5) +
  geom_vline(xintercept = c(0.03, 0.20, 0.50, 0.65), colour = "grey80",
             linetype = "22", linewidth = 0.3) +
  geom_line(linewidth = 0.7) +
  scale_colour_manual(values = c("#D4D8DD", "#2166AC", "#D98C00", "#1D9E75"), name = NULL) +
  coord_cartesian(xlim = c(0, 0.9)) +
  labs(x = "Aridity Index (AI)", y = "% of classified area (per AI bin)",
       title = "Dominance composition vs aridity (exploratory, fine bins)") +
  theme_classic(base_size = 9) +
  theme(legend.position = "bottom")
ggsave(out_plot, p, width = 150, height = 95, units = "mm", dpi = 200, bg = "white")
message("  Saved plot: ", out_plot)

# ------------------------------------------------------------------------------
# File-protection confirmation
# ------------------------------------------------------------------------------
message("\n=== File protection confirmation (inputs READ-ONLY) ===")
mt_after <- file.info(c(file_nc, file_adj, file_ai))$mtime
for (i in seq_along(mt_after))
  message(sprintf("  %-45s : %s", basename(c(file_nc, file_adj, file_ai)[i]),
                  if (identical(mt_before[i], mt_after[i])) "UNCHANGED" else "*** CHANGED ***"))

message("\nOutputs written:")
message("  ", out_comp)
message("  ", out_posi)
message("  ", out_plot)
message("\nDone (exploratory cross-tab; nothing overwritten).")
