#!/usr/bin/env Rscript

# ==============================================================================
# Script: 04_analysis/04_control_ordering.R
#
# Project: Aridity and timescale bound the separability of vegetation water stress
#
# Purpose (NUMBERS ONLY — exploratory/defensive, not a paper figure):
#   Quantify the ORDERING / LOCATION of the dominant control along the aridity
#   gradient, from the FDR dominance classification (same as 42_dominance_aridity
#   _crosstab.R). Focus on the SOIL vs JOINT pair; AED reported separately.
#
#   *** Conceptual note: this is EXTENT / LOCATION of the dominant control by
#   aridity, NOT correlation STRENGTH. It is a different metric from the aridity
#   figure (|rho*|), where AED peaks in semi-arid. Do not conflate the two. ***
#
#   A) Per category (soil, AED, joint): area-weighted AI distribution of its
#      cell x month members — median, mean, p10/p25/p75/p90 — plus a simple
#      overlap diagnostic between soil and joint.
#   B) Per Zomer aridity class: % area of each category, the most frequent
#      control per class, and the gradient crossover where JOINT overtakes SOIL
#      (fine AI bins).
#   C) AED alone: its AI distribution (from A) and the aridity classes where its
#      UNIQUE dominant control (code 3) concentrates.
#
#   Aridity classes + AI raster = Zomer (2022), identical to 42_* and Fig. b.
#   Area-weighted by cellSize; pooled over the 12 months.
#
# Inputs (READ-ONLY; mtimes confirmed unchanged at the end):
#   outputs/varpart_global/varpart_signif_global_2blocks.nc   (subds: shared)
#   outputs/varpart_global/fdr_adjusted_pvalues.nc            (p_full/soil/demand_adj)
#   data/processed/aridity/ai_1982_2022_period.nc             (continuous AI)
#
# Outputs:
#   outputs/tables/control_ordering_ai_distribution.csv
#   outputs/tables/control_ordering_by_aridity_class.csv
#
# Run from repo root:
#   Rscript 04_analysis/04_control_ordering.R
#
# Dependencies: terra
# ==============================================================================

suppressPackageStartupMessages({
  library(terra)
})

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
dir.create(tables_dir, recursive = TRUE, showWarnings = FALSE)
out_dist <- file.path(tables_dir, "control_ordering_ai_distribution.csv")
out_cls  <- file.path(tables_dir, "control_ordering_by_aridity_class.csv")
for (f in c(out_dist, out_cls))
  if (file.exists(f)) stop("Output already exists (not overwriting): ", f, call. = FALSE)
for (f in c(file_nc, file_adj, file_ai))
  if (!file.exists(f)) stop("Required input not found: ", f, call. = FALSE)
mt_before <- file.info(c(file_nc, file_adj, file_ai))$mtime

ALPHA <- 0.05
# Dominance codes: 1 = no signal, 2 = soil, 3 = AED (demand), 4 = joint (coupled)
cat_codes  <- c(soil = 2L, AED = 3L, joint = 4L)
ai_breaks  <- c(-Inf, 0.03, 0.20, 0.50, 0.65, Inf)
ai_labels  <- c("Hyper-arid", "Arid", "Semi-arid", "Dry sub-humid", "Humid")

classify_set_na <- function(r) { r[r == -9999] <- NA; r }
classify_month <- function(pf, ps, pd, sh) {
  out <- rep(NA_integer_, length(sh)); valid <- !is.na(sh)
  full_sig <- valid & !is.na(pf) & pf < ALPHA
  soil_sig <- !is.na(ps) & ps < ALPHA
  dem_sig  <- !is.na(pd) & pd < ALPHA
  out[valid & !full_sig]               <- 1L
  out[full_sig &  soil_sig & !dem_sig] <- 2L
  out[full_sig & !soil_sig &  dem_sig] <- 3L
  out[full_sig & ((!soil_sig & !dem_sig) | (soil_sig & dem_sig))] <- 4L
  out
}
wquantile <- function(x, w, probs) {
  ok <- !is.na(x) & !is.na(w) & w > 0; x <- x[ok]; w <- w[ok]
  if (length(x) == 0) return(rep(NA_real_, length(probs)))
  o <- order(x); x <- x[o]; w <- w[o]
  cw <- (cumsum(w) - 0.5 * w) / sum(w)
  approx(cw, x, xout = probs, rule = 2, ties = "ordered")$y
}
wmean <- function(x, w) { ok <- !is.na(x) & !is.na(w); sum(x[ok] * w[ok]) / sum(w[ok]) }

# ------------------------------------------------------------------------------
# Build monthly dominance + AI vectors (pooled over 12 months, area-weighted)
# ------------------------------------------------------------------------------
message("Classifying dominance (FDR) and pooling with AI...")
r_sh <- classify_set_na(rast(file_nc,  subds = "shared"))
r_pf <- classify_set_na(rast(file_adj, subds = "p_full_adj"))
r_ps <- classify_set_na(rast(file_adj, subds = "p_soil_adj"))
r_pd <- classify_set_na(rast(file_adj, subds = "p_demand_adj"))

cls_stack <- rast(r_sh); values(cls_stack) <- NA
for (m in 1:12)
  values(cls_stack[[m]]) <- classify_month(values(r_pf[[m]]), values(r_ps[[m]]),
                                            values(r_pd[[m]]), values(r_sh[[m]]))
if (is.na(crs(cls_stack)) || crs(cls_stack) == "") crs(cls_stack) <- "EPSG:4326"

ai <- rast(file_ai)[[1]]; if (is.na(crs(ai))) crs(ai) <- "EPSG:4326"
if (xmax(ai) > 180 + 1e-6) ai <- rotate(ai)
ai_r  <- resample(ai, cls_stack[[1]], method = "bilinear")
ai_cl <- classify(ai_r, matrix(c(head(ai_breaks, -1), tail(ai_breaks, -1), 1:5), ncol = 3),
                  right = FALSE)

area_v <- as.vector(values(cellSize(cls_stack[[1]], unit = "km", mask = FALSE)))
ai_vec <- as.vector(values(ai_r))
cl_vec <- as.vector(values(ai_cl))

code_all <- unlist(lapply(1:12, function(m) as.vector(values(cls_stack[[m]]))))
area_all <- rep(area_v, 12); ai_all <- rep(ai_vec, 12); cl_all <- rep(cl_vec, 12)
keep <- !is.na(code_all) & !is.na(ai_all)
code_all <- code_all[keep]; area_all <- area_all[keep]
ai_all <- ai_all[keep]; cl_all <- cl_all[keep]

# ------------------------------------------------------------------------------
# A) AI distribution per category (area-weighted)
# ------------------------------------------------------------------------------
cat("\n=== A) AI distribution per dominant control (area-weighted, 12 months) ===\n")
distr <- do.call(rbind, lapply(names(cat_codes), function(nm) {
  k <- cat_codes[[nm]]; sel <- code_all == k
  q <- wquantile(ai_all[sel], area_all[sel], c(.10, .25, .50, .75, .90))
  data.frame(control = nm, n_cellmonths = sum(sel),
             area_Mkm2 = round(sum(area_all[sel]) / 12 / 1e6, 3),
             mean = round(wmean(ai_all[sel], area_all[sel]), 3),
             p10 = round(q[1], 3), p25 = round(q[2], 3), median = round(q[3], 3),
             p75 = round(q[4], 3), p90 = round(q[5], 3))
}))
print(distr, row.names = FALSE)

# Soil vs joint ordering + overlap diagnostics
soil_sel <- code_all == 2L; joint_sel <- code_all == 4L
med_soil  <- wquantile(ai_all[soil_sel],  area_all[soil_sel],  0.5)
med_joint <- wquantile(ai_all[joint_sel], area_all[joint_sel], 0.5)
# fraction of JOINT area drier than soil median, and SOIL area wetter than joint median
f_joint_below_soilmed <- 100 * sum(area_all[joint_sel & ai_all < med_soil]) / sum(area_all[joint_sel])
f_soil_above_jointmed <- 100 * sum(area_all[soil_sel  & ai_all > med_joint]) / sum(area_all[soil_sel])
cat(sprintf("\n  Ordering: soil median AI = %.3f  vs  joint median AI = %.3f  (soil %s joint)\n",
            med_soil, med_joint, if (med_soil < med_joint) "DRIER than" else "NOT drier than"))
cat(sprintf("  Overlap: %.1f%% of JOINT area is drier than the soil median; %.1f%% of SOIL area is wetter than the joint median.\n",
            f_joint_below_soilmed, f_soil_above_jointmed))
write.csv(distr, out_dist, row.names = FALSE)

# ------------------------------------------------------------------------------
# B) Per Zomer aridity class: % area of each control + most frequent + crossover
# ------------------------------------------------------------------------------
cat("\n=== B) Control composition per aridity class (% of the class analysed area) ===\n")
area_mat <- tapply(area_all, list(factor(cl_all, levels = 1:5), factor(code_all, levels = 1:4)), sum)
area_mat[is.na(area_mat)] <- 0
row_tot <- rowSums(area_mat)
pct <- sweep(area_mat, 1, ifelse(row_tot > 0, row_tot, NA), "/") * 100

cls_tab <- data.frame(
  ai_class = ai_labels,
  ai_range = c("<0.03", "0.03-0.20", "0.20-0.50", "0.50-0.65", ">=0.65"),
  no_signal = round(pct[, 1], 1),
  soil = round(pct[, 2], 1), AED = round(pct[, 3], 1), joint = round(pct[, 4], 1),
  area_Mkm2 = round(row_tot / 12 / 1e6, 3), row.names = NULL
)
# Most frequent control among the three (soil/AED/joint) per class
ctrl_only <- pct[, 2:4, drop = FALSE]; colnames(ctrl_only) <- c("soil", "AED", "joint")
cls_tab$top_control <- c("soil", "AED", "joint")[apply(ctrl_only, 1, which.max)]
print(cls_tab, row.names = FALSE)
write.csv(cls_tab, out_cls, row.names = FALSE)

# Crossover soil -> joint along the continuous gradient (fine AI bins)
cat("\n--- Soil -> Joint crossover along the AI gradient (fine bins, width 0.02) ---\n")
fb <- seq(0, 1, by = 0.02)
mids <- head(fb, -1) + 0.01
bin  <- cut(ai_all, breaks = fb, right = FALSE, include.lowest = TRUE)
soil_b  <- tapply(area_all[soil_sel],  bin[soil_sel],  sum)
joint_b <- tapply(area_all[joint_sel], bin[joint_sel], sum)
soil_b[is.na(soil_b)] <- 0; joint_b[is.na(joint_b)] <- 0
dif <- as.numeric(joint_b) - as.numeric(soil_b)              # >0 where joint leads
tot_b <- as.numeric(soil_b) + as.numeric(joint_b)
valid_b <- tot_b > 0
cross_idx <- which(dif > 0 & c(FALSE, head(dif, -1) <= 0) & valid_b)
if (length(cross_idx) > 0) {
  cat(sprintf("  Joint first overtakes soil (area) at AI ~ %.2f\n", mids[cross_idx[1]]))
} else {
  cat("  No clean soil->joint crossover found in the AI bins.\n")
}
# also which control leads per Zomer class (soil vs joint only)
lead_sj <- ifelse(pct[, 2] >= pct[, 4], "soil", "joint")
cat("  Leading of {soil,joint} per Zomer class: ",
    paste(sprintf("%s=%s", ai_labels, lead_sj), collapse = " | "), "\n")

# ------------------------------------------------------------------------------
# C) AED alone: AI distribution (from A) + where its unique control concentrates
# ------------------------------------------------------------------------------
cat("\n=== C) AED (unique dominant, code 3) — where it falls on the gradient ===\n")
cat("  (Location only; driven by boreal-spring seasonality. NOT the |rho*| strength peak.)\n")
aed_sel <- code_all == 3L
aed_area_by_cls <- tapply(area_all[aed_sel], factor(cl_all[aed_sel], levels = 1:5), sum)
aed_area_by_cls[is.na(aed_area_by_cls)] <- 0
aed_share <- 100 * as.numeric(aed_area_by_cls) / sum(aed_area_by_cls)   # % of AED area per class
aed_tab <- data.frame(ai_class = ai_labels,
                      pct_of_AED_area = round(aed_share, 1),
                      pct_within_class = round(pct[, 3], 1))
print(aed_tab, row.names = FALSE)
q_aed <- wquantile(ai_all[aed_sel], area_all[aed_sel], c(.25, .5, .75))
cat(sprintf("  AED AI: median %.3f (p25 %.3f, p75 %.3f); concentrates in class(es): %s\n",
            q_aed[2], q_aed[1], q_aed[3],
            paste(ai_labels[order(aed_share, decreasing = TRUE)][1:2], collapse = ", ")))

# ------------------------------------------------------------------------------
# File-protection confirmation
# ------------------------------------------------------------------------------
cat("\n=== File protection confirmation (inputs READ-ONLY) ===\n")
inputs_all <- c(file_nc, file_adj, file_ai)
mt_after <- file.info(inputs_all)$mtime
for (i in seq_along(inputs_all))
  cat(sprintf("  %-40s : %s\n", basename(inputs_all[i]),
              if (identical(mt_before[i], mt_after[i])) "UNCHANGED" else "*** CHANGED ***"))
cat("\nWrote:\n  ", out_dist, "\n  ", out_cls, "\n")
