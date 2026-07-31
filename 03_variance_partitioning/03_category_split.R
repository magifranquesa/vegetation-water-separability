#!/usr/bin/env Rscript

# ==============================================================================
# Script: 03_variance_partitioning/03_category_split.R
#
# Project: Aridity and timescale bound the separability of vegetation water stress
# Purpose:
#   Classify each cell x month, from the FDR-adjusted two-block variance
#   partitioning, into five mutually exclusive categories of hydroclimatic
#   control, area-weighted and stratified by aridity.
#
#   The distinction between categories 4 and 5 is essential: category 4 (neither
#   unique fraction significant) is genuine non-separability of the two drivers,
#   whereas category 5 (both unique fractions significant) is the opposite — each
#   driver carries independently resolvable variance. Reporting them separately
#   avoids merging non-separable and separable cells under one label.
#
# Categories produced (5, mutually exclusive; FDR-adjusted, area-weighted):
#   1 no_signal      full model not significant
#   2 soil           full sig, ONLY soil unique significant
#   3 demand         full sig, ONLY demand unique significant
#   4 coupled_pure   full sig, NEITHER unique significant   <- non-separable
#   5 both_unique    full sig, BOTH unique significant       <- separable
#
# Inputs (READ-ONLY; mtimes confirmed unchanged at the end):
#   outputs/varpart_global/varpart_signif_global_2blocks.nc  (shared, total_r2)
#   outputs/varpart_global/fdr_adjusted_pvalues.nc           (p_full/soil/demand_adj)
#   data/processed/aridity/ai_1982_2022_period.nc            (continuous AI)
#
# Outputs (refuses to overwrite):
#   outputs/tables/joint_split.csv          (5-cat area-weighted shares)
#   outputs/tables/joint_split_by_aridity.csv
#
# Run from repo root (AFTER 03_variance_partitioning/02_apply_fdr.R):
#   Rscript 03_variance_partitioning/03_category_split.R
#
# Dependencies: terra
# ==============================================================================

suppressPackageStartupMessages({
  library(terra)
})

config_file <- file.path("R", "config.R")
if (file.exists(config_file)) source(config_file)

# ── Period window (see R/period.R) ───────────────────────────────────────────
# Must match the window used in 01_varpart.R and 02_apply_fdr.R:
#   VWS_PERIOD=1982-2001 Rscript 03_variance_partitioning/03_category_split.R
source(file.path("R", "period.R"))

tables_dir <- if (exists("paths") && !is.null(paths$tables)) paths$tables else
  file.path("outputs", "tables")
dir.create(tables_dir, recursive = TRUE, showWarnings = FALSE)

file_nc  <- file.path("outputs", "varpart_global",
                      paste0("varpart_signif_global_2blocks", PERIOD_SFX, ".nc"))
file_adj <- file.path("outputs", "varpart_global",
                      paste0("fdr_adjusted_pvalues", PERIOD_SFX, ".nc"))
file_ai  <- if (exists("paths") && !is.null(paths$aridity_index)) paths$aridity_index else
  file.path("data", "processed", "aridity", "ai_1982_2022_period.nc")

out_main <- file.path(tables_dir, paste0("joint_split", PERIOD_SFX, ".csv"))
out_arid <- file.path(tables_dir, paste0("joint_split_by_aridity", PERIOD_SFX, ".csv"))

for (f in c(file_nc, file_adj, file_ai))
  if (!file.exists(f)) stop("Required input not found: ", f, call. = FALSE)
for (f in c(out_main, out_arid))
  if (file.exists(f))
    stop("REFUSING TO OVERWRITE: ", f, "\n  Delete/rename it manually to regenerate.", call. = FALSE)

mt_before <- file.info(c(file_nc, file_adj, file_ai))$mtime

FV    <- -9999
ALPHA <- 0.05

months_index <- c(1, 4, 7, 10)
months_names <- c("January", "April", "July", "October")

CATS <- c("no_signal", "soil", "demand", "coupled_pure", "both_unique")

ai_breaks <- c(-Inf, 0.03, 0.20, 0.50, 0.65, Inf)
ai_labels <- c("Hyper-arid", "Arid", "Semi-arid", "Dry sub-humid", "Humid")

set_na <- function(v) { v[abs(v - FV) < 1] <- NA; v }

wquantile <- function(x, w, probs) {
  ok <- !is.na(x) & !is.na(w) & w > 0; x <- x[ok]; w <- w[ok]
  if (!length(x)) return(rep(NA_real_, length(probs)))
  o <- order(x); x <- x[o]; w <- w[o]
  cw <- (cumsum(w) - 0.5 * w) / sum(w)
  approx(cw, x, xout = probs, rule = 2, ties = "ordered")$y
}

# ==============================================================================
# Load (read-only)
# ==============================================================================

message("Loading varpart + FDR-adjusted p-values...")
r_sh  <- rast(file_nc,  subds = "shared")
r_tot <- rast(file_nc,  subds = "total_r2")
r_pf  <- rast(file_adj, subds = "p_full_adj")
r_ps  <- rast(file_adj, subds = "p_soil_adj")
r_pd  <- rast(file_adj, subds = "p_demand_adj")

if (is.na(crs(r_sh)) || crs(r_sh) == "") crs(r_sh) <- "EPSG:4326"

message("Computing geodesic cell area (terra::cellSize, EPSG:4326)...")
area_v <- as.vector(values(cellSize(r_sh[[1]], unit = "km", mask = FALSE)))

message("Aligning aridity index to the analysis grid...")
ai <- rast(file_ai)[[1]]
if (is.na(crs(ai))) crs(ai) <- "EPSG:4326"
if (xmax(ai) > 180 + 1e-6) ai <- rotate(ai)
ai_r   <- resample(ai, r_sh[[1]], method = "bilinear")
ai_v   <- as.vector(values(ai_r))
ai_cls <- cut(ai_v, breaks = ai_breaks, labels = ai_labels, right = FALSE)

# ==============================================================================
# Classify: 5 mutually exclusive categories, FDR-adjusted
# ==============================================================================

classify_month <- function(m) {
  sh <- set_na(values(r_sh[[m]])[, 1])
  tt <- set_na(values(r_tot[[m]])[, 1])
  pf <- set_na(values(r_pf[[m]])[, 1])
  ps <- set_na(values(r_ps[[m]])[, 1])
  pd <- set_na(values(r_pd[[m]])[, 1])

  valid    <- !is.na(tt)
  full_sig <- valid & !is.na(pf) & pf < ALPHA
  soil_sig <- !is.na(ps) & ps < ALPHA
  dem_sig  <- !is.na(pd) & pd < ALPHA

  code <- rep(NA_integer_, length(sh))
  code[valid & !full_sig]                       <- 1L  # no_signal
  code[full_sig &  soil_sig & !dem_sig]         <- 2L  # soil
  code[full_sig & !soil_sig &  dem_sig]         <- 3L  # demand
  code[full_sig & !soil_sig & !dem_sig]         <- 4L  # coupled_pure  <- TRUE non-separable
  code[full_sig &  soil_sig &  dem_sig]         <- 5L  # both_unique   <- SEPARABLE
  list(code = code, sh = sh, tot = tt)
}

message("Classifying 12 months...")
code_all <- integer(0); area_all <- numeric(0); sh_all <- numeric(0)
tot_all  <- numeric(0); ai_all   <- numeric(0); cls_all <- character(0)
mon_all  <- integer(0)

for (m in 1:12) {
  cm  <- classify_month(m)
  keep <- !is.na(cm$code)
  code_all <- c(code_all, cm$code[keep])
  area_all <- c(area_all, area_v[keep])
  sh_all   <- c(sh_all,   cm$sh[keep])
  tot_all  <- c(tot_all,  cm$tot[keep])
  ai_all   <- c(ai_all,   ai_v[keep])
  cls_all  <- c(cls_all,  as.character(ai_cls[keep]))
  mon_all  <- c(mon_all,  rep(m, sum(keep)))
  cat(sprintf("  month %2d: %d valid cells\n", m, sum(keep)))
}

N_area <- sum(area_all)

# ==============================================================================
# A. Split of the joint control area (categories 4 + 5)
# ==============================================================================

cat("\n================================================================\n")
cat("=== Split of the joint control area (categories 4 + 5)\n")
cat("===   FDR-adjusted (p_adj < 0.05), area-weighted (geodesic km2)\n")
cat("================================================================\n")

is_acop <- code_all %in% c(4L, 5L)          # joint control (categories 4 + 5)
a_pure  <- sum(area_all[code_all == 4L])
a_both  <- sum(area_all[code_all == 5L])
a_acop  <- a_pure + a_both

cat(sprintf("\n  Joint control (categories 4 + 5) = %.1f x 10^6 km2 (cell-months)\n",
            a_acop / 1e6))
cat(sprintf("    (4) coupled_pure  neither unique significant  (non-separable) : %6.2f%%\n",
            a_pure / a_acop * 100))
cat(sprintf("    (5) both_unique   both unique significant     (separable)     : %6.2f%%\n",
            a_both / a_acop * 100))

# ==============================================================================
# B. 5-category area-weighted shares: whole domain and signal-only
# ==============================================================================

csv_rows <- list()
add_rows <- function(scope, domain, cats, area) {
  tot <- sum(area)
  for (k in seq_along(cats))
    csv_rows[[length(csv_rows) + 1L]] <<- data.frame(
      scope = scope, domain = domain, category = cats[k],
      area_km2 = area[k], pct_area = if (tot > 0) area[k] / tot * 100 else NA_real_,
      stringsAsFactors = FALSE)
}

area_by_cat <- function(sel) sapply(1:5, function(k) sum(area_all[sel & code_all == k]))

print_cats <- function(title, area, cats = CATS) {
  tot <- sum(area)
  cat(sprintf("\n  %s\n", title))
  cat(sprintf("  %-14s | %16s | %8s\n", "category", "area_km2", "% area"))
  cat(sprintf("  %s\n", strrep("-", 46)))
  for (k in seq_along(cats))
    cat(sprintf("  %-14s | %16.1f | %7.2f%%\n", cats[k], area[k], area[k] / tot * 100))
}

cat("\n----------------------------------------------------------------\n")
cat("B. Five-category split (area-weighted, all 12 months)\n")
cat("----------------------------------------------------------------\n")

a_glob <- area_by_cat(rep(TRUE, length(code_all)))
print_cats("WHOLE vegetated domain:", a_glob)
add_rows("global_12mo", "whole", CATS, a_glob)

print_cats("SIGNAL-ONLY (excl. no_signal):",
           a_glob[2:5], CATS[2:5])
add_rows("global_12mo", "signal_only", CATS[2:5], a_glob[2:5])

cat("\n  Merged joint bucket (categories 4 + 5 combined), signal-only shares:\n")
sig_tot <- sum(a_glob[2:5])
cat(sprintf("    soil   : %6.2f%%\n", a_glob[2] / sig_tot * 100))
cat(sprintf("    demand : %6.2f%%\n", a_glob[3] / sig_tot * 100))
cat(sprintf("    joint  : %6.2f%%  = coupled_pure %.2f%% + both_unique %.2f%%\n",
            (a_glob[4] + a_glob[5]) / sig_tot * 100,
            a_glob[4] / sig_tot * 100, a_glob[5] / sig_tot * 100))

# per displayed month
cat("\n----- per month (signal-only shares) -----\n")
for (mi in seq_along(months_index)) {
  sel <- mon_all == months_index[mi]
  a_m <- area_by_cat(sel)
  print_cats(sprintf("%s -- SIGNAL-ONLY:", months_names[mi]), a_m[2:5], CATS[2:5])
  add_rows(months_names[mi], "signal_only", CATS[2:5], a_m[2:5])
}

# ==============================================================================
# C. Explained variance and shared sign, per category
# ==============================================================================

cat("\n----------------------------------------------------------------\n")
cat("C. Full-model adjusted R2 per category (area-weighted quantiles)\n")
cat("   Splitting the joint bucket shows whether the joint-category R2 reflects the\n")
cat("   genuinely coupled cells or is lowered by merging them with separable cells.\n")
cat("----------------------------------------------------------------\n")
cat(sprintf("  %-14s | %7s %7s %7s | %10s | %10s\n",
            "category", "p25", "median", "p75", "% shared>0", "area_Mkm2"))
cat(sprintf("  %s\n", strrep("-", 68)))
for (k in 2:5) {
  sel <- code_all == k
  q   <- wquantile(tot_all[sel], area_all[sel], c(.25, .50, .75))
  psh <- 100 * sum(area_all[sel & sh_all > 0], na.rm = TRUE) / sum(area_all[sel])
  cat(sprintf("  %-14s | %7.3f %7.3f %7.3f | %9.1f%% | %10.2f\n",
              CATS[k], q[1], q[2], q[3], psh, sum(area_all[sel]) / 12 / 1e6))
  csv_rows[[length(csv_rows) + 1L]] <- data.frame(
    scope = "r2_and_shared", domain = "global_12mo", category = CATS[k],
    area_km2 = sum(area_all[sel]), pct_area = NA_real_, stringsAsFactors = FALSE)
}

# ==============================================================================
# D. Along the aridity gradient + the corrected soil -> coupled crossover
# ==============================================================================

cat("\n----------------------------------------------------------------\n")
cat("D. Composition per Zomer aridity class (% of the class analysed area)\n")
cat("----------------------------------------------------------------\n")

arid_rows <- list()
for (cl in ai_labels) {
  sel_cl <- cls_all == cl & !is.na(cls_all)
  if (!any(sel_cl)) next
  a_cl  <- sapply(1:5, function(k) sum(area_all[sel_cl & code_all == k]))
  tot   <- sum(a_cl)
  arid_rows[[length(arid_rows) + 1L]] <- data.frame(
    ai_class = cl,
    area_Mkm2 = round(tot / 12 / 1e6, 3),
    no_signal    = round(a_cl[1] / tot * 100, 1),
    soil         = round(a_cl[2] / tot * 100, 1),
    demand       = round(a_cl[3] / tot * 100, 1),
    coupled_pure = round(a_cl[4] / tot * 100, 1),
    both_unique  = round(a_cl[5] / tot * 100, 1),
    stringsAsFactors = FALSE)
}
arid_df <- do.call(rbind, arid_rows)
print(arid_df, row.names = FALSE)
write.csv(arid_df, out_arid, row.names = FALSE)

cat("\n--- Soil -> coupled_pure crossover along the AI gradient ---\n")
cat("    Aridity index at which coupled_pure (category 4) first overtakes soil\n")
cat("    control (category 2), computed in 0.02-wide AI bins.\n")

fb   <- seq(0, 1, by = 0.02)
mids <- head(fb, -1) + 0.01
bin  <- cut(ai_all, breaks = fb, right = FALSE, include.lowest = TRUE)

soil_b <- tapply(area_all[code_all == 2L], bin[code_all == 2L], sum)
pure_b <- tapply(area_all[code_all == 4L], bin[code_all == 4L], sum)
merg_b <- tapply(area_all[is_acop],        bin[is_acop],        sum)
soil_b[is.na(soil_b)] <- 0; pure_b[is.na(pure_b)] <- 0; merg_b[is.na(merg_b)] <- 0

first_cross <- function(a, b) {              # first bin where b overtakes a
  d  <- as.numeric(b) - as.numeric(a)
  ok <- (as.numeric(a) + as.numeric(b)) > 0
  i  <- which(d > 0 & c(FALSE, head(d, -1) <= 0) & ok)
  if (length(i)) mids[i[1]] else NA_real_
}

x_merged <- first_cross(soil_b, merg_b)      # soil -> merged joint (cats 4+5)
x_pure   <- first_cross(soil_b, pure_b)      # soil -> coupled_pure (cat 4)

cat(sprintf("\n  soil -> joint (cats 4+5)     crossover at AI ~ %s\n",
            ifelse(is.na(x_merged), "none", sprintf("%.2f", x_merged))))
cat(sprintf("  soil -> coupled_pure (cat 4) crossover at AI ~ %s\n",
            ifelse(is.na(x_pure), "none", sprintf("%.2f", x_pure))))

csv_df <- do.call(rbind, csv_rows)
write.csv(csv_df, out_main, row.names = FALSE)

# ==============================================================================
# File protection
# ==============================================================================

cat("\n=== File protection (inputs READ-ONLY) ===\n")
inputs <- c(file_nc, file_adj, file_ai)
mt_after <- file.info(inputs)$mtime
for (i in seq_along(inputs))
  cat(sprintf("  %-42s : %s\n", basename(inputs[i]),
              if (identical(mt_before[i], mt_after[i])) "UNCHANGED" else "*** CHANGED ***"))
cat("\n  Created:\n")
cat(sprintf("    %s\n    %s\n", out_main, out_arid))
cat("\n=== DONE ===\n")
