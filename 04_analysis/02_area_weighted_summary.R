#!/usr/bin/env Rscript

# ==============================================================================
# Script: 04_analysis/02_area_weighted_summary.R
#
# Project: Aridity and timescale bound the separability of vegetation water stress
#
# Purpose:
#   FINAL paper numbers: AREA-WEIGHTED breakdown of the 4 final categories,
#   using BH-FDR significance (p_adj < 0.05). Geodesic cell area from
#   terra::cellSize() (EPSG:4326). Both area-weighted and pixel-count figures
#   are reported side by side.
#
#   4 final categories (FDR):
#     sin señal : p_full_adj not significant
#     suelo     : p_full_adj sig AND only p_soil_adj sig
#     demanda   : p_full_adj sig AND only p_demand_adj sig
#     acoplado  : p_full_adj sig AND (neither unique sig OR both unique sig)
#
#   Reports the 4-category split (area-weighted) three ways:
#     (a) global — all 12 months aggregated
#     (b) per displayed month — Jan / Apr / Jul / Oct (seasonality)
#     (c) over the WHOLE domain AND over SIGNAL-ONLY cells (excl. sin señal)
#   Plus, within "acoplado": % redundant (shared>0) vs compensatory (shared<0),
#   area-weighted, global and per month.
#
#   *** FILE PROTECTION ***
#   - Both NetCDFs are opened READ-ONLY (terra reads do not write). They are
#     NEVER modified.
#   - The output CSV uses a NEW name; file.exists() is checked and the script
#     STOPS rather than overwrite.
#   - Source mtimes are recorded before/after and confirmed unchanged.
#
# Output:
#   outputs/varpart_global/area_weighted_summary.csv
#
# Run from repo root (AFTER 03_variance_partitioning/02_apply_fdr.R):
#   Rscript 04_analysis/02_area_weighted_summary.R
#
# Dependencies: terra
# ==============================================================================

suppressPackageStartupMessages(library(terra))

file_nc  <- file.path("outputs", "varpart_global", "varpart_signif_global_2blocks.nc")
file_adj <- file.path("outputs", "varpart_global", "fdr_adjusted_pvalues.nc")
file_csv <- file.path("outputs", "varpart_global", "area_weighted_summary.csv")

FV    <- -9999
ALPHA <- 0.05

months_index <- c(1, 4, 7, 10)
months_names <- c("January", "April", "July", "October")

CATS <- c("sin señal", "suelo", "demanda", "acoplado")

# ── Guards ────────────────────────────────────────────────────────────────────
if (!file.exists(file_nc))  stop("Source NetCDF not found: ", file_nc)
if (!file.exists(file_adj)) stop("Adjusted-p NetCDF not found (run 03_variance_partitioning/02_apply_fdr.R first): ", file_adj)
if (file.exists(file_csv))
  stop("REFUSING TO OVERWRITE — output already exists: ", file_csv,
       "\n  Delete/rename it manually if you really want to regenerate.")

mt_nc_before  <- file.info(file_nc)$mtime
mt_adj_before <- file.info(file_adj)$mtime
cat(sprintf("Inputs (READ-ONLY):\n  %s  (mtime %s)\n  %s  (mtime %s)\n",
            file_nc, mt_nc_before, file_adj, mt_adj_before))

# ==============================================================================
# Load rasters (read-only) and compute geodesic cell area
# ==============================================================================

r_sh  <- rast(file_nc,  subds = "shared")
r_tot <- rast(file_nc,  subds = "total_r2")
r_pf  <- rast(file_adj, subds = "p_full_adj")
r_ps  <- rast(file_adj, subds = "p_soil_adj")
r_pd  <- rast(file_adj, subds = "p_demand_adj")

set_na <- function(v) { v[abs(v - FV) < 1] <- NA; v }

# Geodesic area per cell (km^2). Depends on latitude only for a geographic CRS.
cat("\nComputing geodesic cell area (terra::cellSize, EPSG:4326)...\n")
area_r   <- cellSize(r_sh[[1]], unit = "km", mask = FALSE)
area_vec <- values(area_r)[, 1]
ncol_r   <- ncol(r_sh); nrow_r <- nrow(r_sh)
lat_row  <- yFromRow(r_sh, seq_len(nrow_r))

# ── Sanity: area should fall off ~cos(lat) toward the poles ───────────────────
cat("\n=== Area sanity check: area(lat) vs cos(lat) ===\n")
row_area <- area_vec[((seq_len(nrow_r) - 1L) * ncol_r) + 1L]   # col 1 of each row
i_eq     <- which.min(abs(lat_row))            # row nearest equator
a_eq     <- row_area[i_eq]
cat(sprintf("  Equator-most row: lat=%.3f  area=%.3f km^2 (reference)\n",
            lat_row[i_eq], a_eq))
cat(sprintf("  %8s  %12s  %12s  %12s  %10s\n",
            "lat", "area_km2", "area/eq", "cos(lat)", "ratio/cos"))
for (tgt in c(0, 20, 40, 60, 75, 83)) {
  ir <- which.min(abs(abs(lat_row) - tgt))
  ar <- row_area[ir]
  cl <- cos(lat_row[ir] * pi / 180)
  cat(sprintf("  %8.2f  %12.4f  %12.4f  %12.4f  %10.4f\n",
              lat_row[ir], ar, ar / a_eq, cl, (ar / a_eq) / cl))
}
cat("  (area/eq should track cos(lat); ratio/cos ~ 1 confirms geodesic scaling)\n")

# ==============================================================================
# Classify each cell×month and accumulate area + pixel counts
# ==============================================================================
#
# code: 1 sin señal | 2 suelo | 3 demanda | 4 acoplado

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
  code[valid & !full_sig]                          <- 1L  # sin señal
  code[full_sig &  soil_sig & !dem_sig]            <- 2L  # suelo
  code[full_sig & !soil_sig &  dem_sig]            <- 3L  # demanda
  code[full_sig & (( !soil_sig & !dem_sig) |
                   (  soil_sig &  dem_sig))]        <- 4L  # acoplado
  list(code = code, valid = valid, sh = sh)
}

# Accumulators: matrix [category, metric]
acc_global_area <- numeric(4); acc_global_pix <- numeric(4)
acc_month_area  <- matrix(0, 4, length(months_index),
                          dimnames = list(CATS, months_names))
acc_month_pix   <- matrix(0, 4, length(months_index),
                          dimnames = list(CATS, months_names))

# Acoplado redundant/compensatory split: [c(redun, comp, zero), metric]
acc_coup_global <- list(area = c(redun = 0, comp = 0, zero = 0),
                        pix  = c(redun = 0, comp = 0, zero = 0))
acc_coup_month  <- list(
  area = matrix(0, 3, length(months_index),
                dimnames = list(c("redun","comp","zero"), months_names)),
  pix  = matrix(0, 3, length(months_index),
                dimnames = list(c("redun","comp","zero"), months_names)))

cat("\nClassifying 12 months and accumulating area/pixels...\n")
for (m in 1:12) {
  cm   <- classify_month(m)
  code <- cm$code; sh <- cm$sh
  for (k in 1:4) {
    sel <- which(code == k)
    acc_global_area[k] <- acc_global_area[k] + sum(area_vec[sel])
    acc_global_pix[k]  <- acc_global_pix[k]  + length(sel)
  }
  # acoplado split (global)
  sel4 <- which(code == 4L)
  if (length(sel4)) {
    a4 <- area_vec[sel4]; s4 <- sh[sel4]
    acc_coup_global$area["redun"] <- acc_coup_global$area["redun"] + sum(a4[s4 > 0], na.rm = TRUE)
    acc_coup_global$area["comp"]  <- acc_coup_global$area["comp"]  + sum(a4[s4 < 0], na.rm = TRUE)
    acc_coup_global$area["zero"]  <- acc_coup_global$area["zero"]  + sum(a4[s4 == 0], na.rm = TRUE)
    acc_coup_global$pix["redun"]  <- acc_coup_global$pix["redun"]  + sum(s4 > 0, na.rm = TRUE)
    acc_coup_global$pix["comp"]   <- acc_coup_global$pix["comp"]   + sum(s4 < 0, na.rm = TRUE)
    acc_coup_global$pix["zero"]   <- acc_coup_global$pix["zero"]   + sum(s4 == 0, na.rm = TRUE)
  }

  # per displayed month
  mi <- match(m, months_index)
  if (!is.na(mi)) {
    for (k in 1:4) {
      sel <- which(code == k)
      acc_month_area[k, mi] <- sum(area_vec[sel])
      acc_month_pix[k, mi]  <- length(sel)
    }
    if (length(sel4)) {
      a4 <- area_vec[sel4]; s4 <- sh[sel4]
      acc_coup_month$area["redun", mi] <- sum(a4[s4 > 0], na.rm = TRUE)
      acc_coup_month$area["comp",  mi] <- sum(a4[s4 < 0], na.rm = TRUE)
      acc_coup_month$area["zero",  mi] <- sum(a4[s4 == 0], na.rm = TRUE)
      acc_coup_month$pix["redun", mi]  <- sum(s4 > 0, na.rm = TRUE)
      acc_coup_month$pix["comp",  mi]  <- sum(s4 < 0, na.rm = TRUE)
      acc_coup_month$pix["zero",  mi]  <- sum(s4 == 0, na.rm = TRUE)
    }
  }
}

# ==============================================================================
# Report helpers
# ==============================================================================

csv_rows <- list()
add_rows <- function(table, scope, domain, cats, pix, area) {
  tot_pix <- sum(pix); tot_area <- sum(area)
  for (k in seq_along(cats)) {
    csv_rows[[length(csv_rows) + 1L]] <<- data.frame(
      table = table, scope = scope, domain = domain, category = cats[k],
      n_pixels  = as.numeric(pix[k]),
      pct_pixels= if (tot_pix  > 0) pix[k]  / tot_pix  * 100 else NA,
      area_km2  = area[k],
      pct_area  = if (tot_area > 0) area[k] / tot_area * 100 else NA,
      stringsAsFactors = FALSE)
  }
}

print_4cat <- function(title, pix, area, cats = CATS) {
  tot_pix <- sum(pix); tot_area <- sum(area)
  cat(sprintf("\n  %s\n", title))
  cat(sprintf("  %-12s | %12s %8s | %16s %8s\n",
              "Categoría", "n_pix", "% pix", "area_km2", "% area"))
  cat(sprintf("  %s\n", strrep("-", 64)))
  for (k in seq_along(cats))
    cat(sprintf("  %-12s | %12.0f %7.2f%% | %16.1f %7.2f%%\n",
                cats[k], pix[k], pix[k] / tot_pix * 100,
                area[k], area[k] / tot_area * 100))
  cat(sprintf("  %s\n", strrep("-", 64)))
  cat(sprintf("  %-12s | %12.0f          | %16.1f\n", "TOTAL", tot_pix, tot_area))
}

# ==============================================================================
# (a) GLOBAL — all 12 months aggregated
# ==============================================================================

cat("\n\n================================================================\n")
cat("=== FINAL NUMBERS (FDR, area-weighted)  —  PAPER\n")
cat("================================================================\n")

cat("\n----- (a) GLOBAL, all 12 months aggregated -----")
print_4cat("WHOLE domain:", acc_global_pix, acc_global_area)
add_rows("4cat", "global_12mo", "whole", CATS, acc_global_pix, acc_global_area)

sig_idx <- 2:4
print_4cat("SIGNAL-ONLY (excl. sin señal):",
           acc_global_pix[sig_idx], acc_global_area[sig_idx], CATS[sig_idx])
add_rows("4cat", "global_12mo", "signal_only", CATS[sig_idx],
         acc_global_pix[sig_idx], acc_global_area[sig_idx])

# ==============================================================================
# (b) PER MONTH — Jan/Apr/Jul/Oct
# ==============================================================================

cat("\n\n----- (b) PER MONTH (seasonality) -----")
for (mi in seq_along(months_index)) {
  print_4cat(sprintf("%s — WHOLE domain:", months_names[mi]),
             acc_month_pix[, mi], acc_month_area[, mi])
  add_rows("4cat", months_names[mi], "whole", CATS,
           acc_month_pix[, mi], acc_month_area[, mi])
  print_4cat(sprintf("%s — SIGNAL-ONLY:", months_names[mi]),
             acc_month_pix[sig_idx, mi], acc_month_area[sig_idx, mi], CATS[sig_idx])
  add_rows("4cat", months_names[mi], "signal_only", CATS[sig_idx],
           acc_month_pix[sig_idx, mi], acc_month_area[sig_idx, mi])
}

# ==============================================================================
# (c) WITHIN ACOPLADO — redundant vs compensatory
# ==============================================================================

cat("\n\n----- (c) WITHIN 'acoplado': redundant (sh>0) vs compensatory (sh<0) -----")

print_coup <- function(title, pix3, area3) {
  # report on redun+comp base (zero is negligible; shown if present)
  labs <- c(redun = "redundante (sh>0)", comp = "compensatorio (sh<0)", zero = "shared==0")
  use  <- if (pix3["zero"] > 0) c("redun","comp","zero") else c("redun","comp")
  tot_pix <- sum(pix3[use]); tot_area <- sum(area3[use])
  cat(sprintf("\n  %s\n", title))
  cat(sprintf("  %-22s | %12s %8s | %16s %8s\n",
              "Sub-categoría", "n_pix", "% pix", "area_km2", "% area"))
  cat(sprintf("  %s\n", strrep("-", 74)))
  for (u in use)
    cat(sprintf("  %-22s | %12.0f %7.2f%% | %16.1f %7.2f%%\n",
                labs[u], pix3[u], pix3[u] / tot_pix * 100,
                area3[u], area3[u] / tot_area * 100))
}

print_coup("GLOBAL (all 12 months):", acc_coup_global$pix, acc_coup_global$area)
add_rows("acoplado_split", "global_12mo", "acoplado",
         c("redundante (sh>0)", "compensatorio (sh<0)", "shared==0"),
         acc_coup_global$pix, acc_coup_global$area)

for (mi in seq_along(months_index)) {
  print_coup(sprintf("%s:", months_names[mi]),
             acc_coup_month$pix[, mi], acc_coup_month$area[, mi])
  add_rows("acoplado_split", months_names[mi], "acoplado",
           c("redundante (sh>0)", "compensatorio (sh<0)", "shared==0"),
           acc_coup_month$pix[, mi], acc_coup_month$area[, mi])
}

# ==============================================================================
# Write CSV (guarded) + file-protection confirmation
# ==============================================================================

csv_df <- do.call(rbind, csv_rows)
write.csv(csv_df, file_csv, row.names = FALSE, fileEncoding = "UTF-8")
cat(sprintf("\n\nCSV written: %s  (%d rows)\n", file_csv, nrow(csv_df)))

cat("\n=== File protection confirmation ===\n")
mt_nc_after  <- file.info(file_nc)$mtime
mt_adj_after <- file.info(file_adj)$mtime
cat(sprintf("  %s : %s\n", basename(file_nc),
            if (identical(mt_nc_before,  mt_nc_after))  "UNCHANGED ✓" else "*** CHANGED ***"))
cat(sprintf("  %s : %s\n", basename(file_adj),
            if (identical(mt_adj_before, mt_adj_after)) "UNCHANGED ✓" else "*** CHANGED ***"))
cat("\n  Files CREATED by this script:\n")
cat(sprintf("    %s\n", file_csv))

cat("\n=== DONE — these are the paper numbers (FDR, area-weighted) ===\n")
