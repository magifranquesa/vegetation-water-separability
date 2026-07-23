#!/usr/bin/env Rscript

# ==============================================================================
# Script: 04_analysis/08_vegetated_area_proportions_fdr.R
#
# Project: Aridity and timescale bound the separability of vegetation water stress
# Purpose:
#
#   Compute the proportions of the vegetated land surface in the five rho*
#   significance classes, using Benjamini-Hochberg FDR significance
#   (p_adj < 0.05), area-weighted, per indicator x calendar month. This is the
#   table that Supplementary Fig. S2 draws; the figure script only reads the CSV,
#   so the figure can be re-rendered without recomputing anything.
#
#   Classes (sum to 100% of the vegetated domain, area-weighted):
#     Significant positive     rho* > 0 and p_adj <  0.05
#     Non-significant positive rho* > 0 and p_adj >= 0.05
#     Non-significant negative rho* < 0 and p_adj >= 0.05
#     Significant negative     rho* < 0 and p_adj <  0.05
#     Not observed             vegetated cell with no valid rho*
#
#   FDR is applied per indicator over all vegetated cell x month p-values pooled,
#   identical to 02_correlations/04_fdr_significance_mask.R (so the result matches
#   the fdr_padj_*.tif masks). It is recomputed here from abs_max_p_value, so this
#   script depends only on the archived correlation NetCDFs and the vegetation
#   mask -- not on the intermediate FDR rasters.
#
# Inputs (READ-ONLY):
#   outputs/intermediate/absmax_spearman/abs_max_correlation_kndvi_{Ep,Et,ED,SMrz,SMs}.nc
#     (variables abs_max_correlation and abs_max_p_value)
#   outputs/intermediate/vegetation_mask_c1.tif
#
# Output:
#   outputs/tables/proportions_vegetated_fdr.csv
#   Columns: variable, month, month_id, sig_pos, sig_neg, nosig_pos, nosig_neg,
#            not_observed, domain_area_mill_km2
#     (the five class columns are fractions of the vegetated domain, per month)
#
# Run from repo root (feeds 05_figures/figS2.R):
#   Rscript 04_analysis/08_vegetated_area_proportions_fdr.R
#
# Dependencies: terra
# ==============================================================================

suppressPackageStartupMessages(library(terra))

config_file <- file.path("R", "config.R")
if (file.exists(config_file)) source(config_file)

dir_absmax <- if (exists("paths") && !is.null(paths$absmax_spearman)) paths$absmax_spearman else
  file.path("outputs", "intermediate", "absmax_spearman")
file_mask <- if (exists("paths") && !is.null(paths$veg_mask_c1)) paths$veg_mask_c1 else
  file.path("outputs", "intermediate", "vegetation_mask_c1.tif")
tables_dir <- if (exists("paths") && !is.null(paths$tables)) paths$tables else
  file.path("outputs", "tables")
dir.create(tables_dir, recursive = TRUE, showWarnings = FALSE)
out_csv <- file.path(tables_dir, "proportions_vegetated_fdr.csv")
if (file.exists(out_csv)) stop("REFUSING TO OVERWRITE: ", out_csv, call. = FALSE)

ALPHA <- 0.05
indicators   <- c(Ep = "AED", Et = "Et", ED = "ED", SMrz = "SMrz", SMs = "SMs")
months_names <- month.name
ncfile <- function(v) file.path(dir_absmax, sprintf("abs_max_correlation_kndvi_%s.nc", v))

for (v in names(indicators))
  if (!file.exists(ncfile(v))) stop("Missing: ", ncfile(v), call. = FALSE)
if (!file.exists(file_mask)) stop("Missing veg mask: ", file_mask, call. = FALSE)
mt_before <- file.info(c(sapply(names(indicators), ncfile), file_mask))$mtime

# ------------------------------------------------------------------------------
# Grid, geodesic cell area, vegetated mask (aligned to the absmax grid)
# ------------------------------------------------------------------------------

message("Setting up grid, geodesic area and vegetated mask...")
r_ref <- rast(ncfile("Ep"), subds = "abs_max_correlation")[[1]]
if (is.na(crs(r_ref)) || crs(r_ref) == "") crs(r_ref) <- "EPSG:4326"
area_v <- as.vector(values(cellSize(r_ref, unit = "km", mask = FALSE)))
veg <- as.vector(values(resample(rast(file_mask), r_ref, method = "near")))
veg <- !is.na(veg) & veg > 0
vidx <- which(veg)
veg_area <- sum(area_v[veg])
veg_area_Mkm2 <- veg_area / 1e6
ncell <- length(area_v)
message(sprintf("  vegetated cells: %d   vegetated area: %.1f x 10^6 km^2",
                length(vidx), veg_area_Mkm2))

# ------------------------------------------------------------------------------
# Per indicator: BH-FDR (pooled veg cell x month) then area-weighted class shares
# ------------------------------------------------------------------------------

rows <- list()
for (v in names(indicators)) {
  lab <- indicators[[v]]
  message(sprintf("=== %s (%s) ===", lab, v))
  r_cor <- rast(ncfile(v), subds = "abs_max_correlation")   # 12 layers
  r_p   <- rast(ncfile(v), subds = "abs_max_p_value")       # 12 layers

  # FDR-adjusted p-values, pooled over vegetated cell x month (as in script 04).
  P <- matrix(NA_real_, nrow = ncell, ncol = 12)
  for (m in 1:12) P[, m] <- as.vector(values(r_p[[m]]))
  pv <- as.vector(P[vidx, ]); ok <- is.finite(pv)
  padj_pooled <- rep(NA_real_, length(pv))
  padj_pooled[ok] <- p.adjust(pv[ok], method = "BH")
  Padj <- matrix(NA_real_, nrow = ncell, ncol = 12)
  Padj[vidx, ] <- padj_pooled

  for (m in 1:12) {
    cor  <- as.vector(values(r_cor[[m]]))
    padj <- Padj[, m]
    obs  <- veg & is.finite(cor)
    sig  <- obs & is.finite(padj) & padj < ALPHA
    a <- function(sel) sum(area_v[sel]) / veg_area
    rows[[length(rows) + 1L]] <- data.frame(
      variable = lab, month = months_names[m], month_id = m,
      sig_pos      = a(sig & cor > 0),
      sig_neg      = a(sig & cor < 0),
      nosig_pos    = a(obs & !sig & cor > 0),
      nosig_neg    = a(obs & !sig & cor < 0),
      not_observed = a(veg & !obs),
      domain_area_mill_km2 = round(veg_area_Mkm2, 3),
      stringsAsFactors = FALSE)
  }
  rm(r_cor, r_p, P, pv, padj_pooled, Padj); gc(verbose = FALSE)
}
df <- do.call(rbind, rows)
write.csv(df, out_csv, row.names = FALSE)
message("CSV written: ", out_csv)

# ------------------------------------------------------------------------------
# File protection
# ------------------------------------------------------------------------------

mt_after <- file.info(c(sapply(names(indicators), ncfile), file_mask))$mtime
cat("\n=== File protection (inputs READ-ONLY) ===\n")
cat(sprintf("  inputs unchanged: %s\n", if (identical(mt_before, mt_after)) "YES" else "*** NO ***"))
cat("=== DONE ===\n")
