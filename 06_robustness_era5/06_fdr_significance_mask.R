#!/usr/bin/env Rscript

# ==============================================================================
# Script: 06_robustness_era5/06_fdr_significance_mask.R
#
# Project: Aridity and timescale bound the separability of vegetation water stress
#
# ERA5-Land counterpart of 02_correlations/04_fdr_significance_mask.R.
# Same method (BH-FDR pooled per indicator over all vegetated cell x month
# p-values); the inputs are rho* computed against ERA5-Land instead of GLEAM.
#
#   Compute and save the Benjamini-Hochberg FDR-adjusted p-value of rho* for each
#   indicator, applying the same significance control as the GLEAM branch to the
#   ERA5-Land inputs. BH is applied per indicator over all vegetated cell x month
#   p-values pooled. Non-vegetated cells are set to NA. Adjusted p-values are
#   written as a 12-layer GeoTIFF per indicator (layer = calendar month).
#
# Inputs (READ-ONLY):
#   outputs/intermediate/absmax_spearman_era5/abs_max_correlation_kndvi_{Ep,Et,ED,SMrz,SMs}_era5.nc
#     (variable abs_max_p_value)
#   outputs/intermediate/vegetation_mask_c1.tif   (shared with the GLEAM branch)
#
# Outputs:
#   outputs/intermediate/absmax_fdr_era5/fdr_padj_{Ep,Et,ED,SMrz,SMs}_era5.tif
#     (12 layers each)
#
# NOTE — the ERA5 NetCDFs differ from the GLEAM ones in ways that matter here:
#   - time dimension is named "month" (GLEAM: "time")   -> irrelevant, subds used
#   - values are double (GLEAM: float)                  -> irrelevant
#   - _FillValue is -9999 (GLEAM: NaN)                  -> MUST become NA
# If -9999 ever reached p.adjust() it would be treated as a valid p-value: it is
# finite and < 0.05, so every masked cell would count as "significant" and the
# figure would be silently wrong. The guard below turns anything outside [0, 1]
# into NA and reports how many cells that affected. It should normally be 0,
# because terra/GDAL honours _FillValue on read; it exists so that a change in
# that behaviour fails loudly instead of quietly.
#
# Run from repo root:
#   Rscript 06_robustness_era5/06_fdr_significance_mask.R
#
# Dependencies: terra
# ==============================================================================

suppressPackageStartupMessages(library(terra))

config_file <- file.path("R", "config.R")
if (file.exists(config_file)) source(config_file)

# ERA5 inputs live in their own directory; do NOT reuse paths$absmax_spearman,
# which points at the GLEAM rho* rasters.
dir_absmax <- file.path("outputs", "intermediate", "absmax_spearman_era5")

# The vegetation mask is grid-based and shared with the GLEAM branch.
file_mask <- if (exists("paths") && !is.null(paths$veg_mask_c1)) paths$veg_mask_c1 else
  file.path("outputs", "intermediate", "vegetation_mask_c1.tif")

dir_out <- file.path("outputs", "intermediate", "absmax_fdr_era5")
dir.create(dir_out, recursive = TRUE, showWarnings = FALSE)

indicators <- c("Ep", "Et", "ED", "SMrz", "SMs")
ncfile  <- function(v) file.path(dir_absmax, sprintf("abs_max_correlation_kndvi_%s_era5.nc", v))
outfile <- function(v) file.path(dir_out, sprintf("fdr_padj_%s_era5.tif", v))

for (v in indicators) if (!file.exists(ncfile(v))) stop("Missing: ", ncfile(v), call. = FALSE)
if (!file.exists(file_mask)) stop("Missing veg mask: ", file_mask, call. = FALSE)
for (v in indicators) if (file.exists(outfile(v)))
  stop("REFUSING TO OVERWRITE: ", outfile(v), call. = FALSE)
mt_before <- file.info(c(sapply(indicators, ncfile), file_mask))$mtime

# ------------------------------------------------------------------------------
# Vegetated mask on the absmax grid
# ------------------------------------------------------------------------------

message("Aligning vegetated mask to the ERA5 absmax grid...")
r_ref <- rast(ncfile("Ep"), subds = "abs_max_p_value")[[1]]
if (is.na(crs(r_ref)) || crs(r_ref) == "") crs(r_ref) <- "EPSG:4326"
veg <- as.vector(values(resample(rast(file_mask), r_ref, method = "near")))
veg <- !is.na(veg) & veg > 0
vidx <- which(veg)
message(sprintf("  vegetated cells: %d", length(vidx)))

# ------------------------------------------------------------------------------
# Per indicator: BH-FDR pooled over vegetated cell x month; write adjusted p
# ------------------------------------------------------------------------------

for (v in indicators) {
  message(sprintf("\n=== %s ===", v))
  r_p <- rast(ncfile(v), subds = "abs_max_p_value")            # 12 layers
  ncell <- ncell(r_p)
  P <- matrix(NA_real_, nrow = ncell, ncol = 12)
  for (m in 1:12) P[, m] <- as.vector(values(r_p[[m]]))

  # Guard against unconverted fill values (see NOTE in the header).
  n_bad <- sum(is.finite(P) & (P < 0 | P > 1))
  if (n_bad > 0) {
    message(sprintf("  WARNING: %d p-values outside [0, 1] (fill values not", n_bad))
    message("           honoured on read?) -> set to NA before BH adjustment.")
    P[is.finite(P) & (P < 0 | P > 1)] <- NA_real_
  }

  # pool vegetated cell x month, BH-adjust, map back
  pv <- as.vector(P[vidx, ])
  ok <- is.finite(pv)
  padj <- rep(NA_real_, length(pv))
  padj[ok] <- p.adjust(pv[ok], method = "BH")
  Padj <- matrix(NA_real_, nrow = ncell, ncol = 12)
  Padj[vidx, ] <- padj

  r_out <- rast(r_ref, nlyrs = 12)
  values(r_out) <- Padj
  names(r_out) <- month.abb
  writeRaster(r_out, outfile(v), overwrite = FALSE,
              gdal = c("COMPRESS=DEFLATE", "TILED=YES"))

  n_fdr <- sum(padj[ok] < 0.05, na.rm = TRUE)
  message(sprintf("  cell-months significant (FDR): %d", n_fdr))
  message("  wrote: ", outfile(v))
  rm(r_p, P, pv, padj, Padj, r_out); gc(verbose = FALSE)
}

# ------------------------------------------------------------------------------
# File protection
# ------------------------------------------------------------------------------
mt_after <- file.info(c(sapply(indicators, ncfile), file_mask))$mtime
cat("\n=== File protection (inputs READ-ONLY) ===\n")
cat(sprintf("  inputs unchanged: %s\n", if (identical(mt_before, mt_after)) "YES" else "*** NO ***"))
cat("\n  Created (ERA5 FDR-adjusted p-value rasters, 12 layers each):\n")
for (v in indicators) cat(sprintf("    %s\n", outfile(v)))
cat("\n=== DONE ===\n")
