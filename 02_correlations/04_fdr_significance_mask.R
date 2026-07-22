#!/usr/bin/env Rscript

# ==============================================================================
# Script: 02_correlations/04_fdr_significance_mask.R
#
# Project: Aridity and timescale bound the separability of vegetation water stress
# Purpose:
#
#   Compute and SAVE the Benjamini-Hochberg FDR-adjusted p-value of rho* for each
#   indicator, so that Fig. 1, Supplementary Fig. S1 and the area summaries can be
#   redrawn with the SAME multiple-comparison control used in the variance
#   partitioning, instead of the nominal p < 0.05. Reusable output; nothing else
#   recomputes the adjustment.
#
#   BH is applied per indicator over all vegetated cell x month p-values pooled
#   (matching "all cell-month tests, by test type" in Methods). Non-vegetated
#   cells are set to NA. Adjusted p-values are written as a 12-layer GeoTIFF per
#   indicator (layer = calendar month).
#
#   Caveat (documented, not fixed here): FDR corrects MULTIPLICITY, not the
#   SELECTION of rho* as max |rho| over five timescales; the latter would require
#   a permutation null for the maximum.
#
# Inputs (READ-ONLY):
#   outputs/intermediate/absmax_spearman/abs_max_correlation_kndvi_{Ep,Et,ED,SMrz,SMs}.nc
#     (variable abs_max_p_value)
#   outputs/intermediate/vegetation_mask_c1.tif
#
# Outputs:
#   outputs/intermediate/absmax_fdr/fdr_padj_{Ep,Et,ED,SMrz,SMs}.tif   (12 layers)
#
# Run from repo root:
#   Rscript 02_correlations/04_fdr_significance_mask.R
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

dir_out <- file.path("outputs", "intermediate", "absmax_fdr")
dir.create(dir_out, recursive = TRUE, showWarnings = FALSE)

indicators <- c("Ep", "Et", "ED", "SMrz", "SMs")
ncfile  <- function(v) file.path(dir_absmax, sprintf("abs_max_correlation_kndvi_%s.nc", v))
outfile <- function(v) file.path(dir_out, sprintf("fdr_padj_%s.tif", v))

for (v in indicators) if (!file.exists(ncfile(v))) stop("Missing: ", ncfile(v), call. = FALSE)
if (!file.exists(file_mask)) stop("Missing veg mask: ", file_mask, call. = FALSE)
for (v in indicators) if (file.exists(outfile(v)))
  stop("REFUSING TO OVERWRITE: ", outfile(v), call. = FALSE)
mt_before <- file.info(c(sapply(indicators, ncfile), file_mask))$mtime

# ------------------------------------------------------------------------------
# Vegetated mask on the absmax grid
# ------------------------------------------------------------------------------

message("Aligning vegetated mask to the absmax grid...")
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

  n_nom <- sum(pv[ok] < 0.05); n_fdr <- sum(padj[ok] < 0.05, na.rm = TRUE)
  message(sprintf("  cell-months: sig nominal=%d  sig FDR=%d  (%.1f%% retained)",
                  n_nom, n_fdr, 100 * n_fdr / n_nom))
  message("  wrote: ", outfile(v))
  rm(r_p, P, pv, padj, Padj, r_out); gc(verbose = FALSE)
}

# ------------------------------------------------------------------------------
# File protection
# ------------------------------------------------------------------------------
mt_after <- file.info(c(sapply(indicators, ncfile), file_mask))$mtime
cat("\n=== File protection (inputs READ-ONLY) ===\n")
cat(sprintf("  inputs unchanged: %s\n", if (identical(mt_before, mt_after)) "YES" else "*** NO ***"))
cat("\n  Created (FDR-adjusted p-value rasters, 12 layers each):\n")
for (v in indicators) cat(sprintf("    %s\n", outfile(v)))
cat("\n=== DONE ===\n")
