#!/usr/bin/env Rscript

# ==============================================================================
# Script: 03_variance_partitioning/02b_apply_fdr_pooled.R
#
# Project: Aridity and timescale bound the separability of vegetation water stress
#
# Purpose:
#   Benjamini-Hochberg FDR for a TWO-WINDOW comparison, pooling both windows
#   into a SINGLE set of tests per test type, so that both are thresholded with
#   the same critical value.
#
#   Why this exists (do not replace it with two runs of 02_apply_fdr.R):
#     BH sets its cutoff from the p-value distribution of the pool it is given.
#     Run per window, each window gets its OWN cutoff, and the window with more
#     small p-values also gets the more permissive one. On the 1982-2001 vs
#     2002-2022 split that turned a modest raw difference into a cliff:
#
#       p_soil, 1982-2001 : 244332 cells at the minimum p (0.001)
#                           BH needed 255001 to cross -> 0.00% significant
#       p_soil, 2002-2022 : 385710 cells at the minimum p
#                           BH cutoff comfortably cleared -> 4.14% significant
#
#     A 4.4% shortfall in one window became "0% vs 4.14%", which reads as a
#     temporal collapse of soil control and is nothing of the sort. Pooling both
#     windows (m doubles, but so does the mass at the minimum p) puts the
#     minimum-p group above the critical value in BOTH windows.
#
#   Note the resolution floor: with NPERM=999 the smallest attainable p is
#   0.001, and the BH critical value only reaches 0.001 around rank 255000 for
#   a single window. Significance therefore depends on how much mass sits at
#   the minimum p, not on how extreme the best cell is.
#
# Inputs (READ-ONLY; mtimes confirmed unchanged):
#   outputs/varpart_global/varpart_signif_global_2blocks_{window}.nc
#     (p_full, p_soil, p_demand, total_r2)  — one per window
#
# Outputs (refuses to overwrite; one per window):
#   outputs/varpart_global/fdr_adjusted_pvalues_pooled_{window}.nc
#     (p_full_adj, p_soil_adj, p_demand_adj)
#
# Run from repo root, AFTER 01_varpart.R has been run for BOTH windows:
#   Rscript 03_variance_partitioning/02b_apply_fdr_pooled.R
#
# This script does NOT read VWS_PERIOD: it reads both windows by construction.
#
# Dependencies: ncdf4
# ==============================================================================

suppressPackageStartupMessages(library(ncdf4))

# ------------------------------------------------------------------------------
# Configuration
# ------------------------------------------------------------------------------

WINDOWS <- c("1982_2001", "2002_2022")

FV    <- -9999.0
ALPHA <- 0.05
TESTS <- c("p_full", "p_soil", "p_demand")

varpart_dir <- file.path("outputs", "varpart_global")

file_src <- setNames(
  file.path(varpart_dir, sprintf("varpart_signif_global_2blocks_%s.nc", WINDOWS)),
  WINDOWS)
file_out <- setNames(
  file.path(varpart_dir, sprintf("fdr_adjusted_pvalues_pooled_%s.nc", WINDOWS)),
  WINDOWS)

for (w in WINDOWS) {
  if (!file.exists(file_src[[w]]))
    stop("Source NetCDF not found: ", file_src[[w]], call. = FALSE)
  if (file.exists(file_out[[w]]))
    stop("REFUSING TO OVERWRITE — output already exists: ", file_out[[w]],
         "\n  Delete/rename it manually if you really want to regenerate.",
         call. = FALSE)
}

mtime_before <- sapply(file_src, function(f) as.character(file.info(f)$mtime))

mask_fv <- function(a) { a[abs(a - FV) < 1] <- NA_real_; a }

# ==============================================================================
# Read both windows (READ-ONLY)
# ==============================================================================

cat("Pooled BH-FDR across windows: ", paste(WINDOWS, collapse = "  +  "), "\n\n")

raw <- list(); axes <- list(); dims <- list()

for (w in WINDOWS) {
  cat(sprintf("Opening READ-ONLY: %s\n", file_src[[w]]))
  nc <- nc_open(file_src[[w]], write = FALSE)
  axes[[w]] <- list(lon   = ncvar_get(nc, "lon"),
                    lat   = ncvar_get(nc, "lat"),
                    month = ncvar_get(nc, "month"))
  raw[[w]] <- lapply(setNames(TESTS, TESTS), function(v) mask_fv(ncvar_get(nc, v)))
  nc_close(nc)
  dims[[w]] <- dim(raw[[w]][["p_full"]])
  cat(sprintf("  dims %s   mtime %s\n",
              paste(dims[[w]], collapse = " × "), mtime_before[[w]]))
}

if (!identical(dims[[WINDOWS[1]]], dims[[WINDOWS[2]]]))
  stop("Windows have different grids; cannot pool.", call. = FALSE)
if (!isTRUE(all.equal(axes[[WINDOWS[1]]]$lat, axes[[WINDOWS[2]]]$lat)))
  stop("Windows have different latitude axes; cannot pool.", call. = FALSE)

# ==============================================================================
# Pooled BH per test type
# ==============================================================================

cat("\n================================================================\n")
cat("=== Benjamini-Hochberg FDR, ONE pool per test type across windows\n")
cat("================================================================\n\n")
cat(sprintf("  Threshold: %.2f\n\n", ALPHA))

adj <- setNames(vector("list", length(WINDOWS)), WINDOWS)
for (w in WINDOWS) adj[[w]] <- list()

for (v in TESTS) {

  v1 <- as.vector(raw[[WINDOWS[1]]][[v]])
  v2 <- as.vector(raw[[WINDOWS[2]]][[v]])
  n1 <- length(v1)

  pooled <- c(v1, v2)
  ok     <- !is.na(pooled)
  a      <- pooled
  a[ok]  <- p.adjust(pooled[ok], method = "BH")

  # Split back, preserving each window's array shape.
  a1 <- a[seq_len(n1)]
  a2 <- a[(n1 + 1L):length(a)]
  dim(a1) <- dims[[WINDOWS[1]]]
  dim(a2) <- dims[[WINDOWS[2]]]
  adj[[WINDOWS[1]]][[v]] <- a1
  adj[[WINDOWS[2]]][[v]] <- a2

  # Report, per window, raw vs pooled-FDR, plus what a per-window FDR would give.
  cat(sprintf("  %s   (pooled m = %d)\n", v, sum(ok)))
  for (w in WINDOWS) {
    r  <- as.vector(raw[[w]][[v]]); r <- r[!is.na(r)]
    ap <- as.vector(adj[[w]][[v]]); ap <- ap[!is.na(ap)]
    solo <- p.adjust(r, method = "BH")
    cat(sprintf("    %-10s n=%9d | RAW %6.2f%% | FDR per-window %6.2f%% | FDR pooled %6.2f%%\n",
                w, length(r), 100 * mean(r < ALPHA),
                100 * mean(solo < ALPHA), 100 * mean(ap < ALPHA)))
  }
  cat("\n")
}

# ==============================================================================
# Write one NetCDF per window
# ==============================================================================

for (w in WINDOWS) {

  cat("Writing: ", file_out[[w]], "\n")

  dim_lon   <- ncdim_def("lon",   "degrees_east",  axes[[w]]$lon)
  dim_lat   <- ncdim_def("lat",   "degrees_north", axes[[w]]$lat)
  dim_month <- ncdim_def("month", "1",             axes[[w]]$month)

  vdef <- function(nm, ln) ncvar_def(nm, "1", list(dim_lon, dim_lat, dim_month),
                                     FV, longname = ln, prec = "float")
  v_full <- vdef("p_full_adj",   "BH-FDR adjusted p-value, full model (pooled across windows)")
  v_soil <- vdef("p_soil_adj",   "BH-FDR adjusted p-value, unique SUMINISTRO (pooled across windows)")
  v_dem  <- vdef("p_demand_adj", "BH-FDR adjusted p-value, unique DEMANDA (pooled across windows)")

  nc_o <- nc_create(file_out[[w]], list(v_full, v_soil, v_dem), force_v4 = TRUE)
  put_fv <- function(nc, nm, a) { a[is.na(a)] <- FV; ncvar_put(nc, nm, a) }
  put_fv(nc_o, "p_full_adj",   adj[[w]][["p_full"]])
  put_fv(nc_o, "p_soil_adj",   adj[[w]][["p_soil"]])
  put_fv(nc_o, "p_demand_adj", adj[[w]][["p_demand"]])

  ncatt_put(nc_o, 0, "title",  "BH-FDR adjusted p-values, pooled across comparison windows")
  ncatt_put(nc_o, 0, "method", paste0("p.adjust(method='BH'); ONE pool per test type spanning ",
                                      paste(WINDOWS, collapse = " + "),
                                      "; identical critical value in both windows"))
  ncatt_put(nc_o, 0, "window", w)
  ncatt_put(nc_o, 0, "pooled_with", paste(setdiff(WINDOWS, w), collapse = ", "))
  ncatt_put(nc_o, 0, "source", basename(file_src[[w]]))
  ncatt_put(nc_o, 0, "note",
            paste("Adjusted p-values only. Source NetCDF NOT modified.",
                  "NOT comparable with the full-record fdr_adjusted_pvalues.nc,",
                  "which is thresholded over a different pool."))
  ncatt_put(nc_o, 0, "date_created", as.character(Sys.time()))
  ncatt_put(nc_o, 0, "created_by", "03_variance_partitioning/02b_apply_fdr_pooled.R")
  nc_close(nc_o)
  cat("  written.\n")
}

# ==============================================================================
# File protection
# ==============================================================================

mtime_after <- sapply(file_src, function(f) as.character(file.info(f)$mtime))
cat("\n=== File protection confirmation ===\n")
for (w in WINDOWS)
  cat(sprintf("  %-46s : %s\n", basename(file_src[[w]]),
              if (identical(mtime_before[[w]], mtime_after[[w]]))
                "UNCHANGED (read-only respected)" else "*** CHANGED ***"))

cat("\n  Files CREATED by this script:\n")
for (w in WINDOWS) cat("    ", file_out[[w]], "\n")
cat("\n=== DONE ===\n")
cat(">>> Next: Rscript 04_analysis/10_window_comparison.R\n")
