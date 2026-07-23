#!/usr/bin/env Rscript

# ==============================================================================
# Script: 04_analysis/01_rhostar_areas.R
#
# Project: Aridity and timescale bound the separability of vegetation water stress
# Purpose:
#
#   Quantify the significant vegetation-hydroclimate correlation area: the share
#   of the vegetated domain where rho* is significant under BH-FDR (p_adj < 0.05),
#   per indicator and calendar month, plus the union across indicators. Reports
#   the area numbers; produces no figure.
#
#   FDR granularity: BH applied per indicator over all vegetated cell x month
#   p-values pooled (matching "all cell-month tests, by test type" in Methods).
#
# Inputs (READ-ONLY):
#   outputs/intermediate/absmax_spearman/abs_max_correlation_kndvi_{Ep,Et,ED,SMrz,SMs}.nc
#     (variable abs_max_p_value)
#   outputs/intermediate/vegetation_mask_c1.tif
#
# Output:
#   outputs/tables/rhostar_fdr_area.csv
#
# Run from repo root:
#   Rscript 04_analysis/01_rhostar_areas.R
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
out_csv <- file.path(tables_dir, "rhostar_fdr_area.csv")
if (file.exists(out_csv)) stop("REFUSING TO OVERWRITE: ", out_csv, call. = FALSE)

ALPHA <- 0.05
indicators <- c(Ep = "AED", Et = "Et", ED = "ED", SMrz = "SMrz", SMs = "SMs")
months_names <- month.name

ncfile <- function(v) file.path(dir_absmax, sprintf("abs_max_correlation_kndvi_%s.nc", v))
for (v in names(indicators))
  if (!file.exists(ncfile(v))) stop("Missing: ", ncfile(v), call. = FALSE)
if (!file.exists(file_mask)) stop("Missing veg mask: ", file_mask, call. = FALSE)
mt_before <- file.info(c(sapply(names(indicators), ncfile), file_mask))$mtime

# ------------------------------------------------------------------------------
# Grid, area, vegetated mask (aligned to the absmax grid via terra)
# ------------------------------------------------------------------------------

message("Setting up grid, geodesic area and vegetated mask...")
r_ref <- rast(ncfile("Ep"), subds = "abs_max_p_value")[[1]]
if (is.na(crs(r_ref)) || crs(r_ref) == "") crs(r_ref) <- "EPSG:4326"

area_v <- as.vector(values(cellSize(r_ref, unit = "km", mask = FALSE)))
msk <- rast(file_mask)
veg <- as.vector(values(resample(msk, r_ref, method = "near")))
veg <- !is.na(veg) & veg > 0
veg_area_Mkm2 <- sum(area_v[veg]) / 1e6
message(sprintf("  vegetated cells: %d   vegetated area: %.1f x 10^6 km^2",
                sum(veg), veg_area_Mkm2))

ncell <- length(area_v)
vidx  <- which(veg)

# union accumulator over indicators, per cell x month
any_fdr <- matrix(FALSE, nrow = ncell, ncol = 12)

rows <- list()

# ------------------------------------------------------------------------------
# Per indicator: BH-FDR significant area (pooled cell x month)
# ------------------------------------------------------------------------------

for (v in names(indicators)) {
  lab <- indicators[[v]]
  message(sprintf("\n=== %s (%s) ===", lab, v))
  r_p <- rast(ncfile(v), subds = "abs_max_p_value")            # 12 layers
  P   <- matrix(NA_real_, nrow = ncell, ncol = 12)
  for (m in 1:12) P[, m] <- as.vector(values(r_p[[m]]))
  rm(r_p); gc(verbose = FALSE)

  # restrict to vegetated cells; pool all months for this indicator's BH
  pv    <- as.vector(P[vidx, ])                                 # veg cells x 12 months
  ok    <- is.finite(pv)
  padj  <- rep(NA_real_, length(pv))
  padj[ok] <- p.adjust(pv[ok], method = "BH")

  fdr_flag <- is.finite(padj) & padj < ALPHA
  fdr_mat  <- matrix(fdr_flag, nrow = length(vidx), ncol = 12)

  for (m in 1:12) {
    a_fdr <- sum(area_v[vidx][fdr_mat[, m]])
    rows[[length(rows) + 1L]] <- data.frame(
      indicator = lab, month = months_names[m],
      area_Mkm2 = round(a_fdr / 1e6, 2),
      pct       = round(100 * a_fdr / sum(area_v[veg]), 1),
      stringsAsFactors = FALSE)
    any_fdr[vidx[fdr_mat[, m]], m] <- TRUE
  }
  sub <- do.call(rbind, rows); sub <- sub[sub$indicator == lab, ]
  message(sprintf("  significant (FDR): %.1f-%.1f %% of vegetated area",
                  min(sub$pct), max(sub$pct)))
  rm(P, pv, padj, fdr_flag, fdr_mat); gc(verbose = FALSE)
}

# ------------------------------------------------------------------------------
# Union across indicators: "significant with at least one indicator"
# ------------------------------------------------------------------------------

message("\n=== UNION (significant with >=1 indicator) ===")
for (m in 1:12) {
  a_fdr <- sum(area_v[any_fdr[, m]])
  rows[[length(rows) + 1L]] <- data.frame(
    indicator = "ANY (>=1)", month = months_names[m],
    area_Mkm2 = round(a_fdr / 1e6, 2),
    pct       = round(100 * a_fdr / sum(area_v[veg]), 1),
    stringsAsFactors = FALSE)
}

res <- do.call(rbind, rows)
write.csv(res, out_csv, row.names = FALSE)

# ------------------------------------------------------------------------------
# Headline console report
# ------------------------------------------------------------------------------

u <- res[res$indicator == "ANY (>=1)", ]
cat("\n================================================================\n")
cat("=== rho* significant area (BH-FDR)\n")
cat("================================================================\n")
cat(sprintf("  Vegetated domain: %.1f x 10^6 km^2\n\n", veg_area_Mkm2))
cat(sprintf("  UNION (>=1 indicator): %.1f-%.1f%% of veg area (%.0f-%.0f x 10^6 km^2)\n",
            min(u$pct), max(u$pct), min(u$area_Mkm2), max(u$area_Mkm2)))

cat("\n  Per indicator (range across months, %% of veg area):\n")
for (lab in indicators) {
  s <- res[res$indicator == lab, ]
  cat(sprintf("    %-6s | %5.1f - %-6.1f\n", lab, min(s$pct), max(s$pct)))
}

cat(sprintf("\nCSV written: %s\n", out_csv))
mt_after <- file.info(c(sapply(names(indicators), ncfile), file_mask))$mtime
cat("\n=== File protection ===\n")
cat(sprintf("  inputs unchanged: %s\n", if (identical(mt_before, mt_after)) "YES" else "*** NO ***"))
cat("=== DONE ===\n")
