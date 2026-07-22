#!/usr/bin/env Rscript

# ==============================================================================
# Script: 06_robustness_era5/08_apply_fdr.R
#
# Project: Aridity and timescale bound the separability of vegetation water stress
#
# Purpose:
#   STEP 1 — Apply Benjamini-Hochberg (BH) FDR correction to the raw p-values
#            of the significance NetCDF, SEPARATELY per test type, pooling ALL
#            months and cells within each type (global FDR per test type, NOT
#            per month):
#              all p_full   -> p_full_adj
#              all p_soil    -> p_soil_adj
#              all p_demand  -> p_demand_adj
#            Adjusted p-values are written to a NEW, separate NetCDF.
#            Reports, per type: % significant RAW (p<0.05) vs FDR (p_adj<0.05).
#
#   STEP 3 — 5-category significance table (no_signal / soil / demand /
#            both_significant / coupled_pure), RAW vs FDR side by side
#            (pixel counts, pre-area-weighting).
#
#   *** FILE PROTECTION ***
#   - The source NetCDF varpart_signif_global_2blocks_era5.nc is the primary
#     result. It is opened READ-ONLY (write = FALSE) and is
#     NEVER modified, appended to, or overwritten.
#   - Every output is written to a NEW path; before writing, file.exists() is
#     checked and the script STOPS if the file already exists.
#   - The source file's modification time is recorded before and after and
#     confirmed unchanged at the end.
#
# Output:
#   outputs/varpart_global/fdr_adjusted_pvalues_era5.nc
#     (p_full_adj, p_soil_adj, p_demand_adj + lon/lat/month coords)
#
# Run from repo root:
#   Rscript 06_robustness_era5/08_apply_fdr.R
#
# Dependencies: ncdf4
# ==============================================================================

suppressPackageStartupMessages(library(ncdf4))

file_src <- file.path("outputs", "varpart_global", "varpart_signif_global_2blocks_era5.nc")
file_out <- file.path("outputs", "varpart_global", "fdr_adjusted_pvalues_era5.nc")

FV    <- -9999.0
ALPHA <- 0.05

# ── Guards ────────────────────────────────────────────────────────────────────
if (!file.exists(file_src)) stop("Source NetCDF not found: ", file_src)
if (file.exists(file_out))
  stop("REFUSING TO OVERWRITE — output already exists: ", file_out,
       "\n  Delete/rename it manually if you really want to regenerate.")

# Record source mtime BEFORE touching anything (read-only sanity).
src_mtime_before <- file.info(file_src)$mtime
cat(sprintf("Source: %s\n  mtime (before): %s\n", file_src, src_mtime_before))

mask_fv <- function(a) { a[abs(a - FV) < 1] <- NA_real_; a }

# ==============================================================================
# Read raw p-values + axes (READ-ONLY)
# ==============================================================================

cat("\nOpening source READ-ONLY (write = FALSE)...\n")
nc <- nc_open(file_src, write = FALSE)

lon   <- ncvar_get(nc, "lon")
lat   <- ncvar_get(nc, "lat")
month <- ncvar_get(nc, "month")

p_full   <- mask_fv(ncvar_get(nc, "p_full"))
p_soil   <- mask_fv(ncvar_get(nc, "p_soil"))
p_demand <- mask_fv(ncvar_get(nc, "p_demand"))
total_r2 <- mask_fv(ncvar_get(nc, "total_r2"))   # to define the valid domain

nc_close(nc)
cat("  Read p_full, p_soil, p_demand, total_r2 + coords. Source closed.\n")

dims <- dim(p_full)
cat(sprintf("  Array dims: %s\n", paste(dims, collapse = " × ")))

# ==============================================================================
# STEP 1 — BH FDR per test type (pool all months & cells of that type)
# ==============================================================================

cat("\n================================================================\n")
cat("=== STEP 1 — Benjamini-Hochberg FDR per test type (global pool)\n")
cat("================================================================\n\n")

bh_adjust <- function(a) {
  ok      <- !is.na(a)
  adj     <- a                       # keep NA where raw is NA
  adj[ok] <- p.adjust(a[ok], method = "BH")
  list(adj = adj, n = sum(ok))
}

res_full <- bh_adjust(p_full)
res_soil <- bh_adjust(p_soil)
res_dem  <- bh_adjust(p_demand)

p_full_adj   <- res_full$adj
p_soil_adj   <- res_soil$adj
p_demand_adj <- res_dem$adj

report_drop <- function(name, raw, adj, n) {
  raw_ok <- raw[!is.na(raw)]
  adj_ok <- adj[!is.na(adj)]
  pc_raw <- mean(raw_ok < ALPHA) * 100
  pc_fdr <- mean(adj_ok < ALPHA) * 100
  cat(sprintf("  %-9s  n=%9d  | RAW p<0.05: %6.2f%%  -> FDR p_adj<0.05: %6.2f%%  (Δ %+6.2f pp)\n",
              name, n, pc_raw, pc_fdr, pc_fdr - pc_raw))
  c(raw = pc_raw, fdr = pc_fdr)
}

cat(sprintf("  Method: p.adjust(method=\"BH\"), pooled over all cells×months per type.\n"))
cat(sprintf("  Threshold: %.2f\n\n", ALPHA))
invisible(report_drop("p_full",   p_full,   p_full_adj,   res_full$n))
invisible(report_drop("p_soil",   p_soil,   p_soil_adj,   res_soil$n))
invisible(report_drop("p_demand", p_demand, p_demand_adj, res_dem$n))

# ==============================================================================
# Write NEW NetCDF with the three adjusted layers (guarded above)
# ==============================================================================

cat("\nWriting adjusted p-values to NEW file: ", file_out, "\n")

dim_lon   <- ncdim_def("lon",   "degrees_east",  lon)
dim_lat   <- ncdim_def("lat",   "degrees_north", lat)
dim_month <- ncdim_def("month", "1",             month)

vdef <- function(nm, ln) ncvar_def(nm, "1", list(dim_lon, dim_lat, dim_month),
                                    FV, longname = ln, prec = "float")
v_full <- vdef("p_full_adj",   "BH-FDR adjusted p-value, full model (global pool)")
v_soil <- vdef("p_soil_adj",   "BH-FDR adjusted p-value, unique SUMINISTRO (global pool)")
v_dem  <- vdef("p_demand_adj", "BH-FDR adjusted p-value, unique DEMANDA (global pool)")

nc_o <- nc_create(file_out, list(v_full, v_soil, v_dem), force_v4 = TRUE)
put_fv <- function(nc, nm, a) { a[is.na(a)] <- FV; ncvar_put(nc, nm, a) }
put_fv(nc_o, "p_full_adj",   p_full_adj)
put_fv(nc_o, "p_soil_adj",   p_soil_adj)
put_fv(nc_o, "p_demand_adj", p_demand_adj)

ncatt_put(nc_o, 0, "title", "BH-FDR adjusted p-values for varpart_signif_global_2blocks_era5 (ERA5-Land)")
ncatt_put(nc_o, 0, "method", "p.adjust(method='BH'); pooled over all cells×months separately per test type")
ncatt_put(nc_o, 0, "source", basename(file_src))
ncatt_put(nc_o, 0, "note", "Adjusted p-values only. Source NetCDF NOT modified.")
ncatt_put(nc_o, 0, "date_created", as.character(Sys.time()))
ncatt_put(nc_o, 0, "created_by", "06_robustness_era5/08_apply_fdr.R")
nc_close(nc_o)
cat("  NetCDF written.\n")

# ==============================================================================
# STEP 3 — 5-category table, RAW vs FDR (pixel counts, exploratory)
# ==============================================================================

cat("\n================================================================\n")
cat("=== STEP 3 — 5-category significance table: RAW vs FDR\n")
cat("===   PRELIMINAR — pixel counts, no area weighting, NOT final\n")
cat("================================================================\n\n")

valid <- !is.na(total_r2)
N     <- sum(valid)

classify5 <- function(pf, ps, pd) {
  full_sig <- valid & !is.na(pf) & pf < ALPHA
  soil_sig <- !is.na(ps) & ps < ALPHA
  dem_sig  <- !is.na(pd) & pd < ALPHA
  out <- array(NA_character_, dim(valid))
  out[valid & !full_sig]               <- "sin señal"
  out[full_sig &  soil_sig & !dem_sig] <- "solo suelo"
  out[full_sig & !soil_sig &  dem_sig] <- "solo demanda"
  out[full_sig &  soil_sig &  dem_sig] <- "ambos sig"
  out[full_sig & !soil_sig & !dem_sig] <- "acoplado puro"
  out
}

CAT_ORDER <- c("sin señal", "solo suelo", "solo demanda", "ambos sig", "acoplado puro")
cat_raw <- classify5(p_full,     p_soil,     p_demand)
cat_fdr <- classify5(p_full_adj, p_soil_adj, p_demand_adj)

t_raw <- table(factor(cat_raw[valid], levels = CAT_ORDER))
t_fdr <- table(factor(cat_fdr[valid], levels = CAT_ORDER))

cat(sprintf("  Valid cell×months: %d\n\n", N))
cat(sprintf("  %-16s | %12s %8s | %12s %8s | %10s\n",
            "Categoría", "n RAW", "% RAW", "n FDR", "% FDR", "Δ pp"))
cat(sprintf("  %s\n", strrep("-", 78)))
for (nm in CAT_ORDER) {
  nr <- as.integer(t_raw[nm]); nf <- as.integer(t_fdr[nm])
  pr <- nr / N * 100; pf <- nf / N * 100
  cat(sprintf("  %-16s | %12d %7.2f%% | %12d %7.2f%% | %+9.2f\n",
              nm, nr, pr, nf, pf, pf - pr))
}
cat(sprintf("  %s\n", strrep("-", 78)))
cat(sprintf("  %-16s | %12d %7.2f%% | %12d %7.2f%% |\n",
            "TOTAL", sum(t_raw), sum(t_raw)/N*100, sum(t_fdr), sum(t_fdr)/N*100))

# ==============================================================================
# File-protection confirmation
# ==============================================================================

cat("\n================================================================\n")
cat("=== File protection confirmation\n")
cat("================================================================\n\n")
src_mtime_after <- file.info(file_src)$mtime
cat(sprintf("  Source mtime before : %s\n", src_mtime_before))
cat(sprintf("  Source mtime after  : %s\n", src_mtime_after))
cat(sprintf("  Source UNCHANGED    : %s\n",
            if (identical(src_mtime_before, src_mtime_after)) "YES ✓ (read-only respected)"
            else "*** NO — investigate! ***"))

cat("\n  Files CREATED by this script:\n")
cat(sprintf("    %s\n", file_out))

cat("\n=== DONE (Step 1 + Step 3) ===\n")
cat(">>> Next: run the figure scripts in 05_figures/.\n")
