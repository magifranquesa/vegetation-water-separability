#!/usr/bin/env Rscript

# ==============================================================================
# Script: 04_analysis/06_supply_demand_coupling.R
#
# Project: Aridity and timescale bound the separability of vegetation water stress
# Purpose:
#   Quantify the coupling between the soil-water supply block and the atmospheric
#   demand block as a function of accumulation timescale. For every grid cell,
#   calendar month and accumulation timescale k in {1, 3, 6, 9, 12} months:
#
#       rho_couple(k) = Spearman( SM_k , AED_k )     [both linearly detrended]
#
#   using the same series, detrending and rank correlation as the vegetation
#   analysis, so the values are directly comparable. Area-weighted means of
#   |rho_couple| are summarised by aridity class and timescale.
#
#   Also reports the variance inflation factor, VIF = 1 / (1 - rho^2), as a
#   measure of the collinearity between the two blocks at each timescale.
#
# Inputs (READ-ONLY):
#   outputs/intermediate/accumulated/SMrz_GLEAM_v4.2a_MO_1982-2022_scale_{k}.nc
#   outputs/intermediate/accumulated/SMs_GLEAM_v4.2a_MO_1982-2022_scale_{k}.nc
#   outputs/intermediate/accumulated/Ep_GLEAM_v4.2a_MO_1982-2022_scale_{k}.nc
#   outputs/intermediate/vegetation_mask_c1.tif
#   data/processed/aridity/ai_1982_2022_period.nc
#
# Outputs (refuses to overwrite):
#   outputs/intermediate/coupling/coupling_{SUPPLY}_Ep_scale_{k}.tif   (12 layers)
#   outputs/tables/supply_demand_coupling.csv
#
# Run from repo root:
#   Rscript 04_analysis/06_supply_demand_coupling.R
#
# Runtime note: this streams the accumulated NetCDFs in latitude blocks, so memory
# stays bounded. Expect tens of minutes, dominated by disk I/O. Reduce LAT_BLOCK if
# memory is tight.
#
# Dependencies: ncdf4, terra
# ==============================================================================

suppressPackageStartupMessages({
  library(ncdf4)
  library(terra)
})

config_file <- file.path("R", "config.R")
if (file.exists(config_file)) source(config_file)

dir_acc <- if (exists("paths") && !is.null(paths$accumulated)) paths$accumulated else
  file.path("outputs", "intermediate", "accumulated")
file_mask <- if (exists("paths") && !is.null(paths$veg_mask_c1)) paths$veg_mask_c1 else
  file.path("outputs", "intermediate", "vegetation_mask_c1.tif")
file_ai <- if (exists("paths") && !is.null(paths$aridity_index)) paths$aridity_index else
  file.path("data", "processed", "aridity", "ai_1982_2022_period.nc")
tables_dir <- if (exists("paths") && !is.null(paths$tables)) paths$tables else
  file.path("outputs", "tables")

dir_out <- file.path("outputs", "intermediate", "coupling")
dir.create(dir_out,    recursive = TRUE, showWarnings = FALSE)
dir.create(tables_dir, recursive = TRUE, showWarnings = FALSE)

out_csv <- file.path(tables_dir, "supply_demand_coupling.csv")
if (file.exists(out_csv))
  stop("REFUSING TO OVERWRITE: ", out_csv, "\n  Delete/rename it manually.", call. = FALSE)

# ------------------------------------------------------------------------------
# Configuration
# ------------------------------------------------------------------------------

SCALES   <- c(1, 3, 6, 9, 12)
DEMAND   <- "Ep"                    # GLEAM potential evaporation == AED
SUPPLIES <- c("SMrz", "SMs")        # root-zone (primary) and surface

LAT_BLOCK <- 30                     # rows of latitude read at a time
NYEAR     <- 41                     # 1982-2022
ALPHA     <- 0.05

ai_breaks <- c(-Inf, 0.03, 0.20, 0.50, 0.65, Inf)
ai_labels <- c("Hyper-arid", "Arid", "Semi-arid", "Dry sub-humid", "Humid")

acc_file <- function(v, k)
  file.path(dir_acc, sprintf("%s_GLEAM_v4.2a_MO_1982-2022_scale_%d.nc", v, k))

for (k in SCALES) for (v in c(DEMAND, SUPPLIES))
  if (!file.exists(acc_file(v, k)))
    stop("Missing accumulated file: ", acc_file(v, k), call. = FALSE)
for (f in c(file_mask, file_ai))
  if (!file.exists(f)) stop("Missing input: ", f, call. = FALSE)

# ------------------------------------------------------------------------------
# Vectorised detrend + Spearman across cells
# ------------------------------------------------------------------------------

# Linear detrend of each ROW of M (cells x nyear). Same convention as
# 01_compute_monthly_spearman_correlations.R: remove the least-squares trend, keep
# the interannual residual.
detrend_rows <- function(M) {
  n  <- ncol(M)
  tt <- seq_len(n) - (n + 1) / 2                 # centred time
  b  <- as.vector(M %*% tt) / sum(tt^2)          # slope per row
  M - rowMeans(M) - outer(b, tt)
}

# TRUE where a row is not constant (a constant series carries no rank information)
row_varies <- function(M) {
  rmin <- do.call(pmin, as.data.frame(M))
  rmax <- do.call(pmax, as.data.frame(M))
  is.finite(rmax - rmin) & (rmax - rmin) > 0
}

# Rank every ROW of an m x n matrix, vectorised.
#
# apply(M, 1, rank) would be prohibitive here (~10^8 rows across the full run), so
# we rank all rows in one radix sort: flatten column-major, sort by (row, value),
# and write ranks 1..n back into the sorted positions of each row.
#
# NOTE: this is ties.method = "first", not "average". On linearly detrended
# floating-point residuals exact ties are vanishingly rare, so the difference from
# stats::cor(method = "spearman") is negligible. Constant rows are excluded upstream
# by row_varies().
rank_rows <- function(M) {
  m <- nrow(M); n <- ncol(M)
  v   <- as.vector(M)                                  # column-major
  rid <- rep.int(seq_len(m), n)                        # row id of each flat element
  o   <- order(rid, v, method = "radix")               # group by row, ascending value
  r   <- numeric(m * n)
  r[o] <- rep.int(seq_len(n), m)                       # ranks 1..n within each row
  dim(r) <- c(m, n)
  r
}

# Spearman between two matched sets of rows: rank each row, then Pearson on ranks.
spearman_rows <- function(A, B) {
  RA <- rank_rows(A); RB <- rank_rows(B)
  ca <- RA - rowMeans(RA)
  cb <- RB - rowMeans(RB)
  out <- rowSums(ca * cb) / sqrt(rowSums(ca^2) * rowSums(cb^2))
  out[!is.finite(out)] <- NA_real_
  out
}

# ------------------------------------------------------------------------------
# Grid geometry (taken from the demand file at scale 1)
# ------------------------------------------------------------------------------

nc0  <- nc_open(acc_file(DEMAND, 1))
lon  <- nc0$dim[[grep("^(lon|longitude|x)$", names(nc0$dim), ignore.case = TRUE)[1]]]$vals
lat  <- nc0$dim[[grep("^(lat|latitude|y)$",  names(nc0$dim), ignore.case = TRUE)[1]]]$vals
ntim <- nc0$dim[[grep("^(time|t)$",          names(nc0$dim), ignore.case = TRUE)[1]]]$len
varname_demand <- names(nc0$var)[1]
nc_close(nc0)

nlon <- length(lon); nlat <- length(lat)
lat_desc <- lat[1] > lat[nlat]          # TRUE if lat runs north -> south

message(sprintf("Grid: %d lon x %d lat x %d time  (lat %s)",
                nlon, nlat, ntim, if (lat_desc) "descending" else "ascending"))
if (ntim < NYEAR * 12)
  stop("Unexpected time length: ", ntim, " (expected >= ", NYEAR * 12, ")", call. = FALSE)

template <- rast(nrows = nlat, ncols = nlon,
                 xmin = min(lon) - 0.05, xmax = max(lon) + 0.05,
                 ymin = min(lat) - 0.05, ymax = max(lat) + 0.05,
                 crs = "EPSG:4326")

# ------------------------------------------------------------------------------
# Main loop: one (supply, scale) pair at a time, streamed in latitude blocks
# ------------------------------------------------------------------------------

# The CDO-accumulated files (scale > 1) carry a 2-D `time_bnds(time, bnds)`
# variable that is listed BEFORE the data variable, so names(nc$var)[1] is not the
# field we want. Pick the first variable that actually has 3 dimensions.
var_first <- function(f) {
  nc <- nc_open(f)
  nd <- vapply(nc$var, function(v) v$ndims, integer(1))
  v3 <- names(nc$var)[nd == 3]
  nc_close(nc)
  if (!length(v3)) stop("No 3-D variable found in ", basename(f), call. = FALSE)
  v3[1]
}

# --- self-test: rank_rows/spearman_rows must reproduce stats::cor exactly --------
local({
  set.seed(1)
  A <- matrix(rnorm(20 * NYEAR), 20, NYEAR)
  B <- matrix(rnorm(20 * NYEAR), 20, NYEAR)
  ours <- spearman_rows(A, B)
  ref  <- vapply(seq_len(20), function(i)
    stats::cor(A[i, ], B[i, ], method = "spearman"), numeric(1))
  if (max(abs(ours - ref)) > 1e-10)
    stop("rank_rows/spearman_rows self-test FAILED (max diff ",
         max(abs(ours - ref)), ") -- refusing to run.", call. = FALSE)
  message("Self-test passed: vectorised Spearman matches stats::cor to 1e-10.")
})

# One pass per accumulation window. AED is read ONCE per block and correlated
# against every supply variable, instead of re-reading it for each.
for (k in SCALES) {

  out_tif <- setNames(
    file.path(dir_out, sprintf("coupling_%s_%s_scale_%d.tif", SUPPLIES, DEMAND, k)),
    SUPPLIES)
  if (all(file.exists(out_tif))) {
    message(sprintf("  scale %2d: all outputs exist, skipping", k)); next
  }

  message(sprintf("\n=== scale %2d months | %s vs {%s} ===",
                  k, DEMAND, paste(SUPPLIES, collapse = ", ")))

  nc_d  <- nc_open(acc_file(DEMAND, k)); v_dem <- var_first(acc_file(DEMAND, k))
  nc_s  <- lapply(SUPPLIES, function(s) nc_open(acc_file(s, k)))
  v_sup <- vapply(SUPPLIES, function(s) var_first(acc_file(s, k)), character(1))
  names(nc_s) <- SUPPLIES

  rho <- lapply(SUPPLIES, function(s) matrix(NA_real_, nrow = nlon * nlat, ncol = 12))
  names(rho) <- SUPPLIES

  y0 <- 1L
  while (y0 <= nlat) {
    nyb <- min(LAT_BLOCK, nlat - y0 + 1L)

    Ad <- ncvar_get(nc_d, v_dem, start = c(1, y0, 1), count = c(nlon, nyb, ntim))
    dim(Ad) <- c(nlon * nyb, ntim)

    for (SUP in SUPPLIES) {
      As <- ncvar_get(nc_s[[SUP]], v_sup[[SUP]],
                      start = c(1, y0, 1), count = c(nlon, nyb, ntim))
      dim(As) <- c(nlon * nyb, ntim)

      for (m in 1:12) {
        idx <- seq(m, by = 12, length.out = NYEAR)
        Xs  <- As[, idx, drop = FALSE]
        Xd  <- Ad[, idx, drop = FALSE]

        ok <- stats::complete.cases(Xs) & stats::complete.cases(Xd)
        if (!any(ok)) next
        ok[ok] <- row_varies(Xs[ok, , drop = FALSE]) & row_varies(Xd[ok, , drop = FALSE])
        if (!any(ok)) next

        r <- spearman_rows(detrend_rows(Xs[ok, , drop = FALSE]),
                           detrend_rows(Xd[ok, , drop = FALSE]))

        # block cell index is lon-major; map back to the global lon-major index
        cb <- which(ok)
        j  <- ((cb - 1L) %/% nlon) + y0        # global lat index
        i  <- ((cb - 1L) %%  nlon) + 1L        # global lon index
        rho[[SUP]][(j - 1L) * nlon + i, m] <- r
      }
      rm(As); gc(verbose = FALSE)
    }

    cat(sprintf("    lat rows %4d-%4d done\n", y0, y0 + nyb - 1L))
    y0 <- y0 + nyb
    rm(Ad); gc(verbose = FALSE)
  }

  nc_close(nc_d); invisible(lapply(nc_s, nc_close))

  # rho is lon-major (lon varies fastest, lat index 1 = first NetCDF row). terra
  # wants row-major from the NORTH-WEST corner, so reshape to [lat, lon], flip if
  # the NetCDF latitudes run south -> north, and flatten.
  for (SUP in SUPPLIES) {
    if (file.exists(out_tif[[SUP]])) next
    vals <- matrix(NA_real_, nrow = nlon * nlat, ncol = 12)
    for (m in 1:12) {
      mat <- matrix(rho[[SUP]][, m], nrow = nlat, ncol = nlon, byrow = TRUE)
      if (!lat_desc) mat <- mat[nlat:1, , drop = FALSE]
      vals[, m] <- as.vector(t(mat))
    }
    r_out <- rast(template, nlyrs = 12)
    values(r_out) <- vals
    names(r_out) <- month.abb
    writeRaster(r_out, out_tif[[SUP]], overwrite = FALSE,
                gdal = c("COMPRESS=DEFLATE", "TILED=YES"))
    message("  wrote: ", out_tif[[SUP]])
    rm(vals, r_out)
  }
  rm(rho); gc(verbose = FALSE)
}

# ==============================================================================
# Summarise: coupling by timescale and aridity class
# ==============================================================================

message("\n=== Summarising coupling by timescale and aridity class ===")

msk <- rast(file_mask)
ai  <- rast(file_ai)[[1]]
if (is.na(crs(ai))) crs(ai) <- "EPSG:4326"
if (xmax(ai) > 180 + 1e-6) ai <- rotate(ai)

ref    <- rast(file.path(dir_out, sprintf("coupling_%s_%s_scale_1.tif", SUPPLIES[1], DEMAND)))[[1]]
msk_v  <- as.vector(values(resample(msk, ref, method = "near")))
ai_v   <- as.vector(values(resample(ai,  ref, method = "bilinear")))
area_v <- as.vector(values(cellSize(ref, unit = "km", mask = FALSE)))
cls_v  <- cut(ai_v, breaks = ai_breaks, labels = ai_labels, right = FALSE)

veg <- !is.na(msk_v) & msk_v > 0 & !is.na(ai_v)

wmean <- function(x, w) { ok <- is.finite(x) & is.finite(w); sum(x[ok] * w[ok]) / sum(w[ok]) }

rows <- list()
for (SUP in SUPPLIES) for (k in SCALES) {
  r <- rast(file.path(dir_out, sprintf("coupling_%s_%s_scale_%d.tif", SUP, DEMAND, k)))
  for (m in 1:12) {
    v   <- as.vector(values(r[[m]]))
    sel <- veg & is.finite(v)

    for (cl in c("ALL", ai_labels)) {
      s <- if (cl == "ALL") sel else sel & !is.na(cls_v) & cls_v == cl
      if (!any(s)) next
      mabs  <- wmean(abs(v[s]), area_v[s])
      msign <- wmean(v[s],      area_v[s])
      rows[[length(rows) + 1L]] <- data.frame(
        supply      = SUP,
        demand      = DEMAND,
        scale_months = k,
        month       = m,
        ai_class    = cl,
        mean_abs_rho = round(mabs, 4),
        mean_rho     = round(msign, 4),
        # VIF = 1/(1-rho^2): the direct answer to "your non-separability is just
        # collinearity". VIF > 5 is the usual red flag; > 10 is severe.
        VIF_at_mean_abs_rho = round(1 / (1 - mabs^2), 3),
        area_Mkm2   = round(sum(area_v[s]) / 1e6, 3),
        stringsAsFactors = FALSE)
    }
  }
  cat(sprintf("  %s scale %2d summarised\n", SUP, k))
}

res <- do.call(rbind, rows)
write.csv(res, out_csv, row.names = FALSE)

# ------------------------------------------------------------------------------
# Console report -- the headline
# ------------------------------------------------------------------------------

cat("\n================================================================\n")
cat("=== DOES SUPPLY-DEMAND COUPLING RISE WITH TIMESCALE?\n")
cat("===   area-weighted mean |rho(SM, AED)| over vegetated land\n")
cat("================================================================\n")

for (SUP in SUPPLIES) {
  cat(sprintf("\n  --- %s vs AED ---\n", SUP))
  cat(sprintf("  %-16s | %s\n", "aridity class",
              paste(sprintf("%9s", paste0(SCALES, "mo")), collapse = " ")))
  cat(sprintf("  %s\n", strrep("-", 68)))
  for (cl in c("ALL", ai_labels)) {
    vals <- sapply(SCALES, function(k) {
      s <- res$supply == SUP & res$scale_months == k & res$ai_class == cl
      if (!any(s)) return(NA_real_)
      mean(res$mean_abs_rho[s], na.rm = TRUE)      # across the 12 months
    })
    cat(sprintf("  %-16s | %s\n", cl,
                paste(sprintf("%9s", sprintf("%.3f", vals)), collapse = " ")))
  }
  vif <- sapply(SCALES, function(k) {
    s <- res$supply == SUP & res$scale_months == k & res$ai_class == "ALL"
    m <- mean(res$mean_abs_rho[s], na.rm = TRUE); 1 / (1 - m^2)
  })
  cat(sprintf("  %-16s | %s\n", "VIF (ALL)",
              paste(sprintf("%9s", sprintf("%.2f", vif)), collapse = " ")))
}

cat("\n  Interpretation: if |rho| rises from 1 to 12 months, supply and demand\n")
cat("  become less separable towards the timescales over which vegetation\n")
cat("  integrates water availability.\n")

cat(sprintf("\nCSV written: %s  (%d rows)\n", out_csv, nrow(res)))
cat("\n=== DONE ===\n")
