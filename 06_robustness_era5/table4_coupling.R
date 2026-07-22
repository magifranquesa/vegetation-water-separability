#!/usr/bin/env Rscript

# ==============================================================================
# Script: 06_robustness_era5/table4_coupling.R
#
# Project: Aridity and timescale bound the separability of vegetation water stress
#
# Supplementary Table 4 - supply-demand coupling by accumulation timescale,
# GLEAM v4.2a vs ERA5-Land.
#
#   rho(SMrz, AED) is computed between the two forcings, with no vegetation index
#   involved, so it is a property of the hydroclimatic products alone. Comparing
#   how it varies with accumulation timescale in GLEAM and in ERA5-Land tests
#   whether the timescale dependence of the coupling depends on the choice of
#   product, independently of the kNDVI analysis and the variance partitioning.
#
# METHOD - mirrors 04_analysis/06_supply_demand_coupling.R:
#   same linear detrending, same vectorised Spearman (ties = "first"), same
#   41 years, same aridity breaks, same area weighting. Only the input product
#   changes. The aridity index is the SAME CRU-derived file used for GLEAM, so
#   the stratification is identical and cannot itself explain a difference.
#
#   Restricted to SMrz, the supply variable behind Fig. 5c. SMs can be added by
#   extending SUPPLIES.
#
# Inputs (READ-ONLY):
#   outputs/intermediate/accumulated_era5/{SMrz,Ep}_era5land_monthly_1982-2022_scale_{k}.nc
#   outputs/intermediate/vegetation_mask_c1.tif
#   data/processed/aridity/ai_1982_2022_period.nc
#   outputs/tables/supply_demand_coupling.csv        (GLEAM, for the comparison)
#
# Outputs (refuses to overwrite):
#   outputs/intermediate/coupling_era5/coupling_SMrz_Ep_scale_{k}_era5.tif
#   outputs/tables/robustness_T4_coupling_era5_full.csv       (per month and class)
#   outputs/tables/robustness_T4_coupling_comparison.csv      (the table for the paper)
#
# Run from repo root:
#   Rscript 06_robustness_era5/table4_coupling.R
#
# Runtime: streams the accumulated NetCDFs in latitude blocks; tens of minutes,
# dominated by disk I/O.
#
# Dependencies: ncdf4, terra
# ==============================================================================

suppressPackageStartupMessages({ library(ncdf4); library(terra) })

config_file <- file.path("R", "config.R")
if (file.exists(config_file)) source(config_file)

dir_acc    <- file.path("outputs", "intermediate", "accumulated_era5")
file_mask  <- file.path("outputs", "intermediate", "vegetation_mask_c1.tif")
file_ai    <- file.path("data", "processed", "aridity", "ai_1982_2022_period.nc")
tables_dir <- if (exists("paths") && !is.null(paths$tables)) paths$tables else
  file.path("outputs", "tables")
file_gleam_csv <- file.path(tables_dir, "supply_demand_coupling.csv")

dir_out <- file.path("outputs", "intermediate", "coupling_era5")
dir.create(dir_out,    recursive = TRUE, showWarnings = FALSE)
dir.create(tables_dir, recursive = TRUE, showWarnings = FALSE)

out_csv_full <- file.path(tables_dir, "robustness_T4_coupling_era5_full.csv")
out_csv_cmp  <- file.path(tables_dir, "robustness_T4_coupling_comparison.csv")
for (f in c(out_csv_full, out_csv_cmp))
  if (file.exists(f)) stop("REFUSING TO OVERWRITE: ", f, call. = FALSE)

# ------------------------------------------------------------------------------
# Configuration (as in 04_analysis/06_supply_demand_coupling.R, ERA5-Land inputs)
# ------------------------------------------------------------------------------

SCALES   <- c(1, 3, 6, 9, 12)
DEMAND   <- "Ep"
SUPPLIES <- c("SMrz")

LAT_BLOCK <- 30
NYEAR     <- 41

ai_breaks <- c(-Inf, 0.03, 0.20, 0.50, 0.65, Inf)
ai_labels <- c("Hyper-arid", "Arid", "Semi-arid", "Dry sub-humid", "Humid")

acc_file <- function(v, k)
  file.path(dir_acc, sprintf("%s_era5land_monthly_1982-2022_scale_%d.nc", v, k))

for (k in SCALES) for (v in c(DEMAND, SUPPLIES))
  if (!file.exists(acc_file(v, k)))
    stop("Missing accumulated file: ", acc_file(v, k), call. = FALSE)
for (f in c(file_mask, file_ai, file_gleam_csv))
  if (!file.exists(f)) stop("Missing input: ", f, call. = FALSE)

# ------------------------------------------------------------------------------
# Vectorised detrend + Spearman
# ------------------------------------------------------------------------------

detrend_rows <- function(M) {
  n  <- ncol(M)
  tt <- seq_len(n) - (n + 1) / 2
  b  <- as.vector(M %*% tt) / sum(tt^2)
  M - rowMeans(M) - outer(b, tt)
}

row_varies <- function(M) {
  rmin <- do.call(pmin, as.data.frame(M))
  rmax <- do.call(pmax, as.data.frame(M))
  is.finite(rmax - rmin) & (rmax - rmin) > 0
}

rank_rows <- function(M) {
  m <- nrow(M); n <- ncol(M)
  v   <- as.vector(M)
  rid <- rep.int(seq_len(m), n)
  o   <- order(rid, v, method = "radix")
  r   <- numeric(m * n)
  r[o] <- rep.int(seq_len(n), m)
  dim(r) <- c(m, n)
  r
}

spearman_rows <- function(A, B) {
  RA <- rank_rows(A); RB <- rank_rows(B)
  ca <- RA - rowMeans(RA)
  cb <- RB - rowMeans(RB)
  out <- rowSums(ca * cb) / sqrt(rowSums(ca^2) * rowSums(cb^2))
  out[!is.finite(out)] <- NA_real_
  out
}

local({
  set.seed(1)
  A <- matrix(rnorm(20 * NYEAR), 20, NYEAR)
  B <- matrix(rnorm(20 * NYEAR), 20, NYEAR)
  ours <- spearman_rows(A, B)
  ref  <- vapply(seq_len(20), function(i)
    stats::cor(A[i, ], B[i, ], method = "spearman"), numeric(1))
  if (max(abs(ours - ref)) > 1e-10)
    stop("Spearman self-test FAILED - refusing to run.", call. = FALSE)
  message("Self-test passed: vectorised Spearman matches stats::cor to 1e-10.")
})

var_first <- function(f) {
  nc <- nc_open(f)
  nd <- vapply(nc$var, function(v) v$ndims, integer(1))
  v3 <- names(nc$var)[nd == 3]
  nc_close(nc)
  if (!length(v3)) stop("No 3-D variable found in ", basename(f), call. = FALSE)
  v3[1]
}

# ------------------------------------------------------------------------------
# Grid geometry
# ------------------------------------------------------------------------------

nc0  <- nc_open(acc_file(DEMAND, 1))
lon  <- nc0$dim[[grep("^(lon|longitude|x)$", names(nc0$dim), ignore.case = TRUE)[1]]]$vals
lat  <- nc0$dim[[grep("^(lat|latitude|y)$",  names(nc0$dim), ignore.case = TRUE)[1]]]$vals
ntim <- nc0$dim[[grep("^(time|t)$",          names(nc0$dim), ignore.case = TRUE)[1]]]$len
nc_close(nc0)

nlon <- length(lon); nlat <- length(lat)
lat_desc <- lat[1] > lat[nlat]

message(sprintf("ERA5-Land grid: %d lon x %d lat x %d time (lat %s)",
                nlon, nlat, ntim, if (lat_desc) "descending" else "ascending"))
if (ntim < NYEAR * 12)
  stop("Unexpected time length: ", ntim, call. = FALSE)

template <- rast(nrows = nlat, ncols = nlon,
                 xmin = min(lon) - 0.05, xmax = max(lon) + 0.05,
                 ymin = min(lat) - 0.05, ymax = max(lat) + 0.05,
                 crs = "EPSG:4326")

# ------------------------------------------------------------------------------
# Coupling per scale
# ------------------------------------------------------------------------------

for (k in SCALES) {

  out_tif <- setNames(
    file.path(dir_out, sprintf("coupling_%s_%s_scale_%d_era5.tif", SUPPLIES, DEMAND, k)),
    SUPPLIES)
  if (all(file.exists(out_tif))) { message(sprintf("  scale %2d: exists, skipping", k)); next }

  message(sprintf("\n=== scale %2d months | %s vs %s (ERA5-Land) ===",
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

        cb <- which(ok)
        j  <- ((cb - 1L) %/% nlon) + y0
        i  <- ((cb - 1L) %%  nlon) + 1L
        rho[[SUP]][(j - 1L) * nlon + i, m] <- r
      }
      rm(As); gc(verbose = FALSE)
    }

    cat(sprintf("    lat rows %4d-%4d done\n", y0, y0 + nyb - 1L))
    y0 <- y0 + nyb
    rm(Ad); gc(verbose = FALSE)
  }

  nc_close(nc_d); invisible(lapply(nc_s, nc_close))

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

# ------------------------------------------------------------------------------
# Summarise by aridity class
# ------------------------------------------------------------------------------

message("\n=== Summarising ERA5-Land coupling by timescale and aridity class ===")

msk <- rast(file_mask)
ai  <- rast(file_ai)[[1]]
if (is.na(crs(ai))) crs(ai) <- "EPSG:4326"
if (xmax(ai) > 180 + 1e-6) ai <- rotate(ai)

ref    <- rast(file.path(dir_out, sprintf("coupling_%s_%s_scale_1_era5.tif", SUPPLIES[1], DEMAND)))[[1]]
msk_v  <- as.vector(values(resample(msk, ref, method = "near")))
ai_v   <- as.vector(values(resample(ai,  ref, method = "bilinear")))
area_v <- as.vector(values(cellSize(ref, unit = "km", mask = FALSE)))
cls_v  <- cut(ai_v, breaks = ai_breaks, labels = ai_labels, right = FALSE)

veg <- !is.na(msk_v) & msk_v > 0 & !is.na(ai_v)
message(sprintf("  vegetated cells with valid AI: %d (%.2f x 10^6 km2)",
                sum(veg), sum(area_v[veg], na.rm = TRUE) / 1e6))

wmean <- function(x, w) { ok <- is.finite(x) & is.finite(w); sum(x[ok] * w[ok]) / sum(w[ok]) }

rows <- list()
for (SUP in SUPPLIES) for (k in SCALES) {
  r <- rast(file.path(dir_out, sprintf("coupling_%s_%s_scale_%d_era5.tif", SUP, DEMAND, k)))
  for (m in 1:12) {
    v   <- as.vector(values(r[[m]]))
    sel <- veg & is.finite(v)
    for (cl in c("ALL", ai_labels)) {
      s <- if (cl == "ALL") sel else sel & !is.na(cls_v) & cls_v == cl
      if (!any(s)) next
      mabs <- wmean(abs(v[s]), area_v[s])
      rows[[length(rows) + 1L]] <- data.frame(
        product = "ERA5-Land", supply = SUP, demand = DEMAND,
        scale_months = k, month = m, ai_class = cl,
        mean_abs_rho = round(mabs, 4),
        mean_rho     = round(wmean(v[s], area_v[s]), 4),
        area_Mkm2    = round(sum(area_v[s]) / 1e6, 3),
        stringsAsFactors = FALSE)
    }
  }
  cat(sprintf("  %s scale %2d summarised\n", SUP, k))
}

res_era5 <- do.call(rbind, rows)
write.csv(res_era5, out_csv_full, row.names = FALSE)
message("Full ERA5 table: ", out_csv_full)

# ------------------------------------------------------------------------------
# T4: the comparison table (mean over the 12 calendar months)
# ------------------------------------------------------------------------------

g <- read.csv(file_gleam_csv, stringsAsFactors = FALSE)
g <- g[g$supply == "SMrz" & g$demand == "Ep", ]

collapse <- function(df, label) {
  out <- aggregate(mean_abs_rho ~ ai_class + scale_months, data = df, FUN = mean)
  names(out)[names(out) == "mean_abs_rho"] <- label
  out
}

cmp <- merge(collapse(g, "GLEAM"), collapse(res_era5, "ERA5_Land"),
             by = c("ai_class", "scale_months"), all = TRUE)
cmp$GLEAM     <- round(cmp$GLEAM, 3)
cmp$ERA5_Land <- round(cmp$ERA5_Land, 3)
cmp$diff      <- round(cmp$ERA5_Land - cmp$GLEAM, 3)
cmp$ai_class  <- factor(cmp$ai_class, levels = c("ALL", ai_labels))
cmp <- cmp[order(cmp$ai_class, cmp$scale_months), ]

write.csv(cmp, out_csv_cmp, row.names = FALSE)

cat("\n================================================================\n")
cat("=== T4: |rho(SMrz, AED)| by accumulation timescale, both products\n")
cat("===   area-weighted, averaged over the 12 calendar months\n")
cat("================================================================\n\n")
print(cmp, row.names = FALSE)

cat("\n  >> The claim under test: |rho| RISES with timescale in BOTH products.\n")
cat("     If it does, the coupling is not an artefact of GLEAM's construction.\n")
cat("     Pay attention to Arid and Semi-arid: those carry the argument.\n")
cat(sprintf("\nComparison table: %s\n", out_csv_cmp))
cat("\n=== DONE ===\n")
