#!/usr/bin/env Rscript

# ==============================================================================
# Script: 06_robustness_era5/checks/verify_tmin_tmax.R
#
# Project: Aridity and timescale bound the separability of vegetation water stress
#
# Purpose (QC — read-only, blocking check before running the FAO-56 ET0 module):
#   The two temperature files both declare the SAME variable ('t2m', GRIB paramId
#   167 = mean 2 m temperature) and the same long_name. Their file names claim one
#   is the minimum and the other the maximum. If both actually held the mean, FAO-56
#   would silently receive Tmax = Tmin = Tmean: zero diurnal range, wrong saturation
#   vapour pressure and wrong net longwave radiation, with no error raised anywhere.
#   Metadata cannot settle this, so the test is physical.
#
#   Also verifies the TIME OFFSET. These files span 1980-01..2022-12 (516 months)
#   while every other input starts in 1981-01 (504 months). Our month t therefore
#   lives at index t + 12. Reading them with start = c(1,1,t) would silently take the
#   temperature of the SAME calendar month ONE YEAR EARLIER: the seasonal cycle would
#   look perfect and only the interannual variability -- the very thing the analysis
#   uses -- would be wrong.
#
#   Expected diurnal temperature range (Tmax - Tmin), monthly means, July:
#     Sahara ~18-22 K | Iberia ~13-16 K | Siberia ~10-13 K | Amazon ~9-11 K
#   Near-zero everywhere would mean both files hold the mean.
#   ~25-35 K would mean they hold MONTHLY extremes, not means of DAILY extremes,
#   which is not what FAO-56 expects either.
#
# Inputs (READ-ONLY):
#   data/raw/ERA5_land_monthly/mnt2m_monthly_mean_1980-1_2022-12.nc   (var t2m)
#   data/raw/ERA5_land_monthly/mxt2m_monthly_mean_1980-1_2022-12.nc   (var t2m)
#   data/raw/ERA5_land_monthly/era5-land_2m_dewpoint_temperature.nc   (var d2m)
#
# Run from repo root:
#   Rscript 06_robustness_era5/checks/verify_tmin_tmax.R
#
# Dependencies: ncdf4
# ==============================================================================

suppressPackageStartupMessages(library(ncdf4))

source(file.path("R", "config_era5.R"))   # canonical paths

f_mn <- era5$tmin
f_mx <- era5$tmax
f_td <- era5$dewpoint
for (f in c(f_mn, f_mx, f_td)) {
  if (!file.exists(f)) stop("Missing: ", f, call. = FALSE)
}

T_OFFSET <- era5_t_offset$tmin        # 0 since the leading 1980 year was trimmed
t_of <- function(y, m) (y - 1981) * 12 + m
probe <- c(Jan2000 = t_of(2000, 1), Jul2000 = t_of(2000, 7))

# Boxes in the files' 0-360 longitude convention
boxes <- list(
  "Sahara (hyper-arid)" = list(lon = c(10, 30),   lat = c(20, 28)),
  "Iberia"              = list(lon = c(352, 358), lat = c(37, 43)),
  "Siberia (boreal)"    = list(lon = c(90, 110),  lat = c(55, 65)),
  "Amazon (forest)"     = list(lon = c(290, 305), lat = c(-8, 0))
)

get_slice <- function(f, v, t) {
  nc <- nc_open(f); on.exit(nc_close(nc))
  ncvar_get(nc, v, start = c(1, 1, t), count = c(-1, -1, 1))
}

cat("\n######## Are mnt2m / mxt2m really Tmin and Tmax? ########\n")

# ------------------------------------------------------------------------------
# 1. Time axis and the +12 offset
# ------------------------------------------------------------------------------
cat("\n=== 1. Time axis ===\n")
nc <- nc_open(f_mn)
tv <- as.numeric(ncvar_get(nc, "valid_time"))
tu <- ncatt_get(nc, "valid_time", "units")$value
lon <- as.numeric(ncvar_get(nc, "longitude")); lat <- as.numeric(ncvar_get(nc, "latitude"))
nc_close(nc)
origin <- sub("^days since ", "", tu)
dates <- as.Date(tv, origin = as.Date(substr(origin, 1, 10)))
cat(sprintf("  n = %d | first = %s | last = %s\n", length(tv), dates[1], dates[length(dates)]))
cat(sprintf("  index 1 -> %s   (must be 1981-01, i.e. our month 1)\n", format(dates[1], "%Y-%m")))
stopifnot(length(tv) == 504, format(dates[1], "%Y-%m") == "1981-01")
cat(sprintf("  -> 504 months, aligned with every other input. T_OFFSET = %d.\n", T_OFFSET))
cat(sprintf("  grid: %d lon [%.1f, %.1f] x %d lat [%.1f, %.1f]\n",
            length(lon), min(lon), max(lon), length(lat), min(lat), max(lat)))

# ------------------------------------------------------------------------------
# 2. Diurnal temperature range
# ------------------------------------------------------------------------------
cat("\n=== 2. Diurnal temperature range, Tmax - Tmin (K) ===\n")
for (pn in names(probe)) {
  t <- probe[[pn]] + T_OFFSET
  mn <- get_slice(f_mn, "t2m", t)
  mx <- get_slice(f_mx, "t2m", t)
  dtr <- mx - mn
  ok <- is.finite(dtr)
  q <- quantile(dtr[ok], c(0.01, 0.5, 0.99), names = FALSE)
  cat(sprintf("\n  -- %s (file index %d)\n", pn, t))
  cat(sprintf("     DTR: min %6.2f | p1 %6.2f | med %6.2f | p99 %6.2f | max %6.2f\n",
              min(dtr[ok]), q[1], q[2], q[3], max(dtr[ok])))
  cat(sprintf("     Tmax <= Tmin in %.3f%% of land cells (should be ~0)\n", 100 * mean(dtr[ok] <= 0)))
  cat(sprintf("     mean Tmin %6.2f degC | mean Tmax %6.2f degC\n",
              mean(mn[ok]) - 273.15, mean(mx[ok]) - 273.15))
  if (pn == "Jul2000") {
    cat("     by region (July):\n")
    for (nm in names(boxes)) {
      b <- boxes[[nm]]
      i <- which(lon >= b$lon[1] & lon <= b$lon[2]); j <- which(lat >= b$lat[1] & lat <= b$lat[2])
      x <- dtr[i, j]; x <- x[is.finite(x)]
      cat(sprintf("       %-22s DTR = %5.2f K\n", nm, mean(x)))
    }
  }
}

# ------------------------------------------------------------------------------
# 3. Cross-check against dewpoint (Tdew must not exceed Tmax)
# ------------------------------------------------------------------------------
cat("\n=== 3. Consistency with dewpoint, Jul 2000 ===\n")
t <- probe[["Jul2000"]]
mn <- get_slice(f_mn, "t2m", t + T_OFFSET)
mx <- get_slice(f_mx, "t2m", t + T_OFFSET)
td <- get_slice(f_td, "d2m", t)                 # d2m file starts in 1981 -> no offset
ok <- is.finite(mn) & is.finite(mx) & is.finite(td)
cat(sprintf("  Tdew > Tmax : %.3f%% of cells (physically almost impossible)\n", 100 * mean(td[ok] > mx[ok])))
cat(sprintf("  Tdew > Tmean: %.2f%% of cells (some is normal in very humid places)\n",
            100 * mean(td[ok] > (mn[ok] + mx[ok]) / 2)))
cat(sprintf("  Tdew > Tmin : %.2f%% of cells (common; nights reach saturation)\n", 100 * mean(td[ok] > mn[ok])))

# ------------------------------------------------------------------------------
# 4. Verdict
# ------------------------------------------------------------------------------
cat("\n=== Verdict ===\n")
cat("  DTR ~ 0 everywhere        -> BOTH files hold the mean. Do NOT use them.\n")
cat("  DTR ~ 8-20 K, region-wise -> they are means of DAILY min/max. Correct for FAO-56.\n")
cat("  DTR ~ 25-35 K             -> they are MONTHLY extremes, not means of daily extremes.\n")
cat("                               FAO-56 would overestimate the vapour-pressure deficit.\n")
cat("  A large 'Tdew > Tmax' fraction would mean the two files are swapped.\n")
cat("\nDone (read-only).\n")
