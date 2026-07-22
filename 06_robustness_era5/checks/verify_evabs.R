#!/usr/bin/env Rscript

# ==============================================================================
# Script: 06_robustness_era5/checks/verify_evabs.R
#
# Project: Aridity and timescale bound the separability of vegetation water stress
#
# Purpose (EXPLORATORY QC — read-only, fast):
#   Test empirically WHICH physical quantity lives in
#   "Era5land_Evaporation from bare soil.nc" (GRIB paramId 228101).
#
#   ERA5-Land documented known issue: the values of three evaporation variables
#   are ROTATED among each other:
#       evabs (228101, "bare soil")   holds  evavt (transpiration)
#       evaow (228102, "open water")  holds  evabs (bare soil)
#       evavt (228103, "transpiration") holds evaow (open water)
#   Metadata cannot reveal this (paramId/long_name are unchanged), so the test
#   must be PHYSICAL.
#
#   Discriminating signatures (annual-ish means, positive fluxes, mm/day):
#     TRANSPIRATION : large over dense forest (Amazon/Congo ~2.5-3.5), ~0 over
#                     Sahara; ratio to total E over forest ~0.6-0.8.
#     BARE SOIL EVAP: small over dense forest (~0.1-0.5, ratio ~0.05-0.15),
#                     it is the dominant component of E over deserts.
#     OPEN WATER    : ~0 almost everywhere on land.
#
#   The same numbers also settle the ACCUMULATION BASIS: this file is stamped on
#   the LAST day of each month while pev/e are stamped on the 2nd. If the values
#   are ~30x larger than the expected mm/day, the file holds a MONTHLY total.
#
# Inputs (READ-ONLY; mtimes confirmed unchanged at the end):
#   data/raw/ERA5_land_monthly_by_hour/Era5land_Evaporation from bare soil.nc
#   data/raw/ERA5_land_monthly_by_hour/Era5land_total_Evaporation.nc
#
# Run from repo root:
#   Rscript 06_robustness_era5/checks/verify_evabs.R
#
# Dependencies: ncdf4
# ==============================================================================

suppressPackageStartupMessages(library(ncdf4))

source(file.path("R", "config_era5.R"))   # canonical paths

f_x <- era5$transpiration   # variable under test: labelled 'evabs' (see MANIFEST)
f_e <- era5$totalevap       # total evaporation
for (f in c(f_x, f_e)) {
  if (!file.exists(f)) stop("Required input not found: ", f, call. = FALSE)
}
mt_before <- file.info(c(f_x, f_e))$mtime

t_of <- function(year, month) (year - 1981) * 12 + month
probe_t <- c(t_of(2000, 1), t_of(2000, 4), t_of(2000, 7), t_of(2000, 10))

# Boxes (lon in the file's 0-360 convention)
boxes <- list(
  "Amazon (dense forest)" = list(lon = c(290, 305), lat = c(-8, 0)),
  "Congo  (dense forest)" = list(lon = c(15, 28),   lat = c(-4, 4)),
  "Sahara (hyper-arid)"   = list(lon = c(10, 30),   lat = c(20, 28)),
  "Sahel  (semi-arid)"    = list(lon = c(0, 15),    lat = c(12, 16)),
  "Siberia (boreal)"      = list(lon = c(90, 110),  lat = c(55, 65))
)

nc_x <- nc_open(f_x); nc_e <- nc_open(f_e)
vx <- names(nc_x$var)[vapply(names(nc_x$var), function(v) nc_x$var[[v]]$ndims == 3, logical(1))][1]
ve <- names(nc_e$var)[vapply(names(nc_e$var), function(v) nc_e$var[[v]]$ndims == 3, logical(1))][1]
lon <- as.numeric(ncvar_get(nc_x, "longitude"))
lat <- as.numeric(ncvar_get(nc_x, "latitude"))

idx_range <- function(v, lo, hi) which(v >= lo & v <= hi)

# Accumulate the 4 probe months (mm/day, sign-flipped to positive fluxes)
Xsum <- NULL; Esum <- NULL
for (t in probe_t) {
  X <- -ncvar_get(nc_x, vx, start = c(1, 1, t), count = c(-1, -1, 1)) * 1000
  E <- -ncvar_get(nc_e, ve, start = c(1, 1, t), count = c(-1, -1, 1)) * 1000
  if (is.null(Xsum)) { Xsum <- X; Esum <- E } else { Xsum <- Xsum + X; Esum <- Esum + E }
}
nc_close(nc_x); nc_close(nc_e)
Xm <- Xsum / length(probe_t); Em <- Esum / length(probe_t)   # 4-month mean

cat("\n############ Which variable is in the 'bare soil' file? ############\n")
cat(sprintf("  Variable read: '%s'   Total evaporation: '%s'\n", vx, ve))
cat("  Values are the mean of Jan/Apr/Jul/Oct 2000, sign-flipped, in mm/day.\n\n")
cat(sprintf("  %-24s %10s %10s %8s\n", "Region", "X", "E", "X/E"))
res <- list()
for (nm in names(boxes)) {
  b <- boxes[[nm]]
  i <- idx_range(lon, b$lon[1], b$lon[2]); j <- idx_range(lat, b$lat[1], b$lat[2])
  x <- Xm[i, j]; e <- Em[i, j]
  ok <- is.finite(x) & is.finite(e)
  mx <- mean(x[ok]); me <- mean(e[ok])
  res[[nm]] <- c(X = mx, E = me, ratio = mx / me)
  cat(sprintf("  %-24s %10.3f %10.3f %8.3f\n", nm, mx, me, mx / me))
}

ok <- is.finite(Xm) & is.finite(Em) & Em > 0.2
cat(sprintf("\n  Global land (E > 0.2 mm/day): median X/E = %.3f   mean X = %.3f mm/day\n",
            median(Xm[ok] / Em[ok]), mean(Xm[is.finite(Xm)])))

# ------------------------------------------------------------------------------
# Verdict
# ------------------------------------------------------------------------------
amz <- res[["Amazon (dense forest)"]]; sah <- res[["Sahara (hyper-arid)"]]
cgo <- res[["Congo  (dense forest)"]]
glob_ratio <- median(Xm[ok] / Em[ok])
cat("\n=== Verdict ===\n")
# The DESERT is the clean discriminator: with no vegetation, transpiration is exactly 0,
# whereas bare-soil evaporation would account for essentially all of E.
if (max(amz["X"], cgo["X"]) > 20) {
  cat("  !! Forest X is far above a plausible daily flux -> the file looks like a MONTHLY\n")
  cat("     accumulation (consistent with its last-day-of-month timestamp). Divide by the\n")
  cat("     number of days, or re-download with the same request as pev/e.\n")
} else if (sah["ratio"] < 0.10 && glob_ratio > 0.35 && max(amz["X"], cgo["X"]) > 1.0) {
  cat("  -> TRANSPIRATION (evavt). Zero over the Sahara (no vegetation), dominant component\n")
  cat("     of E over forest, global median X/E ~0.5. Magnitudes are daily fluxes.\n")
  cat("     Confirms the documented ERA5-Land variable rotation. USE THIS FILE AS Et.\n")
} else if (sah["ratio"] > 0.5) {
  cat("  -> BARE SOIL EVAPORATION (evabs, as labelled): it accounts for most of E in the\n")
  cat("     desert. The rotation issue would NOT affect this download, and a separate\n")
  cat("     variable would be needed for Et.\n")
} else {
  cat("  -> Ambiguous. Inspect the printed numbers and the maps before deciding.\n")
}
cat("\n  Reference expectations (mm/day, annual-ish):\n")
cat("    transpiration : Amazon 2.5-3.5, ratio 0.6-0.8 | Sahara ~0.00\n")
cat("    bare soil     : Amazon 0.1-0.5, ratio 0.05-0.15 | dominant fraction of E in deserts\n")
cat("    open water    : ~0 nearly everywhere on land\n")

cat("\n=== File protection confirmation (inputs READ-ONLY) ===\n")
mt_after <- file.info(c(f_x, f_e))$mtime
for (i in 1:2)
  cat(sprintf("  %-42s : %s\n", basename(c(f_x, f_e)[i]),
              if (identical(mt_before[i], mt_after[i])) "UNCHANGED" else "*** CHANGED ***"))
cat("\nDone.\n")
