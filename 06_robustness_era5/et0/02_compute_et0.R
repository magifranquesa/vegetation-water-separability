# ==============================================================================
# 06_robustness_era5/et0/02_compute_et0.R
#
# Monthly FAO-56 Penman-Monteith reference evapotranspiration (ET0, the
# atmospheric evaporative demand) from ERA5-Land monthly-means.
#
# Inputs are configured in R/config_era5.R (paths in `era5`, variable names in
# `era5_var`); all share the ERA5-Land 0.1 deg global grid (valid_time = 504
# months, 1981-2022). The variables used here are:
#   Tmin, Tmax  -> t2m       [K]     (monthly means of daily min/max 2 m temp.)
#   dewpoint    -> d2m       [K]
#   ssrd        -> ssrd      [J m-2, mean daily accumulation]
#   wind        -> u10, v10  [m s-1]
#   elevation   -> z         [m]     (from 01_prepare_elevation.R, via the .rds)
#
# Output:
#   data/processed/et0_era5land/et0_era5land_monthly_1981-2022.nc  [mm month-1]
#
# Variable names inside each file are auto-detected. Run from the repository root:
#   Rscript 06_robustness_era5/et0/02_compute_et0.R
# ==============================================================================

library(ncdf4)

here <- file.path("06_robustness_era5", "et0")
source(file.path(here, "fao56_functions.R"))

# ------------------------------------------------------------------------------
# 1. Configuration
# ------------------------------------------------------------------------------
# The mix of CDS products below is DELIBERATE:
#   - ssrd is an accumulated variable and is taken from 'monthly averaged reanalysis
#     by hour of day' at 00:00, because the plain monthly product is corrupted from
#     September 2022 onwards (ECMWF known issue, 2024-01) for all accumulated fluxes.
#   - dewpoint and wind are instantaneous, so the plain monthly product is required:
#     the by-hour product at 00:00 would give only the 00 UTC snapshot.
source(file.path("R", "config_era5.R"))   # canonical paths, offsets, metadata

files <- list(
  tmin = era5$tmin,       # var 't2m'; series starts 1980 -> offset below
  tmax = era5$tmax,       # var 't2m'
  tdew = era5$dewpoint,
  wind = era5$wind,
  ssrd = era5$ssrd
)

# All inputs now start in January 1981, so month t is at index t in every file. The
# temperature files used to start in 1980 and were trimmed with CDO. The offsets are
# still read from the config and asserted below, because reading a file one year out
# raises no exception: the seasonal cycle would look perfect while every interannual
# anomaly came from the wrong year -- and that is exactly what the analysis measures.
t_off <- list(tmin = era5_t_offset$tmin, tmax = era5_t_offset$tmax,
              tdew = era5_t_offset$dewpoint, wind = era5_t_offset$wind,
              ssrd = era5_t_offset$ssrd)

elevation_rds <- file.path(here, "elevation_era5land.rds")

out_dir <- file.path("data", "processed", "et0_era5land")
out_nc  <- file.path(out_dir, "et0_era5land_monthly_1981-2022.nc")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# ssrd is the mean DAILY accumulation in J m-2 (by-hour product at 00:00, step 24 h,
# stamped on the 2nd of the month). Confirmed empirically: the July maximum reaches
# 31.5 MJ m-2, exactly the clear-sky ceiling (0.75 x Ra). FAO-56 wants MJ m-2 day-1.
SSRD_TO_MJ_M2_DAY <- 1e6

# ------------------------------------------------------------------------------
# 2. Open inputs, read coordinates and time axis
# ------------------------------------------------------------------------------
stopifnot(all(vapply(files, file.exists, logical(1))))
nc <- lapply(files, nc_open)

var <- lapply(nc, data_vars)                       # detected variable name(s) per file
message("Detected variables:")
for (k in names(var)) message("  ", k, ": ", paste(var[[k]], collapse = ", "))

lon <- ncvar_get(nc$tdew, "longitude")             # 3600, 0-360
lat <- ncvar_get(nc$tdew, "latitude")              # 1801, +90 -> -90
tsec <- ncvar_get(nc$tdew, "valid_time")           # seconds since 1970-01-01
tunits <- ncatt_get(nc$tdew, "valid_time", "units")$value

dates  <- as.Date(as.POSIXct(tsec, origin = "1970-01-01", tz = "UTC"))
years  <- as.integer(format(dates, "%Y"))
months <- as.integer(format(dates, "%m"))
ndays  <- mapply(days_in_month, years, months)
# Mid-month day-of-year, the representative J for monthly extraterrestrial radiation
midJ <- as.integer(format(as.Date(sprintf("%04d-%02d-15", years, months)), "%j"))

nlon <- length(lon); nlat <- length(lat); ntime <- length(tsec)

# Safety guard: the index-based reads below assume every input shares the
# reference grid (0-360 ERA5-Land, taken from d2m). Stop rather than silently
# misaligning variables if a file arrives on a different grid.
for (k in names(nc)) {
  lo <- ncvar_get(nc[[k]], "longitude")
  la <- ncvar_get(nc[[k]], "latitude")
  if (length(lo) != nlon || length(la) != nlat ||
      max(abs(lo - lon)) > 1e-2 || max(abs(la - lat)) > 1e-2) {
    stop(sprintf("Grid mismatch in '%s' (%s): expected the 0-360 ERA5-Land grid. Reproject/rotate it to match before running.",
                 k, files[[k]]))
  }
  # Each file must still reach month ntime once its own index offset is applied.
  have <- nc[[k]]$dim$valid_time$len
  need <- ntime + t_off[[k]]
  if (have < need) {
    stop(sprintf("'%s' has %d months but %d are needed (offset %d). Check the period.",
                 k, have, need, t_off[[k]]))
  }
}

# Decode a CF time axis ('days since ...' or 'seconds since ...') to Date.
nc_dates <- function(ncobj) {
  tv <- as.numeric(ncvar_get(ncobj, "valid_time"))
  u  <- ncatt_get(ncobj, "valid_time", "units")$value
  org <- substr(sub("^\\w+ since ", "", u), 1, 10)
  if (grepl("^seconds", u)) {
    as.Date(as.POSIXct(tv, origin = paste(org, "00:00:00"), tz = "UTC"))
  } else {
    as.Date(org) + tv
  }
}

# Every input must resolve to the SAME calendar month once its offset is applied.
# Assert it rather than trust it: a one-year slip raises no exception and would leave
# the seasonal cycle looking perfect while every interannual anomaly came from the
# wrong year. Also check the LAST month, which catches a file of the right length but
# the wrong period.
d_ref <- nc_dates(nc$tdew)
for (k in names(nc)) {
  d_k <- nc_dates(nc[[k]])
  for (i in c(1L, ntime)) {
    got <- format(d_k[i + t_off[[k]]], "%Y-%m")
    exp <- format(d_ref[i], "%Y-%m")
    if (got != exp) {
      stop(sprintf("Time misalignment in '%s': month %d resolves to %s, expected %s.",
                   k, i, got, exp))
    }
  }
}
message(sprintf("Time alignment OK: %s .. %s (%d months) in all %d inputs.",
                format(d_ref[1], "%Y-%m"), format(d_ref[ntime], "%Y-%m"),
                ntime, length(nc)))

# Elevation and per-latitude fields, broadcast to the [lon x lat] grid
z <- readRDS(elevation_rds)
stopifnot(nrow(z) == nlon, ncol(z) == nlat)
lat_rad <- lat * pi / 180

# ------------------------------------------------------------------------------
# 3. Define output NetCDF (compressed, written month by month)
# ------------------------------------------------------------------------------
londim  <- ncdim_def("longitude", "degrees_east",  lon)
latdim  <- ncdim_def("latitude",  "degrees_north", lat)
timedim <- ncdim_def("valid_time", tunits, tsec, unlim = TRUE)

FILL <- 1e20
et0_var <- ncvar_def(
  name       = "et0",
  units      = "mm month-1",
  dim        = list(londim, latdim, timedim),
  missval    = FILL,
  longname   = "FAO-56 Penman-Monteith reference evapotranspiration",
  prec       = "float",
  compression = 4
)
ncout <- nc_create(out_nc, et0_var, force_v4 = TRUE)
ncatt_put(ncout, "et0", "method", "FAO-56 Penman-Monteith (Allen et al. 1998)")
ncatt_put(ncout, "et0", "source", "ERA5-Land monthly-means (reanalysis-era5-land-monthly-means)")
ncatt_put(ncout, 0, "history", paste("Created", Sys.time(),
                                      "by 06_robustness_era5/et0/02_compute_et0.R"))

# ------------------------------------------------------------------------------
# 4. Loop over months
# ------------------------------------------------------------------------------
read_month <- function(ncobj, vname, t) {
  ncvar_get(ncobj, vname, start = c(1, 1, t), count = c(-1, -1, 1))
}

for (t in seq_len(ntime)) {
  Tmin <- read_month(nc$tmin, var$tmin[1], t) - 273.15
  Tmax <- read_month(nc$tmax, var$tmax[1], t) - 273.15
  Tdew <- read_month(nc$tdew, var$tdew[1], t) - 273.15
  Rs   <- read_month(nc$ssrd, var$ssrd[1], t) / SSRD_TO_MJ_M2_DAY
  u10  <- read_month(nc$wind, var$wind[1], t)
  v10  <- read_month(nc$wind, var$wind[2], t)

  u2 <- wind_10m_to_2m(u10, v10)

  # Extraterrestrial radiation varies with latitude and month only: build a
  # per-latitude vector and broadcast across longitude.
  Ra_vec <- extraterrestrial_radiation(lat_rad, midJ[t])
  Ra <- matrix(Ra_vec, nrow = nlon, ncol = nlat, byrow = TRUE)

  et0_day <- et0_penman_monteith(Tmax, Tmin, Tdew, u2, Rs, Ra, z)  # mm day-1
  et0_mon <- et0_day * ndays[t]                                    # mm month-1

  # Ocean / missing cells arrive as NaN; write them as the fill value via NA.
  et0_mon[!is.finite(et0_mon)] <- NA

  # Diagnose BEFORE writing: ncvar_put substitutes NA by the fill value (1e20) and,
  # because R may modify the vector in place, reading et0_mon afterwards reports the
  # fill value as if it were the maximum ET0.
  if (t == 1) {
    ssrd_raw <- read_month(nc$ssrd, var$ssrd[1], t)
    message(sprintf("Diagnostics month 1 (%s):", format(dates[1], "%Y-%m")))
    message(sprintf("  raw ssrd [%.3g, %.3g] J m-2  (mean daily accumulation)",
                    min(ssrd_raw, na.rm = TRUE), max(ssrd_raw, na.rm = TRUE)))
    # Rs must not exceed the clear-sky ceiling Rso = (0.75 + 2e-5 z) Ra. In January
    # the maximum sits on the Antarctic plateau: 24 h of sun, ~3000 m, so ~36 is right.
    message(sprintf("  Rs  [%.1f, %.1f] | Ra [%.1f, %.1f] MJ m-2 day-1  (Rs <= ~0.81 Ra)",
                    min(Rs, na.rm = TRUE), max(Rs, na.rm = TRUE),
                    min(Ra, na.rm = TRUE), max(Ra, na.rm = TRUE)))
    message(sprintf("  ET0 [%.1f, %.1f] mm/month  (expect a max of ~150-250, never >400)",
                    min(et0_mon, na.rm = TRUE), max(et0_mon, na.rm = TRUE)))
    message(sprintf("  land cells: %.1f%%", 100 * mean(!is.na(et0_mon))))
  }

  ncvar_put(ncout, "et0", et0_mon, start = c(1, 1, t), count = c(-1, -1, 1))
  if (t %% 24 == 0 || t == ntime)
    message(sprintf("  %d / %d months (%s)", t, ntime, format(dates[t], "%Y-%m")))
}

# ------------------------------------------------------------------------------
# 5. Clean up
# ------------------------------------------------------------------------------
nc_close(ncout)
lapply(nc, nc_close)
message("Done. Wrote ", out_nc)
