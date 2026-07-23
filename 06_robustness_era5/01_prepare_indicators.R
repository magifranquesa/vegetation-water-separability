#!/usr/bin/env Rscript

# ==============================================================================
# Script: 06_robustness_era5/01_prepare_indicators.R
#
# Project: Aridity and timescale bound the separability of vegetation water stress
#
# Purpose:
#   Turn the raw ERA5-Land monthly download into the five project indicators, on
#   the GLEAM/kNDVI grid, ready to feed the correlation pipeline as a robustness
#   replication of the GLEAM-based analysis.
#
#   Ep (AED)  <- pev            (or the FAO-56 ET0 from et0_era5land/, recommended)
#   E         <- e              (total evaporation, only used to build ED)
#   Et        <- "evabs" file   (ERA5-Land known issue: 228101 holds TRANSPIRATION;
#                                verified in 06_robustness_era5/checks/verify_evabs.R)
#   ED        =  E - Ep         (project definition, 01_compute_evaporation_deficit.sh)
#   SMs       <- swvl1          (0-7 cm)
#   SMrz      <- depth-weighted mean of swvl1..3 (7/21/72 cm ~ 0-100 cm, GLEAM root zone)
#
#   Conversions applied to the accumulated fluxes:
#     ERA5 sign is positive DOWNWARD, so evaporation is negative -> flip sign.
#     Values are the mean DAILY accumulation in m (stream="mnth" at 00:00, step 24 h),
#     so:  mm/month = -raw * 1000 * days_in_month.
#   Soil moisture is instantaneous (m3 m-3): no sign flip, no time scaling.
#
#   REGRID (exact, no interpolation library):
#     ERA5-Land is node-registered  (lon 0..359.9, lat 90..-90, 1801 x 3600).
#     GLEAM/kNDVI is cell-centred   (lon -179.95..179.95, lat 89.95..-89.95, 1800 x 3600).
#     Every GLEAM cell centre falls EXACTLY on the midpoint of an ERA5-Land 2x2 node
#     block, so bilinear resampling reduces to a plain 2x2 average after rotating the
#     longitude to -180/180. The grid assumption is asserted at runtime.
#
# Inputs (READ-ONLY):
#   data/raw/ERA5_land_monthly_by_hour/Era5land_Potential evaporation.nc
#   data/raw/ERA5_land_monthly_by_hour/Era5land_total_Evaporation.nc
#   data/raw/ERA5_land_monthly_by_hour/Era5land_Evaporation from bare soil.nc  (= Et)
#   data/raw/ERA5_land_monthly/era5-land_soil_layers.nc   (swvl1..3; product 'moda')
#   outputs/intermediate/kndvi/kndvi.nc              (target-grid template)
#   outputs/intermediate/vegetation_mask_c1.tif      (optional mask)
#
# Outputs (one NetCDF per indicator):
#   data/processed/era5land_indicators/<VAR>_era5land_monthly_1981-2022.nc
#
# Run from repo root:
#   Rscript 06_robustness_era5/01_prepare_indicators.R
#
# Dependencies: ncdf4 (+ terra only if APPLY_VEG_MASK = TRUE)
# ==============================================================================

suppressPackageStartupMessages(library(ncdf4))

# ------------------------------------------------------------------------------
# Configuration
# ------------------------------------------------------------------------------
source(file.path("R", "config_era5.R"))   # canonical paths, offsets, metadata

F_E  <- era5$totalevap
F_ET <- era5$transpiration        # variable inside is 'evabs'; it holds transpiration
F_PEV <- NA_character_            # ERA5-Land pev was deleted: biased, and AED = FAO-56 ET0

# Soil moisture comes from product 'monthly_averaged_reanalysis' (stream='moda'), the
# true monthly mean for instantaneous variables. The by-hour copy was the 00 UTC
# snapshot only and has been deleted.
SOIL_FILE <- era5$soillayers

# AED source. ERA5-Land's own 'pev' is biased high (Iberia 9.8 mm/day in July against
# ~6-7 from FAO-56, with extremes of 72 mm/day) and has been deleted, so the only
# option is the FAO-56 Penman-Monteith product built in et0_era5land/.
AED_SOURCE <- "et0"
F_ET0 <- file.path("data", "processed", "et0_era5land", "et0_era5land_monthly_1981-2022.nc")

# Which indicators to build. Et/E/SMs/SMrz are final regardless of the AED choice and
# are already done; Ep and ED depend on AED_SOURCE, so they come last.
# Already built: c("Et", "E", "SMs", "SMrz")
VARS <- c("Ep", "ED")
# The GLEAM inputs feed the correlation pipeline UNMASKED (02_compute_correlations/
# applies no mask; C1 is applied downstream, in the figures and area summaries). To keep
# the ERA5 replication a clean drop-in swap, leave this FALSE so both sources share the
# same domain. Set TRUE only if you deliberately want pre-masked indicator files.
APPLY_VEG_MASK <- FALSE
F_TEMPLATE <- file.path("outputs", "intermediate", "kndvi", "kndvi.nc")
F_VEG      <- file.path("outputs", "intermediate", "vegetation_mask_c1.tif")

out_dir <- file.path("data", "processed", "era5land_indicators")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# ------------------------------------------------------------------------------
# Time axis (canonical: 1981-01 .. 2022-12, stamped on the 1st)
# ------------------------------------------------------------------------------
dates  <- seq(as.Date("1981-01-01"), by = "month", length.out = 504)
tsec   <- as.numeric(as.POSIXct(paste(dates, "00:00:00"), tz = "UTC"))
ndays  <- as.integer(diff(c(dates, seq(dates[504], by = "month", length.out = 2)[2])))
ntime  <- length(dates)

# ------------------------------------------------------------------------------
# Grid: source (ERA5-Land) -> target (GLEAM), verified against the template
# ------------------------------------------------------------------------------
NLON_S <- 3600L; NLAT_S <- 1801L
NLON_T <- 3600L; NLAT_T <- 1800L

lon_src <- (0:(NLON_S - 1)) * 0.1                  # 0 .. 359.9
lat_src <- 90 - (0:(NLAT_S - 1)) * 0.1             # 90 .. -90
rot     <- c(1801:3600, 1:1800)                    # rotate 0-360 -> -180/180
lon_rot <- ifelse(lon_src >= 180, lon_src - 360, lon_src)[rot]   # -180 .. 179.9

lon_tgt <- lon_rot + 0.05                          # -179.95 .. 179.95
lat_tgt <- lat_src[1:NLAT_T] - 0.05                # 89.95 .. -89.95

nc_t <- nc_open(F_TEMPLATE)
lon_ref <- as.numeric(ncvar_get(nc_t, "lon")); lat_ref <- as.numeric(ncvar_get(nc_t, "lat"))
nc_close(nc_t)
stopifnot(length(lon_ref) == NLON_T, length(lat_ref) == NLAT_T)
if (max(abs(lon_tgt - lon_ref)) > 1e-3 || max(abs(lat_tgt - lat_ref)) > 1e-3) {
  stop("Grid assumption failed: derived target grid does not match ", F_TEMPLATE, call. = FALSE)
}
message("Grid check OK: ERA5-Land 2x2 node average lands exactly on the GLEAM cell centres.")

# Exact regrid: rotate longitude, then average the four surrounding nodes.
to_gleam <- function(m) {                          # m: [3600 lon(0-360), 1801 lat]
  R  <- m[rot, ]
  sh <- c(2:NLON_S, 1L)                            # circular +1 in longitude
  A <- R[, 1:NLAT_T];      B <- R[sh, 1:NLAT_T]
  C <- R[, 2:NLAT_S];      D <- R[sh, 2:NLAT_S]
  n_ok <- (!is.na(A)) + (!is.na(B)) + (!is.na(C)) + (!is.na(D))
  A[is.na(A)] <- 0; B[is.na(B)] <- 0; C[is.na(C)] <- 0; D[is.na(D)] <- 0
  out <- (A + B + C + D) / n_ok
  out[n_ok == 0] <- NA
  out                                              # [3600 lon, 1800 lat], lat north-first
}

# ------------------------------------------------------------------------------
# Optional C1 vegetation mask, brought to the same [lon, lat] orientation
# ------------------------------------------------------------------------------
veg_mat <- NULL
if (APPLY_VEG_MASK) {
  suppressPackageStartupMessages(library(terra))
  v <- rast(F_VEG)
  e <- as.vector(ext(v)); r <- res(v)               # e = xmin, xmax, ymin, ymax
  stopifnot(ncol(v) == NLON_T,
            abs(r[1] - 0.1) < 1e-6, abs(r[2] - 0.1) < 1e-6,
            abs(e[1] - (-180)) < 1e-6)
  # The C1 mask sits on the GLEAM grid but is cropped in latitude (it stops at 65 S,
  # excluding Antarctica), so it covers only the top rows of the 1800-row target grid.
  # Drop it into place; everything outside it stays masked out.
  row0 <- as.integer(round((90 - e[4]) / 0.1))     # target rows skipped above the mask
  rows <- (row0 + 1):(row0 + nrow(v))
  stopifnot(min(rows) >= 1, max(rows) <= NLAT_T)
  m <- t(as.matrix(v, wide = TRUE))                # [lon, lat], lat north-first
  veg_mat <- matrix(FALSE, NLON_T, NLAT_T)
  veg_mat[, rows] <- !is.na(m) & m > 0
  message(sprintf("Vegetation mask C1: %d x %d -> target rows %d-%d (%.1f%% of grid vegetated).",
                  ncol(v), nrow(v), min(rows), max(rows), 100 * mean(veg_mat)))
}

# ------------------------------------------------------------------------------
# Helpers
# ------------------------------------------------------------------------------
main_var <- function(nc) {
  v <- names(nc$var)
  v[vapply(v, function(x) nc$var[[x]]$ndims == 3, logical(1))][1]
}
read_slice <- function(nc, v, t) ncvar_get(nc, v, start = c(1, 1, t), count = c(-1, -1, 1))

FILL <- 1e20
make_out <- function(varname, units, longname) {
  f <- file.path(out_dir, sprintf("%s_era5land_monthly_1981-2022.nc", varname))
  if (file.exists(f)) stop("Output exists (not overwriting): ", f, call. = FALSE)
  d_lon <- ncdim_def("lon", "degrees_east", lon_tgt)
  d_lat <- ncdim_def("lat", "degrees_north", lat_tgt)
  d_t   <- ncdim_def("time", "seconds since 1970-01-01", tsec, unlim = TRUE)
  v <- ncvar_def(varname, units, list(d_lon, d_lat, d_t), FILL, longname, "float", compression = 4)
  nc <- nc_create(f, v, force_v4 = TRUE)
  ncatt_put(nc, 0, "source", "ERA5-Land monthly (stream=mnth for fluxes, moda for soil)")
  ncatt_put(nc, 0, "regrid", "exact 2x2 node average onto the GLEAM/kNDVI cell centres")
  list(nc = nc, name = varname, path = f)
}
put_month <- function(o, mat, t) {
  if (!is.null(veg_mat)) mat[!veg_mat] <- NA
  mat[!is.finite(mat)] <- NA
  ncvar_put(o$nc, o$name, mat, start = c(1, 1, t), count = c(-1, -1, 1))
}

# ------------------------------------------------------------------------------
# Open sources / create outputs
# ------------------------------------------------------------------------------
need_flux <- any(c("Ep", "E", "Et", "ED") %in% VARS)
need_soil <- any(c("SMs", "SMrz") %in% VARS)
need_Ep   <- any(c("Ep", "ED") %in% VARS)
need_E    <- any(c("E", "ED") %in% VARS)

src <- list()
if (need_Ep) {
  if (AED_SOURCE == "et0") {
    if (!file.exists(F_ET0)) stop("AED_SOURCE='et0' but not found: ", F_ET0, call. = FALSE)
    src$ep <- nc_open(F_ET0)
    message("AED source: FAO-56 ET0 (already mm/month).")
  } else {
    stop("AED_SOURCE='pev' is not supported: ERA5-Land potential evaporation is biased ",
         "high. Use AED_SOURCE='et0' and build the FAO-56 product first with ",
         "06_robustness_era5/et0/02_compute_et0.R.", call. = FALSE)
  }
}
if (need_E)  { stopifnot(file.exists(F_E));  src$e  <- nc_open(F_E) }
if ("Et" %in% VARS) { stopifnot(file.exists(F_ET)); src$et <- nc_open(F_ET) }
if (need_soil) {
  if (!file.exists(SOIL_FILE)) {
    stop("Soil file not found (download it with product 'monthly_averaged_reanalysis'): ",
         SOIL_FILE, call. = FALSE)
  }
  src$soil <- nc_open(SOIL_FILE)
}

out <- list()
if ("Ep"   %in% VARS) out$Ep   <- make_out("Ep",   "mm month-1", "Atmospheric evaporative demand")
if ("E"    %in% VARS) out$E    <- make_out("E",    "mm month-1", "Total evaporation")
if ("Et"   %in% VARS) out$Et   <- make_out("Et",   "mm month-1", "Transpiration (ERA5-Land rotated evabs)")
if ("ED"   %in% VARS) out$ED   <- make_out("ED",   "mm month-1", "Evaporative deficit (E - Ep)")
if ("SMs"  %in% VARS) out$SMs  <- make_out("SMs",  "m3 m-3",     "Surface soil moisture (swvl1, 0-7 cm)")
if ("SMrz" %in% VARS) out$SMrz <- make_out("SMrz", "m3 m-3",     "Root-zone soil moisture (swvl1-3, 0-100 cm)")

# ------------------------------------------------------------------------------
# Month loop
# ------------------------------------------------------------------------------
message(sprintf("Processing %d months for: %s", ntime, paste(VARS, collapse = ", ")))
for (t in seq_len(ntime)) {
  Ep <- NULL; E <- NULL
  if (need_Ep) {
    v <- main_var(src$ep); raw <- read_slice(src$ep, v, t)
    Ep <- if (AED_SOURCE == "et0") raw else -raw * 1000 * ndays[t]
  }
  if (need_E) {
    raw <- read_slice(src$e, main_var(src$e), t); E <- -raw * 1000 * ndays[t]
  }
  if ("Ep" %in% VARS) put_month(out$Ep, to_gleam(Ep), t)
  if ("E"  %in% VARS) put_month(out$E,  to_gleam(E),  t)
  if ("ED" %in% VARS) put_month(out$ED, to_gleam(E - Ep), t)
  if ("Et" %in% VARS) {
    raw <- read_slice(src$et, main_var(src$et), t)
    put_month(out$Et, to_gleam(-raw * 1000 * ndays[t]), t)
  }
  if (need_soil) {
    s1 <- read_slice(src$soil, "swvl1", t)
    if ("SMs" %in% VARS) put_month(out$SMs, to_gleam(s1), t)
    if ("SMrz" %in% VARS) {
      s2 <- read_slice(src$soil, "swvl2", t); s3 <- read_slice(src$soil, "swvl3", t)
      put_month(out$SMrz, to_gleam((7 * s1 + 21 * s2 + 72 * s3) / 100), t)
    }
  }
  if (t %% 24 == 0 || t == ntime) message(sprintf("  %3d/%d (%s)", t, ntime, format(dates[t], "%Y-%m")))
}

for (o in out) nc_close(o$nc)
for (s in src) nc_close(s)

message("\nWritten:")
for (o in out) message("  ", o$path)
message("\nDone.")
