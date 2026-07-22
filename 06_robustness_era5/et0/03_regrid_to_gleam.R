# ==============================================================================
# 06_robustness_era5/et0/03_regrid_to_gleam.R
#
# Bring the native ET0 product (ERA5-Land 0-360, node-registered, 1801x3600)
# onto the project analysis grid = GLEAM / kNDVI grid (-180/180, cell-centred,
# 1800x3600). Needed so ET0 can be crossed cell-by-cell with kNDVI in the
# correlation pipeline.
#
# Two mismatches are handled: (1) longitude convention 0-360 -> -180/180 via
# rotate(); (2) half-cell registration + latitude count (1801 vs 1800) via a
# bilinear resample onto a template read from an actual kNDVI/GLEAM file.
#
# Writes a separate file from the native 0-360 output of 02_compute_et0.R.
# Run from the repository root, after 02:
#   Rscript 06_robustness_era5/et0/03_regrid_to_gleam.R
# ==============================================================================

library(terra)

# --- Paths --------------------------------------------------------------------
native_nc   <- file.path("data", "processed", "et0_era5land",
                         "et0_era5land_monthly_1981-2022.nc")           # from 02
template_nc <- file.path("outputs", "intermediate", "kndvi", "kndvi.nc")  # GLEAM grid
out_nc      <- file.path("data", "processed", "et0_era5land",
                         "et0_era5land_monthly_1981-2022_gleamgrid.nc")

stopifnot(file.exists(native_nc), file.exists(template_nc))
if (file.exists(out_nc))
  stop("Output already exists (not overwriting): ", out_nc)

# --- Load ---------------------------------------------------------------------
r <- rast(native_nc, subds = "et0")
crs(r) <- "EPSG:4326"
message(sprintf("Native ET0: %d layers | xrange [%.1f, %.1f]",
                nlyr(r), xmin(r), xmax(r)))

tmpl <- rast(template_nc)[[1]]
crs(tmpl) <- "EPSG:4326"
message(sprintf("Template (GLEAM) grid: %d x %d | xrange [%.2f, %.2f]",
                nrow(tmpl), ncol(tmpl), xmin(tmpl), xmax(tmpl)))

# --- 0-360 -> -180/180, then resample onto the GLEAM grid ---------------------
if (xmax(r) > 180 + 1e-6) r <- rotate(r)                 # left=TRUE -> -180/180
r <- resample(r, tmpl, method = "bilinear")

# --- Write (separate file; native kept) ---------------------------------------
writeCDF(r, out_nc, varname = "et0", unit = "mm month-1",
         longname = "FAO-56 Penman-Monteith reference evapotranspiration (GLEAM grid)",
         compression = 4, missval = 1e20, overwrite = FALSE)
message("Wrote: ", out_nc)
