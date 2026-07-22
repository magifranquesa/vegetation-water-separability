#!/usr/bin/env Rscript

# Build and save the C1 vegetated domain mask.
# C1: Et > 0 in at least one month (1981-2022) AND mean kNDVI > 0.025
#     (approx. NDVI > 0.16; ~108.6 M km²)
#
# Spatial extent: -180 to 180 longitude, -65 to 90 latitude (full data domain).
#
# Output:
#   outputs/intermediate/vegetation_mask_c1.tif  <- binary mask (INT1U) for reuse

suppressPackageStartupMessages(library(terra))

config_file <- file.path("R", "config.R")
if (file.exists(config_file)) source(config_file)

et_file <- if (exists("paths") && !is.null(paths$et_gleam)) paths$et_gleam else
  file.path("outputs", "intermediate", "concatenated", "Et_GLEAM_v4.2a_MO_1981-2022_NCO.nc")

kndvi_file <- if (exists("paths") && !is.null(paths$kndvi)) paths$kndvi else
  file.path("outputs", "intermediate", "kndvi", "kndvi.nc")

mask_out <- if (exists("paths") && !is.null(paths$veg_mask_c1)) paths$veg_mask_c1 else
  file.path("outputs", "intermediate", "vegetation_mask_c1.tif")

dir.create(dirname(mask_out), recursive = TRUE, showWarnings = FALSE)

e_ll <- ext(-180, 180, -65, 90)

message("Building C1 mask (Et > 0 any month AND mean kNDVI > 0.025)...")

et_nc   <- rast(et_file)
et_crop <- crop(et_nc, e_ll)
m1      <- app(et_crop, function(x) as.integer(any(x > 0, na.rm = TRUE)))

kndvi_nc   <- rast(kndvi_file)
kndvi_crop <- crop(kndvi_nc, e_ll)
k1         <- app(kndvi_crop, function(x) {
  mn <- mean(x, na.rm = TRUE)
  as.integer(!is.na(mn) && mn > 0.025)
})

c1 <- m1 & k1
crs(c1) <- "EPSG:4326"

cell_area  <- cellSize(c1, unit = "km")
domain_km2 <- global(cell_area * c1, "sum", na.rm = TRUE)$sum
message("Vegetated domain area (C1): ", round(domain_km2 / 1e6, 2), " million km²")

writeRaster(c1, mask_out, datatype = "INT1U", overwrite = TRUE)
message("Mask saved: ", mask_out)
