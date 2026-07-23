#!/usr/bin/env Rscript
# -*- coding: utf-8 -*-

# ==============================================================================
# Script: 05_figures/figS4.R
#
# Project: Aridity and timescale bound the separability of vegetation water stress
#
# Purpose:
#   Two-panel supplementary map, arranged side by side (a | b), in Equal Earth
#   projection:
#
#     a) vegetated land mask used in the study
#     b) Zomer aridity-index classes
#
#   The aridity raster is expected to be ALREADY CLASSIFIED with integer values:
#     1 = Hyper arid
#     2 = Arid
#     3 = Semi-arid
#     4 = Dry sub-humid
#     5 = Humid
#
# Inputs (READ-ONLY):
#   ./data/processed/aridity/AI_1982_2022_Zomer_classes.tif
#   ./data/external/admin_equal_earth_clean.shp
#   ./outputs/intermediate/vegetation_mask_c1.tif
#
# Outputs:
#   ./outputs/figures/figS4_vegetated_land_aridity.tif
#   ./outputs/figures/figS4_vegetated_land_aridity.pdf
#
# Run from repo root or anywhere:
#   Rscript 05_figures/figS4.R
#
# Dependencies: terra, ggplot2, patchwork, sf, grid
# ==============================================================================

suppressPackageStartupMessages({
  library(terra)
  library(ggplot2)
  library(patchwork)
  library(sf)
  library(grid)
})

# ------------------------------------------------------------------------------
# Paths
# ------------------------------------------------------------------------------

config_file <- file.path("R", "config.R")
if (file.exists(config_file)) source(config_file)

project_root <- "."

file_ai_classes <- file.path(project_root, "data", "processed", "aridity",
                             "AI_1982_2022_Zomer_classes.tif")

file_world_equal_earth <- if (exists("paths") && !is.null(paths$world_equal_earth)) {
  paths$world_equal_earth
} else {
  file.path(project_root, "data", "external", "admin_equal_earth_clean.shp")
}

file_veg_mask <- if (exists("paths") && !is.null(paths$veg_mask_c1)) {
  paths$veg_mask_c1
} else {
  file.path(project_root, "outputs", "intermediate", "vegetation_mask_c1.tif")
}

out_file_tif <- if (exists("paths") && !is.null(paths$figures)) {
  file.path(paths$figures, "figS4_vegetated_land_aridity.tif")
} else {
  file.path(project_root, "outputs", "figures", "figS4_vegetated_land_aridity.tif")
}
out_file_pdf <- sub("\\.tif$", ".pdf", out_file_tif)
dir.create(dirname(out_file_tif), recursive = TRUE, showWarnings = FALSE)

mt_ai_before    <- file.info(file_ai_classes)$mtime
mt_world_before <- file.info(file_world_equal_earth)$mtime
mt_veg_before   <- file.info(file_veg_mask)$mtime

# ------------------------------------------------------------------------------
# Figure options
# ------------------------------------------------------------------------------

font_family <- "Arial"

fig_width_mm  <- 180
fig_height_mm <- 50    # single row of two world maps -> short canvas (tune this)
fig_dpi       <- 600

panel_label_size_pt <- 8.0
axis_text_size_pt   <- 5.8
legend_title_pt     <- 5.8    # shrunk for half-column maps
legend_text_pt      <- 5.2    # shrunk for half-column maps (>= 5 pt Nature min)
annot_title_mm      <- 1.80   # geom_text units; approx. 5.1 pt (Nature >= 5 pt)
annot_text_mm       <- 1.78   # geom_text units; approx. 5.1 pt (Nature >= 5 pt)

line_grid_mm  <- 0.12         # approx. 0.34 pt
line_coast_mm <- 0.09         # approx. 0.25 pt

# Manual graticule labels. Using manual labels avoids the misplaced/left-aligned
# axis labels sometimes produced by coord_sf(label_graticule = ...) in Equal Earth.
coord_label_size_mm <- 1.85   # geom_text units; approx. 5.3 pt
# Half-column maps: 30° longitude labels collide, so use 60° spacing here
# (this also thins the meridian graticule, which was too dense at this size).
lon_breaks <- seq(-180, 180, by = 60)
lat_breaks <- seq(-60,  60, by = 30)

# The ArcGIS supplementary overview used the full longitude frame. If you want
# to match the tighter framing of the main dominance maps, set this to TRUE.
USE_PROJECT_COMPACT_FRAME <- FALSE
extent_longlat <- if (USE_PROJECT_COMPACT_FRAME) {
  ext(-170, 174, -60, 85)
} else {
  ext(-180, 180, -60, 85)
}

# Keep the period label used in the existing ArcGIS panel annotation. Change to
# "1982–2022" only if the manuscript/caption should be strictly raster-period based.
veg_period_label <- "1981–2022"

# Colours aligned with the ArcGIS version and the project map style.
vegetation_colour <- "#64A900"
ocean_colour      <- "#FFFFFF"
land_base_colour  <- "#F1F1EF"
coastline_colour  <- "#9A9A9A"
grid_colour       <- "#BDBDBD"
text_colour       <- "#2A2A2A"
axis_text_colour  <- "#6E6E6E"
legend_bg_colour  <- grDevices::adjustcolor("#FFFFFF", alpha.f = 0.86)
annot_bg_colour   <- grDevices::adjustcolor("#FFFFFF", alpha.f = 0.82)

ai_levels <- c("Hyper Arid", "Arid", "Semi-Arid", "Dry sub-humid", "Humid")
ai_values <- setNames(1:5, ai_levels)
ai_colours <- c(
  "Hyper Arid"    = "#9C5B04",
  "Arid"          = "#D9B15F",
  "Semi-Arid"     = "#F3E5B8",
  "Dry sub-humid" = "#CDEDEA",
  "Humid"         = "#006D63"
)

# ------------------------------------------------------------------------------
# Helpers
# ------------------------------------------------------------------------------

check_file_exists <- function(file_path, description = "file") {
  if (!file.exists(file_path)) {
    stop("Required ", description, " not found: ", file_path, call. = FALSE)
  }
}

set_longlat_crs_if_missing <- function(r) {
  if (is.na(crs(r))) crs(r) <- "EPSG:4326"
  r
}

build_land_outline <- function(world_sf) {
  world_sf <- st_make_valid(world_sf)
  if (is.na(st_crs(world_sf))) st_crs(world_sf) <- 8857
  world_sf <- st_transform(world_sf, 8857)
  st_as_sf(st_sfc(suppressWarnings(st_union(st_geometry(world_sf))), crs = st_crs(world_sf)))
}

warp_to_equal_earth_cat <- function(raster_layer,
                                    write_datatype = "INT2S",
                                    gdal_datatype  = "Int16",
                                    nodata_value   = -9999) {
  raster_layer <- set_longlat_crs_if_missing(raster_layer)
  tmp_in  <- tempfile(fileext = ".tif")
  tmp_out <- tempfile(fileext = ".tif")
  on.exit({
    if (file.exists(tmp_in))  file.remove(tmp_in)
    if (file.exists(tmp_out)) file.remove(tmp_out)
  }, add = TRUE)

  # terra::writeRaster() does not use GDAL's "Byte" name for datatype.
  # INT2S keeps the categorical values unchanged and can safely store -9999
  # as nodata, avoiding both the Byte warning and GDAL nodata clamping.
  writeRaster(
    raster_layer,
    tmp_in,
    overwrite = TRUE,
    datatype = write_datatype,
    NAflag = nodata_value
  )

  sf::gdal_utils(
    util = "warp",
    source = tmp_in,
    destination = tmp_out,
    options = c(
      "-overwrite",
      "-s_srs", "EPSG:4326",
      "-t_srs", "EPSG:8857",
      "-r", "near",
      "-ot", gdal_datatype,
      "-srcnodata", as.character(nodata_value),
      "-dstnodata", as.character(nodata_value)
    )
  )
  if (!file.exists(tmp_out)) {
    stop("GDAL warp failed: temporary projected raster was not created.", call. = FALSE)
  }
  r <- rast(tmp_out)
  r[r == nodata_value] <- NA
  r
}

ee_x <- function(lon, lat) {
  sf::st_coordinates(sf::st_transform(
    sf::st_sfc(sf::st_point(c(lon, lat)), crs = sf::st_crs(4326)),
    sf::st_crs(8857)
  ))[1L, "X"]
}

format_lon <- function(lon) {
  if (abs(lon) == 180) return("180°")
  if (lon == 0) return("0°")
  paste0(abs(lon), "°", ifelse(lon < 0, "W", "E"))
}

format_lat <- function(lat) {
  if (lat == 0) return("0°")
  paste0(abs(lat), "°", ifelse(lat < 0, "S", "N"))
}

project_xy <- function(lon, lat) {
  xy <- sf::st_coordinates(sf::st_transform(
    sf::st_sfc(sf::st_point(c(lon, lat)), crs = sf::st_crs(4326)),
    sf::st_crs(8857)
  ))
  c(x = unname(xy[1L, "X"]), y = unname(xy[1L, "Y"]))
}

make_graticule_lines <- function() {
  lat_seq <- seq(extent_longlat[3], extent_longlat[4], length.out = 260)
  lon_seq <- seq(extent_longlat[1], extent_longlat[2], length.out = 520)

  meridians <- lapply(lon_breaks, function(lon) {
    sf::st_linestring(cbind(rep(lon, length(lat_seq)), lat_seq))
  })
  parallels <- lapply(lat_breaks, function(lat) {
    sf::st_linestring(cbind(lon_seq, rep(lat, length(lon_seq))))
  })

  gr <- sf::st_sf(
    type = c(rep("meridian", length(meridians)), rep("parallel", length(parallels))),
    geometry = sf::st_sfc(c(meridians, parallels), crs = sf::st_crs(4326))
  )
  sf::st_transform(gr, 8857)
}

make_coord_labels <- function() {
  x_rng <- diff(ee_xlim)
  y_rng <- diff(ee_ylim)
  nudge_x <- 0.010 * x_rng     # outward offset for latitude labels (sides)
  nudge_y <- 0.016 * y_rng     # outward offset for longitude labels (top/bottom)

  # Equal Earth is an oval, not a rectangle. Each label is anchored to where its
  # graticule line meets the projected boundary, so it sits ON its grid line:
  #   longitude -> meridian crossing the top / bottom edge (follows the curve),
  #   latitude  -> parallel endpoint on the left / right edge.
  # All are placed just OUTSIDE that point. Meridians converge poleward, so the
  # top longitude row is naturally narrower than the bottom — that is the grid.
  lon_top <- do.call(rbind, lapply(lon_breaks, function(lon) {
    xy <- project_xy(lon, extent_longlat[4])
    data.frame(x = xy[["x"]], y = xy[["y"]] + nudge_y, label = format_lon(lon),
               hjust = 0.5, vjust = 0, stringsAsFactors = FALSE)
  }))
  lon_bottom <- do.call(rbind, lapply(lon_breaks, function(lon) {
    xy <- project_xy(lon, extent_longlat[3])
    data.frame(x = xy[["x"]], y = xy[["y"]] - nudge_y, label = format_lon(lon),
               hjust = 0.5, vjust = 1, stringsAsFactors = FALSE)
  }))

  lat_left <- do.call(rbind, lapply(lat_breaks, function(lat) {
    xy <- project_xy(extent_longlat[1], lat)
    data.frame(
      x = xy[["x"]] - nudge_x, y = xy[["y"]],
      label = format_lat(lat), hjust = 1, vjust = 0.5,
      stringsAsFactors = FALSE
    )
  }))

  lat_right <- do.call(rbind, lapply(lat_breaks, function(lat) {
    xy <- project_xy(extent_longlat[2], lat)
    data.frame(
      x = xy[["x"]] + nudge_x, y = xy[["y"]],
      label = format_lat(lat), hjust = 0, vjust = 0.5,
      stringsAsFactors = FALSE
    )
  }))

  d <- rbind(lon_top, lon_bottom, lat_left, lat_right)
  d <- d[is.finite(d$x) & is.finite(d$y), , drop = FALSE]
  rownames(d) <- NULL
  d
}

prepare_vegetation_mask <- function(r) {
  r <- crop(set_longlat_crs_if_missing(r), extent_longlat)
  # vegetation_mask_c1 is expected to be 1 for vegetated land and 0/NA otherwise.
  r[r <= 0] <- NA
  r[r > 0]  <- 1
  r
}

prepare_ai_class_raster <- function(r) {
  r <- crop(set_longlat_crs_if_missing(r), extent_longlat)
  r[!(r %in% ai_values)] <- NA
  r
}

raster_to_df <- function(r, value_name) {
  d <- as.data.frame(r, xy = TRUE, na.rm = TRUE)
  names(d) <- c("x", "y", value_name)
  d
}

make_base_map <- function(letter) {
  suppressMessages({
    ggplot() +
    coord_sf(
      crs = st_crs(8857),
      xlim = ee_xlim,
      ylim = ee_ylim,
      datum = NA,
      expand = FALSE,
      clip = "off"
    ) +
    geom_sf(data = world_land, fill = land_base_colour, colour = NA, inherit.aes = FALSE) +
    labs(tag = letter) +
    theme_void(base_family = font_family) +
    theme(
      plot.tag = element_text(
        family = font_family,
        face = "bold",
        size = panel_label_size_pt,
        colour = text_colour
      ),
      plot.tag.position = c(0.045, 0.92),
      panel.background = element_rect(fill = ocean_colour, colour = NA),
      plot.background  = element_rect(fill = "white", colour = NA),
      axis.title = element_blank(),
      axis.text = element_blank(),
      axis.ticks = element_blank(),
      # Room OUTSIDE the oval boundary for the coordinate labels (sides:
      # latitude; top & bottom: longitude). clip = "off" lets them render here.
      plot.margin = margin(4.0, 8.0, 4.0, 8.0, unit = "mm")
    )
  })
}

add_map_frame <- function(p) {
  suppressMessages({
    p +
    geom_sf(data = graticule_lines, fill = NA, colour = grid_colour,
            linewidth = line_grid_mm, inherit.aes = FALSE) +
    geom_sf(data = world_land, fill = NA, colour = coastline_colour,
            linewidth = line_coast_mm, inherit.aes = FALSE)
  })
}

add_coord_labels <- function(p) {
  p +
    geom_text(
      data = coord_labels,
      aes(x = x, y = y, label = label, hjust = hjust, vjust = vjust),
      inherit.aes = FALSE,
      family = font_family,
      size = coord_label_size_mm,
      colour = axis_text_colour
    )
}

make_vegetation_panel <- function(veg_df, veg_area_mkm2) {
  x_rng <- diff(ee_xlim)
  y_rng <- diff(ee_ylim)

  # Annotation placed low over the Southern Ocean (white) so it covers no land;
  # no background box is needed there.
  annot_x  <- mean(ee_xlim) + 0.105 * x_rng
  annot_y1 <- ee_ylim[1] + 0.115 * y_rng
  annot_y2 <- ee_ylim[1] + 0.053 * y_rng

  p <- make_base_map("a") +
    geom_tile(data = veg_df, aes(x, y), fill = vegetation_colour, linewidth = 0)

  p <- add_map_frame(p) +
    annotate(
      "text",
      x = annot_x,
      y = annot_y1,
      label = paste0("Vegetated land (", veg_period_label, "; ~",
                     round(veg_area_mkm2), " × 10⁶ km²)"),
      family = font_family,
      fontface = "bold",
      size = annot_title_mm,
      colour = text_colour,
      hjust = 0.5
    ) +
    annotate(
      "text",
      x = annot_x,
      y = annot_y2,
      label = "Et > 0 in ≥1 month and mean kNDVI > 0.025",
      family = font_family,
      fontface = "plain",
      size = annot_text_mm,
      colour = text_colour,
      hjust = 0.5
    )

  add_coord_labels(p)
}

make_aridity_panel <- function(ai_df) {
  p <- make_base_map("b") +
    geom_tile(data = ai_df, aes(x, y, fill = aridity_class), linewidth = 0) +
    scale_fill_manual(
      values = ai_colours,
      limits = ai_levels,
      drop = FALSE,
      name = "Aridity Index"
    )

  p <- add_map_frame(p) +
    guides(fill = guide_legend(
      ncol = 1,
      byrow = TRUE,
      title.position = "top",
      override.aes = list(colour = NA)
    )) +
    theme(
      legend.position = c(0.115, 0.350),
      legend.justification = c(0, 0.5),
      legend.background = element_blank(),
      legend.key = element_blank(),
      legend.title = element_text(family = font_family, size = legend_title_pt,
                                  face = "bold", colour = text_colour),
      legend.text = element_text(family = font_family, size = legend_text_pt,
                                 colour = text_colour),
      legend.key.width = unit(2.6, "mm"),
      legend.key.height = unit(2.6, "mm"),
      legend.spacing.y = unit(0.1, "mm"),
      legend.margin = margin(0.8, 1.0, 0.8, 1.0, unit = "mm")
    )

  add_coord_labels(p)
}

save_outputs <- function(fig) {
  ggsave(out_file_tif, fig, width = fig_width_mm, height = fig_height_mm,
         units = "mm", dpi = fig_dpi, bg = "white",
         device = "tiff", compression = "lzw")
  ggsave(out_file_pdf, fig, width = fig_width_mm, height = fig_height_mm,
         units = "mm", dpi = fig_dpi, bg = "white", device = cairo_pdf)
}

# ------------------------------------------------------------------------------
# Read inputs
# ------------------------------------------------------------------------------

message("Generating supplementary figure (side-by-side): vegetated land + Zomer aridity classes.")
check_file_exists(file_ai_classes, "Zomer aridity-class raster")
check_file_exists(file_world_equal_earth, "Equal Earth world boundary shapefile")
check_file_exists(file_veg_mask, "vegetation mask raster")

sf::sf_use_s2(FALSE)
world_sf <- st_read(file_world_equal_earth, quiet = TRUE)
st_crs(world_sf) <- 8857
world_land <- build_land_outline(world_sf)

r_ai  <- rast(file_ai_classes)
r_veg <- rast(file_veg_mask)

# ------------------------------------------------------------------------------
# Prepare raster layers
# ------------------------------------------------------------------------------

message("Preparing vegetation mask...")
veg_ll <- prepare_vegetation_mask(r_veg)
veg_area <- cellSize(veg_ll, unit = "km", mask = FALSE)
veg_area_vec <- as.vector(values(veg_area))
veg_vec <- as.vector(values(veg_ll))
veg_area_mkm2 <- sum(veg_area_vec[!is.na(veg_vec)], na.rm = TRUE) / 1e6
message(sprintf("  Vegetated area: %.2f × 10^6 km²", veg_area_mkm2))

message("Preparing aridity classes (expected values 1--5)...")
ai_ll <- prepare_ai_class_raster(r_ai)
ai_vals <- sort(unique(as.vector(values(ai_ll))))
ai_vals <- ai_vals[!is.na(ai_vals)]
message("  Aridity class values found after masking/crop: ", paste(ai_vals, collapse = ", "))
missing_vals <- setdiff(ai_values, ai_vals)
if (length(missing_vals) > 0) {
  message("  Note: class value(s) absent in cropped raster: ", paste(missing_vals, collapse = ", "))
}

message("Warping rasters to Equal Earth (EPSG:8857, nearest-neighbour)...")
veg_ee <- warp_to_equal_earth_cat(veg_ll)
ai_ee  <- warp_to_equal_earth_cat(ai_ll)

if (USE_PROJECT_COMPACT_FRAME) {
  ee_xlim <- c(ee_x(extent_longlat[1], 60), ee_x(extent_longlat[2], -40))
} else {
  ee_xlim <- c(ext(veg_ee)[1], ext(veg_ee)[2])
}
ee_ylim <- c(ext(veg_ee)[3], ext(veg_ee)[4])

graticule_lines <- make_graticule_lines()
coord_labels    <- make_coord_labels()

veg_df <- raster_to_df(veg_ee, "veg")
veg_df <- veg_df[!is.na(veg_df$veg), , drop = FALSE]

ai_df <- raster_to_df(ai_ee, "class_value")
ai_df <- ai_df[!is.na(ai_df$class_value), , drop = FALSE]
ai_df$aridity_class <- names(ai_values)[match(ai_df$class_value, ai_values)]
ai_df <- ai_df[!is.na(ai_df$aridity_class), , drop = FALSE]
ai_df$aridity_class <- factor(ai_df$aridity_class, levels = ai_levels)

if (nrow(veg_df) == 0) stop("Vegetation mask has no non-NA cells after crop/warp.", call. = FALSE)
if (nrow(ai_df) == 0) stop("No aridity class values 1--5 found after crop/warp.", call. = FALSE)

# ------------------------------------------------------------------------------
# Build and save figure
# ------------------------------------------------------------------------------

message("Building panels...")
p_veg <- make_vegetation_panel(veg_df, veg_area_mkm2)
p_ai  <- make_aridity_panel(ai_df)

# ONLY change vs the stacked script: one row, two columns (a | b).
# Parenthesised because in R `+` binds tighter than `|`.
fig <- (p_veg | p_ai) +
  plot_layout(widths = c(1, 1)) &
  theme(
    plot.background = element_rect(fill = "white", colour = NA),
    plot.margin = margin(0.5, 0.5, 0.5, 0.5, unit = "mm")
  )

save_outputs(fig)

message("\nFigure saved to:")
message("  ", out_file_tif)
message("  ", out_file_pdf)
message(sprintf("  (canvas %.0f × %.0f mm; %d dpi TIFF)", fig_width_mm, fig_height_mm, fig_dpi))

# ------------------------------------------------------------------------------
# File-protection confirmation + created files
# ------------------------------------------------------------------------------

message("\n=== File protection confirmation ===")
message(sprintf("  %s : %s", basename(file_ai_classes),
                if (identical(mt_ai_before, file.info(file_ai_classes)$mtime))
                  "UNCHANGED (read-only respected)" else "*** CHANGED ***"))
message(sprintf("  %s : %s", basename(file_world_equal_earth),
                if (identical(mt_world_before, file.info(file_world_equal_earth)$mtime))
                  "UNCHANGED (read-only respected)" else "*** CHANGED ***"))
message(sprintf("  %s : %s", basename(file_veg_mask),
                if (identical(mt_veg_before, file.info(file_veg_mask)$mtime))
                  "UNCHANGED (read-only respected)" else "*** CHANGED ***"))
message("\n  Files CREATED by this script:")
message("    ", out_file_tif)
message("    ", out_file_pdf)
message("\nDone.")
