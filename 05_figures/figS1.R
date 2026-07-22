#!/usr/bin/env Rscript

# ==============================================================================
# Script: 05_figures/figS1.R
#
# Project: Aridity and timescale bound the separability of vegetation water stress
#
# Supplementary Fig. S1: maps of rho* for the eight calendar months not shown in
# Fig. 1 (same projection, colours and 8x5 landscape layout as Fig. 1). rho* is
# shown where the BH-FDR-adjusted p-value < 0.05. Maps only (no ridgelines).
#
# Inputs (READ-ONLY):
#   outputs/intermediate/absmax_spearman/abs_max_correlation_kndvi_{Ep,Et,ED,SMrz,SMs}.nc
#   outputs/intermediate/absmax_fdr/fdr_padj_{Ep,Et,ED,SMrz,SMs}.tif
#   outputs/intermediate/vegetation_mask_c1.tif
#   data/external/admin_equal_earth_clean.shp
#
# Output:
#   outputs/figures/figS1_remaining_months_vegmask_grey_FDR.tif
#   outputs/figures/figS1_remaining_months_vegmask_grey_FDR.pdf
#
# Run from repo root (AFTER 02_correlations/04_fdr_significance_mask.R):
#   Rscript 05_figures/figS1.R
#
# Dependencies: terra, ggplot2, scico, patchwork, sf, grid
# ==============================================================================

suppressPackageStartupMessages({
  library(terra); library(ggplot2); library(scico); library(patchwork); library(sf); library(grid)
})

config_file <- file.path("R", "config.R")
if (file.exists(config_file)) source(config_file)

dir_absmax <- if (exists("paths") && !is.null(paths$absmax_spearman)) paths$absmax_spearman else
  file.path("outputs", "intermediate", "absmax_spearman")
dir_fdr <- file.path("outputs", "intermediate", "absmax_fdr")
file_world_equal_earth <- if (exists("paths") && !is.null(paths$world_equal_earth))
  paths$world_equal_earth else file.path("data", "external", "admin_equal_earth_clean.shp")
file_veg_mask <- if (exists("paths") && !is.null(paths$veg_mask_c1)) paths$veg_mask_c1 else
  file.path("outputs", "intermediate", "vegetation_mask_c1.tif")

out_file_tif <- file.path(if (exists("paths") && !is.null(paths$figures)) paths$figures else
  file.path("outputs", "figures"), "figS1_remaining_months_vegmask_grey_FDR.tif")
out_file_pdf <- sub("\\.tif$", ".pdf", out_file_tif)
for (f in c(out_file_tif, out_file_pdf))
  if (file.exists(f)) stop("REFUSING TO OVERWRITE: ", f, call. = FALSE)
dir.create(dirname(out_file_tif), recursive = TRUE, showWarnings = FALSE)

# ------------------------------------------------------------------------------
# Constants (8 remaining months, landscape)
# ------------------------------------------------------------------------------

indicators <- c("Ep", "Et", "ED", "SMrz", "SMs")
indicator_labels <- c(Ep = "AED", Et = "Et", ED = "ED", SMrz = "SMrz", SMs = "SMs")
months_index <- c(2, 3, 5, 6, 8, 9, 11, 12)
months_names <- c("February", "March", "May", "June", "August", "September", "November", "December")

font_family   <- "sans"
fig_width_mm  <- 340
fig_height_mm <- 150
fig_dpi       <- 600
base_font_size      <- 5.6
month_title_size    <- 6.8
row_label_size      <- 2.45
colourbar_base_size <- 5.6

coastline_colour  <- "#9A9A9A"
row_label_colour  <- "#1A1A1A"
veg_bg_colour     <- "#D4D8DD"
excluded_land_col <- "#FFFFFF"
ocean_colour      <- "#FFFFFF"

rho_colours    <- scico(50, palette = "bam")
colour_limits  <- c(-1, 1)
extent_longlat <- ext(-170, 170, -65, 90)

ALPHA <- 0.05
fdrfile <- function(v) file.path(dir_fdr, sprintf("fdr_padj_%s.tif", v))

# ------------------------------------------------------------------------------
# Helpers (significance = FDR p_adj < ALPHA)
# ------------------------------------------------------------------------------

check_file_exists <- function(fp, d = "file")
  if (!file.exists(fp)) stop("Required ", d, " not found: ", fp, call. = FALSE)

build_land_outline <- function(world_sf) {
  world_sf <- st_make_valid(world_sf)
  if (is.na(st_crs(world_sf))) st_crs(world_sf) <- 8857
  world_sf <- st_transform(world_sf, 8857)
  st_as_sf(st_sfc(suppressWarnings(st_union(st_geometry(world_sf))), crs = st_crs(world_sf)))
}

warp_to_equal_earth <- function(raster_layer, method = "bilinear") {
  if (is.na(crs(raster_layer))) crs(raster_layer) <- "EPSG:4326"
  tmp_in <- tempfile(fileext = ".tif"); tmp_out <- tempfile(fileext = ".tif")
  on.exit({ if (file.exists(tmp_in)) file.remove(tmp_in) }, add = TRUE)
  writeRaster(raster_layer, tmp_in, overwrite = TRUE)
  sf::gdal_utils(util = "warp", source = tmp_in, destination = tmp_out,
    options = c("-overwrite", "-s_srs", "EPSG:4326", "-t_srs", "EPSG:8857",
                "-r", method, "-ot", "Float32", "-dstnodata", "-9999"))
  if (!file.exists(tmp_out)) stop("GDAL warp failed.")
  rast(tmp_out)
}

make_map_panel <- function(correlation_layer, padj_layer, veg_mask, world_sf,
                           month_name, indicator_label,
                           show_column_title = FALSE, show_row_label = FALSE) {
  correlation_layer <- crop(correlation_layer, extent_longlat)
  padj_layer        <- crop(padj_layer,        extent_longlat)
  veg_mask          <- crop(veg_mask,          extent_longlat)

  correlation_layer[veg_mask == 0] <- NA
  sig_ok <- !is.na(padj_layer) & (padj_layer < ALPHA)
  correlation_layer[!sig_ok] <- NA

  veg_layer <- veg_mask; veg_layer[veg_layer == 0] <- NA

  correlation_ee <- warp_to_equal_earth(correlation_layer, method = "bilinear")
  veg_ee         <- warp_to_equal_earth(veg_layer,         method = "near")
  ext_ee         <- ext(veg_ee)

  ee_xlim_west <- sf::st_coordinates(sf::st_transform(
    sf::st_sfc(sf::st_point(c(extent_longlat[1], 60)), crs = sf::st_crs(4326)), sf::st_crs(8857)))[1L, "X"]
  ee_xlim_east <- sf::st_coordinates(sf::st_transform(
    sf::st_sfc(sf::st_point(c(extent_longlat[2], -40)), crs = sf::st_crs(4326)), sf::st_crs(8857)))[1L, "X"]
  x_rng <- ee_xlim_east - ee_xlim_west

  veg_data  <- as.data.frame(veg_ee, xy = TRUE, na.rm = TRUE); names(veg_data)  <- c("x", "y", "veg")
  plot_data <- as.data.frame(correlation_ee, xy = TRUE, na.rm = TRUE); names(plot_data) <- c("x", "y", "rho_star")

  p <- ggplot() +
    geom_sf(data = world_sf, fill = excluded_land_col, colour = NA, inherit.aes = FALSE) +
    geom_tile(data = veg_data, aes(x, y), fill = veg_bg_colour) +
    geom_tile(data = plot_data, aes(x, y, fill = rho_star)) +
    scale_fill_gradientn(colours = rho_colours, limits = colour_limits, na.value = NA, guide = "none") +
    geom_sf(data = world_sf, fill = NA, colour = coastline_colour, linewidth = 0.09, inherit.aes = FALSE) +
    coord_sf(crs = st_crs(8857), xlim = c(ee_xlim_west, ee_xlim_east),
             ylim = c(ext_ee[3], ext_ee[4]), expand = FALSE, datum = NA, clip = "off") +
    theme_void(base_family = font_family) +
    theme(plot.margin = margin(0, 0, 0, 0, unit = "mm"),
          panel.background = element_rect(fill = ocean_colour, colour = NA))

  if (show_column_title)
    p <- p + ggtitle(month_name) +
      theme(plot.title = element_text(hjust = 0.5, size = month_title_size, face = "bold",
                                      family = font_family, margin = margin(b = 0.8, unit = "mm")))
  if (show_row_label && nrow(veg_data) > 0)
    p <- p + annotate("text", x = ee_xlim_west - 0.018 * x_rng, y = mean(veg_data$y, na.rm = TRUE),
                      label = indicator_label, angle = 90, hjust = 0.5, vjust = 0.5,
                      size = row_label_size, family = font_family, fontface = "bold",
                      colour = row_label_colour) +
      theme(plot.margin = margin(0, 0, 0, 2.2, unit = "mm"))
  p
}

# ------------------------------------------------------------------------------
# Read inputs + build 8x5 map grid
# ------------------------------------------------------------------------------

message("Generating Supplementary Fig. S1 (FDR): 8 remaining months, rho* under BH-FDR.")
check_file_exists(file_world_equal_earth, "Equal Earth world boundary shapefile")
check_file_exists(file_veg_mask, "vegetation mask C1 TIF")
for (v in indicators) check_file_exists(fdrfile(v), paste0("FDR mask ", v))

sf::sf_use_s2(FALSE)
world_sf <- st_make_valid(st_read(file_world_equal_earth, quiet = TRUE))
st_crs(world_sf) <- 8857
world_land <- build_land_outline(world_sf)

reference_raster <- rast(file.path(dir_absmax, paste0("abs_max_correlation_kndvi_", indicators[1], ".nc")))
veg_mask <- resample(rast(file_veg_mask), reference_raster[[1]], method = "near")

plots_all <- list()
for (indicator in indicators) {
  message("  Reading indicator: ", indicator)
  cor_stack  <- rast(file.path(dir_absmax, paste0("abs_max_correlation_kndvi_", indicator, ".nc")))
  cor_stack  <- cor_stack[[grep("abs_max_correlation", names(cor_stack))]]
  padj_stack <- rast(fdrfile(indicator))
  for (i in seq_along(months_index)) {
    mi <- months_index[i]
    plots_all[[paste(indicator, i, sep = "_")]] <- make_map_panel(
      correlation_layer = cor_stack[[mi]], padj_layer = padj_stack[[mi]],
      veg_mask = veg_mask, world_sf = world_land,
      month_name = months_names[i], indicator_label = indicator_labels[[indicator]],
      show_column_title = indicator == indicators[1], show_row_label = i == 1)
  }
}

# Horizontal colourbar legend (matched to Fig 1)
colorbar_row <- ggplot(
  data.frame(x = 0, y = 0, rho_star = 0),
  aes(x = x, y = y, fill = rho_star)
) +
  geom_tile(alpha = 0) +
  scale_fill_gradientn(
    colours = rho_colours,
    limits  = colour_limits,
    breaks  = c(-1, -0.5, 0, 0.5, 1),
    name    = expression("Spearman's " * rho * "*"),
    guide   = guide_colourbar(
      direction = "horizontal",
      title.position = "top",
      title.hjust = 0.5,
      barwidth = unit(54, "mm"),
      barheight = unit(2.0, "mm")
    )
  ) +
  theme_void(base_family = font_family, base_size = base_font_size) +
  theme(
    legend.position = "bottom",
    legend.title = element_text(size = 6.0, colour = "black"),
    legend.text = element_text(size = colourbar_base_size, colour = "black"),
    legend.key.width = unit(3, "mm"),
    legend.key.height = unit(2, "mm"),
    legend.box.margin = margin(-4.5, 0, 0, 0, unit = "mm"),
    legend.margin = margin(0, 0, 0, 0, unit = "mm"),
    legend.box.spacing = unit(0, "mm"),
    plot.margin = margin(-4.5, 1.0, -0.6, 1.0, unit = "mm")
  ) +
  labs(x = NULL, y = NULL)

colorbar_centered <- (plot_spacer() | colorbar_row | plot_spacer()) +
  plot_layout(widths = c(3, 2, 3))

map_grid <- wrap_plots(plots_all, ncol = length(months_index), byrow = TRUE) &
  theme(plot.margin = margin(0, 0, 0, 0, unit = "mm"))

fig <- map_grid / colorbar_centered +
  plot_layout(heights = c(11, 0.32)) &
  theme(
    plot.margin = margin(0.5, 0.5, 0.5, 2.5, unit = "mm"),
    plot.background = element_rect(fill = "white", colour = NA)
  )

ggsave(out_file_tif, fig, width = fig_width_mm, height = fig_height_mm, units = "mm",
       dpi = fig_dpi, bg = "white", device = "tiff", compression = "lzw")
ggsave(out_file_pdf, fig, width = fig_width_mm, height = fig_height_mm, units = "mm",
       dpi = fig_dpi, bg = "white", device = cairo_pdf)
message("Supplementary Fig. S1 (FDR) saved to:\n  ", out_file_tif, "\n  ", out_file_pdf)
