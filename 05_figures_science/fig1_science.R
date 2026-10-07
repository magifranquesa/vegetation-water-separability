#!/usr/bin/env Rscript

# ==============================================================================
# Script: 05_figures_science/fig1_science.R
#
# Project: Aridity and timescale bound the separability of vegetation water stress
#
# Figure 1: maximum vegetation-hydroclimate correlation (rho*) for five
# hydroclimatic indicators across four representative months. rho* is shown where
# the Benjamini-Hochberg FDR-adjusted p-value < 0.05 (from
# 02_correlations/04_fdr_significance_mask.R); vegetated cells without a
# significant correlation are shown in grey and non-vegetated cells in white.
# Lower panels show the distribution of rho* across vegetated cells.
#
# Inputs (READ-ONLY):
#   outputs/intermediate/absmax_spearman/abs_max_correlation_kndvi_{Ep,Et,ED,SMrz,SMs}.nc
#     (variable abs_max_correlation)
#   outputs/intermediate/absmax_fdr/fdr_padj_{Ep,Et,ED,SMrz,SMs}.tif
#   outputs/intermediate/vegetation_mask_c1.tif
#   data/external/admin_equal_earth_clean.shp
#
# Output:
#   outputs/figures_science/fig1_spatial_absmax_correlations_vegmask_grey_FDR_science.tif
#   (PDF output commented out)
#
# Run from repo root (AFTER 02_correlations/04_fdr_significance_mask.R):
#   Rscript 05_figures_science/fig1_science.R
#
# Dependencies: terra, ggplot2, ggridges, scico, patchwork, sf, grid, dplyr
# ==============================================================================

suppressPackageStartupMessages({
  library(terra); library(ggplot2); library(ggridges); library(scico)
  library(patchwork); library(sf); library(grid); library(dplyr)
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

out_file_tif <- file.path("outputs", "figures_science",
  "fig1_spatial_absmax_correlations_vegmask_grey_FDR_science.tif")
# out_file_pdf <- sub("\\.tif$", ".pdf", out_file_tif)   # PDF output disabled
for (f in c(out_file_tif))
  if (file.exists(f)) stop("REFUSING TO OVERWRITE: ", f, call. = FALSE)
dir.create(dirname(out_file_tif), recursive = TRUE, showWarnings = FALSE)

# ------------------------------------------------------------------------------
# Constants
# ------------------------------------------------------------------------------

indicators <- c("Ep", "Et", "ED", "SMrz", "SMs")
indicator_labels <- c(Ep = "AED", Et = "Et", ED = "ED", SMrz = "SMrz", SMs = "SMs")
months_index <- c(1, 4, 7, 10)
months_names <- c("January", "April", "July", "October")

font_family   <- "sans"
fig_width_mm  <- 180
fig_height_mm <- 172
fig_dpi       <- 300
base_font_size   <- 5.6
month_title_size <- 6.8
row_label_size   <- 2.45
ridge_base_size  <- 5.4
n_samp <- 5000

coastline_colour  <- "#9A9A9A"
row_label_colour  <- "#1A1A1A"
veg_bg_colour     <- "#D4D8DD"
excluded_land_col <- "#FFFFFF"
ocean_colour      <- "#FFFFFF"

colour_palette <- "bam"
colour_limits  <- c(-1, 1)
rho_colours    <- scico(50, palette = colour_palette)
extent_longlat <- ext(-170, 170, -65, 90)

ALPHA <- 0.05
fdrfile <- function(v) file.path(dir_fdr, sprintf("fdr_padj_%s.tif", v))

# ------------------------------------------------------------------------------
# Helpers (only the significance source changes: FDR p_adj < ALPHA)
# ------------------------------------------------------------------------------

sample_density_with_sig <- function(month_i, n, veg_mask) {
  bind_rows(lapply(indicators, function(v) {
    f_nc    <- file.path(dir_absmax, paste0("abs_max_correlation_kndvi_", v, ".nc"))
    cor_lyr <- rast(f_nc)[[grep("abs_max_correlation", names(rast(f_nc)))]][[month_i]]
    padj_lyr <- rast(fdrfile(v))[[month_i]]                 # FDR-adjusted p
    cor_lyr  <- crop(cor_lyr,  extent_longlat)
    padj_lyr <- crop(padj_lyr, extent_longlat)
    vm       <- crop(veg_mask, extent_longlat)
    cor_lyr[vm == 0] <- NA
    both <- c(cor_lyr, padj_lyr)
    s    <- terra::spatSample(both, size = n, method = "random", na.rm = TRUE, as.df = TRUE)
    colnames(s) <- c("rho", "padj")
    data.frame(indicator = indicator_labels[[v]], rho = s$rho, sig = s$padj < ALPHA)
  }))
}

make_ridgeline <- function(df, title = "", show_y_labels = TRUE) {
  fct_levels <- c("SMs", "SMrz", "ED", "Et", "AED")
  dens_df <- bind_rows(lapply(fct_levels, function(lbl) {
    d <- df[df$indicator == lbl, ]; d_sig <- d[d$sig, ]
    if (nrow(d) < 5) return(NULL)
    p_sig <- nrow(d_sig) / nrow(d)
    da <- density(d$rho, from = -1, to = 1, n = 512)
    if (nrow(d_sig) > 5) {
      ds <- density(d_sig$rho, from = -1, to = 1, n = 512, bw = da$bw)
      y_sig <- pmin(ds$y * p_sig, da$y)
    } else y_sig <- rep(0, 512)
    y_max <- max(da$y)
    data.frame(indicator = lbl, x = da$x, y_all = da$y / y_max, y_sig = y_sig / y_max)
  })) %>% mutate(indicator = factor(indicator, levels = fct_levels))

  ggplot(dens_df, aes(x = x, y = indicator, group = indicator)) +
    geom_ridgeline(aes(height = y_all), fill = veg_bg_colour, colour = "grey35",
                   alpha = 0.78, linewidth = 0.22, scale = 0.90) +
    geom_ridgeline_gradient(aes(height = y_sig, fill = x), colour = "grey35",
                            linewidth = 0.22, scale = 0.90) +
    scale_fill_gradientn(colours = rho_colours, limits = colour_limits, guide = "none") +
    geom_vline(xintercept = 0, linetype = "dashed", linewidth = 0.28, colour = "grey45") +
    scale_x_continuous(limits = c(-1, 1), breaks = c(-1, -0.5, 0, 0.5, 1)) +
    labs(title = title, x = NULL, y = NULL) +
    theme_minimal(base_size = ridge_base_size, base_family = font_family) +
    theme(plot.title = element_text(face = "bold", hjust = 0.5, size = ridge_base_size),
          panel.grid = element_blank(), axis.ticks.x = element_line(linewidth = 0.22),
          axis.ticks.y = element_line(color = if (show_y_labels) "black" else "transparent"),
          axis.text.x = element_text(size = ridge_base_size - 0.4, colour = "#3A3A3A"),
          axis.text.y = element_text(face = "bold", size = ridge_base_size,
                                     color = if (show_y_labels) "black" else "transparent"),
          plot.margin = margin(0, 2, 1.2, 2, unit = "mm"))
}

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
  # keep only rho* significant under FDR (p_adj < ALPHA)
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
# Read inputs
# ------------------------------------------------------------------------------

message("Generating Figure 1 (FDR): rho* maps significant under BH-FDR.")
check_file_exists(file_world_equal_earth, "Equal Earth world boundary shapefile")
check_file_exists(file_veg_mask, "vegetation mask C1 TIF")
for (v in indicators) check_file_exists(fdrfile(v), paste0("FDR mask ", v))

sf::sf_use_s2(FALSE)
world_sf <- st_read(file_world_equal_earth, quiet = TRUE)
st_crs(world_sf) <- 8857
world_sf <- st_make_valid(world_sf)
world_land <- build_land_outline(world_sf)

reference_raster <- rast(file.path(dir_absmax, paste0("abs_max_correlation_kndvi_", indicators[1], ".nc")))
veg_mask <- resample(rast(file_veg_mask), reference_raster[[1]], method = "near")

# ------------------------------------------------------------------------------
# Map panels
# ------------------------------------------------------------------------------

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

# ------------------------------------------------------------------------------
# Ridgelines
# ------------------------------------------------------------------------------

cache_dir  <- file.path("outputs", "intermediate", "fig1_cache")
cache_file <- file.path(cache_dir, sprintf("density_jan_apr_jul_oct_n%d_vegmask_FDR.rds", n_samp))
if (file.exists(cache_file)) {
  message("Loading cached FDR density samples: ", cache_file)
  dc <- readRDS(cache_file); samp_jan <- dc$jan; samp_apr <- dc$apr; samp_jul <- dc$jul; samp_oct <- dc$oct
} else {
  message("Sampling rho* (vegetated) for density plots (FDR significance)...")
  samp_jan <- sample_density_with_sig(1,  n_samp, veg_mask)
  samp_apr <- sample_density_with_sig(4,  n_samp, veg_mask)
  samp_jul <- sample_density_with_sig(7,  n_samp, veg_mask)
  samp_oct <- sample_density_with_sig(10, n_samp, veg_mask)
  dir.create(cache_dir, recursive = TRUE, showWarnings = FALSE)
  saveRDS(list(jan = samp_jan, apr = samp_apr, jul = samp_jul, oct = samp_oct), cache_file)
}

ridge_jan <- make_ridgeline(samp_jan, show_y_labels = TRUE)
ridge_apr <- make_ridgeline(samp_apr, show_y_labels = FALSE)
ridge_jul <- make_ridgeline(samp_jul, show_y_labels = FALSE)
ridge_oct <- make_ridgeline(samp_oct, show_y_labels = FALSE)

colorbar_row <- ggplot(data.frame(x = 0, y = 0, rho_star = 0), aes(x, y, fill = rho_star)) +
  geom_tile(alpha = 0) +
  scale_fill_gradientn(colours = rho_colours, limits = colour_limits,
    breaks = c(-1, -0.5, 0, 0.5, 1), name = expression("Spearman's " * rho * "*"),
    guide = guide_colourbar(direction = "horizontal", title.position = "top", title.hjust = 0.5,
                            barwidth = unit(54, "mm"), barheight = unit(2.0, "mm"))) +
  theme_void(base_family = font_family, base_size = base_font_size) +
  theme(legend.position = "bottom", legend.title = element_text(size = 6.0, colour = "black"),
        legend.text = element_text(size = 5.6, colour = "black"),
        legend.key.width = unit(3, "mm"), legend.key.height = unit(2, "mm"),
        legend.box.margin = margin(-1.2, 0, 0, 0, unit = "mm"),
        plot.margin = margin(-1.0, 1.0, -0.6, 1.0, unit = "mm")) + labs(x = NULL, y = NULL)

colorbar_centered <- (plot_spacer() | colorbar_row | plot_spacer()) + plot_layout(widths = c(3, 2, 3))

map_grid <- wrap_plots(plots_all, ncol = length(months_index), byrow = TRUE) &
  theme(plot.margin = margin(0, 0, 0, 0, unit = "mm"))
density_panel <- (ridge_jan | ridge_apr | ridge_jul | ridge_oct) / colorbar_centered +
  plot_layout(heights = c(9.1, 0.42))
density_row <- (plot_spacer() | density_panel | plot_spacer()) + plot_layout(widths = c(0.85, 18.3, 0.85))

fig <- map_grid / density_row + plot_layout(heights = c(5.15, 1.15)) &
  theme(plot.margin = margin(0.5, 0.5, 0.5, 0.5, unit = "mm"),
        plot.background = element_rect(fill = "white", colour = NA))

ggsave(out_file_tif, fig, width = fig_width_mm, height = fig_height_mm, units = "mm",
       dpi = fig_dpi, bg = "white", device = "tiff", compression = "lzw")
#ggsave(out_file_pdf, fig, width = fig_width_mm, height = fig_height_mm, units = "mm",
#       dpi = fig_dpi, bg = "white", device = cairo_pdf)
message("Figure 1 (FDR) saved to:\n  ", out_file_tif)
