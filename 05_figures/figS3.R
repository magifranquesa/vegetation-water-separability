#!/usr/bin/env Rscript

# ==============================================================================
# Script: 05_figures/figS3.R
#
# Project: Aridity and timescale bound the separability of vegetation water stress
#
# Purpose:
#   SUPPLEMENTARY figure: signed shared-variance map (redundant vs compensatory
#   coupling). Decomposes the "coupled" category of the main figure into its two
#   modes. Reuses the cartographic machinery / aesthetics of the main figure
#   (Equal Earth EPSG:8857, GDAL warp, dissolved borderless land outline, 2x2
#   Jan/Apr/Jul/Oct layout, typography, dimensions, export). The signed shared
#   layer is CONTINUOUS, so a bilinear warp is added; everything else is shared
#   with the main figure.
#
#   Content:
#     - 2x2 maps of `shared` for Jan/Apr/Jul/Oct, ONLY where the full model is
#       significant (p_full_adj < 0.05). Non-signal cells use the main figure's
#       background grey (#D4D8DD).
#     - Diverging colour scale centred at 0: blue = shared > 0 (redundant,
#       positive supply–demand co-variation); red = shared < 0 (compensatory,
#       anticorrelated suppression); white ~ 0.
#     - Symmetric limit = ±p98 of |shared| over signal cells, SAME for all 4
#       panels (≈ ±0.37).
#     - Integrated month labels, panel letters a–d.
#     - Foot legend: diverging colourbar, "compensatory (shared < 0)" at the red
#       end, "redundant (shared > 0)" at the blue end; title "Shared variance
#       fraction".
#
#   *** FILE PROTECTION ***
#   - Inputs opened READ-ONLY (terra reads do not write); never modified.
#     mtimes recorded and confirmed unchanged at end.
#
# Input (READ-ONLY):
#   outputs/varpart_global/varpart_signif_global_2blocks.nc  (shared)
#   outputs/varpart_global/fdr_adjusted_pvalues.nc           (p_full_adj)
#   data/external/admin_equal_earth_clean.shp
#
# Output:
#   outputs/figures/figS_shared_signed_coupling.tif
#   outputs/figures/figS_shared_signed_coupling.pdf
#
# Run from repo root (AFTER 03_variance_partitioning/02_apply_fdr.R):
#   Rscript 05_figures/figS3.R
#
# Dependencies: terra, ggplot2, patchwork, sf, grid, scales
# ==============================================================================

suppressPackageStartupMessages({
  library(terra)
  library(ggplot2)
  library(patchwork)
  library(sf)
  library(grid)
  library(scales)
})

# ------------------------------------------------------------------------------
# Paths
# ------------------------------------------------------------------------------

config_file <- file.path("R", "config.R")
if (file.exists(config_file)) source(config_file)

file_nc <- if (exists("paths") && !is.null(paths$varpart_signif)) {
  paths$varpart_signif
} else {
  file.path("outputs", "varpart_global", "varpart_signif_global_2blocks.nc")
}
file_adj <- file.path("outputs", "varpart_global", "fdr_adjusted_pvalues.nc")
file_world_equal_earth <- file.path("data", "external", "admin_equal_earth_clean.shp")

out_file_tif <- if (exists("paths") && !is.null(paths$figures)) {
  file.path(paths$figures, "figS_shared_signed_coupling.tif")
} else {
  file.path("outputs", "figures", "figS_shared_signed_coupling.tif")
}
out_file_pdf <- sub("\\.tif$", ".pdf", out_file_tif)
dir.create(dirname(out_file_tif), recursive = TRUE, showWarnings = FALSE)

mt_nc_before  <- file.info(file_nc)$mtime
mt_adj_before <- file.info(file_adj)$mtime

# ------------------------------------------------------------------------------
# Options  (maps geometry/typography UNCHANGED from the main figure)
# ------------------------------------------------------------------------------

months_index  <- c(1, 4, 7, 10)
months_names  <- c("January", "April", "July", "October")
panel_letters <- letters[seq_along(months_index)]

ALPHA <- 0.05
FV    <- -9999

extent_longlat <- ext(-170, 174, -60, 85)
font_family    <- "Arial"

fig_width_mm  <- 180
fig_dpi       <- 600
panel_margin_mm       <- 0.03
outer_margin_top_mm   <- 0.55
outer_margin_right_mm <- 0.05
outer_margin_bottom_mm<- 0.05
outer_margin_left_mm  <- 0.05

# Maps row keeps the main figure's height; the colourbar legend needs a little
# more room than the discrete one (two-line end labels).
maps_unit_mm <- 106 / 1.032     # = main figure maps-row height
legend_mm    <- 13

# Colours
excluded_land_col  <- "#F1F1EF"   # context/analysed land without value
ocean_colour       <- "#FFFFFF"
coastline_colour   <- "#9A9A9A"
nonsig_col         <- "#D4D8DD"   # no-signal cells (main figure background grey)
panel_label_colour <- "#1A1A1A"

col_comp <- "#B2182B"   # shared < 0  (compensatory)
col_zero <- "#F7F7F7"
col_red  <- "#2166AC"   # shared > 0  (redundant)   [blue]

# ------------------------------------------------------------------------------
# Helpers
# ------------------------------------------------------------------------------

check_file_exists <- function(fp, d = "file")
  if (!file.exists(fp)) stop("Required ", d, " not found: ", fp, call. = FALSE)

set_fv_na <- function(r) { r[r == FV] <- NA; r }

warp_ee_cont <- function(r) {
  if (is.na(crs(r))) crs(r) <- "EPSG:4326"
  ti <- tempfile(fileext = ".tif"); to <- tempfile(fileext = ".tif")
  on.exit({ if (file.exists(ti)) file.remove(ti); if (file.exists(to)) file.remove(to) }, add = TRUE)
  writeRaster(r, ti, overwrite = TRUE)
  sf::gdal_utils("warp", ti, to, options = c(
    "-overwrite", "-s_srs", "EPSG:4326", "-t_srs", "EPSG:8857",
    "-r", "bilinear", "-ot", "Float32", "-dstnodata", "-9999"))
  if (!file.exists(to)) stop("GDAL warp (cont) failed.")
  x <- rast(to); x[x == -9999] <- NA; x
}
warp_ee_cat <- function(r) {
  if (is.na(crs(r))) crs(r) <- "EPSG:4326"
  ti <- tempfile(fileext = ".tif"); to <- tempfile(fileext = ".tif")
  on.exit({ if (file.exists(ti)) file.remove(ti); if (file.exists(to)) file.remove(to) }, add = TRUE)
  writeRaster(r, ti, overwrite = TRUE)
  sf::gdal_utils("warp", ti, to, options = c(
    "-overwrite", "-s_srs", "EPSG:4326", "-t_srs", "EPSG:8857",
    "-r", "near", "-ot", "Int16", "-dstnodata", "-9999"))
  if (!file.exists(to)) stop("GDAL warp (cat) failed.")
  x <- rast(to); x[x == -9999] <- NA; x
}

build_land_outline <- function(world_sf) {
  world_sf <- st_make_valid(world_sf)
  if (is.na(st_crs(world_sf))) st_crs(world_sf) <- 8857
  world_sf <- st_transform(world_sf, 8857)
  st_as_sf(st_sfc(suppressWarnings(st_union(st_geometry(world_sf))), crs = st_crs(world_sf)))
}

# Equal Earth framing (same trick as the main figure)
ee_x <- function(lon, lat)
  sf::st_coordinates(sf::st_transform(
    sf::st_sfc(sf::st_point(c(lon, lat)), crs = sf::st_crs(4326)), sf::st_crs(8857)))[1L, "X"]
ee_xlim_west <- ee_x(extent_longlat[1], 60)
ee_xlim_east <- ee_x(extent_longlat[2], -40)

# ------------------------------------------------------------------------------
# Read inputs
# ------------------------------------------------------------------------------

message("Generating SUPPLEMENTARY signed-shared figure (2x2, Equal Earth).")
check_file_exists(file_nc,  "significance varpart NetCDF (shared)")
check_file_exists(file_adj, "FDR adjusted-p NetCDF")
check_file_exists(file_world_equal_earth, "Equal Earth world boundary shapefile")

sf::sf_use_s2(FALSE)
world_sf <- st_read(file_world_equal_earth, quiet = TRUE)
st_crs(world_sf) <- 8857
world_land <- build_land_outline(world_sf)

r_sh <- set_fv_na(rast(file_nc,  subds = "shared"))
r_pf <- set_fv_na(rast(file_adj, subds = "p_full_adj"))

# Symmetric colour limit: p98 of |shared| over signal cells (4 displayed months).
sh_abs <- c()
for (m in months_index) {
  shv <- values(r_sh[[m]]); pfv <- values(r_pf[[m]])
  sig <- !is.na(pfv) & pfv < ALPHA & !is.na(shv)
  sh_abs <- c(sh_abs, abs(shv[sig]))
}
A_lim <- as.numeric(quantile(sh_abs, 0.98, na.rm = TRUE))
if (!is.finite(A_lim) || A_lim <= 0) A_lim <- max(sh_abs, na.rm = TRUE)
message(sprintf("  Symmetric colour limit |shared| p98 (signal cells): %.4f", A_lim))

# Console: exact area-weighted redundant vs compensatory over the signal cells
# actually shown here (all 12 months; uses only shared + p_full_adj).
area_vec <- as.vector(values(cellSize(r_sh[[1]], unit = "km", mask = FALSE)))
a_red <- 0; a_comp <- 0
for (m in 1:12) {
  shv <- as.vector(values(r_sh[[m]])); pfv <- as.vector(values(r_pf[[m]]))
  sig <- !is.na(shv) & !is.na(pfv) & pfv < ALPHA
  a_red  <- a_red  + sum(area_vec[sig & shv > 0], na.rm = TRUE)
  a_comp <- a_comp + sum(area_vec[sig & shv < 0], na.rm = TRUE)
}
tot_a <- a_red + a_comp
message(sprintf(paste0("  [Info] Over ALL significant cells (12 mo, area-weighted): ",
                       "redundant %.1f%% / compensatory %.1f%%"),
                a_red / tot_a * 100, a_comp / tot_a * 100))
message("  (Caption annotation uses the WITHIN-COUPLED split ~86%/14% from the paper.)")

# ------------------------------------------------------------------------------
# Panel builder (continuous signed shared)
# ------------------------------------------------------------------------------

make_signed_panel <- function(m, month_name, letter) {
  sh_l <- crop(r_sh[[m]], extent_longlat)
  pf_l <- crop(r_pf[[m]], extent_longlat)
  shv  <- values(sh_l); pfv <- values(pf_l)
  valid <- !is.na(shv)
  sig   <- valid & !is.na(pfv) & pfv < ALPHA

  sh_sig <- sh_l; v <- values(sh_sig); v[!sig] <- NA; values(sh_sig) <- v
  ns <- sh_l;     w <- rep(NA_integer_, length(shv)); w[valid & !sig] <- 1L; values(ns) <- w

  sh_ee <- warp_ee_cont(sh_sig)
  ns_ee <- warp_ee_cat(ns)
  ext_ee <- ext(sh_ee)

  df_v  <- as.data.frame(sh_ee, xy = TRUE, na.rm = TRUE); names(df_v) <- c("x","y","val")
  df_ns <- as.data.frame(ns_ee, xy = TRUE, na.rm = TRUE); names(df_ns) <- c("x","y","z")
  df_v$val <- pmax(pmin(df_v$val, A_lim), -A_lim)

  x_rng <- ee_xlim_east - ee_xlim_west
  y_rng <- ext_ee[4] - ext_ee[3]

  base_map <- ggplot() +
    geom_sf(data = world_land, fill = excluded_land_col, colour = NA, inherit.aes = FALSE) +
    geom_tile(data = df_ns, aes(x, y), fill = nonsig_col) +
    geom_tile(data = df_v, aes(x, y, fill = val)) +
    geom_sf(data = world_land, fill = NA, colour = coastline_colour,
            linewidth = 0.09, inherit.aes = FALSE) +
    scale_fill_gradient2(low = col_comp, mid = col_zero, high = col_red, midpoint = 0,
                         limits = c(-A_lim, A_lim), oob = scales::squish, guide = "none") +
    coord_sf(crs = st_crs(8857), xlim = c(ee_xlim_west, ee_xlim_east),
             ylim = c(ext_ee[3], ext_ee[4]), expand = FALSE, datum = NA) +
    labs(tag = letter) +
    theme_void(base_family = font_family) +
    theme(legend.position = "none",
          plot.tag = element_text(family = font_family, face = "bold", size = 8.0, colour = panel_label_colour),
          plot.tag.position = c(0, 1),
          plot.margin = margin(1.6, panel_margin_mm, panel_margin_mm, 1.1, unit = "mm"),
          panel.background = element_rect(fill = ocean_colour, colour = NA),
          plot.background = element_rect(fill = "white", colour = NA))

  # Integrated month label (lower-left region, as in the main figure)
  base_map +
    annotate("text",
             x = ee_xlim_west + 0.39 * x_rng,
             y = ext_ee[3] + 0.10 * y_rng,
             label = month_name, hjust = 0.5, vjust = 0,
             family = font_family, fontface = "plain",
             size = 1.9, colour = panel_label_colour)
}

# Diverging colourbar legend (foot), built once.
make_cbar_legend <- function() {
  dd <- data.frame(x = seq(-A_lim, A_lim, length.out = 200), y = 1)
  src <- ggplot(dd, aes(x, y, fill = x)) + geom_tile() +
    scale_fill_gradient2(low = col_comp, mid = col_zero, high = col_red, midpoint = 0,
      limits = c(-A_lim, A_lim),
      breaks = c(-A_lim, 0, A_lim),
      labels = c("compensatory\n(shared < 0)", "0", "redundant\n(shared > 0)"),
      name = "Shared variance fraction",
      guide = guide_colourbar(direction = "horizontal", title.position = "top",
                              title.hjust = 0.5, label.position = "bottom",
                              barwidth = unit(48, "mm"), barheight = unit(2.0, "mm"),
                              ticks.colour = "grey30")) +
    theme_void(base_family = font_family) +
    theme(legend.position = "bottom",
          legend.title = element_text(size = 5.9, face = "plain"),
          legend.text  = element_text(size = 5.6),
          legend.margin = margin(0, 0, 0, 0, unit = "mm"))
  g <- ggplotGrob(src)
  g$grobs[[which(sapply(g$grobs, function(x) x$name) == "guide-box")]]
}

# ------------------------------------------------------------------------------
# Build + assemble
# ------------------------------------------------------------------------------

message("\nBuilding panels...")
panels <- lapply(seq_along(months_index), function(i)
  make_signed_panel(months_index[i], months_names[i], panel_letters[i]))

cbar <- make_cbar_legend()
map_grid <- wrap_plots(panels, ncol = 2, byrow = TRUE)

fig_total_mm <- maps_unit_mm + legend_mm
fig <- map_grid / wrap_elements(cbar) +
  plot_layout(heights = c(1, legend_mm / maps_unit_mm)) &
  theme(plot.margin = margin(outer_margin_top_mm, outer_margin_right_mm,
                             outer_margin_bottom_mm, outer_margin_left_mm, unit = "mm"),
        plot.background = element_rect(fill = "white", colour = NA))

ggsave(out_file_tif, fig, width = fig_width_mm, height = fig_total_mm, units = "mm",
       dpi = fig_dpi, bg = "white", device = "tiff", compression = "lzw")
ggsave(out_file_pdf, fig, width = fig_width_mm, height = fig_total_mm, units = "mm",
       dpi = fig_dpi, bg = "white", device = cairo_pdf)

message("\nFigure saved to:")
message("  ", out_file_tif)
message("  ", out_file_pdf)
message(sprintf("  (canvas %.0f × %.1f mm)", fig_width_mm, fig_total_mm))

# ------------------------------------------------------------------------------
# File-protection confirmation + created files
# ------------------------------------------------------------------------------

message("\n=== File protection confirmation ===")
message(sprintf("  %s : %s", basename(file_nc),
                if (identical(mt_nc_before,  file.info(file_nc)$mtime))  "UNCHANGED (read-only respected)" else "*** CHANGED ***"))
message(sprintf("  %s : %s", basename(file_adj),
                if (identical(mt_adj_before, file.info(file_adj)$mtime)) "UNCHANGED (read-only respected)" else "*** CHANGED ***"))
message("\n  Files CREATED by this script:")
message("    ", out_file_tif)
message("    ", out_file_pdf)
message("\nDone (supplementary signed-shared map).")
