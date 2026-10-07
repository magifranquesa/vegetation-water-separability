#!/usr/bin/env Rscript

# ==============================================================================
# Script: 05_figures/fig2.R
#
# Project: Aridity and timescale bound the separability of vegetation water stress
#
# Figure 2: classification of the hydroclimatic control on interannual vegetation
# activity into five variance-partitioning categories, mapped for four
# representative months, with the seasonal cycle of control shares (e,f) and the
# distribution of full-model adjusted R2 per category (g).
#
#   The five categories separate the two significance patterns that a merged
#   "joint" class would combine:
#       (iv) coupled_pure  full_sig & !soil & !dem  -> non-separable
#       (v)  both_unique   full_sig &  soil &  dem  -> both resolvable
#   The R2 distribution (g) shows these are different populations: (v) has the
#   highest R2 (separable), while (iv) is the largest class at moderate R2.
#
#   Everything else -- Equal Earth projection, GDAL warp, land outline, 2x2 map
#   layout, mini-bars, seasonal panels (e,f), violin machinery, typography, file
#   protection.
#
# Input (READ-ONLY; mtimes confirmed unchanged):
#   outputs/varpart_global/varpart_signif_global_2blocks.nc  (shared, total_r2)
#   outputs/varpart_global/fdr_adjusted_pvalues.nc           (p_*_adj)
#   data/external/admin_equal_earth_clean.shp
#
# Output:
#   outputs/figures/fig2_dominance_split_5cat.tif
#   (PDF output commented out)
#
# Run from repo root (AFTER 03_variance_partitioning/02_apply_fdr.R):
#   Rscript 05_figures/fig2.R
#
# Dependencies: terra, ggplot2, patchwork, sf, grid, dplyr, tibble
# ==============================================================================

suppressPackageStartupMessages({
  library(terra); library(ggplot2); library(patchwork)
  library(sf); library(grid); library(dplyr)
})

config_file <- file.path("R", "config.R")
if (file.exists(config_file)) source(config_file)

file_nc <- if (exists("paths") && !is.null(paths$varpart_signif)) paths$varpart_signif else
  file.path("outputs", "varpart_global", "varpart_signif_global_2blocks.nc")
file_adj <- file.path("outputs", "varpart_global", "fdr_adjusted_pvalues.nc")
file_world_equal_earth <- file.path("data", "external", "admin_equal_earth_clean.shp")

out_file_tif <- file.path("outputs", "figures",
                          "fig2_dominance_split_5cat.tif")
# out_file_pdf <- sub("\\.tif$", ".pdf", out_file_tif)   # PDF output disabled
dir.create(dirname(out_file_tif), recursive = TRUE, showWarnings = FALSE)
for (f in c(out_file_tif))
  if (file.exists(f)) stop("REFUSING TO OVERWRITE: ", f, call. = FALSE)

mt_nc_before  <- file.info(file_nc)$mtime
mt_adj_before <- file.info(file_adj)$mtime

# ------------------------------------------------------------------------------
# Figure options
# ------------------------------------------------------------------------------

months_index  <- c(1, 4, 7, 10)
months_names  <- c("January", "April", "July", "October")
panel_letters <- LETTERS[seq_along(months_index)]

ALPHA <- 0.05
extent_longlat <- ext(-170, 174, -60, 85)
font_family <- "Arial"

fig_width_mm  <- 180
fig_height_mm <- 106
fig_dpi       <- 300

panel_margin_mm        <- 0.03
outer_margin_top_mm    <- 0.55
outer_margin_right_mm  <- 0.05
outer_margin_bottom_mm <- 0.05
outer_margin_left_mm   <- 0.05
legend_height_ratio    <- 0.032
band_height_mm         <- 34

# --- 5-level palette: coupled (JOINT, non-separable) + NEW both (resolvable) ---
cat_levels <- c("sin_senal", "soil", "demand", "coupled", "both")
cat_labels <- c(
  sin_senal = "No significant signal",
  soil      = "Soil moisture",
  demand    = "AED",
  coupled   = "Joint (non-separable)",
  both      = "Both resolvable")
cat_colours <- c(
  sin_senal = "#D4D8DD",
  soil      = "#2166AC",
  demand    = "#D98C00",
  coupled   = "#1D9E75",
  both      = "#762A83")   # NEW: muted purple for the separable minority

ocean_colour       <- "#FFFFFF"
excluded_land_col  <- "#F1F1EF"
coastline_colour   <- "#9A9A9A"
panel_label_colour <- "#1A1A1A"

NCAT <- 5L

# Band / violin cover the four significant classes now
band_cats    <- c("soil", "demand", "coupled", "both")
band_colours <- cat_colours[band_cats]
sz_axis_title <- 5.5; sz_axis_text <- 5; sz_strip <- 6
lw_line <- 0.30; lw_guide <- 0.18; pt_size <- 0.60
mon_labs  <- c("J","F","M","A","M","J","J","A","S","O","N","D")
viol_labs <- c(soil = "Soil", demand = "AED", coupled = "Joint", both = "Both")
SAMPLE_N  <- 40000L
set.seed(42L)

# ------------------------------------------------------------------------------
# Helpers
# ------------------------------------------------------------------------------

check_file_exists <- function(fp, d = "file")
  if (!file.exists(fp)) stop("Required ", d, " not found: ", fp, call. = FALSE)
classify_set_na <- function(r) { r[r == -9999] <- NA; r }

warp_to_equal_earth_cat <- function(raster_layer) {
  if (is.na(crs(raster_layer))) crs(raster_layer) <- "EPSG:4326"
  tmp_in <- tempfile(fileext = ".tif"); tmp_out <- tempfile(fileext = ".tif")
  on.exit({ if (file.exists(tmp_in)) file.remove(tmp_in); if (file.exists(tmp_out)) file.remove(tmp_out) }, add = TRUE)
  writeRaster(raster_layer, tmp_in, overwrite = TRUE)
  sf::gdal_utils(util = "warp", source = tmp_in, destination = tmp_out,
    options = c("-overwrite", "-s_srs", "EPSG:4326", "-t_srs", "EPSG:8857",
                "-r", "near", "-ot", "Int16", "-dstnodata", "-9999"))
  if (!file.exists(tmp_out)) stop("GDAL warp failed.")
  r <- rast(tmp_out); r[r == -9999] <- NA; r
}

# Split coupled_pure (4) from both_unique (5)
classify_month <- function(pf, ps, pd, sh) {
  out   <- rep(NA_integer_, length(sh)); valid <- !is.na(sh)
  full_sig <- valid & !is.na(pf) & pf < ALPHA
  soil_sig <- !is.na(ps) & ps < ALPHA
  dem_sig  <- !is.na(pd) & pd < ALPHA
  out[valid & !full_sig]                 <- 1L  # no signal
  out[full_sig &  soil_sig & !dem_sig]   <- 2L  # soil
  out[full_sig & !soil_sig &  dem_sig]   <- 3L  # demand
  out[full_sig & !soil_sig & !dem_sig]   <- 4L  # coupled_pure : JOINT, non-separable
  out[full_sig &  soil_sig &  dem_sig]   <- 5L  # both_unique  : BOTH resolvable
  out
}

build_land_outline <- function(world_sf) {
  world_sf <- st_make_valid(world_sf)
  if (is.na(st_crs(world_sf))) st_crs(world_sf) <- 8857
  world_sf <- st_transform(world_sf, 8857)
  st_as_sf(st_sfc(suppressWarnings(st_union(st_geometry(world_sf))), crs = st_crs(world_sf)))
}

category_summary <- function(vals, area = NULL) {
  vals <- as.vector(vals); keep <- !is.na(vals)
  w <- if (!is.null(area)) as.vector(area)[keep] else rep(1, sum(keep))
  vk <- vals[keep]; fac <- factor(vk, levels = 1:NCAT)
  wsum <- tapply(w, fac, sum); wsum[is.na(wsum)] <- 0
  ntab <- table(fac); total <- sum(wsum)
  pct <- if (total > 0) as.numeric(wsum) / total * 100 else rep(0, NCAT)
  tibble::tibble(
    code = 1:NCAT, cat = cat_levels, cat_f = factor(cat_levels, levels = cat_levels),
    n = as.integer(ntab), pct = pct,
    xmin = c(0, head(cumsum(pct), -1)), xmax = cumsum(pct),
    xmid = (c(0, head(cumsum(pct), -1)) + cumsum(pct)) / 2,
    pct_label = ifelse(pct >= 7, paste0(round(pct), "%"), ""))
}

print_pct <- function(vals, area, label) {
  s <- category_summary(vals, area)
  message("  [", label, "]  (area-weighted)  n_pix=", sum(s$n))
  for (k in seq_len(nrow(s)))
    message(sprintf("    %-12s : %10d px  (%5.2f%% area)", s$cat[k], s$n[k], s$pct[k]))
}

wquantile <- function(x, w, probs) {
  ok <- !is.na(x) & !is.na(w) & w > 0; x <- x[ok]; w <- w[ok]
  if (length(x) == 0) return(rep(NA_real_, length(probs)))
  o <- order(x); x <- x[o]; w <- w[o]
  cw <- (cumsum(w) - 0.5 * w) / sum(w)
  approx(cw, x, xout = probs, rule = 2, ties = "ordered")$y
}

make_mini_bar <- function(summary_df) {
  no_sig_pct <- summary_df$pct[summary_df$cat == "sin_senal"]
  sig_pct    <- sum(summary_df$pct[summary_df$cat %in% band_cats])
  raw_x <- c(no_sig_pct / 2, no_sig_pct + sig_pct / 2)
  label_df <- data.frame(
    x = pmin(pmax(raw_x, 12), 88), y = c(0.74, 0.74),
    label = c(paste0("No sig. ", round(no_sig_pct), "%"),
              paste0("Sig. ", round(sig_pct), "%")))
  ggplot(summary_df) +
    geom_rect(aes(xmin = xmin, xmax = xmax, ymin = 0.18, ymax = 0.46, fill = cat_f), colour = NA) +
    geom_text(data = label_df, aes(x = x, y = y, label = label), inherit.aes = FALSE,
              family = font_family, size = 1.55, colour = "grey10") +
    scale_fill_manual(values = cat_colours, limits = cat_levels, drop = FALSE) +
    coord_cartesian(xlim = c(0, 100), ylim = c(0, 1), expand = FALSE, clip = "off") +
    theme_void(base_family = font_family) +
    theme(legend.position = "none", plot.margin = margin(0, 0, 0, 0, unit = "mm"),
          panel.background = element_rect(fill = grDevices::adjustcolor("#FFFFFF", alpha.f = 0), colour = NA),
          plot.background  = element_rect(fill = grDevices::adjustcolor("#FFFFFF", alpha.f = 0), colour = NA))
}

make_dominance_panel <- function(cls_layer, world_land, month_name, panel_letter) {
  cls_layer <- crop(cls_layer, extent_longlat)
  area_layer <- cellSize(cls_layer, unit = "km", mask = FALSE)
  stats_df   <- category_summary(values(cls_layer), values(area_layer))

  cls_ee <- warp_to_equal_earth_cat(cls_layer); ext_ee <- ext(cls_ee)
  ee_xlim_west <- sf::st_coordinates(sf::st_transform(
    sf::st_sfc(sf::st_point(c(extent_longlat[1], 60)), crs = sf::st_crs(4326)), sf::st_crs(8857)))[1L, "X"]
  ee_xlim_east <- sf::st_coordinates(sf::st_transform(
    sf::st_sfc(sf::st_point(c(extent_longlat[2], -40)), crs = sf::st_crs(4326)), sf::st_crs(8857)))[1L, "X"]

  pd <- as.data.frame(cls_ee, xy = TRUE, na.rm = TRUE); names(pd) <- c("x", "y", "cat")
  pd <- pd[!is.na(pd$cat) & pd$cat %in% 1:NCAT, , drop = FALSE]
  pd$cat_f <- factor(cat_levels[pd$cat], levels = cat_levels)
  x_rng <- ee_xlim_east - ee_xlim_west; y_rng <- ext_ee[4] - ext_ee[3]

  base_map <- ggplot() +
    geom_sf(data = world_land, fill = excluded_land_col, colour = NA, inherit.aes = FALSE) +
    geom_tile(data = pd, aes(x, y, fill = cat_f), linewidth = 0) +
    geom_sf(data = world_land, fill = NA, colour = coastline_colour, linewidth = 0.09, inherit.aes = FALSE) +
    scale_fill_manual(values = cat_colours, labels = cat_labels, name = "Dominant control:",
                      limits = cat_levels, drop = FALSE, na.value = NA) +
    coord_sf(crs = st_crs(8857), xlim = c(ee_xlim_west, ee_xlim_east),
             ylim = c(ext_ee[3], ext_ee[4]), expand = FALSE, datum = NA) +
    labs(tag = panel_letter) +
    theme_void(base_family = font_family) +
    theme(legend.position = "none",
          plot.tag = element_text(family = font_family, face = "bold", size = 8.0, colour = panel_label_colour),
          plot.tag.position = c(0, 1),
          plot.margin = margin(1.6, panel_margin_mm, panel_margin_mm, 1.1, unit = "mm"),
          panel.background = element_rect(fill = ocean_colour, colour = NA),
          plot.background = element_rect(fill = "white", colour = NA))

  bar_grob <- ggplotGrob(make_mini_bar(stats_df))
  bar_xmin <- ee_xlim_west + 0.269  * x_rng; bar_xmax <- ee_xlim_west + 0.5155 * x_rng
  bar_ymin <- ext_ee[3] + 0.018 * y_rng;     bar_ymax <- ext_ee[3] + 0.125 * y_rng
  bar_xmid <- (bar_xmin + bar_xmax) / 2;     month_x  <- bar_xmid - 0.015 * x_rng

  base_map +
    annotate("text", x = month_x, y = bar_ymax + 0.020 * y_rng, label = month_name,
             hjust = 0.5, vjust = 0, family = font_family, fontface = "plain",
             size = 1.9, colour = panel_label_colour) +
    annotation_custom(grob = bar_grob, xmin = bar_xmin, xmax = bar_xmax, ymin = bar_ymin, ymax = bar_ymax)
}

make_legend <- function() {
  legend_df <- data.frame(x = seq_along(cat_levels), y = 1, cat_f = factor(cat_levels, levels = cat_levels))
  legend_src <- ggplot(legend_df, aes(x, y, fill = cat_f)) +
    geom_tile(width = 0.7, height = 0.7) +
    scale_fill_manual(values = cat_colours, labels = cat_labels, name = "Dominant control:",
                      limits = cat_levels, drop = FALSE) +
    theme_void(base_family = font_family, base_size = 5.8) +
    theme(legend.position = "bottom", legend.direction = "horizontal",
          legend.title = element_text(size = 5.9, face = "plain"),
          legend.text = element_text(size = 5.9, margin = margin(l = 1.2, r = 1.8, unit = "mm")),
          legend.key.width = unit(2.6, "mm"), legend.key.height = unit(2.6, "mm"),
          legend.spacing.x = unit(1.0, "mm"), legend.margin = margin(0, 0, 0, 0, unit = "mm"),
          legend.box.margin = margin(-0.8, 0, 0, 0, unit = "mm"), legend.box.spacing = unit(0, "mm"),
          plot.margin = margin(-0.8, 0, 0, 0, unit = "mm")) +
    guides(fill = guide_legend(nrow = 1, byrow = TRUE, title.position = "left",
                               title.hjust = 0, override.aes = list(colour = NA)))
  g <- ggplotGrob(legend_src)
  g$grobs[[which(sapply(g$grobs, function(x) x$name) == "guide-box")]]
}

smooth_series <- function(s) {
  s  <- s[order(s$month), ]; xo <- seq(1, 12, length.out = 240)
  yy <- spline(x = c(s$month - 12, s$month, s$month + 12), y = rep(s$pct, 3), xout = xo)$y
  data.frame(month = xo, pct = pmax(yy, 0), category = s$category[1])
}

make_band_panel <- function(d, letter, title, show_y) {
  d$category <- factor(d$category, levels = band_cats)
  dsm <- do.call(rbind, lapply(split(d, d$category), smooth_series))
  dsm$category <- factor(dsm$category, levels = band_cats)
  ggplot() +
    geom_vline(xintercept = months_index, colour = "grey85", linewidth = lw_guide, linetype = "22") +
    geom_line(data = dsm, aes(month, pct, colour = category), linewidth = lw_line) +
    geom_point(data = d, aes(month, pct, colour = category), size = pt_size) +
    scale_colour_manual(values = band_colours, limits = band_cats, name = NULL) +
    scale_x_continuous(breaks = 1:12, labels = mon_labs, expand = expansion(0.02)) +
    scale_y_continuous(breaks = seq(0, 60, 20), expand = expansion(0.02)) +
    coord_cartesian(ylim = c(0, 60)) +
    labs(tag = letter, title = title, x = NULL, y = if (show_y) "% of significant area" else NULL) +
    theme_classic(base_family = font_family) +
    theme(plot.tag = element_text(family = font_family, face = "bold", size = 8),
          plot.tag.position = c(0.004, 0.99),
          plot.title = element_text(size = sz_strip, hjust = 0.5, margin = margin(b = 1.0)),
          axis.title.y = element_text(size = sz_axis_title),
          axis.text = element_text(size = sz_axis_text, colour = "black"),
          axis.text.y = if (show_y) element_text(size = sz_axis_text, colour = "black") else element_blank(),
          axis.ticks.y = if (show_y) element_line(linewidth = lw_guide) else element_blank(),
          axis.line = element_line(linewidth = lw_guide, colour = "black"),
          axis.ticks = element_line(linewidth = lw_guide, colour = "black"),
          legend.position = "none", plot.margin = margin(1.5, 2, 0.5, 1.5, unit = "mm"))
}

make_violin_panel <- function(viol_df, stat_df, letter) {
  ggplot() +
    geom_violin(data = viol_df, aes(category, r2, fill = category),
                colour = NA, width = 0.85, scale = "area", alpha = 0.85) +
    geom_linerange(data = stat_df, aes(x = category, ymin = p5, ymax = p95),
                   linewidth = 0.3, colour = "grey20") +
    geom_linerange(data = stat_df, aes(x = category, ymin = p25, ymax = p75),
                   linewidth = 1.0, colour = "grey15") +
    geom_point(data = stat_df, aes(x = category, y = median),
               size = 1.0, shape = 21, fill = "white", colour = "grey15", stroke = 0.35) +
    scale_fill_manual(values = band_colours, guide = "none") +
    scale_x_discrete(labels = viol_labs) +
    labs(tag = letter, x = NULL, y = expression("Total adj." * R^2)) +
    theme_classic(base_family = font_family) +
    theme(plot.tag = element_text(family = font_family, face = "bold", size = 8),
          plot.tag.position = c(0.004, 0.99),
          axis.title.y = element_text(size = sz_axis_title),
          axis.text = element_text(size = sz_axis_text, colour = "black"),
          axis.line = element_line(linewidth = lw_guide, colour = "black"),
          axis.ticks = element_line(linewidth = lw_guide, colour = "black"),
          legend.position = "none", plot.margin = margin(1.5, 2, 0.5, 1.5, unit = "mm"))
}

save_outputs <- function(fig, height_mm) {
  ggsave(out_file_tif, fig, width = fig_width_mm, height = height_mm, units = "mm",
         dpi = fig_dpi, bg = "white", device = "tiff", compression = "lzw")
  # ggsave(out_file_pdf, fig, width = fig_width_mm, height = height_mm, units = "mm",
  #        dpi = fig_dpi, bg = "white", device = cairo_pdf)
}

# ------------------------------------------------------------------------------
# Read inputs
# ------------------------------------------------------------------------------

message("Generating Fig. 2: maps + seasonal band (e,f) + R2 violin (g).")
check_file_exists(file_nc,  "significance varpart NetCDF (shared/total_r2)")
check_file_exists(file_adj, "FDR adjusted-p NetCDF")
check_file_exists(file_world_equal_earth, "Equal Earth world boundary shapefile")

sf::sf_use_s2(FALSE)
world_sf <- st_read(file_world_equal_earth, quiet = TRUE)
st_crs(world_sf) <- 8857
world_land <- build_land_outline(world_sf)

r_sh  <- classify_set_na(rast(file_nc,  subds = "shared"))
r_tot <- classify_set_na(rast(file_nc,  subds = "total_r2"))
r_pf  <- classify_set_na(rast(file_adj, subds = "p_full_adj"))
r_ps  <- classify_set_na(rast(file_adj, subds = "p_soil_adj"))
r_pd  <- classify_set_na(rast(file_adj, subds = "p_demand_adj"))

# ------------------------------------------------------------------------------
# Classify + summaries (5 categories)
# ------------------------------------------------------------------------------

message("Classifying significance dominance per month (FDR, 5 categories)...")
cls_stack <- rast(r_sh); values(cls_stack) <- NA
for (m in 1:12)
  values(cls_stack[[m]]) <- classify_month(values(r_pf[[m]]), values(r_ps[[m]]),
                                           values(r_pd[[m]]), values(r_sh[[m]]))

area_full <- cellSize(cls_stack[[1]], unit = "km", mask = FALSE)
area_vec  <- as.vector(values(area_full))
lat_vec   <- as.vector(values(init(cls_stack[[1]], "y")))
hemiN     <- lat_vec > 0; hemiS <- lat_vec < 0

message("\n=== Control summary: % per category (AREA-weighted) ===")
all_vals <- unlist(lapply(1:12, function(m) values(cls_stack[[m]])))
all_area <- rep(area_vec, 12)
print_pct(all_vals, all_area, "GLOBAL all months")
for (i in seq_along(months_index))
  print_pct(values(cls_stack[[months_index[i]]]), area_vec, months_names[i])

# Seasonal band data (e/f): area-weighted, signal-only, per hemisphere x month
message("\nComputing seasonal band data (area-weighted, signal-only)...")
band_rows <- list()
for (m in 1:12) {
  code <- as.vector(values(cls_stack[[m]]))
  for (h in c("N", "S")) {
    hm  <- if (h == "N") hemiN else hemiS
    tot <- sum(area_vec[hm & code %in% 2:5], na.rm = TRUE)
    for (k in 2:5) {
      cn  <- band_cats[k - 1L]
      a_k <- sum(area_vec[hm & code == k], na.rm = TRUE)
      band_rows[[length(band_rows) + 1L]] <- data.frame(
        month = m, hemisphere = h, category = cn,
        pct = if (tot > 0) a_k / tot * 100 else NA_real_, stringsAsFactors = FALSE)
    }
  }
}
band_df <- do.call(rbind, band_rows)

# Violin data (g): pool total_r2 by category, area-weighted stats + subsample
message("Pooling total_r2 by category for the R2 violin (g)...")
pool <- setNames(lapply(band_cats, function(x) list(r2 = numeric(0), ar = numeric(0))), band_cats)
for (m in 1:12) {
  code <- as.vector(values(cls_stack[[m]])); tt <- as.vector(values(r_tot[[m]]))
  for (k in 2:5) {
    cn <- band_cats[k - 1L]; sel <- which(code == k)
    pool[[cn]]$r2 <- c(pool[[cn]]$r2, tt[sel]); pool[[cn]]$ar <- c(pool[[cn]]$ar, area_vec[sel])
  }
}

message("\n=== total_r2 by category (AREA-weighted) ===")
message(sprintf("  %-10s %12s %16s %8s %8s %8s %8s",
                "category", "n_cellmon", "area_km2", "median", "p25", "p75", "p95"))
stat_rows <- list(); viol_rows <- list()
for (cn in band_cats) {
  r2 <- pool[[cn]]$r2; ar <- pool[[cn]]$ar
  q  <- wquantile(r2, ar, c(0.05, 0.25, 0.50, 0.75, 0.95))
  message(sprintf("  %-10s %12d %16.1f %8.4f %8.4f %8.4f %8.4f",
                  cn, length(r2), sum(ar, na.rm = TRUE), q[3], q[2], q[4], q[5]))
  stat_rows[[cn]] <- data.frame(category = cn, n = length(r2), area_km2 = sum(ar, na.rm = TRUE),
                                p5 = q[1], p25 = q[2], median = q[3], p75 = q[4], p95 = q[5],
                                stringsAsFactors = FALSE)
  ok <- !is.na(r2) & !is.na(ar) & ar > 0
  r2o <- r2[ok]; aro <- ar[ok]; n_take <- min(SAMPLE_N, length(r2o))
  idx <- sample.int(length(r2o), n_take, replace = FALSE, prob = aro / sum(aro))
  viol_rows[[cn]] <- data.frame(category = cn, r2 = r2o[idx], stringsAsFactors = FALSE)
}
stat_df <- do.call(rbind, stat_rows); stat_df$category <- factor(stat_df$category, levels = band_cats)
viol_df <- do.call(rbind, viol_rows); viol_df$category <- factor(viol_df$category, levels = band_cats)

# ------------------------------------------------------------------------------
# Build panels + assemble
# ------------------------------------------------------------------------------

message("\nBuilding map panels (a-d)...")
panels <- vector("list", length(months_index))
for (i in seq_along(months_index))
  panels[[i]] <- make_dominance_panel(cls_stack[[months_index[i]]], world_land,
                                      months_names[i], panel_letters[i])

message("Building band panels (e, f, g)...")
pN <- make_band_panel(band_df[band_df$hemisphere == "N", ], "E", "Northern Hemisphere", show_y = TRUE)
pS <- make_band_panel(band_df[band_df$hemisphere == "S", ], "F", "Southern Hemisphere", show_y = FALSE)
pg <- make_violin_panel(viol_df, stat_df, "G")

band_grob     <- patchwork::patchworkGrob(pN | pS | pg)
shared_legend <- make_legend()
map_grid      <- wrap_plots(panels, ncol = 2, byrow = TRUE)

maps_unit_mm <- fig_height_mm / (1 + legend_height_ratio)
band_ratio   <- band_height_mm / maps_unit_mm
fig_total_mm <- maps_unit_mm * (1 + band_ratio + legend_height_ratio)

fig <- map_grid / wrap_elements(full = band_grob) / wrap_elements(shared_legend) +
  plot_layout(heights = c(1, band_ratio, legend_height_ratio)) &
  theme(plot.margin = margin(outer_margin_top_mm, outer_margin_right_mm,
                             outer_margin_bottom_mm, outer_margin_left_mm, unit = "mm"),
        plot.background = element_rect(fill = "white", colour = NA))

save_outputs(fig, fig_total_mm)

message("\nFigure saved to:\n  ", out_file_tif)
message(sprintf("  (canvas %.0f x %.1f mm; band %.0f mm)", fig_width_mm, fig_total_mm, band_height_mm))

message("\n=== File protection confirmation ===")
message(sprintf("  %s : %s", basename(file_nc),
                if (identical(mt_nc_before,  file.info(file_nc)$mtime))  "UNCHANGED" else "*** CHANGED ***"))
message(sprintf("  %s : %s", basename(file_adj),
                if (identical(mt_adj_before, file.info(file_adj)$mtime)) "UNCHANGED" else "*** CHANGED ***"))
message("\nDone: Fig. 2 (five-category dominance split).")
