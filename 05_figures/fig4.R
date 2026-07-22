#!/usr/bin/env Rscript

# ==============================================================================
# Script: 05_figures/fig4.R
#
# Project: Aridity and timescale bound the separability of vegetation water stress
# Design: Nature-style, coherent with Fig. 1 and dominance maps
#
# Purpose:
#   Spatial agreement between rho_max maps across indicators.
#   Panel a: July pairwise scatterplots (hex-bin density).
#   Panel b: monthly Spearman heatmap.
#
# Design principles in this version:
#   - Same final figure logic as the redesigned dominance maps: 180 mm width,
#     compact margins, Arial/Helvetica-like typography, lower-case panel tags.
#   - Heatmap uses the same signed-rho colour scale as Fig. 1: scico palette
#     "bam", with limits [-1, 1]. This keeps a single visual grammar for
#     signed correlation values across all panels.
#   - Scatterplot density is intentionally neutral/blue rather than using the
#     rho colour scale, because hex-bin colour represents pixel count, not rho.
#   - Red regression lines have been removed; fitted trends are shown with a
#     quiet dark grey line, and reference lines with a light dashed line.
#   - In-panel rho annotations use semi-transparent white labels without borders,
#     matching the visual language used in the dominance maps.
#
# Input:
#   outputs/intermediate/max_spearman/max_correlation_kndvi_<v>.nc
#   outputs/tables/spatial_correlations_rho_max.csv
#   outputs/intermediate/vegetation_mask_c1.tif
#
# Output:
#   outputs/figures/fig2_spatial_agreement_complex_nature.tif
#   outputs/figures/fig2_spatial_agreement_complex_nature.pdf
#   outputs/figures/fig2_spatial_agreement_complex_nature.png
#
# Notes:
#   - Run summarise_spatial_correlations_rho_max.R first to generate the CSV.
#
# Dependencies:
#   terra, ggplot2, hexbin, dplyr, tidyr, patchwork, readr, scico, grid
# ==============================================================================
suppressPackageStartupMessages({
  library(terra)
  library(ggplot2)
  library(hexbin)
  library(dplyr)
  library(tidyr)
  library(patchwork)
  library(readr)
  library(scico)
  library(grid)
})

# ------------------------------------------------------------------------------
# Paths
# ------------------------------------------------------------------------------

config_file <- file.path("R", "config.R")
if (file.exists(config_file)) source(config_file)

dir_max <- if (exists("paths") && !is.null(paths$max_spearman)) {
  paths$max_spearman
} else {
  file.path("outputs", "intermediate", "max_spearman")
}

heatmap_csv <- if (exists("paths") && !is.null(paths$tables)) {
  file.path(paths$tables, "spatial_correlations_rho_max.csv")
} else {
  file.path("outputs", "tables", "spatial_correlations_rho_max.csv")
}

file_veg_mask <- if (exists("paths") && !is.null(paths$veg_mask_c1)) {
  paths$veg_mask_c1
} else {
  file.path("outputs", "intermediate", "vegetation_mask_c1.tif")
}

out_file_tif <- if (exists("paths") && !is.null(paths$figures)) {
  file.path(paths$figures, "fig2_spatial_agreement_complex_nature.tif")
} else {
  file.path("outputs", "figures", "fig2_spatial_agreement_complex_nature.tif")
}

out_file_pdf <- sub("\\.tif$", ".pdf", out_file_tif)
out_file_png <- sub("\\.tif$", ".png", out_file_tif)
dir.create(dirname(out_file_tif), recursive = TRUE, showWarnings = FALSE)

# ------------------------------------------------------------------------------
# Figure options
# ------------------------------------------------------------------------------

n_samp    <- 50000
month_idx <- 7
e_ll      <- ext(-170, 170, -65, 90)

variables  <- c("Ep", "Et", "ED", "SMrz", "SMs")
var_labels <- c(Ep = "AED", Et = "Et", ED = "ED", SMrz = "SMrz", SMs = "SMs")

month_abbrev <- c("Jan", "Feb", "Mar", "Apr", "May", "Jun",
                  "Jul", "Aug", "Sep", "Oct", "Nov", "Dec")
month_full   <- c("January", "February", "March", "April", "May", "June",
                  "July", "August", "September", "October", "November", "December")
month_lookup <- setNames(month_abbrev, month_full)

# Explicit top-to-bottom heatmap row order. Keeping this as a named vector makes
# the editorial ordering transparent and avoids hidden reordering logic inside
# the data-processing pipeline.
pair_levels_topdown <- c(
  "AED vs ED",
  "AED vs SMrz",
  "AED vs SMs",
  "AED vs Et",
  "Et vs ED",
  "Et vs SMrz",
  "Et vs SMs",
  "ED vs SMrz",
  "ED vs SMs",
  "SMrz vs SMs"
)


# Same signed-rho palette used in Fig. 1.
colour_palette <- "bam"
colour_limits  <- c(-1, 1)
rho_colours    <- scico(50, palette = colour_palette)

# Indicator colours are retained for consistency with other scripts, although
# this figure does not directly map indicator identity with colour.
ind_colours <- c(
  AED  = "#E6A817",
  Et   = "#4DAF4A",
  ED   = "#377EB8",
  SMrz = "#8B4513",
  SMs  = "#D2691E"
)

# Typography and export size. 180 mm is a typical double-column figure width.
font_family  <- "Arial"
fig_width_mm <- 180
fig_height_mm <- 114
fig_dpi      <- 600

base_font_size       <- 6.1
axis_font_size       <- 5.3
axis_title_size      <- 5.7
scatter_title_size   <- 6.0
scatter_label_size   <- 1.95
heatmap_text_size    <- 2.15
legend_title_size    <- 6.0
legend_text_size     <- 5.6
tag_font_size        <- 8

scatter_fill_low  <- "#F2F2F2"
scatter_fill_high <- "#2166AC"
trend_colour      <- "#242424"
ref_line_colour   <- "#B8B8B8"
label_fill        <- grDevices::adjustcolor("#FFFFFF", alpha.f = 0.6)
label_text_colour <- "#1A1A1A"

# ------------------------------------------------------------------------------
# Helpers
# ------------------------------------------------------------------------------

check_file_exists <- function(file_path, description = "file") {
  if (!file.exists(file_path)) {
    stop("Required ", description, " not found: ", file_path, call. = FALSE)
  }
}

# ------------------------------------------------------------------------------
# Load and resample vegetation mask once
# ------------------------------------------------------------------------------

check_file_exists(file_veg_mask, "vegetation mask")
check_file_exists(heatmap_csv, "heatmap CSV")

message("Loading vegetation mask...")
reference_raster <- rast(
  file.path(dir_max, paste0("max_correlation_kndvi_", variables[1], ".nc"))
)[[1]]
veg_mask_raw <- rast(file_veg_mask)
veg_mask_res <- crop(resample(veg_mask_raw, reference_raster, method = "near"), e_ll)

# ------------------------------------------------------------------------------
# Helper: load rho_max stack with vegetation mask applied
# ------------------------------------------------------------------------------

load_rho_stack <- function(month_i) {
  stack <- list()
  for (v in variables) {
    f_nc <- file.path(dir_max, paste0("max_correlation_kndvi_", v, ".nc"))
    check_file_exists(f_nc, paste0("max correlation NetCDF for ", v))

    nc      <- rast(f_nc)
    cor_lyr <- nc[[grep("max_correlation", names(nc))]][[month_i]]
    cor_lyr <- crop(cor_lyr, e_ll)

    # Exclude non-vegetated pixels (C1 mask replaces aridity filter).
    cor_lyr[veg_mask_res == 0] <- NA

    names(cor_lyr) <- v
    stack[[v]] <- cor_lyr
  }
  rast(stack)
}

sample_pair_rho <- function(R_stack, xvar, yvar, n) {
  R_pair <- R_stack[[c(xvar, yvar)]]
  s <- terra::spatSample(R_pair, size = n, method = "random",
                         na.rm = TRUE, as.df = TRUE)
  s <- tidyr::drop_na(s)
  colnames(s) <- var_labels[c(xvar, yvar)]
  s
}

# ------------------------------------------------------------------------------
# Load heatmap table
# ------------------------------------------------------------------------------

message("Loading vegmask heatmap CSV...")
df_corr <- read_csv(heatmap_csv, show_col_types = FALSE) %>%
  mutate(
    Month = ifelse(Month %in% names(month_lookup), month_lookup[Month], Month),
    Month = factor(Month, levels = month_abbrev),
    Pair  = factor(Pair, levels = rev(pair_levels_topdown)),
    sig   = pval < 0.05
  )

get_heatmap_rho <- function(pair_label, month_label = "Jul") {
  val <- df_corr %>%
    filter(as.character(Pair) == pair_label,
           as.character(Month) == month_label) %>%
    pull(rho)

  if (length(val) != 1) {
    warning("Could not find unique heatmap rho for pair: ", pair_label,
            " and month: ", month_label)
    return(NA_real_)
  }
  round(val, 2)
}

# ------------------------------------------------------------------------------
# Panel a: scatterplots
# ------------------------------------------------------------------------------

message("Building panel a: scatterplots...")

make_scatter <- function(df, xvar, yvar, title, pair_label,
                         ref_slope = 1, month_label = "Jul",
                         show_ref_line = TRUE) {
  map_rho <- get_heatmap_rho(pair_label, month_label)

  rho_label <- paste0("\u03c1 = ", sprintf("%.2f", map_rho))

  p <- ggplot(df, aes(x = .data[[xvar]], y = .data[[yvar]])) +
    geom_hex(bins = 70) +
    scale_fill_gradient(
      low = scatter_fill_low,
      high = scatter_fill_high,
      trans = "sqrt",
      name = "Pixel count"
    )

  if (show_ref_line) {
    p <- p + geom_abline(
      slope = ref_slope,
      intercept = 0,
      linetype = "dashed",
      linewidth = 0.25,
      colour = ref_line_colour
    )
  }

  p <- p +
    geom_smooth(
      method = "lm",
      se = FALSE,
      colour = trend_colour,
      linewidth = 0.28
    ) +
    # Annotation box: draw the semi-transparent label background first,
    # with no visible text or border. A second text layer is then placed
    # on top. This avoids the device-dependent black outline produced by
    # geom_label()/annotate("label") when colour is also used for text.
    annotate(
      "label",
      x = -0.965,
      y = 0.965,
      label = rho_label,
      hjust = 0,
      vjust = 1,
      family = font_family,
      fontface = "bold",
      size = scatter_label_size,
      label.size = 0,
      label.r = unit(0.55, "mm"),
      label.padding = unit(0.55, "mm"),
      fill = label_fill,
      colour = NA
    ) +
    annotate(
      "text",
      x = -0.965,
      y = 0.965,
      label = rho_label,
      hjust = 0,
      vjust = 1,
      family = font_family,
      fontface = "bold",
      size = scatter_label_size,
      lineheight = 0.95,
      colour = label_text_colour
    ) +
    coord_fixed(xlim = c(-1, 1), ylim = c(-1, 1), expand = FALSE) +
    scale_x_continuous(breaks = c(-1, -0.5, 0, 0.5, 1)) +
    scale_y_continuous(breaks = c(-1, -0.5, 0, 0.5, 1)) +
    labs(
      title = title,
      x = paste0("\u03c1max(kNDVI, ", xvar, ")"),
      y = paste0("\u03c1max(kNDVI, ", yvar, ")")
    ) +
    theme_minimal(base_family = font_family, base_size = base_font_size) +
    theme(
      plot.title = element_text(
        face = "bold",
        hjust = 0.5,
        size = scatter_title_size,
        margin = margin(b = 1.0, unit = "mm")
      ),
      axis.title = element_text(size = axis_title_size, colour = "black"),
      axis.text = element_text(size = axis_font_size, colour = "#3A3A3A"),
      axis.ticks = element_line(linewidth = 0.22, colour = "#333333"),
      axis.ticks.length = unit(1.0, "mm"),
      panel.grid = element_blank(),
      legend.position = "none",
      plot.margin = margin(0.6, 1.1, 0.6, 1.1, unit = "mm")
    )

  p
}

cache_dir  <- file.path("outputs", "intermediate", "fig2_cache")
cache_file <- file.path(cache_dir, sprintf("scatter_jul_n%d_vegmask.rds", n_samp))

if (file.exists(cache_file)) {
  message("Loading cached scatter samples from: ", cache_file)
  sc             <- readRDS(cache_file)
  samp_smrz_sms <- sc$smrz_sms
  samp_aed_ed   <- sc$aed_ed
  samp_ed_sms   <- sc$ed_sms
  samp_aed_et   <- sc$aed_et
} else {
  message("Sampling scatter data (n = ", n_samp, ")...")
  R_jul          <- load_rho_stack(month_idx)
  samp_smrz_sms <- sample_pair_rho(R_jul, "SMrz", "SMs", n_samp)
  samp_aed_ed   <- sample_pair_rho(R_jul, "Ep",   "ED",  n_samp)
  samp_ed_sms   <- sample_pair_rho(R_jul, "ED",   "SMs", n_samp)
  samp_aed_et   <- sample_pair_rho(R_jul, "Ep",   "Et",  n_samp)

  dir.create(cache_dir, recursive = TRUE, showWarnings = FALSE)
  saveRDS(
    list(smrz_sms = samp_smrz_sms,
         aed_ed   = samp_aed_ed,
         ed_sms   = samp_ed_sms,
         aed_et   = samp_aed_et),
    cache_file
  )
  message("Cache saved to: ", cache_file)
}

scatter_pos <- make_scatter(
  samp_smrz_sms, "SMrz", "SMs",
  title = "SMrz vs SMs",
  pair_label = "SMrz vs SMs",
  ref_slope = 1
)

scatter_neg <- make_scatter(
  samp_aed_ed, "AED", "ED",
  title = "AED vs ED",
  pair_label = "AED vs ED",
  ref_slope = -1,
  show_ref_line = FALSE
)

scatter_mid <- make_scatter(
  samp_ed_sms, "ED", "SMs",
  title = "ED vs SMs",
  pair_label = "ED vs SMs",
  ref_slope = 1
)

scatter_low <- make_scatter(
  samp_aed_et, "AED", "Et",
  title = "AED vs Et",
  pair_label = "AED vs Et",
  ref_slope = 1
)

panel_a <- (scatter_pos | scatter_neg | scatter_mid | scatter_low) +
  plot_annotation(
    subtitle = "July",
    theme = theme(
      plot.subtitle = element_text(
        hjust = 0.5,
        size = 6.2,
        face = "italic",
        family = font_family,
        margin = margin(b = 1.2, unit = "mm")
      ),
      plot.margin = margin(0, 0, 0, 0, unit = "mm")
    )
  )

# ------------------------------------------------------------------------------
# Panel b: heatmap
# ------------------------------------------------------------------------------

message("Building panel b: heatmap...")

panel_b <- ggplot(df_corr, aes(x = Month, y = Pair, fill = rho)) +
  geom_tile(color = "white", linewidth = 0.18) +
  geom_text(
    aes(label = sprintf("%.2f", rho), color = sig),
    size = heatmap_text_size,
    family = font_family
  ) +
  scale_color_manual(values = c("TRUE" = "black", "FALSE" = "#8C8C8C"),
                     guide = "none") +
  scale_fill_gradientn(
    colours = rho_colours,
    limits = colour_limits,
    name = "Spatial Spearman's \u03c1 between \u03c1max maps",
    guide = guide_colourbar(
      direction = "horizontal",
      title.position = "top",
      title.hjust = 0.5,
      barwidth = unit(54, "mm"),
      barheight = unit(2.0, "mm")
    )
  ) +
  theme_minimal(base_family = font_family, base_size = base_font_size) +
  theme(
    axis.text.x = element_text(size = axis_font_size, angle = 0,
                               hjust = 0.5, colour = "black"),
    axis.text.y = element_text(size = axis_font_size, colour = "black"),
    axis.ticks = element_blank(),
    panel.grid = element_blank(),
    legend.position = "bottom",
    legend.title = element_text(size = legend_title_size, colour = "black"),
    legend.text = element_text(size = legend_text_size, colour = "black"),
    legend.key.width = unit(3, "mm"),
    legend.key.height = unit(2, "mm"),
    legend.box.margin = margin(0.5, 0, 0, 0, unit = "mm"),
    plot.margin = margin(1.2, 1.0, 0.0, 1.0, unit = "mm")
  ) +
  labs(x = NULL, y = NULL)

# ------------------------------------------------------------------------------
# Combine and save
# ------------------------------------------------------------------------------

message("Combining panels...")

fig <- wrap_elements(panel_a) / wrap_elements(panel_b) +
  plot_layout(heights = c(0.92, 1.18)) +
  plot_annotation(
    tag_levels = "a",
    theme = theme(
      plot.margin = margin(0.3, 0.5, 0.2, 0.5, unit = "mm"),
      plot.background = element_rect(fill = "white", colour = NA)
    )
  ) &
  theme(
    plot.tag = element_text(
      face = "bold",
      size = tag_font_size,
      family = font_family,
      colour = "black"
    )
  )

save_outputs <- function(fig) {
  ggsave(
    filename = out_file_tif,
    plot = fig,
    width = fig_width_mm,
    height = fig_height_mm,
    units = "mm",
    dpi = fig_dpi,
    bg = "white",
    device = "tiff",
    compression = "lzw"
  )

  ggsave(
     filename = out_file_pdf,
     plot = fig,
     width = fig_width_mm,
     height = fig_height_mm,
     units = "mm",
     dpi = fig_dpi,
     bg = "white",
     device = cairo_pdf
  )

  # ggsave(
  #   filename = out_file_png,
  #   plot = fig,
  #   width = fig_width_mm,
  #   height = fig_height_mm,
  #   units = "mm",
  #   dpi = 300,
  #   bg = "white"
  # )
}

save_outputs(fig)

message("Figure saved to:")
message("  ", out_file_tif)
# message("  ", out_file_pdf)
# message("  ", out_file_png)
message("Done.")
