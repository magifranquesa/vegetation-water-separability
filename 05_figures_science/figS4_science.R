#!/usr/bin/env Rscript

# ==============================================================================
# Script: 05_figures_science/figS4_science.R
#
# Project: Aridity and timescale bound the separability of vegetation water stress
# Design: Nature-style, coherent with Fig. 1
#
# Fig. S4 (Fig. 4 in the original numbering): spatial agreement between the
# rho_max maps across indicators.
#   Heatmap of the Spearman correlation between the rho_max maps of every
#   indicator pair, across calendar months, computed over vegetated grid cells.
#   Same signed-rho colour scale as Fig. 1 (scico "bam", limits [-1, 1]), so a
#   single visual grammar for signed correlation values is kept across figures.
#
# Input:
#   outputs/tables/spatial_correlations_rho_max.csv
#
# Output:
#   outputs/figures_science/figS4_spatial_agreement_science.tif
#   (PDF output commented out)
#
# Dependencies: ggplot2, dplyr, readr, scico, grid
# ==============================================================================
suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(readr)
  library(scico)
  library(grid)
})

# ------------------------------------------------------------------------------
# Paths
# ------------------------------------------------------------------------------

config_file <- file.path("R", "config.R")
if (file.exists(config_file)) source(config_file)

heatmap_csv <- if (exists("paths") && !is.null(paths$tables)) {
  file.path(paths$tables, "spatial_correlations_rho_max.csv")
} else {
  file.path("outputs", "tables", "spatial_correlations_rho_max.csv")
}

out_file_tif <- file.path("outputs", "figures_science",
                          "figS4_spatial_agreement_science.tif")
# out_file_pdf <- sub("\\.tif$", ".pdf", out_file_tif)   # PDF output disabled
dir.create(dirname(out_file_tif), recursive = TRUE, showWarnings = FALSE)

# ------------------------------------------------------------------------------
# Figure options
# ------------------------------------------------------------------------------

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
rho_colours   <- scico(50, palette = "bam")
colour_limits <- c(-1, 1)

# Typography and export size. 180 mm is a typical double-column figure width;
# tune fig_height_mm to taste after the first render.
font_family   <- "Arial"
fig_width_mm  <- 180
fig_height_mm <- 55
fig_dpi       <- 300

base_font_size    <- 6.1
axis_font_size    <- 5.3
heatmap_text_size <- 2.15
legend_title_size <- 6.0
legend_text_size  <- 5.6

# ------------------------------------------------------------------------------
# Load heatmap table
# ------------------------------------------------------------------------------

if (!file.exists(heatmap_csv)) {
  stop("Required heatmap CSV not found: ", heatmap_csv, call. = FALSE)
}

message("Loading heatmap CSV...")
df_corr <- read_csv(heatmap_csv, show_col_types = FALSE) %>%
  mutate(
    Month = ifelse(Month %in% names(month_lookup), month_lookup[Month], Month),
    Month = factor(Month, levels = month_abbrev),
    Pair  = factor(Pair, levels = rev(pair_levels_topdown)),
    sig   = pval < 0.05
  )

# ------------------------------------------------------------------------------
# Heatmap
# ------------------------------------------------------------------------------

message("Building heatmap...")

fig <- ggplot(df_corr, aes(x = Month, y = Pair, fill = rho)) +
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
    name = "Spatial Spearman's ρ between ρmax maps",
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
    plot.margin = margin(1.5, 1.5, 1.5, 1.5, unit = "mm"),
    plot.background = element_rect(fill = "white", colour = NA)
  ) +
  labs(x = NULL, y = NULL)

# ------------------------------------------------------------------------------
# Save
# ------------------------------------------------------------------------------

ggsave(out_file_tif, fig, width = fig_width_mm, height = fig_height_mm, units = "mm",
       dpi = fig_dpi, bg = "white", device = "tiff", compression = "lzw")
#ggsave(out_file_pdf, fig, width = fig_width_mm, height = fig_height_mm, units = "mm",
#       dpi = fig_dpi, bg = "white", device = cairo_pdf)

message("Figure saved to:\n  ", out_file_tif)
message("Done.")
