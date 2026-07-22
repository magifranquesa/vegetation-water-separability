#!/usr/bin/env Rscript

# ==============================================================================
# Script: 05_figures/figS2.R
#
# Project: Aridity and timescale bound the separability of vegetation water stress
#
# Supplementary Fig. S2: proportions of the vegetated land surface in the five
# significance classes, using BH-FDR significance (p_adj < 0.05). Computes the
# proportions and draws the stacked-bar figure.
#
#   Classes (sum to 100% of the vegetated domain, area-weighted):
#     Significant positive     rho* > 0 and p_adj < 0.05
#     Non-significant positive rho* > 0 and p_adj >= 0.05
#     Non-significant negative rho* < 0 and p_adj >= 0.05
#     Significant negative     rho* < 0 and p_adj < 0.05
#     Not observed             vegetated cell with no valid rho*
#
# Inputs (READ-ONLY):
#   outputs/intermediate/absmax_spearman/abs_max_correlation_kndvi_{Ep,Et,ED,SMrz,SMs}.nc
#   outputs/intermediate/absmax_fdr/fdr_padj_{Ep,Et,ED,SMrz,SMs}.tif
#   outputs/intermediate/vegetation_mask_c1.tif
#
# Outputs:
#   outputs/tables/proportions_vegetated_fdr.csv
#   outputs/figures/figS2_significance_vegetated_proportions_bars_FDR.tif
#   outputs/figures/figS2_significance_vegetated_proportions_bars_FDR.pdf
#
# Run from repo root (AFTER 02_correlations/04_fdr_significance_mask.R):
#   Rscript 05_figures/figS2.R
#
# Dependencies: terra, ggplot2, dplyr, tidyr, grid
# ==============================================================================

suppressPackageStartupMessages({
  library(terra); library(ggplot2); library(dplyr); library(tidyr); library(grid)
})

config_file <- file.path("R", "config.R")
if (file.exists(config_file)) source(config_file)

dir_absmax <- if (exists("paths") && !is.null(paths$absmax_spearman)) paths$absmax_spearman else
  file.path("outputs", "intermediate", "absmax_spearman")
dir_fdr <- file.path("outputs", "intermediate", "absmax_fdr")
file_veg_mask <- if (exists("paths") && !is.null(paths$veg_mask_c1)) paths$veg_mask_c1 else
  file.path("outputs", "intermediate", "vegetation_mask_c1.tif")
tables_dir  <- if (exists("paths") && !is.null(paths$tables))  paths$tables  else file.path("outputs", "tables")
figures_dir <- if (exists("paths") && !is.null(paths$figures)) paths$figures else file.path("outputs", "figures")

out_csv <- file.path(tables_dir,  "proportions_vegetated_fdr.csv")
out_tif <- file.path(figures_dir, "figS2_significance_vegetated_proportions_bars_FDR.tif")
out_pdf <- sub("\\.tif$", ".pdf", out_tif)
for (f in c(out_csv, out_tif, out_pdf))
  if (file.exists(f)) stop("REFUSING TO OVERWRITE: ", f, call. = FALSE)

ALPHA <- 0.05
indic_files <- c(Ep = "AED", Et = "Et", ED = "ED", SMrz = "SMrz", SMs = "SMs")
fdrfile <- function(v) file.path(dir_fdr, sprintf("fdr_padj_%s.tif", v))
ncfile  <- function(v) file.path(dir_absmax, sprintf("abs_max_correlation_kndvi_%s.nc", v))

# ==============================================================================
# Compute FDR proportions per indicator x month
# ==============================================================================

message("Setting up grid, area and vegetated mask...")
r_ref <- rast(ncfile("Ep"), subds = "abs_max_correlation")[[1]]
if (is.na(crs(r_ref)) || crs(r_ref) == "") crs(r_ref) <- "EPSG:4326"
area_v <- as.vector(values(cellSize(r_ref, unit = "km", mask = FALSE)))
veg <- as.vector(values(resample(rast(file_veg_mask), r_ref, method = "near")))
veg <- !is.na(veg) & veg > 0
veg_area <- sum(area_v[veg])
veg_area_Mkm2 <- veg_area / 1e6

month_levels <- month.name
rows <- list()
for (v in names(indic_files)) {
  lab <- indic_files[[v]]
  message(sprintf("=== %s (%s) ===", lab, v))
  r_cor  <- rast(ncfile(v), subds = "abs_max_correlation")
  r_padj <- rast(fdrfile(v))
  for (m in 1:12) {
    cor <- as.vector(values(r_cor[[m]]))
    padj <- as.vector(values(r_padj[[m]]))
    obs  <- veg & is.finite(cor)
    sig  <- obs & is.finite(padj) & padj < ALPHA
    a <- function(sel) sum(area_v[sel]) / veg_area
    rows[[length(rows) + 1L]] <- data.frame(
      variable = lab, month = month_levels[m], month_id = m,
      sig_pos   = a(sig & cor > 0),
      sig_neg   = a(sig & cor < 0),
      nosig_pos = a(obs & !sig & cor > 0),
      nosig_neg = a(obs & !sig & cor < 0),
      not_observed = a(veg & !obs),
      domain_area_mill_km2 = round(veg_area_Mkm2, 3),
      stringsAsFactors = FALSE)
  }
}
df <- do.call(rbind, rows)
write.csv(df, out_csv, row.names = FALSE)
message("CSV written: ", out_csv)

# ==============================================================================
# Plot (reproduces figS2 style)
# ==============================================================================

font_family <- "sans"
fig_width_mm <- 180; fig_height_mm <- 95; fig_dpi <- 600
base_font_size <- 6.0; axis_font_size <- 5.4; axis_title_size <- 6.0
strip_font_size <- 6.8; legend_title_size <- 6.0; legend_text_size <- 5.8
bar_label_size <- 2.1

indicator_order <- c("AED", "Et", "ED", "SMrz", "SMs")
month_labels <- c("Jan","Feb","Mar","Apr","May","Jun","Jul","Aug","Sep","Oct","Nov","Dec")
correlation_classes <- c("Significant positive", "Non-significant positive",
                         "Non-significant negative", "Significant negative", "Not observed")
class_colours <- c("Significant positive" = "#114F02", "Non-significant positive" = "#ADCC89",
                   "Non-significant negative" = "#D88FC6", "Significant negative" = "#6C0B52",
                   "Not observed" = "#D3D3D3")
text_colours <- c("Significant positive" = "white", "Non-significant positive" = "black",
                  "Non-significant negative" = "black", "Significant negative" = "white",
                  "Not observed" = "grey40")

df_long <- df %>%
  mutate(sig_pos = sig_pos * 100, sig_neg = sig_neg * 100, nosig_pos = nosig_pos * 100,
         nosig_neg = nosig_neg * 100, no_data = not_observed * 100) %>%
  select(variable, month, month_id, sig_pos, sig_neg, nosig_pos, nosig_neg, no_data) %>%
  pivot_longer(cols = c(sig_pos, nosig_pos, nosig_neg, sig_neg, no_data),
               names_to = "class_code", values_to = "percent") %>%
  mutate(correlation_class = recode(class_code,
           "sig_pos" = "Significant positive", "nosig_pos" = "Non-significant positive",
           "nosig_neg" = "Non-significant negative", "sig_neg" = "Significant negative",
           "no_data" = "Not observed"),
         month = factor(month, levels = month.name, labels = month_labels),
         variable = factor(variable, levels = indicator_order),
         correlation_class = factor(correlation_class, levels = correlation_classes))

fig <- ggplot(df_long, aes(x = month, y = percent, fill = correlation_class)) +
  geom_col(width = 0.85, colour = "white", linewidth = 0.15) +
  geom_text(aes(label = ifelse(percent > 5, sprintf("%.0f", percent), ""), colour = correlation_class),
            position = position_stack(vjust = 0.5), size = bar_label_size, family = font_family) +
  facet_wrap(~ variable, ncol = 3) +
  scale_fill_manual(values = class_colours, name = "Correlation class", drop = FALSE) +
  scale_colour_manual(values = text_colours, guide = "none", drop = FALSE) +
  guides(fill = guide_legend(nrow = 2, byrow = TRUE)) +
  scale_y_continuous(expand = c(0, 0), breaks = seq(0, 100, 20), labels = function(x) paste0(x, "%")) +
  labs(x = NULL, y = "Vegetated land surface") +
  theme_minimal(base_size = base_font_size, base_family = font_family) +
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5, size = axis_font_size, colour = "#3A3A3A"),
        axis.text.y = element_text(size = axis_font_size, colour = "#3A3A3A"),
        axis.title.y = element_text(size = axis_title_size, colour = "black"),
        axis.ticks.x = element_line(linewidth = 0.22, colour = "#333333"),
        axis.ticks.y = element_line(linewidth = 0.22, colour = "#333333"),
        axis.ticks.length = unit(1.0, "mm"), panel.grid = element_blank(),
        strip.text = element_text(face = "bold", size = strip_font_size, colour = "black"),
        legend.position = "bottom", legend.title = element_text(size = legend_title_size, colour = "black"),
        legend.text = element_text(size = legend_text_size, colour = "black"),
        legend.key.width = unit(3.0, "mm"), legend.key.height = unit(2.6, "mm"),
        legend.spacing.x = unit(1.4, "mm"), legend.box.margin = margin(-0.6, 0, 0, 0, unit = "mm"),
        panel.spacing = unit(1.4, "mm"), plot.margin = margin(0.8, 0.8, 0.6, 0.8, unit = "mm"),
        plot.background = element_rect(fill = "white", colour = NA),
        panel.background = element_rect(fill = "white", colour = NA))

ggsave(out_tif, fig, width = fig_width_mm, height = fig_height_mm, units = "mm",
       dpi = fig_dpi, bg = "white", device = "tiff", compression = "lzw")
ggsave(out_pdf, fig, width = fig_width_mm, height = fig_height_mm, units = "mm",
       dpi = fig_dpi, bg = "white", device = cairo_pdf)
message("Supplementary Fig. S2 (FDR) saved to:\n  ", out_tif, "\n  ", out_pdf)
cat("\n=== DONE ===\n")
