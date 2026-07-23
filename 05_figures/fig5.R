#!/usr/bin/env Rscript

# ==============================================================================
# Script: 05_figures/fig5.R
#
# Project: Aridity and timescale bound the separability of vegetation water stress
#
# Figure 5: response timescales of vegetation-hydroclimate associations by
# aridity class, and supply-demand coupling as a function of accumulation
# timescale. Panels a/b use BH-FDR significance for the timescale curves and the
# significance-area bars; panel c is the coupling correlation (not thresholded).
#
#   a  Water-limited regime: mean response timescale by aridity class (Et+, SMrz+,
#      SMs+, ED+, AED-), four month facets, significance-area bars below.
#   b  Opposite-sign regime: AED+, ED- (dashed).
#   c  Coupling |rho(SMrz, AED)| vs accumulation timescale, by aridity class;
#      it rises with timescale, most in the semi-arid and dry sub-humid classes.
#
# Input (READ-ONLY; mtimes confirmed unchanged):
#   outputs/tables/timescale_bars_stress_fdr.csv       (panels a,b; from 04_analysis/09_timescale_bars.R)
#   outputs/tables/timescale_bars_nonstress_fdr.csv    (panels a,b; from 04_analysis/09_timescale_bars.R)
#   outputs/tables/supply_demand_coupling.csv          (panel c;   from 04_analysis/06_supply_demand_coupling.R)
#
# Output:
#   outputs/figures/fig5_memory_coupling_combined_FDR.tif
#   outputs/figures/fig5_memory_coupling_combined_FDR.pdf
#
# Run from repo root:
#   Rscript 05_figures/fig5.R
#
# Dependencies: dplyr, readr, ggplot2, patchwork, grid, scales, tibble
# ==============================================================================

suppressPackageStartupMessages({
  library(dplyr); library(readr); library(ggplot2)
  library(patchwork); library(grid); library(scales)
})

config_file <- file.path("R", "config.R")
if (file.exists(config_file)) source(config_file)

tables_dir  <- if (exists("paths") && !is.null(paths$tables))  paths$tables  else file.path("outputs", "tables")
figures_dir <- if (exists("paths") && !is.null(paths$figures)) paths$figures else file.path("outputs", "figures")
dir.create(figures_dir, recursive = TRUE, showWarnings = FALSE)

f_stress   <- file.path(tables_dir, "timescale_bars_stress_fdr.csv")      # FDR version
f_nonstr   <- file.path(tables_dir, "timescale_bars_nonstress_fdr.csv")   # FDR version
f_coupling <- file.path(tables_dir, "supply_demand_coupling.csv")
out_tif    <- file.path(figures_dir, "fig5_memory_coupling_combined_FDR.tif")
out_pdf    <- sub("\\.tif$", ".pdf", out_tif)

for (f in c(f_stress, f_nonstr, f_coupling))
  if (!file.exists(f)) stop("Required input not found: ", f, call. = FALSE)
for (f in c(out_tif, out_pdf))
  if (file.exists(f)) stop("REFUSING TO OVERWRITE: ", f, call. = FALSE)
mt_before <- file.info(c(f_stress, f_nonstr, f_coupling))$mtime

# ------------------------------------------------------------------------------
# Style
# ------------------------------------------------------------------------------

font_family   <- "sans"
fig_width_mm  <- 180
fig_height_mm <- 138          # a/b block + per-panel legends + faceted panel c
                              # (raised from 132 to fit aridity-class labels under panel a's bars)
fig_dpi       <- 600

base_font_size   <- 5.8
strip_font_size  <- 6.8
axis_font_size   <- 5.4
axis_title_size  <- 6.0
tag_font_size    <- 8
line_width_main  <- 0.50
point_size_main  <- 1.05
bar_label_size   <- 1.8
legend_text_size <- 5.8

ind_colours <- c(AED = "#D98C00", Et = "#4DAF4A", ED = "#9E1F63",
                 SMrz = "#084C8D", SMs = "#6BAED6")
face_colours <- c("Et+" = ind_colours[["Et"]], "SMrz+" = ind_colours[["SMrz"]],
                  "SMs+" = ind_colours[["SMs"]], "ED+" = ind_colours[["ED"]],
                  "AED-" = ind_colours[["AED"]], "AED+" = ind_colours[["AED"]],
                  "ED-" = ind_colours[["ED"]])
face_labels <- c("Et+" = "Et⁺", "SMrz+" = "SMrz⁺", "SMs+" = "SMs⁺",
                 "ED+" = "ED⁺", "AED-" = "AED⁻", "AED+" = "AED⁺", "ED-" = "ED⁻")

months_kept    <- c("January", "April", "July", "October")
months_index   <- c(1, 4, 7, 10)   # integer month codes in the coupling CSV
aridity_levels <- c("Arid", "Semi-arid", "Dry sub-humid", "Humid")
aridity_labels <- c("Arid" = "Arid", "Semi-arid" = "Semi-arid",
                    "Dry sub-humid" = "Dry\nsub-humid", "Humid" = "Humid")

bar_max_width <- 0.28; bar_h <- 0.32; y_bar_top <- 0.48; bar_x_offset <- 0.16

# panel c specifics -- aridity-class colours taken from panel b of figS4.R
# for cross-figure consistency.
# (That ramp is designed for map fills; the pale mid-classes are anchored here
#  with dark-bordered points so the lines remain legible.)
aridity_colours <- c("Arid" = "#D9B15F", "Semi-arid" = "#F3E5B8",
                     "Dry sub-humid" = "#CDEDEA", "Humid" = "#006D63")
aridity_lw <- c("Arid" = line_width_main * 1.2, "Semi-arid" = line_width_main * 1.5,
                "Dry sub-humid" = line_width_main * 1.2, "Humid" = line_width_main * 1.2)
scales_kept <- c(1, 3, 6, 9, 12)

# ------------------------------------------------------------------------------
# a/b data
# ------------------------------------------------------------------------------

message("Loading a/b CSVs and coupling CSV...")
df_stress    <- read_csv(f_stress, show_col_types = FALSE)
df_nonstress <- read_csv(f_nonstr, show_col_types = FALSE)

compute_monthly <- function(df, sign_char) {
  df %>%
    filter(month %in% months_kept, aridity_class %in% aridity_levels) %>%
    group_by(indicator, month, aridity_class) %>%
    summarise(mean_scale = sum(scale * perc) / sum(perc),
              sig_perc_total = mean(sig_perc_total), .groups = "drop") %>%
    mutate(face = paste0(indicator, sign_char))
}

mean_all <- bind_rows(
  compute_monthly(df_stress    %>% filter(sign == "pos"), "+"),
  compute_monthly(df_stress    %>% filter(sign == "neg"), "-"),
  compute_monthly(df_nonstress %>% filter(sign == "pos"), "+"),
  compute_monthly(df_nonstress %>% filter(sign == "neg"), "-")
) %>% mutate(indicator = sub("[+-]$", "", face)) %>%
  mutate(aridity_class = factor(aridity_class, levels = aridity_levels),
         month = factor(month, levels = months_kept))

make_bar_data <- function(df_mean, face_vec) {
  df_mean %>% filter(face %in% face_vec) %>%
    mutate(face_idx = match(as.character(face), face_vec),
           xpos = match(as.character(aridity_class), aridity_levels),
           bar_w = sig_perc_total / 100 * bar_max_width,
           ymax = y_bar_top - (face_idx - 1) * bar_h, ymin = ymax - bar_h,
           xmin = xpos - bar_x_offset, xmax = xmin + bar_w, ycenter = (ymin + ymax) / 2,
           bar_label = ifelse(sig_perc_total == 0, "0%",
                       ifelse(round(sig_perc_total) == 0, "<1%", paste0(round(sig_perc_total), "%"))),
           indicator = sub("[+-]$", "", face))
}

y_scale_line <- scale_y_continuous(limits = c(1.6, 8.2), breaks = c(2, 4, 6, 8), oob = scales::squish)
bar_strip_ylim <- function(n_faces) c(y_bar_top - n_faces * bar_h - 0.05, y_bar_top + 0.05)
x_scale_classes <- scale_x_continuous(breaks = seq_along(aridity_levels), labels = aridity_labels,
                                      limits = c(0.55, 4.55), expand = expansion(mult = c(0, 0)))

# ── Row A ─────────────────────────────────────────────────────────────────────
stress_order <- c("Et+", "SMrz+", "SMs+", "ED+", "AED-")
df_a <- mean_all %>% filter(face %in% stress_order) %>%
  mutate(face = factor(face, levels = stress_order), indicator = sub("[+-]$", "", face),
         xpos = match(as.character(aridity_class), aridity_levels))
bars_a <- make_bar_data(mean_all, stress_order) %>% mutate(face = factor(face, levels = stress_order))

row_a <- ggplot(df_a, aes(xpos, mean_scale, colour = face, group = face)) +
  geom_line(linewidth = line_width_main, show.legend = FALSE) +
  geom_point(size = point_size_main, show.legend = FALSE) +
  facet_wrap(~month, nrow = 1) +
  scale_colour_manual(values = face_colours[stress_order], guide = "none") +
  x_scale_classes + y_scale_line +
  labs(tag = "a", x = NULL, y = "Mean timescale (months)") +
  theme_minimal(base_family = font_family, base_size = base_font_size) +
  theme(strip.text = element_text(face = "bold", size = strip_font_size, family = font_family),
        panel.grid.major = element_line(linewidth = 0.25, colour = "#E8E8E8"),
        panel.grid.minor = element_blank(),
        axis.title.y = element_text(size = axis_title_size, colour = "black"),
        axis.text.y = element_text(size = axis_font_size, colour = "#3A3A3A"),
        axis.ticks.y = element_line(linewidth = 0.22, colour = "#333333"),
        axis.text.x = element_blank(), axis.ticks.x = element_blank(), legend.position = "none",
        plot.tag = element_text(face = "bold", size = tag_font_size, family = font_family,
                                colour = "black", margin = margin(r = 2.0, unit = "mm")),
        plot.margin = margin(0.4, 1.0, -0.4, 1.0, unit = "mm"))

bars_a_strip <- ggplot(bars_a) +
  geom_rect(aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, fill = face), show.legend = FALSE) +
  geom_text(aes(x = xmax + 0.025, y = ycenter, label = bar_label), size = bar_label_size,
            hjust = 0, colour = "#1A1A1A", family = font_family, show.legend = FALSE) +
  facet_wrap(~month, nrow = 1) +
  scale_fill_manual(values = face_colours[stress_order], guide = "none") +
  x_scale_classes + coord_cartesian(ylim = bar_strip_ylim(length(stress_order)), clip = "off") +
  theme_void(base_family = font_family, base_size = base_font_size) +
  theme(strip.text = element_blank(),
        axis.text.x = element_text(size = axis_font_size, colour = "#3A3A3A"),   # aridity-class labels under panel a's bars
        axis.ticks.x = element_line(linewidth = 0.22, colour = "#333333"),
        axis.ticks.length = unit(1.0, "mm"),
        plot.margin = margin(-1.6, 1.0, 1.2, 1.0, unit = "mm"))

# ── Row B ─────────────────────────────────────────────────────────────────────
mirror_order <- c("AED+", "ED-")
df_b <- mean_all %>% filter(face %in% mirror_order) %>%
  mutate(face = factor(face, levels = mirror_order), month = factor(month, levels = months_kept),
         indicator = sub("[+-]$", "", face), xpos = match(as.character(aridity_class), aridity_levels))
bars_b <- make_bar_data(mean_all, mirror_order) %>% mutate(face = factor(face, levels = mirror_order))

row_b <- ggplot(df_b, aes(xpos, mean_scale, colour = face, group = face)) +
  geom_line(linewidth = line_width_main, linetype = "dashed", show.legend = FALSE) +
  geom_point(size = point_size_main, show.legend = FALSE) +
  facet_wrap(~month, nrow = 1) +
  scale_colour_manual(values = face_colours[mirror_order], guide = "none") +
  x_scale_classes + y_scale_line +
  labs(tag = "b", x = NULL, y = "Mean timescale (months)") +
  theme_minimal(base_family = font_family, base_size = base_font_size) +
  theme(strip.text = element_blank(), strip.background = element_blank(),
        panel.grid.major = element_line(linewidth = 0.25, colour = "#E8E8E8"),
        panel.grid.minor = element_blank(),
        axis.title.y = element_text(size = axis_title_size, colour = "black"),
        axis.text.y = element_text(size = axis_font_size, colour = "#3A3A3A"),
        axis.ticks = element_line(linewidth = 0.22, colour = "#333333"),
        axis.ticks.length = unit(1.0, "mm"), axis.text.x = element_blank(),
        axis.ticks.x = element_blank(), legend.position = "none",
        plot.tag = element_text(face = "bold", size = tag_font_size, family = font_family,
                                colour = "black", margin = margin(r = 2.0, b = 0.6, unit = "mm")),
        plot.margin = margin(0.4, 1.0, -0.4, 1.5, unit = "mm"))

bars_b_strip <- ggplot(bars_b) +
  geom_rect(aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, fill = face), show.legend = FALSE) +
  geom_text(aes(x = xmax + 0.025, y = ycenter, label = bar_label), size = bar_label_size,
            hjust = 0, colour = "#1A1A1A", family = font_family, show.legend = FALSE) +
  facet_wrap(~month, nrow = 1) +
  scale_fill_manual(values = face_colours[mirror_order], guide = "none") +
  x_scale_classes + coord_cartesian(ylim = bar_strip_ylim(length(mirror_order)), clip = "off") +
  theme_void(base_family = font_family, base_size = base_font_size) +
  theme(strip.text = element_blank(),
        axis.text.x = element_text(size = axis_font_size, colour = "#3A3A3A"),
        axis.ticks.x = element_line(linewidth = 0.22, colour = "#333333"),
        axis.ticks.length = unit(1.0, "mm"), plot.margin = margin(-1.5, 1.0, 0.0, 1.0, unit = "mm"))

# ── per-panel legends: panel a's series below a, panel b's below b ────────────
make_face_legend <- function(faces, ltys, xs, xlim) {
  df <- tibble::tibble(face = factor(faces, levels = faces),
                       label = unname(face_labels[faces]), x = xs, lty = ltys)
  ggplot(df) +
    geom_segment(aes(x = x - 0.22, xend = x + 0.18, y = 1, yend = 1, colour = face, linetype = lty),
                 linewidth = line_width_main, lineend = "round") +
    geom_point(aes(x = x, y = 1, colour = face), size = point_size_main) +
    geom_text(aes(x = x + 0.30, y = 1, label = label), hjust = 0, vjust = 0.5,
              size = 1.8, family = font_family, colour = "black") +
    scale_colour_manual(values = face_colours[faces], guide = "none") +
    scale_linetype_manual(values = c(solid = "solid", dashed = "dashed"), guide = "none") +
    coord_cartesian(xlim = xlim, ylim = c(0.84, 1.16), clip = "off") +
    theme_void(base_family = font_family) +
    theme(plot.margin = margin(-0.2, 0, 0.2, 0, unit = "mm"))
}

# panel a: 5 solid series, centred under the full width
legend_a <- make_face_legend(stress_order, rep("solid", length(stress_order)),
                             xs = c(1.6, 3.1, 4.55, 5.85, 7.05), xlim = c(0.55, 8.85))
# panel b: 2 dashed series, centred
legend_b <- make_face_legend(mirror_order, rep("dashed", length(mirror_order)),
                             xs = c(3.6, 5.3), xlim = c(0.55, 8.85))

# ------------------------------------------------------------------------------
# Panel c — SM-AED coupling vs timescale
# ------------------------------------------------------------------------------

message("Building panel c (coupling vs timescale, by month)...")
# Faceted by the SAME four representative months as panels a/b (seasonality kept).
dc <- read.csv(f_coupling, stringsAsFactors = FALSE) %>%
  filter(supply == "SMrz", demand == "Ep",
         ai_class %in% aridity_levels, scale_months %in% scales_kept,
         month %in% months_index) %>%
  mutate(rho = mean_abs_rho,
         ai_class = factor(ai_class, levels = aridity_levels),
         month = factor(months_kept[match(month, months_index)], levels = months_kept))

y_lo <- floor(min(dc$rho) * 20) / 20; y_hi <- ceiling(max(dc$rho) * 20) / 20

panel_c <- ggplot(dc, aes(scale_months, rho, colour = ai_class, group = ai_class)) +
  geom_line(aes(linewidth = ai_class)) +
  geom_point(aes(fill = ai_class), shape = 21, colour = "#3A3A3A",
             stroke = 0.28, size = point_size_main) +              # dark border, size as in panel a
  facet_wrap(~month, nrow = 1) +
  scale_colour_manual(values = aridity_colours, guide = "none") +  # lines coloured, no legend
  scale_fill_manual(values = aridity_colours, name = NULL) +       # points drive the legend
  scale_linewidth_manual(values = aridity_lw, guide = "none") +
  scale_x_continuous(breaks = scales_kept, expand = expansion(mult = c(0.03, 0.04))) +
  scale_y_continuous(limits = c(y_lo, y_hi), breaks = seq(0, 1, 0.05),
                     labels = number_format(accuracy = 0.01), expand = expansion(mult = c(0.03, 0.05))) +
  labs(tag = "c", x = "Accumulation timescale (months)", y = expression("|" * rho * "(SMrz, AED)|")) +
  guides(fill = guide_legend(nrow = 1, override.aes = list(size = 2.2, stroke = 0.3))) +
  theme_minimal(base_family = font_family, base_size = base_font_size) +
  theme(strip.text = element_blank(),           # columns already labelled by panel a
        plot.tag = element_text(face = "bold", size = tag_font_size, family = font_family, colour = "black"),
        plot.tag.position = c(0.004, 0.99),
        panel.grid.major = element_line(linewidth = 0.25, colour = "#E8E8E8"),
        panel.grid.minor = element_blank(),
        axis.title = element_text(size = axis_title_size, colour = "black"),
        axis.text = element_text(size = axis_font_size, colour = "#3A3A3A"),
        axis.ticks = element_line(linewidth = 0.22, colour = "#333333"),
        axis.ticks.length = unit(1.0, "mm"),
        legend.position = "bottom",
        legend.text = element_text(size = legend_text_size - 0.2, colour = "black"),
        legend.key.width = unit(4.5, "mm"), legend.key.height = unit(2.5, "mm"),
        legend.box.margin = margin(-1.5, 0, 0, 0, unit = "mm"),
        plot.margin = margin(1.5, 1.0, 0.5, 1.5, unit = "mm"))

# ------------------------------------------------------------------------------
# Assemble and save
# ------------------------------------------------------------------------------

message("Assembling combined Fig. 5 (a / b / c)...")
fig <- row_a / bars_a_strip / legend_a / row_b / bars_b_strip / legend_b / panel_c +
  plot_layout(heights = c(0.90, 0.48, 0.12, 0.70, 0.14, 0.09, 1.05)) &
  theme(plot.margin = margin(0.5, 0.5, 0.5, 0.5, unit = "mm"),
        plot.background = element_rect(fill = "white", colour = NA))

ggsave(out_tif, fig, width = fig_width_mm, height = fig_height_mm, units = "mm",
       dpi = fig_dpi, bg = "white", device = "tiff", compression = "lzw")
#ggsave(out_pdf, fig, width = fig_width_mm, height = fig_height_mm, units = "mm",
#       dpi = fig_dpi, bg = "white", device = cairo_pdf)

message("\nFigure saved to:\n  ", out_tif, "\n  ", out_pdf)

mt_after <- file.info(c(f_stress, f_nonstr, f_coupling))$mtime
cat("\n=== File protection (inputs READ-ONLY) ===\n")
inputs <- c(f_stress, f_nonstr, f_coupling)
for (i in seq_along(inputs))
  cat(sprintf("  %-46s : %s\n", basename(inputs[i]),
              if (identical(mt_before[i], mt_after[i])) "UNCHANGED" else "*** CHANGED ***"))
cat("\n=== DONE ===\n")
