#!/usr/bin/env Rscript
# F01 Main — Training Volume (A) + DXA LBM (B) + VL Thickness (C)
# Produces single-column + double-column composites + xlsx

withr::local_dir(rprojroot::find_root(rprojroot::has_file("setup.R")))

library(dplyr)
library(ggplot2)
library(ggsignif)
library(patchwork)
library(cowplot)

source("04_Figures_v2/shared/style.R")
source("04_Figures_v2/shared/figure_supplement_helpers.R")

BASE <- "04_Figures_v2/F01"
RPT_PNG <- file.path(BASE, "b_reports", "main", "png")
RPT_PDF <- file.path(BASE, "b_reports", "main", "pdf")
PNL_PNG <- file.path(RPT_PNG, "panels")
PNL_PDF <- file.path(RPT_PDF, "panels")
DAT <- file.path(BASE, "c_data")
for (d in c(PNL_PNG, PNL_PDF, DAT)) dir.create(d, recursive = TRUE, showWarnings = FALSE)

TMPL <- "04_Figures_v2/F01/a_script/_prepost_template.R"

# Panels B and C via the shared template, which loads and coerces meta

cfg <- list(
  dv_col = "DXA_LBM_kg", y_label = "DXA LBM (kg)",
  delta_label = expression(bold(Delta ~ "DXA LBM (kg)")),
  title = "DXA Lean Body Mass", tag = "b", output_prefix = "pB",
  file_tag = "panel_B_dxa_lbm", audit_file = "panel_B_dxa_lbm.csv",
  file_prefix = "MAIN", rpt_png = PNL_PNG, rpt_pdf = PNL_PDF, dat = DAT,
  coerce_cols = TRUE, use_plotmath_subtitle = TRUE
)
source(TMPL)

cfg <- list(
  dv_col = "VL_thick_cm", y_label = "VL thickness (cm)",
  delta_label = expression(bold(Delta ~ "VL thickness (cm)")),
  title = "VL Thickness", tag = "c", output_prefix = "pC",
  file_tag = "panel_C_vl_thickness", audit_file = "panel_C_vl_thickness.csv",
  file_prefix = "MAIN", rpt_png = PNL_PNG, rpt_pdf = PNL_PDF, dat = DAT,
  coerce_cols = TRUE, use_plotmath_subtitle = TRUE,
  y_breaks = c(0, 0.5, 1.0), y_labels = c("0", ".5", "1")
)
source(TMPL)

# Panel A: Training Volume — reuses meta loaded and coerced by the template
tv_df <- meta |>
  filter(Timepoint == "Post") |>
  select(subject_key, Group, tv = Total_Training_Volume_kg) |>
  filter(!is.na(tv))

stats_A <- t.test(tv ~ Group, data = tv_df)
norm_sub <- sprintf("Welch t %s", fmt_p(stats_A$p.value))

audit_A <- tv_df |>
  summarise(
    n = n(), mean = mean(tv), sd = sd(tv), sem = sd / sqrt(n),
    shapiro_p = shapiro.test(tv)$p.value, .by = Group
  ) |>
  mutate(t_test_p = stats_A$p.value)
write.csv(audit_A, file.path(DAT, "panel_A_training_volume.csv"), row.names = FALSE)

tv_df$tv_scaled <- tv_df$tv / 1e5

pA <- ggplot(tv_df, aes(Group, tv_scaled, fill = Group)) +
  add_age_bands(1.5, 2.5) +
  geom_bar(stat = "summary", fun = mean, width = 0.6, color = "grey30", linewidth = 0.3) +
  geom_errorbar(stat = "summary", fun.data = mean_se, width = 0.2, linewidth = 0.4) +
  geom_jitter(width = 0.15, size = 1, alpha = 0.35, shape = 16, color = "grey30") +
  geom_signif(
    comparisons = list(c("Young", "Old")),
    annotations = fmt_p_plot(stats_A$p.value),
    parse = TRUE, textsize = 1.5, size = 0.3, tip_length = 0.02,
    y_position = max(tv_df$tv_scaled) * 1.15
  ) +
  scale_fill_manual(values = delta_bar_colors) +
  scale_x_discrete(expand = expansion(add = 0.3)) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.18))) +
  labs(
    title = "Training Volume", subtitle = norm_sub,
    y = expression(bold("Volume (" %*% 10^5 ~ "kg)")), x = NULL, tag = "a"
  ) +
  FIG_THEME +
  theme(legend.position = "none")

ggsave(file.path(PNL_PNG, "MAIN_panel_A_training_volume.png"), pA,
  width = 90, height = 150, units = "mm", dpi = 300
)
pA_title <- "Training Volume"
pA_subtitle <- paste0('italic("', norm_sub, '")')
pA <- strip_for_composite(pA)

# Single-column layout (85 × 125 mm) — for journal column width
sc_cfg <- list(
  w = 85, h = 125,
  tag_x = 0.02,
  ttl_x = 0.08,
  sub_off = 0.022,
  y = c(A = 0.979, B = 0.614, C = 0.314) # panel top positions (normalized)
)

# Double-column layout (178 × 75 mm) — for full-width display
dc_cfg <- list(
  w = 178, h = 75,
  tag_sz = 7, # pt — fixed (not via composite_text_sizes)
  ttl_sz = 6,
  sub_sz = 4,
  x_a = 0.010, # tag x for panel A column
  x_bc = 0.360, # tag x for panels B/C column
  ttl_off = 0.030, # title x offset from tag
  sub_off = 0.028,
  y_top = 0.984, # top row y
  y_mid = 0.516 # mid row y (panel C)
)

pB_comp <- (pB_left | pB_right) + plot_layout(widths = c(0.65, 0.35))
pC_comp <- (pC_left | pC_right) + plot_layout(widths = c(0.65, 0.35))

sc <- (pA / pB_comp / pC_comp) +
  plot_layout(heights = c(1.0, 0.8, 0.8)) &
  theme(plot.margin = margin_auto(10, 2, 2))

txt <- composite_text_sizes(sc_cfg$h)

sc <- ggdraw(sc) +
  draw_label("A", x = sc_cfg$tag_x, y = sc_cfg$y["A"], size = txt$tag, fontface = "bold", hjust = 0, vjust = 1) +
  draw_label(pA_title, x = sc_cfg$ttl_x, y = sc_cfg$y["A"], size = txt$title, fontface = "bold", hjust = 0, vjust = 1) +
  annotate("text",
    x = sc_cfg$ttl_x, y = sc_cfg$y["A"] - sc_cfg$sub_off, label = pA_subtitle,
    parse = TRUE, hjust = 0, vjust = 1, size = txt$subtitle / .pt, colour = "grey30"
  ) +
  draw_label("B", x = sc_cfg$tag_x, y = sc_cfg$y["B"], size = txt$tag, fontface = "bold", hjust = 0, vjust = 1) +
  draw_label(pB_title, x = sc_cfg$ttl_x, y = sc_cfg$y["B"], size = txt$title, fontface = "bold", hjust = 0, vjust = 1) +
  annotate("text",
    x = sc_cfg$ttl_x, y = sc_cfg$y["B"] - sc_cfg$sub_off, label = pB_subtitle,
    parse = TRUE, hjust = 0, vjust = 1, size = txt$subtitle / .pt, colour = "grey30"
  ) +
  draw_label("C", x = sc_cfg$tag_x, y = sc_cfg$y["C"], size = txt$tag, fontface = "bold", hjust = 0, vjust = 1) +
  draw_label(pC_title, x = sc_cfg$ttl_x, y = sc_cfg$y["C"], size = txt$title, fontface = "bold", hjust = 0, vjust = 1) +
  annotate("text",
    x = sc_cfg$ttl_x, y = sc_cfg$y["C"] - sc_cfg$sub_off, label = pC_subtitle,
    parse = TRUE, hjust = 0, vjust = 1, size = txt$subtitle / .pt, colour = "grey30"
  )

ggsave(file.path(RPT_PDF, "MAIN_F01_composite_single_col.pdf"), sc,
  width = sc_cfg$w, height = sc_cfg$h, units = "mm", device = get_pdf_device()
)
ggsave(file.path(RPT_PNG, "MAIN_F01_composite_single_col.png"), sc,
  width = sc_cfg$w, height = sc_cfg$h, units = "mm", dpi = 300
)

dc <- wrap_elements(full = pA + theme(plot.margin = margin_auto(8, 2, 10))) +
  wrap_elements(full = pB_comp & theme(plot.margin = margin_auto(4, 2, 2))) +
  wrap_elements(full = pC_comp & theme(plot.margin = margin_auto(4, 2, 4))) +
  plot_layout(design = "AB\nAC", widths = c(0.35, 0.65), heights = c(1, 1))

dc <- ggdraw(dc) +
  draw_label("A",
    x = dc_cfg$x_a, y = dc_cfg$y_top,
    size = dc_cfg$tag_sz, fontface = "bold", hjust = 0, vjust = 1
  ) +
  draw_label(pA_title,
    x = dc_cfg$x_a + dc_cfg$ttl_off, y = dc_cfg$y_top,
    size = dc_cfg$ttl_sz, fontface = "bold", hjust = 0, vjust = 1
  ) +
  annotate("text",
    x = dc_cfg$x_a + dc_cfg$ttl_off, y = dc_cfg$y_top - dc_cfg$sub_off,
    label = pA_subtitle,
    parse = TRUE, hjust = 0, vjust = 1, size = dc_cfg$sub_sz / .pt, colour = "grey30"
  ) +
  draw_label("B",
    x = dc_cfg$x_bc, y = dc_cfg$y_top,
    size = dc_cfg$tag_sz, fontface = "bold", hjust = 0, vjust = 1
  ) +
  draw_label(pB_title,
    x = dc_cfg$x_bc + dc_cfg$ttl_off, y = dc_cfg$y_top,
    size = dc_cfg$ttl_sz, fontface = "bold", hjust = 0, vjust = 1
  ) +
  annotate("text",
    x = dc_cfg$x_bc + dc_cfg$ttl_off, y = dc_cfg$y_top - dc_cfg$sub_off,
    label = pB_subtitle,
    parse = TRUE, hjust = 0, vjust = 1, size = dc_cfg$sub_sz / .pt, colour = "grey30"
  ) +
  draw_label("C",
    x = dc_cfg$x_bc, y = dc_cfg$y_mid,
    size = dc_cfg$tag_sz, fontface = "bold", hjust = 0, vjust = 1
  ) +
  draw_label(pC_title,
    x = dc_cfg$x_bc + dc_cfg$ttl_off, y = dc_cfg$y_mid,
    size = dc_cfg$ttl_sz, fontface = "bold", hjust = 0, vjust = 1
  ) +
  annotate("text",
    x = dc_cfg$x_bc + dc_cfg$ttl_off, y = dc_cfg$y_mid - dc_cfg$sub_off,
    label = pC_subtitle,
    parse = TRUE, hjust = 0, vjust = 1, size = dc_cfg$sub_sz / .pt, colour = "grey30"
  )

ggsave(file.path(RPT_PDF, "MAIN_F01_composite.pdf"), dc,
  width = dc_cfg$w, height = dc_cfg$h, units = "mm", device = get_pdf_device()
)
ggsave(file.path(RPT_PNG, "MAIN_F01_composite.png"), dc,
  width = dc_cfg$w, height = dc_cfg$h, units = "mm", dpi = 300
)

message("F01 main composites done")
