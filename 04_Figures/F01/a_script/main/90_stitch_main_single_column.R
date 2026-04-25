# F01 — Main Figure Composite (Panels A + B + C)
# J Physiol single-column (85 mm): A / B / C stacked vertically.
# B and C each have left (pre/post) + right (delta) sub-panels.
# Titles, subtitles, and tags placed at composite level via cowplot.

setwd(rprojroot::find_rstudio_root_file())

suppressPackageStartupMessages({
  library(patchwork)
  library(cowplot)
})

source("04_Figures/F01/a_script/main/panels/panel_A.R")
source("04_Figures/F01/a_script/main/panels/panel_B.R")
source("04_Figures/F01/a_script/main/panels/panel_C.R")

RPT_PNG <- "04_Figures/F01/b_reports/main/png"
RPT_PDF <- "04_Figures/F01/b_reports/main/pdf"

# Compose B and C sub-panels (clean — no titles/tags/legends)
pB_comp <- (pB_left | pB_right) + plot_layout(widths = c(0.65, 0.35))
pC_comp <- (pC_left | pC_right) + plot_layout(widths = c(0.65, 0.35))

# Stack A / B / C vertically with tight spacing
composite <- (pA / pB_comp / pC_comp) +
  plot_layout(heights = c(1.0, 0.8, 0.8)) &
  theme(plot.margin = margin(10, 2, 2, 2))

COMP_W <- 85    # J Physiol single-column
COMP_H <- 125
TAG_SZ <- composite_text_sizes(COMP_H)$tag
TTL_SZ <- composite_text_sizes(COMP_H)$title
SUB_SZ <- composite_text_sizes(COMP_H)$subtitle

# Row top positions (normalised to COMP_H):
#   Panel A top ≈ 1.0, Panel B top ≈ 1 - 1.0/2.6 ≈ 0.615, Panel C top ≈ 1 - 1.8/2.6 ≈ 0.308
Y_A <- 0.979;  Y_B <- 0.614;  Y_C <- 0.314
X_TAG <- 0.02
X_TTL <- 0.08   # title starts right of tag
X_SUB <- 0.08   # subtitle same x as title
TTL_OFFSET <- 0.000  # title just below tag line
TAG_DY_A   <- -0.003  # per-panel tag vertical offset
TAG_DY_BC  <- -0.002
SUB_OFFSET <- 0.022  # subtitle below title

composite <- ggdraw(composite) +
  # Panel A
  draw_label("A",          x = X_TAG, y = Y_A - TAG_DY_A,   size = TAG_SZ, fontface = "bold", hjust = 0, vjust = 1) +
  draw_label(pA_title,     x = X_TTL, y = Y_A - TTL_OFFSET, size = TTL_SZ, fontface = "bold", hjust = 0, vjust = 1) +
  annotate("text", x = X_SUB, y = Y_A - SUB_OFFSET, label = pA_subtitle, parse = TRUE, hjust = 0, vjust = 1, size = SUB_SZ / .pt, colour = "grey30") +
  # Panel B
  draw_label("B",          x = X_TAG, y = Y_B - TAG_DY_BC,  size = TAG_SZ, fontface = "bold", hjust = 0, vjust = 1) +
  draw_label(pB_title,     x = X_TTL, y = Y_B - TTL_OFFSET, size = TTL_SZ, fontface = "bold", hjust = 0, vjust = 1) +
  annotate("text", x = X_SUB, y = Y_B - SUB_OFFSET, label = pB_subtitle, parse = TRUE, hjust = 0, vjust = 1, size = SUB_SZ / .pt, colour = "grey30") +
  # Panel C
  draw_label("C",          x = X_TAG, y = Y_C - TAG_DY_BC,  size = TAG_SZ, fontface = "bold", hjust = 0, vjust = 1) +
  draw_label(pC_title,     x = X_TTL, y = Y_C - TTL_OFFSET, size = TTL_SZ, fontface = "bold", hjust = 0, vjust = 1) +
  annotate("text", x = X_SUB, y = Y_C - SUB_OFFSET, label = pC_subtitle, parse = TRUE, hjust = 0, vjust = 1, size = SUB_SZ / .pt, colour = "grey30")

ggsave(file.path(RPT_PDF, "MAIN_F01_composite.pdf"), composite,
       width = COMP_W, height = COMP_H, units = "mm", device = get_pdf_device())
ggsave(file.path(RPT_PNG, "MAIN_F01_composite.png"), composite,
       width = COMP_W, height = COMP_H, units = "mm", dpi = 300)
cat("F01 composite done\n")

# --- Supplementary Excel: one sheet per panel ---
source("04_Figures/shared/figure_supplement_helpers.R")

cat("=== F01 supplementary workbook ===\n")
f01_specs <- list(
  list(name = "panel_A",         path = "04_Figures/F01/c_data/panel_A_training_volume.csv"),
  list(name = "panel_B",         path = "04_Figures/F01/c_data/panel_B_dxa_lbm.csv"),
  list(name = "panel_C",         path = "04_Figures/F01/c_data/panel_C_vl_thickness.csv"),
  list(name = "SUPP_panel_A",    path = "04_Figures/F01/c_data/supp/panel_A_deadlift_1rm.csv"),
  list(name = "SUPP_panel_B",    path = "04_Figures/F01/c_data/supp/panel_B_type_II_fcsa.csv"),
  list(name = "SUPP_panel_C",    path = "04_Figures/F01/c_data/supp/panel_C_type_I_fcsa.csv")
)
build_workbook(
  "04_Figures/F01/c_data/F01_supplementary.xlsx",
  title = "F01 \u2014 Figure 1 source data",
  description = "Phenotype panels A\u2013C (main) and SUPP panels A\u2013C.",
  overview_df = data.frame(
    Sheet = c("panel_A", "panel_B", "panel_C",
              "SUPP_panel_A", "SUPP_panel_B", "SUPP_panel_C"),
    Description = c(
      "Panel A: Training volume summary stats by age group + t-test",
      "Panel B: DXA lean body mass paired t-test by age group",
      "Panel C: VL ultrasound thickness paired t-test by age group",
      "SUPP Panel A: Deadlift 1RM pre/post comparison",
      "SUPP Panel B: Type II fiber cross-sectional area",
      "SUPP Panel C: Type I fiber cross-sectional area"),
    stringsAsFactors = FALSE),
  sheet_specs = f01_specs
)
cleanup_after_workbook(f01_specs, extra_subdirs = c("04_Figures/F01/c_data/supp"))
