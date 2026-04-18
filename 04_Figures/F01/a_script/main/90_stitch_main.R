# F01 — Main Figure Composite (Panels A + B + C)
# Flat 2×3 grid for true plot-area alignment across columns.
# Col 1: Panel A (spanning both rows)
# Col 2: B_left (row 1), C_left (row 2)   — pre/post bar charts
# Col 3: B_right (row 1), C_right (row 2) — delta bar charts
#
# Panel scripts export pA, pB_left, pB_right, pC_left, pC_right as
# individual ggplot objects (composed into pB/pC for standalone export).

setwd(rprojroot::find_rstudio_root_file())

suppressPackageStartupMessages(library(patchwork))

source("04_Figures/F01/a_script/main/panels/panel_A.R")
source("04_Figures/F01/a_script/main/panels/panel_B.R")
source("04_Figures/F01/a_script/main/panels/panel_C.R")

RPT_PNG <- "04_Figures/F01/b_reports/main/png"
RPT_PDF <- "04_Figures/F01/b_reports/main/pdf"

# Add invisible title/subtitle spacers to right sub-panels so their
# plot areas align vertically with the titled left sub-panels.
# B_left/C_left have 1-line title + 2-line subtitle; match that height.
pB_right_s <- pB_right + labs(title = " ", subtitle = " \n ")
pC_right_s <- pC_right + labs(title = " ", subtitle = " \n ")

# Trim bottom margins on B/C left panels to tighten the gap between rows
# while keeping B+C top/bottom aligned with panel A
pB_left  <- pB_left  + theme(plot.margin = margin(5, 5, 2, 5))
pC_left  <- pC_left  + theme(plot.margin = margin(2, 5, 5, 5))
pB_right_s <- pB_right_s + theme(plot.margin = margin(5, 5, 2, 5))
pC_right_s <- pC_right_s + theme(plot.margin = margin(2, 5, 5, 5))

# Flat design: A spans rows, B/C decomposed into left+right columns.
# Patchwork assigns plots to design chars in ALPHABETICAL order,
# so use ABC (row 1: A, B_left, B_right) / ADE (row 2: A, C_left, C_right).
design <- "ABC\nADE"

composite <- pA + pB_left + pB_right_s + pC_left + pC_right_s +
  plot_layout(design = design,
              widths  = c(90, 78, 42),
              heights = c(1, 1)) &
  theme(plot.tag      = element_text(face = "bold", size = 15,
                                      margin = margin(t = 3, r = 0, b = 0, l = 4)),
        plot.title    = element_text(margin = margin(b = 1)),
        plot.subtitle = element_text(margin = margin(t = 0)))

COMP_W <- 215   # 90 + 78 + 42 + margins
COMP_H <- 155

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
