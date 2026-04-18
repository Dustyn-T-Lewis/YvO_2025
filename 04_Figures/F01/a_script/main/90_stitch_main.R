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

source("04_Figures/F01/a_script/panels/panel_A.R")
source("04_Figures/F01/a_script/panels/panel_B.R")
source("04_Figures/F01/a_script/panels/panel_C.R")

RPT_DIR <- "04_Figures/F01/b_reports"
WRITING_DIR <- "/Users/dtl0018/Library/CloudStorage/OneDrive-AuburnUniversity/YvO_writing/Figures/F01"
BOX_DIR     <- "/Users/dtl0018/Library/CloudStorage/Box-Box/YvO_proteomics_manuscript_2026-04-16/02_Figures/F01_phenotype"
dir.create(WRITING_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(BOX_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(BOX_DIR, "pdf"), recursive = TRUE, showWarnings = FALSE)

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
  theme(plot.tag = element_text(face = "bold", size = 15))

COMP_W <- 215   # 90 + 78 + 42 + margins
COMP_H <- 155

ggsave(file.path(RPT_DIR, "MAIN_F01_composite.pdf"), composite,
       width = COMP_W, height = COMP_H, units = "mm", device = get_pdf_device())
ggsave(file.path(RPT_DIR, "MAIN_F01_composite.png"), composite,
       width = COMP_W, height = COMP_H, units = "mm", dpi = 300)
file.copy(file.path(RPT_DIR, "MAIN_F01_composite.pdf"), file.path(WRITING_DIR, "01_main_composite_f01.pdf"), overwrite = TRUE)
file.copy(file.path(RPT_DIR, "MAIN_F01_composite.png"), file.path(WRITING_DIR, "01_main_composite_f01.png"), overwrite = TRUE)
file.copy(file.path(RPT_DIR, "MAIN_F01_composite.png"), file.path(BOX_DIR, "F01_main.png"), overwrite = TRUE)
file.copy(file.path(RPT_DIR, "MAIN_F01_composite.pdf"), file.path(BOX_DIR, "pdf", "F01_main.pdf"), overwrite = TRUE)

cat("F01 composite done\n")

# --- Supplementary Excel: one sheet per panel ---
source("04_Figures/shared/figure_supplement_helpers.R")

cat("=== F01 supplementary workbook ===\n")
f01_specs <- list(
  list(name = "panel_A",         path = "04_Figures/F01/c_data/panel_A_training_volume.csv"),
  list(name = "panel_B",         path = "04_Figures/F01/c_data/panel_B_dxa_lbm.csv"),
  list(name = "panel_C",         path = "04_Figures/F01/c_data/panel_C_vl_thickness.csv"),
  list(name = "SUPP_panel_D",    path = "04_Figures/F01/c_data/supp/panel_D_deadlift_1rm.csv"),
  list(name = "SUPP_panel_E",    path = "04_Figures/F01/c_data/supp/panel_E_type_II_fcsa.csv"),
  list(name = "SUPP_panel_F",    path = "04_Figures/F01/c_data/supp/panel_F_type_I_fcsa.csv"),
  list(name = "data_dictionary", path = "04_Figures/F01/c_data/00_data_dictionary.csv")
)
build_workbook(
  "04_Figures/F01/c_data/F01_supplementary.xlsx",
  title = "F01 \u2014 Figure 1 source data",
  description = "Phenotype panels A\u2013C (main) and SUPP panels D\u2013F.",
  overview_df = data.frame(
    Sheet = c("panel_A", "panel_B", "panel_C",
              "SUPP_panel_D", "SUPP_panel_E", "SUPP_panel_F",
              "data_dictionary"),
    Description = c(
      "Panel A: Training volume summary stats by age group + t-test",
      "Panel B: DXA lean body mass paired t-test by age group",
      "Panel C: VL ultrasound thickness paired t-test by age group",
      "SUPP Panel D: Deadlift 1RM pre/post comparison",
      "SUPP Panel E: Type II fiber cross-sectional area",
      "SUPP Panel F: Type I fiber cross-sectional area",
      "Column definitions for each panel CSV"),
    stringsAsFactors = FALSE),
  sheet_specs = f01_specs
)
file.copy("04_Figures/F01/c_data/F01_supplementary.xlsx", file.path(WRITING_DIR, "04_supp_data_f01.xlsx"), overwrite = TRUE)
file.copy("04_Figures/F01/c_data/F01_supplementary.xlsx", file.path(BOX_DIR, "F01_source_data.xlsx"), overwrite = TRUE)
cleanup_after_workbook(f01_specs, extra_subdirs = c("04_Figures/F01/c_data/supp"))
