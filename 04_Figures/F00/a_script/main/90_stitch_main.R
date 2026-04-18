# F00_QC — Pipeline Quality Control Composite
# 6-panel 3x2 grid summarizing stages 00_input -> 01_normalization -> 02_Imputation.
# Panels:
#   A Protein filter cascade        B Per-sample missingness      C Pre-normalization PCA
#   D Post-normalization PCA        E MAR/MNAR classification     F Observed vs imputed density
# Output: 04_Figures/F00/b_reports/MAIN_F00_QC_composite.{pdf,png} (300 x 200 mm)

setwd(rprojroot::find_rstudio_root_file())

suppressPackageStartupMessages({
  library(patchwork)
  library(ggplot2)
})

# Source panels — each creates a named plot object and saves a standalone PNG
source("04_Figures/F00/a_script/panels/panel_A.R")  # -> pA (filter cascade)
source("04_Figures/F00/a_script/panels/panel_B.R")  # -> pB (sample missingness)
source("04_Figures/F00/a_script/panels/panel_C.R")  # -> pC (pre-norm PCA)
source("04_Figures/F00/a_script/panels/panel_D.R")  # -> pD (post-norm PCA)
source("04_Figures/F00/a_script/panels/panel_E.R")  # -> pE (MAR/MNAR classification)
source("04_Figures/F00/a_script/panels/panel_F.R")  # -> pF (imputation density)

RPT_DIR <- "04_Figures/F00/b_reports"
DAT_DIR <- "04_Figures/F00/c_data"
dir.create(RPT_DIR, recursive = TRUE, showWarnings = FALSE)
WRITING_DIR <- "/Users/dtl0018/Library/CloudStorage/OneDrive-AuburnUniversity/YvO_writing/Figures/F00_QC"
BOX_DIR     <- "/Users/dtl0018/Library/CloudStorage/Box-Box/YvO_proteomics_manuscript_2026-04-16/04_Pipeline_QC"
dir.create(WRITING_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(BOX_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(BOX_DIR, "pdf"), recursive = TRUE, showWarnings = FALSE)
pdf_device <- get_pdf_device()

# Re-read imputation RDS for composite subtitle (matrix dimensions)
int_imp <- readRDS("02_Imputation/c_data/00_report_intermediates.rds")

# --- Composite ---
composite <- (pA | pB | pC) / (pD | pE | pF) +
  plot_annotation(
    title = "YvO Proteomics Pipeline \u2014 Quality Control Summary",
    subtitle = sprintf(
      "00_input -> 01_normalization (cyclic loess) -> 02_Imputation (missForest) | %s proteins x %s samples | %s",
      format(nrow(int_imp$mat), big.mark = ","),
      format(ncol(int_imp$mat), big.mark = ","),
      format(Sys.Date(), "%Y-%m-%d")),
    theme = theme(plot.title    = element_text(face = "bold", size = 14),
                  plot.subtitle = element_text(face = "italic", size = 10,
                                               color = "grey30"))
  ) &
  theme(plot.tag = element_text(face = "bold", size = 15))

COMP_W <- 300
COMP_H <- 200

ggsave(file.path(RPT_DIR, "MAIN_F00_QC_composite.pdf"), composite,
       width = COMP_W, height = COMP_H, units = "mm", device = pdf_device)
ggsave(file.path(RPT_DIR, "MAIN_F00_QC_composite.png"), composite,
       width = COMP_W, height = COMP_H, units = "mm", dpi = 300)
file.copy(file.path(RPT_DIR, "MAIN_F00_QC_composite.pdf"), file.path(WRITING_DIR, "01_main_composite_f00.pdf"), overwrite = TRUE)
file.copy(file.path(RPT_DIR, "MAIN_F00_QC_composite.png"), file.path(WRITING_DIR, "01_main_composite_f00.png"), overwrite = TRUE)
file.copy(file.path(RPT_DIR, "MAIN_F00_QC_composite.png"), file.path(BOX_DIR, "F00_main.png"), overwrite = TRUE)
file.copy(file.path(RPT_DIR, "MAIN_F00_QC_composite.pdf"), file.path(BOX_DIR, "pdf", "F00_main.pdf"), overwrite = TRUE)

message(sprintf("F00_QC composite saved -> %s", RPT_DIR))

# --- Supplementary Excel: one sheet per panel ---
source("04_Figures/shared/figure_supplement_helpers.R")

cat("=== F00 supplementary workbook ===\n")
f00_specs <- list(
  list(name = "panel_A", path = "04_Figures/F00/c_data/panel_A_filter_cascade.csv"),
  list(name = "panel_B", path = "04_Figures/F00/c_data/panel_B_sample_missingness.csv"),
  list(name = "panel_C", path = "04_Figures/F00/c_data/panel_C_pca_pre.csv"),
  list(name = "panel_D", path = "04_Figures/F00/c_data/panel_D_pca_post.csv"),
  list(name = "panel_E", path = "04_Figures/F00/c_data/panel_E_miss_classification.csv"),
  list(name = "panel_F", path = "04_Figures/F00/c_data/panel_F_imputation_summary.csv")
)
build_workbook(
  "04_Figures/F00/c_data/F00_supplementary.xlsx",
  title = "F00 \u2014 QC figure source data",
  description = "Panels A\u2013F source data for pipeline QC composite (normalization + imputation stages).",
  overview_df = data.frame(
    Sheet = c("panel_A", "panel_B", "panel_C",
              "panel_D", "panel_E", "panel_F"),
    Description = c(
      "Panel A: Protein filter cascade (retained vs removed counts per step)",
      "Panel B: Per-sample missingness (missing protein counts; outlier flag)",
      "Panel C: Pre-normalization PCA scores with Mahalanobis outlier flag",
      "Panel D: Post-normalization (cyclic loess) PCA scores with group/time factors",
      "Panel E: MAR/MNAR/Complete classification counts and percentages",
      "Panel F: Observed vs imputed intensity distribution summary stats + missForest OOB error"),
    stringsAsFactors = FALSE),
  sheet_specs = f00_specs
)
file.copy("04_Figures/F00/c_data/F00_supplementary.xlsx", file.path(WRITING_DIR, "04_supp_data_f00.xlsx"), overwrite = TRUE)
file.copy("04_Figures/F00/c_data/F00_supplementary.xlsx", file.path(BOX_DIR, "F00_source_data.xlsx"), overwrite = TRUE)
cleanup_after_workbook(f00_specs)
