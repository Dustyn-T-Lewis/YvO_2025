#!/usr/bin/env Rscript
# Build Supplementary Data S1: Pipeline Excel (raw -> normalized -> imputed -> DEP)
# Multi-sheet workbook with bold header, frozen pane, auto-width.

library(openxlsx)
library(readxl)
library(readr)

setwd(rprojroot::find_rstudio_root_file())

OUT_FILE <- "04_Figures/shared/S1_pipeline.xlsx"

add_sheet <- function(wb, name, data, description = NULL) {
  addWorksheet(wb, name)
  writeData(wb, name, data)
  hs <- createStyle(textDecoration = "bold")
  addStyle(wb, name, hs, rows = 1, cols = seq_len(ncol(data)), gridExpand = TRUE)
  freezePane(wb, name, firstRow = TRUE)
  setColWidths(wb, name, cols = seq_len(ncol(data)), widths = "auto")
  cat(sprintf("  + %s: %d rows x %d cols\n", name, nrow(data), ncol(data)))
}

wb <- createWorkbook()

# --- README sheet ---
readme <- data.frame(
  Sheet = c("S1a_metadata", "S1b_raw_intensities", "S1c_filter_log",
            "S1d_normalized", "S1e_mar_mnar", "S1f_imputed",
            "S1g_benchmark", "S1h_dep_all_contrasts", "S1i_dep_summary",
            "S1j_sensitivity"),
  Description = c(
    "Sample-level metadata: Col_ID, Group, Timepoint, phenotypes",
    "Raw DIA-NN log2 protein intensities (2,824 proteins x 62 samples)",
    "Filter cascade: HPA tissue -> blood contaminants -> dedup -> missingness",
    "Cycloess-normalized log2 intensities (2,113 proteins x 62 samples)",
    "Per-protein MAR/MNAR missingness classification + reliability flag",
    "missForest-imputed log2 intensities (2,113 proteins x 62 samples)",
    "23-method imputation benchmark: composite score, NRMSE, FC_rho, NES_rho",
    "All proteins x 5 contrasts: logFC, CI, t, p, FDR, Pi-score, sig_pi",
    "DEP counts per contrast at p<0.05, FDR<0.05, FDR<0.10, Pi<0.05",
    "Imputation + outlier sensitivity analyses"
  ),
  Source = c(
    "00_input/YvO_meta.xlsx",
    "00_input/YvO_raw.xlsx",
    "01_normalization/c_data/00_report_intermediates.rds (filter_log)",
    "01_normalization/c_data/02_normalized.csv",
    "02_Imputation/c_data/02_mar_mnar_classification.csv",
    "02_Imputation/c_data/01_imputed.csv",
    "02_Imputation/c_data/benchmark/04_composite_ranking.csv",
    "03_DEP/c_data/03_combined_results.csv",
    "03_DEP/c_data/02_DA_summary.csv",
    "03_DEP/c_data/09_imputation_sensitivity.csv + 11_outlier_sensitivity.csv"
  ),
  stringsAsFactors = FALSE
)
add_sheet(wb, "README", readme)

# --- S1a: metadata ---
meta <- as.data.frame(read_excel("00_input/YvO_meta.xlsx"))
add_sheet(wb, "S1a_metadata", meta)

# --- S1b: raw intensities ---
raw <- as.data.frame(read_excel("00_input/YvO_raw.xlsx"))
add_sheet(wb, "S1b_raw_intensities", raw)

# --- S1c: filter log ---
intermediates <- readRDS("01_normalization/c_data/00_report_intermediates.rds")
flog <- as.data.frame(intermediates$filter_log)
add_sheet(wb, "S1c_filter_log", flog)

# --- S1d: normalized ---
norm <- read_csv("01_normalization/c_data/02_normalized.csv", show_col_types = FALSE)
add_sheet(wb, "S1d_normalized", as.data.frame(norm))

# --- S1e: MAR/MNAR classification ---
marmnar <- read_csv("02_Imputation/c_data/02_mar_mnar_classification.csv",
                     show_col_types = FALSE)
add_sheet(wb, "S1e_mar_mnar", as.data.frame(marmnar))

# --- S1f: imputed ---
imp <- read_csv("02_Imputation/c_data/01_imputed.csv", show_col_types = FALSE)
add_sheet(wb, "S1f_imputed", as.data.frame(imp))

# --- S1g: benchmark ---
bench <- read_csv("02_Imputation/c_data/benchmark/04_composite_ranking.csv",
                   show_col_types = FALSE)
add_sheet(wb, "S1g_benchmark", as.data.frame(bench))

# --- S1h: DEP all contrasts ---
dep <- read_csv("03_DEP/c_data/03_combined_results.csv", show_col_types = FALSE)
add_sheet(wb, "S1h_dep_all_contrasts", as.data.frame(dep))

# --- S1i: DEP summary ---
da_summ <- read_csv("03_DEP/c_data/02_DA_summary.csv", show_col_types = FALSE)
add_sheet(wb, "S1i_dep_summary", as.data.frame(da_summ))

# --- S1j: sensitivity ---
imp_sens <- read_csv("03_DEP/c_data/09_imputation_sensitivity.csv",
                      show_col_types = FALSE)
out_sens <- read_csv("03_DEP/c_data/11_outlier_sensitivity.csv",
                      show_col_types = FALSE)
add_sheet(wb, "S1j_imputation_sens", as.data.frame(imp_sens))
add_sheet(wb, "S1j_outlier_sens", as.data.frame(out_sens))

# --- Save ---
saveWorkbook(wb, OUT_FILE, overwrite = TRUE)
cat(sprintf("\nSaved: %s (%.1f MB)\n", OUT_FILE, file.size(OUT_FILE) / 1e6))

# Copy to QC_reports (legacy location)
qc_dir <- "/Users/dtl0018/Library/CloudStorage/OneDrive-AuburnUniversity/YvO_writing/QC_reports"
if (dir.exists(qc_dir)) {
  file.copy(OUT_FILE, file.path(qc_dir, basename(OUT_FILE)), overwrite = TRUE)
  cat(sprintf("Copied to: %s\n", qc_dir))
}

# Copy to manuscript Figures/F00_QC (primary location for reviewers)
f00_qc_dir <- "/Users/dtl0018/Library/CloudStorage/OneDrive-AuburnUniversity/YvO_writing/Figures/F00_QC"
if (dir.exists(f00_qc_dir)) {
  file.copy(OUT_FILE, file.path(f00_qc_dir, basename(OUT_FILE)), overwrite = TRUE)
  cat(sprintf("Copied to: %s\n", f00_qc_dir))
}
