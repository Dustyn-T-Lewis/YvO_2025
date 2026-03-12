#!/usr/bin/env Rscript
# 03_imputation_supp.R — Supplementary Excel workbook
# Input:  CSV files from c_data_v2/ (from 01_run_imputation.R)
# Output: c_data_v2/10_imputation_supp.xlsx

# --- Libraries ---------------------------------------------------------------
library(openxlsx)
library(readr)
library(tibble)

# --- Configuration -----------------------------------------------------------
setwd(rprojroot::find_rstudio_root_file())

DATA_DIR <- "02_Imputation/c_data_v2"

# --- Load data ---------------------------------------------------------------
bench_sum  <- read_csv(file.path(DATA_DIR, "03_benchmark_summary.csv"),
                       show_col_types = FALSE)
bench_df   <- read_csv(file.path(DATA_DIR, "04_benchmark_raw_iterations.csv"),
                       show_col_types = FALSE)
ext_sum    <- read_csv(file.path(DATA_DIR, "05_benchmark_extended.csv"),
                       show_col_types = FALSE)
bin_sum    <- read_csv(file.path(DATA_DIR, "06_benchmark_per_intensity.csv"),
                       show_col_types = FALSE)
miss_class <- read_csv(file.path(DATA_DIR, "02_mar_mnar_classification.csv"),
                       show_col_types = FALSE)
mnar_audit <- read_csv(file.path(DATA_DIR, "08_mnar_imputation_audit.csv"),
                       show_col_types = FALSE)
info_lines <- readLines(file.path(DATA_DIR, "09_imputation_summary.txt"))
info_df    <- tibble(
  Parameter = sub(" = .*", "", info_lines),
  Value     = sub(".* = ", "", info_lines))

# --- Helper ------------------------------------------------------------------
add_sheet <- function(wb, name, df) {
  hs <- createStyle(textDecoration = "bold", border = "Bottom", fgFill = "#DCE6F1")
  addWorksheet(wb, name)
  writeData(wb, name, df, headerStyle = hs)
  freezePane(wb, name, firstActiveRow = 2)
  setColWidths(wb, name, cols = seq_len(ncol(df)), widths = "auto")
}

# --- Build workbook ----------------------------------------------------------
wb <- createWorkbook()

readme <- tibble(
  Sheet = c("Benchmark", "Extended", "Per_Intensity", "Classification",
            "MNAR_Audit", "Iterations", "Summary"),
  Description = c(
    "Method ranking by NRMSE and PSS (12 methods x 20 iterations)",
    "Top 5 methods: NRMSE + PSS + per-intensity-tertile breakdown",
    "Per-intensity bin NRMSE for all methods (low/mid/high abundance)",
    "Per-protein Complete/MAR/MNAR classification with reliability flag",
    "MNAR pre/post means, shift, Cohen's d",
    "Raw per-iteration NRMSE and PSS for all methods",
    "Pipeline summary statistics"))
add_sheet(wb, "README", readme)
add_sheet(wb, "Benchmark", bench_sum)
add_sheet(wb, "Extended", ext_sum)
add_sheet(wb, "Per_Intensity", bin_sum)
add_sheet(wb, "Classification", miss_class)
add_sheet(wb, "MNAR_Audit", mnar_audit)
add_sheet(wb, "Iterations", bench_df)
add_sheet(wb, "Summary", info_df)

saveWorkbook(wb, file.path(DATA_DIR, "10_imputation_supp.xlsx"), overwrite = TRUE)

cat(sprintf("Supplementary workbook written: %s/10_imputation_supp.xlsx\n", DATA_DIR))
