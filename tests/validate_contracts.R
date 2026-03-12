#!/usr/bin/env Rscript
# Contract validation: checks that stage outputs exist and have expected properties.
# Run after any pipeline change: Rscript tests/validate_contracts.R
# Exit code 0 = all pass, 1 = failures found.

setwd(rprojroot::find_rstudio_root_file())

failures <- character(0)
check <- function(condition, msg) {
  if (!condition) failures <<- c(failures, msg)
  invisible(condition)
}

cat("=== Stage 01: Normalization ===\n")
check(file.exists("01_normalization/c_data/02_normalized.csv"),
      "01: normalized CSV missing")
check(file.exists("01_normalization/c_data/03_DAList_normalized.rds"),
      "01: normalized DAList missing")
if (file.exists("01_normalization/c_data/02_normalized.csv")) {
  norm <- read.csv("01_normalization/c_data/02_normalized.csv", nrows = 5)
  check(ncol(norm) > 50, "01: normalized CSV has too few columns (expected 60+)")
}
if (file.exists("01_normalization/c_data/03_DAList_normalized.rds")) {
  dal <- readRDS("01_normalization/c_data/03_DAList_normalized.rds")
  check("metadata" %in% names(dal), "01: DAList missing metadata slot")
  check(all(c("Col_ID", "Group", "Timepoint") %in% names(dal$metadata)),
        "01: DAList metadata missing required columns")
}

cat("=== Stage 02: Imputation ===\n")
check(file.exists("02_Imputation/c_data/01_imputed.csv"),
      "02: imputed CSV missing")
check(file.exists("02_Imputation/c_data/02_mar_mnar_classification.csv"),
      "02: MAR/MNAR classification missing")
if (file.exists("02_Imputation/c_data/01_imputed.csv")) {
  imp <- read.csv("02_Imputation/c_data/01_imputed.csv", nrows = 5)
  check(!any(is.na(imp[, -1])), "02: imputed CSV still has NAs in first 5 rows")
}

cat("=== Stage 03: DEP ===\n")
check(file.exists("03_DEP/c_data/03_combined_results.csv"),
      "03: combined results CSV missing")
check(file.exists("03_DEP/c_data/01_limma_DAList.rds"),
      "03: limma DAList missing")
if (file.exists("03_DEP/c_data/03_combined_results.csv")) {
  dep <- read.csv("03_DEP/c_data/03_combined_results.csv", nrows = 5)
  # Wide format: columns are suffixed by contrast (e.g., logFC_Aging, P.Value_Aging)
  expected_cols <- c("gene", "logFC_Aging", "P.Value_Aging", "adj.P.Val_Aging",
                     "logFC_Training_Young", "logFC_Training_Old",
                     "pi_score_Aging", "pi_score_Reversal")
  present <- expected_cols %in% names(dep)
  check(all(present),
        paste("03: combined results missing columns:",
              paste(expected_cols[!present], collapse = ", ")))
}

cat("=== Stage 04: Figures ===\n")
for (fig in paste0("F", 0:6)) {
  pdf_path <- file.path("04_Figures", fig, "b_reports",
                         paste0("Figure_", sub("F", "", fig), ".pdf"))
  check(file.exists(pdf_path),
        paste0("04: ", fig, " PDF missing at ", pdf_path))
}

# Summary
cat("\n=== RESULTS ===\n")
if (length(failures) == 0) {
  cat("All contract checks PASSED.\n")
  quit(status = 0)
} else {
  cat(length(failures), "FAILURES:\n")
  for (f in failures) cat("  FAIL:", f, "\n")
  quit(status = 1)
}
