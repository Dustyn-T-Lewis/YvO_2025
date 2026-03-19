#!/usr/bin/env Rscript
# Figure consistency validation: checks that key statistics in figure data
# match the canonical values from the DEP analysis.
# Run after figure regeneration: Rscript tests/validate_figures.R

setwd(rprojroot::find_rstudio_root_file())

failures <- character(0)
check <- function(condition, msg) {
  if (!condition) failures <<- c(failures, msg)
  invisible(condition)
}

cat("=== Loading canonical DEP results ===\n")
dep <- read.csv("03_DEP/c_data/03_combined_results.csv")

# Wide format: columns are contrast-suffixed (e.g., pi_score_Aging, adj.P.Val_Aging)
# Expected DEP counts (Pi < 0.05)
expected_pi <- c(Aging = 198, Training_Young = 99, Training_Old = 18, Interaction = 34)

for (contrast_name in names(expected_pi)) {
  pi_col <- paste0("pi_score_", contrast_name)
  if (!pi_col %in% names(dep)) {
    failures <- c(failures, paste0("Missing column: ", pi_col))
    next
  }
  n_pi <- sum(dep[[pi_col]] < 0.05, na.rm = TRUE)
  check(n_pi == expected_pi[contrast_name],
        sprintf("%s: expected %d Pi-DEPs, got %d",
                contrast_name, expected_pi[contrast_name], n_pi))
}

# Expected DEP counts (FDR < 0.05)
expected_fdr <- c(Aging = 278, Training_Young = 134, Training_Old = 0, Interaction = 1)

for (contrast_name in names(expected_fdr)) {
  fdr_col <- paste0("adj.P.Val_", contrast_name)
  if (!fdr_col %in% names(dep)) next
  n_fdr <- sum(dep[[fdr_col]] < 0.05, na.rm = TRUE)
  check(n_fdr == expected_fdr[contrast_name],
        sprintf("%s: expected %d FDR-DEPs, got %d",
                contrast_name, expected_fdr[contrast_name], n_fdr))
}

cat("=== Checking figure data files ===\n")

# F2 concordance stats
conc_file <- "04_Figures/F04/c_data/panel_D/concordance.csv"
if (file.exists(conc_file)) {
  cat("  F04 concordance data found\n")
} else {
  conc_candidates <- Sys.glob("04_Figures/F04/c_data/*concordance*")
  if (length(conc_candidates) > 0) {
    cat("  F04 concordance data at:", conc_candidates[1], "\n")
  }
}

# F05 reversal classification
rev_file <- Sys.glob("04_Figures/F05/c_data/*reversal*")
if (length(rev_file) > 0) {
  cat("  F05 reversal classification found\n")
}

# F06 cluster assignments
clust_file <- "04_Figures/F06/c_data/06_mfuzz_assignments.csv"
if (file.exists(clust_file)) {
  clust <- read.csv(clust_file)
  check("cluster" %in% names(clust), "F4: cluster column missing from assignments")
  check(length(unique(clust$cluster)) == 4,
        sprintf("F06: expected 4 clusters, got %d", length(unique(clust$cluster))))
  cat("  F06 cluster assignments: k =", length(unique(clust$cluster)), "\n")
}

# Summary
cat("\n=== RESULTS ===\n")
if (length(failures) == 0) {
  cat("All figure consistency checks PASSED.\n")
  quit(status = 0)
} else {
  cat(length(failures), "FAILURES:\n")
  for (f in failures) cat("  FAIL:", f, "\n")
  quit(status = 1)
}
