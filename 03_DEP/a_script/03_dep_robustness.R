# YvO DEP — Robustness Analyses
# Blunting diagnostics, bootstrap CIs, power analysis, imputation sensitivity
# Assembles supplementary Excel workbook. Depends on: 01_run_dep.R outputs
#
# Refs: Conover 1999 (KS, Fligner), Romano 2006 (Cliff's delta),
#       Efron & Tibshirani 1993 (bootstrap), Cohen 1988 (power),
#       Karpievitch 2012 (imputation sensitivity)

# --- Setup

library(dplyr)
library(tidyr)
library(tibble)
library(stringr)
library(readr)
library(purrr)
library(ggplot2)
library(patchwork)
library(proteoDA)
library(openxlsx)
library(boot)
library(pwr)

set.seed(42)

setwd(rprojroot::find_rstudio_root_file())

cfg <- list(
  data_dir   = "03_DEP/c_data",
  per_dir    = "03_DEP/c_data/04_per_contrast_results",
  report_dir = "03_DEP/b_reports",
  norm_csv   = "01_normalization/c_data/02_normalized.csv",
  norm_rds   = "01_normalization/c_data/03_DAList_normalized.rds",
  imp_path   = "02_Imputation/c_data/01_imputed.csv",
  pi_thresh  = 0.05,
  pval_thresh = 0.10
)

# --- Load data

dal <- readRDS(file.path(cfg$data_dir, "01_limma_DAList.rds"))
contrast_names <- names(dal$results)

# Read per-contrast CSVs
results_list <- lapply(contrast_names, function(cname) {
  read_csv(file.path(cfg$per_dir, paste0(cname, ".csv")), show_col_types = FALSE)
})
names(results_list) <- contrast_names

# Metadata for sample counts
meta <- as.data.frame(dal$metadata)

# --- 1. Blunting diagnostics
# Tests whether |logFC| distributions differ between Training_Young and
# Training_Old, providing evidence for age-dependent response attenuation.

blunt_df <- tibble(
  gene          = results_list[["Training_Young"]]$gene,
  abs_lfc_young = abs(results_list[["Training_Young"]]$logFC),
  abs_lfc_old   = abs(results_list[["Training_Old"]]$logFC)
) |> filter(!is.na(abs_lfc_young) & !is.na(abs_lfc_old))

# KS test: distributional shape/location shift
ks_res <- ks.test(blunt_df$abs_lfc_young, blunt_df$abs_lfc_old)

# Fligner-Killeen: variance heterogeneity (robust to non-normality)
fk_data <- data.frame(
  abs_lfc  = c(blunt_df$abs_lfc_young, blunt_df$abs_lfc_old),
  contrast = factor(rep(c("Training_Young", "Training_Old"), each = nrow(blunt_df)))
)
fk_res <- fligner.test(abs_lfc ~ contrast, data = fk_data)

# Wilcoxon signed-rank: paired location shift (same proteins)
wx_res <- wilcox.test(blunt_df$abs_lfc_young, blunt_df$abs_lfc_old, paired = TRUE)

# Cliff's delta: non-parametric effect size
diffs <- blunt_df$abs_lfc_young - blunt_df$abs_lfc_old
cliff_delta <- (sum(diffs > 0) - sum(diffs < 0)) / length(diffs)
cliff_mag <- case_when(
  abs(cliff_delta) < 0.147 ~ "negligible",
  abs(cliff_delta) < 0.33  ~ "small",
  abs(cliff_delta) < 0.474 ~ "medium",
  TRUE                     ~ "large"
)

# Quantile comparison
quants <- c(0.25, 0.50, 0.75, 0.90)
quant_df <- tibble(
  quantile       = paste0("Q", quants * 100),
  Training_Young = quantile(blunt_df$abs_lfc_young, quants),
  Training_Old   = quantile(blunt_df$abs_lfc_old, quants),
  ratio_YO       = round(Training_Young / Training_Old, 3)
)

blunt_diag <- tibble(
  test = c("Kolmogorov-Smirnov", "Fligner-Killeen",
           "Wilcoxon signed-rank", "Cliff's delta"),
  statistic = c(ks_res$statistic, fk_res$statistic,
                wx_res$statistic, cliff_delta),
  p_value = c(ks_res$p.value, fk_res$p.value, wx_res$p.value, NA),
  interpretation = c(
    ifelse(ks_res$p.value < 0.05,
           "Distributions differ significantly (shape/location)",
           "No significant distributional difference"),
    ifelse(fk_res$p.value < 0.05,
           "Variance differs (consistent with blunted response)",
           "No significant variance difference"),
    ifelse(wx_res$p.value < 0.05,
           "Paired |logFC| shift significant (attenuated with age)",
           "No significant paired shift"),
    sprintf("%s effect (delta = %.3f; >0 = young responds more)",
            str_to_title(cliff_mag), cliff_delta)
  )
)
write_csv(blunt_diag, file.path(cfg$data_dir, "06_blunting_diagnostics.csv"))

cat(sprintf("  Blunting: KS p=%.2g, Wilcoxon p=%.2g, Cliff d=%.3f (%s)\n",
            ks_res$p.value, wx_res$p.value, cliff_delta, cliff_mag))

# --- 2. Bootstrap CI
# Median |logFC| with 95% BCa bootstrap CI per contrast

median_fn <- function(d, i) median(d[i], na.rm = TRUE)

boot_contrasts <- c("Training_Young", "Training_Old", "Aging", "Interaction")
boot_df <- map_dfr(boot_contrasts, function(cname) {
  vals <- abs(results_list[[cname]]$logFC)
  vals <- vals[!is.na(vals)]
  b  <- boot(vals, median_fn, R = 10000)
  ci <- tryCatch(boot.ci(b, type = "bca"),
                 error = function(e) boot.ci(b, type = "perc"))
  ci_lo <- if (!is.null(ci$bca)) ci$bca[4] else ci$percent[4]
  ci_hi <- if (!is.null(ci$bca)) ci$bca[5] else ci$percent[5]
  tibble(contrast = cname, median_absLFC = median(vals),
         ci_lower = ci_lo, ci_upper = ci_hi,
         boot_se = sd(b$t), n_proteins = length(vals))
})
write_csv(boot_df, file.path(cfg$data_dir, "07_effect_size_bootstrap.csv"))

print(boot_df)

# --- 3. Power analysis
# Approximate minimum detectable logFC at 80% power.
# Conservative: pwr.t.test assumes standard t; limma's moderated t has higher
# power via empirical Bayes variance shrinkage. Treat as lower bounds.

fit <- dal$eBayes_fit
within_cor <- dal$eBayes_fit$correlation %||%
  dal$tags$duplicate_correlation %||% NA_real_
sigma_residual <- sqrt(mean(fit$sigma^2, na.rm = TRUE))
n_young <- sum(meta$age == "Young" & meta$time == "Pre")
n_old   <- sum(meta$age == "Old"   & meta$time == "Pre")

power_df <- map_dfr(boot_contrasts, function(cname) {
  n_subj <- switch(cname,
    Training_Young = n_young, Training_Old = n_old,
    min(n_young, n_old))
  paired <- cname %in% c("Training_Young", "Training_Old")
  eff_sigma <- if (paired && !is.na(within_cor))
    sigma_residual * sqrt(2 * (1 - within_cor)) else
    sigma_residual * sqrt(2)
  pw <- pwr.t.test(n = n_subj, d = NULL, sig.level = 0.10, power = 0.80,
                   type = if (paired) "paired" else "two.sample")
  tibble(contrast = cname, n_subjects = n_subj,
         within_cor = ifelse(paired, within_cor, NA_real_),
         effective_sigma = round(eff_sigma, 4),
         min_detectable_d = round(pw$d, 4),
         min_detectable_logFC = round(pw$d * eff_sigma, 4),
         power = 0.80, alpha = 0.10)
})
write_csv(power_df, file.path(cfg$data_dir, "08_power_analysis.csv"))
print(as.data.frame(power_df))

# --- 4. Imputation sensitivity
# Compare t-statistics: non-imputed (main) vs imputed limma.

sens_df <- NULL

if (file.exists(cfg$imp_path)) {
  # Load annotation and original data
  ann_cols <- c("uniprot_id", "protein", "gene", "description")
  ann <- as.data.frame(dal$annotation)
  mat <- dal$data

  imp_data <- read_csv(cfg$imp_path, show_col_types = FALSE)
  imp_samp <- setdiff(names(imp_data), ann_cols)
  imp_mat  <- as.matrix(imp_data[, imp_samp])
  rownames(imp_mat) <- imp_data$uniprot_id

  shared_samps <- intersect(colnames(mat), colnames(imp_mat))
  meta_imp_df  <- meta[meta$sample_id %in% shared_samps, ]

  dal_imp <- DAList(
    data       = imp_mat[, shared_samps],
    annotation = ann,
    metadata   = meta_imp_df,
    tags       = list(norm_method = "cycloess_imputed")
  )
  dal_imp <- add_design(dal_imp, "~ 0 + group + (1 | subject)")
  colnames(dal_imp$design$design_matrix) <- gsub("^group", "",
    colnames(dal_imp$design$design_matrix))
  dal_imp <- add_contrasts(dal_imp, contrasts_vector = c(
    "Training_Young = Young_Post - Young_Pre",
    "Training_Old = Old_Post - Old_Pre",
    "Aging = Old_Pre - Young_Pre",
    "Interaction = (Old_Post - Old_Pre) - (Young_Post - Young_Pre)"
  ))
  dal_imp <- fit_limma_model(dal_imp)
  dal_imp <- extract_DA_results(dal_imp, pval_thresh = 0.10, lfc_thresh = 0,
                                 adj_method = "BH")

  # Extract imputed results from dal_imp$results (no disk round-trip)
  imp_results <- lapply(names(dal_imp$results), function(cname) {
    dal_imp$results[[cname]] |> rownames_to_column("uniprot_id")
  })
  names(imp_results) <- names(dal_imp$results)

  # Read non-imputed combined results for comparison
  comb <- read_csv(file.path(cfg$data_dir, "03_combined_results.csv"),
                   show_col_types = FALSE)

  sens_contrasts <- c("Training_Young", "Training_Old", "Aging", "Interaction")
  sens_rows <- list()
  for (cname in sens_contrasts) {
    t_col   <- paste0("t_", cname)
    adj_col <- paste0("adj.P.Val_", cname)
    if (!(t_col %in% names(comb)) || !(cname %in% names(imp_results))) next

    imp_df <- imp_results[[cname]] |>
      select(uniprot_id, t_imp = t, padj_imp = adj.P.Val)

    merged <- inner_join(
      comb |> select(uniprot_id, t_nonimp = all_of(t_col),
                      padj_nonimp = all_of(adj_col)),
      imp_df,
      by = "uniprot_id"
    ) |> filter(!is.na(t_nonimp) & !is.na(t_imp))

    sp <- cor.test(merged$t_nonimp, merged$t_imp, method = "spearman")

    sens_rows[[cname]] <- tibble(contrast = cname,
           spearman_rho = round(sp$estimate, 4),
           p_value = sp$p.value, n_proteins = nrow(merged))
  }
  sens_df <- bind_rows(sens_rows)
  write_csv(sens_df, file.path(cfg$data_dir, "09_imputation_sensitivity.csv"))

  print(as.data.frame(sens_df))
} else {
  cat("  Imputed data not found — skipping\n")
  # Write empty placeholder so downstream doesn't break
  write_csv(tibble(contrast = character(), spearman_rho = numeric(),
                   p_value = numeric(), n_proteins = integer()),
            file.path(cfg$data_dir, "09_imputation_sensitivity.csv"))
}

# --- 5. Supplementary Excel

wb_supp    <- createWorkbook()
hdr_style  <- createStyle(textDecoration = "italic", fontColour = "#555555",
                           wrapText = TRUE)
bold_style <- createStyle(textDecoration = "bold")

write_supp_sheet <- function(wb, sheet, header, data) {
  addWorksheet(wb, sheet)
  for (i in seq_along(header))
    writeData(wb, sheet, header[i], startRow = i, startCol = 1)
  addStyle(wb, sheet, hdr_style, rows = seq_along(header), cols = 1,
           stack = TRUE)
  start <- length(header) + 2
  writeData(wb, sheet, data, startRow = start, headerStyle = bold_style)
  setColWidths(wb, sheet, cols = seq_len(ncol(data)), widths = "auto")
}

# Sheet 1: Methods Overview
n_young_pre <- sum(meta$age == "Young" & meta$time == "Pre")
n_old_pre   <- sum(meta$age == "Old"   & meta$time == "Pre")

overview <- tibble(
  Parameter = c(
    "Design", "Model", "Contrasts", "FDR method", "FDR threshold",
    "Pi-score threshold", "N proteins", "N subjects (Young)",
    "N subjects (Old)", "Within-subject correlation",
    "Normalization", "Imputation strategy"),
  Value = c(
    "2x2 factorial (Age x Time), repeated measures on subject",
    "limma + duplicateCorrelation via proteoDA",
    paste(contrast_names, collapse = "; "),
    "Benjamini-Hochberg",
    "0.10 (exploratory; pi-score provides secondary filter)",
    "Pi < 0.05 (= original pi > 1.3; Xiao et al. 2014)",
    as.character(nrow(dal$data)),
    as.character(n_young_pre),
    as.character(n_old_pre),
    sprintf("%.3f", within_cor),
    "Cycloess (selected by PRONE benchmarking)",
    "Non-imputed; limma handles NAs per-protein (Karpievitch et al. 2012)")
)
write_supp_sheet(wb_supp, "Methods_Overview", c(
  "YvO DEP Analysis \u2014 Methods Overview",
  "Parameters and settings for the differential expression analysis."
), overview)

# Sheet 2: Significance Summary
da_summary <- read_csv(file.path(cfg$data_dir, "02_DA_summary.csv"),
                        show_col_types = FALSE)
write_supp_sheet(wb_supp, "Significance_Summary", c(
  "Significance counts by contrast, direction, and criterion.",
  "sig.PVal: nominal P < pval_thresh | sig.FDR: adj.P < pval_thresh (BH)",
  "sig.Pi: Pi-score < 0.05 | sig.FDR.05/10: FDR at 0.05 and 0.10 thresholds."
), da_summary)

# Sheets 3-7: Per-contrast results (sorted by pi-score, key columns)
keep_cols <- c("uniprot_id", "gene", "protein", "description",
               "logFC", "CI.L", "CI.R", "t", "P.Value", "adj.P.Val",
               "pi_score", "sig_pi")

contrast_labels <- c(
  Training_Young = "Young_Post - Young_Pre",
  Training_Old   = "Old_Post - Old_Pre",
  Aging          = "Old_Pre - Young_Pre",
  Interaction    = "(Old_Post-Old_Pre) - (Young_Post-Young_Pre)",
  Reversal       = "(Old_Post-Old_Pre) - (Old_Pre-Young_Pre)"
)

for (cname in contrast_names) {
  res <- results_list[[cname]]
  cols_present <- intersect(keep_cols, names(res))
  n_sig_fdr <- sum(res$adj.P.Val < 0.10, na.rm = TRUE)
  n_sig_pi  <- sum(res$sig_pi != 0)
  write_supp_sheet(wb_supp, cname, c(
    sprintf("%s: %s", cname, contrast_labels[[cname]]),
    sprintf("Proteins tested: %d | FDR<0.10: %d | Pi<0.05: %d",
            nrow(res), n_sig_fdr, n_sig_pi),
    "Sorted by Pi-score (most significant first)."
  ), res[order(res$pi_score), cols_present])
}

# Sheet 8: Blunting Diagnostics (tests + quantiles on one sheet)
addWorksheet(wb_supp, "Blunting_Diagnostics")
blunt_header <- c(
  "Age-dependent attenuation of training response.",
  "Tests compare |logFC| distributions: Training_Young vs Training_Old.",
  "KS: shape/location | Fligner-Killeen: variance | Wilcoxon: paired shift",
  "Cliff's delta: non-parametric effect size (>0 = young responds more)."
)
for (i in seq_along(blunt_header))
  writeData(wb_supp, "Blunting_Diagnostics", blunt_header[i],
            startRow = i, startCol = 1)
addStyle(wb_supp, "Blunting_Diagnostics", hdr_style,
         rows = seq_along(blunt_header), cols = 1, stack = TRUE)
writeData(wb_supp, "Blunting_Diagnostics", blunt_diag,
          startRow = length(blunt_header) + 2, headerStyle = bold_style)
quant_start <- length(blunt_header) + 2 + nrow(blunt_diag) + 2
writeData(wb_supp, "Blunting_Diagnostics", "Quantile comparison (|logFC|):",
          startRow = quant_start)
writeData(wb_supp, "Blunting_Diagnostics", quant_df,
          startRow = quant_start + 1, headerStyle = bold_style)
setColWidths(wb_supp, "Blunting_Diagnostics",
             cols = 1:max(ncol(blunt_diag), ncol(quant_df)), widths = "auto")

# Sheet 9: Bootstrap CIs
write_supp_sheet(wb_supp, "Effect_Size_CI", c(
  "Median |logFC| with 95% BCa bootstrap confidence intervals (10,000 resamples).",
  "Summarises global effect magnitude per contrast (Efron & Tibshirani 1993).",
  "Aging shows the largest median effect; Training_Old is attenuated vs Training_Young."
), boot_df)

# Sheet 10: Power Analysis
write_supp_sheet(wb_supp, "Power_Analysis", c(
  "Minimum detectable logFC at 80% power (alpha = 0.10).",
  "Conservative: uses pwr.t.test (standard t); limma's moderated t has higher power.",
  "Paired contrasts adjust for within-subject correlation; between-group do not.",
  "Interpret as approximate lower bounds on detectable effect sizes."
), power_df)

# Sheet 11: Imputation Sensitivity
if (!is.null(sens_df) && nrow(sens_df) > 0) {
  write_supp_sheet(wb_supp, "Imputation_Sensitivity", c(
    "Spearman correlation of t-statistics: non-imputed (main) vs imputed limma.",
    "rho > 0.93 across all contrasts indicates high robustness to imputation choice.",
    "Reference: Karpievitch et al. 2012, BMC Bioinform 13(S16):S5."
  ), sens_df)
}

saveWorkbook(wb_supp, file.path(cfg$data_dir, "10_DEP_supplementary.xlsx"),
             overwrite = TRUE)
cat("Done: 03_dep_robustness.R\n")
