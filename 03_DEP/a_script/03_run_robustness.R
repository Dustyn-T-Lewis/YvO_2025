#!/usr/bin/env Rscript
# Stage 03: Robustness analyses
# Blunting, bootstrap CIs, power analysis, imputation sensitivity
# Adds sheets to 03_DEP_results.xlsx

setwd(rprojroot::find_root(rprojroot::has_file("setup.R")))

library(dplyr)
library(tibble)
library(purrr)
library(readxl)
library(openxlsx)
library(boot)
library(pwr)
library(proteoDA)

set.seed(42)

DAT  <- "03_DEP/c_data"
XLSX <- file.path(DAT, "03_DEP_results.xlsx")

dal <- readRDS(file.path(DAT, "01_limma_DAList.rds"))
contrast_names <- names(dal$results)
meta <- as.data.frame(dal$metadata)

results_list <- lapply(contrast_names, \(cn) {
  as.data.frame(read_excel(XLSX, sheet = cn))
})
names(results_list) <- contrast_names

# Blunting diagnostics

blunt_df <- tibble(
  gene      = results_list[["Training_Young"]]$gene,
  abs_lfc_y = abs(results_list[["Training_Young"]]$logFC),
  abs_lfc_o = abs(results_list[["Training_Old"]]$logFC)
) |> filter(!is.na(abs_lfc_y) & !is.na(abs_lfc_o))

ks_res <- ks.test(blunt_df$abs_lfc_y, blunt_df$abs_lfc_o)
fk_res <- fligner.test(
  abs_lfc ~ contrast,
  data = data.frame(
    abs_lfc  = c(blunt_df$abs_lfc_y, blunt_df$abs_lfc_o),
    contrast = rep(c("Young", "Old"), each = nrow(blunt_df))))
wx_res <- wilcox.test(blunt_df$abs_lfc_y, blunt_df$abs_lfc_o, paired = TRUE)
diffs  <- blunt_df$abs_lfc_y - blunt_df$abs_lfc_o
cliff  <- (sum(diffs > 0) - sum(diffs < 0)) / length(diffs)
cliff_mag <- case_when(
  abs(cliff) < 0.147 ~ "negligible", abs(cliff) < 0.33 ~ "small",
  abs(cliff) < 0.474 ~ "medium",     TRUE ~ "large")

blunt_diag <- tibble(
  test = c("KS", "Fligner-Killeen", "Wilcoxon signed-rank", "Cliff's delta"),
  statistic = c(ks_res$statistic, fk_res$statistic, wx_res$statistic, cliff),
  p_value = c(ks_res$p.value, fk_res$p.value, wx_res$p.value, NA),
  interpretation = c(
    ifelse(ks_res$p.value < 0.05, "Distributions differ", "No difference"),
    ifelse(fk_res$p.value < 0.05, "Variance differs", "No difference"),
    ifelse(wx_res$p.value < 0.05, "Paired shift significant", "No shift"),
    sprintf("%s (d=%.3f; >0 = young responds more)",
            stringr::str_to_title(cliff_mag), cliff)))

message(sprintf("Blunting: KS p=%.2g, Cliff d=%.3f (%s)",
                ks_res$p.value, cliff, cliff_mag))

# Bootstrap CI (median |logFC|, BCa, 10k reps)

boot_df <- list_rbind(lapply(contrast_names, \(cname) {
  vals <- abs(results_list[[cname]]$logFC)
  vals <- vals[!is.na(vals)]
  b  <- boot(vals, \(d, i) median(d[i]), R = 10000)
  ci <- tryCatch(boot.ci(b, type = "bca"),
                 error = \(e) boot.ci(b, type = "perc"))
  ci_lo <- if (!is.null(ci$bca)) ci$bca[4] else ci$percent[4]
  ci_hi <- if (!is.null(ci$bca)) ci$bca[5] else ci$percent[5]
  tibble(contrast = cname, median_absLFC = median(vals),
         ci_lower = ci_lo, ci_upper = ci_hi,
         boot_se = sd(b$t), n_proteins = length(vals))
}))

# Power analysis (min detectable logFC at 80% power)

fit <- dal$eBayes_fit
within_cor <- fit$correlation %||% dal$tags$duplicate_correlation %||% NA_real_
sigma_res  <- sqrt(mean(fit$sigma^2, na.rm = TRUE))
n_young    <- sum(meta$age == "Young" & meta$time == "Pre")
n_old      <- sum(meta$age == "Old"   & meta$time == "Pre")

power_df <- list_rbind(lapply(contrast_names, \(cname) {
  n_subj <- switch(cname,
    Training_Young = n_young, Training_Old = n_old, min(n_young, n_old))
  paired  <- cname %in% c("Training_Young", "Training_Old")
  eff_sig <- if (paired && !is.na(within_cor))
    sigma_res * sqrt(2 * (1 - within_cor)) else sigma_res * sqrt(2)
  pw <- pwr.t.test(n = n_subj, d = NULL, sig.level = 0.10, power = 0.80,
                   type = if (paired) "paired" else "two.sample")
  tibble(contrast = cname, n_subjects = n_subj,
         within_cor = ifelse(paired, within_cor, NA_real_),
         effective_sigma = round(eff_sig, 4),
         min_detectable_d = round(pw$d, 4),
         min_detectable_logFC = round(pw$d * eff_sig, 4),
         power = 0.80, alpha = 0.10)
}))

# Imputation sensitivity

IMP_RDS <- "02_Imputation/c_data/01_DAList_imputed.rds"
sens_df <- tibble(contrast = character(), spearman_rho = numeric(),
                  p_value = numeric(), n_proteins = integer())

if (file.exists(IMP_RDS)) {
  dal_imp_raw <- readRDS(IMP_RDS)
  imp_mat <- as.matrix(dal_imp_raw$data)

  # Align imp_mat rows to the DEP annotation row order. The imputed RDS row
  # order (gene_order alphabetical) differs from the DEP DAList annotation
  # order (normalization order). DAList() silently re-labels rownames(data)
  # to match annotation while leaving values in place — without this match()
  # protein labels detach from intensities.
  ann_dep <- as.data.frame(dal$annotation)
  ord <- match(ann_dep$uniprot_id, rownames(imp_mat))
  if (any(is.na(ord))) {
    stop("Imputed DAList is missing proteins from the DEP annotation.")
  }
  imp_mat <- imp_mat[ord, , drop = FALSE]
  stopifnot(identical(rownames(imp_mat), ann_dep$uniprot_id))

  shared_samps <- intersect(colnames(dal$data), colnames(imp_mat))

  dal_imp <- DAList(
    data       = imp_mat[, shared_samps],
    annotation = ann_dep,
    metadata   = meta[meta$sample_id %in% shared_samps, ],
    tags       = list(norm_method = "cycloess_imputed"))
  dal_imp <- add_design(dal_imp, "~ 0 + group + (1 | subject)")
  colnames(dal_imp$design$design_matrix) <- gsub("^group", "",
    colnames(dal_imp$design$design_matrix))
  dal_imp <- add_contrasts(dal_imp, contrasts_vector = c(
    "Training_Young = Young_Post - Young_Pre",
    "Training_Old   = Old_Post - Old_Pre",
    "Aging          = Old_Pre - Young_Pre",
    "Interaction    = (Old_Post - Old_Pre) - (Young_Post - Young_Pre)"))
  dal_imp <- fit_limma_model(dal_imp)
  dal_imp <- extract_DA_results(dal_imp, pval_thresh = 0.10,
                                 lfc_thresh = 0, adj_method = "BH")

  comb <- readr::read_csv(file.path(DAT, "03_combined_results.csv"),
                           show_col_types = FALSE)

  sens_df <- list_rbind(lapply(contrast_names, \(cname) {
    t_col <- paste0("t_", cname)
    if (!(t_col %in% names(comb)) || !(cname %in% names(dal_imp$results))) {
      return(NULL)
    }
    imp_t <- dal_imp$results[[cname]] |>
      rownames_to_column("uniprot_id") |>
      select(uniprot_id, t_imp = t)
    merged <- inner_join(
      comb |> select(uniprot_id, t_nonimp = all_of(t_col)),
      imp_t, by = join_by(uniprot_id)) |>
      filter(!is.na(t_nonimp) & !is.na(t_imp))
    sp <- cor.test(merged$t_nonimp, merged$t_imp, method = "spearman")
    tibble(contrast = cname, spearman_rho = round(sp$estimate, 4),
           p_value = sp$p.value, n_proteins = nrow(merged))
  }))
  message("Imputation sensitivity: ", paste(sprintf("%s rho=%.3f",
    sens_df$contrast, sens_df$spearman_rho), collapse = " | "))
} else {
  message("Imputed DAList not found — skipping sensitivity")
}

# Add robustness sheets to xlsx

write_sheet <- function(wb, name, data) {
  addWorksheet(wb, name)
  writeData(wb, name, data,
    headerStyle = createStyle(textDecoration = "bold", fgFill = "#DCE6F1"))
  freezePane(wb, name, firstRow = TRUE)
  setColWidths(wb, name, cols = seq_len(ncol(data)), widths = "auto")
}

wb <- loadWorkbook(XLSX)
robustness_sheets <- c("blunting", "bootstrap_ci", "power_analysis", "imputation_sensitivity")
for (s in intersect(robustness_sheets, names(wb))) removeWorksheet(wb, s)
write_sheet(wb, "blunting",               blunt_diag)
write_sheet(wb, "bootstrap_ci",           boot_df)
write_sheet(wb, "power_analysis",         power_df)
if (nrow(sens_df) > 0) {
  write_sheet(wb, "imputation_sensitivity", sens_df)
}
saveWorkbook(wb, XLSX, overwrite = TRUE)

# Box copy

BOX <- Sys.getenv("YVO_BOX_SUPP", unset = "")
if (nzchar(BOX) && dir.exists(BOX)) {
  box_tbl <- file.path(BOX, "tables")
  dir.create(box_tbl, recursive = TRUE, showWarnings = FALSE)
  file.copy(XLSX, file.path(box_tbl, "S03_Table_DEP.xlsx"), overwrite = TRUE)
  message("Copied to Box: tables/S03_Table_DEP.xlsx")
}

message("Done: robustness analyses added to ", basename(XLSX))
