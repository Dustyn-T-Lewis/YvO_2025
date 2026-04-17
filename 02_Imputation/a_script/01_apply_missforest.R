#!/usr/bin/env Rscript
# 01_apply_missforest.R — missForest imputation with MAR/MNAR audit
#
# Loads cycloess-normalized matrix, classifies missingness as MAR/MNAR for
# downstream auditing, then applies missForest. The classification is used
# to flag and audit MNAR proteins, not to drive class-specific imputation.
#
# Method provenance + benchmark numbers: see 02_Imputation/README.md
#
# Outputs (c_data/):
#   01_imputed.csv               — imputed matrix with annotation
#   01_DAList_imputed.rds        — DAList with missingness annotations
#   02_mar_mnar_classification.csv
#   03_imputation_mask.csv
#   04_mnar_imputation_audit.csv
#   05_imputation_summary.txt
#   00_report_intermediates.rds

library(missForest)
library(msImpute)
library(dplyr)
library(tidyr)
library(tibble)
library(readr)
library(stringr)

setwd(rprojroot::find_rstudio_root_file())

# --- Configuration
cfg <- list(
  NORM_CSV        = "01_normalization/c_data/02_normalized.csv",
  NORM_RDS        = "01_normalization/c_data/03_DAList_normalized.rds",
  DATA_DIR        = "02_Imputation/c_data",
  MISS_UNRELIABLE = 50
)
dir.create(cfg$DATA_DIR, showWarnings = FALSE, recursive = TRUE)

# --- Palettes
PAL_GT    <- c(Young_Pre = scales::alpha("#4393C3", 0.5), Young_Post = "#4393C3",
               Old_Pre   = scales::alpha("#D6604D", 0.5), Old_Post   = "#D6604D")
PAL_MAR   <- c(MAR = "#4393C3", MNAR = "#D6604D")
PAL_CLASS <- c(Complete = "#4DAF4A", MAR = "#4393C3", MNAR = "#D6604D")

# --- 1. Load data
df  <- read_csv(cfg$NORM_CSV, show_col_types = FALSE)
ann <- df |> select(uniprot_id, gene, protein, description)
mat <- as.matrix(df[, -(1:4)])
rownames(mat) <- df$gene

if (any(duplicated(df$gene))) {
  warning("Duplicate gene names; using uniprot_id as rownames")
  rownames(mat) <- df$uniprot_id
}

dal_norm <- readRDS(cfg$NORM_RDS)
meta <- tibble(
  Col_ID     = dal_norm$metadata$Col_ID,
  Group      = dal_norm$metadata$Group,
  Timepoint  = dal_norm$metadata$Timepoint,
  Group_Time = dal_norm$metadata$Group_Time
)
stopifnot(setequal(meta$Col_ID, colnames(mat)))
cat(sprintf("Loaded: %d proteins x %d samples\n", nrow(mat), ncol(mat)))

# --- 2. Missingness profiling
prot_miss <- rowSums(is.na(mat))
prot_pct  <- prot_miss / ncol(mat) * 100
obs_means <- rowMeans(mat, na.rm = TRUE)
pct_miss  <- round(sum(is.na(mat)) / length(mat) * 100, 2)
cat(sprintf("Missing: %d / %d (%.2f%%)\n", sum(is.na(mat)), length(mat), pct_miss))

miss_by_group <- sapply(unique(meta$Group_Time), function(g) {
  cols <- meta$Col_ID[meta$Group_Time == g]
  rowSums(is.na(mat[, cols, drop = FALSE])) / length(cols) * 100
})

# --- 3. MAR/MNAR classification
has_na <- which(prot_miss > 0 & prot_miss < ncol(mat))
miss_class <- tibble(gene = rownames(mat), n_miss = prot_miss,
                     pct_miss = prot_pct, mean_intensity = obs_means)

mar_result <- tryCatch({
  feat <- msImpute::selectFeatures(mat[has_na, ], method = "ebm",
                                   group = meta$Group_Time)
  mar_names <- feat$name[feat$msImpute_feature == TRUE]
  n_incomplete <- nrow(feat)
  if (length(mar_names) < 0.05 * n_incomplete) {
    cat("  EBM degenerate (<5% MAR) -- falling back to k-means\n")
    NULL
  } else {
    cat(sprintf("EBM: %d MAR / %d MNAR\n", length(mar_names), n_incomplete - length(mar_names)))
    list(mar_genes = mar_names, method = "msImpute_ebm")
  }
}, error = function(e) {
  cat(sprintf("EBM failed (%s) -- falling back to k-means\n", conditionMessage(e)))
  NULL
})

if (is.null(mar_result)) {
  set.seed(42)
  mc_sub <- miss_class |> filter(n_miss > 0, n_miss < ncol(mat))
  km <- kmeans(scale(cbind(mc_sub$mean_intensity, mc_sub$pct_miss)),
               centers = 2, nstart = 25)
  cl_means <- tapply(mc_sub$mean_intensity, km$cluster, mean)
  mnar_cl <- which.min(cl_means)
  mar_result <- list(
    mar_genes = mc_sub$gene[km$cluster != mnar_cl],
    method = sprintf("k-means (cluster means: %.1f vs %.1f)",
                     cl_means[mnar_cl], cl_means[-mnar_cl]))
}

miss_class <- miss_class |>
  mutate(
    classification = case_when(
      n_miss == 0 ~ "Complete",
      gene %in% mar_result$mar_genes ~ "MAR",
      TRUE ~ "MNAR"),
    imputation_reliable = classification == "Complete" | pct_miss < cfg$MISS_UNRELIABLE
  )

mnar_genes <- miss_class$gene[miss_class$classification == "MNAR"]

# Group-stratified missingness (Fisher test per MNAR protein)
group_miss_pval <- setNames(rep(NA_real_, nrow(miss_class)), miss_class$gene)
for (g in mnar_genes) {
  ct <- sapply(unique(meta$Group_Time), function(gt) {
    cols <- meta$Col_ID[meta$Group_Time == gt]
    c(missing = sum(is.na(mat[g, cols])), observed = sum(!is.na(mat[g, cols])))
  })
  group_miss_pval[g] <- tryCatch(
    fisher.test(ct, simulate.p.value = TRUE, B = 2000)$p.value,
    error = function(e) NA_real_)
}
miss_class$group_miss_pval <- group_miss_pval[miss_class$gene]

n_mar_prots  <- sum(miss_class$classification == "MAR")
n_mnar_prots <- length(mnar_genes)
n_comp_prots <- sum(miss_class$classification == "Complete")
mar_miss_vals  <- sum(miss_class$n_miss[miss_class$classification == "MAR"])
mnar_miss_vals <- sum(miss_class$n_miss[miss_class$classification == "MNAR"])
total_miss_vals <- mar_miss_vals + mnar_miss_vals

cat(sprintf("Classification (%s): MAR %d | MNAR %d | Complete %d\n",
            mar_result$method, n_mar_prots, n_mnar_prots, n_comp_prots))

# --- 4. Apply missForest
# Note: missForest treats all missing values as MAR (ignorable mechanism).
# The MAR/MNAR classification above is retained for downstream flagging and
# reporting; the benchmark (a_script/benchmark/) validates that MAR-only
# imputation performs well on this dataset's missingness pattern.
cat("Imputing with missForest...\n")
set.seed(42)
mf_result <- missForest::missForest(t(mat), maxiter = 10, ntree = 100, verbose = TRUE)
mat_imp <- t(mf_result$ximp)
rownames(mat_imp) <- rownames(mat)
colnames(mat_imp) <- colnames(mat)
stopifnot(sum(is.na(mat_imp)) == 0)
cat(sprintf("missForest OOB error: %.4f\n", mf_result$OOBerror[1]))

# --- 5. Export
was_na <- is.na(mat)
best <- "missForest"

# MNAR audit
mnar_audit <- tibble(
  gene      = mnar_genes,
  pre_mean  = rowMeans(mat[mnar_genes, ], na.rm = TRUE),
  post_mean = rowMeans(mat_imp[mnar_genes, ]),
  pre_sd    = apply(mat[mnar_genes, ], 1, sd, na.rm = TRUE),
  pct_miss  = prot_pct[mnar_genes],
  shift     = rowMeans(mat_imp[mnar_genes, ]) - rowMeans(mat[mnar_genes, ], na.rm = TRUE),
  effect_d  = (rowMeans(mat_imp[mnar_genes, ]) - rowMeans(mat[mnar_genes, ], na.rm = TRUE)) /
              apply(mat[mnar_genes, ], 1, sd, na.rm = TRUE),
  imputation_reliable = prot_pct[mnar_genes] < cfg$MISS_UNRELIABLE
)

# 01_imputed.csv
write_csv(bind_cols(ann, as_tibble(mat_imp)),
          file.path(cfg$DATA_DIR, "01_imputed.csv"))

# 01_DAList_imputed.rds
dal <- readRDS(cfg$NORM_RDS)
mat_imp_uid <- mat_imp
rownames(mat_imp_uid) <- ann$uniprot_id
dal$data <- mat_imp_uid
n_ann <- nrow(dal$annotation)
dal$annotation <- merge(
  dal$annotation,
  miss_class |> select(gene, n_miss, pct_miss,
                       miss_classification = classification,
                       imputation_reliable),
  by = "gene", all.x = TRUE, sort = FALSE)
stopifnot(nrow(dal$annotation) == n_ann)
saveRDS(dal, file.path(cfg$DATA_DIR, "01_DAList_imputed.rds"))

# Supporting tables
write_csv(miss_class, file.path(cfg$DATA_DIR, "02_mar_mnar_classification.csv"))
write_csv(bind_cols(tibble(gene = rownames(was_na)), as_tibble(was_na)),
          file.path(cfg$DATA_DIR, "03_imputation_mask.csv"))
write_csv(mnar_audit, file.path(cfg$DATA_DIR, "04_mnar_imputation_audit.csv"))

# Summary
info <- list(
  n_proteins   = nrow(mat), n_samples = ncol(mat), pct_missing = pct_miss,
  n_complete   = n_comp_prots, n_mar_proteins = n_mar_prots,
  n_mnar_proteins = n_mnar_prots, n_mar_values = mar_miss_vals,
  n_mnar_values = mnar_miss_vals,
  pct_mar_values = round(mar_miss_vals / total_miss_vals * 100, 1),
  classification_method = mar_result$method,
  n_unreliable = sum(!miss_class$imputation_reliable),
  best_method  = best,
  oob_error    = as.numeric(round(mf_result$OOBerror[1], 4)))
writeLines(paste(names(info), info, sep = " = "),
           file.path(cfg$DATA_DIR, "05_imputation_summary.txt"))

# Report intermediates
saveRDS(list(
  mat = mat, mat_imp = mat_imp, was_na = was_na, ann = ann, meta = meta,
  miss_class = miss_class, miss_by_group = miss_by_group,
  prot_pct = prot_pct, pct_miss = pct_miss,
  mnar_genes = mnar_genes, mnar_audit = mnar_audit,
  best = best, n_mar_prots = n_mar_prots, n_mnar_prots = n_mnar_prots,
  mar_miss_vals = mar_miss_vals, mnar_miss_vals = mnar_miss_vals,
  total_miss_vals = total_miss_vals,
  classification_method = mar_result$method,
  oob_error = as.numeric(mf_result$OOBerror[1]),
  PAL_GT = PAL_GT, PAL_MAR = PAL_MAR, PAL_CLASS = PAL_CLASS
), file.path(cfg$DATA_DIR, "00_report_intermediates.rds"))

cat(sprintf("Done: missForest | %d proteins x %d samples | OOB=%.4f\n",
            nrow(mat_imp), ncol(mat_imp), mf_result$OOBerror[1]))

writeLines(capture.output(sessionInfo()), file.path(cfg$DATA_DIR, "sessionInfo.txt"))
