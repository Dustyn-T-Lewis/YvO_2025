#!/usr/bin/env Rscript
# YvO Normalization — DIA-MS skeletal muscle proteomics (Young vs Old × Pre/Post)
#
# Pipeline: HPA tissue filter → dedup → 66% group missingness →
#           consensus outlier detection → cycloess normalization (proteoDA)
#
# Outputs (c_data/): 00-03 DAList/CSV artifacts
# Reports (b_reports/): 01_norm_comparison, 02_qc_pre, 03_qc_post (.pdf)
#
# Refs: Thurman 2023 (proteoDA), Bolstad 2003 (cycloess),
#       Brenes 2024 (CV on linear scale), Huang 2024 (SEAOP outlier)

library(proteoDA)
library(readxl)
library(readr)
library(dplyr)
library(tidyr)
library(stringr)

set.seed(42)
setwd(rprojroot::find_rstudio_root_file())

# --- Configuration -----------------------------------------------------------

cfg <- list(
  # Input files
  raw_file  = "00_input/YvO_raw.xlsx",
  meta_file = "00_input/YvO_meta.xlsx",
  hpa_file  = "00_input/HPA_skeletal_muscle_annotations.tsv",

  # Output directories
  report_dir = "01_normalization/b_reports",
  data_dir   = "01_normalization/c_data",

  # Thresholds
  min_reps    = 10,          # min non-missing samples in ≥1 Group_Time
  outlier_k   = 3,          # methods that must agree for consensus
  mahal_p     = 0.01,       # PCA Mahalanobis chi-sq cutoff
  mad_k       = 3,          # MAD multiplier for median intensity
  norm_method = "cycloess"
)

dir.create(cfg$report_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(cfg$data_dir,   recursive = TRUE, showWarnings = FALSE)

# --- Helper ------------------------------------------------------------------

run_pca <- function(mat, metadata, log_transform = TRUE) {
  # Median-impute for PCA only (imputed values never exported)
  for (j in seq_len(ncol(mat)))
    mat[is.na(mat[, j]), j] <- median(mat[, j], na.rm = TRUE)
  if (log_transform) mat <- log2(mat + 1)
  pca <- prcomp(t(mat), center = TRUE, scale. = TRUE)
  ve  <- round(summary(pca)$importance[2, 1:3] * 100, 1)
  pc  <- as.data.frame(pca$x[, 1:3]) |>
    mutate(Col_ID = rownames(pca$x)) |>
    left_join(metadata, by = "Col_ID")
  list(pca = pca, scores = pc, var_exp = ve)
}

# =============================================================================
# 1. LOAD DATA
# =============================================================================

raw <- read_excel(cfg$raw_file)
annot_cols <- c("uniprot_id", "protein", "gene", "description", "n_seq")
annotation <- raw[, annot_cols]
intensity  <- raw[, setdiff(names(raw), annot_cols)]

metadata <- as.data.frame(read_excel(cfg$meta_file))
rownames(metadata) <- metadata$Col_ID
stopifnot("Sample mismatch" = setequal(colnames(intensity), metadata$Col_ID))
intensity <- intensity[, metadata$Col_ID]

n_raw <- nrow(annotation)
filter_log <- tibble(step = "Raw input", n_before = NA_integer_,
                     n_after = n_raw, n_removed = NA_integer_)
cat(sprintf("Raw: %d proteins x %d samples\n", n_raw, ncol(intensity)))

# =============================================================================
# 2. HPA TISSUE FILTER
# =============================================================================

hpa <- read_tsv(cfg$hpa_file, show_col_types = FALSE) |>
  select(Gene, Ensembl, Evidence,
         Protein_class    = `Protein class`,
         Subcellular_main = `Subcellular main location`,
         Interactions) |>
  distinct(Gene, .keep_all = TRUE)

n_before <- nrow(annotation)
keep_hpa   <- annotation$gene %in% hpa$Gene
intensity  <- intensity[keep_hpa, ]
annotation <- annotation[keep_hpa, ] |> left_join(hpa, by = c("gene" = "Gene"))
removed_genes <- setdiff(raw$gene, annotation$gene)

filter_log <- bind_rows(filter_log, tibble(
  step = "HPA tissue filter", n_before = n_before,
  n_after = nrow(annotation), n_removed = n_before - nrow(annotation)))
cat(sprintf("HPA: %d -> %d (-%d)\n",
            n_before, nrow(annotation), n_before - nrow(annotation)))

# =============================================================================
# 3. DEDUPLICATE BY UNIPROT ID
# =============================================================================

if (any(duplicated(annotation$uniprot_id))) {
  n_before_dup <- nrow(annotation)
  annotation$row_mean <- rowMeans(data.matrix(intensity), na.rm = TRUE)
  keep_idx <- annotation |>
    mutate(row_idx = row_number()) |>
    group_by(uniprot_id) |>
    slice_max(row_mean, n = 1, with_ties = FALSE) |>
    pull(row_idx)
  annotation <- annotation[keep_idx, ]
  intensity  <- intensity[keep_idx, ]
  annotation$row_mean <- NULL
  filter_log <- bind_rows(filter_log, tibble(
    step = "Deduplication", n_before = n_before_dup,
    n_after = nrow(annotation), n_removed = n_before_dup - nrow(annotation)))
  cat(sprintf("Deduplicated: %d proteins\n", nrow(annotation)))
}

# =============================================================================
# 4. ASSEMBLE DAList & MISSINGNESS FILTER
# =============================================================================

int_mat <- as.data.frame(data.matrix(intensity))
rownames(int_mat) <- annotation$uniprot_id
annot_df <- as.data.frame(annotation); rownames(annot_df) <- annotation$uniprot_id
meta_df  <- as.data.frame(metadata);   rownames(meta_df)  <- metadata$Col_ID

dal <- DAList(data = int_mat, annotation = annot_df, metadata = meta_df)
dal <- zero_to_missing(dal)

n_before <- nrow(dal$data)
dal <- filter_proteins_by_group(dal, min_reps = cfg$min_reps, min_groups = 1,
                                grouping_column = "Group_Time")

filter_log <- bind_rows(filter_log, tibble(
  step = sprintf("Missingness (>=%d reps in >=1 group)", cfg$min_reps),
  n_before = n_before, n_after = nrow(dal$data),
  n_removed = n_before - nrow(dal$data)))
filter_log <- filter_log |> mutate(pct_of_raw = round(n_after / n_raw * 100, 1))
cat(sprintf("Missingness: %d -> %d (-%d)\n",
            n_before, nrow(dal$data), n_before - nrow(dal$data)))

filtered_proteins <- bind_rows(
  tibble(uniprot_id = raw$uniprot_id, gene = raw$gene,
         description = raw$description) |>
    filter(gene %in% removed_genes) |>
    mutate(removal_step = "HPA tissue filter"),
  annot_df |>
    filter(!uniprot_id %in% rownames(dal$data)) |>
    select(uniprot_id, gene, description) |>
    mutate(removal_step = sprintf("Missingness (<%d reps in all groups)", cfg$min_reps))
) |> distinct(uniprot_id, .keep_all = TRUE)

# =============================================================================
# 5. OUTLIER DETECTION (4-method consensus, flag if ≥3/4)
# =============================================================================

# Method 1: Sample missingness (pooled threshold, flags individual samples)
pct_missing <- colMeans(is.na(dal$data)) * 100

miss_info <- dal$metadata |>
  select(Col_ID, Subject_ID, Group, Timepoint) |>
  mutate(pct_missing = pct_missing[Col_ID],
         prefix = str_remove(Col_ID, "_(Pre|Post)$")) |>
  group_by(prefix) |>
  mutate(delta_missing = abs(diff(pct_missing))) |>
  ungroup()

miss_thresh  <- quantile(pct_missing, 0.75) + 1.5 * IQR(pct_missing)
delta_thresh <- quantile(miss_info$delta_missing, 0.75) +
  1.5 * IQR(miss_info$delta_missing)

miss_info$miss_flag <- miss_info$pct_missing > miss_thresh |
  miss_info$delta_missing > delta_thresh

# Method 2: PCA Mahalanobis distance
complete_mat <- dal$data[rowSums(is.na(dal$data)) == 0, ]
pca_pre  <- run_pca(complete_mat, dal$metadata, log_transform = TRUE)
pc3      <- pca_pre$pca$x[, 1:3]
mahal    <- mahalanobis(pc3, colMeans(pc3), cov(pc3))
pca_flags <- tibble(Col_ID = colnames(dal$data), mahal_dist = mahal,
                    pca_flag = mahal > qchisq(1 - cfg$mahal_p, df = 3))

# Method 3: MAD-based median intensity
samp_med   <- apply(log2(dal$data + 1), 2, median, na.rm = TRUE)
global_med <- median(samp_med)
mad_val    <- mad(samp_med)
mad_flags  <- tibble(Col_ID = names(samp_med), sample_median = samp_med,
                     mad_flag = abs(samp_med - global_med) > cfg$mad_k * mad_val)

# Method 4: Inter-sample correlation (median pairwise Pearson r)
cor_mat     <- cor(log2(dal$data + 1), use = "pairwise.complete.obs")
med_cor     <- apply(cor_mat, 2, function(x) median(x[x < 1], na.rm = TRUE))
cor_med_all <- median(med_cor)
cor_mad_all <- mad(med_cor)
cor_flags   <- tibble(Col_ID = names(med_cor), median_cor = med_cor,
                      cor_flag = med_cor < cor_med_all - cfg$mad_k * cor_mad_all)

# Consensus (≥3 of 4 methods)
outlier_diag <- miss_info |>
  left_join(pca_flags, by = "Col_ID") |>
  left_join(mad_flags, by = "Col_ID") |>
  left_join(cor_flags, by = "Col_ID") |>
  mutate(n_flags = miss_flag + pca_flag + mad_flag + cor_flag,
         consensus_outlier = n_flags >= cfg$outlier_k)

n_outliers <- sum(outlier_diag$consensus_outlier)
cat(sprintf("Outliers: %d sample(s) flagged (%d/4 consensus)\n",
            n_outliers, cfg$outlier_k))

# Snapshot pre-outlier state for diagnostics (detection plot includes outliers)
data_pre_outlier <- dal$data
meta_pre_outlier <- dal$metadata

# Remove outlier samples (individual timepoints, not paired)
outlier_ids <- outlier_diag |> filter(consensus_outlier) |> pull(Col_ID)
if (n_outliers > 0) {
  dal <- filter_samples(dal, !(Col_ID %in% outlier_ids))
  cat(sprintf("Removed: %s (%d remain)\n", paste(outlier_ids, collapse = ", "), ncol(dal$data)))
}

# =============================================================================
# 6. NORMALIZE (cycloess via proteoDA)
# =============================================================================

write_norm_report(dal, grouping_column = "Group_Time",
                  output_dir = cfg$report_dir,
                  filename = "01_norm_comparison.pdf", overwrite = TRUE)

write_qc_report(dal, color_column = "Group_Time",
                output_dir = cfg$report_dir,
                filename = "02_qc_pre.pdf", overwrite = TRUE)

dal_pre <- dal
saveRDS(dal_pre, file.path(cfg$data_dir, "01_DAList_prenorm.rds"))

dal <- normalize_data(dal, norm_method = cfg$norm_method)
cat(sprintf("Normalized (%s): %d proteins x %d samples\n",
            cfg$norm_method, nrow(dal$data), ncol(dal$data)))

write_qc_report(dal, color_column = "Group_Time",
                output_dir = cfg$report_dir,
                filename = "03_qc_post.pdf", overwrite = TRUE)

# =============================================================================
# 7. EXPORT
# =============================================================================

export_df <- bind_cols(
  as_tibble(dal$annotation) |> select(uniprot_id, protein, gene, description),
  as_tibble(dal$data))

write_csv(export_df, file.path(cfg$data_dir, "02_normalized.csv"))
saveRDS(dal, file.path(cfg$data_dir, "03_DAList_normalized.rds"))

# =============================================================================
# 8. COMPUTE PLOT DATA & SAVE INTERMEDIATES for 02_norm_reports.R
# =============================================================================

filter_bar_data <- filter_log |>
  filter(!is.na(n_removed)) |>
  mutate(step = factor(step, levels = step)) |>
  pivot_longer(c(n_after, n_removed), names_to = "status", values_to = "n") |>
  mutate(status = recode(status, n_after = "Retained", n_removed = "Removed"))

miss_bar_data <- meta_pre_outlier |>
  select(Col_ID, Group_Time) |>
  mutate(detected = colSums(!is.na(data_pre_outlier[, Col_ID])),
         missing  = nrow(data_pre_outlier) - detected,
         is_outlier = Col_ID %in% outlier_ids) |>
  pivot_longer(c(detected, missing), names_to = "status", values_to = "n") |>
  mutate(status = str_to_title(status))

subj_var <- dal$metadata |>
  mutate(iqr = apply(log2(dal$data[, Col_ID] + 1), 2, IQR, na.rm = TRUE),
         Subject_ID = str_remove(Col_ID, "_(Pre|Post)$")) |>
  select(Col_ID, Subject_ID, Group, Timepoint, Group_Time, iqr)

log_dat <- log2(dal$data + 1)
grp_vec <- dal$metadata$Group_Time[match(colnames(log_dat), dal$metadata$Col_ID)]
eta2_vals <- apply(log_dat, 1, function(x) {
  ok <- !is.na(x)
  if (sum(ok) < 4) return(NA_real_)
  xk <- x[ok]; gk <- grp_vec[ok]
  ss_b <- sum(tapply(xk, gk, length) * (tapply(xk, gk, mean) - mean(xk))^2)
  ss_t <- sum((xk - mean(xk))^2)
  if (ss_t > 0) ss_b / ss_t else NA_real_
})

pca_post <- run_pca(dal$data, dal$metadata, log_transform = FALSE)

intermediates <- list(
  cfg              = cfg,
  filter_log       = filter_log,
  filter_bar_data  = filter_bar_data,
  miss_bar_data    = miss_bar_data,
  n_raw            = n_raw,
  n_outliers       = n_outliers,
  outlier_diag     = outlier_diag,
  outlier_ids      = outlier_ids,
  miss_thresh      = miss_thresh,
  delta_thresh     = delta_thresh,
  pca_pre          = pca_pre,
  pca_post         = pca_post,
  global_med       = global_med,
  mad_val          = mad_val,
  subj_var         = subj_var,
  eta2_vals        = eta2_vals,
  filtered_proteins = filtered_proteins,
  data_pre_outlier = data_pre_outlier,
  meta_pre_outlier = meta_pre_outlier,
  dal_nrow         = nrow(dal$data),
  dal_ncol         = ncol(dal$data)
)

saveRDS(intermediates, file.path(cfg$data_dir, "00_report_intermediates.rds"))

cat(sprintf("Done: %d proteins x %d samples -> %s/\n",
            nrow(dal$data), ncol(dal$data), cfg$data_dir))
