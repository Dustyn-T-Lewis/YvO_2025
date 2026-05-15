#!/usr/bin/env Rscript
# Stage 01: Normalize
# HPA filter → blood removal → dedup → missingness → outlier detection → cycloess
#
# Outputs:
#   c_data/01_normalization.xlsx      8-sheet supplement
#   c_data/03_DAList_normalized.rds   proteoDA object for stages 02–04
#   c_data/00_report_intermediates.rds diagnostic data for reports + F00
#   b_reports/01–03*.pdf              proteoDA QC reports

setwd(rprojroot::find_rstudio_root_file())

library(proteoDA)
library(readxl)
library(readr)
library(dplyr)
library(tidyr)
library(stringr)
library(openxlsx)

set.seed(42)

RPT <- "01_normalization/b_reports"
DAT <- "01_normalization/c_data"
dir.create(RPT, recursive = TRUE, showWarnings = FALSE)
dir.create(DAT, recursive = TRUE, showWarnings = FALSE)

MIN_REPS  <- 10L
OUTLIER_K <- 3L
MAHAL_P   <- 0.01
MAD_K     <- 3

run_pca <- function(mat, metadata, log_transform = TRUE) {
  for (j in seq_len(ncol(mat)))
    mat[is.na(mat[, j]), j] <- median(mat[, j], na.rm = TRUE)
  if (log_transform) mat <- log2(mat)
  pca <- prcomp(t(mat), center = TRUE, scale. = TRUE)
  ve  <- round(summary(pca)$importance[2, 1:3] * 100, 1)
  list(pca = pca,
       scores = as.data.frame(pca$x[, 1:3]) |>
         mutate(Col_ID = rownames(pca$x)) |>
         left_join(metadata, by = join_by(Col_ID)),
       var_exp = ve)
}

write_sheet <- function(wb, name, data) {
  addWorksheet(wb, name)
  writeData(wb, name, data,
    headerStyle = createStyle(textDecoration = "bold", fgFill = "#DCE6F1"))
  freezePane(wb, name, firstRow = TRUE)
  setColWidths(wb, name, cols = seq_len(ncol(data)), widths = "auto")
}

# 1. Load

raw <- read_excel("00_input/YvO_raw.xlsx")
annot_cols <- c("uniprot_id", "protein", "gene", "description", "n_seq")
annotation <- raw[, annot_cols]
intensity  <- raw[, setdiff(names(raw), annot_cols)]

metadata <- as.data.frame(read_excel("00_input/YvO_meta.xlsx"))
rownames(metadata) <- metadata$Col_ID
stopifnot("Sample mismatch" = setequal(colnames(intensity), metadata$Col_ID))
intensity <- intensity[, metadata$Col_ID]

n_raw <- nrow(annotation)
filter_log <- tibble(step = "Raw input", n_before = NA_integer_,
                     n_after = n_raw, n_removed = NA_integer_)
message(sprintf("Raw: %d proteins x %d samples", n_raw, ncol(intensity)))

# 2. HPA tissue filter

hpa <- read_tsv("00_input/HPA_skeletal_muscle_annotations.tsv",
                show_col_types = FALSE) |>
  select(Gene, Ensembl, Evidence,
         Protein_class    = `Protein class`,
         Subcellular_main = `Subcellular main location`,
         Interactions) |>
  distinct(Gene, .keep_all = TRUE)

n_before <- nrow(annotation)
keep_hpa <- annotation$gene %in% hpa$Gene
intensity  <- intensity[keep_hpa, ]
annotation <- annotation[keep_hpa, ] |>
  left_join(hpa, by = join_by(gene == Gene))
removed_genes <- setdiff(raw$gene, annotation$gene)

filter_log <- bind_rows(filter_log, tibble(
  step = "HPA tissue filter", n_before = n_before,
  n_after = nrow(annotation), n_removed = n_before - nrow(annotation)))

# 3. Blood contaminant removal (Geyer 2016 + HPA Ig)

BLOOD_CONTAMINANTS <- c(
  "HBA1", "HBA2", "HBB",
  "ALB", "TF", "HP", "HPX", "GC",
  "APOA1", "APOA2", "APOB", "APOC1", "APOC2", "APOC3",
  "FGA", "FGB", "FGG", "F2", "PLG",
  "C3", "C4A", "C4B", "C5", "C6", "C7", "C8A", "C8B", "C8G", "C9",
  "CFB", "CFH", "CFI", "C1QB", "C1QC", "C1R", "C1S", "C2",
  "SERPINA1", "SERPINA3", "A2M", "ORM1", "ORM2", "AHSG", "ITIH4",
  "AGT", "AMBP", "KNG1", "HRG", "VTN")
hpa_ig <- hpa$Gene[grepl("Immunoglobulin genes", hpa$Protein_class, fixed = TRUE)]
blood_genes <- unique(c(BLOOD_CONTAMINANTS, hpa_ig))

n_before <- nrow(annotation)
keep <- !annotation$gene %in% blood_genes
intensity  <- intensity[keep, ]
annotation <- annotation[keep, ]

filter_log <- bind_rows(filter_log, tibble(
  step = "Blood contaminant removal", n_before = n_before,
  n_after = nrow(annotation), n_removed = n_before - nrow(annotation)))

# 4. Deduplicate by UniProt ID (keep max mean intensity)

if (any(duplicated(annotation$uniprot_id))) {
  n_before <- nrow(annotation)
  annotation$row_mean <- rowMeans(data.matrix(intensity), na.rm = TRUE)
  keep_idx <- annotation |>
    mutate(row_idx = row_number()) |>
    slice_max(row_mean, n = 1, with_ties = FALSE, by = uniprot_id) |>
    pull(row_idx)
  annotation <- annotation[keep_idx, ]
  intensity  <- intensity[keep_idx, ]
  annotation$row_mean <- NULL
  filter_log <- bind_rows(filter_log, tibble(
    step = "Deduplication", n_before = n_before,
    n_after = nrow(annotation), n_removed = n_before - nrow(annotation)))
}

# 5. Assemble DAList + missingness filter

int_mat <- as.data.frame(data.matrix(intensity))
rownames(int_mat) <- annotation$uniprot_id
annot_df <- as.data.frame(annotation); rownames(annot_df) <- annotation$uniprot_id
meta_df  <- as.data.frame(metadata);   rownames(meta_df)  <- metadata$Col_ID

dal <- DAList(data = int_mat, annotation = annot_df, metadata = meta_df) |>
  zero_to_missing()

n_before <- nrow(dal$data)
dal <- filter_proteins_by_group(dal, min_reps = MIN_REPS, min_groups = 1,
                                grouping_column = "Group_Time")

filter_log <- bind_rows(filter_log, tibble(
  step = sprintf("Missingness (>=%d in >=1 group)", MIN_REPS),
  n_before = n_before, n_after = nrow(dal$data),
  n_removed = n_before - nrow(dal$data))) |>
  mutate(pct_of_raw = round(n_after / n_raw * 100, 1))

removed_blood <- raw |>
  select(uniprot_id, gene, description) |>
  filter(gene %in% blood_genes, !gene %in% removed_genes)

filtered_proteins <- bind_rows(
  tibble(uniprot_id = raw$uniprot_id, gene = raw$gene,
         description = raw$description) |>
    filter(gene %in% removed_genes) |>
    mutate(removal_step = "HPA tissue filter"),
  removed_blood |> mutate(removal_step = "Blood contaminant"),
  as_tibble(annot_df) |>
    filter(!uniprot_id %in% rownames(dal$data)) |>
    select(uniprot_id, gene, description) |>
    mutate(removal_step = "Missingness")) |>
  distinct(uniprot_id, .keep_all = TRUE)

# 6. Outlier detection (4-method consensus, >=3/4)

pct_missing <- colMeans(is.na(dal$data)) * 100

miss_info <- dal$metadata |>
  select(Col_ID, Subject_ID, Group, Timepoint) |>
  mutate(pct_missing = pct_missing[Col_ID],
         prefix = str_remove(Col_ID, "_(Pre|Post)$")) |>
  mutate(delta_missing = if (n() == 2) abs(diff(pct_missing)) else NA_real_,
         .by = prefix)

miss_thresh  <- quantile(pct_missing, 0.75) + 1.5 * IQR(pct_missing)
delta_thresh <- quantile(miss_info$delta_missing, 0.75, na.rm = TRUE) +
  1.5 * IQR(miss_info$delta_missing, na.rm = TRUE)
miss_info$miss_flag <- miss_info$pct_missing > miss_thresh |
  coalesce(miss_info$delta_missing > delta_thresh, FALSE)

complete_mat <- dal$data[rowSums(is.na(dal$data)) == 0, ]
pca_pre <- run_pca(complete_mat, dal$metadata, log_transform = TRUE)
pc3     <- pca_pre$pca$x[, 1:3]
mahal   <- mahalanobis(pc3, colMeans(pc3), cov(pc3))
pca_flags <- tibble(Col_ID = colnames(dal$data), mahal_dist = mahal,
                    pca_flag = mahal > qchisq(1 - MAHAL_P, df = 3))

samp_med   <- apply(log2(dal$data), 2, median, na.rm = TRUE)
global_med <- median(samp_med)
mad_val    <- mad(samp_med)
mad_flags  <- tibble(Col_ID = names(samp_med), sample_median = samp_med,
                     mad_flag = abs(samp_med - global_med) > MAD_K * mad_val)

cor_mat <- cor(log2(dal$data), use = "pairwise.complete.obs")
med_cor <- apply(cor_mat, 2, \(x) median(x[x < 1], na.rm = TRUE))
cor_flags <- tibble(Col_ID = names(med_cor), median_cor = med_cor,
                    cor_flag = med_cor < median(med_cor) - MAD_K * mad(med_cor))

outlier_diag <- miss_info |>
  left_join(pca_flags, by = join_by(Col_ID)) |>
  left_join(mad_flags, by = join_by(Col_ID)) |>
  left_join(cor_flags, by = join_by(Col_ID)) |>
  mutate(n_flags = miss_flag + pca_flag + mad_flag + cor_flag,
         consensus_outlier = n_flags >= OUTLIER_K)

outlier_ids <- outlier_diag |> filter(consensus_outlier) |> pull(Col_ID)
n_outliers <- length(outlier_ids)
message(sprintf("Outliers: %d flagged (>=%d/4)", n_outliers, OUTLIER_K))

data_pre_outlier <- dal$data
meta_pre_outlier <- dal$metadata

if (n_outliers > 0) {
  dal <- filter_samples(dal, !(Col_ID %in% outlier_ids))
  message(sprintf("Removed: %s (%d remain)",
                  paste(outlier_ids, collapse = ", "), ncol(dal$data)))
}

# 7. Normalize (cycloess via proteoDA)

write_norm_report(dal, grouping_column = "Group_Time",
                  output_dir = RPT, filename = "01_norm_comparison.pdf",
                  overwrite = TRUE)
write_qc_report(dal, color_column = "Group_Time",
                output_dir = RPT, filename = "02_qc_pre.pdf", overwrite = TRUE)

dal <- normalize_data(dal, norm_method = "cycloess")

write_qc_report(dal, color_column = "Group_Time",
                output_dir = RPT, filename = "03_qc_post.pdf", overwrite = TRUE)

message(sprintf("Normalized: %d proteins x %d samples",
                nrow(dal$data), ncol(dal$data)))

# 8. Build xlsx

norm_df <- bind_cols(
  as_tibble(dal$annotation) |> select(uniprot_id, protein, gene, description),
  as_tibble(dal$data))

meta_out <- as.data.frame(metadata) |>
  mutate(QC_Status = if_else(Col_ID %in% outlier_ids, "Excluded", "Retained"))

pheno <- as.data.frame(read_excel("00_input/YvO_pheno_calc.xlsx"))
pheno <- pheno[rowSums(!is.na(pheno) & pheno != "") > 0, , drop = FALSE]

sample_cols <- dal$metadata$Col_ID
prot_miss <- tibble(
  gene    = norm_df$gene,
  n_miss  = rowSums(is.na(norm_df[, sample_cols])),
  pct_miss = round(100 * n_miss / length(sample_cols), 2))

samp_miss <- tibble(
  Col_ID   = sample_cols,
  n_miss   = colSums(is.na(dal$data)),
  pct_miss = round(100 * n_miss / nrow(dal$data), 2)) |>
  left_join(dal$metadata |> select(Col_ID, Group, Timepoint),
            by = join_by(Col_ID))

wb <- createWorkbook()
write_sheet(wb, "sample_metadata",     meta_out)
write_sheet(wb, "phenotype",           pheno)
write_sheet(wb, "normalized_matrix",   norm_df)
write_sheet(wb, "filter_log",          filter_log)
write_sheet(wb, "outlier_diagnostics", outlier_diag)
write_sheet(wb, "filtered_proteins",   filtered_proteins)
write_sheet(wb, "protein_missingness", prot_miss)
write_sheet(wb, "sample_missingness",  samp_miss)
saveWorkbook(wb, file.path(DAT, "01_normalization.xlsx"), overwrite = TRUE)

# CSV for downstream stages (benchmark, figures)
readr::write_csv(norm_df, file.path(DAT, "02_normalized.csv"))

# 9. Save R objects

saveRDS(dal, file.path(DAT, "03_DAList_normalized.rds"))

# Compute report/F00 diagnostic data
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
  mutate(iqr = apply(dal$data[, Col_ID], 2, IQR, na.rm = TRUE),
         Subject_ID = str_remove(Col_ID, "_(Pre|Post)$")) |>
  select(Col_ID, Subject_ID, Group, Timepoint, Group_Time, iqr)

grp_vec <- dal$metadata$Group_Time[match(colnames(dal$data), dal$metadata$Col_ID)]
eta2_vals <- apply(dal$data, 1, function(x) {
  ok <- !is.na(x)
  if (sum(ok) < 4) return(NA_real_)
  xk <- x[ok]; gk <- grp_vec[ok]
  ss_b <- sum(tapply(xk, gk, length) * (tapply(xk, gk, mean) - mean(xk))^2)
  ss_t <- sum((xk - mean(xk))^2)
  if (ss_t > 0) ss_b / ss_t else NA_real_
})

pca_post <- run_pca(dal$data, dal$metadata, log_transform = FALSE)

saveRDS(list(
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
  dal_ncol         = ncol(dal$data),
  mahal_p          = MAHAL_P,
  mad_k            = MAD_K
), file.path(DAT, "00_report_intermediates.rds"))

message(sprintf("Done: %d proteins x %d samples -> %s/",
                nrow(dal$data), ncol(dal$data), DAT))
