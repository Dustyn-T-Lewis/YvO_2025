#!/usr/bin/env Rscript
###############################################################################
#   YvO Normalization — DIA-MS proteomics (Young vs Old, skeletal muscle)
#   HPA tissue filter, 66% group missingness, consensus outlier detection,
#   cycloess normalization via proteoDA
#
#   Reports:  01_norm_comparison, 02_qc_pre, 03_qc_post, 04_diagnostics
#
#   Key references:
#     Thurman et al. 2023 JOSS 8:5184 (proteoDA)
#     Bolstad et al. 2003 Bioinformatics 19:185 (cyclic loess)
#     Valikangas et al. 2018 Brief Bioinform 19:1 (norm evaluation)
#     Arend et al. 2025 Brief Bioinform 26:bbaf201 (PRONE)
#     Brenes 2024 J Proteome Res 23:5274 (CV on linear scale)
#     Huang et al. 2024 Brief Bioinform 25:bbae129 (SEAOP outlier)
###############################################################################

if (!requireNamespace("pacman", quietly = TRUE)) install.packages("pacman")
pacman::p_load(
  proteoDA,
  readxl, readr, dplyr, tidyr, stringr,
  ggplot2, ggrepel, patchwork, openxlsx
)

set.seed(42)

setwd(rprojroot::find_rstudio_root_file())
input_dir  <- "00_input"
report_dir <- "01_normalization/b_reports"
data_dir   <- "01_normalization/c_data"

dir.create(report_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(data_dir,   recursive = TRUE, showWarnings = FALSE)

# Colors (self-contained; no dependency on figure palettes)
pal_group_time <- c(
  Young_Pre  = scales::alpha("#4393C3", 0.5), Young_Post = "#4393C3",
  Old_Pre    = scales::alpha("#D6604D", 0.5), Old_Post   = "#D6604D"
)
shape_tp <- c(Pre = 16, Post = 17)

# --- Helpers -----------------------------------------------------------------

impute_median <- function(mat) {
  for (j in seq_len(ncol(mat)))
    mat[is.na(mat[, j]), j] <- median(mat[, j], na.rm = TRUE)
  mat
}

run_pca <- function(mat, metadata, log_transform = TRUE) {
  if (log_transform) mat <- log2(impute_median(mat) + 1)
  else mat <- impute_median(mat)
  pca <- prcomp(t(mat), center = TRUE, scale. = TRUE)
  var_exp <- round(summary(pca)$importance[2, 1:3] * 100, 1)
  pc_df <- as.data.frame(pca$x[, 1:3]) %>%
    mutate(Col_ID = rownames(.)) %>%
    left_join(metadata, by = "Col_ID")
  list(pca = pca, pc_df = pc_df, var_exp = var_exp)
}

# =============================================================================
# 1. LOAD DATA
# =============================================================================

raw <- read_excel(file.path(input_dir, "YvO_raw.xlsx"))
annot_cols <- c("uniprot_id", "protein", "gene", "description", "n_seq")
annotation <- raw[, annot_cols]
intensity  <- raw[, setdiff(names(raw), annot_cols)]
cat(sprintf("Raw: %d proteins x %d samples\n", nrow(raw), ncol(intensity)))

metadata <- as.data.frame(read_excel(file.path(input_dir, "YvO_meta.xlsx")))
rownames(metadata) <- metadata$Col_ID
stopifnot("Sample mismatch" = setequal(colnames(intensity), metadata$Col_ID))
intensity <- intensity[, metadata$Col_ID]

n_raw <- nrow(annotation)
filter_log <- tibble(step = "Raw input", n_before = NA_integer_,
                     n_after = n_raw, n_removed = NA_integer_)

# =============================================================================
# 2. HPA TISSUE FILTER
# =============================================================================

hpa <- read_tsv(file.path(input_dir, "HPA_skeletal_muscle_annotations.tsv"),
                show_col_types = FALSE) %>%
  dplyr::select(Gene, Ensembl, Evidence,
                Protein_class    = `Protein class`,
                Subcellular_main = `Subcellular main location`,
                Interactions) %>%
  distinct(Gene, .keep_all = TRUE)

n_before <- nrow(annotation)
keep_mask <- annotation$gene %in% hpa$Gene
intensity  <- intensity[keep_mask, ]
annotation <- annotation[keep_mask, ] %>%
  left_join(hpa, by = c("gene" = "Gene"))
n_after <- nrow(annotation)

filter_log <- bind_rows(filter_log, tibble(
  step = "HPA non-muscle removal", n_before = n_before,
  n_after = n_after, n_removed = n_before - n_after))
cat(sprintf("HPA filter: %d -> %d (%d removed)\n", n_before, n_after, n_before - n_after))
removed_genes <- setdiff(raw$gene, annotation$gene)

# =============================================================================
# 3. ASSEMBLE DAList
# =============================================================================

if (any(duplicated(annotation$uniprot_id))) {
  n_before_dup <- nrow(annotation)
  annotation$row_mean <- rowMeans(data.matrix(intensity), na.rm = TRUE)
  keep_idx <- annotation %>%
    mutate(row_idx = row_number()) %>%
    group_by(uniprot_id) %>%
    slice_max(row_mean, n = 1, with_ties = FALSE) %>%
    pull(row_idx)
  annotation <- annotation[keep_idx, ]; intensity <- intensity[keep_idx, ]
  annotation$row_mean <- NULL
  filter_log <- bind_rows(filter_log, tibble(
    step = "Deduplication", n_before = n_before_dup,
    n_after = nrow(annotation), n_removed = n_before_dup - nrow(annotation)))
  cat(sprintf("Deduplicated: %d unique proteins\n", nrow(annotation)))
}

intensity_mat <- as.data.frame(data.matrix(intensity))
rownames(intensity_mat) <- annotation$uniprot_id
annot_df <- as.data.frame(annotation); rownames(annot_df) <- annotation$uniprot_id
meta_df  <- as.data.frame(metadata);   rownames(meta_df)  <- metadata$Col_ID

dal <- DAList(data = intensity_mat, annotation = annot_df, metadata = meta_df)
cat(sprintf("DAList: %d proteins x %d samples\n", nrow(dal$data), ncol(dal$data)))

# =============================================================================
# 4. QUALITY FILTERING
# =============================================================================

dal <- zero_to_missing(dal)
n_before <- nrow(dal$data)

group_prop <- dal$metadata %>%
  split(.$Group_Time) %>%
  lapply(function(g) rowMeans(!is.na(dal$data[, g$Col_ID, drop = FALSE]))) %>%
  bind_cols()
keep <- apply(group_prop >= 0.66, 1, any)
dal <- filter_proteins_by_annotation(dal, keep)

n_after <- nrow(dal$data)
filter_log <- bind_rows(filter_log, tibble(
  step = "Missingness (>=66% in >=1 group)", n_before = n_before,
  n_after = n_after, n_removed = n_before - n_after))
cat(sprintf("Missingness filter: %d -> %d (%d removed)\n",
            n_before, n_after, n_before - n_after))

filter_log <- filter_log %>% mutate(pct_of_raw = round(n_after / n_raw * 100, 1))

filtered_proteins <- bind_rows(
  tibble(uniprot_id = raw$uniprot_id, gene = raw$gene,
         description = raw$description) %>%
    filter(gene %in% removed_genes) %>%
    mutate(removal_step = "HPA tissue filter"),
  annot_df %>%
    filter(!uniprot_id %in% rownames(dal$data)) %>%
    dplyr::select(uniprot_id, gene, description) %>%
    mutate(removal_step = "Missingness (<66% in all groups)")
) %>% distinct(uniprot_id, .keep_all = TRUE)
print(filter_log)

# Filtering stacked bar
filter_plot_data <- filter_log %>%
  filter(!is.na(n_removed)) %>%
  mutate(step = factor(step, levels = step)) %>%
  pivot_longer(c(n_after, n_removed), names_to = "status", values_to = "n") %>%
  mutate(status = recode(status, n_after = "Retained", n_removed = "Removed"))

p_filter_bar <- ggplot(filter_plot_data, aes(x = step, y = n, fill = status)) +
  geom_col(width = 0.7) +
  geom_text(aes(label = n), position = position_stack(vjust = 0.5), size = 4.5) +
  scale_fill_manual(values = c(Retained = "#2166AC", Removed = "#B2182B")) +
  labs(x = NULL, y = "Proteins", fill = NULL,
       title = "Protein retention") +
  theme_minimal(base_size = 14) +
  theme(axis.text.x = element_text(angle = 25, hjust = 1, size = 11),
        legend.position = "top",
        plot.title = element_text(size = 14, face = "bold"),
        plot.margin = margin(2, 2, 2, 2))

# Per-sample missingness bars
miss_profile <- dal$metadata %>%
  dplyr::select(Col_ID, Group_Time) %>%
  mutate(n_detected = colSums(!is.na(dal$data[, Col_ID])),
         n_missing  = nrow(dal$data) - n_detected) %>%
  pivot_longer(c(n_detected, n_missing), names_to = "status", values_to = "n") %>%
  mutate(status = recode(status, n_detected = "Detected", n_missing = "Missing"))

p_miss_bars <- ggplot(miss_profile, aes(x = reorder(Col_ID, -n * (status == "Detected")),
                                         y = n, fill = status)) +
  geom_col(width = 0.8) +
  scale_fill_manual(values = c(Detected = "#2166AC", Missing = "#D6604D")) +
  facet_grid(~ Group_Time, scales = "free_x", space = "free_x") +
  labs(x = NULL, y = "Proteins", fill = NULL,
       title = "Per-sample detection vs missingness") +
  theme_minimal(base_size = 13) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 6),
        legend.position = "top",
        strip.text = element_text(face = "bold", size = 12),
        plot.title = element_text(size = 14, face = "bold"),
        plot.margin = margin(2, 2, 2, 2))

# =============================================================================
# 5. CV VARIABILITY
# =============================================================================

# CV on linear-scale intensities (Brenes 2024)
lin_dat <- dal$data

cv_by_group <- dal$metadata %>%
  split(.$Group_Time) %>%
  lapply(function(g) {
    sub <- lin_dat[, g$Col_ID, drop = FALSE]
    tibble(uniprot_id = rownames(sub),
           cv = apply(sub, 1, function(x) {
             x <- x[!is.na(x)]
             if (length(x) < 2) return(NA_real_)
             sd(x) / mean(x)
           }),
           Group_Time = g$Group_Time[1])
  }) %>% bind_rows()

p_cv_violin <- ggplot(cv_by_group, aes(x = Group_Time, y = cv, fill = Group_Time)) +
  geom_violin(alpha = 0.7, draw_quantiles = c(0.25, 0.5, 0.75)) +
  scale_fill_manual(values = pal_group_time, guide = "none") +
  coord_cartesian(ylim = c(0, quantile(cv_by_group$cv, 0.99, na.rm = TRUE))) +
  labs(x = NULL, y = "CV (linear scale)",
       title = "Inter-individual protein variability by group") +
  theme_minimal(base_size = 12)

cv_wide <- cv_by_group %>% pivot_wider(names_from = Group_Time, values_from = cv)

plot_cv_scatter <- function(df, x_col, y_col) {
  ggplot(df, aes(x = .data[[x_col]], y = .data[[y_col]])) +
    geom_point(alpha = 0.3, size = 1, color = "gray30") +
    geom_abline(slope = 1, linetype = "dashed", color = "red", alpha = 0.5) +
    geom_smooth(method = "lm", se = FALSE, linewidth = 0.6, color = "#2166AC") +
    coord_cartesian(xlim = c(0, quantile(df[[x_col]], 0.99, na.rm = TRUE)),
                    ylim = c(0, quantile(df[[y_col]], 0.99, na.rm = TRUE))) +
    labs(x = paste("CV:", x_col), y = paste("CV:", y_col)) +
    theme_minimal(base_size = 11)
}

p_cv1 <- plot_cv_scatter(cv_wide, "Young_Pre", "Young_Post")
p_cv2 <- plot_cv_scatter(cv_wide, "Old_Pre",   "Old_Post")
p_cv3 <- plot_cv_scatter(cv_wide, "Young_Pre",  "Old_Pre")

# =============================================================================
# 6. OUTLIER DETECTION & REMOVAL
# =============================================================================

pct_missing <- colMeans(is.na(dal$data)) * 100

delta_missing <- dal$metadata %>%
  dplyr::select(Col_ID, Subject_ID, Group, Timepoint) %>%
  mutate(pct_missing = pct_missing[Col_ID],
         prefix = str_remove(Col_ID, "_(Pre|Post)$")) %>%
  arrange(prefix, match(Timepoint, c("Pre", "Post"))) %>%
  group_by(prefix) %>%
  mutate(delta_missing = abs(pct_missing - dplyr::first(pct_missing))) %>%
  ungroup()

miss_threshold  <- quantile(pct_missing, 0.75) + 1.5 * IQR(pct_missing)
delta_vals      <- delta_missing$delta_missing[delta_missing$Timepoint != "Pre"]
delta_threshold <- quantile(delta_vals, 0.75) + 1.5 * IQR(delta_vals)

delta_missing <- delta_missing %>%
  mutate(miss_flag = pct_missing > miss_threshold |
           (Timepoint != "Pre" & delta_missing > delta_threshold))

zero_miss_mat <- dal$data[rowSums(is.na(dal$data)) == 0, ]
cat(sprintf("Pre-norm PCA: %d complete proteins (of %d)\n",
            nrow(zero_miss_mat), nrow(dal$data)))
pca_out   <- run_pca(zero_miss_mat, dal$metadata, log_transform = TRUE)
pc_scores <- pca_out$pca$x[, 1:3]
mahal_dist <- mahalanobis(pc_scores, colMeans(pc_scores), cov(pc_scores))
pca_flags <- tibble(Col_ID = colnames(dal$data), mahal_dist = mahal_dist,
                    pca_flag = mahal_dist > qchisq(0.99, df = 3))

sample_medians <- apply(log2(dal$data + 1), 2, median, na.rm = TRUE)
global_median  <- median(sample_medians)
mad_val        <- mad(sample_medians)
mad_flags <- tibble(Col_ID = names(sample_medians), sample_median = sample_medians,
                    mad_flag = abs(sample_medians - global_median) > 3 * mad_val)

outlier_diag <- delta_missing %>%
  left_join(pca_flags, by = "Col_ID") %>%
  left_join(mad_flags, by = "Col_ID") %>%
  mutate(n_flags = miss_flag + pca_flag + mad_flag,
         consensus_outlier = n_flags >= 3)

n_outliers <- sum(outlier_diag$consensus_outlier)
cat(sprintf("Outlier consensus: %d sample(s) flagged (3/3 methods)\n", n_outliers))

p_out_a <- ggplot(outlier_diag, aes(pct_missing, delta_missing)) +
  geom_point(aes(color = consensus_outlier, shape = Timepoint), size = 3) +
  geom_text_repel(data = . %>% filter(consensus_outlier),
                  aes(label = prefix), size = 2.5, color = "red", max.overlaps = 20) +
  geom_hline(yintercept = delta_threshold, linetype = "dashed", color = "red", alpha = 0.5) +
  geom_vline(xintercept = miss_threshold, linetype = "dashed", color = "red", alpha = 0.5) +
  scale_color_manual(values = c("FALSE" = "gray40", "TRUE" = "red")) +
  scale_shape_manual(values = shape_tp) +
  labs(x = "% Missing", y = "|Delta Missingness|", title = "A: Paired Missingness") +
  theme_minimal(base_size = 11) + theme(legend.position = "bottom")

pca_outlier_df <- pca_out$pc_df %>%
  left_join(outlier_diag %>% dplyr::select(Col_ID, consensus_outlier, prefix),
            by = "Col_ID")
p_out_b <- ggplot(pca_outlier_df, aes(PC1, PC2, shape = Timepoint)) +
  geom_point(aes(color = ifelse(consensus_outlier, "Outlier", as.character(Group_Time))),
             size = 3.5, alpha = 0.85) +
  geom_text_repel(data = . %>% filter(consensus_outlier),
                  aes(label = prefix), size = 2.5, color = "red", max.overlaps = 20) +
  scale_color_manual(values = c(pal_group_time, Outlier = "red"), name = "Group") +
  scale_shape_manual(values = shape_tp) +
  labs(x = sprintf("PC1 (%.1f%%)", pca_out$var_exp[1]),
       y = sprintf("PC2 (%.1f%%)", pca_out$var_exp[2]),
       title = "B: PCA (pre-normalization)") +
  theme_minimal(base_size = 11) + theme(legend.position = "bottom")

p_out_c <- ggplot(outlier_diag, aes(reorder(prefix, sample_median), sample_median)) +
  geom_point(aes(color = consensus_outlier), size = 2.5) +
  geom_text_repel(data = . %>% filter(consensus_outlier),
                  aes(label = prefix), size = 2.5, color = "red") +
  geom_hline(yintercept = global_median) +
  geom_hline(yintercept = global_median + c(-3, 3) * mad_val,
             linetype = "dashed", color = "red", alpha = 0.5) +
  scale_color_manual(values = c("FALSE" = "gray40", "TRUE" = "red")) +
  labs(x = "Sample", y = "Median log2 intensity", title = "C: MAD Outliers") +
  theme_minimal(base_size = 11) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 4),
        legend.position = "none")

# Remove outliers (paired)
if (n_outliers > 0) {
  flagged_ids <- outlier_diag %>% filter(consensus_outlier) %>% pull(Col_ID)
  flagged_prefixes <- unique(str_remove(flagged_ids, "_(Pre|Post)$"))
  remove_ids <- dal$metadata$Col_ID[
    str_remove(dal$metadata$Col_ID, "_(Pre|Post)$") %in% flagged_prefixes]
  cat(sprintf("Removing (with pairs): %s\n", paste(flagged_prefixes, collapse = ", ")))
  dal <- filter_samples(dal, !(Col_ID %in% remove_ids))
  cat(sprintf("%d samples remain\n", ncol(dal$data)))
}

# =============================================================================
# 7. NORMALIZATION & proteoDA REPORTS
# =============================================================================

write_norm_report(dal, grouping_column = "Group_Time", output_dir = report_dir,
                  filename = "01_norm_comparison.pdf", overwrite = TRUE)
dal_pre <- dal
saveRDS(dal_pre, file.path(data_dir, "01_DAList_prenorm.rds"))

write_qc_report(dal, color_column = "Group_Time", label_column = "Col_ID",
                output_dir = report_dir, filename = "02_qc_pre.pdf",
                overwrite = TRUE)

dal <- normalize_data(dal, norm_method = "cycloess")
cat(sprintf("Cycloess: %d proteins x %d samples\n", nrow(dal$data), ncol(dal$data)))

write_qc_report(dal, color_column = "Group_Time", label_column = "Col_ID",
                output_dir = report_dir, filename = "03_qc_post.pdf",
                overwrite = TRUE)

# =============================================================================
# 8. REPORT 4: DIAGNOSTICS (3 pages)
# =============================================================================

# Variance decomposition: eta-squared per protein (Group_Time ANOVA)
log_dat <- log2(dal$data + 1)
grp_vec <- dal$metadata$Group_Time[match(colnames(log_dat), dal$metadata$Col_ID)]

var_decomp <- apply(log_dat, 1, function(x) {
  keep <- !is.na(x)
  if (sum(keep) < 4) return(c(within = NA_real_, between = NA_real_, eta2 = NA_real_))
  xk <- x[keep]; gk <- grp_vec[keep]
  grand_mean <- mean(xk)
  grp_means  <- tapply(xk, gk, mean)
  grp_n      <- tapply(xk, gk, length)
  ss_between <- sum(grp_n * (grp_means - grand_mean)^2)
  ss_total   <- sum((xk - grand_mean)^2)
  ss_within  <- ss_total - ss_between
  df_b <- length(grp_means) - 1
  df_w <- sum(keep) - length(grp_means)
  ms_w <- if (df_w > 0) ss_within / df_w else NA_real_
  ms_b <- if (df_b > 0) ss_between / df_b else NA_real_
  eta2 <- if (ss_total > 0) ss_between / ss_total else NA_real_
  c(within = ms_w, between = ms_b, eta2 = eta2)
})
var_df <- as.data.frame(t(var_decomp)) %>%
  mutate(uniprot_id = rownames(log_dat)) %>%
  filter(!is.na(eta2))

p_eta2 <- ggplot(var_df, aes(x = eta2)) +
  geom_histogram(bins = 50, fill = "#2166AC", color = "white", alpha = 0.8) +
  geom_vline(xintercept = median(var_df$eta2), linetype = "dashed", color = "red") +
  annotate("text", x = median(var_df$eta2) + 0.02, y = Inf, vjust = 2,
           label = sprintf("median = %.2f", median(var_df$eta2)),
           size = 3.5, color = "red") +
  labs(x = expression(eta^2 ~ "(between-group / total variance)"),
       y = "Proteins",
       title = "Variance partition by experimental group") +
  theme_minimal(base_size = 12)

p_var_scatter <- ggplot(var_df, aes(x = within, y = between)) +
  geom_point(alpha = 0.3, size = 1, color = "gray30") +
  geom_abline(slope = 1, linetype = "dashed", color = "red", alpha = 0.5) +
  geom_smooth(method = "lm", se = FALSE, linewidth = 0.6, color = "#2166AC") +
  scale_x_log10() + scale_y_log10() +
  labs(x = "Within-group MS (log10)", y = "Between-group MS (log10)",
       title = "Within vs between-group variance") +
  theme_minimal(base_size = 11)

# --- Post-norm PCA -----------------------------------------------------------
pca_post <- run_pca(dal$data, dal$metadata, log_transform = FALSE)

p_pca_gt <- ggplot(pca_post$pc_df, aes(PC1, PC2, color = Group_Time, shape = Timepoint)) +
  geom_point(size = 3.5, alpha = 0.85) +
  stat_ellipse(aes(group = Group_Time), type = "norm", level = 0.68, linewidth = 0.7) +
  scale_color_manual(values = pal_group_time) + scale_shape_manual(values = shape_tp) +
  labs(x = sprintf("PC1 (%.1f%%)", pca_post$var_exp[1]),
       y = sprintf("PC2 (%.1f%%)", pca_post$var_exp[2]),
       title = "Post-normalization PCA") +
  theme_minimal(base_size = 12) + theme(legend.position = "bottom")

# --- Build PDF ---------------------------------------------------------------
pdf(file.path(report_dir, "04_diagnostics.pdf"), width = 20, height = 13)

# Page 1: Filtering
print(
  p_filter_bar / p_miss_bars +
    plot_layout(heights = c(1, 1.2)) +
    plot_annotation(title = "Protein Filtering",
                    subtitle = sprintf("%d raw -> %d retained | %d samples",
                                       n_raw, nrow(dal$data), ncol(dal$data)),
                    theme = theme(plot.title = element_text(size = 18, face = "bold"),
                                  plot.subtitle = element_text(size = 14)))
)

# Page 2: Variance
print(
  (p_cv_violin | p_eta2 | p_var_scatter) /
    (p_cv1 | p_cv2 | p_cv3 | p_pca_gt) +
    plot_layout(heights = c(1, 1)) +
    plot_annotation(title = "Variance")
)

# Page 3: Outlier Diagnostics
print(
  (p_out_a | p_out_b | p_out_c) +
    plot_annotation(title = "Outlier Diagnostics (3-method consensus)")
)

dev.off()

# =============================================================================
# 9. NORMALIZATION METHOD COMPARISON (PRONE-style ranking)
# =============================================================================

norm_metric <- function(mat, groups, metric = "cv") {
  grp_list <- split(seq_len(ncol(mat)), groups)
  if (metric == "cor") {
    vals <- unlist(lapply(grp_list, function(idx) {
      sub <- mat[, idx, drop = FALSE]
      if (ncol(sub) < 2) return(numeric(0))
      cm <- cor(sub, use = "pairwise.complete.obs")
      cm[lower.tri(cm)]
    }))
    return(mean(vals, na.rm = TRUE))
  }
  vals <- unlist(lapply(grp_list, function(idx) {
    sub <- mat[, idx, drop = FALSE]
    apply(sub, 1, function(x) {
      x <- x[!is.na(x)]
      if (length(x) < 2) return(NA_real_)
      if (metric == "cv") sd(x) / abs(mean(x)) else mad(x, constant = 1)
    })
  }))
  median(vals, na.rm = TRUE)
}

dal_pre$metadata$group <- factor(dal_pre$metadata$Group_Time,
                                  levels = c("Young_Pre","Young_Post","Old_Pre","Old_Post"))
norm_methods <- c("log2", "median", "mean", "vsn", "quantile", "cycloess", "rlr", "gi")
norm_scores <- lapply(norm_methods, function(nm) {
  dal_n <- tryCatch(normalize_data(dal_pre, norm_method = nm), error = function(e) NULL)
  if (is.null(dal_n)) return(NULL)
  mat_n <- as.matrix(dal_n$data); grps <- dal_n$metadata$group
  tibble(norm = nm,
         PCV  = round(norm_metric(mat_n, grps, "cv"),  4),
         PMAD = round(norm_metric(mat_n, grps, "mad"), 4),
         COR  = round(norm_metric(mat_n, grps, "cor"), 4))
}) |> bind_rows()
norm_scores <- norm_scores %>%
  mutate(PCV_rank  = rank(PCV), PMAD_rank = rank(PMAD), COR_rank = rank(-COR),
         composite = (PCV_rank + PMAD_rank + COR_rank) / 3) %>%
  arrange(composite)
write_csv(norm_scores, file.path(data_dir, "04_norm_quality_scores.csv"))
print(norm_scores %>% dplyr::select(norm, PCV, PMAD, COR, composite))

# =============================================================================
# 10. EXPORT
# =============================================================================

export_df <- bind_cols(
  as_tibble(dal$annotation) %>% dplyr::select(uniprot_id, protein, gene, description),
  as_tibble(dal$data))

write_csv(export_df, file.path(data_dir, "02_normalized.csv"))
saveRDS(dal, file.path(data_dir, "03_DAList_normalized.rds"))

# --- Supplementary workbook ---------------------------------------------------
add_sheet <- function(wb, sheet_name, title_text, notes_text, data_df) {
  addWorksheet(wb, sheet_name)
  title_style <- createStyle(fontSize = 13, textDecoration = "bold")
  notes_style <- createStyle(fontSize = 10, fontColour = "#555555", wrapText = TRUE)
  header_style <- createStyle(textDecoration = "bold", border = "Bottom",
                              fgFill = "#DCE6F1")
  writeData(wb, sheet_name, x = title_text, startRow = 1, startCol = 1)
  mergeCells(wb, sheet_name, cols = 1:ncol(data_df), rows = 1)
  addStyle(wb, sheet_name, title_style, rows = 1, cols = 1)
  writeData(wb, sheet_name, x = notes_text, startRow = 2, startCol = 1)
  mergeCells(wb, sheet_name, cols = 1:ncol(data_df), rows = 2)
  addStyle(wb, sheet_name, notes_style, rows = 2, cols = 1)
  writeData(wb, sheet_name, x = data_df, startRow = 4, headerStyle = header_style)
  freezePane(wb, sheet_name, firstActiveRow = 5, firstActiveCol = 1)
  setColWidths(wb, sheet_name, cols = 1:ncol(data_df), widths = "auto")
}

wb <- createWorkbook()

add_sheet(wb, "Pipeline_Summary",
  "Protein Filtering Pipeline Summary",
  "step: filtering stage | n_before/n_after: protein counts | n_removed: proteins lost | pct_of_raw: cumulative retention",
  filter_log)

add_sheet(wb, "Normalization_Ranking",
  "Normalization Method Comparison (PRONE framework)",
  "PCV: pooled CV (lower=better) | PMAD: pooled MAD (lower=better) | COR: intraclass correlation (higher=better) | composite: mean rank (Arend et al. 2025, Brief Bioinform 26:bbaf201)",
  norm_scores)

add_sheet(wb, "Outlier_Diagnostics",
  "Per-Sample Outlier Diagnostics (3-method consensus)",
  "miss_flag: paired missingness outlier | pca_flag: Mahalanobis PC1-3 (p<0.01) | mad_flag: median intensity >3 MAD | consensus_outlier: TRUE if all 3 agree",
  outlier_diag)

add_sheet(wb, "Filtered_Proteins",
  "All Proteins Removed by Filtering",
  "removal_step: HPA tissue filter or Missingness (<66% in all groups) | uniprot_id, gene, description: protein identifiers",
  filtered_proteins)

saveWorkbook(wb, file.path(data_dir, "05_normalization_supp.xlsx"), overwrite = TRUE)

cat(sprintf("Exported: %d proteins x %d samples\n", nrow(dal$data), ncol(dal$data)))
