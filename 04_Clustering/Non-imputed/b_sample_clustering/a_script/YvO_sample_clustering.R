################################################################################
#
#   YvO Sample Clustering — Non-imputed Data
#
#   Workflow:
#     1.  Setup & data loading + metadata
#     2.  Handle NAs (column median imputation for PCA/t-SNE/UMAP)
#     3.  PCA (scree, PC1xPC2, PC2xPC3 biplots)
#     4.  Hierarchical clustering (Euclidean, dendrogram colored by group)
#     5.  t-SNE (perplexity: 5, 10, 15, 30)
#     6.  UMAP (default + tuned)
#     7.  Outlier detection (Mahalanobis in PC space)
#     8.  Sample correlation heatmap
#     9.  Export
#
#   Input: 01_normalization/c_data/01_normalized.csv (2124 proteins, NAs present)
#   Note: NAs imputed with column medians for dimensionality reduction only.
#
################################################################################

# ==============================================================================
# 1.  SETUP & DATA LOADING
# ==============================================================================

cat("=== YvO Sample Clustering (Non-imputed) ===\n\n")
cat(">> 1 -- Setup & Data Loading\n")

suppressPackageStartupMessages({
  library(tidyverse)
  library(Rtsne)
  library(umap)
  library(factoextra)
  library(ggrepel)
  library(patchwork)
  library(pheatmap)
  library(dendextend)
})

base_dir   <- normalizePath(file.path(dirname(getwd()), "..", "..", ".."), mustWork = TRUE)
DATA_FILE  <- file.path(base_dir, "01_normalization", "c_data", "01_normalized.csv")
REPORT_DIR <- file.path(base_dir, "04_Clustering", "Non-imputed", "b_sample_clustering", "b_reports")
DATA_DIR   <- file.path(base_dir, "04_Clustering", "Non-imputed", "b_sample_clustering", "c_data")

dir.create(REPORT_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(DATA_DIR,   recursive = TRUE, showWarnings = FALSE)

df <- read_csv(DATA_FILE, show_col_types = FALSE)
ann_cols   <- c("uniprot_id", "protein", "gene", "description")
samp_names <- setdiff(names(df), ann_cols)
mat        <- as.matrix(df[, samp_names])
rownames(mat) <- df$uniprot_id

n_missing <- sum(is.na(mat))
pct_missing <- round(100 * n_missing / length(mat), 1)
cat(sprintf("   Loaded: %d proteins x %d samples (%d NAs, %.1f%%)\n",
            nrow(mat), ncol(mat), n_missing, pct_missing))

# Build metadata
meta <- tibble(sample_id = samp_names) %>%
  mutate(
    prefix   = str_extract(sample_id, "^[A-Z]+"),
    subj_num = str_extract(sample_id, "S\\d+"),
    time     = str_extract(sample_id, "(Pre|Post)$"),
    age      = if_else(str_detect(prefix, "^O"), "Old", "Young"),
    subject  = paste0(prefix, "_", subj_num),
    group    = paste(age, time, sep = "_")
  )
meta$group <- factor(meta$group, levels = c("Young_Pre", "Young_Post", "Old_Pre", "Old_Post"))

group_colors <- c("Young_Pre" = "#4A90D9", "Young_Post" = "#2E5FA1",
                   "Old_Pre" = "#D94A4A", "Old_Post" = "#A12E2E")

# ==============================================================================
# 2.  HANDLE NAs (Column Median Imputation for Dimensionality Reduction)
# ==============================================================================

cat("\n>> 2 -- Handle Missing Values\n")

# Transpose: samples as rows, proteins as columns
mat_t <- t(mat)

# Impute NAs with column medians (each column = one protein across samples)
cat("   Imputing NAs with column medians for PCA/t-SNE/UMAP...\n")
mat_t_imp <- mat_t
for (j in 1:ncol(mat_t_imp)) {
  na_idx <- is.na(mat_t_imp[, j])
  if (any(na_idx)) {
    col_median <- median(mat_t_imp[, j], na.rm = TRUE)
    mat_t_imp[na_idx, j] <- col_median
  }
}

# Remove any remaining all-NA columns
all_na_cols <- apply(mat_t_imp, 2, function(x) all(is.na(x)))
if (any(all_na_cols)) {
  cat(sprintf("   Removed %d all-NA proteins\n", sum(all_na_cols)))
  mat_t_imp <- mat_t_imp[, !all_na_cols]
}

# Remove zero-variance columns
col_vars <- apply(mat_t_imp, 2, var, na.rm = TRUE)
zero_var <- col_vars == 0 | is.na(col_vars)
if (any(zero_var)) {
  cat(sprintf("   Removed %d zero-variance proteins\n", sum(zero_var)))
  mat_t_imp <- mat_t_imp[, !zero_var]
}

cat(sprintf("   Matrix for DR: %d samples x %d proteins\n",
            nrow(mat_t_imp), ncol(mat_t_imp)))
cat("   NOTE: Column median imputation used ONLY for dimensionality reduction.\n")

# ==============================================================================
# 3.  PCA
# ==============================================================================

cat("\n>> 3 -- PCA\n")

pca_res <- prcomp(mat_t_imp, center = TRUE, scale. = TRUE)
pca_df <- as.data.frame(pca_res$x[, 1:10])
pca_df$sample_id <- rownames(mat_t)
pca_df <- left_join(pca_df, meta, by = "sample_id")

var_explained <- summary(pca_res)$importance[2, ] * 100

p_scree <- ggplot(tibble(PC = 1:20, Var = var_explained[1:20]),
                   aes(x = PC, y = Var)) +
  geom_col(fill = "#D94A4A") +
  geom_line(linewidth = 0.8) +
  geom_point(size = 2) +
  labs(x = "Principal Component", y = "Variance Explained (%)",
       title = "PCA Scree Plot") +
  theme_bw()

p_pc12 <- ggplot(pca_df, aes(x = PC1, y = PC2, color = group)) +
  geom_point(size = 3) +
  geom_text_repel(aes(label = sample_id), size = 2, max.overlaps = 10) +
  scale_color_manual(values = group_colors) +
  labs(x = sprintf("PC1 (%.1f%%)", var_explained[1]),
       y = sprintf("PC2 (%.1f%%)", var_explained[2]),
       title = "PCA: PC1 vs PC2") +
  theme_bw()

p_pc23 <- ggplot(pca_df, aes(x = PC2, y = PC3, color = group)) +
  geom_point(size = 3) +
  geom_text_repel(aes(label = sample_id), size = 2, max.overlaps = 10) +
  scale_color_manual(values = group_colors) +
  labs(x = sprintf("PC2 (%.1f%%)", var_explained[2]),
       y = sprintf("PC3 (%.1f%%)", var_explained[3]),
       title = "PCA: PC2 vs PC3") +
  theme_bw()

p_pca <- p_scree / (p_pc12 | p_pc23) +
  plot_annotation(title = "PCA — Sample Clustering (Non-imputed, median-imputed for DR)")
ggsave(file.path(REPORT_DIR, "01_pca.pdf"), p_pca, width = 14, height = 10)

cat("   Saved: 01_pca.pdf\n")

# ==============================================================================
# 4.  HIERARCHICAL CLUSTERING
# ==============================================================================

cat("\n>> 4 -- Hierarchical Clustering\n")

dist_samples <- dist(mat_t_imp, method = "euclidean")
hc_samples   <- hclust(dist_samples, method = "ward.D2")

dend <- as.dendrogram(hc_samples)
group_vec <- meta$group[match(labels(dend), meta$sample_id)]
dend <- color_branches(dend, clusters = as.numeric(factor(group_vec)),
                        col = group_colors[levels(meta$group)])

pdf(file.path(REPORT_DIR, "02_dendrogram.pdf"), width = 14, height = 6)
par(mar = c(8, 4, 3, 1))
plot(dend, main = "Sample Dendrogram (Euclidean, Ward.D2, Non-imputed)")
legend("topright", legend = levels(meta$group), fill = group_colors, cex = 0.8)
dev.off()

cat("   Saved: 02_dendrogram.pdf\n")

# ==============================================================================
# 5.  t-SNE
# ==============================================================================

cat("\n>> 5 -- t-SNE\n")

perplexities <- c(5, 10, 15, 30)
tsne_plots <- list()

for (perp in perplexities) {
  if (perp >= nrow(mat_t_imp)) next
  set.seed(42)
  tsne_res <- Rtsne(mat_t_imp, perplexity = perp, dims = 2, check_duplicates = FALSE)
  tsne_df <- tibble(
    tSNE1 = tsne_res$Y[, 1],
    tSNE2 = tsne_res$Y[, 2],
    sample_id = rownames(mat_t)
  ) %>% left_join(meta, by = "sample_id")

  tsne_plots[[as.character(perp)]] <- ggplot(tsne_df, aes(x = tSNE1, y = tSNE2, color = group)) +
    geom_point(size = 3) +
    scale_color_manual(values = group_colors) +
    labs(title = sprintf("perplexity = %d", perp)) +
    theme_bw(base_size = 9) +
    theme(legend.position = "none")
}

p_tsne <- wrap_plots(tsne_plots, ncol = 2) +
  plot_annotation(title = "t-SNE — Sample Clustering (Non-imputed)")
ggsave(file.path(REPORT_DIR, "03_tsne.pdf"), p_tsne, width = 10, height = 8)

cat("   Saved: 03_tsne.pdf\n")

# ==============================================================================
# 6.  UMAP
# ==============================================================================

cat("\n>> 6 -- UMAP\n")

set.seed(42)
umap_default <- umap(mat_t_imp)
umap_df <- tibble(
  UMAP1 = umap_default$layout[, 1],
  UMAP2 = umap_default$layout[, 2],
  sample_id = rownames(mat_t)
) %>% left_join(meta, by = "sample_id")

p_umap1 <- ggplot(umap_df, aes(x = UMAP1, y = UMAP2, color = group)) +
  geom_point(size = 3) +
  geom_text_repel(aes(label = sample_id), size = 2, max.overlaps = 10) +
  scale_color_manual(values = group_colors) +
  labs(title = "UMAP (default)") +
  theme_bw()

umap_config <- umap.defaults
umap_config$n_neighbors <- 10
umap_config$min_dist <- 0.01
set.seed(42)
umap_tuned <- umap(mat_t_imp, config = umap_config)
umap_df2 <- tibble(
  UMAP1 = umap_tuned$layout[, 1],
  UMAP2 = umap_tuned$layout[, 2],
  sample_id = rownames(mat_t)
) %>% left_join(meta, by = "sample_id")

p_umap2 <- ggplot(umap_df2, aes(x = UMAP1, y = UMAP2, color = group)) +
  geom_point(size = 3) +
  geom_text_repel(aes(label = sample_id), size = 2, max.overlaps = 10) +
  scale_color_manual(values = group_colors) +
  labs(title = "UMAP (n_neighbors=10, min_dist=0.01)") +
  theme_bw()

p_umap <- p_umap1 | p_umap2
p_umap <- p_umap + plot_annotation(title = "UMAP — Sample Clustering (Non-imputed)")
ggsave(file.path(REPORT_DIR, "04_umap.pdf"), p_umap, width = 14, height = 6)

cat("   Saved: 04_umap.pdf\n")

# ==============================================================================
# 7.  OUTLIER DETECTION
# ==============================================================================

cat("\n>> 7 -- Outlier Detection\n")

pc_scores <- pca_res$x[, 1:10]
center <- colMeans(pc_scores)
cov_mat <- cov(pc_scores)
mahal_dist <- mahalanobis(pc_scores, center, cov_mat)

outlier_df <- tibble(
  sample_id = rownames(mat_t),
  mahal_dist = mahal_dist,
  p_value = pchisq(mahal_dist, df = 10, lower.tail = FALSE),
  outlier = p_value < 0.01
) %>% left_join(meta, by = "sample_id")

n_outliers <- sum(outlier_df$outlier)
cat(sprintf("   Outliers (Mahalanobis, p < 0.01): %d\n", n_outliers))
if (n_outliers > 0) {
  cat("   Flagged samples:\n")
  print(outlier_df %>% filter(outlier) %>% dplyr::select(sample_id, group, mahal_dist, p_value))
}

# ==============================================================================
# 8.  SAMPLE CORRELATION HEATMAP
# ==============================================================================

cat("\n>> 8 -- Sample Correlation Heatmap\n")

# Use pairwise complete obs for correlation with NAs
cor_samples <- cor(mat, method = "pearson", use = "pairwise.complete.obs")

col_ann <- data.frame(
  Age   = meta$age,
  Time  = meta$time,
  Group = meta$group,
  row.names = meta$sample_id
)

pdf(file.path(REPORT_DIR, "05_sample_correlation.pdf"), width = 12, height = 10)
pheatmap(cor_samples,
         annotation_col = col_ann,
         annotation_row = col_ann,
         show_rownames = FALSE,
         show_colnames = FALSE,
         main = "Sample Pearson correlation (Non-imputed, pairwise complete)",
         fontsize = 8)
dev.off()

cat("   Saved: 05_sample_correlation.pdf\n")

# ==============================================================================
# 9.  EXPORT
# ==============================================================================

cat("\n>> 9 -- Export\n")

pca_export <- pca_df %>% dplyr::select(sample_id, PC1:PC10, age, time, group)
write_csv(pca_export, file.path(DATA_DIR, "sample_pca_scores.csv"))
write_csv(outlier_df, file.path(DATA_DIR, "sample_outlier_flags.csv"))

cat("   Saved: sample_pca_scores.csv, sample_outlier_flags.csv\n")
cat("   NOTE: Column median imputation used only for dimensionality reduction.\n")

cat("\n")
sessionInfo()

cat("\n=== YvO Sample Clustering (Non-imputed) Complete ===\n")
