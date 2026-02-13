################################################################################
#
#   YvO Sample Clustering — Imputed Data
#
#   Workflow:
#     1.  Setup & data loading + metadata
#     2.  PCA (scree, PC1xPC2, PC2xPC3 biplots)
#     3.  Hierarchical clustering (Euclidean, dendrogram colored by group)
#     4.  t-SNE (perplexity: 5, 10, 15, 30)
#     5.  UMAP (default + tuned)
#     6.  Outlier detection (Mahalanobis in PC space)
#     7.  Sample correlation heatmap
#     8.  Export
#
#   Input: 02_Imputation/c_data/01_imputed.csv (2124 proteins x 62 samples)
#
################################################################################

# ==============================================================================
# 1.  SETUP & DATA LOADING
# ==============================================================================

cat("=== YvO Sample Clustering (Imputed) ===\n\n")
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
DATA_FILE  <- file.path(base_dir, "02_Imputation", "c_data", "01_imputed.csv")
REPORT_DIR <- file.path(base_dir, "04_Clustering", "Imputed", "b_sample_clustering", "b_reports")
DATA_DIR   <- file.path(base_dir, "04_Clustering", "Imputed", "b_sample_clustering", "c_data")

dir.create(REPORT_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(DATA_DIR,   recursive = TRUE, showWarnings = FALSE)

df <- read_csv(DATA_FILE, show_col_types = FALSE)
ann_cols   <- c("uniprot_id", "protein", "gene", "description")
samp_names <- setdiff(names(df), ann_cols)
mat        <- as.matrix(df[, samp_names])
rownames(mat) <- df$uniprot_id

cat(sprintf("   Loaded: %d proteins x %d samples\n", nrow(mat), ncol(mat)))

# Transpose: samples as rows
mat_t <- t(mat)

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
# 2.  PCA
# ==============================================================================

cat("\n>> 2 -- PCA\n")

pca_res <- prcomp(mat_t, center = TRUE, scale. = TRUE)
pca_df <- as.data.frame(pca_res$x[, 1:10])
pca_df$sample_id <- rownames(mat_t)
pca_df <- left_join(pca_df, meta, by = "sample_id")

var_explained <- summary(pca_res)$importance[2, ] * 100

# Scree plot
p_scree <- ggplot(tibble(PC = 1:20, Var = var_explained[1:20]),
                   aes(x = PC, y = Var)) +
  geom_col(fill = "#4A90D9") +
  geom_line(linewidth = 0.8) +
  geom_point(size = 2) +
  labs(x = "Principal Component", y = "Variance Explained (%)",
       title = "PCA Scree Plot") +
  theme_bw()

# PC1 x PC2
p_pc12 <- ggplot(pca_df, aes(x = PC1, y = PC2, color = group)) +
  geom_point(size = 3) +
  geom_text_repel(aes(label = sample_id), size = 2, max.overlaps = 10) +
  scale_color_manual(values = group_colors) +
  labs(x = sprintf("PC1 (%.1f%%)", var_explained[1]),
       y = sprintf("PC2 (%.1f%%)", var_explained[2]),
       title = "PCA: PC1 vs PC2") +
  theme_bw()

# PC2 x PC3
p_pc23 <- ggplot(pca_df, aes(x = PC2, y = PC3, color = group)) +
  geom_point(size = 3) +
  geom_text_repel(aes(label = sample_id), size = 2, max.overlaps = 10) +
  scale_color_manual(values = group_colors) +
  labs(x = sprintf("PC2 (%.1f%%)", var_explained[2]),
       y = sprintf("PC3 (%.1f%%)", var_explained[3]),
       title = "PCA: PC2 vs PC3") +
  theme_bw()

p_pca <- p_scree / (p_pc12 | p_pc23) +
  plot_annotation(title = "PCA — Sample Clustering (Imputed)")
ggsave(file.path(REPORT_DIR, "01_pca.pdf"), p_pca, width = 14, height = 10)

cat("   Saved: 01_pca.pdf\n")

# ==============================================================================
# 3.  HIERARCHICAL CLUSTERING
# ==============================================================================

cat("\n>> 3 -- Hierarchical Clustering\n")

dist_samples <- dist(mat_t, method = "euclidean")
hc_samples   <- hclust(dist_samples, method = "ward.D2")

dend <- as.dendrogram(hc_samples)
group_vec <- meta$group[match(labels(dend), meta$sample_id)]
dend <- color_branches(dend, clusters = as.numeric(factor(group_vec)),
                        col = group_colors[levels(meta$group)])

pdf(file.path(REPORT_DIR, "02_dendrogram.pdf"), width = 14, height = 6)
par(mar = c(8, 4, 3, 1))
plot(dend, main = "Sample Dendrogram (Euclidean, Ward.D2, Imputed)")
legend("topright", legend = levels(meta$group), fill = group_colors, cex = 0.8)
dev.off()

cat("   Saved: 02_dendrogram.pdf\n")

# ==============================================================================
# 4.  t-SNE
# ==============================================================================

cat("\n>> 4 -- t-SNE\n")

perplexities <- c(5, 10, 15, 30)
tsne_plots <- list()

for (perp in perplexities) {
  if (perp >= nrow(mat_t)) next
  set.seed(42)
  tsne_res <- Rtsne(mat_t, perplexity = perp, dims = 2, check_duplicates = FALSE)
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
  plot_annotation(title = "t-SNE — Sample Clustering (Imputed)")
ggsave(file.path(REPORT_DIR, "03_tsne.pdf"), p_tsne, width = 10, height = 8)

cat("   Saved: 03_tsne.pdf\n")

# ==============================================================================
# 5.  UMAP
# ==============================================================================

cat("\n>> 5 -- UMAP\n")

# Default UMAP
set.seed(42)
umap_default <- umap(mat_t)
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

# Tuned UMAP (n_neighbors = 10, min_dist = 0.01)
umap_config <- umap.defaults
umap_config$n_neighbors <- 10
umap_config$min_dist <- 0.01
set.seed(42)
umap_tuned <- umap(mat_t, config = umap_config)
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
p_umap <- p_umap + plot_annotation(title = "UMAP — Sample Clustering (Imputed)")
ggsave(file.path(REPORT_DIR, "04_umap.pdf"), p_umap, width = 14, height = 6)

cat("   Saved: 04_umap.pdf\n")

# ==============================================================================
# 6.  OUTLIER DETECTION
# ==============================================================================

cat("\n>> 6 -- Outlier Detection\n")

# Mahalanobis distance in PC space (first 10 PCs)
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
# 7.  SAMPLE CORRELATION HEATMAP
# ==============================================================================

cat("\n>> 7 -- Sample Correlation Heatmap\n")

cor_samples <- cor(mat, method = "pearson")

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
         main = "Sample Pearson correlation (Imputed)",
         fontsize = 8)
dev.off()

cat("   Saved: 05_sample_correlation.pdf\n")

# ==============================================================================
# 8.  EXPORT
# ==============================================================================

cat("\n>> 8 -- Export\n")

# PCA scores
pca_export <- pca_df %>%
  dplyr::select(sample_id, PC1:PC10, age, time, group)
write_csv(pca_export, file.path(DATA_DIR, "sample_pca_scores.csv"))

# Outlier flags
write_csv(outlier_df, file.path(DATA_DIR, "sample_outlier_flags.csv"))

cat("   Saved: sample_pca_scores.csv, sample_outlier_flags.csv\n")

cat("\n")
sessionInfo()

cat("\n=== YvO Sample Clustering (Imputed) Complete ===\n")
