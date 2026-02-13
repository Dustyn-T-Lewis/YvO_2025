################################################################################
#
#   YvO Protein Clustering — Imputed Data
#
#   Workflow:
#     1.  Setup & data loading
#     2.  Hierarchical clustering (1 - Pearson, Ward.D2)
#     3.  Optimal k (silhouette, gap statistic, elbow)
#     4.  K-means at optimal k
#     5.  Co-expression heatmap with dendrogram
#     6.  Cluster characterization (mean logFC per contrast)
#     7.  Export
#
#   Input: 02_Imputation/c_data/01_imputed.csv (2124 proteins x 62 samples)
#   DEP results: 03_DEP/Imputed/a_limma/c_data/per_contrast_results/*.csv
#
################################################################################

# ==============================================================================
# 1.  SETUP & DATA LOADING
# ==============================================================================

cat("=== YvO Protein Clustering (Imputed) ===\n\n")
cat(">> 1 -- Setup & Data Loading\n")

suppressPackageStartupMessages({
  library(tidyverse)
  library(pheatmap)
  library(factoextra)
  library(cluster)
  library(dendextend)
  library(patchwork)
  library(RColorBrewer)
})

base_dir   <- normalizePath(file.path(dirname(getwd()), "..", "..", ".."), mustWork = TRUE)
DATA_FILE  <- file.path(base_dir, "02_Imputation", "c_data", "01_imputed.csv")
DEP_DIR    <- file.path(base_dir, "03_DEP", "Imputed", "a_limma", "c_data", "per_contrast_results")
REPORT_DIR <- file.path(base_dir, "04_Clustering", "Imputed", "a_protein_clustering", "b_reports")
DATA_DIR   <- file.path(base_dir, "04_Clustering", "Imputed", "a_protein_clustering", "c_data")

dir.create(REPORT_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(DATA_DIR,   recursive = TRUE, showWarnings = FALSE)

df <- read_csv(DATA_FILE, show_col_types = FALSE)
cat(sprintf("   Loaded: %d proteins x %d columns\n", nrow(df), ncol(df)))

ann_cols   <- c("uniprot_id", "protein", "gene", "description")
ann        <- df[, ann_cols]
samp_names <- setdiff(names(df), ann_cols)
mat        <- as.matrix(df[, samp_names])
rownames(mat) <- ann$uniprot_id

cat(sprintf("   Intensity matrix: %d proteins x %d samples\n", nrow(mat), ncol(mat)))

# Parse metadata for group annotation
meta <- tibble(sample_id = samp_names) %>%
  mutate(
    prefix = str_extract(sample_id, "^[A-Z]+"),
    time   = str_extract(sample_id, "(Pre|Post)$"),
    age    = if_else(str_detect(prefix, "^O"), "Old", "Young"),
    group  = paste(age, time, sep = "_")
  )

# ==============================================================================
# 2.  HIERARCHICAL CLUSTERING
# ==============================================================================

cat("\n>> 2 -- Hierarchical Clustering\n")

# Pearson correlation distance
cor_mat  <- cor(t(mat), method = "pearson")
dist_mat <- as.dist(1 - cor_mat)

hc <- hclust(dist_mat, method = "ward.D2")
cat(sprintf("   Hierarchical clustering: %d proteins, Ward.D2\n", nrow(mat)))

# ==============================================================================
# 3.  OPTIMAL k
# ==============================================================================

cat("\n>> 3 -- Optimal k Selection\n")

# Silhouette method
p_sil <- fviz_nbclust(mat, kmeans, method = "silhouette", k.max = 15,
                       nstart = 25) +
  labs(title = "Silhouette method") +
  theme_bw()

# Elbow method
p_elbow <- fviz_nbclust(mat, kmeans, method = "wss", k.max = 15,
                          nstart = 25) +
  labs(title = "Elbow method") +
  theme_bw()

# Gap statistic
set.seed(42)
gap_stat <- clusGap(mat, FUN = kmeans, nstart = 25, K.max = 15, B = 50)
p_gap <- fviz_gap_stat(gap_stat) +
  labs(title = "Gap statistic") +
  theme_bw()

p_optk <- p_sil / p_elbow / p_gap +
  plot_annotation(title = "Optimal k selection — Protein Clustering (Imputed)")

ggsave(file.path(REPORT_DIR, "01_optimal_k.pdf"), p_optk, width = 8, height = 12)

# Determine optimal k from silhouette
sil_vals <- sapply(2:15, function(k) {
  km <- kmeans(mat, centers = k, nstart = 25)
  ss <- silhouette(km$cluster, dist_mat)
  mean(ss[, 3])
})
optimal_k <- which.max(sil_vals) + 1
cat(sprintf("   Optimal k (silhouette): %d\n", optimal_k))

# ==============================================================================
# 4.  K-MEANS AT OPTIMAL k
# ==============================================================================

cat("\n>> 4 -- K-means Clustering\n")

set.seed(42)
km_fit <- kmeans(mat, centers = optimal_k, nstart = 50, iter.max = 100)

cluster_assignments <- tibble(
  uniprot_id = rownames(mat),
  cluster    = km_fit$cluster
) %>% left_join(ann %>% dplyr::select(uniprot_id, gene), by = "uniprot_id")

cat(sprintf("   K-means: k=%d, tot.withinss=%.1f\n", optimal_k, km_fit$tot.withinss))
cat("   Cluster sizes:\n")
print(table(km_fit$cluster))

# Cluster visualization
p_clust <- fviz_cluster(km_fit, data = mat, geom = "point", ellipse.type = "convex",
                         palette = brewer.pal(min(optimal_k, 8), "Set2"),
                         ggtheme = theme_bw()) +
  labs(title = sprintf("K-means clustering (k=%d, Imputed)", optimal_k))

ggsave(file.path(REPORT_DIR, "02_kmeans_pca.pdf"), p_clust, width = 8, height = 6)

# ==============================================================================
# 5.  CO-EXPRESSION HEATMAP
# ==============================================================================

cat("\n>> 5 -- Co-expression Heatmap\n")

# Annotation for columns (samples)
col_ann <- data.frame(
  Group = meta$group,
  row.names = meta$sample_id
)

# Annotation for rows (proteins) - cluster assignment
row_ann <- data.frame(
  Cluster = factor(km_fit$cluster),
  row.names = rownames(mat)
)

# Scale rows for visualization
mat_scaled <- t(scale(t(mat)))

# Select top variable proteins for readable heatmap
row_vars <- apply(mat, 1, var)
top_var_idx <- order(row_vars, decreasing = TRUE)[1:min(200, nrow(mat))]

pdf(file.path(REPORT_DIR, "03_coexpression_heatmap.pdf"), width = 12, height = 14)
pheatmap(mat_scaled[top_var_idx, ],
         clustering_distance_rows = "correlation",
         clustering_method = "ward.D2",
         annotation_col = col_ann,
         annotation_row = row_ann[top_var_idx, , drop = FALSE],
         show_rownames = FALSE,
         show_colnames = FALSE,
         main = "Co-expression heatmap (top 200 variable proteins, Imputed)",
         fontsize = 8)
dev.off()

cat("   Saved: 03_coexpression_heatmap.pdf\n")

# ==============================================================================
# 6.  CLUSTER CHARACTERIZATION
# ==============================================================================

cat("\n>> 6 -- Cluster Characterization\n")

contrast_names <- c("Training_Young", "Training_Old", "Aging", "Interaction")

char_rows <- list()
for (cname in contrast_names) {
  dep_path <- file.path(DEP_DIR, paste0(cname, ".csv"))
  if (!file.exists(dep_path)) next
  dep <- read_csv(dep_path, show_col_types = FALSE)

  dep_clust <- dep %>%
    inner_join(cluster_assignments %>% dplyr::select(uniprot_id, cluster),
               by = "uniprot_id")

  summary <- dep_clust %>%
    group_by(cluster) %>%
    summarise(
      mean_logFC = mean(logFC, na.rm = TRUE),
      median_logFC = median(logFC, na.rm = TRUE),
      n_sig_FDR05 = sum(adj.P.Val < 0.05, na.rm = TRUE),
      n_total = n(),
      .groups = "drop"
    ) %>%
    mutate(contrast = cname)

  char_rows <- c(char_rows, list(summary))
}

char_df <- bind_rows(char_rows)

p_char <- ggplot(char_df, aes(x = factor(cluster), y = mean_logFC, fill = contrast)) +
  geom_col(position = position_dodge(width = 0.8), width = 0.7) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  labs(x = "Cluster", y = "Mean log2FC", fill = "Contrast",
       title = "Cluster characterization by DEP logFC (Imputed)") +
  theme_bw(base_size = 11)

ggsave(file.path(REPORT_DIR, "04_cluster_characterization.pdf"), p_char,
       width = 10, height = 5)

cat("   Saved: 04_cluster_characterization.pdf\n")

# ==============================================================================
# 7.  EXPORT
# ==============================================================================

cat("\n>> 7 -- Export\n")

write_csv(cluster_assignments, file.path(DATA_DIR, "protein_clusters.csv"))
write_csv(char_df, file.path(DATA_DIR, "cluster_characterization.csv"))

cat("   Saved: protein_clusters.csv, cluster_characterization.csv\n")

cat("\n")
sessionInfo()

cat("\n=== YvO Protein Clustering (Imputed) Complete ===\n")
