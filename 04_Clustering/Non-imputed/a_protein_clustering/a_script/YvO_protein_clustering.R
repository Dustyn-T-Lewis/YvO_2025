################################################################################
#
#   YvO Protein Clustering — Non-imputed Data
#
#   Workflow:
#     1.  Setup & data loading
#     2.  Handle missing values (na.omit)
#     3.  Hierarchical clustering (1 - Pearson, Ward.D2)
#     4.  Optimal k (silhouette, gap statistic, elbow)
#     5.  K-means at optimal k
#     6.  Co-expression heatmap with dendrogram
#     7.  Cluster characterization (mean logFC per contrast)
#     8.  Export
#
#   Input: 01_normalization/c_data/01_normalized.csv (2124 proteins, NAs present)
#   DEP results: 03_DEP/Non-imputed/a_limma/c_data/per_contrast_results/*.csv
#
################################################################################

# ==============================================================================
# 1.  SETUP & DATA LOADING
# ==============================================================================

cat("=== YvO Protein Clustering (Non-imputed) ===\n\n")
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
DATA_FILE  <- file.path(base_dir, "01_normalization", "c_data", "01_normalized.csv")
DEP_DIR    <- file.path(base_dir, "03_DEP", "Non-imputed", "a_limma", "c_data", "per_contrast_results")
REPORT_DIR <- file.path(base_dir, "04_Clustering", "Non-imputed", "a_protein_clustering", "b_reports")
DATA_DIR   <- file.path(base_dir, "04_Clustering", "Non-imputed", "a_protein_clustering", "c_data")

dir.create(REPORT_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(DATA_DIR,   recursive = TRUE, showWarnings = FALSE)

df <- read_csv(DATA_FILE, show_col_types = FALSE)
cat(sprintf("   Loaded: %d proteins x %d columns\n", nrow(df), ncol(df)))

ann_cols   <- c("uniprot_id", "protein", "gene", "description")
ann        <- df[, ann_cols]
samp_names <- setdiff(names(df), ann_cols)
mat        <- as.matrix(df[, samp_names])
rownames(mat) <- ann$uniprot_id

n_missing <- sum(is.na(mat))
cat(sprintf("   Intensity matrix: %d proteins x %d samples (%d NAs)\n",
            nrow(mat), ncol(mat), n_missing))

# Parse metadata for group annotation
meta <- tibble(sample_id = samp_names) %>%
  mutate(
    prefix = str_extract(sample_id, "^[A-Z]+"),
    time   = str_extract(sample_id, "(Pre|Post)$"),
    age    = if_else(str_detect(prefix, "^O"), "Old", "Young"),
    group  = paste(age, time, sep = "_")
  )

# ==============================================================================
# 2.  HANDLE MISSING VALUES
# ==============================================================================

cat("\n>> 2 -- Handle Missing Values\n")

n_before <- nrow(mat)
mat_complete <- na.omit(mat)
n_after <- nrow(mat_complete)
n_removed <- n_before - n_after

cat(sprintf("   Proteins before na.omit: %d\n", n_before))
cat(sprintf("   Proteins after na.omit:  %d\n", n_after))
cat(sprintf("   Proteins removed (had NAs): %d (%.1f%%)\n",
            n_removed, 100 * n_removed / n_before))

# Filter annotation to match
ann_complete <- ann %>% filter(uniprot_id %in% rownames(mat_complete))

# ==============================================================================
# 3.  HIERARCHICAL CLUSTERING
# ==============================================================================

cat("\n>> 3 -- Hierarchical Clustering\n")

cor_mat  <- cor(t(mat_complete), method = "pearson")
dist_mat <- as.dist(1 - cor_mat)

hc <- hclust(dist_mat, method = "ward.D2")
cat(sprintf("   Hierarchical clustering: %d proteins, Ward.D2\n", nrow(mat_complete)))

# ==============================================================================
# 4.  OPTIMAL k
# ==============================================================================

cat("\n>> 4 -- Optimal k Selection\n")

p_sil <- fviz_nbclust(mat_complete, kmeans, method = "silhouette", k.max = 15,
                       nstart = 25) +
  labs(title = "Silhouette method") +
  theme_bw()

p_elbow <- fviz_nbclust(mat_complete, kmeans, method = "wss", k.max = 15,
                          nstart = 25) +
  labs(title = "Elbow method") +
  theme_bw()

set.seed(42)
gap_stat <- clusGap(mat_complete, FUN = kmeans, nstart = 25, K.max = 15, B = 50)
p_gap <- fviz_gap_stat(gap_stat) +
  labs(title = "Gap statistic") +
  theme_bw()

p_optk <- p_sil / p_elbow / p_gap +
  plot_annotation(title = "Optimal k selection — Protein Clustering (Non-imputed)")

ggsave(file.path(REPORT_DIR, "01_optimal_k.pdf"), p_optk, width = 8, height = 12)

sil_vals <- sapply(2:15, function(k) {
  km <- kmeans(mat_complete, centers = k, nstart = 25)
  ss <- silhouette(km$cluster, dist_mat)
  mean(ss[, 3])
})
optimal_k <- which.max(sil_vals) + 1
cat(sprintf("   Optimal k (silhouette): %d\n", optimal_k))

# ==============================================================================
# 5.  K-MEANS AT OPTIMAL k
# ==============================================================================

cat("\n>> 5 -- K-means Clustering\n")

set.seed(42)
km_fit <- kmeans(mat_complete, centers = optimal_k, nstart = 50, iter.max = 100)

cluster_assignments <- tibble(
  uniprot_id = rownames(mat_complete),
  cluster    = km_fit$cluster
) %>% left_join(ann_complete %>% dplyr::select(uniprot_id, gene), by = "uniprot_id")

cat(sprintf("   K-means: k=%d, tot.withinss=%.1f\n", optimal_k, km_fit$tot.withinss))
cat("   Cluster sizes:\n")
print(table(km_fit$cluster))

p_clust <- fviz_cluster(km_fit, data = mat_complete, geom = "point",
                         ellipse.type = "convex",
                         palette = brewer.pal(min(optimal_k, 8), "Set2"),
                         ggtheme = theme_bw()) +
  labs(title = sprintf("K-means clustering (k=%d, Non-imputed, %d proteins)", optimal_k, n_after))

ggsave(file.path(REPORT_DIR, "02_kmeans_pca.pdf"), p_clust, width = 8, height = 6)

# ==============================================================================
# 6.  CO-EXPRESSION HEATMAP
# ==============================================================================

cat("\n>> 6 -- Co-expression Heatmap\n")

col_ann <- data.frame(Group = meta$group, row.names = meta$sample_id)
row_ann <- data.frame(Cluster = factor(km_fit$cluster), row.names = rownames(mat_complete))

mat_scaled <- t(scale(t(mat_complete)))
row_vars <- apply(mat_complete, 1, var)
top_var_idx <- order(row_vars, decreasing = TRUE)[1:min(200, nrow(mat_complete))]

pdf(file.path(REPORT_DIR, "03_coexpression_heatmap.pdf"), width = 12, height = 14)
pheatmap(mat_scaled[top_var_idx, ],
         clustering_distance_rows = "correlation",
         clustering_method = "ward.D2",
         annotation_col = col_ann,
         annotation_row = row_ann[top_var_idx, , drop = FALSE],
         show_rownames = FALSE,
         show_colnames = FALSE,
         main = sprintf("Co-expression heatmap (top 200 variable, Non-imputed, %d proteins)", n_after),
         fontsize = 8)
dev.off()

cat("   Saved: 03_coexpression_heatmap.pdf\n")

# ==============================================================================
# 7.  CLUSTER CHARACTERIZATION
# ==============================================================================

cat("\n>> 7 -- Cluster Characterization\n")

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
       title = "Cluster characterization by DEP logFC (Non-imputed)") +
  theme_bw(base_size = 11)

ggsave(file.path(REPORT_DIR, "04_cluster_characterization.pdf"), p_char,
       width = 10, height = 5)

# ==============================================================================
# 8.  EXPORT
# ==============================================================================

cat("\n>> 8 -- Export\n")

write_csv(cluster_assignments, file.path(DATA_DIR, "protein_clusters.csv"))
write_csv(char_df, file.path(DATA_DIR, "cluster_characterization.csv"))

cat("   Saved: protein_clusters.csv, cluster_characterization.csv\n")
cat(sprintf("   NOTE: %d proteins excluded due to missing values\n", n_removed))

cat("\n")
sessionInfo()

cat("\n=== YvO Protein Clustering (Non-imputed) Complete ===\n")
