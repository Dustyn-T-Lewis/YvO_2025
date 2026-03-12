################################################################################
#   Option C: Kappa-Based Multi-Database Pathway Classification
#
#   Replicates the Metascape kappa clustering approach (Zhou et al. 2019,
#   Nature Communications) programmatically in R:
#
#   1. Collect all MSigDB C2:CP gene sets (Reactome + KEGG + WikiPathways +
#      BioCarta + PID)
#   2. For foreground proteins, identify which gene sets they belong to
#   3. Build binary gene membership matrix (genes × terms)
#   4. Compute pairwise Cohen's kappa between terms (gene overlap statistic)
#   5. Hierarchical clustering (average linkage), cut at distance 0.7
#      (kappa >= 0.3)
#   6. Select most significant term per cluster as representative
#   7. Assign proteins to clusters via their most-specific member term
#   8. Label clusters from representative term names
#
#   Dependencies: msigdbr (CRAN), org.Hs.eg.db (Bioconductor)
#
#   Exports:
#     assign_kappa_classification(fg_genes, all_genes, min_cat_size = 2)
#       → tibble(gene, cluster_id, cluster_label)
#
#   Output files (when run standalone):
#     data/kappa_classification_f2.csv
#     data/kappa_classification_f3.csv
#     data/kappa_coverage_stats.csv
#     data/kappa_cluster_terms.csv
################################################################################

suppressPackageStartupMessages({
  library(msigdbr)
  library(org.Hs.eg.db)
  library(AnnotationDbi)
  library(dplyr)
  library(tidyr)
  library(readr)
  library(stats)
})

# === Cohen's kappa for binary vectors ========================================
# Vectorized computation of kappa between all pairs of columns in a binary matrix

compute_kappa_matrix <- function(bin_mat) {
  # bin_mat: genes (rows) x terms (columns), binary (0/1)
  n <- nrow(bin_mat)
  p <- ncol(bin_mat)

  # Pairwise agreement: both 1 (a), both 0 (d)
  # a = t(X) %*% X = number of genes in both terms
  a <- crossprod(bin_mat)
  # d = t(1-X) %*% (1-X) = genes in neither term
  d <- crossprod(1 - bin_mat)
  # b = t(X) %*% (1-X) = in term i but not j
  b <- crossprod(bin_mat, 1 - bin_mat)
  # c = t(1-X) %*% X = in term j but not i
  cc <- crossprod(1 - bin_mat, bin_mat)

  # Cohen's kappa = (p_o - p_e) / (1 - p_e)
  p_o <- (a + d) / n
  p_e <- ((a + b) * (a + cc) + (cc + d) * (b + d)) / n^2
  kappa <- (p_o - p_e) / (1 - p_e)

  # Fix NaN (when p_e == 1, both vectors identical)
  kappa[is.nan(kappa)] <- 1
  diag(kappa) <- 1

  kappa
}

# === Main classification function =============================================

assign_kappa_classification <- function(fg_genes, all_genes, min_cat_size = 2,
                                         kappa_threshold = 0.3,
                                         max_cluster_label_len = 50) {

  # 1. Get C2:CP gene sets from msigdbr
  c2cp <- msigdbr(species = "Homo sapiens", collection = "C2",
                   subcollection = "CP") %>%
    bind_rows(
      msigdbr(species = "Homo sapiens", collection = "C2", subcollection = "CP:REACTOME"),
      msigdbr(species = "Homo sapiens", collection = "C2", subcollection = "CP:KEGG_MEDICUS"),
      msigdbr(species = "Homo sapiens", collection = "C2", subcollection = "CP:KEGG_LEGACY"),
      msigdbr(species = "Homo sapiens", collection = "C2", subcollection = "CP:WIKIPATHWAYS"),
      msigdbr(species = "Homo sapiens", collection = "C2", subcollection = "CP:BIOCARTA"),
      msigdbr(species = "Homo sapiens", collection = "C2", subcollection = "CP:PID")
    ) %>%
    distinct(gs_name, gene_symbol) %>%
    filter(gene_symbol %in% all_genes)

  # Filter out disease/cancer-specific terms
  disease_patterns <- "DISEASE|CANCER|TUMOR|CARCINOMA|LEUKEMIA|LYMPHOMA|MELANOMA|GLIOMA|HEPATITIS|HIV|INFECTION|VIRAL|BACTERIAL"
  c2cp <- c2cp %>%
    filter(!grepl(disease_patterns, gs_name, ignore.case = TRUE))

  message(sprintf("C2:CP after disease filtering: %d terms, %d gene-term pairs",
                  n_distinct(c2cp$gs_name), nrow(c2cp)))

  # 2. Filter to terms that contain at least 2 foreground genes
  fg_terms <- c2cp %>%
    filter(gene_symbol %in% fg_genes) %>%
    count(gs_name, name = "n_fg") %>%
    filter(n_fg >= 2) %>%
    pull(gs_name)

  c2cp_fg <- c2cp %>% filter(gs_name %in% fg_terms)
  message(sprintf("Terms with >=2 foreground genes: %d", length(fg_terms)))

  if (length(fg_terms) < 3) {
    warning("Too few terms for kappa clustering")
    return(tibble(gene = fg_genes, cluster_id = 0L, cluster_label = "Other"))
  }

  # 3. Build binary gene membership matrix (universe genes × terms)
  # Use only genes in the universe that appear in at least one term
  relevant_genes <- unique(c2cp_fg$gene_symbol)
  term_list <- split(c2cp_fg$gene_symbol, c2cp_fg$gs_name)

  bin_mat <- matrix(0L, nrow = length(relevant_genes), ncol = length(term_list),
                    dimnames = list(relevant_genes, names(term_list)))
  for (term in names(term_list)) {
    bin_mat[term_list[[term]], term] <- 1L
  }

  message(sprintf("Binary matrix: %d genes × %d terms", nrow(bin_mat), ncol(bin_mat)))

  # 4. Compute pairwise kappa
  message("Computing pairwise kappa matrix...")
  kappa_mat <- compute_kappa_matrix(bin_mat)

  # 5. Hierarchical clustering on kappa distance
  kappa_dist <- as.dist(1 - kappa_mat)
  hc <- hclust(kappa_dist, method = "average")
  clusters <- cutree(hc, h = 1 - kappa_threshold)  # distance = 1 - kappa
  message(sprintf("Kappa clustering: %d clusters at threshold %.2f",
                  max(clusters), kappa_threshold))

  # 6. Select representative term per cluster (most fg genes)
  cluster_df <- tibble(
    gs_name = names(clusters),
    cluster_id = as.integer(clusters)
  )

  fg_term_counts <- c2cp_fg %>%
    filter(gene_symbol %in% fg_genes) %>%
    count(gs_name, name = "n_fg_genes")

  representatives <- cluster_df %>%
    left_join(fg_term_counts, by = "gs_name") %>%
    mutate(n_fg_genes = coalesce(n_fg_genes, 0L)) %>%
    group_by(cluster_id) %>%
    slice_max(n_fg_genes, n = 1, with_ties = FALSE) %>%
    ungroup() %>%
    mutate(
      # Clean up representative name for labeling
      cluster_label = gs_name %>%
        gsub("^(REACTOME_|KEGG_|WP_|BIOCARTA_|PID_)", "", .) %>%
        gsub("_", " ", .) %>%
        substr(1, max_cluster_label_len)
    )

  # 7. Assign foreground genes to clusters via most-specific member term
  fg_gene_terms <- c2cp_fg %>%
    filter(gene_symbol %in% fg_genes) %>%
    inner_join(cluster_df, by = "gs_name")

  # Count fg genes per cluster (for specificity weighting)
  cluster_sizes <- fg_gene_terms %>%
    distinct(gene_symbol, cluster_id) %>%
    count(cluster_id, name = "n_fg_cluster")

  best_assignment <- fg_gene_terms %>%
    distinct(gene_symbol, cluster_id) %>%
    left_join(cluster_sizes, by = "cluster_id") %>%
    arrange(n_fg_cluster) %>%  # prefer smaller (more specific) clusters
    group_by(gene_symbol) %>%
    slice_head(n = 1) %>%
    ungroup() %>%
    left_join(representatives %>% select(cluster_id, cluster_label), by = "cluster_id")

  # 8. Unmapped → "Other"
  unmapped <- setdiff(fg_genes, best_assignment$gene_symbol)
  if (length(unmapped) > 0) {
    best_assignment <- bind_rows(best_assignment,
      tibble(gene_symbol = unmapped, cluster_id = 0L,
             cluster_label = "Other", n_fg_cluster = NA_integer_))
  }

  # 9. Merge small clusters
  small_clusters <- best_assignment %>%
    count(cluster_label) %>%
    filter(n < min_cat_size, cluster_label != "Other") %>%
    pull(cluster_label)
  if (length(small_clusters) > 0) {
    best_assignment <- best_assignment %>%
      mutate(cluster_label = ifelse(cluster_label %in% small_clusters,
                                     "Other", cluster_label))
  }

  list(
    classification = best_assignment %>%
      transmute(gene = gene_symbol, cluster_id, cluster_label) %>%
      arrange(cluster_label, gene),
    representatives = representatives,
    cluster_df = cluster_df
  )
}

# === Standalone execution =====================================================

if (sys.nframe() == 0 || identical(environment(), globalenv())) {
  setwd(rprojroot::find_rstudio_root_file())

  message("\n=== Running Kappa multi-database classification (Option C) ===\n")

  dep_df <- read_csv("03_DEP/c_data/03_combined_results.csv", show_col_types = FALSE)
  all_genes <- unique(dep_df$gene)

  # F2 foreground
  f2_fg <- dep_df %>%
    filter(pi_score_Training_Young < 0.05 | pi_score_Training_Old < 0.05 |
           pi_score_Interaction < 0.05) %>%
    pull(gene) %>% unique()

  # F3 foreground
  f3_fg <- dep_df %>%
    filter(pi_score_Aging < 0.05 | pi_score_Reversal < 0.05) %>%
    pull(gene) %>% unique()

  message(sprintf("F2 foreground: %d proteins", length(f2_fg)))
  message(sprintf("F3 foreground: %d proteins", length(f3_fg)))

  # Classify
  f2_kappa <- assign_kappa_classification(f2_fg, all_genes)
  f3_kappa <- assign_kappa_classification(f3_fg, all_genes)

  # Save
  out_dir <- "04_Figures_v1/shared/data"
  write_csv(f2_kappa$classification, file.path(out_dir, "kappa_classification_f2.csv"))
  write_csv(f3_kappa$classification, file.path(out_dir, "kappa_classification_f3.csv"))

  # Save cluster term mappings
  write_csv(
    bind_rows(
      f2_kappa$representatives %>% mutate(figure = "F2"),
      f3_kappa$representatives %>% mutate(figure = "F3")
    ),
    file.path(out_dir, "kappa_cluster_terms.csv")
  )

  # Coverage stats
  compute_kappa_stats <- function(cls, label) {
    df <- cls$classification
    n_total  <- nrow(df)
    n_mapped <- sum(df$cluster_label != "Other")
    cat_counts <- df %>% filter(cluster_label != "Other") %>% count(cluster_label)
    tibble(
      figure = label,
      n_total = n_total,
      n_mapped = n_mapped,
      pct_mapped = round(100 * n_mapped / n_total, 1),
      n_categories = nrow(cat_counts),
      largest_cat = if (nrow(cat_counts) > 0) cat_counts$cluster_label[which.max(cat_counts$n)] else NA_character_,
      largest_pct = if (nrow(cat_counts) > 0) round(100 * max(cat_counts$n) / n_total, 1) else NA_real_,
      smallest_cat_n = if (nrow(cat_counts) > 0) min(cat_counts$n) else NA_integer_
    )
  }

  stats <- bind_rows(
    compute_kappa_stats(f2_kappa, "F2"),
    compute_kappa_stats(f3_kappa, "F3")
  )
  write_csv(stats, file.path(out_dir, "kappa_coverage_stats.csv"))

  message("\n--- F2 Kappa clusters ---")
  print(table(f2_kappa$classification$cluster_label))
  message("\n--- F3 Kappa clusters ---")
  print(table(f3_kappa$classification$cluster_label))
  message(sprintf("\nF2 coverage: %.1f%% (%d/%d), %d clusters",
                  stats$pct_mapped[1], stats$n_mapped[1], stats$n_total[1],
                  stats$n_categories[1]))
  message(sprintf("F3 coverage: %.1f%% (%d/%d), %d clusters",
                  stats$pct_mapped[2], stats$n_mapped[2], stats$n_total[2],
                  stats$n_categories[2]))
}
