################################################################################
#   Compare Pathway Classifications: GO Slim (B) vs Reactome (A) vs Kappa (C)
#
#   Loads classification outputs from all three approaches, computes comparison
#   metrics, and produces:
#     - Summary comparison table (CSV)
#     - Cross-tabulation (pairwise agreement)
#     - Sanity-check heatmaps (category overlap)
#
#   Run this AFTER running each of the three classification scripts:
#     test_reactome_classification.R  (Option A)
#     go_slim_categories.R export     (Option B)
#     test_kappa_classification.R     (Option C)
#
#   Output:
#     data/comparison_summary.csv
#     data/cross_tabulation_f2.csv
#     data/cross_tabulation_f3.csv
################################################################################

suppressPackageStartupMessages({
  library(tidyverse)
  library(readr)
})

setwd(rprojroot::find_rstudio_root_file())

DATA_DIR <- "04_Figures_v1/shared/data"

# === 1. Load DEP data (foreground definitions) ================================

dep_df <- read_csv("03_DEP/c_data/03_combined_results.csv", show_col_types = FALSE)
all_genes <- unique(dep_df$gene)

f2_fg <- dep_df %>%
  filter(pi_score_Training_Young < 0.05 | pi_score_Training_Old < 0.05 |
         pi_score_Interaction < 0.05) %>%
  pull(gene) %>% unique()

f3_fg <- dep_df %>%
  filter(pi_score_Aging < 0.05 | pi_score_Reversal < 0.05) %>%
  pull(gene) %>% unique()

# === 2. Run GO Slim classification (Option B) =================================

# go_slim_categories.R loads AnnotationDbi which masks dplyr::select.
# Source it then re-attach dplyr to restore priority.
source("04_Figures_v1/shared/go_slim_categories.R")
library(dplyr)

message("\n=== Option B: GO Slim classification ===")
f2_goslim <- assign_go_slim_consolidated(f2_fg, all_genes) %>%
  transmute(gene, go_slim_category = as.character(consolidated))
f3_goslim <- assign_go_slim_consolidated(f3_fg, all_genes) %>%
  transmute(gene, go_slim_category = as.character(consolidated))

# Also export the mapping
export_slim_mapping(f2_fg, all_genes, file.path(DATA_DIR, "go_slim_f2"))
export_slim_mapping(f3_fg, all_genes, file.path(DATA_DIR, "go_slim_f3"))

# === 3. Load Reactome classification (Option A) ===============================

message("\n=== Option A: Reactome classification ===")
f2_react_file <- file.path(DATA_DIR, "reactome_classification_f2.csv")
f3_react_file <- file.path(DATA_DIR, "reactome_classification_f3.csv")

if (!file.exists(f2_react_file) || !file.exists(f3_react_file)) {
  message("Running Reactome classification...")
  source("04_Figures_v1/shared/test_reactome_classification.R")
}

f2_react <- read_csv(f2_react_file, show_col_types = FALSE) %>%
  dplyr::rename(reactome_category = category)
f3_react <- read_csv(f3_react_file, show_col_types = FALSE) %>%
  dplyr::rename(reactome_category = category)

# === 4. Load Kappa classification (Option C) ==================================

message("\n=== Option C: Kappa classification ===")
f2_kappa_file <- file.path(DATA_DIR, "kappa_classification_f2.csv")
f3_kappa_file <- file.path(DATA_DIR, "kappa_classification_f3.csv")

if (!file.exists(f2_kappa_file) || !file.exists(f3_kappa_file)) {
  message("Running Kappa classification...")
  source("04_Figures_v1/shared/test_kappa_classification.R")
}

f2_kappa <- read_csv(f2_kappa_file, show_col_types = FALSE) %>%
  dplyr::rename(kappa_category = cluster_label)
f3_kappa <- read_csv(f3_kappa_file, show_col_types = FALSE) %>%
  dplyr::rename(kappa_category = cluster_label) %>%
  dplyr::select(gene, kappa_category)

# === 5. Optional: Load PathCards (Option D) ===================================

f2_pathcards <- NULL
f3_pathcards <- NULL
pathcards_file <- file.path(DATA_DIR, "pathcards_manual.csv")
if (file.exists(pathcards_file)) {
  message("Loading PathCards manual classification (Option D)...")
  pathcards <- read_csv(pathcards_file, show_col_types = FALSE)
  f2_pathcards <- pathcards %>% filter(gene_symbol %in% f2_fg) %>%
    dplyr::rename(gene = gene_symbol, pathcards_category = superpath_name)
  f3_pathcards <- pathcards %>% filter(gene_symbol %in% f3_fg) %>%
    dplyr::rename(gene = gene_symbol, pathcards_category = superpath_name)
}

# === 6. Merge all classifications =============================================

merge_classifications <- function(goslim, react, kappa, pathcards = NULL) {
  merged <- goslim %>%
    full_join(react, by = "gene") %>%
    full_join(kappa %>% dplyr::select(gene, kappa_category), by = "gene")

  if (!is.null(pathcards)) {
    merged <- merged %>%
      full_join(pathcards, by = "gene")
  }

  merged
}

f2_merged <- merge_classifications(f2_goslim, f2_react, f2_kappa, f2_pathcards)
f3_merged <- merge_classifications(f3_goslim, f3_react, f3_kappa, f3_pathcards)

# === 7. Compute comparison metrics ============================================

compute_method_stats <- function(categories, method_name, figure_name) {
  cats <- categories[!is.na(categories)]
  n_total <- length(cats)
  n_mapped <- sum(cats != "Other")
  cat_table <- table(cats[cats != "Other"])

  tibble(
    figure = figure_name,
    method = method_name,
    n_total = n_total,
    n_mapped = n_mapped,
    pct_mapped = round(100 * n_mapped / n_total, 1),
    n_categories = length(cat_table),
    largest_category = if (length(cat_table) > 0) names(which.max(cat_table)) else NA_character_,
    largest_pct = if (length(cat_table) > 0) round(100 * max(cat_table) / n_total, 1) else NA_real_,
    smallest_cat_n = if (length(cat_table) > 0) min(cat_table) else NA_integer_,
    cats_under_5 = sum(cat_table < 5)
  )
}

summary_stats <- bind_rows(
  compute_method_stats(f2_merged$go_slim_category, "GO Slim (B)", "F2"),
  compute_method_stats(f2_merged$reactome_category, "Reactome (A)", "F2"),
  compute_method_stats(f2_merged$kappa_category, "Kappa Multi-DB (C)", "F2"),
  compute_method_stats(f3_merged$go_slim_category, "GO Slim (B)", "F3"),
  compute_method_stats(f3_merged$reactome_category, "Reactome (A)", "F3"),
  compute_method_stats(f3_merged$kappa_category, "Kappa Multi-DB (C)", "F3")
)

if (!is.null(f2_pathcards)) {
  summary_stats <- bind_rows(summary_stats,
    compute_method_stats(f2_merged$pathcards_category, "PathCards (D)", "F2"),
    compute_method_stats(f3_merged$pathcards_category, "PathCards (D)", "F3")
  )
}

write_csv(summary_stats, file.path(DATA_DIR, "comparison_summary.csv"))

# === 8. Cross-tabulations =====================================================

cross_tab <- function(merged, col_a, col_b, label_a, label_b) {
  df <- merged %>%
    filter(!is.na(.data[[col_a]]), !is.na(.data[[col_b]])) %>%
    count(.data[[col_a]], .data[[col_b]], name = "n_genes")

  df %>%
    dplyr::rename(!!label_a := !!col_a, !!label_b := !!col_b)
}

# F2 cross-tabs
f2_gs_react <- cross_tab(f2_merged, "go_slim_category", "reactome_category",
                          "GO_Slim", "Reactome")
f2_gs_kappa <- cross_tab(f2_merged, "go_slim_category", "kappa_category",
                          "GO_Slim", "Kappa")
f2_react_kappa <- cross_tab(f2_merged, "reactome_category", "kappa_category",
                             "Reactome", "Kappa")

f2_cross <- bind_rows(
  f2_gs_react %>% mutate(comparison = "GO_Slim vs Reactome"),
  f2_gs_kappa %>%
    dplyr::rename(GO_Slim = GO_Slim, Kappa = Kappa) %>%
    pivot_longer(c(GO_Slim, Kappa), names_to = "method", values_to = "category") %>%
    mutate(comparison = "GO_Slim vs Kappa")
)

# Save simplified cross-tab: per-gene all three labels
write_csv(f2_merged, file.path(DATA_DIR, "cross_tabulation_f2.csv"))
write_csv(f3_merged, file.path(DATA_DIR, "cross_tabulation_f3.csv"))

# === 9. Pairwise agreement (Jaccard-like) =====================================

# For each pair of methods: what fraction of proteins get a "semantically
# similar" category? We use a simple metric: for proteins classified by both
# methods as non-"Other", do they share any keyword?

keyword_overlap <- function(cat_a, cat_b) {
  # Simple: check if any word in cat_a appears in cat_b (case-insensitive)
  words_a <- tolower(unlist(strsplit(cat_a, "[ &,/]+")))
  words_b <- tolower(unlist(strsplit(cat_b, "[ &,/]+")))
  # Remove common stop words
  stop <- c("of", "the", "and", "in", "a", "to", "by", "with", "for", "from",
            "other", "st:", "regulation", "process")
  words_a <- setdiff(words_a, stop)
  words_b <- setdiff(words_b, stop)
  length(intersect(words_a, words_b)) > 0
}

pairwise_agreement <- function(merged, col_a, col_b) {
  both <- merged %>%
    filter(!is.na(.data[[col_a]]), !is.na(.data[[col_b]]),
           .data[[col_a]] != "Other", .data[[col_b]] != "Other")
  if (nrow(both) == 0) return(NA_real_)

  agrees <- mapply(keyword_overlap, both[[col_a]], both[[col_b]])
  round(100 * mean(agrees), 1)
}

agreement <- tibble(
  comparison = c("GO Slim vs Reactome", "GO Slim vs Kappa", "Reactome vs Kappa"),
  f2_agreement_pct = c(
    pairwise_agreement(f2_merged, "go_slim_category", "reactome_category"),
    pairwise_agreement(f2_merged, "go_slim_category", "kappa_category"),
    pairwise_agreement(f2_merged, "reactome_category", "kappa_category")
  ),
  f3_agreement_pct = c(
    pairwise_agreement(f3_merged, "go_slim_category", "reactome_category"),
    pairwise_agreement(f3_merged, "go_slim_category", "kappa_category"),
    pairwise_agreement(f3_merged, "reactome_category", "kappa_category")
  )
)
write_csv(agreement, file.path(DATA_DIR, "pairwise_agreement.csv"))

# === 10. Print summary ========================================================

message("\n")
message("=" %>% strrep(70))
message("PATHWAY CLASSIFICATION COMPARISON")
message("=" %>% strrep(70))

message("\n--- Summary Statistics ---")
print(as.data.frame(summary_stats), row.names = FALSE)

message("\n--- Pairwise Keyword Agreement ---")
print(as.data.frame(agreement), row.names = FALSE)

message("\n--- Files written ---")
message(sprintf("  %s/comparison_summary.csv", DATA_DIR))
message(sprintf("  %s/cross_tabulation_f2.csv", DATA_DIR))
message(sprintf("  %s/cross_tabulation_f3.csv", DATA_DIR))
message(sprintf("  %s/pairwise_agreement.csv", DATA_DIR))

message("\nReview cross_tabulation_*.csv for per-gene category assignments.")
message("Use comparison_summary.csv for the manuscript methods comparison table.")
