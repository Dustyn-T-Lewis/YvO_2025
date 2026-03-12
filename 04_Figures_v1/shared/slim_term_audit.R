################################################################################
# Audit: Which GO Slim BP terms actually have hits in F2/F3 foregrounds?
# Shows raw term-level mapping before any super-category grouping.
################################################################################

setwd("/Users/dtl0018/Desktop/A_Proteomics_Analysis/A_YvO_2025/")

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(GO.db)
  library(org.Hs.eg.db)
  library(AnnotationDbi)
})

source("04_Figures/shared/go_slim_categories.R")
source("04_Figures/F2/a_script/YvO_F2_setup.R")

# --- Define foregrounds ---
f2_fg <- dep_df %>%
  filter(!is.na(logFC_Training_Young), !is.na(logFC_Training_Old)) %>%
  filter(pi_score_Training_Young < 0.05 | pi_score_Training_Old < 0.05 |
         pi_score_Interaction < 0.05) %>%
  pull(gene) %>% unique()

f3_fg <- dep_df %>%
  filter(!is.na(logFC_Aging), !is.na(logFC_Training_Old)) %>%
  filter(pi_score_Aging < 0.05 | pi_score_Reversal < 0.05) %>%
  pull(gene) %>% unique()

universe <- dep_df$gene

cat(sprintf("F2 foreground: %d genes\nF3 foreground: %d genes\nUniverse: %d genes\n\n",
            length(f2_fg), length(f3_fg), length(universe)))

# --- Map ALL universe genes to slim terms (reuse infrastructure) ---
suppressMessages({
  all_entrez <- AnnotationDbi::mapIds(org.Hs.eg.db, keys = universe,
                      keytype = "SYMBOL", column = "ENTREZID",
                      multiVals = "first")
  all_go <- AnnotationDbi::select(org.Hs.eg.db,
               keys = as.character(na.omit(all_entrez)),
               keytype = "ENTREZID",
               columns = c("SYMBOL", "GO", "ONTOLOGY"))
})

all_bp <- all_go %>%
  filter(ONTOLOGY == "BP", !is.na(GO)) %>%
  distinct(SYMBOL, GO)

ancestors <- as.list(GOBPANCESTOR)
all_go_ids <- unique(all_bp$GO)

go_to_slim <- setNames(
  lapply(all_go_ids, function(go_id) {
    hits <- character(0)
    if (go_id %in% bp_slim) hits <- go_id
    anc <- ancestors[[go_id]]
    if (!is.null(anc)) hits <- c(hits, intersect(anc, bp_slim))
    unique(hits)
  }),
  all_go_ids
)

# Gene -> slim term (many-to-many before 1:1 assignment)
gene_slim <- all_bp %>%
  mutate(slim_list = go_to_slim[GO]) %>%
  unnest(slim_list) %>%
  select(SYMBOL, slim = slim_list) %>%
  distinct()

# --- Get GO term names ---
slim_names <- AnnotationDbi::Term(bp_slim)
slim_name_df <- tibble(slim = bp_slim, term_name = slim_names)

# --- Count genes per slim term for each set ---
univ_counts <- gene_slim %>%
  count(slim, name = "n_universe")

f2_counts <- gene_slim %>%
  filter(SYMBOL %in% f2_fg) %>%
  count(slim, name = "n_f2")

f3_counts <- gene_slim %>%
  filter(SYMBOL %in% f3_fg) %>%
  count(slim, name = "n_f3")

# --- Combine ---
slim_summary <- slim_name_df %>%
  left_join(univ_counts, by = "slim") %>%
  left_join(f2_counts, by = "slim") %>%
  left_join(f3_counts, by = "slim") %>%
  mutate(across(starts_with("n_"), ~replace_na(.x, 0L))) %>%
  # Add current super-category mapping
  mutate(current_super = ifelse(slim %in% names(SLIM_SUPER),
                                 SLIM_SUPER[slim], "(unmapped)")) %>%
  arrange(desc(n_f3)) %>%
  filter(n_f2 > 0 | n_f3 > 0)  # only show terms with hits

cat("=" %>% rep(120) %>% paste(collapse = ""), "\n")
cat("GO SLIM BP TERMS WITH HITS IN F2 OR F3 FOREGROUNDS\n")
cat("(many-to-many: a gene can map to multiple slim terms before 1:1 assignment)\n")
cat("=" %>% rep(120) %>% paste(collapse = ""), "\n\n")

# Print as formatted table
cat(sprintf("%-12s  %-45s  %5s  %5s  %5s  %-25s\n",
            "GO_ID", "Term Name", "Univ", "F2", "F3", "Current Super"))
cat(paste(rep("-", 120), collapse = ""), "\n")

for (i in seq_len(nrow(slim_summary))) {
  r <- slim_summary[i, ]
  cat(sprintf("%-12s  %-45s  %5d  %5d  %5d  %-25s\n",
              r$slim,
              substr(r$term_name, 1, 45),
              r$n_universe, r$n_f2, r$n_f3,
              r$current_super))
}

cat(sprintf("\nTotal slim terms with hits: %d / 62\n", nrow(slim_summary)))

# --- Also show unmapped slim terms (in bp_slim but NOT in SLIM_SUPER) ---
unmapped_slims <- setdiff(bp_slim, names(SLIM_SUPER))
unmapped_with_hits <- slim_summary %>% filter(slim %in% unmapped_slims)

cat(sprintf("\n\nSLIM TERMS NOT MAPPED TO ANY SUPER-CATEGORY (%d total, %d with hits):\n",
            length(unmapped_slims), nrow(unmapped_with_hits)))
if (nrow(unmapped_with_hits) > 0) {
  for (i in seq_len(nrow(unmapped_with_hits))) {
    r <- unmapped_with_hits[i, ]
    cat(sprintf("  %-12s  %-45s  F2=%d  F3=%d\n",
                r$slim, r$term_name, r$n_f2, r$n_f3))
  }
}

# --- Coverage at 1:1 level (using assign_go_slim_super) ---
cat("\n\n1:1 ASSIGNMENT COVERAGE (after specificity weighting):\n")
f2_assigned <- assign_go_slim_super(f2_fg, universe)
f3_assigned <- assign_go_slim_super(f3_fg, universe)

f2_other <- sum(f2_assigned$super == "Other", na.rm = TRUE)
f3_other <- sum(f3_assigned$super == "Other", na.rm = TRUE)

cat(sprintf("  F2: %d/%d mapped to a super-category, %d in Other (%.1f%%)\n",
            nrow(f2_assigned) - f2_other, nrow(f2_assigned),
            f2_other, 100 * f2_other / nrow(f2_assigned)))
cat(sprintf("  F3: %d/%d mapped to a super-category, %d in Other (%.1f%%)\n",
            nrow(f3_assigned) - f3_other, nrow(f3_assigned),
            f3_other, 100 * f3_other / nrow(f3_assigned)))

# Show which genes are in "Other" and what slim terms they DO have
cat("\n\nGENES IN 'OTHER' — what slim terms do they have?\n")
other_genes_f2 <- f2_assigned %>% filter(super == "Other") %>% pull(gene)
other_genes_f3 <- f3_assigned %>% filter(super == "Other") %>% pull(gene)
other_genes <- union(other_genes_f2, other_genes_f3)

if (length(other_genes) > 0) {
  other_slims <- gene_slim %>%
    filter(SYMBOL %in% other_genes) %>%
    left_join(slim_name_df, by = "slim") %>%
    mutate(in_super = ifelse(slim %in% names(SLIM_SUPER), "YES", "no"))
  
  for (g in other_genes) {
    g_slims <- other_slims %>% filter(SYMBOL == g)
    if (nrow(g_slims) > 0) {
      terms_str <- paste(sprintf("%s (%s) [%s]", g_slims$slim,
                                  substr(g_slims$term_name, 1, 30),
                                  g_slims$in_super), collapse = "; ")
    } else {
      terms_str <- "(no GO:BP slim annotation)"
    }
    in_f2 <- ifelse(g %in% other_genes_f2, "F2", "")
    in_f3 <- ifelse(g %in% other_genes_f3, "F3", "")
    cat(sprintf("  %-10s [%s%s]: %s\n", g, in_f2, in_f3, terms_str))
  }
}

# --- Show per-super-category counts for both ---
cat("\n\nCURRENT SUPER-CATEGORY DISTRIBUTION:\n")
cat(sprintf("%-30s  %5s  %5s\n", "Super-Category", "F2", "F3"))
cat(paste(rep("-", 45), collapse = ""), "\n")
all_supers <- union(levels(f2_assigned$super), levels(f3_assigned$super))
for (s in all_supers) {
  n2 <- sum(f2_assigned$super == s, na.rm = TRUE)
  n3 <- sum(f3_assigned$super == s, na.rm = TRUE)
  if (n2 > 0 | n3 > 0) {
    cat(sprintf("%-30s  %5d  %5d\n", s, n2, n3))
  }
}
