################################################################################
# Test: GO:BP → 20 consolidated categories via semantic similarity clustering
#
# Approach:
#   1. Full GO:BP ORA on union of F2+F3 foreground genes
#   2. Cluster significant terms by semantic similarity (Wang method)
#   3. Pick 20 representative terms (highest IC per cluster)
#   4. Assign each gene to the most specific of those 20 terms
#   5. Report coverage
################################################################################

setwd("/Users/dtl0018/Desktop/A_Proteomics_Analysis/A_YvO_2025/")

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(clusterProfiler)
  library(org.Hs.eg.db)
  library(GO.db)
  library(GOSemSim)
  library(AnnotationDbi)
})

# Check if rrvgo is available
has_rrvgo <- requireNamespace("rrvgo", quietly = TRUE)
if (has_rrvgo) library(rrvgo)

source("04_Figures/F2/a_script/YvO_F2_setup.R")

# --- Define gene lists ---
f2_fg <- dep_df %>%
  filter(!is.na(logFC_Training_Young), !is.na(logFC_Training_Old)) %>%
  filter(pi_score_Training_Young < 0.05 | pi_score_Training_Old < 0.05 |
         pi_score_Interaction < 0.05) %>%
  pull(gene) %>% unique()

f3_fg <- dep_df %>%
  filter(!is.na(logFC_Aging), !is.na(logFC_Training_Old)) %>%
  filter(pi_score_Aging < 0.05 | pi_score_Reversal < 0.05) %>%
  pull(gene) %>% unique()

combined_fg <- union(f2_fg, f3_fg)
universe <- dep_df$gene

cat(sprintf("F2: %d genes, F3: %d genes, Combined: %d genes, Universe: %d\n\n",
            length(f2_fg), length(f3_fg), length(combined_fg), length(universe)))

# --- Map to Entrez IDs ---
suppressMessages({
  fg_entrez <- AnnotationDbi::mapIds(org.Hs.eg.db, keys = combined_fg,
                      keytype = "SYMBOL", column = "ENTREZID",
                      multiVals = "first")
  univ_entrez <- AnnotationDbi::mapIds(org.Hs.eg.db, keys = universe,
                        keytype = "SYMBOL", column = "ENTREZID",
                        multiVals = "first")
})

fg_entrez <- na.omit(fg_entrez)
univ_entrez <- na.omit(univ_entrez)

cat(sprintf("Entrez mapped: %d/%d foreground, %d/%d universe\n\n",
            length(fg_entrez), length(combined_fg),
            length(univ_entrez), length(universe)))

# --- Run GO:BP ORA ---
cat("Running GO:BP ORA...\n")
ego <- enrichGO(
  gene     = as.character(fg_entrez),
  universe = as.character(univ_entrez),
  OrgDb    = org.Hs.eg.db,
  ont      = "BP",
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.05,
  qvalueCutoff  = 0.2,
  minGSSize = 5,
  maxGSSize = 500,   # exclude ultra-broad terms
  readable  = TRUE   # convert Entrez back to symbols
)

ego_df <- as.data.frame(ego)
cat(sprintf("Significant GO:BP terms (padj < 0.05, 5-500 genes): %d\n\n", nrow(ego_df)))

if (nrow(ego_df) < 20) {
  cat("WARNING: fewer than 20 significant terms. Relaxing thresholds...\n")
  ego <- enrichGO(
    gene     = as.character(fg_entrez),
    universe = as.character(univ_entrez),
    OrgDb    = org.Hs.eg.db,
    ont      = "BP",
    pAdjustMethod = "BH",
    pvalueCutoff  = 0.2,
    qvalueCutoff  = 0.5,
    minGSSize = 3,
    maxGSSize = 800,
    readable  = TRUE
  )
  ego_df <- as.data.frame(ego)
  cat(sprintf("After relaxing: %d terms\n\n", nrow(ego_df)))
}

# --- Cluster by semantic similarity ---
cat("Clustering by semantic similarity...\n")

if (has_rrvgo && nrow(ego_df) >= 20) {
  # Use rrvgo for semantic similarity reduction
  sim_matrix <- rrvgo::calculateSimMatrix(
    ego_df$ID,
    orgdb = "org.Hs.eg.db",
    ont = "BP",
    method = "Wang"
  )
  
  # Reduce to clusters — aim for ~20 by adjusting threshold
  # Try different thresholds to get close to 20
  for (thresh in seq(0.5, 0.9, by = 0.05)) {
    reduced <- rrvgo::reduceSimMatrix(
      sim_matrix,
      setNames(ego_df$p.adjust, ego_df$ID),
      threshold = thresh,
      orgdb = "org.Hs.eg.db"
    )
    n_parents <- length(unique(reduced$parentTerm))
    if (n_parents >= 18 && n_parents <= 22) break
  }
  
  cat(sprintf("  Threshold=%.2f → %d parent clusters\n", thresh, n_parents))
  
  # Get the top 20 parent terms (by cluster size, then significance)
  parent_summary <- reduced %>%
    group_by(parent, parentTerm) %>%
    summarise(
      n_children = n(),
      best_pval = min(score),
      .groups = "drop"
    ) %>%
    arrange(best_pval) %>%
    head(20)
  
  # For each parent cluster, get the representative term (the parent itself)
  rep_terms <- parent_summary$parent
  rep_names <- parent_summary$parentTerm
  
} else {
  # Fallback: use simplify() from clusterProfiler
  cat("  Using clusterProfiler::simplify() fallback...\n")
  
  # Simplify with decreasing cutoff until we get ~20 terms
  for (cutoff in seq(0.5, 0.9, by = 0.05)) {
    ego_simp <- simplify(ego, cutoff = cutoff, by = "p.adjust", select_fun = min)
    n_terms <- nrow(as.data.frame(ego_simp))
    if (n_terms >= 18 && n_terms <= 25) break
  }
  cat(sprintf("  Cutoff=%.2f → %d terms\n", cutoff, n_terms))
  
  simp_df <- as.data.frame(ego_simp) %>%
    arrange(p.adjust) %>%
    head(20)
  
  rep_terms <- simp_df$ID
  rep_names <- simp_df$Description
}

cat(sprintf("\n=== 20 CONSOLIDATED GO:BP CATEGORIES ===\n\n"))
cat(sprintf("%-12s  %-50s  %6s  %6s\n", "GO_ID", "Term", "padj", "Size"))
cat(paste(rep("-", 85), collapse = ""), "\n")

for (i in seq_along(rep_terms)) {
  term_row <- ego_df %>% filter(ID == rep_terms[i])
  if (nrow(term_row) > 0) {
    cat(sprintf("%-12s  %-50s  %.1e  %d/%d\n",
                rep_terms[i],
                substr(term_row$Description[1], 1, 50),
                term_row$p.adjust[1],
                term_row$Count[1],
                as.integer(sub("/.*", "", term_row$GeneRatio[1]))))
  }
}

# --- Now assign each gene to the most specific of these 20 terms ---
cat("\n\n=== GENE ASSIGNMENT (1:1, most specific term wins) ===\n\n")

# Get all GO:BP annotations for foreground genes
suppressMessages({
  all_go <- AnnotationDbi::select(org.Hs.eg.db,
               keys = as.character(fg_entrez),
               keytype = "ENTREZID",
               columns = c("SYMBOL", "GO", "ONTOLOGY"))
})

gene_bp <- all_go %>%
  filter(ONTOLOGY == "BP", !is.na(GO)) %>%
  distinct(SYMBOL, GO)

# Map each gene's GO terms to ancestors, find which of the 20 rep terms they match
ancestors <- as.list(GOBPANCESTOR)

gene_to_rep <- gene_bp %>%
  rowwise() %>%
  mutate(rep_hits = list({
    ancs <- c(GO, ancestors[[GO]])
    ancs <- ancs[!is.na(ancs)]
    intersect(ancs, rep_terms)
  })) %>%
  ungroup() %>%
  unnest(rep_hits) %>%
  select(SYMBOL, rep_term = rep_hits) %>%
  distinct()

# Information content for specificity: use term size (smaller = more specific)
term_sizes <- ego_df %>%
  filter(ID %in% rep_terms) %>%
  select(ID, Count) %>%
  deframe()

# 1:1 assignment: pick the most specific (smallest Count) rep term per gene
gene_best <- gene_to_rep %>%
  mutate(term_size = term_sizes[rep_term]) %>%
  filter(!is.na(term_size)) %>%
  arrange(term_size) %>%
  group_by(SYMBOL) %>%
  slice_head(n = 1) %>%
  ungroup()

# Map back to gene symbols for F2/F3 foreground
f2_assigned <- gene_best %>% filter(SYMBOL %in% f2_fg)
f3_assigned <- gene_best %>% filter(SYMBOL %in% f3_fg)

# Coverage
f2_covered <- n_distinct(f2_assigned$SYMBOL)
f3_covered <- n_distinct(f3_assigned$SYMBOL)
combined_covered <- n_distinct(gene_best$SYMBOL)

cat(sprintf("Coverage:\n"))
cat(sprintf("  F2: %d/%d (%.1f%%)\n", f2_covered, length(f2_fg),
            100 * f2_covered / length(f2_fg)))
cat(sprintf("  F3: %d/%d (%.1f%%)\n", f3_covered, length(f3_fg),
            100 * f3_covered / length(f3_fg)))
cat(sprintf("  Combined: %d/%d (%.1f%%)\n", combined_covered, length(combined_fg),
            100 * combined_covered / length(combined_fg)))

# Category distribution
cat(sprintf("\n\nCATEGORY DISTRIBUTION:\n"))
cat(sprintf("%-50s  %5s  %5s  %5s\n", "Category", "F2", "F3", "Both"))
cat(paste(rep("-", 70), collapse = ""), "\n")

term_name_map <- setNames(ego_df$Description[match(rep_terms, ego_df$ID)], rep_terms)

for (rt in rep_terms) {
  n2 <- sum(f2_assigned$rep_term == rt, na.rm = TRUE)
  n3 <- sum(f3_assigned$rep_term == rt, na.rm = TRUE)
  nb <- sum(gene_best$rep_term == rt, na.rm = TRUE)
  tname <- substr(term_name_map[rt], 1, 50)
  if (nb > 0)
    cat(sprintf("%-50s  %5d  %5d  %5d\n", tname, n2, n3, nb))
}

# Unmapped genes
f2_unmapped <- setdiff(f2_fg, f2_assigned$SYMBOL)
f3_unmapped <- setdiff(f3_fg, f3_assigned$SYMBOL)
cat(sprintf("\nUnmapped genes:\n  F2: %d — %s\n  F3: %d — %s\n",
            length(f2_unmapped), paste(head(f2_unmapped, 20), collapse = ", "),
            length(f3_unmapped), paste(head(f3_unmapped, 20), collapse = ", ")))

# --- Also test: what if we use the universe GO:BP to define categories? ---
cat("\n\n=== ALTERNATIVE: Top 20 GO:BP terms by gene count in universe (size 20-400) ===\n")

suppressMessages({
  univ_go <- AnnotationDbi::select(org.Hs.eg.db,
                keys = as.character(univ_entrez),
                keytype = "ENTREZID",
                columns = c("SYMBOL", "GO", "ONTOLOGY"))
})

univ_bp <- univ_go %>%
  filter(ONTOLOGY == "BP", !is.na(GO)) %>%
  distinct(SYMBOL, GO)

# Count genes per GO:BP term in universe
term_counts_univ <- univ_bp %>%
  count(GO, name = "n_univ") %>%
  filter(n_univ >= 20, n_univ <= 400)

# Get term names
term_names_all <- AnnotationDbi::Term(term_counts_univ$GO)
term_counts_univ$name <- term_names_all

# Count foreground genes per term
fg_term_counts <- gene_bp %>%
  filter(GO %in% term_counts_univ$GO) %>%
  count(GO, name = "n_fg")

term_counts_univ <- term_counts_univ %>%
  left_join(fg_term_counts, by = "GO") %>%
  mutate(n_fg = replace_na(n_fg, 0L),
         enrichment = (n_fg / length(combined_fg)) / (n_univ / length(universe)))

# Top 20 by enrichment (must have >=3 fg genes)
top_enriched <- term_counts_univ %>%
  filter(n_fg >= 3) %>%
  arrange(desc(enrichment)) %>%
  head(40)

cat(sprintf("\n%-12s  %-45s  %5s  %5s  %6s\n",
            "GO_ID", "Term", "Univ", "FG", "Enrich"))
cat(paste(rep("-", 85), collapse = ""), "\n")
for (i in seq_len(nrow(top_enriched))) {
  r <- top_enriched[i, ]
  cat(sprintf("%-12s  %-45s  %5d  %5d  %5.1fx\n",
              r$GO, substr(r$name, 1, 45), r$n_univ, r$n_fg, r$enrichment))
}
