################################################################################
# GO:BP hierarchy-based 20-category classification
#
# Strategy: Use depth-2 children of "biological_process" as the natural
# top-level partitioning, then split the largest ones at depth-3 to reach ~20.
# This is structure-driven (from the ontology), not data-driven (no bias).
################################################################################

setwd("/Users/dtl0018/Desktop/A_Proteomics_Analysis/A_YvO_2025/")

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(GO.db)
  library(org.Hs.eg.db)
  library(AnnotationDbi)
})

source("04_Figures/F2/a_script/YvO_F2_setup.R")

# --- Gene lists ---
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

# --- Get ALL GO:BP annotations for universe ---
suppressMessages({
  univ_entrez <- AnnotationDbi::mapIds(org.Hs.eg.db, keys = universe,
                        keytype = "SYMBOL", column = "ENTREZID",
                        multiVals = "first")
  univ_go <- AnnotationDbi::select(org.Hs.eg.db,
               keys = as.character(na.omit(univ_entrez)),
               keytype = "ENTREZID",
               columns = c("SYMBOL", "GO", "ONTOLOGY"))
})

univ_bp <- univ_go %>%
  filter(ONTOLOGY == "BP", !is.na(GO)) %>%
  distinct(SYMBOL, GO)

cat(sprintf("Universe: %d genes, %d with GO:BP annotations\n",
            length(universe), n_distinct(univ_bp$SYMBOL)))

# --- Get GO:BP hierarchy ---
ancestors <- as.list(GOBPANCESTOR)
children  <- as.list(GOBPCHILDREN)
bp_root   <- "GO:0008150"  # biological_process

# Depth-2: direct children of biological_process
depth2 <- children[[bp_root]]
depth2 <- depth2[!is.na(depth2)]

cat(sprintf("\nDepth-2 children of biological_process: %d terms\n", length(depth2)))

# For each depth-2 term, count how many universe genes have it as an ancestor
# (gene annotated to any descendant of the depth-2 term)
gene_to_ancestors <- univ_bp %>%
  rowwise() %>%
  mutate(anc_list = list({
    ancs <- ancestors[[GO]]
    if (is.null(ancs)) character(0) else c(GO, ancs[!is.na(ancs)])
  })) %>%
  ungroup()

# Flatten: gene -> all ancestor terms
gene_anc_flat <- gene_to_ancestors %>%
  unnest(anc_list) %>%
  dplyr::select(SYMBOL, ancestor = anc_list) %>%
  distinct()

# Count per depth-2 term
depth2_counts <- gene_anc_flat %>%
  filter(ancestor %in% depth2) %>%
  count(ancestor, name = "n_univ") %>%
  arrange(desc(n_univ))

# Get term names
depth2_counts$name <- AnnotationDbi::Term(depth2_counts$ancestor)

# Also count F2/F3 foreground hits
f2_anc <- gene_anc_flat %>% filter(SYMBOL %in% f2_fg)
f3_anc <- gene_anc_flat %>% filter(SYMBOL %in% f3_fg)

f2_d2 <- f2_anc %>% filter(ancestor %in% depth2) %>% count(ancestor, name = "n_f2")
f3_d2 <- f3_anc %>% filter(ancestor %in% depth2) %>% count(ancestor, name = "n_f3")

depth2_counts <- depth2_counts %>%
  left_join(f2_d2, by = "ancestor") %>%
  left_join(f3_d2, by = "ancestor") %>%
  mutate(across(starts_with("n_f"), ~replace_na(.x, 0L)))

cat(sprintf("\n%-12s  %-45s  %6s  %5s  %5s\n",
            "GO_ID", "Depth-2 Term", "Univ", "F2", "F3"))
cat(paste(rep("-", 80), collapse = ""), "\n")
for (i in seq_len(nrow(depth2_counts))) {
  r <- depth2_counts[i, ]
  cat(sprintf("%-12s  %-45s  %6d  %5d  %5d\n",
              r$ancestor, substr(r$name, 1, 45), r$n_univ, r$n_f2, r$n_f3))
}

# --- Identify the top 3-4 largest depth-2 terms to split ---
cat("\n\n=== SPLITTING LARGE DEPTH-2 TERMS AT DEPTH-3 ===\n")

large_terms <- depth2_counts %>% filter(n_univ >= 500)
cat(sprintf("Terms with >=500 universe genes (to split): %s\n",
            paste(large_terms$name, collapse = ", ")))

# For each large term, get its depth-3 children and their counts
depth3_results <- list()

for (lt in large_terms$ancestor) {
  lt_children <- children[[lt]]
  lt_children <- lt_children[!is.na(lt_children)]
  
  d3_counts <- gene_anc_flat %>%
    filter(ancestor %in% lt_children) %>%
    count(ancestor, name = "n_univ") %>%
    arrange(desc(n_univ))
  
  d3_counts$name <- AnnotationDbi::Term(d3_counts$ancestor)
  d3_counts$parent <- lt
  d3_counts$parent_name <- AnnotationDbi::Term(lt)
  
  d3_f2 <- f2_anc %>% filter(ancestor %in% lt_children) %>% count(ancestor, name = "n_f2")
  d3_f3 <- f3_anc %>% filter(ancestor %in% lt_children) %>% count(ancestor, name = "n_f3")
  
  d3_counts <- d3_counts %>%
    left_join(d3_f2, by = "ancestor") %>%
    left_join(d3_f3, by = "ancestor") %>%
    mutate(across(starts_with("n_f"), ~replace_na(.x, 0L))) %>%
    filter(n_univ >= 20)  # only substantial terms
  
  depth3_results[[lt]] <- d3_counts
  
  cat(sprintf("\n  %s (%s) — %d depth-3 children with >=20 genes:\n",
              lt, AnnotationDbi::Term(lt), nrow(d3_counts)))
  for (j in seq_len(min(nrow(d3_counts), 15))) {
    r <- d3_counts[j, ]
    cat(sprintf("    %-12s  %-40s  %5d  %4d  %4d\n",
                r$ancestor, substr(r$name, 1, 40), r$n_univ, r$n_f2, r$n_f3))
  }
}

# --- Build the proposed 20-term set ---
# Start with depth-2 terms that are NOT too large
small_d2 <- depth2_counts %>%
  filter(n_univ < 500, n_univ >= 5) %>%
  dplyr::select(go_id = ancestor, name, n_univ, n_f2, n_f3) %>%
  mutate(source = "depth-2")

cat(sprintf("\n\n=== PROPOSED 20 CONSOLIDATED GO:BP CATEGORIES ===\n"))
cat(sprintf("(Depth-2 terms kept as-is if <500 genes; depth-3 splits for large ones)\n\n"))

# Print the small depth-2 terms
cat(sprintf("KEPT AT DEPTH-2 (%d terms):\n", nrow(small_d2)))
cat(sprintf("%-12s  %-45s  %6s  %5s  %5s\n", "GO_ID", "Term", "Univ", "F2", "F3"))
cat(paste(rep("-", 80), collapse = ""), "\n")
for (i in seq_len(nrow(small_d2))) {
  r <- small_d2[i, ]
  cat(sprintf("%-12s  %-45s  %6d  %5d  %5d\n",
              r$go_id, substr(r$name, 1, 45), r$n_univ, r$n_f2, r$n_f3))
}

cat(sprintf("\nSPLIT AT DEPTH-3 (from %d large terms):\n", nrow(large_terms)))
for (lt in large_terms$ancestor) {
  d3 <- depth3_results[[lt]] %>%
    filter(n_univ >= 30) %>%
    head(8)  # top children
  cat(sprintf("\n  From %s:\n", AnnotationDbi::Term(lt)))
  for (j in seq_len(nrow(d3))) {
    r <- d3[j, ]
    cat(sprintf("  %-12s  %-40s  %5d  %4d  %4d\n",
                r$ancestor, substr(r$name, 1, 40), r$n_univ, r$n_f2, r$n_f3))
  }
}

# --- Test 1:1 assignment with a candidate set ---
# Combine: small depth-2 + top depth-3 from large terms
candidate_terms <- small_d2$go_id

for (lt in large_terms$ancestor) {
  d3 <- depth3_results[[lt]] %>%
    filter(n_univ >= 30) %>%
    head(5)
  candidate_terms <- c(candidate_terms, d3$ancestor)
}

cat(sprintf("\n\nTOTAL CANDIDATE TERMS: %d\n", length(candidate_terms)))

# 1:1 assignment: for each gene, find all candidate ancestors, pick most specific
gene_candidates <- gene_anc_flat %>%
  filter(ancestor %in% candidate_terms) %>%
  left_join(tibble(ancestor = candidate_terms,
                   term_univ = sapply(candidate_terms, function(t)
                     sum(gene_anc_flat$ancestor == t))),
            by = "ancestor")

gene_best <- gene_candidates %>%
  arrange(term_univ) %>%  # most specific (smallest universe count) first
  group_by(SYMBOL) %>%
  slice_head(n = 1) %>%
  ungroup()

f2_mapped <- gene_best %>% filter(SYMBOL %in% f2_fg)
f3_mapped <- gene_best %>% filter(SYMBOL %in% f3_fg)

cat(sprintf("\n1:1 COVERAGE:\n"))
cat(sprintf("  F2: %d/%d (%.1f%%)\n", n_distinct(f2_mapped$SYMBOL), length(f2_fg),
            100 * n_distinct(f2_mapped$SYMBOL) / length(f2_fg)))
cat(sprintf("  F3: %d/%d (%.1f%%)\n", n_distinct(f3_mapped$SYMBOL), length(f3_fg),
            100 * n_distinct(f3_mapped$SYMBOL) / length(f3_fg)))

# Distribution
cand_names <- AnnotationDbi::Term(candidate_terms)
cat(sprintf("\nDISTRIBUTION (1:1 assignment):\n"))
cat(sprintf("%-12s  %-40s  %5s  %5s\n", "GO_ID", "Term", "F2", "F3"))
cat(paste(rep("-", 70), collapse = ""), "\n")

dist_df <- gene_best %>%
  mutate(in_f2 = SYMBOL %in% f2_fg, in_f3 = SYMBOL %in% f3_fg) %>%
  group_by(ancestor) %>%
  summarise(n_f2 = sum(in_f2), n_f3 = sum(in_f3), .groups = "drop") %>%
  arrange(desc(n_f2 + n_f3))

for (i in seq_len(nrow(dist_df))) {
  r <- dist_df[i, ]
  tname <- cand_names[r$ancestor]
  if (is.na(tname)) tname <- r$ancestor
  cat(sprintf("%-12s  %-40s  %5d  %5d\n",
              r$ancestor, substr(tname, 1, 40), r$n_f2, r$n_f3))
}

# Unmapped
f2_unmapped <- setdiff(f2_fg, gene_best$SYMBOL[gene_best$SYMBOL %in% f2_fg])
f3_unmapped <- setdiff(f3_fg, gene_best$SYMBOL[gene_best$SYMBOL %in% f3_fg])
cat(sprintf("\nUnmapped → Other:\n"))
cat(sprintf("  F2: %d genes\n  F3: %d genes\n", length(f2_unmapped), length(f3_unmapped)))
if (length(f2_unmapped) > 0) cat(sprintf("  F2: %s\n", paste(f2_unmapped, collapse = ", ")))
if (length(f3_unmapped) > 0) cat(sprintf("  F3: %s\n", paste(head(f3_unmapped, 30), collapse = ", ")))
