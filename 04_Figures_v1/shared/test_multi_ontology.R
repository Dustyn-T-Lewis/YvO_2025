################################################################################
#   Empirical test: does combining GO:BP + GO:MF + GO:CC inflate enrichment
#   compared to BP-only GO Slim classification?
#
#   Three approaches compared:
#     A) BP-only (current system)
#     B) BP-primary with CC/MF rescue for "Other" genes
#     C) Full multi-ontology merged categories
#
#   Uses combined F2+F3 foreground as test case.
################################################################################

suppressPackageStartupMessages({
  library(tidyverse)
  library(GO.db)
  library(org.Hs.eg.db)
  library(AnnotationDbi)
})

# === 1. DATA LOADING ==========================================================

setwd("/Users/dtl0018/Desktop/A_Proteomics_Analysis/A_YvO_2025/")
source("04_Figures/F2/a_script/YvO_F2_setup.R")
# go_slim_categories.R is already sourced by F2 setup, but let's be explicit:
# source("04_Figures/shared/go_slim_categories.R")

message("\n=== MULTI-ONTOLOGY ENRICHMENT TEST ===\n")

# === 2. DEFINE GENE LISTS =====================================================

universe <- dep_df$gene
message(sprintf("Universe: %d genes", length(universe)))

# F2 foreground: pi_score < 0.05 for any of the three F2 contrasts
fg_f2 <- dep_df %>%
  filter(pi_score_Training_Young < 0.05 |
         pi_score_Training_Old < 0.05 |
         pi_score_Interaction < 0.05) %>%
  pull(gene)
message(sprintf("F2 foreground: %d genes", length(fg_f2)))

# F3 foreground: pi_score < 0.05 for Aging or Reversal
fg_f3 <- dep_df %>%
  filter(pi_score_Aging < 0.05 |
         pi_score_Reversal < 0.05) %>%
  pull(gene)
message(sprintf("F3 foreground: %d genes", length(fg_f3)))

# Combined foreground
fg_combined <- union(fg_f2, fg_f3)
fg_genes <- fg_combined
message(sprintf("Combined foreground: %d genes", length(fg_genes)))

# === 3. GET GO ANNOTATIONS (ALL THREE ONTOLOGIES) =============================

message("\nFetching GO annotations...")
suppressMessages({
  all_entrez <- AnnotationDbi::mapIds(org.Hs.eg.db, keys = universe,
                      keytype = "SYMBOL", column = "ENTREZID",
                      multiVals = "first")
  all_go <- AnnotationDbi::select(org.Hs.eg.db,
               keys = as.character(na.omit(all_entrez)),
               keytype = "ENTREZID",
               columns = c("SYMBOL", "GO", "ONTOLOGY"))
})

all_bp <- all_go %>% filter(ONTOLOGY == "BP", !is.na(GO)) %>% distinct(SYMBOL, GO)
all_cc <- all_go %>% filter(ONTOLOGY == "CC", !is.na(GO)) %>% distinct(SYMBOL, GO)
all_mf <- all_go %>% filter(ONTOLOGY == "MF", !is.na(GO)) %>% distinct(SYMBOL, GO)

message(sprintf("  BP annotations: %d gene-term pairs (%d genes)",
                nrow(all_bp), n_distinct(all_bp$SYMBOL)))
message(sprintf("  CC annotations: %d gene-term pairs (%d genes)",
                nrow(all_cc), n_distinct(all_cc$SYMBOL)))
message(sprintf("  MF annotations: %d gene-term pairs (%d genes)",
                nrow(all_mf), n_distinct(all_mf$SYMBOL)))

# === 4. APPROACH A: BP-ONLY (CURRENT SYSTEM) ==================================

message("\n--- APPROACH A: BP-only (current) ---")

# Map GO terms to slim ancestors via GOBPANCESTOR
ancestors_bp <- as.list(GOBPANCESTOR)
all_go_ids <- unique(all_bp$GO)

go_to_slim <- setNames(
  lapply(all_go_ids, function(go_id) {
    hits <- character(0)
    if (go_id %in% bp_slim) hits <- go_id
    anc <- ancestors_bp[[go_id]]
    if (!is.null(anc)) hits <- c(hits, intersect(anc, bp_slim))
    unique(hits)
  }),
  all_go_ids
)

# Map ALL genes to slim -> super (background)
all_gene_slim_bp <- all_bp %>%
  mutate(slim_list = go_to_slim[GO]) %>%
  unnest(slim_list) %>%
  select(SYMBOL, slim = slim_list) %>%
  distinct()

all_gene_super_bp <- all_gene_slim_bp %>%
  mutate(super = SLIM_SUPER[slim]) %>%
  filter(!is.na(super)) %>%
  distinct(SYMBOL, super)

# Foreground: specificity-weighted 1:1 assignment (same as assign_go_slim_super)
fg_gene_slim <- all_gene_slim_bp %>% filter(SYMBOL %in% fg_genes)
fg_gene_super <- fg_gene_slim %>%
  mutate(super = SLIM_SUPER[slim]) %>%
  filter(!is.na(super))

fg_term_counts <- fg_gene_slim %>% count(slim, name = "n_fg")

best_super_A <- fg_gene_super %>%
  left_join(fg_term_counts, by = "slim") %>%
  mutate(priority = ifelse(super == "Development", 2, 1)) %>%
  arrange(priority, n_fg) %>%
  group_by(SYMBOL) %>%
  slice_head(n = 1) %>%
  ungroup()

unmapped_A <- setdiff(fg_genes, best_super_A$SYMBOL)
if (length(unmapped_A) > 0) {
  best_super_A <- bind_rows(best_super_A,
    tibble(SYMBOL = unmapped_A, slim = "OTHER", super = "Other",
           n_fg = NA_integer_, priority = 3L))
}

fg_assign_A <- best_super_A %>% select(gene = SYMBOL, super)

coverage_A <- mean(fg_assign_A$super != "Other") * 100
message(sprintf("  Coverage: %.1f%% (%d/%d classified)",
                coverage_A, sum(fg_assign_A$super != "Other"), nrow(fg_assign_A)))

# === 5. APPROACH B: BP-PRIMARY, CC/MF RESCUE =================================

message("\n--- APPROACH B: BP-primary + CC/MF rescue ---")

# CC/MF rescue mapping
cc_mf_rescue <- tribble(
  ~GO,           ~ontology, ~super,
  "GO:0005739",  "CC",      "Mitochondria & Energy",
  "GO:0005576",  "CC",      "ECM & Adhesion",
  "GO:0005856",  "CC",      "Cytoskeleton & Motility",
  "GO:0005634",  "CC",      "Gene Expression",
  "GO:0005886",  "CC",      "Transport",
  "GO:0003824",  "MF",      "Metabolism",
  "GO:0005198",  "MF",      "Cytoskeleton & Motility"
)

# Get CC/MF ancestors for rescue terms
ancestors_cc <- as.list(GOCCANCESTOR)
ancestors_mf <- as.list(GOMFANCESTOR)

# For each gene, check if any of its CC/MF annotations are (or descend from) rescue terms
rescue_cc_terms <- cc_mf_rescue %>% filter(ontology == "CC") %>% pull(GO)
rescue_mf_terms <- cc_mf_rescue %>% filter(ontology == "MF") %>% pull(GO)

# Map CC annotations to rescue terms (term itself or has rescue term as ancestor)
cc_gene_rescue <- all_cc %>%
  filter(SYMBOL %in% unmapped_A) %>%
  rowwise() %>%
  mutate(rescue_hit = {
    hits <- character(0)
    if (GO %in% rescue_cc_terms) hits <- GO
    anc <- ancestors_cc[[GO]]
    if (!is.null(anc)) hits <- c(hits, intersect(anc, rescue_cc_terms))
    list(unique(hits))
  }) %>%
  ungroup() %>%
  unnest(rescue_hit) %>%
  left_join(cc_mf_rescue %>% select(GO, super), by = c("rescue_hit" = "GO")) %>%
  select(SYMBOL, super) %>%
  distinct()

mf_gene_rescue <- all_mf %>%
  filter(SYMBOL %in% unmapped_A) %>%
  rowwise() %>%
  mutate(rescue_hit = {
    hits <- character(0)
    if (GO %in% rescue_mf_terms) hits <- GO
    anc <- ancestors_mf[[GO]]
    if (!is.null(anc)) hits <- c(hits, intersect(anc, rescue_mf_terms))
    list(unique(hits))
  }) %>%
  ungroup() %>%
  unnest(rescue_hit) %>%
  left_join(cc_mf_rescue %>% select(GO, super), by = c("rescue_hit" = "GO")) %>%
  select(SYMBOL, super) %>%
  distinct()

# Combine rescue assignments (take first hit if multi-category)
rescue_all <- bind_rows(cc_gene_rescue, mf_gene_rescue) %>%
  group_by(SYMBOL) %>%
  slice_head(n = 1) %>%
  ungroup()

message(sprintf("  Rescue recovered %d of %d 'Other' genes",
                nrow(rescue_all), length(unmapped_A)))

# Merge rescue into Approach A assignments
fg_assign_B <- fg_assign_A %>%
  left_join(rescue_all %>% rename(gene = SYMBOL, rescue_super = super), by = "gene") %>%
  mutate(super = ifelse(super == "Other" & !is.na(rescue_super), rescue_super, super)) %>%
  select(gene, super)

coverage_B <- mean(fg_assign_B$super != "Other") * 100
message(sprintf("  Coverage: %.1f%% (%d/%d classified)",
                coverage_B, sum(fg_assign_B$super != "Other"), nrow(fg_assign_B)))

# Background: also rescue universe genes the same way
bg_unmapped <- setdiff(universe, all_gene_super_bp$SYMBOL)

bg_cc_rescue <- all_cc %>%
  filter(SYMBOL %in% bg_unmapped) %>%
  rowwise() %>%
  mutate(rescue_hit = {
    hits <- character(0)
    if (GO %in% rescue_cc_terms) hits <- GO
    anc <- ancestors_cc[[GO]]
    if (!is.null(anc)) hits <- c(hits, intersect(anc, rescue_cc_terms))
    list(unique(hits))
  }) %>%
  ungroup() %>%
  unnest(rescue_hit) %>%
  left_join(cc_mf_rescue %>% select(GO, super), by = c("rescue_hit" = "GO")) %>%
  select(SYMBOL, super) %>%
  distinct()

bg_mf_rescue <- all_mf %>%
  filter(SYMBOL %in% bg_unmapped) %>%
  rowwise() %>%
  mutate(rescue_hit = {
    hits <- character(0)
    if (GO %in% rescue_mf_terms) hits <- GO
    anc <- ancestors_mf[[GO]]
    if (!is.null(anc)) hits <- c(hits, intersect(anc, rescue_mf_terms))
    list(unique(hits))
  }) %>%
  ungroup() %>%
  unnest(rescue_hit) %>%
  left_join(cc_mf_rescue %>% select(GO, super), by = c("rescue_hit" = "GO")) %>%
  select(SYMBOL, super) %>%
  distinct()

bg_rescue_all <- bind_rows(bg_cc_rescue, bg_mf_rescue) %>%
  group_by(SYMBOL) %>%
  slice_head(n = 1) %>%
  ungroup()

# Expanded background for Approach B
all_gene_super_B_expanded <- bind_rows(
  all_gene_super_bp,
  bg_rescue_all %>% rename(SYMBOL = SYMBOL)
) %>% distinct(SYMBOL, super)

# === 6. APPROACH C: MULTI-ONTOLOGY MERGED CATEGORIES =========================

message("\n--- APPROACH C: Multi-ontology merged ---")

# Define each category using BP + CC + MF terms
# We keep the same BP Slim terms AND add CC/MF terms
multi_ont_map <- bind_rows(
  # Keep all existing BP -> super mappings
  tibble(GO = names(SLIM_SUPER), super = SLIM_SUPER, ontology = "BP"),
  # Add CC mappings
  tribble(
    ~GO,           ~super,                     ~ontology,
    "GO:0005739",  "Mitochondria & Energy",    "CC",
    "GO:0005576",  "ECM & Adhesion",           "CC",
    "GO:0005856",  "Cytoskeleton & Motility",  "CC",
    "GO:0005634",  "Gene Expression",          "CC",
    "GO:0005886",  "Transport",                "CC",
    "GO:0005829",  "Metabolism",               "CC",      # cytosol
    "GO:0005730",  "Gene Expression",          "CC",      # nucleolus
    "GO:0015629",  "Cytoskeleton & Motility",  "CC",      # actin cytoskeleton
    "GO:0005764",  "Protein Homeostasis",      "CC",      # lysosome
    "GO:0005783",  "Transport",                "CC",      # ER
    "GO:0005794",  "Transport",                "CC",      # Golgi
    "GO:0030018",  "Muscle & Contractile",     "CC",      # Z disc
    "GO:0016459",  "Muscle & Contractile",     "CC",      # myosin complex
    "GO:0005615",  "ECM & Adhesion",           "CC",      # extracellular space
    "GO:0031012",  "ECM & Adhesion",           "CC",      # extracellular matrix
  ),
  # Add MF mappings
  tribble(
    ~GO,           ~super,                     ~ontology,
    "GO:0003824",  "Metabolism",               "MF",      # catalytic activity
    "GO:0005198",  "Cytoskeleton & Motility",  "MF",      # structural molecule
    "GO:0003700",  "Gene Expression",          "MF",      # DNA-binding TF
    "GO:0003723",  "Gene Expression",          "MF",      # RNA binding
    "GO:0005215",  "Transport",                "MF",      # transporter activity
    "GO:0003779",  "Cytoskeleton & Motility",  "MF",      # actin binding
    "GO:0008233",  "Protein Homeostasis",      "MF",      # peptidase activity
    "GO:0016887",  "Mitochondria & Energy",    "MF",      # ATPase activity
    "GO:0003924",  "Transport",                "MF",      # GTPase activity
    "GO:0042802",  "Protein Homeostasis",      "MF",      # identical protein binding
  )
)

# For Approach C, map genes via all three ontologies using ancestor traversal
# BP: use existing go_to_slim
# CC: traverse GOCCANCESTOR to find CC terms in multi_ont_map
# MF: traverse GOMFANCESTOR to find MF terms in multi_ont_map

cc_target_terms <- multi_ont_map %>% filter(ontology == "CC") %>% pull(GO)
mf_target_terms <- multi_ont_map %>% filter(ontology == "MF") %>% pull(GO)

message("  Mapping CC annotations to multi-ontology categories...")
cc_unique_ids <- unique(all_cc$GO)
go_to_cc_target <- setNames(
  lapply(cc_unique_ids, function(go_id) {
    hits <- character(0)
    if (go_id %in% cc_target_terms) hits <- go_id
    anc <- ancestors_cc[[go_id]]
    if (!is.null(anc)) hits <- c(hits, intersect(anc, cc_target_terms))
    unique(hits)
  }),
  cc_unique_ids
)

message("  Mapping MF annotations to multi-ontology categories...")
mf_unique_ids <- unique(all_mf$GO)
go_to_mf_target <- setNames(
  lapply(mf_unique_ids, function(go_id) {
    hits <- character(0)
    if (go_id %in% mf_target_terms) hits <- go_id
    anc <- ancestors_mf[[go_id]]
    if (!is.null(anc)) hits <- c(hits, intersect(anc, mf_target_terms))
    unique(hits)
  }),
  mf_unique_ids
)

# Map all genes to super-categories via CC
all_gene_super_cc <- all_cc %>%
  mutate(target_list = go_to_cc_target[GO]) %>%
  unnest(target_list) %>%
  left_join(multi_ont_map %>% filter(ontology == "CC") %>% select(GO, super),
            by = c("target_list" = "GO")) %>%
  filter(!is.na(super)) %>%
  select(SYMBOL, super) %>%
  distinct()

# Map all genes to super-categories via MF
all_gene_super_mf <- all_mf %>%
  mutate(target_list = go_to_mf_target[GO]) %>%
  unnest(target_list) %>%
  left_join(multi_ont_map %>% filter(ontology == "MF") %>% select(GO, super),
            by = c("target_list" = "GO")) %>%
  filter(!is.na(super)) %>%
  select(SYMBOL, super) %>%
  distinct()

# Combine all three ontologies for background
all_gene_super_C <- bind_rows(
  all_gene_super_bp,
  all_gene_super_cc,
  all_gene_super_mf
) %>% distinct(SYMBOL, super)

# Foreground assignment: specificity-weighted like Approach A but from multi-ont pool
fg_super_C_pool <- all_gene_super_C %>% filter(SYMBOL %in% fg_genes)

# Count per gene per super to pick most specific
fg_super_C_counts <- fg_super_C_pool %>%
  count(super, name = "n_fg_super")

best_super_C <- fg_super_C_pool %>%
  left_join(fg_super_C_counts, by = "super") %>%
  mutate(priority = ifelse(super == "Development", 2, 1)) %>%
  arrange(priority, n_fg_super) %>%
  group_by(SYMBOL) %>%
  slice_head(n = 1) %>%
  ungroup()

unmapped_C <- setdiff(fg_genes, best_super_C$SYMBOL)
if (length(unmapped_C) > 0) {
  best_super_C <- bind_rows(best_super_C,
    tibble(SYMBOL = unmapped_C, super = "Other",
           n_fg_super = NA_integer_, priority = 3L))
}

fg_assign_C <- best_super_C %>% select(gene = SYMBOL, super)

coverage_C <- mean(fg_assign_C$super != "Other") * 100
message(sprintf("  Coverage: %.1f%% (%d/%d classified)",
                coverage_C, sum(fg_assign_C$super != "Other"), nrow(fg_assign_C)))

# === 7. FISHER'S EXACT TESTS FOR ALL APPROACHES ==============================

message("\n=== RUNNING FISHER'S EXACT TESTS ===\n")

categories <- setdiff(SUPER_CATEGORY_ORDER, "Other")
n_fg <- length(fg_genes)
n_bg <- length(universe)

run_fisher <- function(fg_assignment, bg_mapping, approach_name) {
  results <- map_dfr(categories, function(cat) {
    fg_in <- sum(fg_assignment$super == cat)
    fg_out <- n_fg - fg_in
    
    # Background: genes in this category (from the full universe mapping)
    bg_genes_in_cat <- bg_mapping %>% filter(super == cat) %>% pull(SYMBOL) %>% unique()
    bg_in <- length(bg_genes_in_cat)
    bg_out <- n_bg - bg_in
    
    # Fisher's exact test (over-representation)
    mat <- matrix(c(fg_in, fg_out, bg_in - fg_in, bg_out - fg_out), nrow = 2)
    # Ensure non-negative
    mat[mat < 0] <- 0
    ft <- fisher.test(mat, alternative = "greater")
    
    # Fold enrichment
    fg_frac <- fg_in / n_fg
    bg_frac <- bg_in / n_bg
    fold_enrich <- ifelse(bg_frac > 0, fg_frac / bg_frac, NA_real_)
    
    tibble(
      category = cat,
      approach = approach_name,
      fg_count = fg_in,
      bg_count = bg_in,
      fisher_p = ft$p.value,
      fold_enrichment = round(fold_enrich, 3),
      fg_pct = round(fg_frac * 100, 1),
      bg_pct = round(bg_frac * 100, 1)
    )
  })
  results
}

# Approach A: BP-only
results_A <- run_fisher(fg_assign_A, all_gene_super_bp, "A_BP_only")

# Approach B with BP-only background (rescue FG, keep BG as BP)
results_B_bp_bg <- run_fisher(fg_assign_B, all_gene_super_bp, "B_rescue_bp_bg")

# Approach B with expanded background (rescue both FG and BG)
results_B_exp_bg <- run_fisher(fg_assign_B, all_gene_super_B_expanded, "B_rescue_exp_bg")

# Approach C: full multi-ontology
results_C <- run_fisher(fg_assign_C, all_gene_super_C, "C_multi_ont")

# === 8. SUMMARY OUTPUT ========================================================

message("\n")
message("=" |> strrep(100))
message("COVERAGE COMPARISON")
message("=" |> strrep(100))
message(sprintf("  Approach A (BP-only):     %.1f%% (%d/%d)",
                coverage_A, sum(fg_assign_A$super != "Other"), n_fg))
message(sprintf("  Approach B (BP+rescue):   %.1f%% (%d/%d)",
                coverage_B, sum(fg_assign_B$super != "Other"), n_fg))
message(sprintf("  Approach C (multi-ont):   %.1f%% (%d/%d)",
                coverage_C, sum(fg_assign_C$super != "Other"), n_fg))

# Print approach A detail
message("\n")
message("=" |> strrep(100))
message("APPROACH A: BP-ONLY (current system)")
message("=" |> strrep(100))
cat(sprintf("%-30s  fg  bg   fg%%  bg%%   fold   p-value\n",
            "Category"))
cat(strrep("-", 85), "\n")
for (i in seq_len(nrow(results_A))) {
  r <- results_A[i, ]
  sig <- ifelse(r$fisher_p < 0.05, " *", "")
  cat(sprintf("%-30s  %3d %4d  %4.1f %4.1f  %5.2f  %.4f%s\n",
              r$category, r$fg_count, r$bg_count,
              r$fg_pct, r$bg_pct, r$fold_enrichment, r$fisher_p, sig))
}

# Print approach B detail (expanded bg)
message("\n")
message("=" |> strrep(100))
message("APPROACH B: BP-primary + CC/MF rescue (expanded background)")
message("=" |> strrep(100))
cat(sprintf("%-30s  fg  bg   fg%%  bg%%   fold   p-value\n",
            "Category"))
cat(strrep("-", 85), "\n")
for (i in seq_len(nrow(results_B_exp_bg))) {
  r <- results_B_exp_bg[i, ]
  sig <- ifelse(r$fisher_p < 0.05, " *", "")
  cat(sprintf("%-30s  %3d %4d  %4.1f %4.1f  %5.2f  %.4f%s\n",
              r$category, r$fg_count, r$bg_count,
              r$fg_pct, r$bg_pct, r$fold_enrichment, r$fisher_p, sig))
}

# Print approach C detail
message("\n")
message("=" |> strrep(100))
message("APPROACH C: MULTI-ONTOLOGY MERGED")
message("=" |> strrep(100))
cat(sprintf("%-30s  fg  bg   fg%%  bg%%   fold   p-value\n",
            "Category"))
cat(strrep("-", 85), "\n")
for (i in seq_len(nrow(results_C))) {
  r <- results_C[i, ]
  sig <- ifelse(r$fisher_p < 0.05, " *", "")
  cat(sprintf("%-30s  %3d %4d  %4.1f %4.1f  %5.2f  %.4f%s\n",
              r$category, r$fg_count, r$bg_count,
              r$fg_pct, r$bg_pct, r$fold_enrichment, r$fisher_p, sig))
}

# === 9. INFLATION/DEFLATION COMPARISON TABLE ==================================

message("\n")
message("=" |> strrep(100))
message("INFLATION/DEFLATION COMPARISON: C vs A")
message("=" |> strrep(100))

comparison <- results_A %>%
  select(category, fg_A = fg_count, bg_A = bg_count, fold_A = fold_enrichment, p_A = fisher_p) %>%
  left_join(
    results_B_exp_bg %>%
      select(category, fg_B = fg_count, bg_B = bg_count, fold_B = fold_enrichment, p_B = fisher_p),
    by = "category"
  ) %>%
  left_join(
    results_C %>%
      select(category, fg_C = fg_count, bg_C = bg_count, fold_C = fold_enrichment, p_C = fisher_p),
    by = "category"
  ) %>%
  mutate(
    fold_change_pct = round((fold_C / fold_A - 1) * 100, 1),
    direction = case_when(
      fold_change_pct > 50  ~ "INFLATED >50%",
      fold_change_pct > 20  ~ "inflated 20-50%",
      fold_change_pct < -50 ~ "DEFLATED >50%",
      fold_change_pct < -20 ~ "deflated 20-50%",
      TRUE                  ~ "stable"
    )
  )

cat(sprintf("\n%-30s | ---- Approach A ---- | ---- Approach B ---- | ---- Approach C ---- | Change\n",
            "Category"))
cat(strrep("-", 130), "\n")
cat(sprintf("%-30s | fg  bg  fold  p-val  | fg  bg  fold  p-val  | fg  bg  fold  p-val  | C vs A\n",
            ""))
cat(strrep("-", 130), "\n")

for (i in seq_len(nrow(comparison))) {
  r <- comparison[i, ]
  flag <- ifelse(r$direction != "stable",
                 sprintf(" <-- %s", r$direction), "")
  cat(sprintf("%-30s | %3d %4d %5.2f %.4f | %3d %4d %5.2f %.4f | %3d %4d %5.2f %.4f | %+.1f%%%s\n",
              r$category,
              r$fg_A, r$bg_A, r$fold_A, r$p_A,
              r$fg_B, r$bg_B, r$fold_B, r$p_B,
              r$fg_C, r$bg_C, r$fold_C, r$p_C,
              r$fold_change_pct, flag))
}

# === 10. APPROACH B: RESCUE FG WITH BP-ONLY BACKGROUND (inflation test) =======

message("\n")
message("=" |> strrep(100))
message("APPROACH B INFLATION TEST: rescue FG only, keep BP-only background")
message("=" |> strrep(100))
message("This is the biased version — if FG is rescued but BG is not, enrichment inflates.\n")

comparison_B <- results_A %>%
  select(category, fg_A = fg_count, fold_A = fold_enrichment) %>%
  left_join(
    results_B_bp_bg %>%
      select(category, fg_B_bias = fg_count, fold_B_bias = fold_enrichment, p_B_bias = fisher_p),
    by = "category"
  ) %>%
  left_join(
    results_B_exp_bg %>%
      select(category, fg_B_fair = fg_count, fold_B_fair = fold_enrichment, p_B_fair = fisher_p),
    by = "category"
  ) %>%
  mutate(
    bias_pct = round((fold_B_bias / fold_A - 1) * 100, 1),
    fair_pct = round((fold_B_fair / fold_A - 1) * 100, 1)
  )

cat(sprintf("%-30s | fold_A | B(bp_bg)  bias%%  | B(exp_bg) fair%%  | Inflated?\n", "Category"))
cat(strrep("-", 95), "\n")
for (i in seq_len(nrow(comparison_B))) {
  r <- comparison_B[i, ]
  flag <- ifelse(abs(r$bias_pct) > 20 & abs(r$fair_pct) <= 20,
                 " <-- rescue inflates without bg correction", "")
  cat(sprintf("%-30s | %5.2f  | %5.2f  %+6.1f%% | %5.2f  %+6.1f%% |%s\n",
              r$category, r$fold_A,
              r$fold_B_bias, r$bias_pct,
              r$fold_B_fair, r$fair_pct, flag))
}

# === 11. SUMMARY =============================================================

message("\n")
message("=" |> strrep(100))
message("KEY FINDINGS")
message("=" |> strrep(100))

n_inflated <- sum(abs(comparison$fold_change_pct) > 50, na.rm = TRUE)
n_moderate <- sum(abs(comparison$fold_change_pct) > 20 &
                  abs(comparison$fold_change_pct) <= 50, na.rm = TRUE)

message(sprintf("\n  Categories with >50%% fold change (C vs A): %d", n_inflated))
message(sprintf("  Categories with 20-50%% fold change (C vs A): %d", n_moderate))
message(sprintf("  Categories stable (<20%% change): %d",
                nrow(comparison) - n_inflated - n_moderate))

# Identify the most problematic categories
problem_cats <- comparison %>% filter(abs(fold_change_pct) > 50)
if (nrow(problem_cats) > 0) {
  message("\n  Problematic categories (>50% change):")
  for (i in seq_len(nrow(problem_cats))) {
    r <- problem_cats[i, ]
    message(sprintf("    %s: A=%.2f -> C=%.2f (%+.1f%%)",
                    r$category, r$fold_A, r$fold_C, r$fold_change_pct))
  }
}

# Coverage gain assessment
message(sprintf("\n  Coverage gain from multi-ontology: %.1f%% -> %.1f%% (+%.1f%%)",
                coverage_A, coverage_C, coverage_C - coverage_A))

# Genes rescued by each approach
rescue_B_genes <- fg_assign_B %>%
  filter(super != "Other") %>%
  filter(!(gene %in% (fg_assign_A %>% filter(super != "Other") %>% pull(gene)))) %>%
  pull(gene)
rescue_C_genes <- fg_assign_C %>%
  filter(super != "Other") %>%
  filter(!(gene %in% (fg_assign_A %>% filter(super != "Other") %>% pull(gene)))) %>%
  pull(gene)

message(sprintf("\n  Genes rescued by Approach B: %d", length(rescue_B_genes)))
if (length(rescue_B_genes) > 0 && length(rescue_B_genes) <= 30)
  message(sprintf("    %s", paste(rescue_B_genes, collapse = ", ")))
message(sprintf("  Genes rescued by Approach C: %d", length(rescue_C_genes)))
if (length(rescue_C_genes) > 0 && length(rescue_C_genes) <= 30)
  message(sprintf("    %s", paste(rescue_C_genes, collapse = ", ")))

# Where did rescued genes go?
if (length(rescue_C_genes) > 0) {
  message("\n  Distribution of rescued genes (Approach C):")
  rescue_dist <- fg_assign_C %>%
    filter(gene %in% rescue_C_genes) %>%
    count(super) %>%
    arrange(desc(n))
  for (i in seq_len(nrow(rescue_dist))) {
    message(sprintf("    %s: %d", rescue_dist$super[i], rescue_dist$n[i]))
  }
}

message("\n=== DONE ===\n")
