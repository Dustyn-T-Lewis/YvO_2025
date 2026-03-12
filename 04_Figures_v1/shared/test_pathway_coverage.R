################################################################################
#   Pathway / Ontology Database Coverage Test
#
#   Tests how well different annotation databases cover the protein lists
#   used in F2 Panel F (~131 proteins) and F3 Panel F (~246 proteins),
#   compared against the full proteome universe (~2124 proteins).
#
#   Databases tested:
#     1. GO:BP (all terms)
#     2. GO:CC (all terms)
#     3. GO:MF (all terms)
#     4. GO Slim BP (62 terms)
#     5. GO Slim BP -> super-categories (12 categories)
#     6. Reactome (all human pathways)
#     7. MSigDB Hallmark (50 gene sets)
#     8. KEGG Legacy (via msigdbr)
#
#   Usage: Rscript 04_Figures/shared/test_pathway_coverage.R
################################################################################

suppressPackageStartupMessages({
  library(tidyverse)
  library(org.Hs.eg.db)
  library(AnnotationDbi)
  library(GO.db)
  library(reactome.db)
  library(msigdbr)
})

# Ensure dplyr::select takes precedence over AnnotationDbi::select
select <- dplyr::select
filter <- dplyr::filter

# === 0. SETUP =================================================================

setwd("/Users/dtl0018/Desktop/A_Proteomics_Analysis/A_YvO_2025/")
set.seed(42)

# Source GO Slim definitions
source("04_Figures/shared/go_slim_categories.R")

# === 1. LOAD PROTEIN LISTS ===================================================

message("=== Loading protein data ===")
dep_df <- read_csv("03_DEP/c_data/03_combined_results.csv", show_col_types = FALSE)
message(sprintf("  Total proteins in dep_df: %d", nrow(dep_df)))

# Universe: all genes
universe <- unique(dep_df$gene)
message(sprintf("  Universe (unique genes): %d", length(universe)))

# F2 foreground: any significant training or interaction pi_score
f2_fg <- dep_df %>%
  filter(pi_score_Training_Young < 0.05 |
         pi_score_Training_Old < 0.05 |
         pi_score_Interaction < 0.05) %>%
  pull(gene) %>% unique()
message(sprintf("  F2 foreground: %d genes", length(f2_fg)))

# F3 foreground: any significant aging or reversal pi_score
f3_fg <- dep_df %>%
  filter(pi_score_Aging < 0.05 |
         pi_score_Reversal < 0.05) %>%
  pull(gene) %>% unique()
message(sprintf("  F3 foreground: %d genes", length(f3_fg)))

# === 2. HELPER: MAP GENE SYMBOLS TO ENTREZ ===================================

message("\n=== Mapping gene symbols to Entrez IDs ===")
suppressMessages({
  entrez_map <- AnnotationDbi::mapIds(org.Hs.eg.db,
                                       keys = universe,
                                       keytype = "SYMBOL",
                                       column = "ENTREZID",
                                       multiVals = "first")
})
n_mapped_entrez <- sum(!is.na(entrez_map))
message(sprintf("  %d / %d universe genes mapped to Entrez IDs (%.1f%%)",
                n_mapped_entrez, length(universe),
                100 * n_mapped_entrez / length(universe)))

# === 3. DATABASE COVERAGE FUNCTIONS ==========================================

# Generic coverage reporter
report_coverage <- function(db_name, covered_universe, covered_f2, covered_f3,
                            n_categories_f2 = NA, n_categories_f3 = NA,
                            n_toplevel = NA) {
  tibble(
    database = db_name,
    univ_covered = sum(covered_universe),
    univ_total = length(universe),
    univ_pct = round(100 * sum(covered_universe) / length(universe), 1),
    f2_covered = sum(covered_f2),
    f2_total = length(f2_fg),
    f2_pct = round(100 * sum(covered_f2) / length(f2_fg), 1),
    f3_covered = sum(covered_f3),
    f3_total = length(f3_fg),
    f3_pct = round(100 * sum(covered_f3) / length(f3_fg), 1),
    n_categories_f2 = n_categories_f2,
    n_categories_f3 = n_categories_f3,
    n_toplevel = n_toplevel
  )
}

results <- list()

# === 3a. GO:BP (all terms) ===================================================
message("\n=== Testing GO:BP (all terms) ===")
suppressMessages({
  go_all <- AnnotationDbi::select(org.Hs.eg.db,
                                   keys = universe,
                                   keytype = "SYMBOL",
                                   columns = c("SYMBOL", "GO", "ONTOLOGY"))
})
go_bp_genes <- go_all %>%
  filter(ONTOLOGY == "BP", !is.na(GO)) %>%
  distinct(SYMBOL) %>%
  pull(SYMBOL)

n_bp_terms_f2 <- go_all %>%
  filter(ONTOLOGY == "BP", !is.na(GO), SYMBOL %in% f2_fg) %>%
  distinct(GO) %>% nrow()
n_bp_terms_f3 <- go_all %>%
  filter(ONTOLOGY == "BP", !is.na(GO), SYMBOL %in% f3_fg) %>%
  distinct(GO) %>% nrow()

results[["GO:BP (all)"]] <- report_coverage(
  "GO:BP (all terms)",
  universe %in% go_bp_genes,
  f2_fg %in% go_bp_genes,
  f3_fg %in% go_bp_genes,
  n_bp_terms_f2, n_bp_terms_f3
)
message(sprintf("  Universe: %d/%d (%.1f%%), F2: %d/%d, F3: %d/%d",
                sum(universe %in% go_bp_genes), length(universe),
                100 * mean(universe %in% go_bp_genes),
                sum(f2_fg %in% go_bp_genes), length(f2_fg),
                sum(f3_fg %in% go_bp_genes), length(f3_fg)))

# === 3b. GO:CC ===============================================================
message("\n=== Testing GO:CC (all terms) ===")
go_cc_genes <- go_all %>%
  filter(ONTOLOGY == "CC", !is.na(GO)) %>%
  distinct(SYMBOL) %>%
  pull(SYMBOL)

n_cc_terms_f2 <- go_all %>%
  filter(ONTOLOGY == "CC", !is.na(GO), SYMBOL %in% f2_fg) %>%
  distinct(GO) %>% nrow()
n_cc_terms_f3 <- go_all %>%
  filter(ONTOLOGY == "CC", !is.na(GO), SYMBOL %in% f3_fg) %>%
  distinct(GO) %>% nrow()

results[["GO:CC"]] <- report_coverage(
  "GO:CC (all terms)",
  universe %in% go_cc_genes,
  f2_fg %in% go_cc_genes,
  f3_fg %in% go_cc_genes,
  n_cc_terms_f2, n_cc_terms_f3
)
message(sprintf("  Universe: %d/%d (%.1f%%), F2: %d/%d, F3: %d/%d",
                sum(universe %in% go_cc_genes), length(universe),
                100 * mean(universe %in% go_cc_genes),
                sum(f2_fg %in% go_cc_genes), length(f2_fg),
                sum(f3_fg %in% go_cc_genes), length(f3_fg)))

# === 3c. GO:MF ===============================================================
message("\n=== Testing GO:MF (all terms) ===")
go_mf_genes <- go_all %>%
  filter(ONTOLOGY == "MF", !is.na(GO)) %>%
  distinct(SYMBOL) %>%
  pull(SYMBOL)

n_mf_terms_f2 <- go_all %>%
  filter(ONTOLOGY == "MF", !is.na(GO), SYMBOL %in% f2_fg) %>%
  distinct(GO) %>% nrow()
n_mf_terms_f3 <- go_all %>%
  filter(ONTOLOGY == "MF", !is.na(GO), SYMBOL %in% f3_fg) %>%
  distinct(GO) %>% nrow()

results[["GO:MF"]] <- report_coverage(
  "GO:MF (all terms)",
  universe %in% go_mf_genes,
  f2_fg %in% go_mf_genes,
  f3_fg %in% go_mf_genes,
  n_mf_terms_f2, n_mf_terms_f3
)
message(sprintf("  Universe: %d/%d (%.1f%%), F2: %d/%d, F3: %d/%d",
                sum(universe %in% go_mf_genes), length(universe),
                100 * mean(universe %in% go_mf_genes),
                sum(f2_fg %in% go_mf_genes), length(f2_fg),
                sum(f3_fg %in% go_mf_genes), length(f3_fg)))

# === 3d. GO Slim BP (62 terms) ===============================================
message("\n=== Testing GO Slim BP (62 terms) ===")

go_bp_detailed <- go_all %>%
  filter(ONTOLOGY == "BP", !is.na(GO)) %>%
  distinct(SYMBOL, GO)

ancestors <- as.list(GOBPANCESTOR)
all_go_ids <- unique(go_bp_detailed$GO)

message("  Computing GO Slim ancestor mapping...")
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

gene_slim_map <- go_bp_detailed %>%
  mutate(slim_list = go_to_slim[GO]) %>%
  unnest(slim_list) %>%
  distinct(SYMBOL, slim = slim_list)

slim_covered_genes <- unique(gene_slim_map$SYMBOL)

n_slim_terms_f2 <- gene_slim_map %>%
  filter(SYMBOL %in% f2_fg) %>%
  distinct(slim) %>% nrow()
n_slim_terms_f3 <- gene_slim_map %>%
  filter(SYMBOL %in% f3_fg) %>%
  distinct(slim) %>% nrow()

results[["GO Slim BP"]] <- report_coverage(
  "GO Slim BP (62 terms)",
  universe %in% slim_covered_genes,
  f2_fg %in% slim_covered_genes,
  f3_fg %in% slim_covered_genes,
  n_slim_terms_f2, n_slim_terms_f3,
  n_toplevel = length(bp_slim)
)
message(sprintf("  Universe: %d/%d (%.1f%%), F2: %d/%d, F3: %d/%d",
                sum(universe %in% slim_covered_genes), length(universe),
                100 * mean(universe %in% slim_covered_genes),
                sum(f2_fg %in% slim_covered_genes), length(f2_fg),
                sum(f3_fg %in% slim_covered_genes), length(f3_fg)))

# === 3e. GO Slim BP -> Super-categories (12 categories) ======================
message("\n=== Testing GO Slim BP -> Super-categories (12 groups) ===")

gene_super_map <- gene_slim_map %>%
  mutate(super = SLIM_SUPER[slim]) %>%
  filter(!is.na(super)) %>%
  distinct(SYMBOL, super)

super_covered_genes <- unique(gene_super_map$SYMBOL)

n_super_cats_f2 <- gene_super_map %>%
  filter(SYMBOL %in% f2_fg) %>%
  distinct(super) %>% nrow()
n_super_cats_f3 <- gene_super_map %>%
  filter(SYMBOL %in% f3_fg) %>%
  distinct(super) %>% nrow()

results[["GO Slim Super"]] <- report_coverage(
  "GO Slim -> Super-cats (12)",
  universe %in% super_covered_genes,
  f2_fg %in% super_covered_genes,
  f3_fg %in% super_covered_genes,
  n_super_cats_f2, n_super_cats_f3,
  n_toplevel = length(unique(SLIM_SUPER))
)
message(sprintf("  Universe: %d/%d (%.1f%%), F2: %d/%d, F3: %d/%d",
                sum(universe %in% super_covered_genes), length(universe),
                100 * mean(universe %in% super_covered_genes),
                sum(f2_fg %in% super_covered_genes), length(f2_fg),
                sum(f3_fg %in% super_covered_genes), length(f3_fg)))

f2_uncovered_super <- setdiff(f2_fg, super_covered_genes)
f3_uncovered_super <- setdiff(f3_fg, super_covered_genes)
message(sprintf("  F2 uncovered by super-cats: %d genes: %s",
                length(f2_uncovered_super),
                paste(sort(f2_uncovered_super), collapse = ", ")))
message(sprintf("  F3 uncovered by super-cats: %d genes: %s",
                length(f3_uncovered_super),
                paste(sort(f3_uncovered_super), collapse = ", ")))

# === 3f. Reactome =============================================================
message("\n=== Testing Reactome ===")

entrez_valid <- entrez_map[!is.na(entrez_map)]

# Get Entrez -> Reactome pathway mapping
suppressMessages({
  ext2path <- as.list(reactomeEXTID2PATHID)
})

# Filter to human pathways (R-HSA-)
reactome_covered_entrez <- names(ext2path)[sapply(ext2path, function(pids) {
  any(grepl("^R-HSA-", pids))
})]

reactome_covered_genes <- names(entrez_valid[entrez_valid %in% reactome_covered_entrez])

# Count unique pathways hit
f2_entrez_in_reactome <- entrez_valid[entrez_valid %in% reactome_covered_entrez &
                                      names(entrez_valid) %in% f2_fg]
reactome_paths_f2 <- unique(unlist(ext2path[f2_entrez_in_reactome]))
reactome_paths_f2 <- reactome_paths_f2[grepl("^R-HSA-", reactome_paths_f2)]

f3_entrez_in_reactome <- entrez_valid[entrez_valid %in% reactome_covered_entrez &
                                      names(entrez_valid) %in% f3_fg]
reactome_paths_f3 <- unique(unlist(ext2path[f3_entrez_in_reactome]))
reactome_paths_f3 <- reactome_paths_f3[grepl("^R-HSA-", reactome_paths_f3)]

# Top-level Reactome categories (well-known root pathways for Homo sapiens)
# These are the ~27 root-level categories that all other pathways nest under
REACTOME_TOPLEVEL <- c(
  "R-HSA-1474244" = "Extracellular matrix organization",
  "R-HSA-1500931" = "Cell-Cell communication",
  "R-HSA-162582"  = "Signal Transduction",
  "R-HSA-1643685" = "Disease",
  "R-HSA-168256"  = "Immune System",
  "R-HSA-1852241" = "Organelle biogenesis and maintenance",
  "R-HSA-392499"  = "Metabolism of proteins",
  "R-HSA-397014"  = "Muscle contraction",
  "R-HSA-5357801" = "Programmed Cell Death",
  "R-HSA-74160"   = "Gene expression (Transcription)",
  "R-HSA-8953897" = "Cellular responses to stimuli",
  "R-HSA-9006931" = "Signaling by Nuclear Receptors",
  "R-HSA-9609507" = "Protein localization",
  "R-HSA-1266738" = "Developmental Biology",
  "R-HSA-1430728" = "Metabolism",
  "R-HSA-5653656" = "Vesicle-mediated transport",
  "R-HSA-382551"  = "Transport of small molecules",
  "R-HSA-1640170" = "Cell Cycle",
  "R-HSA-9612973" = "Autophagy",
  "R-HSA-73886"   = "Chromosome Maintenance",
  "R-HSA-109582"  = "Hemostasis",
  "R-HSA-112316"  = "Neuronal System",
  "R-HSA-73894"   = "DNA Repair",
  "R-HSA-69306"   = "DNA Replication",
  "R-HSA-8963743" = "Digestion and absorption",
  "R-HSA-1474165" = "Reproduction"
)

# For each gene, find which Reactome pathways it's in, then check if the
# pathway name suggests a top-level category. Since we lack hierarchy,
# we check direct membership only at the top level.
suppressMessages({
  path2name_df <- AnnotationDbi::toTable(reactomePATHID2NAME)
})

# Check how many top-level pathways have genes from each foreground
toplevel_hit_f2 <- names(REACTOME_TOPLEVEL)[names(REACTOME_TOPLEVEL) %in% reactome_paths_f2]
toplevel_hit_f3 <- names(REACTOME_TOPLEVEL)[names(REACTOME_TOPLEVEL) %in% reactome_paths_f3]

# But more useful: use reactomePATHID2EXTID to check membership in top-level pathways
suppressMessages({
  path2ext <- as.list(reactomePATHID2EXTID)
})

# Count genes mapping to each top-level
toplevel_gene_counts <- tibble(
  pathway_id = names(REACTOME_TOPLEVEL),
  pathway_name = unname(REACTOME_TOPLEVEL)
) %>%
  rowwise() %>%
  mutate(
    members = list(path2ext[[pathway_id]]),
    n_universe = sum(entrez_valid %in% members),
    n_f2 = sum(entrez_valid[names(entrez_valid) %in% f2_fg] %in% members),
    n_f3 = sum(entrez_valid[names(entrez_valid) %in% f3_fg] %in% members)
  ) %>%
  ungroup() %>%
  dplyr::select(-members)

results[["Reactome"]] <- report_coverage(
  "Reactome (all human)",
  universe %in% reactome_covered_genes,
  f2_fg %in% reactome_covered_genes,
  f3_fg %in% reactome_covered_genes,
  length(reactome_paths_f2), length(reactome_paths_f3),
  n_toplevel = length(REACTOME_TOPLEVEL)
)
message(sprintf("  Universe: %d/%d (%.1f%%), F2: %d/%d, F3: %d/%d",
                sum(universe %in% reactome_covered_genes), length(universe),
                100 * mean(universe %in% reactome_covered_genes),
                sum(f2_fg %in% reactome_covered_genes), length(f2_fg),
                sum(f3_fg %in% reactome_covered_genes), length(f3_fg)))
message(sprintf("  Reactome pathways hit -- F2: %d, F3: %d",
                length(reactome_paths_f2), length(reactome_paths_f3)))

# === 3g. MSigDB Hallmark =====================================================
message("\n=== Testing MSigDB Hallmark (50 gene sets) ===")
hallmark <- msigdbr(species = "Homo sapiens", collection = "H")
hallmark_genes <- unique(hallmark$gene_symbol)

n_hallmark_sets_f2 <- hallmark %>%
  filter(gene_symbol %in% f2_fg) %>%
  distinct(gs_name) %>% nrow()
n_hallmark_sets_f3 <- hallmark %>%
  filter(gene_symbol %in% f3_fg) %>%
  distinct(gs_name) %>% nrow()

results[["Hallmark"]] <- report_coverage(
  "MSigDB Hallmark (50)",
  universe %in% hallmark_genes,
  f2_fg %in% hallmark_genes,
  f3_fg %in% hallmark_genes,
  n_hallmark_sets_f2, n_hallmark_sets_f3,
  n_toplevel = 50
)
message(sprintf("  Universe: %d/%d (%.1f%%), F2: %d/%d, F3: %d/%d",
                sum(universe %in% hallmark_genes), length(universe),
                100 * mean(universe %in% hallmark_genes),
                sum(f2_fg %in% hallmark_genes), length(f2_fg),
                sum(f3_fg %in% hallmark_genes), length(f3_fg)))
message(sprintf("  Hallmark sets hit -- F2: %d/50, F3: %d/50",
                n_hallmark_sets_f2, n_hallmark_sets_f3))

# === 3h. KEGG Legacy (via msigdbr) ===========================================
message("\n=== Testing KEGG Legacy (via msigdbr C2:CP:KEGG_LEGACY) ===")
kegg <- msigdbr(species = "Homo sapiens", collection = "C2",
                subcollection = "CP:KEGG_LEGACY")
kegg_genes <- unique(kegg$gene_symbol)

n_kegg_sets_f2 <- kegg %>%
  filter(gene_symbol %in% f2_fg) %>%
  distinct(gs_name) %>% nrow()
n_kegg_sets_f3 <- kegg %>%
  filter(gene_symbol %in% f3_fg) %>%
  distinct(gs_name) %>% nrow()

results[["KEGG"]] <- report_coverage(
  "KEGG Legacy (via msigdbr)",
  universe %in% kegg_genes,
  f2_fg %in% kegg_genes,
  f3_fg %in% kegg_genes,
  n_kegg_sets_f2, n_kegg_sets_f3,
  n_toplevel = length(unique(kegg$gs_name))
)
message(sprintf("  Universe: %d/%d (%.1f%%), F2: %d/%d, F3: %d/%d",
                sum(universe %in% kegg_genes), length(universe),
                100 * mean(universe %in% kegg_genes),
                sum(f2_fg %in% kegg_genes), length(f2_fg),
                sum(f3_fg %in% kegg_genes), length(f3_fg)))
message(sprintf("  KEGG pathways hit -- F2: %d, F3: %d (of %d total)",
                n_kegg_sets_f2, n_kegg_sets_f3, length(unique(kegg$gs_name))))

# === 3i. Reactome via msigdbr (for comparison) ===============================
message("\n=== Testing Reactome via msigdbr (C2:CP:REACTOME) ===")
reactome_msigdb <- msigdbr(species = "Homo sapiens", collection = "C2",
                           subcollection = "CP:REACTOME")
reactome_msigdb_genes <- unique(reactome_msigdb$gene_symbol)

n_reactome_msigdb_f2 <- reactome_msigdb %>%
  filter(gene_symbol %in% f2_fg) %>%
  distinct(gs_name) %>% nrow()
n_reactome_msigdb_f3 <- reactome_msigdb %>%
  filter(gene_symbol %in% f3_fg) %>%
  distinct(gs_name) %>% nrow()

results[["Reactome (msigdbr)"]] <- report_coverage(
  "Reactome (via msigdbr)",
  universe %in% reactome_msigdb_genes,
  f2_fg %in% reactome_msigdb_genes,
  f3_fg %in% reactome_msigdb_genes,
  n_reactome_msigdb_f2, n_reactome_msigdb_f3,
  n_toplevel = length(unique(reactome_msigdb$gs_name))
)
message(sprintf("  Universe: %d/%d (%.1f%%), F2: %d/%d, F3: %d/%d",
                sum(universe %in% reactome_msigdb_genes), length(universe),
                100 * mean(universe %in% reactome_msigdb_genes),
                sum(f2_fg %in% reactome_msigdb_genes), length(f2_fg),
                sum(f3_fg %in% reactome_msigdb_genes), length(f3_fg)))

# === 4. SUMMARY TABLE ========================================================

message("\n")
message(strrep("=", 120))
message("PATHWAY / ONTOLOGY DATABASE COVERAGE SUMMARY")
message(strrep("=", 120))

summary_df <- bind_rows(results) %>%
  arrange(desc(univ_pct))

cat("\n")
cat(sprintf("%-30s | %16s | %16s | %16s | %8s | %8s | %8s\n",
            "Database", "Universe", "F2 Foreground", "F3 Foreground",
            "Cats F2", "Cats F3", "TopLvl"))
cat(strrep("-", 120), "\n")

for (i in seq_len(nrow(summary_df))) {
  row <- summary_df[i, ]
  cat(sprintf("%-30s | %4d/%4d %5.1f%% | %4d/%4d  %5.1f%% | %4d/%4d  %5.1f%% | %8s | %8s | %8s\n",
              row$database,
              row$univ_covered, row$univ_total, row$univ_pct,
              row$f2_covered, row$f2_total, row$f2_pct,
              row$f3_covered, row$f3_total, row$f3_pct,
              ifelse(is.na(row$n_categories_f2), "-",
                     as.character(row$n_categories_f2)),
              ifelse(is.na(row$n_categories_f3), "-",
                     as.character(row$n_categories_f3)),
              ifelse(is.na(row$n_toplevel), "-",
                     as.character(row$n_toplevel))))
}

# === 5. REACTOME TOP-LEVEL CATEGORY BREAKDOWN ================================

cat("\n")
message(strrep("=", 90))
message("REACTOME TOP-LEVEL PATHWAY CATEGORIES (direct membership)")
message(strrep("=", 90))

toplevel_print <- toplevel_gene_counts %>%
  arrange(desc(n_universe)) %>%
  filter(n_universe > 0 | n_f2 > 0 | n_f3 > 0)

cat(sprintf("\n%-45s | %8s | %6s | %6s\n",
            "Reactome Top-Level Pathway", "Universe", "F2", "F3"))
cat(strrep("-", 75), "\n")
for (i in seq_len(nrow(toplevel_print))) {
  row <- toplevel_print[i, ]
  cat(sprintf("%-45s | %8d | %6d | %6d\n",
              substr(row$pathway_name, 1, 45),
              row$n_universe, row$n_f2, row$n_f3))
}

# === 6. HALLMARK SET DETAIL ===================================================

cat("\n")
message(strrep("=", 80))
message("HALLMARK GENE SETS WITH HITS IN FOREGROUND")
message(strrep("=", 80))

hallmark_detail <- hallmark %>%
  group_by(gs_name) %>%
  summarise(
    total_genes = n_distinct(gene_symbol),
    f2_hits = sum(gene_symbol %in% f2_fg),
    f3_hits = sum(gene_symbol %in% f3_fg),
    univ_hits = sum(gene_symbol %in% universe),
    .groups = "drop"
  ) %>%
  mutate(gs_short = sub("^HALLMARK_", "", gs_name)) %>%
  filter(f2_hits > 0 | f3_hits > 0) %>%
  arrange(desc(f2_hits + f3_hits))

cat(sprintf("\n%-45s | %5s | %5s | %5s | %5s\n",
            "Hallmark Set", "Total", "Univ", "F2", "F3"))
cat(strrep("-", 75), "\n")
for (i in seq_len(nrow(hallmark_detail))) {
  row <- hallmark_detail[i, ]
  cat(sprintf("%-45s | %5d | %5d | %5d | %5d\n",
              substr(row$gs_short, 1, 45),
              row$total_genes, row$univ_hits, row$f2_hits, row$f3_hits))
}

# Re-establish dplyr::select after AnnotationDbi loaded it
select <- dplyr::select
# === 7. SUPER-CATEGORY DISTRIBUTION ==========================================

cat("\n")
message(strrep("=", 80))
message("GO SLIM SUPER-CATEGORY DISTRIBUTION IN FOREGROUNDS")
message(strrep("=", 80))

message("  Computing F2 super-category assignment...")
f2_super <- assign_go_slim_super(f2_fg, universe)
message("  Computing F3 super-category assignment...")
f3_super <- assign_go_slim_super(f3_fg, universe)

f2_super_counts <- f2_super %>% count(super, name = "f2_n") %>%
  mutate(f2_pct = round(100 * f2_n / sum(f2_n), 1))
f3_super_counts <- f3_super %>% count(super, name = "f3_n") %>%
  mutate(f3_pct = round(100 * f3_n / sum(f3_n), 1))

super_combined <- full_join(f2_super_counts, f3_super_counts, by = "super") %>%
  mutate(across(where(is.numeric), ~replace_na(., 0))) %>%
  arrange(match(as.character(super), SUPER_CATEGORY_ORDER))

cat(sprintf("\n%-30s | %6s %6s | %6s %6s\n",
            "Super-category", "F2 n", "F2 %", "F3 n", "F3 %"))
cat(strrep("-", 70), "\n")
for (i in seq_len(nrow(super_combined))) {
  row <- super_combined[i, ]
  cat(sprintf("%-30s | %6d %5.1f%% | %6d %5.1f%%\n",
              as.character(row$super),
              row$f2_n, row$f2_pct,
              row$f3_n, row$f3_pct))
}
cat(sprintf("%-30s | %6d %5.1f%% | %6d %5.1f%%\n",
            "TOTAL",
            sum(super_combined$f2_n), 100.0,
            sum(super_combined$f3_n), 100.0))

# === 8. UNCOVERED GENES ======================================================

cat("\n")
message(strrep("=", 80))
message("GENES MISSING FROM EACH DATABASE")
message(strrep("=", 80))

cat("\nF2 foreground genes NOT covered:\n")
cat(sprintf("  GO:BP all:      %d genes\n", sum(!f2_fg %in% go_bp_genes)))
cat(sprintf("  GO Slim BP:     %d genes\n", sum(!f2_fg %in% slim_covered_genes)))
cat(sprintf("  GO Slim Super:  %d genes: %s\n",
            length(f2_uncovered_super),
            paste(sort(f2_uncovered_super), collapse = ", ")))
cat(sprintf("  Reactome:       %d genes: %s\n",
            sum(!f2_fg %in% reactome_covered_genes),
            paste(sort(f2_fg[!f2_fg %in% reactome_covered_genes]), collapse = ", ")))
cat(sprintf("  Hallmark:       %d genes: %s\n",
            sum(!f2_fg %in% hallmark_genes),
            paste(sort(f2_fg[!f2_fg %in% hallmark_genes]), collapse = ", ")))
cat(sprintf("  KEGG:           %d genes: %s\n",
            sum(!f2_fg %in% kegg_genes),
            paste(sort(f2_fg[!f2_fg %in% kegg_genes]), collapse = ", ")))

cat("\nF3 foreground genes NOT covered:\n")
cat(sprintf("  GO:BP all:      %d genes\n", sum(!f3_fg %in% go_bp_genes)))
cat(sprintf("  GO Slim BP:     %d genes\n", sum(!f3_fg %in% slim_covered_genes)))
cat(sprintf("  GO Slim Super:  %d genes: %s\n",
            length(f3_uncovered_super),
            paste(sort(f3_uncovered_super), collapse = ", ")))
cat(sprintf("  Reactome:       %d genes: %s\n",
            sum(!f3_fg %in% reactome_covered_genes),
            paste(sort(f3_fg[!f3_fg %in% reactome_covered_genes]), collapse = ", ")))
cat(sprintf("  Hallmark:       %d genes: %s\n",
            sum(!f3_fg %in% hallmark_genes),
            paste(sort(f3_fg[!f3_fg %in% hallmark_genes]), collapse = ", ")))
cat(sprintf("  KEGG:           %d genes: %s\n",
            sum(!f3_fg %in% kegg_genes),
            paste(sort(f3_fg[!f3_fg %in% kegg_genes]), collapse = ", ")))

# === 9. CROSS-DATABASE OVERLAP ================================================

cat("\n")
message(strrep("=", 80))
message("CROSS-DATABASE OVERLAP: GENES COVERED BY N DATABASES")
message(strrep("=", 80))

all_dbs <- list(
  "GO:BP" = go_bp_genes,
  "GO:CC" = go_cc_genes,
  "GO:MF" = go_mf_genes,
  "GO Slim BP" = slim_covered_genes,
  "GO Slim Super" = super_covered_genes,
  "Reactome" = reactome_covered_genes,
  "Hallmark" = hallmark_genes,
  "KEGG" = kegg_genes
)

for (fg_label in c("F2", "F3")) {
  fg_genes <- if (fg_label == "F2") f2_fg else f3_fg
  db_coverage <- sapply(fg_genes, function(g) {
    sum(sapply(all_dbs, function(db) g %in% db))
  })

  cat(sprintf("\n%s foreground -- number of databases covering each gene:\n", fg_label))
  tbl <- table(db_coverage)
  for (n in sort(as.integer(names(tbl)))) {
    cat(sprintf("  %d databases: %d genes (%.1f%%)\n",
                n, tbl[as.character(n)],
                100 * tbl[as.character(n)] / length(fg_genes)))
  }
  poorly_covered <- names(db_coverage[db_coverage <= 2])
  if (length(poorly_covered) > 0 && length(poorly_covered) <= 30) {
    cat(sprintf("  Genes covered by <= 2 databases: %s\n",
                paste(sort(poorly_covered), collapse = ", ")))
  }
}

message("\n=== Coverage analysis complete ===")
