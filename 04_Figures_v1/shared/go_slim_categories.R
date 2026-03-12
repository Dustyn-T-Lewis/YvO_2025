################################################################################
#   GO Slim Consolidated Pathway Assignment Utility
#   Shared by F2/panel_F.R and F3/panel_F.R
#
#   Maps genes to 12 biologically coherent consolidated pathways via GO Slim
#   Generic BP terms (62 curated) and GOBPANCESTOR hierarchy traversal.
#
#   Exports:
#     bp_slim                    — character vector of 62 GO Slim Generic BP IDs
#     SLIM_CONSOLIDATED          — named vector: GO Slim ID -> consolidated pathway name
#     CONSOLIDATED_COLORS        — named vector: consolidated pathway -> hex color
#     CONSOLIDATED_PATHWAY_ORDER — ordered factor levels for consistent plotting
#     assign_go_slim_consolidated(fg_genes, all_genes) — returns tibble(gene, slim, consolidated)
#
#   Legacy aliases (backwards compatibility):
#     SLIM_SUPER, SUPER_COLORS, SUPER_CATEGORY_ORDER, assign_go_slim_super
################################################################################

suppressPackageStartupMessages({
  library(GO.db)
  library(org.Hs.eg.db)
  library(AnnotationDbi)
  library(dplyr)
  library(tidyr)
})

# === GO Slim Generic BP terms (hard-coded from goslim_generic.obo) ============

bp_slim <- c(
  "GO:0000278", "GO:0000910", "GO:0002181", "GO:0002376", "GO:0003012",
  "GO:0003013", "GO:0003014", "GO:0003016", "GO:0005975", "GO:0006091",
  "GO:0006260", "GO:0006281", "GO:0006310", "GO:0006325", "GO:0006351",
  "GO:0006355", "GO:0006399", "GO:0006457", "GO:0006520", "GO:0006629",
  "GO:0006766", "GO:0006886", "GO:0006913", "GO:0006914", "GO:0006954",
  "GO:0007005", "GO:0007010", "GO:0007018", "GO:0007031", "GO:0007059",
  "GO:0007126", "GO:0007155", "GO:0007163", "GO:0007586", "GO:0009100",
  "GO:0012501", "GO:0016071", "GO:0016192", "GO:0023052", "GO:0030154",
  "GO:0030163", "GO:0030198", "GO:0032200", "GO:0034330", "GO:0042060",
  "GO:0042180", "GO:0042254", "GO:0044782", "GO:0048856", "GO:0048870",
  "GO:0050877", "GO:0051604", "GO:0055085", "GO:0055086", "GO:0061024",
  "GO:0065003", "GO:0071941", "GO:0072659", "GO:0098542", "GO:0098754",
  "GO:0140014", "GO:1901135"
)

# === 12 consolidated pathways: each slim GO ID -> biologically coherent group ==
# Note: GO:0023052 ("signaling") and GO:0050877 ("nervous system process")
# remain in bp_slim but are NOT mapped to any consolidated pathway. They were
# removed because GO:0023052 is too broad — nearly all signaling proteins
# match it, producing a dominant enrichment bar that squashes all others.

SLIM_CONSOLIDATED <- c(
  "GO:0003012" = "Muscle & Contractile",
  "GO:0003013" = "Circulatory System",
  "GO:0030198" = "ECM & Adhesion",
  "GO:0007155" = "ECM & Adhesion",
  "GO:0034330" = "ECM & Adhesion",
  "GO:0042060" = "ECM & Adhesion",
  "GO:0007010" = "Cytoskeleton & Motility",
  "GO:0048870" = "Cytoskeleton & Motility",
  "GO:0007018" = "Cytoskeleton & Motility",
  "GO:0007163" = "Cytoskeleton & Motility",
  "GO:0044782" = "Cytoskeleton & Motility",
  "GO:0002376" = "Immune & Inflammation",
  "GO:0006954" = "Immune & Inflammation",
  "GO:0098542" = "Immune & Inflammation",
  "GO:0006520" = "Metabolism",
  "GO:0005975" = "Metabolism",
  "GO:0006629" = "Metabolism",
  "GO:0006091" = "Metabolism",
  "GO:0042180" = "Metabolism",
  "GO:0055086" = "Metabolism",
  "GO:0006766" = "Metabolism",
  "GO:0071941" = "Metabolism",
  "GO:0098754" = "Metabolism",
  "GO:0007586" = "Metabolism",
  "GO:1901135" = "Metabolism",
  "GO:0007005" = "Mitochondria & Energy",
  "GO:0007031" = "Mitochondria & Energy",
  "GO:0006457" = "Protein Homeostasis",
  "GO:0030163" = "Protein Homeostasis",
  "GO:0006914" = "Protein Homeostasis",
  "GO:0051604" = "Protein Homeostasis",
  "GO:0065003" = "Protein Homeostasis",
  "GO:0009100" = "Protein Homeostasis",
  "GO:0055085" = "Transport",
  "GO:0016192" = "Transport",
  "GO:0006886" = "Transport",
  "GO:0006913" = "Transport",
  "GO:0072659" = "Transport",
  "GO:0061024" = "Transport",
  "GO:0006351" = "Gene Expression",
  "GO:0006355" = "Gene Expression",
  "GO:0016071" = "Gene Expression",
  "GO:0006399" = "Gene Expression",
  "GO:0002181" = "Gene Expression",
  "GO:0042254" = "Gene Expression",
  "GO:0006325" = "Gene Expression",
  "GO:0006281" = "DNA & Cell Cycle",
  "GO:0006260" = "DNA & Cell Cycle",
  "GO:0006310" = "DNA & Cell Cycle",
  "GO:0032200" = "DNA & Cell Cycle",
  "GO:0000278" = "DNA & Cell Cycle",
  "GO:0140014" = "DNA & Cell Cycle",
  "GO:0007059" = "DNA & Cell Cycle",
  "GO:0000910" = "DNA & Cell Cycle",
  "GO:0007126" = "DNA & Cell Cycle",
  "GO:0048856" = "Development",
  "GO:0030154" = "Development",
  "GO:0012501" = "Development",
  "GO:0003014" = "Development",
  "GO:0003016" = "Development"
)

# === Ordered factor levels for consistent plotting ============================

CONSOLIDATED_PATHWAY_ORDER <- c(
  "Muscle & Contractile", "Cytoskeleton & Motility", "ECM & Adhesion",
  "Metabolism", "Mitochondria & Energy", "Protein Homeostasis",
  "Transport", "Gene Expression",
  "Immune & Inflammation", "DNA & Cell Cycle", "Circulatory System",
  "Development", "Other"
)

# === Consolidated pathway colors (muted, distinguishable) =====================

CONSOLIDATED_COLORS <- c(
  "Muscle & Contractile"   = "#E57373",
  "Cytoskeleton & Motility" = "#FFB74D",
  "ECM & Adhesion"          = "#FFF176",
  "Metabolism"               = "#AED581",
  "Mitochondria & Energy"    = "#4DB6AC",
  "Protein Homeostasis"      = "#4FC3F7",
  "Transport"                = "#7986CB",
  "Gene Expression"          = "#BA68C8",
  "Immune & Inflammation"    = "#A1887F",
  "DNA & Cell Cycle"         = "#90A4AE",
  "Circulatory System"       = "#CE93D8",
  "Development"              = "#B0BEC5",
  "Other"                    = "#D0D0D0"
)

# === Main function: assign genes to GO Slim consolidated pathways =============
#
# Args:
#   fg_genes  — character vector of foreground gene symbols
#   all_genes — character vector of universe gene symbols (for background)
#   min_cat_size — merge categories with fewer than this many genes into "Other"
#
# Returns:
#   tibble with columns: gene, slim (GO ID or "OTHER"), consolidated (pathway name)

assign_go_slim_consolidated <- function(fg_genes, all_genes, min_cat_size = 2) {

  # 1. Get GO:BP annotations for all genes
  suppressMessages({
    all_entrez <- AnnotationDbi::mapIds(org.Hs.eg.db, keys = all_genes,
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

  # 2. Map GO terms to slim ancestors via GOBPANCESTOR
  ancestors  <- as.list(GOBPANCESTOR)
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

  # 3. Map ALL genes to slim terms (for background)
  all_gene_slim <- all_bp %>%
    mutate(slim_list = go_to_slim[GO]) %>%
    unnest(slim_list) %>%
    dplyr::select(SYMBOL, slim = slim_list) %>%
    distinct()

  # 4. Map to consolidated pathways
  all_gene_consolidated <- all_gene_slim %>%
    mutate(consolidated = SLIM_CONSOLIDATED[slim]) %>%
    filter(!is.na(consolidated)) %>%
    distinct(SYMBOL, consolidated)

  # 5. Specificity-weighted 1:1 assignment for foreground genes
  fg_gene_slim <- all_gene_slim %>% filter(SYMBOL %in% fg_genes)
  fg_gene_consolidated <- fg_gene_slim %>%
    mutate(consolidated = SLIM_CONSOLIDATED[slim]) %>%
    filter(!is.na(consolidated))

  fg_term_counts <- fg_gene_slim %>% count(slim, name = "n_fg")

  best_consolidated <- fg_gene_consolidated %>%
    left_join(fg_term_counts, by = "slim") %>%
    mutate(priority = ifelse(consolidated == "Development", 2, 1)) %>%
    arrange(priority, n_fg) %>%
    group_by(SYMBOL) %>%
    slice_head(n = 1) %>%
    ungroup()

  # 6. Unmapped genes -> "Other"
  unmapped <- setdiff(fg_genes, best_consolidated$SYMBOL)
  if (length(unmapped) > 0) {
    best_consolidated <- bind_rows(best_consolidated,
      tibble(SYMBOL = unmapped, slim = "OTHER", consolidated = "Other",
             n_fg = NA_integer_, priority = 3L))
  }

  # 7. Merge small pathways into "Other"
  small_cats <- best_consolidated %>% count(consolidated) %>%
    filter(n < min_cat_size, consolidated != "Other") %>% pull(consolidated)
  if (length(small_cats) > 0) {
    best_consolidated <- best_consolidated %>%
      mutate(consolidated = ifelse(consolidated %in% small_cats, "Other", consolidated))
  }

  # 8. Return clean result
  best_consolidated %>%
    transmute(gene = SYMBOL, slim, consolidated) %>%
    mutate(consolidated = factor(consolidated, levels = CONSOLIDATED_PATHWAY_ORDER))
}

# === Legacy aliases (backwards compatibility) =================================
SLIM_SUPER           <- SLIM_CONSOLIDATED
SUPER_COLORS         <- CONSOLIDATED_COLORS
SUPER_CATEGORY_ORDER <- CONSOLIDATED_PATHWAY_ORDER
assign_go_slim_super <- function(fg_genes, all_genes, min_cat_size = 2) {
  res <- assign_go_slim_consolidated(fg_genes, all_genes, min_cat_size)
  res %>% rename(super = consolidated)
}

# === Export function: document the GO Slim mapping (Option B) =================
#
# Writes supplementary tables documenting the 62 GO Slim → 12 consolidated
# pathway mapping, coverage statistics, and excluded terms.
#
# Args:
#   fg_genes  — character vector of foreground gene symbols
#   all_genes — character vector of universe gene symbols
#   outdir    — output directory for CSVs
#
# Writes:
#   {outdir}/go_slim_mapping.csv          — full 62 term → 12 category mapping
#   {outdir}/go_slim_coverage.csv         — coverage statistics
#   {outdir}/go_slim_excluded_terms.csv   — excluded terms with rationale

export_slim_mapping <- function(fg_genes, all_genes, outdir) {
  dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

  # 1. Full mapping table: GO ID → term name → consolidated category
  slim_terms <- AnnotationDbi::Term(bp_slim)
  mapping_df <- tibble(
    go_id = bp_slim,
    go_term = slim_terms[bp_slim],
    consolidated = ifelse(bp_slim %in% names(SLIM_CONSOLIDATED),
                          SLIM_CONSOLIDATED[bp_slim],
                          NA_character_)
  ) %>%
    arrange(consolidated, go_id)

  write_csv(mapping_df, file.path(outdir, "go_slim_mapping.csv"))

  # 2. Excluded terms
  excluded <- mapping_df %>%
    filter(is.na(consolidated)) %>%
    mutate(reason = case_when(
      go_id == "GO:0023052" ~ "Too broad: nearly all signaling proteins match, producing a dominant category that squashes all others",
      go_id == "GO:0050877" ~ "Nervous system process: irrelevant for skeletal muscle proteomics",
      TRUE ~ "Unknown"
    )) %>%
    dplyr::select(go_id, go_term, reason)
  write_csv(excluded, file.path(outdir, "go_slim_excluded_terms.csv"))

  # 3. Run classification and compute coverage stats
  result <- assign_go_slim_consolidated(fg_genes, all_genes)

  cat_counts <- result %>%
    count(consolidated, name = "n_genes") %>%
    arrange(desc(n_genes))

  n_total  <- nrow(result)
  n_mapped <- sum(result$consolidated != "Other")
  n_cats   <- n_distinct(result$consolidated[result$consolidated != "Other"])
  largest  <- cat_counts %>% filter(consolidated != "Other") %>% slice_head(n = 1)
  smallest <- cat_counts %>% filter(consolidated != "Other") %>% slice_tail(n = 1)

  coverage <- tibble(
    metric = c("total_proteins", "mapped_proteins", "pct_mapped",
               "n_categories", "largest_category", "largest_pct",
               "smallest_category", "smallest_n",
               "n_other", "n_slim_terms", "n_excluded_terms"),
    value = c(
      as.character(n_total),
      as.character(n_mapped),
      sprintf("%.1f", 100 * n_mapped / n_total),
      as.character(n_cats),
      largest$consolidated,
      sprintf("%.1f", 100 * largest$n_genes / n_total),
      smallest$consolidated,
      as.character(smallest$n_genes),
      as.character(n_total - n_mapped),
      as.character(length(bp_slim)),
      as.character(nrow(excluded))
    )
  )
  write_csv(coverage, file.path(outdir, "go_slim_coverage.csv"))

  # 4. Per-category counts
  write_csv(cat_counts, file.path(outdir, "go_slim_category_counts.csv"))

  message(sprintf("GO Slim export: %d/%d mapped (%.1f%%), %d categories, %d excluded terms",
                  n_mapped, n_total, 100 * n_mapped / n_total, n_cats, nrow(excluded)))

  invisible(list(
    mapping   = mapping_df,
    excluded  = excluded,
    coverage  = coverage,
    counts    = cat_counts,
    result    = result
  ))
}

# === Keyword-based functional classifier for ORA term IDs ====================
#
# Maps raw msigdbr ID strings (e.g. "HALLMARK_MYOGENESIS",
# "GOBP_SARCOMERE_ORGANIZATION") to the same 12 consolidated categories
# used in Panel F heatmaps. Priority-ordered rules using grepl().
#
# Args:
#   ids — character vector of msigdbr term IDs (gs_name strings)
#
# Returns:
#   character vector of category names (same length as ids)

classify_pathway_func <- function(ids) {
  rules <- list(
    "Muscle & Contractile"   = "MYOGEN|MYOFIBRIL|SARCOMERE|MUSCLE_|CONTRACTILE|ACTOMYOSIN|MYOSIN|I_BAND",
    "Cytoskeleton & Motility" = "CYTOSKELET|ACTIN_BIND|STRUCTURAL_MOLECULE|MOTIL|SUPRAMOLECUL",
    "ECM & Adhesion"          = "EXTRACELLULAR_MATRIX|COLLAGEN|BASEMENT_MEMBRANE|ADHESION|APICAL_JUNCTION|EMT|ENCAPSULATING",
    "Metabolism"               = "METABOL|GLYCOLY|FATTY_ACID|LIPID|AMINO_ACID|ADIPOGEN",
    "Mitochondria & Energy"    = "MITOCHOND|OXIDATIVE_PHOSPH|ELECTRON_TRANSFER|RESPIRATORY|OXIDOREDUCT",
    "Protein Homeostasis"      = "PROTEASOM|UBIQUITIN|AUTOPHAGY|MTORC1|PROTEIN_FOLD",
    "Transport"                = "TRANSPORT(?!.*ELECTRON)|VESICLE|ENDOCYT|SECRETI",
    "Gene Expression"          = "TRANSCRIPT|TRANSLAT|RIBOSOM|SPLICEOSOM|MYC_TARGET|E2F_TARGET",
    "Immune & Inflammation"    = "IMMUN|INFLAMMA|INTERFERON|IL2|IL6|TNFA|NF.KB|COMPLEMENT",
    "DNA & Cell Cycle"         = "DNA_REPAIR|CELL_CYCLE|MITOTIC|P53_PATHWAY",
    "Circulatory System"       = "ANGIOGEN|BLOOD_VESSEL|HYPOXIA",
    "Development"              = "UV_RESPONSE|GROWTH_FACTOR|WNT|HEDGEHOG|NOTCH|TGF_BETA|KRAS"
  )
  vapply(toupper(ids), function(id) {
    for (cat in names(rules)) {
      if (grepl(rules[[cat]], id, perl = TRUE)) return(cat)
    }
    "Other"
  }, character(1), USE.NAMES = FALSE)
}
