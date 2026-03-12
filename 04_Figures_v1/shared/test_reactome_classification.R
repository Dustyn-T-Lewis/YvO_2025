################################################################################
#   Option A: Reactome Top-Level Hierarchy Classification
#
#   Maps proteins to Reactome leaf pathways, then traverses the hierarchy to
#   top-level ancestors. Specificity-weighted 1:1 assignment (same heuristic
#   as GO Slim: prefer smallest foreground category).
#
#   Requires:
#     - reactome.db (Bioconductor)
#     - org.Hs.eg.db (Bioconductor)
#     - data/ReactomePathwaysRelation.txt (downloaded from reactome.org)
#     - data/ReactomePathways.txt (downloaded from reactome.org)
#
#   Exports:
#     assign_reactome_toplevel(fg_genes, all_genes, min_cat_size = 2)
#       → tibble(gene, pathway_id, category)
#
#   Output files (when run as standalone):
#     data/reactome_classification_f2.csv
#     data/reactome_classification_f3.csv
#     data/reactome_coverage_stats.csv
################################################################################

suppressPackageStartupMessages({
  library(reactome.db)
  library(org.Hs.eg.db)
  library(AnnotationDbi)
  library(dplyr)
  library(tidyr)
  library(readr)
})

# === Load Reactome hierarchy and pathway names ================================

SHARED_DIR <- dirname(if (exists("SCRIPT_DIR")) SCRIPT_DIR else "04_Figures_v1/shared/x")

rel_file <- file.path(SHARED_DIR, "data", "ReactomePathwaysRelation.txt")
pw_file  <- file.path(SHARED_DIR, "data", "ReactomePathways.txt")

if (!file.exists(rel_file)) {
  # Try alternative path
  rel_file <- "04_Figures_v1/shared/data/ReactomePathwaysRelation.txt"
  pw_file  <- "04_Figures_v1/shared/data/ReactomePathways.txt"
}

stopifnot(file.exists(rel_file), file.exists(pw_file))

# Parse hierarchy (parent -> child)
hierarchy <- read_tsv(rel_file, col_names = c("parent", "child"),
                      show_col_types = FALSE) %>%
  filter(grepl("^R-HSA-", parent), grepl("^R-HSA-", child))

# Parse pathway names
pathway_names <- read_tsv(pw_file, col_names = c("id", "name", "species"),
                          show_col_types = FALSE) %>%
  filter(species == "Homo sapiens") %>%
  select(id, name)
pw_name_map <- setNames(pathway_names$name, pathway_names$id)

# === Identify top-level Reactome categories ===================================

# Top-level = pathways that appear as parents but never as children
all_parents  <- unique(hierarchy$parent)
all_children <- unique(hierarchy$child)
top_level_ids <- setdiff(all_parents, all_children)

top_level_names <- pw_name_map[top_level_ids]
message(sprintf("Found %d top-level Reactome categories:", length(top_level_ids)))
for (i in seq_along(top_level_ids)) {
  message(sprintf("  %s: %s", top_level_ids[i], top_level_names[i]))
}

# Filter out irrelevant categories
EXCLUDE_CATEGORIES <- c(
  "Drug ADME", "Reproduction", "Sensory Perception",
  "Neuronal System", "Digestion and absorption"
)
top_level_keep <- top_level_ids[!top_level_names %in% EXCLUDE_CATEGORIES]

# === Build ancestor lookup (child → set of top-level ancestors) ===============

# Build adjacency list: child -> parents
child_to_parents <- split(hierarchy$parent, hierarchy$child)

find_top_ancestors <- function(pathway_id, memo = new.env(parent = emptyenv())) {
  if (!is.null(memo[[pathway_id]])) return(memo[[pathway_id]])
  if (pathway_id %in% top_level_ids) {
    memo[[pathway_id]] <- pathway_id
    return(pathway_id)
  }
  parents <- child_to_parents[[pathway_id]]
  if (is.null(parents) || length(parents) == 0) {
    memo[[pathway_id]] <- character(0)
    return(character(0))
  }
  ancestors <- unique(unlist(lapply(parents, find_top_ancestors, memo = memo)))
  memo[[pathway_id]] <- ancestors
  return(ancestors)
}

# Precompute all pathway -> top-level mappings
all_hsa_ids <- unique(c(hierarchy$parent, hierarchy$child))
memo_env <- new.env(parent = emptyenv())
pw_to_toplevel <- setNames(
  lapply(all_hsa_ids, find_top_ancestors, memo = memo_env),
  all_hsa_ids
)

# === Signal Transduction sub-splitting ========================================
# "Signal Transduction" is too broad (~40% of annotated proteins).
# Split it into its immediate children while keeping other top-level intact.

st_id <- top_level_ids[top_level_names == "Signal Transduction"]
if (length(st_id) == 1) {
  st_children <- hierarchy %>%
    filter(parent == st_id) %>%
    pull(child) %>%
    unique()
  message(sprintf("Splitting 'Signal Transduction' into %d sub-categories",
                  length(st_children)))

  # For pathways whose top-level is Signal Transduction, replace with
  # the immediate ST child they descend from
  st_child_to_parents <- split(hierarchy$parent, hierarchy$child)

  find_st_subcat <- function(pathway_id) {
    if (pathway_id %in% st_children) return(pathway_id)
    if (pathway_id == st_id) return(st_id)  # shouldn't happen for leaf
    parents <- child_to_parents[[pathway_id]]
    if (is.null(parents)) return(NA_character_)
    for (p in parents) {
      res <- find_st_subcat(p)
      if (!is.na(res)) return(res)
    }
    NA_character_
  }
}

# === Main classification function =============================================

assign_reactome_toplevel <- function(fg_genes, all_genes, min_cat_size = 2) {

  # 1. Map gene symbols to Entrez IDs
  suppressMessages({
    all_entrez <- AnnotationDbi::mapIds(org.Hs.eg.db, keys = all_genes,
                                        keytype = "SYMBOL", column = "ENTREZID",
                                        multiVals = "first")
  })
  entrez_to_sym <- setNames(names(all_entrez), all_entrez)

  # 2. Get Reactome pathway memberships via reactome.db
  suppressMessages({
    reactome_map <- AnnotationDbi::select(reactome.db,
                                           keys = as.character(na.omit(all_entrez)),
                                           keytype = "ENTREZID",
                                           columns = c("PATHID", "PATHNAME")) %>%
      filter(grepl("^R-HSA-", PATHID)) %>%
      mutate(gene = entrez_to_sym[ENTREZID]) %>%
      filter(!is.na(gene)) %>%
      distinct(gene, PATHID)
  })

  # 3. Map each leaf pathway to top-level category
  gene_categories <- reactome_map %>%
    mutate(top_levels = pw_to_toplevel[PATHID]) %>%
    unnest(top_levels) %>%
    filter(top_levels %in% top_level_keep) %>%
    rename(top_id = top_levels)

  # 4. Handle Signal Transduction splitting
  if (exists("st_id") && length(st_id) == 1) {
    gene_categories <- gene_categories %>%
      mutate(
        category = ifelse(
          top_id == st_id,
          {
            sub_ids <- vapply(PATHID, find_st_subcat, character(1))
            ifelse(is.na(sub_ids) | sub_ids == st_id,
                   "Signal Transduction (Other)",
                   paste0("ST: ", pw_name_map[sub_ids]))
          },
          pw_name_map[top_id]
        )
      )
  } else {
    gene_categories <- gene_categories %>%
      mutate(category = pw_name_map[top_id])
  }

  # 5. Filter to foreground
  fg_cats <- gene_categories %>%
    filter(gene %in% fg_genes) %>%
    distinct(gene, category)

  # 6. Specificity-weighted 1:1 assignment (prefer smallest fg category)
  fg_cat_counts <- fg_cats %>% count(category, name = "n_fg")

  best <- fg_cats %>%
    left_join(fg_cat_counts, by = "category") %>%
    arrange(n_fg) %>%
    group_by(gene) %>%
    slice_head(n = 1) %>%
    ungroup()

  # 7. Unmapped → "Other"
  unmapped <- setdiff(fg_genes, best$gene)
  if (length(unmapped) > 0) {
    best <- bind_rows(best,
      tibble(gene = unmapped, category = "Other", n_fg = NA_integer_))
  }

  # 8. Merge small categories
  small_cats <- best %>% count(category) %>%
    filter(n < min_cat_size, category != "Other") %>% pull(category)
  if (length(small_cats) > 0) {
    best <- best %>%
      mutate(category = ifelse(category %in% small_cats, "Other", category))
  }

  best %>%
    transmute(gene, category) %>%
    arrange(category, gene)
}

# === Standalone execution =====================================================

if (sys.nframe() == 0 || identical(environment(), globalenv())) {
  setwd(rprojroot::find_rstudio_root_file())

  message("\n=== Running Reactome classification (Option A) ===\n")

  # Load foreground proteins from F2 and F3
  dep_df <- read_csv("03_DEP/c_data/03_combined_results.csv", show_col_types = FALSE)
  all_genes <- unique(dep_df$gene)

  # F2 foreground: interaction + any training-significant proteins
  f2_fg <- dep_df %>%
    filter(pi_score_Training_Young < 0.05 | pi_score_Training_Old < 0.05 |
           pi_score_Interaction < 0.05) %>%
    pull(gene) %>% unique()

  # F3 foreground: aging + reversal significant proteins
  f3_fg <- dep_df %>%
    filter(pi_score_Aging < 0.05 | pi_score_Reversal < 0.05) %>%
    pull(gene) %>% unique()

  message(sprintf("F2 foreground: %d proteins", length(f2_fg)))
  message(sprintf("F3 foreground: %d proteins", length(f3_fg)))

  # Classify
  f2_react <- assign_reactome_toplevel(f2_fg, all_genes)
  f3_react <- assign_reactome_toplevel(f3_fg, all_genes)

  # Save
  out_dir <- "04_Figures_v1/shared/data"
  write_csv(f2_react, file.path(out_dir, "reactome_classification_f2.csv"))
  write_csv(f3_react, file.path(out_dir, "reactome_classification_f3.csv"))

  # Coverage stats
  f2_stats <- f2_react %>%
    summarise(
      figure = "F2",
      n_total = n(),
      n_mapped = sum(category != "Other"),
      pct_mapped = round(100 * n_mapped / n_total, 1),
      n_categories = n_distinct(category[category != "Other"]),
      largest_cat = names(which.max(table(category[category != "Other"]))),
      largest_pct = round(100 * max(table(category[category != "Other"])) / n_total, 1),
      smallest_cat_n = min(table(category[category != "Other"]))
    )
  f3_stats <- f3_react %>%
    summarise(
      figure = "F3",
      n_total = n(),
      n_mapped = sum(category != "Other"),
      pct_mapped = round(100 * n_mapped / n_total, 1),
      n_categories = n_distinct(category[category != "Other"]),
      largest_cat = names(which.max(table(category[category != "Other"]))),
      largest_pct = round(100 * max(table(category[category != "Other"])) / n_total, 1),
      smallest_cat_n = min(table(category[category != "Other"]))
    )
  stats <- bind_rows(f2_stats, f3_stats)
  write_csv(stats, file.path(out_dir, "reactome_coverage_stats.csv"))

  message("\n--- F2 Reactome categories ---")
  print(table(f2_react$category))
  message("\n--- F3 Reactome categories ---")
  print(table(f3_react$category))
  message(sprintf("\nF2 coverage: %.1f%% (%d/%d), %d categories",
                  f2_stats$pct_mapped, f2_stats$n_mapped, f2_stats$n_total,
                  f2_stats$n_categories))
  message(sprintf("F3 coverage: %.1f%% (%d/%d), %d categories",
                  f3_stats$pct_mapped, f3_stats$n_mapped, f3_stats$n_total,
                  f3_stats$n_categories))
}
