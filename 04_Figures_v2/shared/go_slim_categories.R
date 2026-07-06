# GO Slim gene-to-category assignment for the F04/F05 Panel B bars.
# The curated terms, consolidation map, order, and colors are the single source
# of truth in go_slim_terms.R; this file only maps genes onto them.

requireNamespace("GO.db", quietly = TRUE)
requireNamespace("org.Hs.eg.db", quietly = TRUE)
requireNamespace("AnnotationDbi", quietly = TRUE)
library(dplyr)
library(tidyr)

source("04_Figures_v2/shared/go_slim_terms.R")

assign_go_slim_consolidated <- function(fg_genes, all_genes, min_cat_size = 2) {
  suppressMessages({
    all_entrez <- AnnotationDbi::mapIds(org.Hs.eg.db::org.Hs.eg.db,
      keys = all_genes,
      keytype = "SYMBOL", column = "ENTREZID",
      multiVals = "first"
    )
    all_go <- AnnotationDbi::select(org.Hs.eg.db::org.Hs.eg.db,
      keys = as.character(na.omit(all_entrez)),
      keytype = "ENTREZID",
      columns = c("SYMBOL", "GO", "ONTOLOGY")
    )
  })
  # re-attach dplyr after AnnotationDbi::select() masks it
  suppressPackageStartupMessages(library(dplyr))
  all_bp <- all_go |>
    filter(ONTOLOGY == "BP", !is.na(GO)) |>
    distinct(SYMBOL, GO)

  ancestors <- as.list(GO.db::GOBPANCESTOR)
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

  all_gene_slim <- all_bp |>
    mutate(slim_list = go_to_slim[GO]) |>
    unnest(slim_list) |>
    select(SYMBOL, slim = slim_list) |>
    distinct()

  all_gene_consolidated <- all_gene_slim |>
    mutate(consolidated = SLIM_CONSOLIDATED[slim]) |>
    filter(!is.na(consolidated)) |>
    distinct(SYMBOL, consolidated)

  fg_gene_slim <- all_gene_slim |> filter(SYMBOL %in% fg_genes)
  fg_gene_consolidated <- fg_gene_slim |>
    mutate(consolidated = SLIM_CONSOLIDATED[slim]) |>
    filter(!is.na(consolidated))

  fg_term_counts <- fg_gene_slim |> count(slim, name = "n_fg")

  best_consolidated <- fg_gene_consolidated |>
    left_join(fg_term_counts, by = "slim") |>
    mutate(priority = ifelse(consolidated == "Development", 2, 1)) |>
    arrange(priority, n_fg) |>
    group_by(SYMBOL) |>
    slice_head(n = 1) |>
    ungroup()

  unmapped <- setdiff(fg_genes, best_consolidated$SYMBOL)
  if (length(unmapped)) {
    best_consolidated <- bind_rows(
      best_consolidated,
      tibble(
        SYMBOL = unmapped, slim = "OTHER", consolidated = "Other",
        n_fg = NA_integer_, priority = 3L
      )
    )
  }

  small_cats <- best_consolidated |>
    count(consolidated) |>
    filter(n < min_cat_size, consolidated != "Other") |>
    pull(consolidated)
  if (length(small_cats)) {
    best_consolidated <- best_consolidated |>
      mutate(consolidated = ifelse(consolidated %in% small_cats, "Other", consolidated))
  }

  best_consolidated |>
    transmute(gene = SYMBOL, slim, consolidated) |>
    mutate(consolidated = factor(consolidated, levels = CONSOLIDATED_PATHWAY_ORDER))
}
