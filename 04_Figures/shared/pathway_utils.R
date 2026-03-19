# Unified Pathway Enrichment Utility
# Post-hoc Jaccard deduplication (Reimand et al. 2019, Nat Protocols)
# Sources: MSigDB Hallmark (H), Canonical Pathways (C2:CP), GO:BP (C5:GO:BP)

#' Greedy Jaccard deduplication of enrichment results
#'
#' Sorts results by padj ascending. For each term, checks Jaccard overlap with
#' all previously kept terms. If Jaccard > cutoff with any kept term, the term
#' is dropped (the more-significant one was already kept).
#'
#' @param results tibble with at least columns: pathway, padj
#' @param pathways named list of gene sets (character vectors)
#' @param jaccard_cutoff numeric, drop if Jaccard > this (default 0.5)
#' @return filtered tibble with redundant terms removed
deduplicate_enrichment <- function(results, pathways, jaccard_cutoff = 0.5) {
  if (nrow(results) == 0) return(results)

  results <- results[order(results$padj), ]
  kept_names <- character(0)
  kept_sets  <- list()
  keep_mask  <- logical(nrow(results))

  for (i in seq_len(nrow(results))) {
    pw_name <- results$pathway[i]
    pw_genes <- pathways[[pw_name]]
    if (is.null(pw_genes)) { keep_mask[i] <- TRUE; next }

    is_redundant <- FALSE
    for (j in seq_along(kept_sets)) {
      inter <- length(intersect(pw_genes, kept_sets[[j]]))
      union <- length(union(pw_genes, kept_sets[[j]]))
      if (union > 0 && (inter / union) > jaccard_cutoff) {
        is_redundant <- TRUE
        break
      }
    }

    if (!is_redundant) {
      keep_mask[i] <- TRUE
      kept_names <- c(kept_names, pw_name)
      kept_sets[[length(kept_sets) + 1]] <- pw_genes
    }
  }

  results[keep_mask, ]
}


#' Build unified pathway collection from MSigDB
#'
#' Combines Hallmark (H), KEGG Medicus, Reactome, WikiPathways, and GO:BP.
#' Optionally includes GO Slim gene sets. Filters disease/cancer terms from
#' curated pathways. Applies size filters.
#'
#' @param species character, default "Homo sapiens"
#' @param min_size integer, minimum gene set size (default 10)
#' @param max_size integer, maximum gene set size (default 500)
#' @param include_goslim logical, include GO Slim gene sets (default TRUE)
#' @return named list of character vectors (pathway name -> gene symbols)
build_pathway_collection <- function(species = "Homo sapiens",
                                     min_size = 10, max_size = 500,
                                     include_goslim = TRUE,
                                     exclude_variants = FALSE,
                                     exclude_wp = FALSE) {
  requireNamespace("msigdbr", quietly = TRUE)

  hallmark <- msigdbr::msigdbr(species = species, collection = "H")
  kegg     <- msigdbr::msigdbr(species = species, collection = "C2",
                                subcollection = "CP:KEGG_MEDICUS")
  reactome <- msigdbr::msigdbr(species = species, collection = "C2",
                                subcollection = "CP:REACTOME")
  wp       <- msigdbr::msigdbr(species = species, collection = "C2",
                                subcollection = "CP:WIKIPATHWAYS")
  gobp     <- msigdbr::msigdbr(species = species, collection = "C5",
                                subcollection = "GO:BP")

  disease_pat <- paste0("DISEASE|CANCER|TUMOR|CARCINOMA|LEUKEMIA|LYMPHOMA|",
                        "MELANOMA|GLIOMA|HEPATITIS|HIV|INFECTION|VIRAL|",
                        "BACTERIAL|PARASIT")
  kegg     <- kegg[!grepl(disease_pat, kegg$gs_name, ignore.case = TRUE), ]
  reactome <- reactome[!grepl(disease_pat, reactome$gs_name, ignore.case = TRUE), ]
  wp       <- wp[!grepl(disease_pat, wp$gs_name, ignore.case = TRUE), ]

  if (exclude_variants) {
    kegg <- kegg[!grepl("_VARIANT_", kegg$gs_name), ]
  }

  cols <- c("gs_name", "gene_symbol")
  sets_list <- list(hallmark[, cols], kegg[, cols], reactome[, cols], gobp[, cols])
  dbs <- c("H", "KEGG", "Reactome", "GO:BP")
  if (!exclude_wp) {
    sets_list <- c(sets_list, list(wp[, cols]))
    dbs <- c(dbs, "WP")
  }
  all_sets <- do.call(rbind, sets_list)

  pw_list <- split(all_sets$gene_symbol, all_sets$gs_name)
  pw_list <- lapply(pw_list, unique)

  if (include_goslim) {
    goslim_sets <- build_goslim_gene_sets(
      species = species, min_size = min_size, max_size = max_size
    )
    pw_list <- c(pw_list, goslim_sets)
    dbs <- c(dbs, "GO Slim")
  }

  sizes <- vapply(pw_list, length, integer(1))
  pw_list <- pw_list[sizes >= min_size & sizes <= max_size]

  message(sprintf("Pathway collection: %d sets (%s), size %d-%d",
                  length(pw_list), paste(dbs, collapse = " + "),
                  min_size, max_size))
  pw_list
}


#' Build GO Slim gene sets from GO.db hierarchy
#'
#' Converts 62 GO Slim Generic BP terms into fgsea-compatible named gene set
#' lists. For each slim term, collects all genes annotated to it or any
#' descendant term via GOBPOFFSPRING.
#'
#' @param species character (unused, included for API consistency)
#' @param min_size integer, minimum gene set size (default 10)
#' @param max_size integer, maximum gene set size (default 500)
#' @return named list of character vectors (GOSLIM_* -> gene symbols)
build_goslim_gene_sets <- function(species = "Homo sapiens",
                                   min_size = 10, max_size = 500) {
  requireNamespace("GO.db", quietly = TRUE)
  requireNamespace("org.Hs.eg.db", quietly = TRUE)
  requireNamespace("AnnotationDbi", quietly = TRUE)

  # 62 GO Slim Generic BP terms (from go_slim_categories.R)
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

  # Get all descendant GO terms for each slim term
  offspring <- as.list(GO.db::GOBPOFFSPRING)

  # Map all BP GO terms -> gene symbols via org.Hs.eg.db
  suppressMessages({
    go_genes <- AnnotationDbi::select(
      org.Hs.eg.db::org.Hs.eg.db,
      keys = AnnotationDbi::keys(org.Hs.eg.db::org.Hs.eg.db, keytype = "GO"),
      keytype = "GO",
      columns = c("SYMBOL", "ONTOLOGY")
    )
  })
  go_bp_genes <- go_genes[!is.na(go_genes$ONTOLOGY) & go_genes$ONTOLOGY == "BP", ]
  go_to_symbols <- split(go_bp_genes$SYMBOL, go_bp_genes$GO)

  # Build gene sets: each slim term + all its descendants
  goslim_sets <- list()
  slim_names <- vapply(bp_slim, function(id) {
    tryCatch(AnnotationDbi::Term(GO.db::GOTERM[[id]]),
             error = function(e) NA_character_)
  }, character(1))

  for (i in seq_along(bp_slim)) {
    go_id <- bp_slim[i]
    go_term <- slim_names[i]
    if (is.na(go_term)) next

    # Collect genes from this term + all offspring
    all_terms <- go_id
    desc <- offspring[[go_id]]
    if (!is.null(desc)) all_terms <- c(all_terms, desc)

    genes <- unique(unlist(go_to_symbols[intersect(all_terms, names(go_to_symbols))],
                           use.names = FALSE))
    genes <- genes[!is.na(genes)]

    if (length(genes) >= min_size && length(genes) <= max_size) {
      set_name <- paste0("GOSLIM_", toupper(gsub(" ", "_", go_term)))
      goslim_sets[[set_name]] <- genes
    }
  }

  message(sprintf("GO Slim: %d/%d terms passed size filter (%d-%d)",
                  length(goslim_sets), length(bp_slim), min_size, max_size))
  goslim_sets
}


#' Run fGSEA on unified pathway collection with post-hoc deduplication
#'
#' @param ranks named numeric vector of gene-level statistics (e.g., t-stats)
#' @param pathways named list from build_pathway_collection()
#' @param jaccard_cutoff Jaccard threshold for dedup (default 0.5)
#' @param nperm integer, nPermSimple for fgseaMultilevel (default 10000)
#' @param min_size integer (default 15)
#' @param max_size integer (default 500)
#' @return tibble with columns: pathway, padj, NES, size, leadingEdge, database
run_fgsea_deduplicated <- function(ranks, pathways, jaccard_cutoff = 0.5,
                                   nperm = 10000, min_size = 15,
                                   max_size = 500) {
  requireNamespace("fgsea", quietly = TRUE)

  res <- fgsea::fgseaMultilevel(
    pathways    = pathways,
    stats       = ranks,
    minSize     = min_size,
    maxSize     = max_size,
    nPermSimple = nperm,
    eps         = 0
  )
  res <- as.data.frame(res)

  res$database <- classify_database(res$pathway)

  res <- tibble::as_tibble(res)

  keep_cols <- c("pathway", "padj", "NES", "size", "leadingEdge",
                 "database", "pval", "ES", "log2err")
  res <- res[, intersect(keep_cols, names(res))]

  sig   <- res[!is.na(res$padj) & res$padj < 0.05, ]
  nonsig <- res[is.na(res$padj) | res$padj >= 0.05, ]

  sig_dedup <- deduplicate_enrichment(sig, pathways, jaccard_cutoff)

  n_removed <- nrow(sig) - nrow(sig_dedup)
  pct <- if (nrow(sig) > 0) round(100 * n_removed / nrow(sig), 1) else 0
  message(sprintf("fGSEA dedup: %d sig -> %d kept (removed %d, %.1f%%)",
                  nrow(sig), nrow(sig_dedup), n_removed, pct))

  rbind(sig_dedup, nonsig)
}


#' Multi-contrast fGSEA with collapsePathways + Jaccard dedup
#'
#' Runs fgseaMultilevel per contrast, applies collapsePathways() per contrast,
#' then Jaccard dedup at a given cutoff. Returns tidy tibble with union of
#' surviving pathways across all contrasts, including NES values for non-significant
#' contrasts (needed for heatmap display).
#'
#' @param stats_list named list of named numeric vectors (contrast -> gene stats)
#' @param pw_list named list of gene sets from build_pathway_collection()
#' @param jaccard_cutoff Jaccard threshold for greedy dedup (default 0.35)
#' @param nperm nPermSimple for fgseaMultilevel (default 10000)
#' @param min_size minimum gene set size (default 15)
#' @param max_size maximum gene set size (default 500)
#' @param padj_cutoff significance threshold (default 0.05)
#' @return list with: long_df (all results, union pathways), sig_union (pathway names)
run_enrichment_pipeline <- function(stats_list, pw_list,
                                    jaccard_cutoff = 0.35,
                                    nperm = 10000,
                                    min_size = 15, max_size = 500,
                                    padj_cutoff = 0.05) {
  requireNamespace("fgsea", quietly = TRUE)

  all_results <- list()

  for (ctr in names(stats_list)) {
    message(sprintf("\n--- %s ---", ctr))
    ranks <- stats_list[[ctr]]

    # Run fgsea (returns data.table)
    res_dt <- fgsea::fgseaMultilevel(
      pathways    = pw_list,
      stats       = ranks,
      minSize     = min_size,
      maxSize     = max_size,
      nPermSimple = nperm,
      eps         = 0
    )

    # collapsePathways needs data.table input
    sig_dt <- res_dt[!is.na(res_dt$padj) & res_dt$padj < padj_cutoff, ]
    if (nrow(sig_dt) > 0) {
      collapsed <- fgsea::collapsePathways(
        fgseaRes     = sig_dt,
        pathways     = pw_list,
        stats        = ranks
      )
      independent <- collapsed$mainPathways
      message(sprintf("collapsePathways: %d sig -> %d independent",
                      nrow(sig_dt), length(independent)))
      # Mark non-independent sig pathways as padj = 1 (effectively removes them)
      drop_pw <- setdiff(sig_dt$pathway, independent)
      if (length(drop_pw) > 0) {
        res_dt$padj[res_dt$pathway %in% drop_pw] <- 1
      }
    }

    # Convert to data.frame for Jaccard dedup
    res <- as.data.frame(res_dt)
    res$database <- classify_database(res$pathway)
    res$contrast <- ctr

    # Jaccard dedup on remaining sig
    sig_after <- res[!is.na(res$padj) & res$padj < padj_cutoff, ]
    sig_dedup <- deduplicate_enrichment(sig_after, pw_list, jaccard_cutoff)
    n_removed <- nrow(sig_after) - nrow(sig_dedup)
    message(sprintf("Jaccard dedup (%.2f): %d -> %d (removed %d)",
                    jaccard_cutoff, nrow(sig_after), nrow(sig_dedup), n_removed))

    # Reset padj for terms that didn't survive dedup
    survived <- sig_dedup$pathway
    dedup_drop <- setdiff(sig_after$pathway, survived)
    if (length(dedup_drop) > 0) {
      res$padj[res$pathway %in% dedup_drop] <- 1
    }

    all_results[[ctr]] <- tibble::as_tibble(res)
  }

  long_df <- dplyr::bind_rows(all_results)

  # Union of surviving sig pathways across all contrasts
  sig_union <- unique(long_df$pathway[!is.na(long_df$padj) & long_df$padj < padj_cutoff])
  message(sprintf("\nUnion of sig pathways: %d", length(sig_union)))

  # Filter to union pathways only
  long_df <- long_df[long_df$pathway %in% sig_union, ]

  # Summary
  for (ctr in names(stats_list)) {
    sub <- long_df[long_df$contrast == ctr, ]
    n_sig <- sum(!is.na(sub$padj) & sub$padj < padj_cutoff)
    n_up  <- sum(!is.na(sub$padj) & sub$padj < padj_cutoff & sub$NES > 0)
    n_dn  <- sum(!is.na(sub$padj) & sub$padj < padj_cutoff & sub$NES < 0)
    message(sprintf("  %s: %d sig (%d up, %d down)", ctr, n_sig, n_up, n_dn))
  }

  list(long_df = long_df, sig_union = sig_union)
}


#' Run over-representation analysis with post-hoc deduplication
#'
#' Uses fgsea::fora() for hypergeometric test across all collections at once,
#' with a single BH correction. Then deduplicates.
#'
#' @param genes character vector of gene symbols (foreground set)
#' @param universe character vector of gene symbols (background)
#' @param pathways named list from build_pathway_collection()
#' @param jaccard_cutoff Jaccard threshold (default 0.5)
#' @param min_size integer (default 10)
#' @param max_size integer (default 500)
#' @param padj_cutoff numeric, significance threshold (default 0.05)
#' @return tibble with columns: pathway, padj, overlap, size, overlapGenes, database
run_ora_deduplicated <- function(genes, universe, pathways,
                                 jaccard_cutoff = 0.5,
                                 min_size = 10, max_size = 500,
                                 padj_cutoff = 0.05) {
  requireNamespace("fgsea", quietly = TRUE)

  genes <- intersect(genes, universe)

  res <- fgsea::fora(
    pathways = pathways,
    genes    = genes,
    universe = universe,
    minSize  = min_size,
    maxSize  = max_size
  )
  res <- as.data.frame(res)
  res$database <- classify_database(res$pathway)

  N <- length(universe)
  K <- length(genes)
  res$odds_ratio <- vapply(seq_len(nrow(res)), function(i) {
    a <- res$overlap[i]          # hits in pathway
    b <- K - a                   # foreground not in pathway
    c <- res$size[i] - a         # pathway not in foreground
    d <- N - K - c               # neither
    if (b == 0 || c == 0) Inf else (a * d) / (b * c)
  }, numeric(1))

  res <- tibble::as_tibble(res)

  sig <- res[!is.na(res$padj) & res$padj < padj_cutoff, ]
  sig_dedup <- deduplicate_enrichment(sig, pathways, jaccard_cutoff)

  n_removed <- nrow(sig) - nrow(sig_dedup)
  pct <- if (nrow(sig) > 0) round(100 * n_removed / nrow(sig), 1) else 0
  message(sprintf("ORA dedup: %d sig -> %d kept (removed %d, %.1f%%)",
                  nrow(sig), nrow(sig_dedup), n_removed, pct))

  sig_dedup
}


#' Assign genes to enrichment-derived pathway categories
#'
#' Runs ORA, deduplicates, takes top N surviving pathways as categories.
#' Each gene is assigned to its most specific (smallest) enriched category.
#' Unassigned genes get "Other".
#'
#' @param genes character vector of gene symbols
#' @param universe character vector of background gene symbols
#' @param pathways named list from build_pathway_collection()
#' @param jaccard_cutoff Jaccard threshold (default 0.5)
#' @param max_categories integer, max pathway categories to use (default 15)
#' @return tibble with columns: gene, pathway_category
assign_enrichment_classes <- function(genes, universe, pathways,
                                      jaccard_cutoff = 0.5,
                                      max_categories = 15,
                                      padj_cutoff = 0.05) {
  ora_res <- run_ora_deduplicated(
    genes    = genes,
    universe = universe,
    pathways = pathways,
    jaccard_cutoff = jaccard_cutoff,
    min_size = 10, max_size = 500,
    padj_cutoff = padj_cutoff
  )

  if (nrow(ora_res) == 0) {
    message("No significant enrichment found; all genes classified as 'Other'")
    return(tibble::tibble(gene = genes, pathway_category = "Other"))
  }

  top_cats <- head(ora_res[order(ora_res$padj), ], max_categories)

  top_cats <- top_cats[order(top_cats$size), ]
  gene_class <- rep(NA_character_, length(genes))
  names(gene_class) <- genes

  for (i in seq_len(nrow(top_cats))) {
    pw_name <- top_cats$pathway[i]
    pw_genes <- pathways[[pw_name]]
    hits <- intersect(genes, pw_genes)
    # Only assign if not already assigned (most-specific-first)
    unassigned_hits <- hits[is.na(gene_class[hits])]
    if (length(unassigned_hits) > 0) {
      gene_class[unassigned_hits] <- pw_name
    }
  }

  gene_class[is.na(gene_class)] <- "Other"

  result <- tibble::tibble(
    gene = names(gene_class),
    pathway_category = unname(gene_class)
  )

  n_assigned <- sum(result$pathway_category != "Other")
  pct <- round(100 * n_assigned / nrow(result), 1)
  message(sprintf("Classification: %d/%d genes assigned (%.1f%%), %d categories + Other",
                  n_assigned, nrow(result), pct, nrow(top_cats)))

  result
}


#' Classify pathway name to database source
#' @param pathway_names character vector of MSigDB gs_name strings
#' @return character vector of database labels
classify_database <- function(pathway_names) {
  dplyr::case_when(
    grepl("^HALLMARK_",       pathway_names) ~ "Hallmark",
    grepl("^REACTOME_",       pathway_names) ~ "Reactome",
    grepl("^KEGG_MEDICUS_",   pathway_names) ~ "KEGG",
    grepl("^KEGG_",           pathway_names) ~ "KEGG",
    grepl("^WP_",             pathway_names) ~ "WikiPathways",
    grepl("^BIOCARTA_",       pathway_names) ~ "BioCarta",
    grepl("^PID_",            pathway_names) ~ "PID",
    grepl("^GOSLIM_",         pathway_names) ~ "GO Slim",
    grepl("^GOBP_",           pathway_names) ~ "GO:BP",
    grepl("^GOCC_",           pathway_names) ~ "GO:CC",
    grepl("^GOMF_",           pathway_names) ~ "GO:MF",
    TRUE ~ "Other"
  )
}


# MSigDB pathway ID -> 15 consolidated categories (keyword rules)
CONSOLIDATED_PATHWAY_ORDER <- c(
  "Muscle & Contractile", "Cytoskeleton & Motility", "ECM & Adhesion",
  "Lipid Metabolism", "Carbohydrate & Energy Metabolism",
  "Amino Acid & Cofactor Metabolism",
  "Mitochondria & Energy", "Protein Homeostasis",
  "Transport", "Translation & Ribosome", "Transcription & Chromatin",
  "Immune & Inflammation", "DNA & Cell Cycle", "Circulatory System",
  "Development", "Other"
)

CONSOLIDATED_COLORS <- c(
  "Muscle & Contractile"              = "#E57373",
  "Cytoskeleton & Motility"           = "#FFB74D",
  "ECM & Adhesion"                    = "#FFF176",
  "Lipid Metabolism"                  = "#AED581",
  "Carbohydrate & Energy Metabolism"  = "#81C784",
  "Amino Acid & Cofactor Metabolism"  = "#66BB6A",
  "Mitochondria & Energy"             = "#4DB6AC",
  "Protein Homeostasis"               = "#4FC3F7",
  "Transport"                         = "#7986CB",
  "Translation & Ribosome"            = "#BA68C8",
  "Transcription & Chromatin"         = "#AB47BC",
  "Immune & Inflammation"             = "#A1887F",
  "DNA & Cell Cycle"                  = "#90A4AE",
  "Circulatory System"                = "#CE93D8",
  "Development"                       = "#B0BEC5",
  "Other"                             = "#D0D0D0"
)

#' Keyword-based classifier: MSigDB pathway ID -> consolidated category
#'
#' @param ids character vector of msigdbr term IDs (gs_name strings)
#' @return character vector of category names (same length as ids)
classify_pathway_func <- function(ids) {
  rules <- list(
    "Muscle & Contractile"              = "MYOGEN|MYOFIBRIL|SARCOMERE|MUSCLE_|CONTRACTILE|ACTOMYOSIN|MYOSIN|I_BAND",
    "Cytoskeleton & Motility"           = "CYTOSKELET|ACTIN_BIND|STRUCTURAL_MOLECULE|MOTIL|SUPRAMOLECUL",
    "ECM & Adhesion"                    = "EXTRACELLULAR_MATRIX|COLLAGEN|BASEMENT_MEMBRANE|ADHESION|APICAL_JUNCTION|EMT|ENCAPSULATING",
    "Lipid Metabolism"                  = "FATTY_ACID|LIPID|ADIPOGEN|STEROID|SPHINGOLIPID|PHOSPHOLIPID|KETONE",
    "Carbohydrate & Energy Metabolism"  = "GLYCOLY|GLUCONEO|CARBOHYDRATE|PENTOSE|PRECURSOR_METABOL",
    "Amino Acid & Cofactor Metabolism"  = "AMINO_ACID|VITAMIN|COFACTOR|NITROGEN|DETOXIF|DIGEST|XENOBIOT",
    "Mitochondria & Energy"             = "MITOCHOND|OXIDATIVE_PHOSPH|ELECTRON_TRANSFER|RESPIRATORY|OXIDOREDUCT",
    "Protein Homeostasis"               = "PROTEASOM|UBIQUITIN|AUTOPHAGY|MTORC1|PROTEIN_FOLD",
    "Transport"                         = "TRANSPORT(?!.*ELECTRON)|VESICLE|ENDOCYT|SECRETI",
    "Translation & Ribosome"            = "TRANSLAT|RIBOSOM|TRNA|MYC_TARGET",
    "Transcription & Chromatin"         = "TRANSCRIPT|SPLICEOSOM|E2F_TARGET|CHROMATIN|MRNA_PROC",
    "Immune & Inflammation"             = "IMMUN|INFLAMMA|INTERFERON|IL2|IL6|TNFA|NF.KB|COMPLEMENT",
    "DNA & Cell Cycle"                  = "DNA_REPAIR|CELL_CYCLE|MITOTIC|P53_PATHWAY",
    "Circulatory System"                = "ANGIOGEN|BLOOD_VESSEL|HYPOXIA",
    "Development"                       = "UV_RESPONSE|GROWTH_FACTOR|WNT|HEDGEHOG|NOTCH|TGF_BETA|KRAS"
  )
  vapply(toupper(ids), function(id) {
    for (cat in names(rules)) {
      if (grepl(rules[[cat]], id, perl = TRUE)) return(cat)
    }
    if (grepl("METABOL", toupper(id), perl = TRUE)) return("Amino Acid & Cofactor Metabolism")
    "Other"
  }, character(1), USE.NAMES = FALSE)
}


assign_consolidated <- function(fg_genes, all_genes, pathways = NULL,
                                max_categories = 25, jaccard_cutoff = 0.5,
                                padj_cutoff = 0.05) {
  if (is.null(pathways)) pathways <- build_pathway_collection()

  gene_classes <- assign_enrichment_classes(
    fg_genes, all_genes, pathways, jaccard_cutoff, max_categories,
    padj_cutoff = padj_cutoff
  )

  gene_classes$consolidated <- classify_pathway_func(gene_classes$pathway_category)
  gene_classes$consolidated[gene_classes$pathway_category == "Other"] <- "Other"
  gene_classes$consolidated <- factor(
    gene_classes$consolidated, levels = CONSOLIDATED_PATHWAY_ORDER
  )

  gene_map <- gene_classes |>
    dplyr::transmute(gene, pathway = consolidated)

  pw_sizes <- vapply(pathways, length, integer(1))
  pw_cats  <- classify_pathway_func(names(pathways))
  names(pw_cats) <- names(pathways)

  bg_best_cat  <- rep("Other", length(all_genes))
  bg_best_size <- rep(Inf, length(all_genes))
  names(bg_best_cat)  <- all_genes
  names(bg_best_size) <- all_genes

  for (pw_name in names(pathways)) {
    hits <- intersect(pathways[[pw_name]], all_genes)
    if (length(hits) == 0) next
    s <- pw_sizes[pw_name]
    better <- hits[bg_best_size[hits] > s]
    if (length(better) > 0) {
      bg_best_cat[better]  <- pw_cats[pw_name]
      bg_best_size[better] <- s
    }
  }

  bg_map <- tibble::tibble(
    gene    = all_genes,
    pathway = factor(unname(bg_best_cat), levels = CONSOLIDATED_PATHWAY_ORDER)
  )

  active_cats <- unique(as.character(gene_map$pathway))
  active_cats <- active_cats[active_cats != "Other"]

  fisher_results <- tibble::tibble(
    pathway_label = active_cats,
    pvalue = vapply(active_cats, function(s) {
      fg_in  <- sum(gene_map$pathway == s)
      fg_out <- nrow(gene_map) - fg_in
      bg_in  <- sum(bg_map$pathway == s)
      bg_out <- nrow(bg_map) - bg_in
      stats::fisher.test(
        matrix(c(fg_in, bg_in, fg_out, bg_out), 2, 2),
        alternative = "greater"
      )$p.value
    }, numeric(1))
  ) |>
    dplyr::mutate(
      p.adjust = stats::p.adjust(pvalue, method = "BH"),
      ID       = pathway_label,
      database = "Consolidated"
    )

  if ("Other" %in% as.character(gene_map$pathway)) {
    fisher_results <- dplyr::bind_rows(fisher_results, tibble::tibble(
      pathway_label = "Other", pvalue = 1, p.adjust = 1,
      ID = "OTHER", database = "Other"
    ))
  }

  list(gene_map = gene_map, bg_map = bg_map, ora = fisher_results)
}


assign_pathways_membership <- function(fg_genes, universe, pathways = NULL,
                                       max_pathways = 12, min_overlap = 2,
                                       jaccard_cutoff = 0.5,
                                       min_size = 10, max_size = 500) {
  requireNamespace("fgsea", quietly = TRUE)

  if (is.null(pathways)) pathways <- build_pathway_collection(min_size = min_size,
                                                               max_size = max_size)

  fg_genes <- intersect(fg_genes, universe)

  res <- fgsea::fora(
    pathways = pathways,
    genes    = fg_genes,
    universe = universe,
    minSize  = min_size,
    maxSize  = max_size
  )
  res <- tibble::as_tibble(as.data.frame(res))
  res$database <- classify_database(res$pathway)

  res <- res[res$overlap >= min_overlap, ]
  res <- res[order(res$pval), ]

  message(sprintf("Membership mapping: %d pathways with overlap >= %d",
                  nrow(res), min_overlap))

  kept_names <- character(0)
  kept_sets  <- list()
  keep_mask  <- logical(nrow(res))

  for (i in seq_len(nrow(res))) {
    pw_name <- res$pathway[i]
    pw_genes <- pathways[[pw_name]]
    if (is.null(pw_genes)) { keep_mask[i] <- TRUE; next }

    is_redundant <- FALSE
    for (j in seq_along(kept_sets)) {
      inter <- length(intersect(pw_genes, kept_sets[[j]]))
      uni   <- length(union(pw_genes, kept_sets[[j]]))
      if (uni > 0 && (inter / uni) > jaccard_cutoff) {
        is_redundant <- TRUE
        break
      }
    }

    if (!is_redundant) {
      keep_mask[i] <- TRUE
      kept_names <- c(kept_names, pw_name)
      kept_sets[[length(kept_sets) + 1]] <- pw_genes
    }
  }

  res_dedup <- res[keep_mask, ]
  message(sprintf("After Jaccard dedup (cutoff=%.2f): %d → %d pathways",
                  jaccard_cutoff, nrow(res), nrow(res_dedup)))

  top_pw <- head(res_dedup, max_pathways)
  message(sprintf("Selected top %d pathways", nrow(top_pw)))

  top_sorted <- top_pw[order(top_pw$size), ]
  gene_class <- setNames(rep(NA_character_, length(fg_genes)), fg_genes)

  for (i in seq_len(nrow(top_sorted))) {
    pw_name <- top_sorted$pathway[i]
    hits <- intersect(fg_genes, pathways[[pw_name]])
    unassigned <- hits[is.na(gene_class[hits])]
    if (length(unassigned) > 0) gene_class[unassigned] <- pw_name
  }
  gene_class[is.na(gene_class)] <- "Other"

  n_mapped <- sum(gene_class != "Other")
  pct_mapped <- round(100 * n_mapped / length(fg_genes), 1)
  message(sprintf("Gene assignment: %d/%d (%.1f%%) mapped",
                  n_mapped, length(fg_genes), pct_mapped))

  gene_map <- tibble::tibble(
    gene       = names(gene_class),
    pathway_id = unname(gene_class),
    pathway    = ifelse(gene_class == "Other", "Other",
                        clean_pathway_name(gene_class)),
    database   = classify_database(gene_class)
  )

  ora_out <- top_pw %>%
    dplyr::mutate(pathway_label = clean_pathway_name(pathway)) %>%
    dplyr::select(pathway, pathway_label, pval, padj, overlap, size,
                  overlapGenes, database)

  list(gene_map = gene_map, ora = ora_out,
       n_mapped = n_mapped, pct_mapped = pct_mapped)
}
