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
#' Combines Hallmark (H), Canonical Pathways (C2:CP), and GO:BP (C5:GO:BP).
#' Filters disease/cancer terms from C2:CP. Applies size filters.
#'
#' @param species character, default "Homo sapiens"
#' @param min_size integer, minimum gene set size (default 10)
#' @param max_size integer, maximum gene set size (default 500)
#' @return named list of character vectors (pathway name -> gene symbols)
build_pathway_collection <- function(species = "Homo sapiens",
                                     min_size = 10, max_size = 500) {
  requireNamespace("msigdbr", quietly = TRUE)

  hallmark <- msigdbr::msigdbr(species = species, collection = "H")
  c2_cp    <- msigdbr::msigdbr(species = species, collection = "C2",
                                subcollection = "CP")
  gobp     <- msigdbr::msigdbr(species = species, collection = "C5",
                                subcollection = "GO:BP")

  disease_pat <- paste0("DISEASE|CANCER|TUMOR|CARCINOMA|LEUKEMIA|LYMPHOMA|",
                        "MELANOMA|GLIOMA|HEPATITIS|HIV|INFECTION|VIRAL|",
                        "BACTERIAL|PARASIT")
  c2_cp <- c2_cp[!grepl(disease_pat, c2_cp$gs_name, ignore.case = TRUE), ]

  all_sets <- rbind(
    hallmark[, c("gs_name", "gene_symbol")],
    c2_cp[, c("gs_name", "gene_symbol")],
    gobp[, c("gs_name", "gene_symbol")]
  )

  pw_list <- split(all_sets$gene_symbol, all_sets$gs_name)
  pw_list <- lapply(pw_list, unique)

  sizes <- vapply(pw_list, length, integer(1))
  pw_list <- pw_list[sizes >= min_size & sizes <= max_size]

  message(sprintf("Pathway collection: %d sets (H + C2:CP + GO:BP), size %d-%d",
                  length(pw_list), min_size, max_size))
  pw_list
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
