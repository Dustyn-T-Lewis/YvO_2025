# ORA Pilot v2: Finding the Sweet Spot for Panel F Enrichment
# Tests database x dedup_level x input_list at Pi<0.05 only
# Key question: does GO Slim ORA give >80% gene coverage WITH statistical significance?
# Output: ora_pilot_v2_summary.csv

setwd(rprojroot::find_rstudio_root_file())

suppressPackageStartupMessages({
  library(tibble)
  library(readr)
})

source("04_Figures/shared/style.R")
source("04_Figures/shared/pathway_utils.R")
source("04_Figures/shared/go_slim_categories.R")

# Re-attach dplyr/tidyr AFTER go_slim_categories.R loads AnnotationDbi
suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
})

dep <- read_csv("03_DEP/c_data/03_combined_results.csv", show_col_types = FALSE)
universe <- dep$gene
set.seed(42)

# --- Section 1: Build All 4 Collections ------------------------------------

cat("Building collections...\n")

# 1. GO Slim gene sets (fora-compatible, ~30-40 sets)
goslim_sets <- build_goslim_gene_sets(min_size = 10, max_size = 500)

# 2. Comprehensive (no WP, no GO Slim)
comprehensive <- build_pathway_collection(
  min_size = 10, max_size = 500,
  include_goslim = FALSE, exclude_variants = TRUE, exclude_wp = TRUE
)

# 3. GO Slim + Hallmark
hallmark_raw <- msigdbr::msigdbr(species = "Homo sapiens", collection = "H")
hallmark_list <- split(hallmark_raw$gene_symbol, hallmark_raw$gs_name)
hallmark_list <- lapply(hallmark_list, unique)
goslim_hallmark <- c(goslim_sets, hallmark_list)

# 4. Full (GO Slim + Comprehensive)
full_collection <- c(goslim_sets, comprehensive)

COLLECTIONS <- list(
  "GO Slim"            = goslim_sets,
  "Comprehensive"      = comprehensive,
  "GO Slim + Hallmark" = goslim_hallmark,
  "Full"               = full_collection
)

for (nm in names(COLLECTIONS)) {
  cat(sprintf("  %s: %d sets\n", nm, length(COLLECTIONS[[nm]])))
}

# --- Section 2: Helper — run_ora_with_dedup() -------------------------------

#' Run ORA with configurable dedup level, returning enrichment + gene coverage
#' @param fg_genes character vector of foreground genes
#' @param universe character vector of background genes
#' @param collection named list of gene sets
#' @param dedup_level "none", "light" (Jaccard 0.35), or "standard" (Jaccard 0.50)
#' @return list with: results (tibble), n_genes_assigned, pct_genes_assigned
run_ora_with_dedup <- function(fg_genes, universe, collection, dedup_level = "none") {
  requireNamespace("fgsea", quietly = TRUE)

  fg_genes <- intersect(fg_genes, universe)
  n_fg <- length(fg_genes)

  # Run fora
  res <- fgsea::fora(
    pathways = collection,
    genes    = fg_genes,
    universe = universe,
    minSize  = 10,
    maxSize  = 500
  )
  res <- tibble::as_tibble(as.data.frame(res))

  # Compute odds ratios
  N <- length(universe)
  K <- n_fg
  res$odds_ratio <- vapply(seq_len(nrow(res)), function(i) {
    a <- res$overlap[i]
    b <- K - a
    cc <- res$size[i] - a
    d <- N - K - cc
    if (b == 0 || cc == 0) Inf else (a * d) / (b * cc)
  }, numeric(1))

  res$database <- classify_database(res$pathway)

  # Filter to significant

  sig <- res[!is.na(res$padj) & res$padj < 0.05, ]

  # Apply dedup
  if (dedup_level == "light" && nrow(sig) > 0) {
    sig <- deduplicate_enrichment(sig, collection, jaccard_cutoff = 0.35)
  } else if (dedup_level == "standard" && nrow(sig) > 0) {
    sig <- deduplicate_enrichment(sig, collection, jaccard_cutoff = 0.50)
  }

  # Gene coverage: how many fg_genes appear in ANY sig term's gene set
  if (nrow(sig) > 0) {
    sig_pw_genes <- unique(unlist(collection[sig$pathway], use.names = FALSE))
    assigned <- intersect(fg_genes, sig_pw_genes)
    n_genes_assigned <- length(assigned)
  } else {
    n_genes_assigned <- 0L
  }
  pct_genes_assigned <- if (n_fg > 0) round(100 * n_genes_assigned / n_fg, 1) else 0

  list(
    results            = sig,
    n_genes_assigned   = n_genes_assigned,
    pct_genes_assigned = pct_genes_assigned
  )
}

# --- Section 3: Build Foreground Lists (Pi < 0.05 only) --------------------

CONTRASTS <- c("Training_Young", "Training_Old", "Interaction")

build_fg_lists <- function(dep, granularity) {
  per_ctr <- setNames(
    lapply(CONTRASTS, function(ctr) {
      dep$gene[dep[[paste0("pi_score_", ctr)]] < 0.05]
    }),
    CONTRASTS
  )

  switch(granularity,
    "per_contrast" = per_ctr,
    "union" = list(Union = unique(unlist(per_ctr))),
    "directional" = {
      sig_TY  <- per_ctr[["Training_Young"]]
      sig_TO  <- per_ctr[["Training_Old"]]
      sig_Int <- per_ctr[["Interaction"]]

      lfc_TY <- setNames(dep$logFC_Training_Young, dep$gene)
      lfc_TO <- setNames(dep$logFC_Training_Old, dep$gene)

      both       <- intersect(sig_TY, sig_TO)
      young_only <- setdiff(sig_TY, c(sig_TO, sig_Int))
      old_only   <- setdiff(sig_TO, c(sig_TY, sig_Int))

      lists <- list(
        "Sig Both Up"    = both[lfc_TY[both] > 0],
        "Sig Both Down"  = both[lfc_TY[both] <= 0],
        "Sig Young Up"   = young_only[lfc_TY[young_only] > 0],
        "Sig Young Down" = young_only[lfc_TY[young_only] <= 0],
        "Sig Old Up"     = old_only[lfc_TO[old_only] > 0],
        "Sig Old Down"   = old_only[lfc_TO[old_only] <= 0],
        "Interaction"    = sig_Int
      )
      lists[vapply(lists, length, integer(1)) > 0]
    }
  )
}

# --- Section 4: Main Loop --------------------------------------------------

DEDUP_LEVELS  <- c("none", "light", "standard")
GRANULARITIES <- c("per_contrast", "union", "directional")

results <- list()
row_idx <- 0L

n_combos <- length(COLLECTIONS) * length(DEDUP_LEVELS) * length(GRANULARITIES)
cat(sprintf("\nRunning %d database x dedup x granularity combos (Pi < 0.05)...\n", n_combos))

for (db_name in names(COLLECTIONS)) {
  collection <- COLLECTIONS[[db_name]]

  for (dedup in DEDUP_LEVELS) {
    for (gran in GRANULARITIES) {
      fg_lists <- build_fg_lists(dep, gran)

      for (list_nm in names(fg_lists)) {
        fg <- fg_lists[[list_nm]]
        n_fg <- length(fg)
        row_idx <- row_idx + 1L

        cat(sprintf("  [%d] %s | %s | %s | %s (n=%d)\n",
                    row_idx, db_name, dedup, gran, list_nm, n_fg))

        # Skip if too few genes
        if (n_fg < 5) {
          results[[row_idx]] <- tibble(
            database = db_name, dedup_level = dedup, granularity = gran,
            list_name = list_nm, n_fg = n_fg,
            n_sig_terms = NA_integer_, n_themes_hit = NA_integer_,
            pct_coverage = NA_real_, n_genes_assigned = NA_integer_,
            pct_genes_assigned = NA_real_,
            top3_terms = NA_character_, median_odds_ratio = NA_real_,
            training_old_viable = NA
          )
          next
        }

        ora <- tryCatch(
          run_ora_with_dedup(fg, universe, collection, dedup),
          error = function(e) {
            message(sprintf("    ERROR: %s", conditionMessage(e)))
            NULL
          }
        )

        if (is.null(ora) || nrow(ora$results) == 0) {
          results[[row_idx]] <- tibble(
            database = db_name, dedup_level = dedup, granularity = gran,
            list_name = list_nm, n_fg = n_fg,
            n_sig_terms = 0L, n_themes_hit = 0L,
            pct_coverage = 0, n_genes_assigned = 0L,
            pct_genes_assigned = 0,
            top3_terms = NA_character_, median_odds_ratio = NA_real_,
            training_old_viable = NA
          )
          next
        }

        res <- ora$results
        n_sig <- nrow(res)

        # Themes via classify_pathway_func
        res$theme <- classify_pathway_func(res$pathway)
        n_themes <- n_distinct(res$theme[res$theme != "Other"])

        # Coverage: genes in any sig pathway / n_fg
        sig_pw_genes <- unique(unlist(collection[res$pathway], use.names = FALSE))
        covered <- length(intersect(fg, sig_pw_genes))
        pct_cov <- round(100 * covered / n_fg, 1)

        # Top 3 terms by odds ratio
        top3 <- res %>%
          arrange(desc(odds_ratio)) %>%
          head(3) %>%
          pull(pathway) %>%
          clean_pathway_name(max_chars = 35) %>%
          paste(collapse = "; ")

        med_or <- median(res$odds_ratio[is.finite(res$odds_ratio)], na.rm = TRUE)

        # Training_Old viability
        to_lists <- c("Training_Old", "Sig Old Up", "Sig Old Down")
        to_viable <- if (list_nm %in% to_lists) n_sig > 0 else NA

        results[[row_idx]] <- tibble(
          database = db_name, dedup_level = dedup, granularity = gran,
          list_name = list_nm, n_fg = n_fg,
          n_sig_terms = n_sig, n_themes_hit = n_themes,
          pct_coverage = pct_cov,
          n_genes_assigned = ora$n_genes_assigned,
          pct_genes_assigned = ora$pct_genes_assigned,
          top3_terms = top3, median_odds_ratio = round(med_or, 2),
          training_old_viable = to_viable
        )
      }
    }
  }
}

# --- Section 5: Output + Recommendation ------------------------------------

summary_df <- bind_rows(results) %>%
  dplyr::select(database, dedup_level, granularity, list_name,
                n_fg, n_sig_terms, n_themes_hit, pct_coverage,
                n_genes_assigned, pct_genes_assigned,
                top3_terms, median_odds_ratio, training_old_viable)

outpath <- "04_Figures/F04/c_data/supp_gallery/ora_pilot_v2_summary.csv"
write_csv(summary_df, outpath)
cat(sprintf("\nSaved: %s (%d rows)\n", outpath, nrow(summary_df)))

# Ranked summary: best combos by pct_genes_assigned * n_themes_hit
cat("\n=== Top 15 combos (ranked by pct_genes_assigned * n_themes_hit) ===\n")
summary_df %>%
  filter(!is.na(n_themes_hit), n_themes_hit > 0,
         !is.na(pct_genes_assigned), pct_genes_assigned > 0) %>%
  mutate(score = pct_genes_assigned * n_themes_hit) %>%
  arrange(desc(score)) %>%
  head(15) %>%
  mutate(label = sprintf("  %s | %s | %s | %s: %d terms, %d themes, %.0f%% genes assigned, score=%.0f",
                         database, dedup_level, granularity, list_name,
                         n_sig_terms, n_themes_hit, pct_genes_assigned, score)) %>%
  pull(label) %>% cat(sep = "\n")

# GO Slim ORA spotlight
cat("\n\n=== GO Slim ORA Results ===\n")
summary_df %>%
  filter(database == "GO Slim", !is.na(n_sig_terms)) %>%
  mutate(label = sprintf("  %s | %s | %s: %d sig terms, %d themes, %.0f%% genes assigned",
                         dedup_level, granularity, list_name,
                         n_sig_terms, n_themes_hit, pct_genes_assigned)) %>%
  pull(label) %>% cat(sep = "\n")

# Training_Old viability
cat("\n\n=== Training_Old Viability ===\n")
to_lists <- c("Training_Old", "Sig Old Up", "Sig Old Down")
summary_df %>%
  filter(list_name %in% to_lists) %>%
  mutate(label = sprintf("  %s | %s | %s | n=%d -> %d sig terms, %.0f%% genes",
                         database, dedup_level, granularity,
                         n_fg, n_sig_terms, pct_genes_assigned)) %>%
  pull(label) %>% cat(sep = "\n")

# Final recommendation
cat("\n\n=== Recommendation ===\n")

# Best GO Slim combo
goslim_best <- summary_df %>%
  filter(database == "GO Slim", !is.na(pct_genes_assigned), pct_genes_assigned > 0) %>%
  mutate(score = pct_genes_assigned * n_themes_hit) %>%
  arrange(desc(score))
if (nrow(goslim_best) > 0) {
  b <- goslim_best[1, ]
  cat(sprintf("Best GO Slim ORA: %s | %s | %s\n", b$dedup_level, b$granularity, b$list_name))
  cat(sprintf("  %d sig terms, %d themes, %.0f%% gene coverage\n",
              b$n_sig_terms, b$n_themes_hit, b$pct_genes_assigned))
}

# Best overall
overall_best <- summary_df %>%
  filter(!is.na(pct_genes_assigned), pct_genes_assigned > 0) %>%
  mutate(score = pct_genes_assigned * n_themes_hit) %>%
  arrange(desc(score))
if (nrow(overall_best) > 0) {
  b <- overall_best[1, ]
  cat(sprintf("\nBest overall: %s | %s | %s | %s\n", b$database, b$dedup_level, b$granularity, b$list_name))
  cat(sprintf("  %d sig terms, %d themes, %.0f%% gene coverage\n",
              b$n_sig_terms, b$n_themes_hit, b$pct_genes_assigned))
}

# TO viability across all databases
to_viable <- summary_df %>%
  filter(training_old_viable == TRUE) %>%
  arrange(desc(pct_genes_assigned))
if (nrow(to_viable) > 0) {
  cat(sprintf("\nTraining_Old viable in %d combos. Best:\n", nrow(to_viable)))
  b <- to_viable[1, ]
  cat(sprintf("  %s | %s | %s | %s: %d terms, %.0f%% genes\n",
              b$database, b$dedup_level, b$granularity, b$list_name,
              b$n_sig_terms, b$pct_genes_assigned))
} else {
  cat("\nNo Training_Old-viable combos found at Pi < 0.05.\n")
}

cat("\nDone.\n")
