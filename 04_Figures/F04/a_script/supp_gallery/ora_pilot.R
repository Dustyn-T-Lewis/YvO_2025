# ORA Pilot Study: Best Enrichment Strategy for Panel F Replacement
# Tests 4 thresholds x 3 collections x 3 granularities across Training contrasts
# Output: ora_pilot_summary.csv with viability metrics per combination

setwd(rprojroot::find_rstudio_root_file())

suppressPackageStartupMessages({
  library(tibble)
  library(readr)
})

source("04_Figures/shared/style.R")
source("04_Figures/shared/pathway_utils.R")
source("04_Figures/shared/go_slim_categories.R")

# Re-attach dplyr/tidyr AFTER go_slim_categories.R loads AnnotationDbi,
# so dplyr::select takes precedence over AnnotationDbi::select
suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
})

dep <- read_csv("03_DEP/c_data/03_combined_results.csv", show_col_types = FALSE)
universe <- dep$gene
set.seed(42)

cat("Building pathway collection...\n")
pw_collection <- build_pathway_collection(
  min_size = 10, max_size = 500,
  include_goslim = FALSE, exclude_variants = TRUE, exclude_wp = TRUE
)
cat(sprintf("  %d pathway sets loaded\n", length(pw_collection)))

# --- Section 2: Define Dimensions ----------------------------------------

CONTRASTS <- c("Training_Young", "Training_Old", "Interaction")

# Threshold functions: dep data.frame, contrast name -> character vector of sig genes
thresholds <- list(
  "Pi < 0.05"  = function(d, ctr) d$gene[d[[paste0("pi_score_", ctr)]] < 0.05],
  "FDR < 0.05" = function(d, ctr) d$gene[d[[paste0("adj.P.Val_", ctr)]] < 0.05],
  "P < 0.05"   = function(d, ctr) d$gene[d[[paste0("P.Value_", ctr)]] < 0.05],
  "P < 0.01"   = function(d, ctr) d$gene[d[[paste0("P.Value_", ctr)]] < 0.01]
)

collections <- c("GO Slim", "Comprehensive", "Comprehensive -> Consolidated")
granularities <- c("per_contrast", "union", "directional")

# --- Section 3: Helper Functions ------------------------------------------

#' Build foreground gene lists given threshold and granularity
#' @return named list of character vectors
build_fg_lists <- function(dep, thresh_fn, granularity) {
  per_ctr <- setNames(
    lapply(CONTRASTS, function(ctr) thresh_fn(dep, ctr)),
    CONTRASTS
  )

  switch(granularity,
    "per_contrast" = per_ctr,
    "union" = list(Union = unique(unlist(per_ctr))),
    "directional" = {
      # Get sig genes per contrast with direction
      sig_TY <- per_ctr[["Training_Young"]]
      sig_TO <- per_ctr[["Training_Old"]]
      sig_Int <- per_ctr[["Interaction"]]

      lfc_TY <- setNames(dep$logFC_Training_Young, dep$gene)
      lfc_TO <- setNames(dep$logFC_Training_Old, dep$gene)

      both <- intersect(sig_TY, sig_TO)
      young_only <- setdiff(sig_TY, c(sig_TO, sig_Int))
      old_only <- setdiff(sig_TO, c(sig_TY, sig_Int))

      lists <- list(
        "Sig Both Up"    = both[lfc_TY[both] > 0],
        "Sig Both Down"  = both[lfc_TY[both] <= 0],
        "Sig Young Up"   = young_only[lfc_TY[young_only] > 0],
        "Sig Young Down" = young_only[lfc_TY[young_only] <= 0],
        "Sig Old Up"     = old_only[lfc_TO[old_only] > 0],
        "Sig Old Down"   = old_only[lfc_TO[old_only] <= 0],
        "Interaction"    = sig_Int
      )
      # Drop empty lists
      lists[vapply(lists, length, integer(1)) > 0]
    }
  )
}


#' Run one enrichment combo and return a 1-row summary tibble
run_one_combo <- function(fg_genes, list_name, universe, collection, pw_collection) {

  n_fg <- length(fg_genes)

  # Skip if too few genes

  if (n_fg < 5) {
    return(tibble(
      list_name = list_name, n_fg = n_fg,
      n_sig_terms = NA_integer_, n_themes_hit = NA_integer_,
      pct_coverage = NA_real_, top3_themes = NA_character_,
      median_odds_ratio = NA_real_
    ))
  }

  if (collection == "GO Slim") {
    # Deterministic assignment, no statistical test
    result <- tryCatch(
      assign_go_slim_consolidated(fg_genes, universe),
      error = function(e) NULL
    )
    if (is.null(result)) {
      return(tibble(
        list_name = list_name, n_fg = n_fg,
        n_sig_terms = NA_integer_, n_themes_hit = NA_integer_,
        pct_coverage = NA_real_, top3_themes = NA_character_,
        median_odds_ratio = NA_real_
      ))
    }
    mapped <- result %>% filter(consolidated != "Other")
    n_themes <- n_distinct(mapped$consolidated)
    pct_cov <- 100 * nrow(mapped) / n_fg
    top3 <- mapped %>% count(consolidated, sort = TRUE) %>%
      head(3) %>% pull(consolidated) %>% paste(collapse = "; ")

    return(tibble(
      list_name = list_name, n_fg = n_fg,
      n_sig_terms = NA_integer_, n_themes_hit = n_themes,
      pct_coverage = round(pct_cov, 1), top3_themes = top3,
      median_odds_ratio = NA_real_
    ))
  }

  if (collection == "Comprehensive") {
    ora_res <- tryCatch(
      run_ora_deduplicated(fg_genes, universe, pw_collection,
                           jaccard_cutoff = 0.5, padj_cutoff = 0.05),
      error = function(e) tibble()
    )
    if (nrow(ora_res) == 0) {
      return(tibble(
        list_name = list_name, n_fg = n_fg,
        n_sig_terms = 0L, n_themes_hit = 0L,
        pct_coverage = 0, top3_themes = NA_character_,
        median_odds_ratio = NA_real_
      ))
    }
    n_sig <- nrow(ora_res)
    # Map to consolidated themes
    ora_res$theme <- classify_pathway_func(ora_res$pathway)
    n_themes <- n_distinct(ora_res$theme[ora_res$theme != "Other"])
    # Coverage: genes in any sig pathway
    sig_pw_genes <- unique(unlist(pw_collection[ora_res$pathway]))
    covered <- length(intersect(fg_genes, sig_pw_genes))
    pct_cov <- 100 * covered / n_fg
    top3 <- ora_res %>% arrange(desc(odds_ratio)) %>%
      head(3) %>% pull(pathway) %>%
      clean_pathway_name(max_chars = 35) %>% paste(collapse = "; ")
    med_or <- median(ora_res$odds_ratio[is.finite(ora_res$odds_ratio)], na.rm = TRUE)

    return(tibble(
      list_name = list_name, n_fg = n_fg,
      n_sig_terms = n_sig, n_themes_hit = n_themes,
      pct_coverage = round(pct_cov, 1), top3_themes = top3,
      median_odds_ratio = round(med_or, 2)
    ))
  }

  if (collection == "Comprehensive -> Consolidated") {
    ora_res <- tryCatch(
      run_ora_deduplicated(fg_genes, universe, pw_collection,
                           jaccard_cutoff = 0.5, padj_cutoff = 0.05),
      error = function(e) tibble()
    )
    if (nrow(ora_res) == 0) {
      return(tibble(
        list_name = list_name, n_fg = n_fg,
        n_sig_terms = 0L, n_themes_hit = 0L,
        pct_coverage = 0, top3_themes = NA_character_,
        median_odds_ratio = NA_real_
      ))
    }
    n_sig <- nrow(ora_res)
    ora_res$theme <- classify_pathway_func(ora_res$pathway)
    themes_hit <- ora_res %>% filter(theme != "Other") %>%
      distinct(theme) %>% nrow()
    # Coverage: genes in any sig pathway
    sig_pw_genes <- unique(unlist(pw_collection[ora_res$pathway]))
    covered <- length(intersect(fg_genes, sig_pw_genes))
    pct_cov <- 100 * covered / n_fg
    top3 <- ora_res %>%
      group_by(theme) %>%
      summarise(best_or = max(odds_ratio[is.finite(odds_ratio)], na.rm = TRUE),
                .groups = "drop") %>%
      filter(theme != "Other") %>%
      arrange(desc(best_or)) %>%
      head(3) %>% pull(theme) %>% paste(collapse = "; ")
    med_or <- median(ora_res$odds_ratio[is.finite(ora_res$odds_ratio)], na.rm = TRUE)

    return(tibble(
      list_name = list_name, n_fg = n_fg,
      n_sig_terms = n_sig, n_themes_hit = themes_hit,
      pct_coverage = round(pct_cov, 1), top3_themes = top3,
      median_odds_ratio = round(med_or, 2)
    ))
  }
}


# --- Section 4: Main Loop ------------------------------------------------

results <- list()
row_idx <- 0L

n_combos <- length(thresholds) * length(collections) * length(granularities)
cat(sprintf("\nRunning %d threshold x collection x granularity combos...\n", n_combos))

for (thresh_name in names(thresholds)) {
  thresh_fn <- thresholds[[thresh_name]]

  for (gran in granularities) {
    fg_lists <- build_fg_lists(dep, thresh_fn, gran)

    for (coll in collections) {
      for (list_nm in names(fg_lists)) {
        fg <- fg_lists[[list_nm]]
        row_idx <- row_idx + 1L

        cat(sprintf("  [%d] %s | %s | %s | %s (n=%d)\n",
                    row_idx, thresh_name, coll, gran, list_nm, length(fg)))

        row <- run_one_combo(fg, list_nm, universe, coll, pw_collection)
        row$threshold  <- thresh_name
        row$collection <- coll
        row$granularity <- gran
        results[[row_idx]] <- row
      }
    }
  }
}

summary_df <- bind_rows(results) %>%
  dplyr::select(threshold, collection, granularity, list_name,
                n_fg, n_sig_terms, n_themes_hit, pct_coverage,
                top3_themes, median_odds_ratio)

# Add training_old_viable flag
to_lists <- c("Training_Old", "Sig Old Up", "Sig Old Down")
summary_df <- summary_df %>%
  mutate(training_old_viable = case_when(
    list_name %in% to_lists & !is.na(n_themes_hit) & n_themes_hit > 0 ~ TRUE,
    list_name %in% to_lists ~ FALSE,
    TRUE ~ NA
  ))


# --- Section 5: Output ---------------------------------------------------

outpath <- "04_Figures/F04/c_data/supp_gallery/ora_pilot_summary.csv"
write_csv(summary_df, outpath)
cat(sprintf("\nSaved: %s (%d rows)\n", outpath, nrow(summary_df)))

# Console summary: best combos ranked by themes * coverage
cat("\n=== Top 10 combos (ranked by n_themes_hit * pct_coverage) ===\n")
summary_df %>%
  filter(!is.na(n_themes_hit), n_themes_hit > 0) %>%
  mutate(score = n_themes_hit * pct_coverage) %>%
  arrange(desc(score)) %>%
  head(10) %>%
  mutate(label = sprintf("  %s | %s | %s | %s: %d themes, %.0f%% cov, score=%.0f",
                         threshold, collection, granularity, list_name,
                         n_themes_hit, pct_coverage, score)) %>%
  pull(label) %>% cat(sep = "\n")

# Training_Old viability
cat("\n\n=== Training_Old viability ===\n")
summary_df %>%
  filter(list_name %in% to_lists) %>%
  mutate(label = sprintf("  %s | %s | %s | n=%d -> %s themes, %s sig terms",
                         threshold, collection, granularity,
                         n_fg,
                         ifelse(is.na(n_themes_hit), "NA", as.character(n_themes_hit)),
                         ifelse(is.na(n_sig_terms), "NA", as.character(n_sig_terms)))) %>%
  pull(label) %>% cat(sep = "\n")

# Recommendation
cat("\n\n=== Recommendation ===\n")
viable <- summary_df %>%
  filter(training_old_viable == TRUE) %>%
  mutate(score = n_themes_hit * pct_coverage) %>%
  arrange(desc(score))
if (nrow(viable) > 0) {
  best <- viable[1, ]
  cat(sprintf("Best Training_Old-viable combo: %s | %s | %s\n",
              best$threshold, best$collection, best$granularity))
  cat(sprintf("  %d themes, %.0f%% coverage, %d sig terms\n",
              best$n_themes_hit, best$pct_coverage,
              ifelse(is.na(best$n_sig_terms), 0L, best$n_sig_terms)))
} else {
  cat("No combo produced viable Training_Old results.\n")
  cat("Consider: P < 0.05 threshold or union granularity.\n")
}

cat("\nDone.\n")
