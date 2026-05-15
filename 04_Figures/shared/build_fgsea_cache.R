#!/usr/bin/env Rscript
# Build fgsea_tstat_all_v2.csv from current Stage 03 t-statistics.
# Source from any figure that needs the cache — skipped if cache is newer
# than 03_DEP/c_data/03_combined_results.csv.

setwd(rprojroot::find_rstudio_root_file())

CACHE_PATH <- "04_Figures/shared/fgsea_tstat_all_v2.csv"
STAGE3_CSV <- "03_DEP/c_data/03_combined_results.csv"

stopifnot("Stage 03 combined_results.csv missing — run 03_DEP/a_script/01_run_dep.R" =
  file.exists(STAGE3_CSV))

if (file.exists(CACHE_PATH) &&
    file.info(CACHE_PATH)$mtime > file.info(STAGE3_CSV)$mtime) {
  message(sprintf("fGSEA cache up to date (%s newer than Stage 03 results) — skipping rebuild",
                  basename(CACHE_PATH)))
} else {
  message("Rebuilding fGSEA cache from current Stage 03 t-statistics...")

  source("04_Figures/shared/pathway_utils.R")
  suppressPackageStartupMessages({
    library(dplyr); library(readr); library(tibble)
  })

  CTRS <- c("Aging", "Training_Young", "Training_Old", "Interaction")
  dep <- read_csv(STAGE3_CSV, show_col_types = FALSE)

  pw_list <- build_pathway_collection(
    min_size = 10, max_size = 500,
    include_goslim   = TRUE,
    exclude_variants = FALSE)

  set.seed(42)
  all_results <- list()
  for (ctr in CTRS) {
    message(sprintf("--- %s ---", ctr))
    ranks <- setNames(dep[[paste0("t_", ctr)]], dep$gene)
    ranks <- sort(ranks[!is.na(ranks)], decreasing = TRUE)

    res <- run_fgsea_deduplicated(
      ranks          = ranks,
      pathways       = pw_list,
      jaccard_cutoff = 0.5,
      nperm          = 10000,
      min_size       = 15,
      max_size       = 500)
    res$contrast <- ctr
    all_results[[ctr]] <- res
  }

  cache_df <- bind_rows(all_results) |>
    mutate(leadingEdge = vapply(leadingEdge, paste, character(1), collapse = ";"))
  cache_df <- cache_df[, c("pathway", "pval", "padj", "log2err", "ES", "NES",
                            "size", "leadingEdge", "database", "contrast")]

  write_csv(cache_df, CACHE_PATH)
  message(sprintf("Wrote %d rows to %s", nrow(cache_df), CACHE_PATH))
}
