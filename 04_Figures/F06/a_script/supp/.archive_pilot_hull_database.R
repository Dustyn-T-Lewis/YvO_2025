# pilot_hull_database.R — Test which pathway database gives best hull groupings
# for hub proteins in each key WGCNA module.
#
# For each key module x 6 database configs:
#   run fora() ORA, report n_sig, coverage, top pathways
#
# Output: 04_Figures/F06/c_data/wgcna/pilot_hull_database.csv

setwd(rprojroot::find_rstudio_root_file())
source("04_Figures/shared/style.R")
source("04_Figures/shared/pathway_utils.R")

suppressPackageStartupMessages({
  library(dplyr)
  library(fgsea)
  library(msigdbr)
})

# --- Pathway name cleaner (inline) ------------------------------------------

clean_pw_name <- function(x) {
  x <- gsub("^(HALLMARK_|REACTOME_|KEGG_MEDICUS_|KEGG_|GOBP_|GOSLIM_)", "", x)
  x <- gsub("_", " ", x)
  x <- tools::toTitleCase(tolower(x))
  ifelse(nchar(x) > 50, paste0(substr(x, 1, 47), "..."), x)
}

# --- Disease exclusion pattern -----------------------------------------------

disease_pat <- paste0("DISEASE|CANCER|TUMOR|CARCINOMA|LEUKEMIA|LYMPHOMA|",
                      "MELANOMA|GLIOMA|HEPATITIS|HIV|INFECTION|VIRAL|",
                      "BACTERIAL|PARASIT")

# --- Size filter helper ------------------------------------------------------

filter_pw_size <- function(pw_list, min_size = 5, max_size = 500) {
  sizes <- vapply(pw_list, length, integer(1))
  pw_list[sizes >= min_size & sizes <= max_size]
}

# --- Load data ---------------------------------------------------------------

hub_df <- read.csv("04_Figures/F06/c_data/wgcna/wgcna_hub_proteins.csv",
                    stringsAsFactors = FALSE)
mod_df <- read.csv("04_Figures/F06/c_data/wgcna/wgcna_module_assignments.csv",
                    stringsAsFactors = FALSE)
key_modules <- trimws(readLines("04_Figures/F06/c_data/wgcna/key_modules.txt"))
key_modules <- key_modules[nzchar(key_modules)]

universe <- mod_df$gene[mod_df$module_color != "grey"]
universe <- unique(universe)

cat(sprintf("Universe: %d non-grey proteins\n", length(universe)))
cat(sprintf("Key modules: %s\n", paste(key_modules, collapse = ", ")))

# --- Build 6 database configs ------------------------------------------------

message("\n=== Building pathway databases ===")

# 1. Hallmark
h_raw <- msigdbr(species = "Homo sapiens", collection = "H")
pw_hallmark <- split(h_raw$gene_symbol, h_raw$gs_name)
pw_hallmark <- lapply(pw_hallmark, unique)
pw_hallmark <- filter_pw_size(pw_hallmark)
message(sprintf("Hallmark: %d sets", length(pw_hallmark)))

# 2. GO Slim
pw_goslim <- build_goslim_gene_sets(
  species = "Homo sapiens", min_size = 5, max_size = 500
)
message(sprintf("GO Slim: %d sets", length(pw_goslim)))

# 3. KEGG only (disease-filtered)
kegg_raw <- msigdbr(species = "Homo sapiens", collection = "C2",
                     subcollection = "CP:KEGG_MEDICUS")
kegg_raw <- kegg_raw[!grepl(disease_pat, kegg_raw$gs_name, ignore.case = TRUE), ]
pw_kegg <- split(kegg_raw$gene_symbol, kegg_raw$gs_name)
pw_kegg <- lapply(pw_kegg, unique)
pw_kegg <- filter_pw_size(pw_kegg)
message(sprintf("KEGG: %d sets", length(pw_kegg)))

# 4. Reactome only (disease-filtered)
react_raw <- msigdbr(species = "Homo sapiens", collection = "C2",
                      subcollection = "CP:REACTOME")
react_raw <- react_raw[!grepl(disease_pat, react_raw$gs_name, ignore.case = TRUE), ]
pw_reactome <- split(react_raw$gene_symbol, react_raw$gs_name)
pw_reactome <- lapply(pw_reactome, unique)
pw_reactome <- filter_pw_size(pw_reactome)
message(sprintf("Reactome: %d sets", length(pw_reactome)))

# 5. GO:BP only
gobp_raw <- msigdbr(species = "Homo sapiens", collection = "C5",
                      subcollection = "GO:BP")
pw_gobp <- split(gobp_raw$gene_symbol, gobp_raw$gs_name)
pw_gobp <- lapply(pw_gobp, unique)
pw_gobp <- filter_pw_size(pw_gobp)
message(sprintf("GO:BP: %d sets", length(pw_gobp)))

# 6. Combined (via build_pathway_collection, min_size=5)
pw_combined <- build_pathway_collection(
  species = "Homo sapiens", min_size = 5, max_size = 500,
  include_goslim = TRUE
)

# --- Collect all DB configs --------------------------------------------------

db_configs <- list(
  Hallmark = pw_hallmark,
  `GO Slim` = pw_goslim,
  KEGG     = pw_kegg,
  Reactome = pw_reactome,
  `GO:BP`  = pw_gobp,
  Combined = pw_combined
)

# --- Run ORA per module x database -------------------------------------------

all_results <- list()

for (mod in key_modules) {
  mod_hubs <- hub_df$gene[hub_df$module == mod]
  mod_hubs <- intersect(mod_hubs, universe)
  n_hubs <- length(mod_hubs)

  cat(sprintf("\n%s\n%s module: %d hub proteins [%s]\n%s\n",
              strrep("=", 70), toupper(mod), n_hubs,
              paste(mod_hubs, collapse = ", "), strrep("=", 70)))

  for (db_name in names(db_configs)) {
    pw_list <- db_configs[[db_name]]

    # Run fora
    ora_res <- tryCatch({
      res <- fgsea::fora(
        pathways = pw_list,
        genes    = mod_hubs,
        universe = universe,
        minSize  = 5,
        maxSize  = 500
      )
      as.data.frame(res)
    }, error = function(e) {
      message(sprintf("  [%s] fora error: %s", db_name, e$message))
      data.frame()
    })

    if (nrow(ora_res) == 0) {
      cat(sprintf("  %-10s: 0 tested pathways\n", db_name))
      all_results[[length(all_results) + 1]] <- data.frame(
        module = mod, database = db_name, n_hubs = n_hubs,
        n_tested = 0, n_sig = 0, n_groups_3plus = 0,
        n_assigned_hubs = 0, pct_assigned = 0,
        top1 = NA_character_, top2 = NA_character_,
        top3 = NA_character_, top4 = NA_character_,
        top1_overlap = NA_integer_, top2_overlap = NA_integer_,
        top3_overlap = NA_integer_, top4_overlap = NA_integer_,
        top1_padj = NA_real_, top2_padj = NA_real_,
        top3_padj = NA_real_, top4_padj = NA_real_,
        stringsAsFactors = FALSE
      )
      next
    }

    n_tested <- nrow(ora_res)
    sig <- ora_res[!is.na(ora_res$padj) & ora_res$padj < 0.05, ]
    sig <- sig[order(sig$padj), ]
    n_sig <- nrow(sig)

    # Top 4 by padj
    top4 <- head(sig, 4)

    # Count stats from top 4
    if (nrow(top4) > 0) {
      # overlapGenes is a list column from fora
      overlap_genes <- top4$overlapGenes
      overlap_counts <- top4$overlap

      n_groups_3plus <- sum(overlap_counts >= 3)

      # Unique hubs covered by top 4
      all_overlap_genes <- unique(unlist(overlap_genes))
      n_assigned <- length(intersect(all_overlap_genes, mod_hubs))
      pct_assigned <- round(100 * n_assigned / n_hubs, 1)

      # Clean pathway names
      pw_names <- clean_pw_name(top4$pathway)
    } else {
      overlap_counts <- integer(0)
      n_groups_3plus <- 0
      n_assigned <- 0
      pct_assigned <- 0
      pw_names <- character(0)
    }

    # Pad to 4
    pw_names_4 <- c(pw_names, rep(NA_character_, 4))[1:4]
    overlap_4  <- c(as.integer(overlap_counts), rep(NA_integer_, 4))[1:4]
    padj_4     <- c(top4$padj, rep(NA_real_, 4))[1:4]

    row <- data.frame(
      module = mod, database = db_name, n_hubs = n_hubs,
      n_tested = n_tested, n_sig = n_sig, n_groups_3plus = n_groups_3plus,
      n_assigned_hubs = n_assigned, pct_assigned = pct_assigned,
      top1 = pw_names_4[1], top2 = pw_names_4[2],
      top3 = pw_names_4[3], top4 = pw_names_4[4],
      top1_overlap = overlap_4[1], top2_overlap = overlap_4[2],
      top3_overlap = overlap_4[3], top4_overlap = overlap_4[4],
      top1_padj = padj_4[1], top2_padj = padj_4[2],
      top3_padj = padj_4[3], top4_padj = padj_4[4],
      stringsAsFactors = FALSE
    )
    all_results[[length(all_results) + 1]] <- row

    # Console summary
    cat(sprintf("  %-10s: %3d tested | %2d sig | %d groups>=3 | %2d/%2d hubs (%.0f%%)",
                db_name, n_tested, n_sig, n_groups_3plus,
                n_assigned, n_hubs, pct_assigned))
    if (nrow(top4) > 0) {
      cat(" | top: ")
      for (j in seq_len(nrow(top4))) {
        cat(sprintf("%s(%d)", pw_names_4[j], overlap_4[j]))
        if (j < nrow(top4)) cat(", ")
      }
    }
    cat("\n")
  }
}

# --- Save results ------------------------------------------------------------

results_df <- do.call(rbind, all_results)
out_path <- "04_Figures/F06/c_data/wgcna/pilot_hull_database.csv"
write.csv(results_df, out_path, row.names = FALSE)
cat(sprintf("\nSaved: %s (%d rows)\n", out_path, nrow(results_df)))

# --- Per-module summary ------------------------------------------------------

cat("\n\n")
cat(strrep("=", 80))
cat("\n  SUMMARY: Best database per module (by pct_assigned, then n_groups_3plus)\n")
cat(strrep("=", 80))
cat("\n\n")

for (mod in key_modules) {
  mod_rows <- results_df[results_df$module == mod, ]
  mod_rows <- mod_rows[order(-mod_rows$pct_assigned, -mod_rows$n_groups_3plus,
                              -mod_rows$n_sig), ]

  cat(sprintf("--- %s (%d hubs) ---\n", toupper(mod), mod_rows$n_hubs[1]))

  for (i in seq_len(nrow(mod_rows))) {
    r <- mod_rows[i, ]
    marker <- if (i == 1) " <-- BEST" else ""
    cat(sprintf("  %d. %-10s: %2d sig, %d grp>=3, %2d/%2d hubs (%.0f%%)%s\n",
                i, r$database, r$n_sig, r$n_groups_3plus,
                r$n_assigned_hubs, r$n_hubs, r$pct_assigned, marker))
    # Show top pathways for this config
    tops <- c(r$top1, r$top2, r$top3, r$top4)
    overlaps <- c(r$top1_overlap, r$top2_overlap, r$top3_overlap, r$top4_overlap)
    padjs <- c(r$top1_padj, r$top2_padj, r$top3_padj, r$top4_padj)
    for (j in seq_along(tops)) {
      if (!is.na(tops[j])) {
        cat(sprintf("       %s (n=%d, padj=%.2e)\n", tops[j], overlaps[j], padjs[j]))
      }
    }
  }
  cat("\n")
}

# --- Green module: ETC complex check -----------------------------------------

cat(strrep("=", 80))
cat("\n  GREEN MODULE: Manual ETC assignments vs database ORA\n")
cat(strrep("=", 80))
cat("\n\n")

green_hubs <- hub_df$gene[hub_df$module == "green"]
etc_patterns <- c("NDUF", "UQCR", "COX", "ATP5")
etc_assigned <- green_hubs[grepl(paste(etc_patterns, collapse = "|"), green_hubs)]
cat(sprintf("Manual ETC pattern (NDUF/UQCR/COX/ATP5): %d/%d hubs matched\n",
            length(etc_assigned), length(green_hubs)))
cat(sprintf("  Matched: %s\n", paste(etc_assigned, collapse = ", ")))
cat(sprintf("  Unmatched: %s\n\n",
            paste(setdiff(green_hubs, etc_assigned), collapse = ", ")))

green_rows <- results_df[results_df$module == "green", ]
cat("Database ORA comparison for green:\n")
for (i in seq_len(nrow(green_rows))) {
  r <- green_rows[i, ]
  cat(sprintf("  %-10s: %2d sig, %.0f%% hubs covered",
              r$database, r$n_sig, r$pct_assigned))
  if (!is.na(r$top1)) cat(sprintf(" | best: %s(%d)", r$top1, r$top1_overlap))
  cat("\n")
}
cat("\n")
cat("NOTE: Manual ETC grouping is domain-specific and may outperform generic\n")
cat("database ORA for the green (ETC/OXPHOS) module. Consider keeping manual\n")
cat("for green and using database ORA for other modules.\n")

cat("\nDone.\n")
