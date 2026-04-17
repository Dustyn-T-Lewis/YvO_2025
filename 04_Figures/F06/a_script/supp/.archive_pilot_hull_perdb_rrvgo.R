# pilot_hull_perdb_rrvgo.R — Compare "Combined single-BH" vs "PerDB BH + rrvgo"
# for hub protein ORA hull groupings across key WGCNA modules.
#
# Approach:
#   1. run_ora_deduplicated() — per-database BH + within-database Jaccard dedup
#   2. rrvgo semantic dedup on GO:BP survivors (Wang similarity, threshold=0.7)
#   3. Combine rrvgo-reduced GO:BP with non-GO:BP survivors, take top 4
#   4. Compare against "Combined single-BH" baseline from pilot_hull_database.csv
#
# Output: 04_Figures/F06/c_data/wgcna/pilot_hull_perdb_rrvgo.csv

setwd(rprojroot::find_rstudio_root_file())
source("04_Figures/shared/style.R")
source("04_Figures/shared/pathway_utils.R")

suppressPackageStartupMessages({
  library(dplyr)
  library(fgsea)
  library(msigdbr)
  library(rrvgo)
})

# --- Pathway name cleaner (inline, matches pilot_hull_database.R) -----------

clean_pw_name <- function(x) {
  x <- gsub("^(HALLMARK_|REACTOME_|KEGG_MEDICUS_|KEGG_|GOBP_|GOSLIM_)", "", x)
  x <- gsub("_", " ", x)
  x <- tools::toTitleCase(tolower(x))
  ifelse(nchar(x) > 50, paste0(substr(x, 1, 47), "..."), x)
}

# --- Load data ---------------------------------------------------------------

hub_df <- read.csv("04_Figures/F06/c_data/wgcna/wgcna_hub_proteins.csv",
                    stringsAsFactors = FALSE)
mod_df <- read.csv("04_Figures/F06/c_data/wgcna/wgcna_module_assignments.csv",
                    stringsAsFactors = FALSE)
key_modules <- trimws(readLines("04_Figures/F06/c_data/wgcna/key_modules.txt"))
key_modules <- key_modules[nzchar(key_modules)]

universe <- unique(mod_df$gene[mod_df$module_color != "grey"])

cat(sprintf("Universe: %d non-grey proteins\n", length(universe)))
cat(sprintf("Key modules: %s\n", paste(key_modules, collapse = ", ")))

# --- Load previous "Combined single-BH" baseline ----------------------------

baseline_df <- read.csv("04_Figures/F06/c_data/wgcna/pilot_hull_database.csv",
                         stringsAsFactors = FALSE)
baseline_combined <- baseline_df[baseline_df$database == "Combined", ]

# --- Build GO:BP gs_name -> GO ID mapping ------------------------------------

message("\n=== Building GO:BP name-to-ID mapping ===")
gobp_raw <- msigdbr(species = "Homo sapiens", collection = "C5",
                     subcollection = "GO:BP")
gobp_name_to_id <- unique(gobp_raw[, c("gs_name", "gs_exact_source")])
gobp_map <- setNames(gobp_name_to_id$gs_exact_source, gobp_name_to_id$gs_name)
message(sprintf("GO:BP mapping: %d terms", length(gobp_map)))

# --- Build pathway collection (min_size=5 for hub ORA) ----------------------

message("\n=== Building pathway collection ===")
pw_collection <- build_pathway_collection(
  species = "Homo sapiens", min_size = 5, max_size = 500,
  include_goslim = TRUE
)

# --- Run per-module analysis ------------------------------------------------

all_results <- list()

for (mod in key_modules) {
  mod_hubs <- hub_df$gene[hub_df$module == mod]
  mod_hubs <- intersect(mod_hubs, universe)
  n_hubs <- length(mod_hubs)

  cat(sprintf("\n%s\n%s module: %d hub proteins [%s]\n%s\n",
              strrep("=", 70), toupper(mod), n_hubs,
              paste(mod_hubs, collapse = ", "), strrep("=", 70)))

  # --- Step 1: run_ora_deduplicated (per-database BH + Jaccard dedup) --------

  message("  Step 1: run_ora_deduplicated (per-DB BH + Jaccard 0.5)")
  ora_dedup <- tryCatch({
    run_ora_deduplicated(
      genes = mod_hubs,
      universe = universe,
      pathways = pw_collection,
      jaccard_cutoff = 0.5,
      min_size = 5,
      max_size = 500,
      padj_cutoff = 0.05
    )
  }, error = function(e) {
    message(sprintf("  ORA error: %s", e$message))
    tibble::tibble()
  })

  n_survivors_pre_rrvgo <- nrow(ora_dedup)
  cat(sprintf("  Survivors after per-DB BH + Jaccard: %d\n", n_survivors_pre_rrvgo))

  if (n_survivors_pre_rrvgo == 0) {
    # No sig results — record empty row
    all_results[[length(all_results) + 1]] <- data.frame(
      module = mod, method = "PerDB_BH_rrvgo", n_hubs = n_hubs,
      n_survivors_pre = 0, n_survivors_post = 0,
      n_groups_3plus = 0, n_assigned_hubs = 0, pct_assigned = 0,
      top1 = NA_character_, top2 = NA_character_,
      top3 = NA_character_, top4 = NA_character_,
      top1_overlap = NA_integer_, top2_overlap = NA_integer_,
      top3_overlap = NA_integer_, top4_overlap = NA_integer_,
      top1_padj = NA_real_, top2_padj = NA_real_,
      top3_padj = NA_real_, top4_padj = NA_real_,
      top1_db = NA_character_, top2_db = NA_character_,
      top3_db = NA_character_, top4_db = NA_character_,
      stringsAsFactors = FALSE
    )
    next
  }

  # --- Step 2: Split into GO:BP vs non-GO:BP survivors -----------------------

  is_gobp <- ora_dedup$database == "GO:BP"
  gobp_survivors <- ora_dedup[is_gobp, ]
  non_gobp_survivors <- ora_dedup[!is_gobp, ]

  cat(sprintf("  GO:BP survivors: %d | Non-GO:BP survivors: %d\n",
              nrow(gobp_survivors), nrow(non_gobp_survivors)))

  # --- Step 3: rrvgo on GO:BP survivors (if >= 2) ----------------------------

  if (nrow(gobp_survivors) >= 2) {
    message("  Step 3: rrvgo semantic reduction on GO:BP survivors")

    # Map pathway names to GO IDs
    go_ids <- gobp_map[gobp_survivors$pathway]
    valid_mask <- !is.na(go_ids)

    if (sum(valid_mask) >= 2) {
      go_ids_valid <- go_ids[valid_mask]
      padj_valid <- gobp_survivors$padj[valid_mask]

      # Calculate similarity matrix (suppress warnings about missing GO IDs)
      sim_matrix <- tryCatch({
        suppressWarnings(
          rrvgo::calculateSimMatrix(
            go_ids_valid,
            orgdb = "org.Hs.eg.db",
            ont = "BP",
            method = "Wang"
          )
        )
      }, error = function(e) {
        message(sprintf("  rrvgo simMatrix error: %s", e$message))
        NULL
      })

      if (!is.null(sim_matrix) && nrow(sim_matrix) >= 2) {
        # Filter to GO IDs that survived simMatrix construction
        sim_go_ids <- rownames(sim_matrix)
        sim_mask <- go_ids_valid %in% sim_go_ids
        scores <- setNames(-log10(padj_valid[sim_mask]),
                           go_ids_valid[sim_mask])

        reduced <- rrvgo::reduceSimMatrix(
          sim_matrix,
          scores = scores,
          threshold = 0.7,
          orgdb = "org.Hs.eg.db"
        )

        # Keep only representative terms (where go == parent)
        representatives <- reduced[reduced$go == reduced$parent, ]
        rep_go_ids <- representatives$go
        cat(sprintf("  rrvgo: %d GO:BP -> %d representatives (threshold=0.7)\n",
                    length(sim_go_ids), length(rep_go_ids)))

        # Map representative GO IDs back to pathway names
        # Build reverse map: GO ID -> gs_name (using the survivors)
        id_to_name <- setNames(gobp_survivors$pathway[valid_mask],
                               go_ids[valid_mask])
        # Filter to survivors whose GO ID is a representative
        rep_names <- id_to_name[names(id_to_name) %in% rep_go_ids]
        gobp_reduced <- gobp_survivors[gobp_survivors$pathway %in% rep_names, ]
      } else {
        message("  rrvgo: simMatrix too small, keeping all GO:BP survivors")
        gobp_reduced <- gobp_survivors
      }
    } else {
      message("  rrvgo: too few valid GO IDs, keeping all GO:BP survivors")
      gobp_reduced <- gobp_survivors
    }
  } else {
    message("  Skipping rrvgo (0-1 GO:BP survivors)")
    gobp_reduced <- gobp_survivors
  }

  cat(sprintf("  GO:BP after rrvgo: %d\n", nrow(gobp_reduced)))

  # --- Step 4: Combine and take top 4 by padj --------------------------------

  final_survivors <- rbind(gobp_reduced, non_gobp_survivors)
  final_survivors <- final_survivors[order(final_survivors$padj), ]
  n_survivors_post <- nrow(final_survivors)

  top4 <- head(final_survivors, 4)

  if (nrow(top4) > 0) {
    overlap_genes_list <- top4$overlapGenes
    overlap_counts <- top4$overlap

    n_groups_3plus <- sum(overlap_counts >= 3)
    all_overlap_genes <- unique(unlist(overlap_genes_list))
    n_assigned <- length(intersect(all_overlap_genes, mod_hubs))
    pct_assigned <- round(100 * n_assigned / n_hubs, 1)

    pw_names <- clean_pw_name(top4$pathway)
    db_names <- top4$database
  } else {
    overlap_counts <- integer(0)
    n_groups_3plus <- 0
    n_assigned <- 0
    pct_assigned <- 0
    pw_names <- character(0)
    db_names <- character(0)
  }

  # Pad to 4
  pw_names_4 <- c(pw_names, rep(NA_character_, 4))[1:4]
  overlap_4  <- c(as.integer(overlap_counts), rep(NA_integer_, 4))[1:4]
  padj_4     <- c(top4$padj, rep(NA_real_, 4))[1:4]
  db_4       <- c(db_names, rep(NA_character_, 4))[1:4]

  row <- data.frame(
    module = mod, method = "PerDB_BH_rrvgo", n_hubs = n_hubs,
    n_survivors_pre = n_survivors_pre_rrvgo,
    n_survivors_post = n_survivors_post,
    n_groups_3plus = n_groups_3plus,
    n_assigned_hubs = n_assigned, pct_assigned = pct_assigned,
    top1 = pw_names_4[1], top2 = pw_names_4[2],
    top3 = pw_names_4[3], top4 = pw_names_4[4],
    top1_overlap = overlap_4[1], top2_overlap = overlap_4[2],
    top3_overlap = overlap_4[3], top4_overlap = overlap_4[4],
    top1_padj = padj_4[1], top2_padj = padj_4[2],
    top3_padj = padj_4[3], top4_padj = padj_4[4],
    top1_db = db_4[1], top2_db = db_4[2],
    top3_db = db_4[3], top4_db = db_4[4],
    stringsAsFactors = FALSE
  )
  all_results[[length(all_results) + 1]] <- row
}

# --- Save results ------------------------------------------------------------

results_df <- do.call(rbind, all_results)
out_path <- "04_Figures/F06/c_data/wgcna/pilot_hull_perdb_rrvgo.csv"
write.csv(results_df, out_path, row.names = FALSE)
cat(sprintf("\nSaved: %s (%d rows)\n", out_path, nrow(results_df)))

# --- Comparison table --------------------------------------------------------

cat("\n")
cat(strrep("=", 90))
cat("\n  COMPARISON: Combined_single_BH vs PerDB_BH_rrvgo (top 4 hull groups)\n")
cat(strrep("=", 90))
cat("\n\n")

for (mod in key_modules) {
  # Baseline: Combined single-BH from previous pilot
  bl <- baseline_combined[baseline_combined$module == mod, ]
  # New: PerDB BH + rrvgo
  nw <- results_df[results_df$module == mod, ]

  n_hubs <- if (nrow(bl) > 0) bl$n_hubs[1] else nw$n_hubs[1]

  cat(sprintf("Module: %s (%d hubs)\n", mod, n_hubs))
  cat(sprintf("  %-22s  Survivors  Groups>=3  Coverage  Top pathways\n", "Method"))
  cat(sprintf("  %s\n", strrep("-", 85)))

  # Baseline row
  if (nrow(bl) > 0) {
    bl_tops <- c(bl$top1, bl$top2, bl$top3, bl$top4)
    bl_tops <- bl_tops[!is.na(bl_tops)]
    bl_tops_str <- if (length(bl_tops) > 0) paste(bl_tops, collapse = ", ") else "(none)"
    cat(sprintf("  %-22s  %3d->4     %d         %3.0f%%      %s\n",
                "Combined_single_BH",
                bl$n_sig, bl$n_groups_3plus, bl$pct_assigned,
                bl_tops_str))
  } else {
    cat(sprintf("  %-22s  (no baseline data)\n", "Combined_single_BH"))
  }

  # New row
  if (nrow(nw) > 0) {
    nw_tops <- c(nw$top1, nw$top2, nw$top3, nw$top4)
    nw_dbs  <- c(nw$top1_db, nw$top2_db, nw$top3_db, nw$top4_db)
    nw_display <- character(0)
    for (k in seq_along(nw_tops)) {
      if (!is.na(nw_tops[k])) {
        nw_display <- c(nw_display, sprintf("%s [%s]", nw_tops[k], nw_dbs[k]))
      }
    }
    nw_tops_str <- if (length(nw_display) > 0) paste(nw_display, collapse = ", ") else "(none)"
    cat(sprintf("  %-22s  %3d->%d     %d         %3.0f%%      %s\n",
                "PerDB_BH_rrvgo",
                nw$n_survivors_pre, nw$n_survivors_post,
                nw$n_groups_3plus, nw$pct_assigned,
                nw_tops_str))
  } else {
    cat(sprintf("  %-22s  (no results)\n", "PerDB_BH_rrvgo"))
  }

  cat("\n")
}

cat("Done.\n")
