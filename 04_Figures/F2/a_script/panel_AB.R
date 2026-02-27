################################################################################
#   Figure 2 — Panels A & B: Volcano Ring Composites
#   Generates: panel_A_volcano.pdf, panel_B_volcano.pdf
#              + c_data/panel_A/volcano_young.csv, c_data/panel_B/volcano_old.csv
################################################################################

if (!exists("dep_df")) source("04_Figures/F2/a_script/YvO_F2_setup.R")

# ==============================================================================
# PANELS A & B — Volcano Ring Composites (volcano + enrichment arcs)
# ==============================================================================

message("Panels A & B: volcano ring composites...")

# Source the ring utility (defines make_volcano_ring_pair)
source("04_Figures/F2/a_script/volcano_ring.R")

# Build the paired ring plots (saves PDFs + ring term CSVs)
pAB_pair <- make_volcano_ring_pair(
  de_df          = dep_df,
  go_df          = fgsea_all,
  contrast_young = "Training_Young",
  contrast_old   = "Training_Old",
  title_young    = "Training Effect (Young)",
  title_old      = "Training Effect (Old)",
  output_dir     = "04_Figures/F2",
  save_outputs   = TRUE
)
pA <- attr(pAB_pair, "p_young")
pB <- attr(pAB_pair, "p_old")

# Also export flat volcano CSVs for tool-ready data access
export_volcano_csv <- function(ctr, panel_dir, filename) {
  col_logFC <- paste0("logFC_", ctr)
  col_pval  <- paste0("P.Value_", ctr)
  col_pi    <- paste0("pi_score_", ctr)
  col_adjp  <- paste0("adj.P.Val_", ctr)

  dep_df %>%
    transmute(
      gene,
      log2_fold_change = round(.data[[col_logFC]], 4),
      neg_log10_pvalue = round(-log10(.data[[col_pval]]), 4),
      pi_score         = round(.data[[col_pi]], 6),
      adjusted_pvalue  = round(.data[[col_adjp]], 6),
      direction = case_when(
        .data[[col_pi]] < 0.05 & .data[[col_logFC]] > 0 ~ "Up",
        .data[[col_pi]] < 0.05 & .data[[col_logFC]] < 0 ~ "Down",
        TRUE ~ "NS"
      )
    ) %>%
    filter(!is.na(log2_fold_change), !is.na(neg_log10_pvalue)) %>%
    arrange(pi_score) %>%
    write_csv(file.path(DAT_DIR, panel_dir, filename))
}

export_volcano_csv("Training_Young", "panel_A", "volcano_young.csv")
export_volcano_csv("Training_Old",   "panel_B", "volcano_old.csv")
message("  Panels A & B saved")
