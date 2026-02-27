################################################################################
#   Figure 2 — Panels A & B: Volcano Ring Composites
#   Generates: panel_A_volcano.pdf, panel_B_volcano.pdf
#              + c_data/panel_A/volcano_young.csv, c_data/panel_B/volcano_old.csv
################################################################################
#
# ── STAT AUDIT (Task 13, 2026-02-27) ─────────────────────────────────────────
# 1. Pi-score threshold (0.05):
#    - Pi-score = -log10(p) * |logFC|, threshold 0.05 applied as a combined
#      significance/effect-size filter. This is a _heuristic_ cutoff; there is
#      no formal Type-I error control on the pi-score itself. The threshold is
#      the same as used in F1 and is justified by the supplementary S1 pi-score
#      distribution analysis (S1_7_pi_score_distributions.pdf). Sensitivity to
#      the 0.05 cutoff was assessed in F1's supplementary; not re-tested here
#      but acknowledged as a potential source of variation.
# 2. fGSEA on ring arcs:
#    - padj (BH correction) is applied within each contrast independently. This
#      is standard for GSEA analyses where each contrast addresses a separate
#      biological question. No cross-contrast correction is applied (documented
#      in YvO_F2_setup.R line 111-113). Top ring terms are filtered at padj<0.05.
# 3. Multiple comparisons:
#    - Two contrasts (Young, Old) are shown as separate panels, not jointly
#      tested. No family-wise correction needed; each panel is self-contained.
# 4. Effect sizes:
#    - NES is the native effect-size metric for GSEA (displayed as arc fill).
#      logFC is the native effect size for protein-level volcano (pi-score
#      integrates logFC with p-value). Both are displayed.
# 5. Sample size: N=15/16 per group; sufficient for limma moderated t-tests
#    with variance shrinkage (Smyth 2004).
# 6. CIs: Not applicable for volcano visualization (CIs on individual protein
#    estimates are available in the underlying limma output but would clutter
#    the display).
# 7. Reproducibility: set.seed(42) in setup; fGSEA uses fgseaMultilevel with
#    deterministic ranking.
# ─────────────────────────────────────────────────────────────────────────────

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
