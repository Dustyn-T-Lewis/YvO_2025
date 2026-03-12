################################################################################
#   Figure 2 — Panels A, B & C: Volcano Ring Composites
#   Generates: panel_A_volcano.pdf/png, panel_B_volcano.pdf/png,
#              panel_C_volcano.pdf/png
#              + c_data/panel_A/volcano_young.csv
#              + c_data/panel_B/volcano_old.csv
#              + c_data/panel_C/volcano_interaction.csv
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
#    - Three contrasts (Young, Old, Interaction) are shown as separate panels,
#      not jointly tested. No family-wise correction needed; each panel is
#      self-contained.
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
# PANELS A, B & C — Volcano Ring Composites (volcano + enrichment arcs)
# ==============================================================================

message("Panels A, B & C: volcano ring composites...")

# Source the ring utility (defines make_volcano_ring, build_panel_legend, etc.)
source("04_Figures/shared/volcano_ring.R")

# ── Direction-balanced ring term selection ────────────────────────────────────

select_ring_terms <- function(go_df, contrast_name, n_each = 6,
                              databases = c("Hallmark", "GO:BP")) {
  sig <- go_df %>%
    filter(contrast == contrast_name, database %in% databases, padj < 0.05)

  pick_direction <- function(sig_df, n_target) {
    pool <- sig_df %>% arrange(padj)
    hm <- pool %>% filter(database == "Hallmark")
    bp <- pool %>% filter(database == "GO:BP")
    n_hm <- min(nrow(hm), ceiling(n_target / 2))
    n_bp <- min(nrow(bp), n_target - n_hm)
    n_hm <- min(nrow(hm), n_target - n_bp)
    bind_rows(hm %>% slice_head(n = n_hm),
              bp %>% slice_head(n = n_bp)) %>%
      slice_head(n = n_target)
  }

  bind_rows(
    pick_direction(sig %>% filter(NES > 0), n_each),
    pick_direction(sig %>% filter(NES < 0), n_each)
  ) %>% slice_head(n = n_each * 2)
}

build_ring_with_gaps <- function(top_terms, contrast_name, go_df,
                                 n_each = 6,
                                 databases = c("Hallmark", "GO:BP")) {
  real_rows <- go_df %>%
    filter(contrast == contrast_name, pathway %in% top_terms$pathway)
  go_subset <- real_rows %>%
    mutate(padj = match(pathway, top_terms$pathway) * 1e-10)

  ring <- prepare_ring_data(
    go_df = go_subset, contrast = contrast_name,
    n_terms = nrow(top_terms), gap_degrees = 3, start_offset = 0,
    databases = databases
  )

  n <- nrow(ring)
  if (n >= 2) {
    gap_normal <- 3; gap_split <- 8
    gaps <- rep(gap_normal, n)
    gaps[min(n_each, n)] <- gap_split
    gaps[n]              <- gap_split
    arc_width_deg <- (360 - sum(gaps)) / n
    cum_offset <- 0
    for (i in seq_len(n)) {
      if (i > 1) cum_offset <- cum_offset + arc_width_deg + gaps[i - 1]
      ring$start_deg[i] <- cum_offset
      ring$end_deg[i]   <- ring$start_deg[i] + arc_width_deg
      ring$mid_deg[i]   <- (ring$start_deg[i] + ring$end_deg[i]) / 2
      ring$start_rad[i] <- ring$start_deg[i] * pi / 180
      ring$end_rad[i]   <- ring$end_deg[i]   * pi / 180
      ring$mid_rad[i]   <- ring$mid_deg[i]   * pi / 180
    }
  }
  ring
}

# ── Build three separate volcano rings ────────────────────────────────────────

# Panel A: Training Effect (Young)
top_terms_A <- select_ring_terms(fgsea_all, "Training_Young", n_each = 6)
ring_A      <- build_ring_with_gaps(top_terms_A, "Training_Young", fgsea_all, n_each = 6)

pA <- make_volcano_ring(
  de_df     = dep_df,
  go_df     = fgsea_all,
  contrast  = "Training_Young",
  contrast_title    = "Training Effect (Young)",
  contrast_subtitle = "Post_Young - Pre_Young",
  title_size    = 16,
  label_size    = TXT_PATHWAY,
  ring_data_override = ring_A
)

# Panel B: Training Effect (Old)
top_terms_B <- select_ring_terms(fgsea_all, "Training_Old", n_each = 6)
ring_B      <- build_ring_with_gaps(top_terms_B, "Training_Old", fgsea_all, n_each = 6)

pB <- make_volcano_ring(
  de_df     = dep_df,
  go_df     = fgsea_all,
  contrast  = "Training_Old",
  contrast_title    = "Training Effect (Old)",
  contrast_subtitle = "Post_Old - Pre_Old",
  title_size    = 16,
  label_size    = TXT_PATHWAY,
  ring_data_override = ring_B
)

# Panel C (volcano): Interaction Effect
top_terms_C <- select_ring_terms(fgsea_all, "Interaction", n_each = 6)
ring_C      <- build_ring_with_gaps(top_terms_C, "Interaction", fgsea_all, n_each = 6)

pC_volcano <- make_volcano_ring(
  de_df     = dep_df,
  go_df     = fgsea_all,
  contrast  = "Interaction",
  contrast_title    = "Interaction Effect",
  contrast_subtitle = "(Old - Young) x (Post - Pre)",
  title_size    = 16,
  label_size    = TXT_PATHWAY,
  ring_data_override = ring_C
)

# ── Save panel figures (160×160mm standalone) ─────────────────────────────────

panel_w <- 160
panel_h <- 160

# Use cairo_pdf if available
pdf_device <- tryCatch(
  { cairo_pdf(tempfile()); dev.off(); cairo_pdf },
  error = function(e) "pdf"
)

ggsave(file.path(RPT_DIR, "panel_A_volcano.pdf"), pA,
       width = panel_w, height = panel_h, units = "mm", device = pdf_device)
ggsave(file.path(RPT_DIR, "panel_A_volcano.png"), pA,
       width = panel_w, height = panel_h, units = "mm", dpi = 300)

ggsave(file.path(RPT_DIR, "panel_B_volcano.pdf"), pB,
       width = panel_w, height = panel_h, units = "mm", device = pdf_device)
ggsave(file.path(RPT_DIR, "panel_B_volcano.png"), pB,
       width = panel_w, height = panel_h, units = "mm", dpi = 300)

ggsave(file.path(RPT_DIR, "panel_C_volcano.pdf"), pC_volcano,
       width = panel_w, height = panel_h, units = "mm", device = pdf_device)
ggsave(file.path(RPT_DIR, "panel_C_volcano.png"), pC_volcano,
       width = panel_w, height = panel_h, units = "mm", dpi = 300)

# ── Export ring term CSVs ─────────────────────────────────────────────────────

ring_data_A <- attr(pA, "ring_data")
if (!is.null(ring_data_A) && nrow(ring_data_A) > 0) {
  write_csv(ring_data_A %>% dplyr::select(-gene_list),
            file.path(DAT_DIR, "panel_A", "ring_terms.csv"))
}

ring_data_B <- attr(pB, "ring_data")
if (!is.null(ring_data_B) && nrow(ring_data_B) > 0) {
  write_csv(ring_data_B %>% dplyr::select(-gene_list),
            file.path(DAT_DIR, "panel_B", "ring_terms.csv"))
}

ring_data_C <- attr(pC_volcano, "ring_data")
if (!is.null(ring_data_C) && nrow(ring_data_C) > 0) {
  dir.create(file.path(DAT_DIR, "panel_C"), recursive = TRUE, showWarnings = FALSE)
  write_csv(ring_data_C %>% dplyr::select(-gene_list),
            file.path(DAT_DIR, "panel_C", "ring_terms.csv"))
}

# ── Export flat volcano CSVs ──────────────────────────────────────────────────

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
export_volcano_csv("Interaction",    "panel_C", "volcano_interaction.csv")

message("  Panels A, B & C (volcano rings) saved")
