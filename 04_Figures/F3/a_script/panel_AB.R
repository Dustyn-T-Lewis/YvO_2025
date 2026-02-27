################################################################################
#   Figure 3 — Panels A & B: Volcano Ring Composites
#   Panel A: Aging Effect (volcano + enrichment arcs)
#   Panel B: Reversal Effect (Aging - Training_Old) (volcano + enrichment arcs)
#
#   Unlike F2 which uses make_volcano_ring_pair() for paired contrasts that
#   share geometry, F3 creates each ring independently since Aging and
#   Reversal are fundamentally different contrasts.
#
#   Generates:
#     b_reports/panel_A_volcano.pdf, panel_A_volcano.png
#     b_reports/panel_B_volcano.pdf, panel_B_volcano.png
#     c_data/panel_A/ring_terms.csv, volcano_aging.csv
#     c_data/panel_B/ring_terms.csv, volcano_reversal.csv
################################################################################

# === 0. SHARED SETUP =========================================================

if (!exists("dep_df")) source("04_Figures/F3/a_script/YvO_F3_setup.R")

# Source the ring utility (defines make_volcano_ring, prepare_ring_data,
# build_panel_legend, clean_ring_label, etc.)
source("04_Figures/F3/a_script/volcano_ring.R")

message("Panels A & B: volcano ring composites...")

# ==============================================================================
# HELPERS — adapted from make_volcano_ring_pair() internals
# ==============================================================================

# --- select_ring_terms(): pick top n_each UP + n_each DOWN pathways,
#     balanced across Hallmark/GO:BP (padj < 0.05)
select_ring_terms <- function(go_df, contrast_name, n_each = 5,
                              databases = c("Hallmark", "GO:BP")) {

  sig <- go_df %>%
    filter(contrast == contrast_name, database %in% databases, padj < 0.05)

  # Interleave databases: aim for ceiling(n/2) Hallmark, rest GO:BP
  pick_direction <- function(sig_df, n_target) {
    pool <- sig_df %>% arrange(padj)
    interleave_db <- function(df, n) {
      hm <- df %>% filter(database == "Hallmark")
      bp <- df %>% filter(database == "GO:BP")
      n_hm <- min(nrow(hm), ceiling(n / 2))
      n_bp <- min(nrow(bp), n - n_hm)
      # Re-balance if one database was short
      n_hm <- min(nrow(hm), n - n_bp)
      bind_rows(hm %>% slice_head(n = n_hm),
                bp %>% slice_head(n = n_bp))
    }
    interleave_db(pool, n_target) %>% slice_head(n = n_target)
  }

  top_terms <- bind_rows(
    pick_direction(sig %>% filter(NES > 0), n_each),
    pick_direction(sig %>% filter(NES < 0), n_each)
  ) %>% slice_head(n = n_each * 2)

  message(sprintf("  select_ring_terms [%s]: %d terms (%d up + %d down)",
                  contrast_name, nrow(top_terms),
                  sum(top_terms$NES > 0), sum(top_terms$NES < 0)))

  if (nrow(top_terms) == 0) {
    stop("No significant terms found for contrast '", contrast_name,
         "' -- cannot build volcano ring")
  }

  top_terms
}

# --- build_ring_with_gaps(): calls prepare_ring_data() then post-processes
#     arc geometry for 8-degree split gaps at UP/DOWN boundary
build_ring_with_gaps <- function(top_terms, contrast_name, go_df,
                                 n_each = 5,
                                 databases = c("Hallmark", "GO:BP")) {

  # Force the selected terms into prepare_ring_data by injecting synthetic
  # padj ordering that preserves the UP-then-DOWN arrangement
  real_rows <- go_df %>%
    filter(contrast == contrast_name, pathway %in% top_terms$pathway)
  go_subset <- real_rows %>%
    mutate(padj = match(pathway, top_terms$pathway) * 1e-10)

  ring <- prepare_ring_data(
    go_df        = go_subset,
    contrast     = contrast_name,
    n_terms      = nrow(top_terms),
    gap_degrees  = 3,
    start_offset = 0,
    databases    = databases
  )

  # Post-process: widen gaps at UP/DOWN hemisphere split points
  n <- nrow(ring)
  if (n >= 2) {
    gap_normal <- 3
    gap_split  <- 8
    start_offset <- 0

    # Gap pattern: split gaps at positions n_each (UP->DOWN) and n (last->first)
    # e.g. for 10 terms: [3, 3, 3, 3, 8, 3, 3, 3, 3, 8]
    gaps <- rep(gap_normal, n)
    gaps[min(n_each, n)] <- gap_split   # between last UP and first DOWN
    gaps[n]              <- gap_split   # between last DOWN and first UP (wraps)

    total_gap <- sum(gaps)
    arc_width_deg <- (360 - total_gap) / n

    # Recalculate start/end for each arc
    cum_offset <- 0
    for (i in seq_len(n)) {
      if (i == 1) {
        cum_offset <- 0
      } else {
        cum_offset <- cum_offset + arc_width_deg + gaps[i - 1]
      }
      ring$start_deg[i] <- start_offset + cum_offset
      ring$end_deg[i]   <- ring$start_deg[i] + arc_width_deg
      ring$mid_deg[i]   <- (ring$start_deg[i] + ring$end_deg[i]) / 2
      ring$start_rad[i] <- ring$start_deg[i] * pi / 180
      ring$end_rad[i]   <- ring$end_deg[i]   * pi / 180
      ring$mid_rad[i]   <- ring$mid_deg[i]   * pi / 180
    }
  }

  ring
}

# ==============================================================================
# PANEL A — Aging Volcano Ring
# ==============================================================================

message("  Building Panel A (Aging)...")

top_terms_aging <- select_ring_terms(fgsea_all, "Aging", n_each = 5)
ring_aging      <- build_ring_with_gaps(top_terms_aging, "Aging", fgsea_all,
                                         n_each = 5)

message("  Ring geometry: Aging ", nrow(ring_aging), " arcs")

pA <- make_volcano_ring(
  de_df              = dep_df,
  go_df              = fgsea_all,
  contrast           = "Aging",
  title              = "Aging Effect",
  ring_data_override = ring_aging
)

# ==============================================================================
# PANEL B — Reversal Volcano Ring (Training_Old - Aging)
# ==============================================================================

message("  Building Panel B (Reversal)...")

top_terms_rev <- select_ring_terms(fgsea_all, "Reversal", n_each = 5)
ring_rev      <- build_ring_with_gaps(top_terms_rev, "Reversal",
                                       fgsea_all, n_each = 5)

message("  Ring geometry: Reversal ", nrow(ring_rev), " arcs")

pB <- make_volcano_ring(
  de_df              = dep_df,
  go_df              = fgsea_all,
  contrast           = "Reversal",
  title              = "Reversal Effect (Training \u2212 Aging)",
  ring_data_override = ring_rev
)

# ==============================================================================
# SAVE — Individual panels with per-panel legend strips
# ==============================================================================

message("  Saving Panel A & B figures...")

# PDF device selection (prefer cairo_pdf for font embedding)
pdf_device <- tryCatch(
  { cairo_pdf(tempfile()); dev.off(); cairo_pdf },
  error = function(e) "pdf"
)

panel_w <- 160   # mm
panel_h <- 180   # mm (extra 20mm for legend strip)

# Build per-panel legends (NES gradient + point legend)
legend_A <- build_panel_legend()
legend_B <- build_panel_legend()

pA_with_key <- pA / legend_A + plot_layout(heights = c(6, 1))
pB_with_key <- pB / legend_B + plot_layout(heights = c(6, 1))

# Panel A
ggsave(file.path(RPT_DIR, "panel_A_volcano.pdf"), pA_with_key,
       width = panel_w, height = panel_h, units = "mm", device = pdf_device)
ggsave(file.path(RPT_DIR, "panel_A_volcano.png"), pA_with_key,
       width = panel_w, height = panel_h, units = "mm", dpi = 300)
message("  Saved: ", file.path(RPT_DIR, "panel_A_volcano.pdf"))
message("  Saved: ", file.path(RPT_DIR, "panel_A_volcano.png"))

# Panel B
ggsave(file.path(RPT_DIR, "panel_B_volcano.pdf"), pB_with_key,
       width = panel_w, height = panel_h, units = "mm", device = pdf_device)
ggsave(file.path(RPT_DIR, "panel_B_volcano.png"), pB_with_key,
       width = panel_w, height = panel_h, units = "mm", dpi = 300)
message("  Saved: ", file.path(RPT_DIR, "panel_B_volcano.pdf"))
message("  Saved: ", file.path(RPT_DIR, "panel_B_volcano.png"))

# ==============================================================================
# EXPORT — Ring term CSVs + flat volcano CSVs
# ==============================================================================

message("  Exporting data CSVs...")

# Ring term CSVs (remove list-column gene_list for flat export)
ring_aging_export <- ring_aging %>% dplyr::select(-gene_list)
ring_rev_export <- ring_rev %>% dplyr::select(-gene_list)

write_csv(ring_aging_export, file.path(DAT_DIR, "panel_A", "ring_terms.csv"))
write_csv(ring_rev_export, file.path(DAT_DIR, "panel_B", "ring_terms.csv"))
message("  Exported: ", file.path(DAT_DIR, "panel_A", "ring_terms.csv"))
message("  Exported: ", file.path(DAT_DIR, "panel_B", "ring_terms.csv"))

# Flat volcano CSVs (gene, logFC, p-value, pi_score, direction)
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

export_volcano_csv("Aging",        "panel_A", "volcano_aging.csv")
export_volcano_csv("Reversal", "panel_B", "volcano_reversal.csv")
message("  Exported: ", file.path(DAT_DIR, "panel_A", "volcano_aging.csv"))
message("  Exported: ", file.path(DAT_DIR, "panel_B", "volcano_reversal.csv"))

message("  Panels A & B complete")
