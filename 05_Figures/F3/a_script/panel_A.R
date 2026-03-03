################################################################################
#   Figure 3 — Panel A: Aging Volcano Ring (wide)
#   Generates: panel_A_volcano.pdf/png, c_data/panel_A/volcano_aging.csv
################################################################################

if (!exists("dep_df")) source("05_Figures/F3/a_script/YvO_F3_setup.R")
source("05_Figures/shared/volcano_ring.R")

message("Panel A: Aging volcano ring...")

# ── Ring term helpers (F3-specific: 8-deg split gaps at UP/DOWN boundary) ─────

select_ring_terms <- function(go_df, contrast_name, n_each = 5,
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

  top_terms <- bind_rows(
    pick_direction(sig %>% filter(NES > 0), n_each),
    pick_direction(sig %>% filter(NES < 0), n_each)
  ) %>% slice_head(n = n_each * 2)

  if (nrow(top_terms) == 0)
    stop("No significant terms for contrast '", contrast_name, "'")
  top_terms
}

build_ring_with_gaps <- function(top_terms, contrast_name, go_df,
                                 n_each = 5,
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

# ── Build volcano ring ───────────────────────────────────────────────────────

top_terms_aging <- select_ring_terms(fgsea_all, "Aging", n_each = 5)
ring_aging      <- build_ring_with_gaps(top_terms_aging, "Aging", fgsea_all, n_each = 5)

pA <- make_volcano_ring(
  de_df = dep_df, go_df = fgsea_all, contrast = "Aging",
  title = "Aging Effect", ring_data_override = ring_aging
) + labs(tag = "A")

# ── Save panel figures (wide: 200 x 200 mm) ─────────────────────────────────

legend_A <- build_panel_legend()
pA_with_key <- pA / legend_A + plot_layout(heights = c(6, 1))

ggsave(file.path(RPT_DIR, "panel_A_volcano.pdf"), pA_with_key,
       width = 200, height = 200, units = "mm", device = cairo_pdf)
ggsave(file.path(RPT_DIR, "panel_A_volcano.png"), pA_with_key,
       width = 200, height = 200, units = "mm", dpi = 300)

# ── Export data CSVs ─────────────────────────────────────────────────────────

ring_data_A <- attr(pA, "ring_data")
if (!is.null(ring_data_A) && nrow(ring_data_A) > 0) {
  write_csv(ring_data_A %>% dplyr::select(-gene_list),
            file.path(DAT_DIR, "panel_A", "ring_terms.csv"))
}

dep_df %>%
  transmute(
    gene,
    log2_fold_change = round(logFC_Aging, 4),
    neg_log10_pvalue = round(-log10(P.Value_Aging), 4),
    pi_score         = round(pi_score_Aging, 6),
    adjusted_pvalue  = round(adj.P.Val_Aging, 6),
    direction = case_when(
      pi_score_Aging < 0.05 & logFC_Aging > 0 ~ "Up",
      pi_score_Aging < 0.05 & logFC_Aging < 0 ~ "Down",
      TRUE ~ "NS"
    )
  ) %>%
  filter(!is.na(log2_fold_change), !is.na(neg_log10_pvalue)) %>%
  arrange(pi_score) %>%
  write_csv(file.path(DAT_DIR, "panel_A", "volcano_aging.csv"))

message("  Panel A saved")
