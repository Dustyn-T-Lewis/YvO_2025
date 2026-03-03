################################################################################
#   Figure 2 — Panel A: Training Response (Young) Volcano Ring
#   Generates: panel_A_volcano.pdf/png, c_data/panel_A/volcano_young.csv
################################################################################

if (!exists("dep_df")) source("05_Figures/F2/a_script/YvO_F2_setup.R")
source("05_Figures/shared/volcano_ring.R")

message("Panel A: Training_Young volcano ring...")

# ── Build volcano ring ───────────────────────────────────────────────────────
pA <- make_volcano_ring(dep_df, fgsea_all,
                        contrast = "Training_Young",
                        title    = "Training Response (Young)")

# ── Save panel figures ───────────────────────────────────────────────────────
legend_A <- build_panel_legend()
pA_with_key <- pA / legend_A + plot_layout(heights = c(6, 1))

ggsave(file.path(RPT_DIR, "panel_A_volcano.pdf"), pA_with_key,
       width = 160, height = 180, units = "mm", device = cairo_pdf)
ggsave(file.path(RPT_DIR, "panel_A_volcano.png"), pA_with_key,
       width = 160, height = 180, units = "mm", dpi = 300)

# ── Export ring data ─────────────────────────────────────────────────────────
ring_data_A <- attr(pA, "ring_data")
if (!is.null(ring_data_A) && nrow(ring_data_A) > 0) {
  write_csv(ring_data_A %>% dplyr::select(-gene_list),
            file.path(DAT_DIR, "panel_A", "ring_terms.csv"))
}

# ── Export volcano CSV ───────────────────────────────────────────────────────
dep_df %>%
  transmute(
    gene,
    log2_fold_change = round(logFC_Training_Young, 4),
    neg_log10_pvalue = round(-log10(P.Value_Training_Young), 4),
    pi_score         = round(pi_score_Training_Young, 6),
    adjusted_pvalue  = round(adj.P.Val_Training_Young, 6),
    direction = case_when(
      pi_score_Training_Young < 0.05 & logFC_Training_Young > 0 ~ "Up",
      pi_score_Training_Young < 0.05 & logFC_Training_Young < 0 ~ "Down",
      TRUE ~ "NS"
    )
  ) %>%
  filter(!is.na(log2_fold_change), !is.na(neg_log10_pvalue)) %>%
  arrange(pi_score) %>%
  write_csv(file.path(DAT_DIR, "panel_A", "volcano_young.csv"))

message("  Panel A saved")
