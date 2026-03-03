################################################################################
#   Figure 2 — Panel C: Age x Training Interaction Volcano Ring
#   Generates: panel_C_volcano.pdf/png, c_data/panel_C/volcano_interaction.csv
################################################################################

if (!exists("dep_df")) source("05_Figures/F2/a_script/YvO_F2_setup.R")
source("05_Figures/shared/volcano_ring.R")

message("Panel C: Interaction volcano ring...")

# ── Build volcano ring ───────────────────────────────────────────────────────
pC <- make_volcano_ring(dep_df, fgsea_all,
                        contrast = "Interaction",
                        title    = "Age \u00d7 Training Interaction")

# ── Save panel figures ───────────────────────────────────────────────────────
legend_C <- build_panel_legend()
pC_with_key <- pC / legend_C + plot_layout(heights = c(6, 1))

ggsave(file.path(RPT_DIR, "panel_C_volcano.pdf"), pC_with_key,
       width = 160, height = 180, units = "mm", device = cairo_pdf)
ggsave(file.path(RPT_DIR, "panel_C_volcano.png"), pC_with_key,
       width = 160, height = 180, units = "mm", dpi = 300)

# ── Export ring data ─────────────────────────────────────────────────────────
ring_data_C <- attr(pC, "ring_data")
if (!is.null(ring_data_C) && nrow(ring_data_C) > 0) {
  write_csv(ring_data_C %>% dplyr::select(-gene_list),
            file.path(DAT_DIR, "panel_C", "ring_terms.csv"))
}

# ── Export volcano CSV ───────────────────────────────────────────────────────
dep_df %>%
  transmute(
    gene,
    log2_fold_change = round(logFC_Interaction, 4),
    neg_log10_pvalue = round(-log10(P.Value_Interaction), 4),
    pi_score         = round(pi_score_Interaction, 6),
    adjusted_pvalue  = round(adj.P.Val_Interaction, 6),
    direction = case_when(
      pi_score_Interaction < 0.05 & logFC_Interaction > 0 ~ "Up",
      pi_score_Interaction < 0.05 & logFC_Interaction < 0 ~ "Down",
      TRUE ~ "NS"
    )
  ) %>%
  filter(!is.na(log2_fold_change), !is.na(neg_log10_pvalue)) %>%
  arrange(pi_score) %>%
  write_csv(file.path(DAT_DIR, "panel_C", "volcano_interaction.csv"))

message("  Panel C saved")
