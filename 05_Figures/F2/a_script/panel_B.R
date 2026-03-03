################################################################################
#   Figure 2 — Panel B: Training Response (Old) Volcano Ring
#   Generates: panel_B_volcano.pdf/png, c_data/panel_B/volcano_old.csv
################################################################################

if (!exists("dep_df")) source("05_Figures/F2/a_script/YvO_F2_setup.R")
source("05_Figures/shared/volcano_ring.R")

message("Panel B: Training_Old volcano ring...")

# ── Build volcano ring ───────────────────────────────────────────────────────
pB <- make_volcano_ring(dep_df, fgsea_all,
                        contrast = "Training_Old",
                        title    = "Training Response (Old)")

# ── Save panel figures ───────────────────────────────────────────────────────
legend_B <- build_panel_legend()
pB_with_key <- pB / legend_B + plot_layout(heights = c(6, 1))

ggsave(file.path(RPT_DIR, "panel_B_volcano.pdf"), pB_with_key,
       width = 160, height = 180, units = "mm", device = cairo_pdf)
ggsave(file.path(RPT_DIR, "panel_B_volcano.png"), pB_with_key,
       width = 160, height = 180, units = "mm", dpi = 300)

# ── Export ring data ─────────────────────────────────────────────────────────
ring_data_B <- attr(pB, "ring_data")
if (!is.null(ring_data_B) && nrow(ring_data_B) > 0) {
  write_csv(ring_data_B %>% dplyr::select(-gene_list),
            file.path(DAT_DIR, "panel_B", "ring_terms.csv"))
}

# ── Export volcano CSV ───────────────────────────────────────────────────────
dep_df %>%
  transmute(
    gene,
    log2_fold_change = round(logFC_Training_Old, 4),
    neg_log10_pvalue = round(-log10(P.Value_Training_Old), 4),
    pi_score         = round(pi_score_Training_Old, 6),
    adjusted_pvalue  = round(adj.P.Val_Training_Old, 6),
    direction = case_when(
      pi_score_Training_Old < 0.05 & logFC_Training_Old > 0 ~ "Up",
      pi_score_Training_Old < 0.05 & logFC_Training_Old < 0 ~ "Down",
      TRUE ~ "NS"
    )
  ) %>%
  filter(!is.na(log2_fold_change), !is.na(neg_log10_pvalue)) %>%
  arrange(pi_score) %>%
  write_csv(file.path(DAT_DIR, "panel_B", "volcano_old.csv"))

message("  Panel B saved")
