################################################################################
#   Figure 1 — Composite Assembly + Supplementary Excel Workbook
#
#   Sources all per-panel scripts, assembles with patchwork into a single
#   Figure_1.pdf/png, and exports F1_supplementary.xlsx.
#
#   Panels:
#     A — CV% violins (inter-individual variability)
#     B — logFC density histograms (effect size distributions)
#     C — PCA biplot + PERMANOVA + legend key
#     D — DEPs per contrast (pseudo-log stacked bar)
#     E — UpSet plot (dual bar + dot matrix)
#     F — fGSEA grouped bar chart (pathway enrichment)
#     S1.7 — Pi-score distributions (supplementary)
################################################################################

# === 0. Source all panel scripts (setup loaded automatically by panel_A) =====

source("04_Figures/F1/a_script/panel_A.R")     # -> pA
source("04_Figures/F1/a_script/panel_B.R")     # -> pB
source("04_Figures/F1/a_script/panel_C.R")     # -> pC, pC_key
source("04_Figures/F1/a_script/panel_D.R")     # -> pD  (+ sig_sets, dir_map for E)
source("04_Figures/F1/a_script/panel_E.R")     # -> pE_bars, pE_dots, pKeys
source("04_Figures/F1/a_script/panel_F.R")     # -> pF
source("04_Figures/F1/a_script/supp_S1_7.R")   # -> s17

message("All F1 panel scripts sourced - assembling composite...")

# === 1. Composite figure (patchwork layout from original monolith) ===========

left_col  <- (pA / pB / pC / pC_key) +
              plot_layout(heights = c(0.20, 0.38, 0.35, 0.07))
mid_col   <- (pD / pE_bars / pE_dots / pKeys) +
              plot_layout(heights = c(0.22, 0.46, 0.18, 0.14))
fig1 <- (left_col | mid_col | pF) +
  plot_layout(widths = c(0.32, 0.43, 0.25)) +
  plot_annotation(
    title = "Proteomic Landscape of Aging and Resistance Training Response",
    subtitle = "Data quality, effect magnitude, contrast overlap, and functional enrichment across the 2 \u00d7 2 factorial design.",
    theme = theme(
      plot.title    = element_text(face = "bold", size = 11, hjust = 0.5),
      plot.subtitle = element_text(size = 8, color = "grey30", hjust = 0.5)
    )
  )

ggsave(file.path(RPT_DIR, "Figure_1.pdf"), fig1,
       width = 380, height = 260, units = "mm", device = pdf)
ggsave(file.path(RPT_DIR, "Figure_1.png"), fig1,
       width = 380, height = 260, units = "mm", dpi = 300)

cat("Saved Figure_1.pdf and Figure_1.png\n")

# === 2. Supplementary Excel workbook =========================================

library(openxlsx)

# Sheet 1: A_variability — CV data from panel A
sheet_A <- cv_df |>
  left_join(cv_med, by = "group")
names(sheet_A)[names(sheet_A) == "cv"]  <- "cv_pct"
names(sheet_A)[names(sheet_A) == "med"] <- "group_median_cv"

# Sheet 2: B_effect_sizes — logFC summary from panel B
sheet_B <- lfc_stats |>
  dplyr::select(contrast, med_abs_lfc, n_above_05)

# Sheet 3: C_pca_scores — PCA scores from panel C
sheet_C <- pca_df |>
  dplyr::select(sample_id, PC1, PC2, age, time, group)

# Sheet 4: D_dep_fractions — DEP fraction data from panel D
sheet_D <- frac_df |>
  dplyr::select(contrast, threshold, n, pct)

# Sheet 5: E_upset_intersections — intersection count data from panel E
#   Build a tidy table of intersection sizes with set membership
upset_tbl <- tibble(
  intersection_id = seq_len(n_int),
  comb_code       = comb_names_ord,
  up_count        = up_ord,
  down_count      = down_ord,
  total           = up_ord + down_ord
)
# Add set membership columns
for (j in seq_along(set_order_ch)) {
  upset_tbl[[set_order_ch[j]]] <- vapply(comb_names_ord, function(code) {
    as.logical(as.integer(strsplit(code, "")[[1]])[j])
  }, logical(1))
}
sheet_E <- upset_tbl

# Sheet 6: F_fgsea_results — full fgsea results from panel F
sheet_F <- fgsea_export

# Sheet 7: _metadata — figure description
sheet_meta <- tibble(
  item = c(
    "figure_title",
    "description",
    "contrasts",
    "statistical_methods",
    "sample_sizes",
    "significance_threshold",
    "software",
    "date_generated"
  ),
  value = c(
    "Figure 1: Proteomic Landscape of Aging and Resistance Training Response",
    "Data quality (CV%), effect magnitude (logFC), PCA, DEP counts, contrast overlap (UpSet), and pathway enrichment (fGSEA) across the 2x2 factorial design.",
    "Aging (Old_Pre - Young_Pre), Training_Young (Young_Post - Young_Pre), Training_Old (Old_Post - Old_Pre), Interaction ((Old_Post-Old_Pre) - (Young_Post-Young_Pre))",
    "Wilcoxon rank-sum (CV%), pairwise Wilcoxon (logFC), PERMANOVA (PCA), pi-score (DEPs), fGSEA multilevel (pathways), rrvgo GO reduction (threshold 0.7)",
    sprintf("%d samples, %d proteins (filtered), %d DEP results",
            nrow(meta), nrow(norm_df), nrow(dep_df)),
    "Pi-score < 0.05 (DEPs), padj < 0.05 (fGSEA)",
    "R, limma, fgsea, rrvgo, patchwork",
    format(Sys.Date(), "%Y-%m-%d")
  )
)

wb <- createWorkbook()
addWorksheet(wb, "A_variability");         writeData(wb, "A_variability",         sheet_A)
addWorksheet(wb, "B_effect_sizes");        writeData(wb, "B_effect_sizes",        sheet_B)
addWorksheet(wb, "C_pca_scores");          writeData(wb, "C_pca_scores",          sheet_C)
addWorksheet(wb, "D_dep_fractions");       writeData(wb, "D_dep_fractions",       sheet_D)
addWorksheet(wb, "E_upset_intersections"); writeData(wb, "E_upset_intersections", sheet_E)
addWorksheet(wb, "F_fgsea_results");       writeData(wb, "F_fgsea_results",       sheet_F)
addWorksheet(wb, "_metadata");             writeData(wb, "_metadata",             sheet_meta)

# Bold headers on every sheet
hs <- createStyle(textDecoration = "bold")
for (sn in names(wb)) {
  addStyle(wb, sn, hs, rows = 1, cols = seq_len(20), gridExpand = TRUE)
  setColWidths(wb, sn, cols = seq_len(20), widths = "auto")
}

xlsx_path <- file.path(DAT_DIR, "F1_supplementary.xlsx")
saveWorkbook(wb, xlsx_path, overwrite = TRUE)
cat(sprintf("Saved %s (%d sheets)\n", xlsx_path, length(names(wb))))

cat("Figure 1 composite assembly complete.\n")
