################################################################################
#   Figure 6 — Composite Assembly
#   Sources all per-panel scripts, builds unified key strip, assembles with
#   patchwork into a single Figure_6.pdf/png.
#
#   Panels:
#     A — PCA of baseline eigengenes + LOOCV ROC (age discrimination)
#     B — Baseline eigengene associates with training response (scatter)
#     C — Training-induced module changes track phenotype (scatter)
#     D — kME x GS scatter (hub proteins)
#     E — Protein-phenotype heatmap
################################################################################

# === 0. Source all panel scripts (F6 setup loaded automatically by each) =====

source("04_Figures/F6/a_script/panel_A.R")     # → pA_pca, pA_roc, pA_combined
source("04_Figures/F6/a_script/panel_B.R")     # → pB
source("04_Figures/F6/a_script/panel_C.R")     # → pC
source("04_Figures/F6/a_script/panel_D.R")     # → pD
source("04_Figures/F6/a_script/panel_E.R")     # → pE

message("All F6 panel scripts sourced — assembling composite...")

# === 1. Build unified key strip =============================================

box_w <- 0.35
box_h <- KEY_BOX_HALF
item_gap <- 0.6

key_items <- list()
x_cursor <- 0

# Age
key_items <- c(key_items, list(tibble(
  x = x_cursor, y = 0, label = "Age:", type = "header", fill = NA_character_)))
x_cursor <- x_cursor + 0.9
for (nm in c("Young", "Old")) {
  key_items <- c(key_items, list(tibble(
    x = x_cursor, y = 0, label = NA_character_, type = "swatch",
    fill = AGE_COLORS[nm])))
  key_items <- c(key_items, list(tibble(
    x = x_cursor + box_w + 0.1, y = 0, label = nm, type = "item",
    fill = NA_character_)))
  x_cursor <- x_cursor + box_w + 0.1 + nchar(nm) * 0.22 + item_gap
}

# Top module eigengenes
x_cursor <- x_cursor + 0.4
key_items <- c(key_items, list(tibble(
  x = x_cursor, y = 0, label = "Top Modules:", type = "header", fill = NA_character_)))
x_cursor <- x_cursor + 2.2
for (mod_me in top3) {
  mod_col <- gsub("^ME", "", mod_me)
  key_items <- c(key_items, list(tibble(
    x = x_cursor, y = 0, label = NA_character_, type = "swatch",
    fill = mod_col)))
  key_items <- c(key_items, list(tibble(
    x = x_cursor + box_w + 0.1, y = 0,
    label = tools::toTitleCase(mod_col),
    type = "item", fill = NA_character_)))
  x_cursor <- x_cursor + box_w + 0.1 + nchar(mod_col) * 0.22 + item_gap
}

# Correlation key
x_cursor <- x_cursor + 0.4
key_items <- c(key_items, list(tibble(
  x = x_cursor, y = 0, label = "Pearson r:", type = "header", fill = NA_character_)))
x_cursor <- x_cursor + 1.8
z_cols <- c("#4393C3", "white", "#D6604D")
z_labs <- c("-1", "0", "+1")
for (zi in seq_along(z_cols)) {
  key_items <- c(key_items, list(tibble(
    x = x_cursor, y = 0, label = NA_character_, type = "swatch",
    fill = z_cols[zi])))
  key_items <- c(key_items, list(tibble(
    x = x_cursor + box_w + 0.05, y = 0, label = z_labs[zi],
    type = "item", fill = NA_character_)))
  x_cursor <- x_cursor + box_w + 0.05 + nchar(z_labs[zi]) * 0.22 + item_gap * 0.6
}

# Phenotype
x_cursor <- x_cursor + 0.4
key_items <- c(key_items, list(tibble(
  x = x_cursor, y = 0, label = "Phenotype:", type = "header", fill = NA_character_)))
x_cursor <- x_cursor + 1.9
for (nm in c("Delta VL", "Delta LBM")) {
  key_items <- c(key_items, list(tibble(
    x = x_cursor, y = 0, label = nm, type = "item", fill = NA_character_)))
  x_cursor <- x_cursor + nchar(nm) * 0.22 + item_gap
}

key_df <- bind_rows(key_items)
sw <- key_df %>% filter(type == "swatch")
hd <- key_df %>% filter(type == "header")
it <- key_df %>% filter(type == "item")

pKey <- ggplot() +
  { if (nrow(sw) > 0)
    geom_rect(data = sw,
              aes(xmin = x, xmax = x + box_w, ymin = -box_h, ymax = box_h),
              fill = sw$fill, color = KEY_BORDER, linewidth = KEY_LW) } +
  geom_text(data = hd, aes(x = x, y = 0, label = label),
            hjust = 0, size = KEY_TITLE, fontface = "bold", color = KEY_HDR_COL) +
  geom_text(data = it, aes(x = x, y = 0, label = label),
            hjust = 0, size = KEY_ITEM, fontface = "bold", color = KEY_ITEM_COL) +
  coord_cartesian(ylim = c(-1, 1), expand = FALSE) +
  theme_void() +
  theme(plot.margin = margin(t = 2, r = 4, b = 2, l = 4))

# === 2. Assemble composite ==================================================

# Row 1: A (PCA+ROC combined) | B (prediction scatter)
row1 <- (pA_combined | pB) + plot_layout(widths = c(0.50, 0.50))

# Row 2: C (plasticity scatter) | D (kME x GS) | E (heatmap)
row2 <- (pC | pD | pE) + plot_layout(widths = c(0.35, 0.30, 0.35))

fig6 <- row1 / row2 / pKey +
  plot_layout(heights = c(0.42, 0.45, 0.13)) +
  plot_annotation(
    title = "Translational Utility of WGCNA Modules",
    subtitle = "Module eigengenes discriminate age groups, associate with training response, and link to phenotypic change.",
    theme = theme(
      plot.title    = element_text(face = "bold", size = 11, hjust = 0.5),
      plot.subtitle = element_text(size = 8, color = "grey30", hjust = 0.5)
    )
  )

# === 3. Save =================================================================

ggsave(file.path(RPT_DIR, "Figure_6.pdf"), fig6,
       width = 380, height = 350, units = "mm",
       device = pdf, limitsize = FALSE)
ggsave(file.path(RPT_DIR, "Figure_6.png"), fig6,
       width = 380, height = 350, units = "mm",
       dpi = 300, limitsize = FALSE)

cat("Saved Figure_6.pdf and Figure_6.png\n")

# === 4. Supplementary Excel Workbook =========================================

library(openxlsx)

wb <- createWorkbook()

addWorksheet(wb, "A_age_discrimination")
addWorksheet(wb, "B_baseline_association")
addWorksheet(wb, "C_plasticity")
addWorksheet(wb, "D_hub_scatter")
addWorksheet(wb, "E_heatmap")
addWorksheet(wb, "_metadata")

# Panel A: PCA scores, ROC curve, permutation, pre-vs-post comparison
writeData(wb, "A_age_discrimination", pca_scores, startRow = 1)
writeData(wb, "A_age_discrimination", roc_df,
          startRow = nrow(pca_scores) + 3)
writeData(wb, "A_age_discrimination", perm_df,
          startRow = nrow(pca_scores) + nrow(roc_df) + 5)
writeData(wb, "A_age_discrimination", pre_post_df,
          startRow = nrow(pca_scores) + nrow(roc_df) + nrow(perm_df) + 7)

# Panel B: Baseline eigengene vs phenotype (long format + stats)
writeData(wb, "B_baseline_association", pred_long, startRow = 1)
writeData(wb, "B_baseline_association", pred_stats,
          startRow = nrow(pred_long) + 3)

# Panel C: Plasticity correlations (long format + stats)
writeData(wb, "C_plasticity", plast_long, startRow = 1)
writeData(wb, "C_plasticity", plast_stats,
          startRow = nrow(plast_long) + 3)

# Panel D: kME x GS hub scatter
writeData(wb, "D_hub_scatter", kme_gs_df, startRow = 1)

# Panel E: Protein-phenotype heatmap
writeData(wb, "E_heatmap", heat_long, startRow = 1)

# Metadata sheet
meta_info <- tibble(
  field = c("figure", "generated", "n_subjects", "n_modules", "n_proteins",
            "top3_modules", "phenotypes"),
  value = c("Figure 6: Translational Utility of WGCNA Modules",
            format(Sys.time(), "%Y-%m-%d %H:%M"),
            as.character(length(common_subj)),
            as.character(ncol(MEs)),
            as.character(ncol(datExpr)),
            paste(gsub("^ME", "", top3), collapse = ", "),
            "Delta VL (cm), Delta LBM (kg)")
)
writeData(wb, "_metadata", meta_info, startRow = 1)

saveWorkbook(wb, file.path(DAT_DIR, "Figure_6_data.xlsx"), overwrite = TRUE)
cat("Saved Figure_6_data.xlsx\n")
