################################################################################
#   Figure 5 — Composite Assembly
#   Sources WGCNA runner (if needed) + all per-panel scripts, builds unified
#   key strip, assembles with patchwork into a single Figure_5.pdf/png.
#
#   Panels:
#     A — Module dendrogram with color bar
#     B — Module-trait correlation heatmap
#     C — Triptych: z-score heatmap | eigengene slopes | enrichment bars
#     D — GO enrichment dotplot (top 5 per module)
#     E — Hub protein network (2x3 grid per module)
#     F — Module preservation (Young vs Old)
################################################################################

# === 0. Ensure WGCNA outputs exist ==========================================

wgcna_cache <- "04_Figures/F5/c_data/wgcna/wgcna_network.rds"
if (!file.exists(wgcna_cache)) {
  message("WGCNA cache not found — running WGCNA pipeline...")
  source("04_Figures/F5/a_script/YvO_WGCNA_run.R")
}

# === 1. Source all panel scripts (F5 setup loaded automatically by each) =====

source("04_Figures/F5/a_script/panel_A.R")     # → pA
source("04_Figures/F5/a_script/panel_B.R")     # → pB
source("04_Figures/F5/a_script/panel_C.R")     # → panel_C
source("04_Figures/F5/a_script/panel_D.R")     # → pD
source("04_Figures/F5/a_script/panel_E.R")     # → panel_E
source("04_Figures/F5/a_script/panel_F.R")     # → pF

message("All F5 panel scripts sourced — assembling composite...")

# === 2. Build unified key strip =============================================

box_w <- 0.35
box_h <- KEY_BOX_HALF
item_gap <- 0.6

key_items <- list()
x_cursor <- 0

# Module colors
mod_names <- sort(unique(module_df$module_color[module_df$module_color != "grey"]))
key_items <- c(key_items, list(tibble(
  x = x_cursor, y = 0, label = "Module:", type = "header", fill = NA_character_)))
x_cursor <- x_cursor + 1.5
for (nm in mod_names) {
  key_items <- c(key_items, list(tibble(
    x = x_cursor, y = 0, label = NA_character_, type = "swatch",
    fill = nm)))
  key_items <- c(key_items, list(tibble(
    x = x_cursor + box_w + 0.1, y = 0, label = tools::toTitleCase(nm),
    type = "item", fill = NA_character_)))
  x_cursor <- x_cursor + box_w + 0.1 + nchar(nm) * 0.22 + item_gap
}

# Age
x_cursor <- x_cursor + 0.4
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

# Trait correlation
x_cursor <- x_cursor + 0.4
key_items <- c(key_items, list(tibble(
  x = x_cursor, y = 0, label = "Trait r:", type = "header", fill = NA_character_)))
x_cursor <- x_cursor + 1.3
z_cols <- c("#2166AC", "white", "#B2182B")
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

key_df <- bind_rows(key_items)
sw <- key_df %>% filter(type == "swatch")
hd <- key_df %>% filter(type == "header")
it <- key_df %>% filter(type == "item")

pKey <- ggplot() +
  geom_rect(data = sw,
            aes(xmin = x, xmax = x + box_w, ymin = -box_h, ymax = box_h),
            fill = sw$fill, color = KEY_BORDER, linewidth = KEY_LW) +
  geom_text(data = hd, aes(x = x, y = 0, label = label),
            hjust = 0, size = KEY_TITLE, fontface = "bold", color = KEY_HDR_COL) +
  geom_text(data = it, aes(x = x, y = 0, label = label),
            hjust = 0, size = KEY_ITEM, fontface = "bold", color = KEY_ITEM_COL) +
  coord_cartesian(ylim = c(-1, 1), expand = FALSE) +
  theme_void() +
  theme(plot.margin = margin(t = 2, r = 4, b = 2, l = 4))

# === 3. Assemble composite ==================================================

# Row 1: A (dendrogram) | B (heatmap)
row1 <- (pA | pB) + plot_layout(widths = c(0.40, 0.60))

# Row 2: C (triptych — spanning full width)
row2 <- wrap_elements(full = panel_C)

# Row 3: D (GO enrichment) | E (hub network) | F (preservation)
row3 <- (pD | wrap_elements(full = panel_E) | pF) +
  plot_layout(widths = c(0.28, 0.42, 0.30))

fig5 <- row1 / row2 / row3 / pKey +
  plot_layout(heights = c(0.15, 0.40, 0.35, 0.10)) +
  plot_annotation(
    title = "WGCNA Co-Expression Network: Module Structure and Functional Characterization",
    subtitle = sprintf("Signed hybrid network | %d modules from %d proteins | %d samples",
                       n_modules, ncol(datExpr), nrow(datExpr)),
    theme = theme(
      plot.title    = element_text(face = "bold", size = 11, hjust = 0.5),
      plot.subtitle = element_text(size = 8, color = "grey30", hjust = 0.5)
    )
  )

# === 4. Save =================================================================

ggsave(file.path(RPT_DIR, "Figure_5.pdf"), fig5,
       width = 600, height = 800, units = "mm",
       device = pdf, limitsize = FALSE)
ggsave(file.path(RPT_DIR, "Figure_5.png"), fig5,
       width = 600, height = 800, units = "mm",
       dpi = 300, limitsize = FALSE)

cat("Saved Figure_5.pdf and Figure_5.png\n")

# === 5. Supplementary Excel Workbook ========================================

library(openxlsx)

wb <- createWorkbook()

addWorksheet(wb, "A_dendrogram")
addWorksheet(wb, "B_heatmap")
addWorksheet(wb, "C_eigengenes")
addWorksheet(wb, "C_enrichment")
addWorksheet(wb, "D_go_enrichment")
addWorksheet(wb, "E_hub_network")
addWorksheet(wb, "F_preservation")
addWorksheet(wb, "_metadata")

# Panel A: dendrogram module assignments
writeData(wb, "A_dendrogram", dendro_data)

# Panel B: module-trait heatmap
writeData(wb, "B_heatmap", heat_df)

# Panel C: module eigengene data
writeData(wb, "C_eigengenes", eigen_all)

# Panel C: triptych GO enrichment
writeData(wb, "C_enrichment", all_enrich)

# Panel D: GO enrichment dotplot
writeData(wb, "D_go_enrichment", go_plot)

# Panel E: hub network nodes (hub + context ring combined)
hub_combined <- bind_rows(
  all_node_data  %>% mutate(node_type = "hub"),
  all_periph_data %>% mutate(node_type = "periphery")
)
writeData(wb, "E_hub_network", hub_combined)

# Panel F: module preservation
writeData(wb, "F_preservation", pres_df)

# Metadata sheet
meta_sheet <- tibble(
  sheet       = c("A_dendrogram", "B_heatmap", "C_eigengenes", "C_enrichment",
                   "D_go_enrichment", "E_hub_network", "F_preservation"),
  description = c(
    "Module assignments from WGCNA dendrogram (unmerged + merged colors)",
    "Module-trait correlation heatmap (Pearson r, BH-corrected p-values)",
    "Module eigengene values per sample for 6 key modules",
    "GO enrichment results for triptych panel (top terms per module)",
    "GO enrichment dotplot data (top 5 per key module, BP preferred)",
    "Hub protein network nodes (top 30 by kME) + periphery ring proteins",
    "Module preservation Zsummary (Pre -> Post training, 200 permutations)"
  ),
  source_panel = c("A", "B", "C", "C", "D", "E", "F"),
  n_rows       = c(nrow(dendro_data), nrow(heat_df), nrow(eigen_all),
                    nrow(all_enrich), nrow(go_plot), nrow(hub_combined),
                    nrow(pres_df))
)
writeData(wb, "_metadata", meta_sheet)

saveWorkbook(wb, file.path(DAT_DIR, "Figure_5_data.xlsx"), overwrite = TRUE)
cat("Saved Figure_5_data.xlsx\n")
