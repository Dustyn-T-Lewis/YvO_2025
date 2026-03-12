################################################################################
#   Figure 2 — Composite Assembly
#   Sources all per-panel scripts, builds unified key strip, assembles with
#   patchwork into a single Figure_2.pdf/png.
#
#   Panels:
#     A — Volcano ring: Training Effect (Young)
#     B — Volcano ring: Training Effect (Old)
#     C — Concordance scatter (logFC x logFC)
#     D — RRHO2 threshold-free concordance map
#     E — fGSEA NES scatter (Hallmark + rrvgo-reduced GO:BP)
#     F — Interaction DEPs: dumbbell | sankey | enrichment bar
################################################################################

# === 0. Source all panel scripts (shared setup loaded automatically) =========

source("04_Figures/F2/a_script/YvO_F2_setup.R")
source("04_Figures/F2/a_script/panel_AB.R")   # → pA, pB
source("04_Figures/F2/a_script/panel_C.R")     # → pC, pC_combined, pC_imp_key_strip
source("04_Figures/F2/a_script/panel_D.R")     # → pD
source("04_Figures/F2/a_script/panel_E.R")     # → pE, pE_key, pE_combined
source("04_Figures/F2/a_script/panel_F.R")     # → pF

message("All F2 panel scripts sourced — assembling composite...")

# === 1. Build unified key strip =============================================

# Key elements: Direction (Up/Down), Significance categories, Database,
# Imputation (MAR/MNAR/Observed), Set size, Interaction pattern

box_w <- 0.35
box_h <- KEY_BOX_HALF
item_gap <- 0.6

key_items <- list()
x_cursor <- 0

# Direction
key_items <- c(key_items, list(tibble(
  x = x_cursor, y = 0, label = "Direction:", type = "header", fill = NA_character_)))
x_cursor <- x_cursor + 1.8
for (nm in c("Up", "Down")) {
  key_items <- c(key_items, list(tibble(
    x = x_cursor, y = 0, label = NA_character_, type = "swatch",
    fill = DIR_COLORS[nm])))
  key_items <- c(key_items, list(tibble(
    x = x_cursor + box_w + 0.1, y = 0, label = nm, type = "item",
    fill = NA_character_)))
  x_cursor <- x_cursor + box_w + 0.1 + nchar(nm) * 0.22 + item_gap
}

# Database
x_cursor <- x_cursor + 0.4
key_items <- c(key_items, list(tibble(
  x = x_cursor, y = 0, label = "Database:", type = "header", fill = NA_character_)))
x_cursor <- x_cursor + 1.8
for (nm in c("Hallmark", "GO:BP")) {
  key_items <- c(key_items, list(tibble(
    x = x_cursor, y = 0, label = NA_character_, type = "swatch",
    fill = if (nm == "Hallmark") "#AA336A" else "#00796B")))
  key_items <- c(key_items, list(tibble(
    x = x_cursor + box_w + 0.1, y = 0, label = nm, type = "item",
    fill = NA_character_)))
  x_cursor <- x_cursor + box_w + 0.1 + nchar(nm) * 0.22 + item_gap
}

# Significance (from Panel E)
x_cursor <- x_cursor + 0.4
key_items <- c(key_items, list(tibble(
  x = x_cursor, y = 0, label = "Significance:", type = "header", fill = NA_character_)))
x_cursor <- x_cursor + 2.2
sig_nms <- c("Interaction", "Sig Both", "Sig Young only", "Sig Old only")
for (nm in sig_nms) {
  key_items <- c(key_items, list(tibble(
    x = x_cursor, y = 0, label = NA_character_, type = "swatch",
    fill = SIG_COLORS[nm])))
  key_items <- c(key_items, list(tibble(
    x = x_cursor + box_w + 0.1, y = 0, label = nm, type = "item",
    fill = NA_character_)))
  x_cursor <- x_cursor + box_w + 0.1 + nchar(nm) * 0.22 + item_gap
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

# === 2. Assemble composite ==================================================

# Row 1: A | B | C (volcanos + concordance scatter)
row1 <- (pA | pB | pC) + plot_layout(widths = c(0.33, 0.33, 0.34))

# Row 2: D | E | F (RRHO + NES scatter + interaction triptych)
row2 <- (pD | pE | pF) + plot_layout(widths = c(0.25, 0.35, 0.40))

fig2 <- row1 / row2 / pKey +
  plot_layout(heights = c(0.44, 0.44, 0.12)) +
  plot_annotation(
    title = "Training-Effect Concordance Across Age Groups",
    subtitle = "Protein-level and pathway-level concordance between young and old training responses, plus interaction DEPs.",
    theme = theme(
      plot.title    = element_text(face = "bold", size = 11, hjust = 0.5),
      plot.subtitle = element_text(size = 8, color = "grey30", hjust = 0.5)
    )
  )

# === 3. Save =================================================================

ggsave(file.path(RPT_DIR, "Figure_2.pdf"), fig2,
       width = 380, height = 350, units = "mm", device = pdf, limitsize = FALSE)
ggsave(file.path(RPT_DIR, "Figure_2.png"), fig2,
       width = 380, height = 350, units = "mm", dpi = 300, limitsize = FALSE)

cat("Saved Figure_2.pdf and Figure_2.png\n")

# === 4. Supplementary Excel workbook =========================================

library(openxlsx)

wb <- createWorkbook()

addWorksheet(wb, "AB_volcano_data")
addWorksheet(wb, "C_concordance")
addWorksheet(wb, "D_rrho2")
addWorksheet(wb, "E_nes_scatter")
addWorksheet(wb, "F_interaction")
addWorksheet(wb, "_metadata")

# Panel AB: volcano data for both contrasts (re-read from saved CSVs)
vol_young <- read_csv(file.path(DAT_DIR, "panel_A", "volcano_young.csv"), show_col_types = FALSE) %>%
  mutate(contrast = "Training_Young")
vol_old <- read_csv(file.path(DAT_DIR, "panel_B", "volcano_old.csv"), show_col_types = FALSE) %>%
  mutate(contrast = "Training_Old")
writeData(wb, "AB_volcano_data", bind_rows(vol_young, vol_old))

# Panel C: concordance scatter
writeData(wb, "C_concordance",
          read_csv(file.path(DAT_DIR, "panel_C", "concordance.csv"), show_col_types = FALSE))

# Panel D: RRHO2 summary
writeData(wb, "D_rrho2",
          read_csv(file.path(DAT_DIR, "panel_D", "rrho2_summary.csv"), show_col_types = FALSE))

# Panel E: NES scatter
writeData(wb, "E_nes_scatter",
          read_csv(file.path(DAT_DIR, "panel_E", "nes_scatter.csv"), show_col_types = FALSE))

# Panel F: interaction classification + sankey links
writeData(wb, "F_interaction",
          read_csv(file.path(DAT_DIR, "panel_F", "sankey_links.csv"), show_col_types = FALSE))

# Metadata sheet
writeData(wb, "_metadata", tibble(
  field = c("figure", "generated", "script", "panels",
            "note_AB", "note_C", "note_D", "note_E", "note_F"),
  value = c("Figure 2", format(Sys.time(), "%Y-%m-%d %H:%M"),
            "YvO_figure2_composite.R",
            "A: Volcano Young, B: Volcano Old, C: Concordance, D: RRHO, E: NES scatter, F: Interaction",
            "Volcano ring data for Training Young and Training Old contrasts",
            "logFC concordance scatter with significance and imputation status",
            "RRHO2 hypergeometric overlap summary (full matrix in panel_D/rrho2_matrix.csv)",
            "fGSEA NES scatter for Hallmark + rrvgo-reduced GO:BP pathways",
            "Interaction DEP gene-pathway links (1:1 mapping) with response classification")
))

saveWorkbook(wb, file.path(DAT_DIR, "F2_supplementary.xlsx"), overwrite = TRUE)
cat("Saved F2_supplementary.xlsx\n")
