################################################################################
#   Figure 2 — Composite Assembly
#   Sources panel scripts, builds unified key strip, assembles 3-row grid with
#   patchwork into a single Figure_2.pdf/png + supplementary Excel.
#
#   Layout:
#     Row 1 (30%): Panel A | Panel B | Panel C   (3 equal-width volcano rings)
#     Row 2 (30%): Panel D (40%) | Panel E (60%)
#     Row 3 (34%): Panel F (full width)
#     Row 4 (6%):  Key strip footer
#
#   Panels:
#     A — Volcano ring: Training Effect (Young)
#     B — Volcano ring: Training Effect (Old)
#     C — Volcano ring: Age x Training Interaction
#     D — Concordance scatter (logFC x logFC, with imputation key)
#     E — RRHO2 threshold-free concordance map + ORA bars
#     F — Heatmap -> Sankey -> pathway bars
################################################################################

# === 0. Source panel scripts (shared setup loaded automatically) ==============

source("05_Figures/F2/a_script/YvO_F2_setup.R")
source("05_Figures/F2/a_script/panel_A.R")   # -> pA
source("05_Figures/F2/a_script/panel_B.R")   # -> pB
source("05_Figures/F2/a_script/panel_C.R")   # -> pC
source("05_Figures/F2/a_script/panel_D.R")   # -> pD (ggplot) + pD_combined (patchwork)
source("05_Figures/F2/a_script/panel_E.R")   # -> pE (patchwork: RRHO + ORA flanks)
source("05_Figures/F2/a_script/panel_F.R")   # -> pF (unified ggplot)

message("All F2 panel scripts sourced -- assembling composite...")

# === 1. Build unified key strip =============================================

KEY_HDR_COL  <- "grey25"
KEY_ITEM_COL <- "grey15"
KEY_BORDER   <- "grey70"

box_w <- 0.30
box_h <- KEY_BOX_HALF
item_gap <- 0.45

key_items <- list()
x_cursor <- 0

# Direction (used in volcano rings A/B/C)
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

# Significance categories (used in Panel D scatter)
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

# Concordance/Discordance (RRHO quadrant shading)
x_cursor <- x_cursor + 0.4
key_items <- c(key_items, list(tibble(
  x = x_cursor, y = 0, label = "Quadrant:", type = "header", fill = NA_character_)))
x_cursor <- x_cursor + 1.6
for (nm in c("Concordant", "Discordant")) {
  fill_col <- if (nm == "Concordant") "#FFE0E0" else "#DCEEFF"
  key_items <- c(key_items, list(tibble(
    x = x_cursor, y = 0, label = NA_character_, type = "swatch",
    fill = fill_col)))
  key_items <- c(key_items, list(tibble(
    x = x_cursor + box_w + 0.1, y = 0, label = nm, type = "item",
    fill = NA_character_)))
  x_cursor <- x_cursor + box_w + 0.1 + nchar(nm) * 0.22 + item_gap
}

key_df_strip <- bind_rows(key_items)
sw <- key_df_strip %>% filter(type == "swatch")
hd <- key_df_strip %>% filter(type == "header")
it <- key_df_strip %>% filter(type == "item")

pKey <- ggplot() +
  geom_rect(data = sw,
            aes(xmin = x, xmax = x + box_w, ymin = -box_h, ymax = box_h),
            fill = sw$fill, color = KEY_BORDER, linewidth = KEY_LW) +
  geom_text(data = hd, aes(x = x, y = 0, label = label),
            hjust = 0, size = KEY_TITLE, fontface = "bold", color = KEY_HDR_COL) +
  geom_text(data = it, aes(x = x, y = 0, label = label),
            hjust = 0, size = KEY_TEXT, fontface = "bold", color = KEY_ITEM_COL) +
  coord_cartesian(xlim = c(-0.5, x_cursor + 0.5), ylim = c(-1, 1), expand = FALSE) +
  theme_void() +
  theme(plot.margin = margin(t = 2, r = 4, b = 2, l = 4))

# === 2. Prepare panels for composite ========================================

tag_theme <- theme(plot.tag = element_text(size = 14, face = "bold"))

# Strip titles from volcano ring panels for composite
pA <- strip_plot_meta(pA) + theme(plot.margin = margin(0, 0, 0, 0))
pB <- strip_plot_meta(pB) + theme(plot.margin = margin(0, 0, 0, 0))
pC <- strip_plot_meta(pC) + theme(plot.margin = margin(0, 0, 0, 0))

# Wrap volcano rings (patchwork objects) and tag
pA <- wrap_elements(full = pA) + labs(tag = "A") + tag_theme
pB <- wrap_elements(full = pB) + labs(tag = "B") + tag_theme
pC <- wrap_elements(full = pC) + labs(tag = "C") + tag_theme

# Panel D: strip title/subtitle, remove legend (key strip handles it)
pD_comp <- pD + labs(tag = "D", title = NULL, subtitle = NULL) +
  tag_theme + theme(legend.position = "none",
                    axis.title = element_text(face = "bold"))

# Panel E: strip inner titles from patchwork composite
pE <- pE +
  plot_annotation(title = NULL, subtitle = NULL) &
  theme(plot.title = element_blank(), plot.subtitle = element_blank())
pE <- wrap_elements(full = pE) + labs(tag = "E") + tag_theme

# Panel F: strip title/subtitle for composite
pF <- strip_plot_meta(pF)
pF <- wrap_elements(full = pF) + labs(tag = "F") + tag_theme

# === 3. Assemble composite (3 rows + key strip) =============================

# Row 1 (30%): A | B | C  (3 equal-width volcano rings)
# Row 2 (30%): D (40%) | E (60%)
# Row 3 (34%): F (full width)
# Row 4 (6%):  Key strip
design <- c(
  area(1,  1,  30, 33),   # A -- volcano Young
  area(1,  34, 30, 66),   # B -- volcano Old
  area(1,  67, 30, 100),  # C -- volcano Interaction
  area(32, 1,  60, 40),   # D -- concordance scatter
  area(32, 41, 60, 100),  # E -- RRHO + ORA bars
  area(62, 1,  94, 100),  # F -- heatmap -> sankey -> pathway bars
  area(96, 1,  100, 100)  # Key strip
)

fig2 <- pA + pB + pC + pD_comp + pE + pF + pKey +
  plot_layout(design = design)

# === 4. Save =================================================================

ggsave(file.path(RPT_DIR, "Figure_2.pdf"), fig2,
       width = 380, height = 400, units = "mm", device = pdf, limitsize = FALSE)
ggsave(file.path(RPT_DIR, "Figure_2.png"), fig2,
       width = 380, height = 400, units = "mm", dpi = 300, limitsize = FALSE)

cat("Saved Figure_2.pdf and Figure_2.png\n")

# === 5. Supplementary Excel workbook =========================================

wb <- createWorkbook()

addWorksheet(wb, "ABC_volcano_data")
addWorksheet(wb, "D_concordance")
addWorksheet(wb, "D_concordance_stats")
addWorksheet(wb, "E_rrho2")
addWorksheet(wb, "E_rrho2_ora")
addWorksheet(wb, "F_classification")
addWorksheet(wb, "_metadata")

# Panels A/B/C: volcano data for all three contrasts
vol_young <- read_csv(file.path(DAT_DIR, "panel_A", "volcano_young.csv"),
                      show_col_types = FALSE) %>% mutate(contrast = "Training_Young")
vol_old <- read_csv(file.path(DAT_DIR, "panel_B", "volcano_old.csv"),
                    show_col_types = FALSE) %>% mutate(contrast = "Training_Old")
vol_int <- read_csv(file.path(DAT_DIR, "panel_C", "volcano_interaction.csv"),
                    show_col_types = FALSE) %>% mutate(contrast = "Interaction")
writeData(wb, "ABC_volcano_data", bind_rows(vol_young, vol_old, vol_int))

# Panel D: concordance scatter + stats
writeData(wb, "D_concordance",
          read_csv(file.path(DAT_DIR, "panel_D", "concordance.csv"), show_col_types = FALSE))
writeData(wb, "D_concordance_stats",
          read_csv(file.path(DAT_DIR, "panel_D", "concordance_stats.csv"), show_col_types = FALSE))

# Panel E: RRHO2 summary
writeData(wb, "E_rrho2",
          read_csv(file.path(DAT_DIR, "panel_E", "rrho2_summary.csv"), show_col_types = FALSE))

# Panel E: RRHO2 ORA results (concordant + discordant)
e_ora <- tryCatch({
  bind_rows(
    read_csv(file.path(DAT_DIR, "panel_E", "rrho2_ora_concordant.csv"), show_col_types = FALSE),
    read_csv(file.path(DAT_DIR, "panel_E", "rrho2_ora_discordant.csv"), show_col_types = FALSE)
  )
}, error = function(e) tibble(note = "No ORA results"))
writeData(wb, "E_rrho2_ora", e_ora)

# Panel F: expanded classification
writeData(wb, "F_classification",
          read_csv(file.path(DAT_DIR, "panel_F", "expanded_classification.csv"),
                   show_col_types = FALSE))

# Metadata sheet
writeData(wb, "_metadata", tibble(
  field = c("figure", "generated", "script", "panels",
            "note_ABC", "note_D", "note_E", "note_E_ora", "note_F"),
  value = c("Figure 2", format(Sys.time(), "%Y-%m-%d %H:%M"),
            "YvO_figure2_composite.R",
            paste("A: Volcano Young, B: Volcano Old, C: Volcano Interaction,",
                  "D: Concordance scatter, E: RRHO + ORA, F: Heatmap-Sankey-Pathway"),
            "Volcano ring data for Training Young, Training Old, and Interaction contrasts",
            "logFC concordance scatter with hexbin NS + category overlays, imputation borders",
            "RRHO2 hypergeometric overlap summary with hotspot gene counts and ORA pathway counts",
            "Per-quadrant ORA results from RRHO hotspot gene extraction (concordant + discordant)",
            "Expanded 7-category classification with GO Slim super-category pathway assignment")
))

saveWorkbook(wb, file.path(DAT_DIR, "F2_supplementary.xlsx"), overwrite = TRUE)
cat("Saved F2_supplementary.xlsx\n")
