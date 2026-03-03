################################################################################
#   Figure 3 — Composite Assembly
#   Sources all per-panel scripts, builds unified key strip, assembles with
#   patchwork into a single Figure_3.pdf/png + supplementary Excel.
#
#   Layout (3 rows + key strip):
#     Row 1 (33%): Panel A (50%) | Panel B (50%)     — 2 wide volcano rings
#     Row 2 (30%): Panel C (40%) | Panel D (60%)     — scatter + RRHO
#     Row 3 (31%): Panel E (full width)              — heatmap + Sankey + bars
#     Row 4 (6%):  Key strip footer
#
#   Panels:
#     A — Volcano ring: Aging Effect
#     B — Volcano ring: Reversal Effect (Training - Aging)
#     C — Reversal scatter (logFC Aging vs logFC Training_Old)
#     D — RRHO2 reversal map + ORA bars
#     E — Reversal classification: heatmap + Sankey + enrichment bars
################################################################################

# === 0. Source panel scripts (shared setup loaded automatically) ==============

source("05_Figures/F3/a_script/YvO_F3_setup.R")
source("05_Figures/F3/a_script/panel_A.R")   # -> pA
source("05_Figures/F3/a_script/panel_B.R")   # -> pB
source("05_Figures/F3/a_script/panel_C.R")   # -> pC (ggplot) + pC_combined (patchwork)
source("05_Figures/F3/a_script/panel_D.R")   # -> pD (patchwork: RRHO + ORA flanks)
source("05_Figures/F3/a_script/panel_E.R")   # -> pE (unified ggplot)

message("All F3 panel scripts sourced -- assembling composite...")

# === 1. Build unified key strip =============================================

KEY_HDR_COL  <- "grey25"
KEY_ITEM_COL <- "grey15"
KEY_BORDER   <- "grey70"

box_w <- 0.30
box_h <- KEY_BOX_HALF
item_gap <- 0.45

key_items <- list()
x_cursor <- 0

# Direction (volcano rings A/B)
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

# Database (ring term colors)
x_cursor <- x_cursor + 0.4
key_items <- c(key_items, list(tibble(
  x = x_cursor, y = 0, label = "Database:", type = "header", fill = NA_character_)))
x_cursor <- x_cursor + 1.8
for (nm in c("Hallmark", "GO:BP")) {
  key_items <- c(key_items, list(tibble(
    x = x_cursor, y = 0, label = NA_character_, type = "swatch",
    fill = DB_COLORS[nm])))
  key_items <- c(key_items, list(tibble(
    x = x_cursor + box_w + 0.1, y = 0, label = nm, type = "item",
    fill = NA_character_)))
  x_cursor <- x_cursor + box_w + 0.1 + nchar(nm) * 0.22 + item_gap
}

# Reversal patterns (Panel E classification)
x_cursor <- x_cursor + 0.4
key_items <- c(key_items, list(tibble(
  x = x_cursor, y = 0, label = "Reversal:", type = "header", fill = NA_character_)))
x_cursor <- x_cursor + 1.6
for (nm in c("Reversed", "Attenuated", "Amplified", "Concordant")) {
  key_items <- c(key_items, list(tibble(
    x = x_cursor, y = 0, label = NA_character_, type = "swatch",
    fill = PATTERN_COLORS[nm])))
  key_items <- c(key_items, list(tibble(
    x = x_cursor + box_w + 0.1, y = 0, label = nm, type = "item",
    fill = NA_character_)))
  x_cursor <- x_cursor + box_w + 0.1 + nchar(nm) * 0.22 + item_gap
}

# Significance (Panel C categories)
x_cursor <- x_cursor + 0.4
key_items <- c(key_items, list(tibble(
  x = x_cursor, y = 0, label = "Significance:", type = "header", fill = NA_character_)))
x_cursor <- x_cursor + 2.2
sig_nms <- c("Reversal", "Sig Both", "Sig Aging only", "Sig Training only")
for (nm in sig_nms) {
  key_items <- c(key_items, list(tibble(
    x = x_cursor, y = 0, label = NA_character_, type = "swatch",
    fill = SIG_COLORS[nm])))
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

# Wrap volcano rings and tag
pA <- wrap_elements(full = pA) + labs(tag = "A") + tag_theme
pB <- wrap_elements(full = pB) + labs(tag = "B") + tag_theme

# Panel C: use the base scatter (not pC_combined); strip title, remove legend
pC_comp <- pC + labs(tag = "C", title = NULL, subtitle = NULL) +
  tag_theme + theme(legend.position = "none",
                    axis.title = element_text(face = "bold"))

# Panel D: patchwork composite (RRHO + ORA bars) -- wrap and tag
pD <- pD +
  plot_annotation(title = NULL, subtitle = NULL) &
  theme(plot.title = element_blank(), plot.subtitle = element_blank())
pD <- wrap_elements(full = pD) + labs(tag = "D") + tag_theme

# Panel E: unified ggplot -- strip title/subtitle, wrap and tag
pE <- strip_plot_meta(pE)
pE <- wrap_elements(full = pE) + labs(tag = "E") + tag_theme

# === 3. Assemble composite (3 rows + key strip) =============================

# Row 1 (33%): A (50%) | B (50%)  — 2 wide volcano rings
# Row 2 (30%): C (40%) | D (60%)  — scatter + RRHO
# Row 3 (31%): E (full width)     — heatmap + Sankey + bars
# Row 4 (6%):  Key strip
design <- c(
  area(1,  1,  33, 50),    # A -- Aging volcano ring
  area(1,  51, 33, 100),   # B -- Reversal volcano ring
  area(35, 1,  64, 40),    # C -- reversal scatter
  area(35, 41, 64, 100),   # D -- RRHO + ORA bars
  area(66, 1,  94, 100),   # E -- heatmap + Sankey + pathway bars
  area(96, 1,  100, 100)   # Key strip
)

fig3 <- pA + pB + pC_comp + pD + pE + pKey +
  plot_layout(design = design)

# === 4. Save =================================================================

ggsave(file.path(RPT_DIR, "Figure_3.pdf"), fig3,
       width = 380, height = 400, units = "mm", device = pdf, limitsize = FALSE)
ggsave(file.path(RPT_DIR, "Figure_3.png"), fig3,
       width = 380, height = 400, units = "mm", dpi = 300, limitsize = FALSE)

cat("Saved Figure_3.pdf and Figure_3.png\n")

# === 5. Supplementary Excel workbook =========================================

wb <- createWorkbook()

addWorksheet(wb, "AB_volcano_data")
addWorksheet(wb, "C_reversal_scatter")
addWorksheet(wb, "C_reversal_scatter_stats")
addWorksheet(wb, "D_rrho2")
addWorksheet(wb, "D_rrho2_ora")
addWorksheet(wb, "E_classification")
addWorksheet(wb, "Reversal_tests")
addWorksheet(wb, "_metadata")

# Panels A/B: volcano data for Aging + Reversal
vol_aging <- read_csv(file.path(DAT_DIR, "panel_A", "volcano_aging.csv"),
                      show_col_types = FALSE) %>% mutate(contrast = "Aging")
vol_rev   <- read_csv(file.path(DAT_DIR, "panel_B", "volcano_reversal.csv"),
                      show_col_types = FALSE) %>% mutate(contrast = "Reversal")
writeData(wb, "AB_volcano_data", bind_rows(vol_aging, vol_rev))

# Panel C: reversal scatter
writeData(wb, "C_reversal_scatter",
          read_csv(file.path(DAT_DIR, "panel_C", "reversal_scatter.csv"),
                   show_col_types = FALSE))

# Panel C: reversal scatter statistics
writeData(wb, "C_reversal_scatter_stats",
          read_csv(file.path(DAT_DIR, "panel_C", "reversal_scatter_stats.csv"),
                   show_col_types = FALSE))

# Panel D: RRHO2 summary
writeData(wb, "D_rrho2",
          read_csv(file.path(DAT_DIR, "panel_D", "rrho2_summary.csv"),
                   show_col_types = FALSE))

# Panel D: RRHO2 ORA results (concordant + discordant combined)
d_ora <- tryCatch({
  bind_rows(
    read_csv(file.path(DAT_DIR, "panel_D", "rrho2_ora_concordant.csv"),
             show_col_types = FALSE),
    read_csv(file.path(DAT_DIR, "panel_D", "rrho2_ora_discordant.csv"),
             show_col_types = FALSE)
  )
}, error = function(e) tibble(note = "No ORA results"))
writeData(wb, "D_rrho2_ora", d_ora)

# Panel E: classification
writeData(wb, "E_classification",
          read_csv(file.path(DAT_DIR, "panel_E", "classification.csv"),
                   show_col_types = FALSE))

# Reversal tests (melov + contingency + signed reversal combined)
rev_tests <- tryCatch({
  bind_rows(
    read_csv(file.path(DAT_DIR, "reversal_tests", "melov_permutation.csv"),
             show_col_types = FALSE) %>% mutate(test = "Melov_permutation"),
    read_csv(file.path(DAT_DIR, "reversal_tests", "reversal_contingency.csv"),
             show_col_types = FALSE) %>% mutate(test = "Contingency_Fisher"),
    read_csv(file.path(DAT_DIR, "reversal_tests", "signed_reversal_score.csv"),
             show_col_types = FALSE) %>% mutate(test = "Signed_reversal")
  )
}, error = function(e) tibble(note = "Could not load reversal test CSVs"))
writeData(wb, "Reversal_tests", rev_tests)

# Metadata sheet
writeData(wb, "_metadata", tibble(
  field = c("figure", "generated", "script", "panels",
            "note_AB", "note_C", "note_C_stats", "note_D", "note_D_ora",
            "note_E", "note_reversal_tests"),
  value = c("Figure 3", format(Sys.time(), "%Y-%m-%d %H:%M"),
            "YvO_figure3_composite.R",
            paste("A: Volcano Aging, B: Volcano Reversal, C: Reversal scatter,",
                  "D: RRHO + ORA, E: Classification heatmap-Sankey-pathway"),
            "Volcano ring data for Aging and Reversal contrasts",
            "logFC reversal scatter (Aging vs Training_Old) with significance and imputation status",
            "Pearson r, Spearman rho, reversal %, and Melov magnitude reversal statistics",
            "RRHO2 hypergeometric overlap summary with hotspot gene counts and ORA pathway counts",
            "Per-quadrant ORA results from RRHO hotspot genes (reversal + exacerbation)",
            "4-pattern reversal classification (Reversed/Attenuated/Amplified/Concordant)",
            "Melov permutation, Fisher's exact contingency, and signed reversal score")
))

saveWorkbook(wb, file.path(DAT_DIR, "F3_supplementary.xlsx"), overwrite = TRUE)
cat("Saved F3_supplementary.xlsx\n")
cat("Figure 3 composite assembly complete.\n")
