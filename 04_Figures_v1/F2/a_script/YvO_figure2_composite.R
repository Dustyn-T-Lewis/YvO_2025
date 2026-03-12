################################################################################
#   Figure 2 — Composite Assembly
#   Sources panel scripts, assembles figure with patchwork + unified key panel.
#   Panels keep their own titles/subtitles (no stripping).
#
#   Layout (100x100 virtual grid):
#     A: Volcano Young    (1-30,  1-25)
#     B: Volcano Old      (31-55, 1-25)
#     C: Volcano Interact (56-85, 1-25)
#     D: Concordance      (1-50,  26-55)
#     E: RRHO + ORA       (51-85, 26-55)
#     F: Heatmap-Sankey   (1-85,  56-100)
#     Key Panel           (86-100, 1-100)
#
#   Target: 420 x 480 mm
################################################################################

# === 0. Source panel scripts (shared setup loaded automatically) ==============

source("04_Figures/F2/a_script/YvO_F2_setup.R")
source("04_Figures/F2/a_script/panel_AB.R")   # -> pA, pB, pC_volcano
source("04_Figures/F2/a_script/panel_C.R")     # -> pC (concordance scatter)
source("04_Figures/F2/a_script/panel_D.R")     # -> pD (RRHO + ORA patchwork)
source("04_Figures/F2/a_script/panel_F.R")     # -> pF (heatmap-sankey)
source("04_Figures/F2/a_script/panel_G.R")     # -> pG (per-contrast ORA)

message("All F2 panel scripts sourced -- assembling composite...")

# === 1. Prepare panels for composite =========================================

# Tag theme for panel letters
tag_theme <- theme(plot.tag = element_text(size = TXT_TAG, face = "bold"))

# Volcano panels: wrap for layout (titles inside ring already)
pA_comp <- wrap_elements(full = pA) + labs(tag = "A") + tag_theme
pB_comp <- wrap_elements(full = pB) + labs(tag = "B") + tag_theme
pC_vol_comp <- wrap_elements(full = pC_volcano) + labs(tag = "C") + tag_theme

# Concordance scatter: keep title/subtitle visible
pD_comp <- pC + labs(tag = "D") + tag_theme

# RRHO + ORA: wrap patchwork composite
pE_comp <- wrap_elements(full = pD) + labs(tag = "E") + tag_theme

# Heatmap-Sankey: wrap for layout
pF_comp <- wrap_elements(full = pF) + labs(tag = "F") + tag_theme

# ORA dot plot: wrap for layout
pG_comp <- wrap_elements(full = pG) + labs(tag = "G") + tag_theme

# === 2. Build Key Panel =======================================================

message("Building unified key panel...")

# Key panel: all legends for the figure in a single theme_void ggplot
# Organized in sections: Volcano Ring | Significance | Scatter | RRHO | ORA | Heatmap | Super-categories

key_x <- 0   # current x cursor
ky <- 1.0    # row height
kpt <- 3.0   # point size
ktxt <- 3.2  # text size
khdr <- 3.8  # header text size

key_annotations <- list()
ka <- function(...) key_annotations[[length(key_annotations) + 1]] <<- list(...)

# --- Section 1: Significance (5-level) ---
ka(type = "header", x = 0, y = 5 * ky, label = "Significance", size = khdr)
sig_names <- c("Interaction", "Sig Both", "Sig Young only", "Sig Old only", "NS")
for (i in seq_along(sig_names)) {
  ka(type = "point", x = 0, y = (5 - i) * ky, fill = SIG_COLORS[sig_names[i]],
     shape = 21, size = kpt, color = "black")
  ka(type = "text", x = 0.4, y = (5 - i) * ky, label = sig_names[i], size = ktxt)
}

# --- Section 2: Scatter Quadrants ---
ka(type = "header", x = 5.5, y = 5 * ky, label = "Scatter Quadrants", size = khdr)
ka(type = "rect", x = 5.5, y = 4 * ky, fill = "#FFE0E0", w = 0.3, h = 0.25)
ka(type = "text", x = 6.0, y = 4 * ky, label = "Concordant", size = ktxt)
ka(type = "rect", x = 5.5, y = 3 * ky, fill = "#DCEEFF", w = 0.3, h = 0.25)
ka(type = "text", x = 6.0, y = 3 * ky, label = "Discordant", size = ktxt)

# --- Section 3: RRHO ---
ka(type = "header", x = 10, y = 5 * ky, label = "RRHO Overlap", size = khdr)
ka(type = "text", x = 10, y = 4 * ky,
   label = "Viridis: -log10(p) hypergeometric", size = ktxt)

# --- Section 4: ORA Quadrants ---
ka(type = "header", x = 10, y = 3 * ky, label = "ORA Quadrants", size = khdr)
ora_names <- names(ORA_QUAD_COLORS)
ora_short <- c("Concordant Up", "Concordant Down", "Y Up / O Down", "Y Down / O Up")
for (i in seq_along(ora_names)) {
  ka(type = "rect", x = 10 + (i - 1) * 3.2, y = 2 * ky,
     fill = ORA_QUAD_COLORS[ora_names[i]], w = 0.3, h = 0.25)
  ka(type = "text", x = 10.4 + (i - 1) * 3.2, y = 2 * ky,
     label = ora_short[i], size = ktxt)
}

# --- Section 5: Heatmap logFC ---
ka(type = "header", x = 23, y = 5 * ky, label = "Heatmap logFC", size = khdr)
ka(type = "text", x = 23, y = 4 * ky,
   label = "Blue (-) | White (0) | Red (+)", size = ktxt)

# --- Section 6: Consolidated pathways ---
ka(type = "header", x = 23, y = 3 * ky, label = "Consolidated Pathways", size = khdr)
sc_names <- setdiff(CONSOLIDATED_PATHWAY_ORDER, "Other")
sc_names <- sc_names[sc_names %in% names(CONSOLIDATED_COLORS)]
n_sc_cols <- 3
for (i in seq_along(sc_names)) {
  col_idx <- (i - 1) %% n_sc_cols
  row_idx <- (i - 1) %/% n_sc_cols
  sx <- 23 + col_idx * 5.5
  sy <- 2 * ky - row_idx * ky
  ka(type = "rect", x = sx, y = sy, fill = CONSOLIDATED_COLORS[sc_names[i]], w = 0.3, h = 0.25)
  ka(type = "text", x = sx + 0.4, y = sy, label = sc_names[i], size = 2.5)
}

# Build key ggplot
pKey <- ggplot() + theme_void() +
  coord_cartesian(xlim = c(-0.5, 40), ylim = c(-3 * ky, 6 * ky), expand = FALSE) +
  theme(plot.margin = margin(2, 2, 2, 2, "mm"))

for (a in key_annotations) {
  if (a$type == "header") {
    pKey <- pKey + annotate("text", x = a$x, y = a$y, label = a$label,
                            hjust = 0, size = a$size, fontface = "bold",
                            color = "grey25")
  } else if (a$type == "text") {
    pKey <- pKey + annotate("text", x = a$x, y = a$y, label = a$label,
                            hjust = 0, size = a$size, color = "grey15")
  } else if (a$type == "point") {
    pKey <- pKey + annotate("point", x = a$x, y = a$y,
                            shape = a$shape, size = a$size,
                            fill = a$fill, color = a$color, stroke = 0.5)
  } else if (a$type == "rect") {
    pKey <- pKey + annotate("rect",
                            xmin = a$x - a$w/2, xmax = a$x + a$w/2,
                            ymin = a$y - a$h/2, ymax = a$y + a$h/2,
                            fill = a$fill, color = "grey50", linewidth = 0.2)
  }
}

pKey_comp <- wrap_elements(full = pKey)

# === 3. Build layout with patchwork::area() ===================================

layout <- c(
  area(t =  1, l =  1, b = 25, r = 23),   # A: Volcano Young
  area(t = 26, l =  1, b = 48, r = 23),   # B: Volcano Old
  area(t = 49, l =  1, b = 72, r = 23),   # C: Volcano Interaction
  area(t =  1, l = 24, b = 40, r = 47),   # D: Concordance scatter
  area(t = 41, l = 24, b = 72, r = 47),   # E: RRHO + ORA
  area(t =  1, l = 48, b = 72, r = 100),  # F: Heatmap-Sankey (wider, 53%)
  area(t = 73, l =  1, b = 90, r = 100),  # G: ORA dot plot
  area(t = 91, l =  1, b = 100, r = 100)  # Key Panel
)

fig2 <- pA_comp + pB_comp + pC_vol_comp + pD_comp + pE_comp + pF_comp + pG_comp + pKey_comp +
  plot_layout(design = layout)

# === 4. Save ==================================================================

ggsave(file.path(RPT_DIR, "Figure_2.pdf"), fig2,
       width = 420, height = 560, units = "mm", device = cairo_pdf, limitsize = FALSE)
ggsave(file.path(RPT_DIR, "Figure_2.png"), fig2,
       width = 420, height = 560, units = "mm", dpi = 300, limitsize = FALSE)

cat("Saved Figure_2.pdf and Figure_2.png\n")

# === 5. Supplementary Excel workbook =========================================

library(openxlsx)

wb <- createWorkbook()

addWorksheet(wb, "A_volcano_young")
addWorksheet(wb, "B_volcano_old")
addWorksheet(wb, "C_volcano_interaction")
addWorksheet(wb, "D_concordance")
addWorksheet(wb, "D_concordance_stats")
addWorksheet(wb, "E_rrho2")
addWorksheet(wb, "E_rrho2_ora")
addWorksheet(wb, "F_classification")
addWorksheet(wb, "F_sankey_links")
addWorksheet(wb, "G_ora_per_contrast")
addWorksheet(wb, "_metadata")

# Panel A: volcano data (Training Young)
vol_young <- tryCatch(
  read_csv(file.path(DAT_DIR, "panel_A", "volcano_young.csv"), show_col_types = FALSE),
  error = function(e) tibble(note = "File not found"))
writeData(wb, "A_volcano_young", vol_young)

# Panel B: volcano data (Training Old)
vol_old <- tryCatch(
  read_csv(file.path(DAT_DIR, "panel_B", "volcano_old.csv"), show_col_types = FALSE),
  error = function(e) tibble(note = "File not found"))
writeData(wb, "B_volcano_old", vol_old)

# Panel C: volcano data (Interaction)
vol_int <- tryCatch(
  read_csv(file.path(DAT_DIR, "panel_C", "ring_terms.csv"), show_col_types = FALSE),
  error = function(e) tibble(note = "Interaction ring terms not available"))
writeData(wb, "C_volcano_interaction", vol_int)

# Panel D: concordance scatter + stats
writeData(wb, "D_concordance",
          tryCatch(read_csv(file.path(DAT_DIR, "panel_C", "concordance.csv"), show_col_types = FALSE),
                   error = function(e) tibble(note = "File not found")))
writeData(wb, "D_concordance_stats",
          tryCatch(read_csv(file.path(DAT_DIR, "panel_C", "concordance_stats.csv"), show_col_types = FALSE),
                   error = function(e) tibble(note = "File not found")))

# Panel E: RRHO2 summary + ORA
writeData(wb, "E_rrho2",
          tryCatch(read_csv(file.path(DAT_DIR, "panel_D", "rrho2_summary.csv"), show_col_types = FALSE),
                   error = function(e) tibble(note = "File not found")))
e_ora <- tryCatch({
  bind_rows(
    read_csv(file.path(DAT_DIR, "panel_D", "rrho2_ora_concordant.csv"), show_col_types = FALSE),
    read_csv(file.path(DAT_DIR, "panel_D", "rrho2_ora_discordant.csv"), show_col_types = FALSE)
  )
}, error = function(e) tibble(note = "No ORA results"))
writeData(wb, "E_rrho2_ora", e_ora)

# Panel F: classification + sankey links
writeData(wb, "F_classification",
          tryCatch(read_csv(file.path(DAT_DIR, "panel_F", "expanded_classification.csv"), show_col_types = FALSE),
                   error = function(e) tibble(note = "File not found")))
writeData(wb, "F_sankey_links",
          tryCatch(read_csv(file.path(DAT_DIR, "panel_F", "sankey_links.csv"), show_col_types = FALSE),
                   error = function(e) tibble(note = "File not found")))

# Panel G: ORA per contrast
writeData(wb, "G_ora_per_contrast",
          tryCatch(read_csv(file.path(DAT_DIR, "panel_G", "ora_per_contrast.csv"), show_col_types = FALSE),
                   error = function(e) tibble(note = "File not found")))

# Metadata sheet
writeData(wb, "_metadata", tibble(
  field = c("figure", "generated", "script", "panels",
            "note_A", "note_B", "note_C", "note_D", "note_E", "note_F", "note_G"),
  value = c("Figure 2", format(Sys.time(), "%Y-%m-%d %H:%M"),
            "YvO_figure2_composite.R",
            "A-C: Volcano rings (Young/Old/Interaction), D: Concordance, E: RRHO+ORA, F: Heatmap-Sankey, G: ORA",
            "Volcano ring: Training Effect (Young) -- Post_Young - Pre_Young",
            "Volcano ring: Training Effect (Old) -- Post_Old - Pre_Old",
            "Volcano ring: Interaction Effect -- (Old-Young) x (Post-Pre)",
            "logFC concordance scatter with symmetric axes, significance categories",
            "RRHO2 hypergeometric overlap + 4-color per-quadrant ORA bars",
            "Expanded 7-category heatmap + consolidated pathway Sankey + gene count bars",
            "Per-contrast ORA bar chart: Hallmark + GO:BP (top 3 per contrast)")
))

saveWorkbook(wb, file.path(DAT_DIR, "F2_supplementary.xlsx"), overwrite = TRUE)
cat("Saved F2_supplementary.xlsx\n")
