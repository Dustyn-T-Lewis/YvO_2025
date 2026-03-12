################################################################################
#   Figure 3 — Composite Assembly
#   Sources panel scripts, assembles figure with patchwork + unified key panel.
#   Panels keep their own titles/subtitles (no stripping).
#   NES scatter (panel_E.R) is supplementary only — not sourced here.
#
#   Layout (100x100 virtual grid):
#     A: Volcano Aging     (1-40,  1-28)
#     B: Volcano Reversal  (41-70, 1-28)
#     C: Reversal scatter  (1-70,  29-58)
#     D: RRHO + ORA        (71-85, 1-58)
#     E: Classification    (1-85,  59-100)
#     Key Panel            (86-100, 1-100)
#
#   Target: 420 x 420 mm
################################################################################

# === 0. Source panel scripts (shared setup loaded automatically) ==============

source("04_Figures/F3/a_script/YvO_F3_setup.R")
source("04_Figures/F3/a_script/panel_AB.R")    # -> pA, pB
source("04_Figures/F3/a_script/panel_C.R")      # -> pC (reversal scatter)
source("04_Figures/F3/a_script/panel_D.R")      # -> pD (RRHO + ORA patchwork)
source("04_Figures/F3/a_script/panel_F.R")      # -> pF (classification triptych)
source("04_Figures/F3/a_script/panel_G.R")      # -> pG (per-contrast ORA)

message("All F3 panel scripts sourced -- assembling composite...")

# === 1. Prepare panels for composite =========================================

# Tag theme for panel letters
tag_theme <- theme(plot.tag = element_text(size = TXT_TAG, face = "bold"))

# Volcano panels: wrap for layout (titles inside ring already)
pA_comp <- wrap_elements(full = pA) + labs(tag = "A") + tag_theme
pB_comp <- wrap_elements(full = pB) + labs(tag = "B") + tag_theme

# Reversal scatter: keep title/subtitle visible
pC_comp <- pC + labs(tag = "C") + tag_theme

# RRHO + ORA: wrap patchwork composite
pD_comp <- wrap_elements(full = pD) + labs(tag = "D") + tag_theme

# Classification triptych: wrap for layout
pE_comp <- wrap_elements(full = pF) + labs(tag = "E") + tag_theme

# ORA dot plot: wrap for layout
pF_ora_comp <- wrap_elements(full = pG) + labs(tag = "F") + tag_theme

# === 2. Build Key Panel =======================================================

message("Building unified key panel...")

# Key panel: all legends for the figure in a single theme_void ggplot
# Organized in sections: Significance | Scatter | RRHO | ORA | Contrasts | Super-categories

ky <- 1.0    # row height
kpt <- 3.0   # point size
ktxt <- 3.2  # text size
khdr <- 3.8  # header text size

key_annotations <- list()
ka <- function(...) key_annotations[[length(key_annotations) + 1]] <<- list(...)

# --- Section 1: Significance (5-level, F3 terminology) ---
ka(type = "header", x = 0, y = 5 * ky, label = "Significance", size = khdr)
sig_names <- c("Reversal", "Sig Both", "Sig Aging only", "Sig Training only", "NS")
for (i in seq_along(sig_names)) {
  ka(type = "point", x = 0, y = (5 - i) * ky, fill = SIG_COLORS[sig_names[i]],
     shape = 21, size = kpt, color = "black")
  ka(type = "text", x = 0.4, y = (5 - i) * ky, label = sig_names[i], size = ktxt)
}

# --- Section 2: Scatter Quadrants ---
ka(type = "header", x = 5.5, y = 5 * ky, label = "Scatter Quadrants", size = khdr)
ka(type = "rect", x = 5.5, y = 4 * ky, fill = "#DCEEFF", w = 0.3, h = 0.25)
ka(type = "text", x = 6.0, y = 4 * ky, label = "Reversal", size = ktxt)
ka(type = "rect", x = 5.5, y = 3 * ky, fill = "#FFE0E0", w = 0.3, h = 0.25)
ka(type = "text", x = 6.0, y = 3 * ky, label = "Exacerbation", size = ktxt)

# --- Section 3: RRHO ---
ka(type = "header", x = 10, y = 5 * ky, label = "RRHO Overlap", size = khdr)
ka(type = "text", x = 10, y = 4 * ky,
   label = "Viridis: -log10(p) hypergeometric", size = ktxt)

# --- Section 4: ORA Quadrants ---
ka(type = "header", x = 10, y = 3 * ky, label = "ORA Quadrants", size = khdr)
ora_names <- names(ORA_QUAD_COLORS)
ora_short <- c("Reversed Ag Up/Tr Down", "Reversed Ag Down/Tr Up",
               "Exacerbated Up", "Exacerbated Down")
for (i in seq_along(ora_names)) {
  ka(type = "rect", x = 10 + (i - 1) * 3.5, y = 2 * ky,
     fill = ORA_QUAD_COLORS[ora_names[i]], w = 0.3, h = 0.25)
  ka(type = "text", x = 10.4 + (i - 1) * 3.5, y = 2 * ky,
     label = ora_short[i], size = ktxt)
}

# --- Section 5: Heatmap Classification ---
ka(type = "header", x = 23, y = 5 * ky, label = "Classification", size = khdr)
ka(type = "rect", x = 23, y = 4 * ky, fill = "#D6604D", w = 0.3, h = 0.25)
ka(type = "text", x = 23.4, y = 4 * ky, label = "Aging Up", size = ktxt)
ka(type = "rect", x = 26, y = 4 * ky, fill = "#4393C3", w = 0.3, h = 0.25)
ka(type = "text", x = 26.4, y = 4 * ky, label = "Aging Down", size = ktxt)
ka(type = "rect", x = 29, y = 4 * ky, fill = "#E65100", w = 0.3, h = 0.25)
ka(type = "text", x = 29.4, y = 4 * ky, label = "Rev. (Age Up)", size = ktxt)
ka(type = "rect", x = 33, y = 4 * ky, fill = "#FFB300", w = 0.3, h = 0.25)
ka(type = "text", x = 33.4, y = 4 * ky, label = "Rev. (Age Down)", size = ktxt)

# --- Section 5b: Heatmap logFC gradient ---
ka(type = "rect", x = 23, y = 3.3 * ky, fill = "#2166AC", w = 0.3, h = 0.25)
ka(type = "text", x = 23.4, y = 3.3 * ky, label = "Down", size = ktxt)
ka(type = "rect", x = 25, y = 3.3 * ky, fill = "white", w = 0.3, h = 0.25)
ka(type = "text", x = 25.4, y = 3.3 * ky, label = "0", size = ktxt)
ka(type = "rect", x = 26.5, y = 3.3 * ky, fill = "#B2182B", w = 0.3, h = 0.25)
ka(type = "text", x = 26.9, y = 3.3 * ky, label = "Up (log2FC)", size = ktxt)

# --- Section 6: Consolidated pathways ---
ka(type = "header", x = 23, y = 2.5 * ky, label = "Consolidated Pathways", size = khdr)
sc_names <- setdiff(CONSOLIDATED_PATHWAY_ORDER, "Other")
sc_names <- sc_names[sc_names %in% names(CONSOLIDATED_COLORS)]
n_sc_cols <- 3
for (i in seq_along(sc_names)) {
  col_idx <- (i - 1) %% n_sc_cols
  row_idx <- (i - 1) %/% n_sc_cols
  sx <- 23 + col_idx * 5.5
  sy <- 1.5 * ky - row_idx * ky
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
  area(t =  1, l =  1, b = 33, r = 24),   # A: Volcano Aging (top-left)
  area(t = 34, l =  1, b = 58, r = 24),   # B: Volcano Reversal (mid-left)
  area(t =  1, l = 25, b = 58, r = 49),   # C: Scatter (top-middle, spans to B)
  area(t = 59, l =  1, b = 72, r = 49),   # D: RRHO + ORA (bottom, left+middle)
  area(t =  1, l = 50, b = 72, r = 100),  # E: Classification (wider, 51%)
  area(t = 73, l =  1, b = 90, r = 100),  # F: ORA dot plot
  area(t = 91, l =  1, b = 100, r = 100)  # Key Panel
)

fig3 <- pA_comp + pB_comp + pC_comp + pD_comp + pE_comp + pF_ora_comp + pKey_comp +
  plot_layout(design = layout)

# === 4. Save ==================================================================

ggsave(file.path(RPT_DIR, "Figure_3.pdf"), fig3,
       width = 420, height = 500, units = "mm", device = cairo_pdf, limitsize = FALSE)
ggsave(file.path(RPT_DIR, "Figure_3.png"), fig3,
       width = 420, height = 500, units = "mm", dpi = 300, limitsize = FALSE)

cat("Saved Figure_3.pdf and Figure_3.png\n")

# === 5. Supplementary Excel workbook =========================================

library(openxlsx)

wb <- createWorkbook()

addWorksheet(wb, "A_volcano_aging")
addWorksheet(wb, "B_volcano_reversal")
addWorksheet(wb, "C_reversal_scatter")
addWorksheet(wb, "D_rrho2")
addWorksheet(wb, "D_rrho2_ora")
addWorksheet(wb, "E_classification")
addWorksheet(wb, "F_ora_per_contrast")
addWorksheet(wb, "_metadata")

# Panel A: volcano data (Aging)
writeData(wb, "A_volcano_aging",
          tryCatch(read_csv(file.path(DAT_DIR, "panel_A", "volcano_aging.csv"), show_col_types = FALSE),
                   error = function(e) tibble(note = "File not found")))

# Panel B: volcano data (Reversal)
writeData(wb, "B_volcano_reversal",
          tryCatch(read_csv(file.path(DAT_DIR, "panel_B", "volcano_reversal.csv"), show_col_types = FALSE),
                   error = function(e) tibble(note = "File not found")))

# Panel C: reversal scatter
writeData(wb, "C_reversal_scatter",
          tryCatch(read_csv(file.path(DAT_DIR, "panel_C", "reversal_scatter.csv"), show_col_types = FALSE),
                   error = function(e) tibble(note = "File not found")))

# Panel D: RRHO2 summary + ORA
writeData(wb, "D_rrho2",
          tryCatch(read_csv(file.path(DAT_DIR, "panel_D", "rrho2_summary.csv"), show_col_types = FALSE),
                   error = function(e) tibble(note = "File not found")))
d_ora <- tryCatch({
  bind_rows(
    read_csv(file.path(DAT_DIR, "panel_D", "rrho2_ora_concordant.csv"), show_col_types = FALSE),
    read_csv(file.path(DAT_DIR, "panel_D", "rrho2_ora_discordant.csv"), show_col_types = FALSE)
  )
}, error = function(e) tibble(note = "No ORA results"))
writeData(wb, "D_rrho2_ora", d_ora)

# Panel E: classification (sankey links)
writeData(wb, "E_classification",
          tryCatch(read_csv(file.path(DAT_DIR, "panel_F", "sankey_links.csv"), show_col_types = FALSE),
                   error = function(e) tibble(note = "File not found")))

# Panel F: ORA per contrast
writeData(wb, "F_ora_per_contrast",
          tryCatch(read_csv(file.path(DAT_DIR, "panel_G", "ora_per_contrast.csv"), show_col_types = FALSE),
                   error = function(e) tibble(note = "File not found")))

# Metadata sheet
writeData(wb, "_metadata", tibble(
  field = c("figure", "generated", "script", "panels",
            "note_A", "note_B", "note_C", "note_D", "note_E", "note_F"),
  value = c("Figure 3", format(Sys.time(), "%Y-%m-%d %H:%M"),
            "YvO_figure3_composite.R",
            "A: Volcano Aging, B: Volcano Reversal, C: Reversal scatter, D: RRHO+ORA, E: Classification, F: ORA",
            "Volcano ring: Aging Effect -- Old_Pre - Young_Pre",
            "Volcano ring: Reversal (Training - Aging)",
            "logFC reversal scatter (Aging vs Training_Old) with significance categories",
            "RRHO2 hypergeometric overlap + 4-color per-quadrant ORA bars (horizontal layout)",
            "Aging & Reversal: heatmap + consolidated pathway Sankey + gene count bars",
            "Per-contrast ORA bar chart: Hallmark + GO:BP (top 3 per contrast)")
))

saveWorkbook(wb, file.path(DAT_DIR, "F3_supplementary.xlsx"), overwrite = TRUE)
cat("Saved F3_supplementary.xlsx\n")
