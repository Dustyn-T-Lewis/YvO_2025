# F04 — Composite: Training Concordance (Panels A–E)
# 3-column geometry-aware layout (panels labeled in visual reading order):
#   Top row:    A (Quadrant ORA)  | B (NES scatter)  | C (Pattern heatmap, full-height)
#   Bottom row: D (fry barcode)   | E (RRHO2)         | C (continues)
# Outputs: MAIN_F04_composite.{pdf,png}
#
# Source order: A, B, D, E, then C. Panel C (pattern heatmap) loads
# AnnotationDbi (via go_slim_categories.R); the S4 select() masking is now
# repaired inside go_slim_categories.R itself (re-attaches dplyr after dispatch).
# C is still sourced last to preserve stat-snapshot variable ordering.

setwd(rprojroot::find_rstudio_root_file())
source("04_Figures/shared/style.R")

suppressPackageStartupMessages({
  library(patchwork)
  library(cowplot)
})

# -- Source panels (C last for stat-snapshot ordering) -----------------------

message("=== F04 Composite: sourcing panels ===")
source("04_Figures/F04/a_script/main/panels/panel_A_ORA.R")        # -> composite (patchwork)
# Snapshot Panel A stats before Panel C overwrites `n_total` and `sig_df`
n_total_A    <- nrow(scatter_df)
n_sig_A      <- n_sig
n_enrich_A   <- n_enrich
r_spear_A    <- r_spear
source("04_Figures/F04/a_script/main/panels/panel_B_nes_scatter.R")  # -> pB (ggplot)
# Snapshot Panel B stats (survive, but alias for clarity)
n_pw_B       <- nrow(fgsea_wide)
n_sig_pw_B   <- n_total_sig
rho_B        <- as.numeric(nes_cor_all$estimate)
rho_lo_B     <- nes_ci_all[1]
rho_hi_B     <- nes_ci_all[2]
pw_conc_B    <- pw_conc_frac
source("04_Figures/F04/a_script/main/panels/panel_D_fry.R")        # -> pD_fry, p1, p2, p_t
# Snapshot Panel D stats
cor_imp_D    <- cor_imp
n_all_D      <- n_all
source("04_Figures/F04/a_script/main/panels/panel_E_rrho2.R")      # -> pE_heat (ggplot)
# Snapshot Panel E stats
n_shared_E   <- n_shared
max_UU_E     <- max_UU
n_UU_E       <- n_UU
source("04_Figures/F04/a_script/main/panels/panel_C_pattern_heatmap.R")  # -> p (ggplot, theme_void)
RPT_PDF <- "04_Figures/F04/b_reports/main/pdf"
RPT_PNG <- "04_Figures/F04/b_reports/main/png"
dir.create(RPT_PDF, recursive = TRUE, showWarnings = FALSE)
dir.create(RPT_PNG, recursive = TRUE, showWarnings = FALSE)
pdf_device <- get_pdf_device()

# -- Layout constants --------------------------------------------------------

COMP_W     <- 420     # mm (wider to fit equal-size D + E panels)
COMP_H     <- 310     # mm
PRINT_SCALE <- 380 / 178  # 380mm source → 178mm print = 2.13x
TAG_SZ     <- round(10 * PRINT_SCALE * 0.85)   # uniform tag size (~18pt)
TTL_SZ     <- round(10 * PRINT_SCALE * 0.85)   # uniform title size (~18pt)
SUB_SZ     <- round(7 * PRINT_SCALE * 0.85)    # uniform subtitle size (~13pt)

# -- Stat-snapshot title strings --------------------------------------------

ttl_A <- "Quadrant ORA (Concordance)"
sub_A <- sprintf("N = %d | %d DEPs (\u03a0) | %d enriched (FDR) | \u03c1 = %.2f",
                 n_total_A, n_sig_A, n_enrich_A, r_spear_A)
ttl_B <- "Protein-to-Pathway"
sub_B <- sprintf("%d proteins | %d pathways", n_total, n_pw)
ttl_C <- "fry: Concordance"
sub_C <- sprintf("n = %d | dupCor = %.3f", n_all_D, cor_imp_D)
ttl_D <- "Pathway Concordance"
sub_D <- sprintf("\u03c1 = %.2f | %.0f%% concordant",
                 rho_B, pw_conc_B * 100)
ttl_E <- "RRHO2 Concordance"
sub_E <- sprintf("%d genes | max %d", n_shared_E, n_UU_E)

# -- Compose: 2-row layout with title spacers, A+B top-aligned --------------
# 15-row, 13-col grid. B extends 2 rows past A for extra height.
# A=6 rows, B=8 rows (top-aligned, B taller). C/D/E = 6 rows.
#
# Patchwork assigns design characters ALPHABETICALLY by addition order.
# Addition order: composite(A), p(B), pD_fry(C), pB(D), pE_heat(E)

layout <- paste(
  "##############",   # row 1:  title spacer (a + b titles)
  "AAAAAAAABBBBBB",   # rows 2-7: A (8 cols) + B (6 cols)
  "AAAAAAAABBBBBB",
  "AAAAAAAABBBBBB",
  "AAAAAAAABBBBBB",
  "AAAAAAAABBBBBB",
  "AAAAAAAABBBBBB",
  "##############",   # row 8: spacer (B same height as A)
  "##############",   # row 9: mid title spacer (c + d + e titles)
  "CCCCCCDDDDEEEE",   # rows 10-15: C=6 (wider fry), D=4, E=4 (D+E equal, adjacent)
  "CCCCCCDDDDEEEE",
  "CCCCCCDDDDEEEE",
  "CCCCCCDDDDEEEE",
  "CCCCCCDDDDEEEE",
  "CCCCCCDDDDEEEE",
  sep = "\n"
)

# Expand Panel A slightly: +5mm height, +2mm width
composite <- composite + plot_annotation(theme = theme(plot.margin = margin(-2.5, -1, -2.5, -1, "mm")))

# Margins: consistent top so C/D/E tops align; spacer row provides title room
pD_fry  <- pD_fry + plot_annotation(theme = theme(plot.margin = margin(3, 5, 0, 0, "mm")))
pB      <- pB + theme(plot.margin = margin(-2.8, 5, 2.8, -5, "mm"))  # shift 5mm left, 2.8mm up
pE_heat <- pE_heat + theme(plot.margin = margin(-2.1, -0.2, 3.4, -3.5, "mm"),
                           axis.title = element_text(face = "bold", size = 8))  # bold + visible axis titles
# Shift Panel B heatmap ~5mm right by shifting coord viewport left
p <- p + coord_cartesian(xlim = c(-0.25, X_BAR_MAX + 1.75),
                         ylim = c(BAR_YMAX + ROW_H * 7.5, -ROW_H * 0.05),
                         expand = FALSE)

fig <- wrap_elements(full = composite) +   # A: top-left (ORA, 8 cols)
       wrap_elements(full = p) +           # B: top-right (heatmap, 6 cols)
       wrap_elements(full = pD_fry) +      # C: bottom-left (fry, 6 cols — wider)
       wrap_elements(full = pB) +          # D: bottom-mid (NES scatter, 4 cols)
       wrap_elements(full = pE_heat) +     # E: bottom-right (RRHO2, 4 cols — equal to D)
       plot_layout(design = layout,
                   widths  = rep(1, 14),
                   heights = c(6.5, rep(10, 6), 4, 4.5, rep(12, 6)))

# -- Panel tags, titles, subtitles via cowplot -------------------------------
# 15-row, 14-col grid. Heights: [4, 10×6, 4, 7, 12×6] = 147 total
#   Top spacer:    y  1.000 → 0.973   (4/147)
#   A rows (6):    y  0.973 → 0.565   (60/147)
#   B extension:   y  0.565 → 0.537   (4/147)   — B only, A empty (taller for heatmap)
#   Mid spacer:    y  0.537 → 0.490   (7/147)   — C/D/E titles + subtitles
#   Bottom panels: y  0.490 → 0.000   (72/147)
X_A <- 0.005;  X_B <- 0.549;  X_C <- 0.012;  X_D <- 0.406;  X_E <- 0.693
X_TTL      <- 0.030        # title/subtitle offset right of tag
TAG_DY     <- -0.002       # raise tag for baseline alignment
SUB_OFFSET <- 0.020        # subtitle gap below title

# Per-panel Y positions — titles in spacer, well above panel tops
Y_A <- 0.984;  Y_B <- 0.984           # a + b: shifted down with panels
Y_C <- 0.512;  Y_D <- 0.511;  Y_E <- 0.511   # C down 0.5mm more; D/E down 2mm more

composite_final <- ggdraw(fig) +
  # Panel a (ORA scatter): top-left
  draw_label("A",    x = X_A,           y = Y_A - TAG_DY,       size = TAG_SZ, fontface = "bold",        hjust = 0, vjust = 1) +
  draw_label(ttl_A,  x = X_A + X_TTL,   y = Y_A,                size = TTL_SZ, fontface = "bold",        hjust = 0, vjust = 1) +
  draw_label(sub_A,  x = X_A + X_TTL,   y = Y_A - SUB_OFFSET,   size = SUB_SZ, fontface = "bold.italic", hjust = 0, vjust = 1, colour = "grey40") +
  # Panel b (pattern heatmap): top-right
  draw_label("B",    x = X_B,           y = Y_B - TAG_DY,       size = TAG_SZ, fontface = "bold",        hjust = 0, vjust = 1) +
  draw_label(ttl_B,  x = X_B + X_TTL,   y = Y_B,                size = TTL_SZ, fontface = "bold",        hjust = 0, vjust = 1) +
  draw_label(sub_B,  x = X_B + X_TTL,   y = Y_B - SUB_OFFSET, size = SUB_SZ, fontface = "bold.italic", hjust = 0, vjust = 1, colour = "grey40") +
  # Panel c (fry barcode): bottom-left
  draw_label("C",    x = X_C,           y = Y_C - TAG_DY,       size = TAG_SZ, fontface = "bold",        hjust = 0, vjust = 1) +
  draw_label(ttl_C,  x = X_C + X_TTL,   y = Y_C,                size = TTL_SZ, fontface = "bold",        hjust = 0, vjust = 1) +
  draw_label(sub_C,  x = X_C + X_TTL,   y = Y_C - SUB_OFFSET,   size = SUB_SZ, fontface = "bold.italic", hjust = 0, vjust = 1, colour = "grey40") +
  # Panel d (NES scatter): bottom-middle
  draw_label("D",    x = X_D,           y = Y_D - TAG_DY,       size = TAG_SZ, fontface = "bold",        hjust = 0, vjust = 1) +
  draw_label(ttl_D,  x = X_D + X_TTL,   y = Y_D,                size = TTL_SZ, fontface = "bold",        hjust = 0, vjust = 1) +
  draw_label(sub_D,  x = X_D + X_TTL,   y = Y_D - SUB_OFFSET,   size = SUB_SZ, fontface = "bold.italic", hjust = 0, vjust = 1, colour = "grey40") +
  # Panel e (RRHO2): bottom-right
  draw_label("E",    x = X_E,           y = Y_E - TAG_DY,       size = TAG_SZ, fontface = "bold",        hjust = 0, vjust = 1) +
  draw_label(ttl_E,  x = X_E + X_TTL,   y = Y_E,                size = TTL_SZ, fontface = "bold",        hjust = 0, vjust = 1) +
  draw_label(sub_E,  x = X_E + X_TTL,   y = Y_E - SUB_OFFSET,   size = SUB_SZ, fontface = "bold.italic", hjust = 0, vjust = 1, colour = "grey40")

# -- Overlay raster legend PNGs below panels --
# (Panel D and E legends are now part of their plot objects — no overlay needed)

# -- Save -------------------------------------------------------------------

ggsave(file.path(RPT_PDF, "MAIN_F04_composite.pdf"), composite_final,
       width = COMP_W, height = COMP_H, units = "mm", device = pdf_device)
ggsave(file.path(RPT_PNG, "MAIN_F04_composite.png"), composite_final,
       width = COMP_W, height = COMP_H, units = "mm", dpi = 300)

message("F04 composite (5-panel, 3-column layout) saved")

# --- Supplementary Excel: one workbook, sheets keyed to figure panels ---
source("04_Figures/shared/figure_supplement_helpers.R")

cat("=== F04 supplementary workbook ===\n")
enrichment_blunting_df <- read.csv("04_Figures/F04/c_data/panel_supp/enrichment_blunting.csv",
                                   stringsAsFactors = FALSE, check.names = FALSE)
f04_specs <- list(
  list(name = "panel_A_ora_quadrant",       path = "04_Figures/F04/c_data/panel_A/ora_quadrant.csv"),
  list(name = "panel_B_pattern_class",      path = "04_Figures/F04/c_data/panel_C_heatmap/pattern_classification.csv"),
  list(name = "panel_B_sankey",             path = "04_Figures/F04/c_data/panel_C_heatmap/sankey_links.csv"),
  list(name = "panel_B_bar",                path = "04_Figures/F04/c_data/panel_C_heatmap/bar_data.csv"),
  list(name = "panel_C_fry_results",        path = "04_Figures/F04/c_data/panel_D_fry/fry_results_all.csv"),
  list(name = "panel_C_fry_driving",        path = "04_Figures/F04/c_data/panel_D_fry/driving_proteins.csv"),
  list(name = "panel_D_nes_scatter",        path = "04_Figures/F04/c_data/panel_B/nes_scatter.csv"),
  list(name = "panel_E_rrho2_summary",      path = "04_Figures/F04/c_data/panel_E/rrho2_summary.csv"),
  list(name = "panel_E_rrho2_hotspot",      path = "04_Figures/F04/c_data/panel_E/rrho2_hotspot_genes.csv"),
  list(name = "panel_E_rrho2_ora_concord",  path = "04_Figures/F04/c_data/panel_E/rrho2_ora_concordant.csv"),
  list(name = "panel_E_rrho2_ora_discord",  path = "04_Figures/F04/c_data/panel_E/rrho2_ora_discordant.csv"),
  list(name = "SUPP_enrichment_blunting",   df = enrichment_blunting_df),
  list(name = "SUPP_ora_dedup",             path = "04_Figures/F04/c_data/panel_supp/SUPP_ora_dedup_sensitivity.csv"),
  list(name = "SUPP_rho_bootstrap",         path = "04_Figures/F04/c_data/panel_supp/SUPP_rho_bootstrap.csv"),
  list(name = "SUPP_threshold_sens",        path = "04_Figures/F04/c_data/panel_supp/SUPP_threshold_sensitivity.csv"),
  list(name = "SUPP_goslim_bars",           path = "04_Figures/F04/c_data/panel_supp/SUPP_goslim_distribution.csv"),
  list(name = "SUPP_fry_leading",           path = "04_Figures/F04/c_data/panel_supp/SUPP_fry_leading_edge.csv")
)
build_workbook(
  "04_Figures/F04/c_data/F04_supplementary.xlsx",
  title = "F04 \u2014 Figure 4 source data",
  description = "Training-concordance diagnostics: quadrant ORA, pathway NES scatter, per-protein pattern classification, fry rotation test, RRHO2.",
  overview_df = data.frame(
    Sheet = c(
      "panel_A_ora_quadrant",
      "panel_B_pattern_class", "panel_B_sankey", "panel_B_bar",
      "panel_C_fry_results", "panel_C_fry_driving",
      "panel_D_nes_scatter",
      "panel_E_rrho2_summary",
      "panel_E_rrho2_hotspot", "panel_E_rrho2_ora_concord", "panel_E_rrho2_ora_discord",
      "SUPP_enrichment_blunting",
      "SUPP_ora_dedup", "SUPP_rho_bootstrap", "SUPP_threshold_sens",
      "SUPP_goslim_bars", "SUPP_fry_leading"),
    Description = c(
      "Panel A: ORA by training-concordance scatter quadrant",
      "Panel B: per-protein concordance pattern classification (heatmap source)",
      "Panel B: pathway-protein sankey links (heatmap source)",
      "Panel B: per-pattern bar chart counts",
      "Panel C: fry rotation test p-values + directions per pathway",
      "Panel C: concordant driving proteins (gene, t-stats, logFC)",
      "Panel D: NES scatter (Training_Young vs Training_Old) per pathway, Spearman + Fisher Z",
      "Panel E: RRHO2 quadrant summary (max -log10p per quadrant)",
      "Panel E: RRHO2 hotspot genes per quadrant",
      "Panel E: ORA on RRHO2 concordant quadrant genes",
      "Panel E: ORA on RRHO2 discordant quadrant genes",
      "SUPP: blunting enrichment (Pi-score-weighted ORA on Training_Old DEPs)",
      "SUPP: ORA dedup sensitivity across Jaccard cutoffs",
      "SUPP: Spearman rho bootstrap (1000 reps, 95% CI)",
      "SUPP: Protein counts per quadrant across significance thresholds",
      "SUPP: GO Slim category distribution by concordance quadrant",
      "SUPP: Top 20 fry driving proteins by |t-stat| in Training Old"),
    stringsAsFactors = FALSE),
  sheet_specs = f04_specs
)
cleanup_after_workbook(f04_specs,
  extra_subdirs = c("04_Figures/F04/c_data/panel_A",
                     "04_Figures/F04/c_data/panel_B",
                     "04_Figures/F04/c_data/panel_C",
                     "04_Figures/F04/c_data/panel_C_heatmap",
                     "04_Figures/F04/c_data/panel_D_fry",
                     "04_Figures/F04/c_data/panel_E",
                     "04_Figures/F04/c_data/panel_supp",
                     "04_Figures/F04/c_data/supp"))
