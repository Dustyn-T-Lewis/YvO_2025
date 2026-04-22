# F05 — Composite: Aging Reversal (Panels A–E)
# 3-column geometry-aware layout (panels labeled in visual reading order):
#   Top row:    A (Quadrant ORA)  | B (NES scatter)  | C (Pattern heatmap, full-height)
#   Bottom row: D (fry barcode)   | E (RRHO2)         | C (continues)
# Outputs: MAIN_F05_composite.{pdf,png}
#
# Source order: A, B, D, E, then C. Panel C (pattern heatmap) loads
# AnnotationDbi (via go_slim_categories.R); the S4 select() masking is now
# repaired inside go_slim_categories.R itself (re-attaches dplyr after dispatch).
# C is still sourced last to preserve stat-snapshot variable ordering.

setwd(rprojroot::find_rstudio_root_file())

suppressPackageStartupMessages({
  library(patchwork)
  library(cowplot)
  library(ggplot2)
})

# -- Source panels (C last for stat-snapshot ordering) -----------------------

message("=== F05 Composite: sourcing panels ===")
source("04_Figures/F05/a_script/main/panels/panel_A_ORA.R")              # -> composite (patchwork)
# Snapshot Panel A stats before Panel C overwrites `n_total` and `sig_df`
n_total_A    <- nrow(scatter_df)
n_sig_A      <- n_sig
n_enrich_A   <- n_enrich
r_pear_A     <- r_pear
source("04_Figures/F05/a_script/main/panels/panel_B_nes_scatter.R")      # -> pB (ggplot)
# Snapshot Panel B stats
n_pw_B       <- nrow(fgsea_wide)
n_sig_pw_B   <- n_total_sig
rho_B        <- as.numeric(nes_cor_all$estimate)
rho_lo_B     <- nes_ci_all[1]
rho_hi_B     <- nes_ci_all[2]
pw_rev_B     <- pw_rev_frac
source("04_Figures/F05/a_script/main/panels/panel_D_fry.R")              # -> pD_fry, p1, p2, p_t
# Snapshot Panel D stats
cor_imp_D    <- cor_imp
n_all_D      <- n_all
circ_r_D     <- circ_r
source("04_Figures/F05/a_script/main/panels/panel_E_rrho2.R")            # -> pE_heat (ggplot)
# Snapshot Panel E stats — F05 strongest signal is the reversed quadrant (UD)
n_shared_E   <- n_shared
max_rev_E    <- max(max_UD, max_DU)
n_rev_E      <- if (max_UD >= max_DU) n_UD else n_DU
source("04_Figures/F05/a_script/main/panels/panel_C_pattern_heatmap.R")  # -> p (ggplot, theme_void)

RPT_PDF <- "04_Figures/F05/b_reports/main/pdf"
RPT_PNG <- "04_Figures/F05/b_reports/main/png"
dir.create(RPT_PDF, recursive = TRUE, showWarnings = FALSE)
dir.create(RPT_PNG, recursive = TRUE, showWarnings = FALSE)
pdf_device <- get_pdf_device()

# -- Layout constants (matched to F04) ---------------------------------------

COMP_W     <- 420     # mm (wider to fit equal-size D + E panels)
COMP_H     <- 310     # mm (taller for deeper bottom row)
PRINT_SCALE <- 380 / 178  # 380mm source → 178mm print = 2.13x
TAG_SZ     <- round(10 * PRINT_SCALE * 0.85)   # uniform tag size (~18pt)
TTL_SZ     <- round(10 * PRINT_SCALE * 0.85)   # uniform title size (~18pt)
SUB_SZ     <- round(7 * PRINT_SCALE * 0.85)    # uniform subtitle size (~13pt)

# -- Stat-snapshot title strings --------------------------------------------

ttl_A <- "Quadrant ORA (Reversal)"
sub_A <- sprintf("N = %d | %d DEPs (\u03a0) | %d enriched (FDR) | r = %.2f",
                 n_total_A, n_sig_A, n_enrich_A, r_pear_A)
ttl_B <- "Protein-to-Pathway"
sub_B <- sprintf("%d proteins | %d pathways", n_total, n_pw)
ttl_C <- "fry: Reversal"
sub_C <- sprintf("n = %d | dupCor = %.3f", n_all_D, cor_imp_D)
ttl_D <- "Pathway Reversal"
sub_D <- sprintf("\u03c1 = %.2f | %.0f%% reversed",
                 rho_B, pw_rev_B * 100)
ttl_E <- "RRHO2 Reversal"
sub_E <- sprintf("%d genes | max %d", n_shared_E, n_rev_E)

# -- Compose: 2-row layout with title spacers, A+B top-aligned (mirrors F04) -
# 16-row, 13-col grid. B extends 2 rows past A for extra height.
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
  "########BBBBBB",   # rows 8-9: B extends past A
  "########BBBBBB",
  "##############",   # row 10: narrow title spacer (c + d + e titles)
  "CCCCCDDDD#EEEE",   # rows 11-16: C=5, D=4, #=1 gap, E=4
  "CCCCCDDDD#EEEE",
  "CCCCCDDDD#EEEE",
  "CCCCCDDDD#EEEE",
  "CCCCCDDDD#EEEE",
  "CCCCCDDDD#EEEE",
  sep = "\n"
)

# Raise bottom panels independently via negative top margins
pD_fry  <- pD_fry + plot_annotation(theme = theme(plot.margin = margin(-9, 0, 0, 0, "mm")))
pB      <- pB + theme(plot.margin = margin(-25, -14, -5, 5, "mm"))
pE_heat <- pE_heat + theme(plot.margin = margin(-24, 0, 0, -10, "mm"))

fig <- wrap_elements(full = composite) +   # A: top-left (ORA, 8 cols)
       wrap_elements(full = p) +           # B: top-right (heatmap, 5 cols)
       wrap_elements(full = pD_fry) +      # C: bottom-left (fry, 5 cols)
       wrap_elements(full = pB) +          # D: bottom-mid (NES scatter, 4 cols)
       wrap_elements(full = pE_heat) +     # E: bottom-right (RRHO2, 4 cols)
       plot_layout(design = layout,
                   widths  = rep(1, 14),
                   heights = c(4, rep(10, 6), 5, 5, 0.5, rep(10, 6)))

# -- Panel tags, titles, subtitles via cowplot -------------------------------
# 13-col, 16-row grid. Heights: [4, 10×6, 10, 10, 3, 10×6] = 147 total
#   Top spacer:    y  1.000 → 0.973   (4/147)
#   A rows (6):    y  0.973 → 0.565   (60/147)
#   B extension:   y  0.565 → 0.429   (20/147)  — B only, A empty
#   Mid spacer:    y  0.429 → 0.408   (3/147)
#   Bottom panels: y  0.408 → 0.000   (60/147)
X_A <- 0.005;  X_B <- 0.539;  X_C <- 0.010;  X_D <- 0.362;  X_E <- 0.675
X_TTL      <- 0.030        # title/subtitle offset right of tag
TAG_DY     <- -0.002       # raise tag for baseline alignment
SUB_OFFSET <- 0.025        # subtitle gap below title

# Per-panel Y positions — a+b vertically aligned, d+e vertically aligned
Y_A <- 0.990;  Y_B <- 0.990           # a + b aligned at top
Y_C <- 0.509;  Y_D <- 0.507;  Y_E <- 0.512   # C down 1mm

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

# -- Save -------------------------------------------------------------------

ggsave(file.path(RPT_PDF, "MAIN_F05_composite.pdf"), composite_final,
       width = COMP_W, height = COMP_H, units = "mm", device = pdf_device)
ggsave(file.path(RPT_PNG, "MAIN_F05_composite.png"), composite_final,
       width = COMP_W, height = COMP_H, units = "mm", dpi = 300)

message("F05 composite (5-panel, 3-column layout) saved")

# --- Supplementary Excel: one workbook, sheets keyed to figure panels ---
source("04_Figures/shared/figure_supplement_helpers.R")

cat("=== F05 supplementary workbook ===\n")
f05_specs <- list(
  list(name = "panel_A_ora_quadrant",       path = "04_Figures/F05/c_data/panel_A/ora_quadrant.csv"),
  list(name = "panel_B_nes_scatter",        path = "04_Figures/F05/c_data/panel_B/nes_scatter.csv"),
  list(name = "panel_C_pattern_class",      path = "04_Figures/F05/c_data/panel_C_heatmap/pattern_classification.csv"),
  list(name = "panel_C_sankey",             path = "04_Figures/F05/c_data/panel_C_heatmap/sankey_links.csv"),
  list(name = "panel_C_bar",                path = "04_Figures/F05/c_data/panel_C_heatmap/bar_data.csv"),
  list(name = "panel_C_pie",                path = "04_Figures/F05/c_data/panel_C_heatmap/pie_data.csv"),
  list(name = "panel_C_circos",             path = "04_Figures/F05/c_data/panel_C_heatmap/circos_links.csv"),
  list(name = "panel_C_union_dep",          path = "04_Figures/F05/c_data/panel_C_heatmap/union_dep_list.csv"),
  list(name = "panel_C_union_ora",          path = "04_Figures/F05/c_data/panel_C_heatmap/union_ora_links.csv"),
  list(name = "panel_D_fry_results",        path = "04_Figures/F05/c_data/panel_D_fry/fry_results_all.csv"),
  list(name = "panel_D_fry_driving",        path = "04_Figures/F05/c_data/panel_D_fry/driving_proteins.csv"),
  list(name = "panel_D_pi_p_overlap",       path = "04_Figures/F05/c_data/panel_D_fry/pi_p_overlap.csv"),
  list(name = "panel_E_rrho2_summary",      path = "04_Figures/F05/c_data/panel_E/rrho2_summary.csv"),
  list(name = "panel_E_rrho2_hotspot",      path = "04_Figures/F05/c_data/panel_E/rrho2_hotspot_genes.csv"),
  list(name = "panel_E_rrho2_ora_conc",     path = "04_Figures/F05/c_data/panel_E/rrho2_ora_concordant.csv"),
  list(name = "panel_E_rrho2_ora_disc",     path = "04_Figures/F05/c_data/panel_E/rrho2_ora_discordant.csv"),
  list(name = "panel_E_chord_combined",     path = "04_Figures/F05/c_data/panel_E/chord_aging_combined.csv"),
  list(name = "panel_E_chord_pathways",     path = "04_Figures/F05/c_data/panel_E/chord_aging_pathways.csv"),
  list(name = "panel_E_rrho2_sensitivity",  path = "04_Figures/F05/c_data/panel_E/rrho2_sensitivity.csv"),
  list(name = "REVERSAL_melov_permutation",     path = "04_Figures/F05/c_data/reversal_tests/melov_permutation.csv"),
  list(name = "REVERSAL_contingency",           path = "04_Figures/F05/c_data/reversal_tests/reversal_contingency.csv"),
  list(name = "REVERSAL_threshold_sensitivity", path = "04_Figures/F05/c_data/reversal_tests/threshold_sensitivity.csv"),
  list(name = "REVERSAL_signed_score",          path = "04_Figures/F05/c_data/reversal_tests/signed_reversal_score.csv"),
  list(name = "SUPP_dir_asym_ora",          path = "04_Figures/F05/c_data/supp/h_directional_asymmetry_ora.csv"),
  list(name = "SUPP_dir_asym_tests",        path = "04_Figures/F05/c_data/supp/h_directional_asymmetry_tests.csv"),
  list(name = "SUPP_nes_scatter_collapsed", path = "04_Figures/F05/c_data/supp/nes_scatter_collapsed.csv"),
  list(name = "SUPP_enrichment_reversal",   path = "04_Figures/F05/c_data/panel_supp/enrichment_reversal.csv"),
  list(name = "SUPP_reversal_pathway_stats", path = "04_Figures/F05/c_data/panel_supp/reversal_pathway_stats.csv"),
  list(name = "SUPP_reversal_vs_effect",    path = "04_Figures/F05/c_data/supp/supp_reversal_vs_effect.csv"),
  list(name = "data_dictionary",            path = "04_Figures/F05/c_data/00_data_dictionary.csv")
)
build_workbook(
  "04_Figures/F05/c_data/F05_supplementary.xlsx",
  title = "F05 \u2014 Figure 5 source data",
  description = "Aging-reversal diagnostics: quadrant ORA, pathway NES scatter, per-protein pattern classification, fry rotation test, RRHO2.",
  overview_df = data.frame(
    Sheet = c(
      "panel_A_ora_quadrant",
      "panel_B_nes_scatter",
      "panel_C_pattern_class", "panel_C_sankey", "panel_C_bar", "panel_C_pie", "panel_C_circos",
      "panel_C_union_dep", "panel_C_union_ora",
      "panel_D_fry_results", "panel_D_fry_driving", "panel_D_pi_p_overlap",
      "panel_E_rrho2_summary", "panel_E_rrho2_hotspot",
      "panel_E_rrho2_ora_conc", "panel_E_rrho2_ora_disc",
      "panel_E_chord_combined", "panel_E_chord_pathways",
      "panel_E_rrho2_sensitivity",
      "REVERSAL_melov_permutation", "REVERSAL_contingency",
      "REVERSAL_threshold_sensitivity", "REVERSAL_signed_score",
      "SUPP_dir_asym_ora", "SUPP_dir_asym_tests",
      "SUPP_nes_scatter_collapsed",
      "SUPP_enrichment_reversal", "SUPP_reversal_pathway_stats",
      "SUPP_reversal_vs_effect",
      "data_dictionary"),
    Description = c(
      "Panel A: ORA by reversal-quadrant scatter (Reversed Up/Down, Exacerbated Up/Down)",
      "Panel B: NES scatter (Aging vs Training_Old) per pathway with Spearman + Fisher Z",
      "Panel C: per-protein reversal pattern classification (heatmap source)",
      "Panel C: pathway-protein sankey links (heatmap source)",
      "Panel C: per-pattern bar chart counts",
      "Panel C: per-pattern pie chart fractions",
      "Panel C: circos protein-pathway links",
      "Panel C: union DEP list (Pi-significant in Aging OR Training_Old)",
      "Panel C: union ORA links (pathway enrichment on union DEP set)",
      "Panel D: fry rotation test for reversal contrasts",
      "Panel D: reversal driving proteins (concordant sign-flip)",
      "Panel D: Pi-score x raw-p overlap audit",
      "Panel E: RRHO2 Aging vs Training_Old quadrant summary",
      "Panel E: RRHO2 hotspot gene lists per quadrant",
      "Panel E: ORA on concordant RRHO2 quadrants",
      "Panel E: ORA on discordant RRHO2 quadrants",
      "Panel E: chord protein-pathway combined edges",
      "Panel E: chord pathway-level edges (NES-weighted)",
      "Panel E: RRHO2 sensitivity (varying overlap thresholds)",
      "REVERSAL: Melov distance permutation test (10k perms, bootstrap CI)",
      "REVERSAL: 2x2 Fisher exact: aging dir x training dir + OR + counts",
      "REVERSAL: reversal % at |logFC| thresholds 0.1/0.2/0.3",
      "REVERSAL: signed reversal score (Pearson r +/- CI for logFC_Aging vs logFC_Training_Old)",
      "SUPP: directional asymmetry ORA per reversal quadrant",
      "SUPP: directional asymmetry tests (binomial + Wilcoxon, BH)",
      "SUPP: NES scatter collapsed (representative pathway per database)",
      "SUPP: pathway-level reversal enrichment",
      "SUPP: per-pathway reversal stats (NES Aging, NES TO, sign agreement)",
      "SUPP: reversal magnitude vs effect-size correlation",
      "Column definitions for each sheet"),
    stringsAsFactors = FALSE),
  sheet_specs = f05_specs
)
cleanup_after_workbook(f05_specs,
  extra_subdirs = c("04_Figures/F05/c_data/panel_A",
                     "04_Figures/F05/c_data/panel_B",
                     "04_Figures/F05/c_data/panel_C",
                     "04_Figures/F05/c_data/panel_C_heatmap",
                     "04_Figures/F05/c_data/panel_D_fry",
                     "04_Figures/F05/c_data/panel_E",
                     "04_Figures/F05/c_data/panel_F",
                     "04_Figures/F05/c_data/panel_F_fry",
                     "04_Figures/F05/c_data/panel_G",
                     "04_Figures/F05/c_data/reversal_tests",
                     "04_Figures/F05/c_data/panel_supp",
                     "04_Figures/F05/c_data/supp"))
