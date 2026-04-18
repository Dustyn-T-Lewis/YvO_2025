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
r_pear_A     <- r_pear
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

COMP_W     <- 720     # mm (geometry pass: wider for panel B horizontal space)
COMP_H     <- 440     # mm (taller so Panel A's square scatter can grow to match Panel B's ~175mm plot panel)
TAG_SZ     <- composite_text_sizes(COMP_H)$tag   # F01-proportional

# -- Title/subtitle styling (F01-proportional to COMP_H) --------------------

TITLE_THEME <- theme(
  plot.title          = element_text(face = "bold", size = composite_text_sizes(COMP_H)$title,
                                     hjust = 0, margin = margin(t = 2, b = 1, unit = "mm")),
  plot.subtitle       = element_text(face = "bold.italic", size = composite_text_sizes(COMP_H)$subtitle,
                                     color = "grey40", hjust = 0,
                                     margin = margin(b = 1, unit = "mm")),
  plot.title.position = "panel",
  plot.margin         = margin(t = 6, r = 2, b = 2, l = 2)
)

# -- Clean and re-label panels ----------------------------------------------

# Panel A: ORA quadrant composite (L=104mm matches individual panel positioning)
pA_clean <- composite +
  plot_annotation(
    title    = "Quadrant ORA (Concordance)",
    subtitle = sprintf("N = %d | %d DEPs (\u03a0) | %d enriched (FDR) | r = %.2f",
                       n_total_A, n_sig_A, n_enrich_A, r_pear_A),
    theme    = TITLE_THEME + theme(
      plot.title    = element_text(margin = margin(t = 4, b = 1, l = 42.1, unit = "mm")),
      plot.subtitle = element_text(margin = margin(b = -3, l = 42.1, unit = "mm")),
      plot.margin   = margin(t = 10, r = 2, b = 2, l = 2)
    )
  )

# Panel B: NES scatter (square) — stats from panel B snapshot
pB_clean <- pB +
  labs(title    = "Pathway Concordance",
       subtitle = sprintf("%d/%d sig. pathways | \u03c1 = %.2f [%.2f, %.2f] | %.0f%% concordant",
                          n_sig_pw_B, n_pw_B, rho_B, rho_lo_B, rho_hi_B, pw_conc_B * 100)) +
  TITLE_THEME +
  theme(plot.margin = margin(t = 14, r = 2, b = 2, l = 2))

# Panel C: pattern heatmap (full-height right anchor, L=6.5mm matches individual panel)
pC_clean <- p +
  labs(title    = "Protein-to-Pathway Structure",
       subtitle = sprintf("%d proteins | %d pathways",
                          n_total, n_pw)) +
  TITLE_THEME +
  theme(plot.title    = element_text(margin = margin(t = 2, b = 1, l = 18, unit = "mm")),
        plot.subtitle = element_text(margin = margin(b = 1, l = 18, unit = "mm")),
        plot.margin   = margin(t = 11, r = 4, b = 15, l = -18))

# Panel D: fry barcode + flanking ORA bars (L=10mm matches individual panel)
pD_clean <- pD_fry +
  plot_annotation(
    title    = "fry: Training Concordance",
    subtitle = sprintf("n = %d | dupCor = %.3f",
                       n_all_D, cor_imp_D),
    theme    = TITLE_THEME + theme(
      plot.title    = element_text(margin = margin(t = 2, b = 1, l = 12, unit = "mm")),
      plot.subtitle = element_text(margin = margin(b = 1, l = 12, unit = "mm")),
      plot.margin   = margin(t = 14, r = 2, b = 2, l = 2)
    )
  )

# Panel E: RRHO2 heatmap (square)
pE_clean <- pE_heat +
  labs(title    = "RRHO2 Concordance",
       subtitle = sprintf("%d shared genes | min P < 1e-%d (%d genes)",
                          n_shared_E, round(max_UU_E), n_UU_E)) +
  TITLE_THEME +
  theme(plot.margin = margin(t = 14, r = -0.5, b = 20, l = 3.5))

# -- Compose: 3-column layout with C spanning both rows ---------------------
# Design string: A spans cols 1-2, B in col 3, C spans cols 4-5 (both rows)
#                D spans cols 1-2, E in col 3
# Patchwork assigns plots to design characters in ALPHABETICAL order: A, B, C, D, E

layout <- "AABCC\nDDECC"

fig <- wrap_elements(full = pA_clean) +  # A: top-left (ORA)
       wrap_elements(full = pB_clean) +  # B: top-middle (NES scatter, square)
       wrap_elements(full = pC_clean) +  # C: right, full-height (pattern heatmap)
       wrap_elements(full = pD_clean) +  # D: bottom-left (fry)
       wrap_elements(full = pE_clean) +  # E: bottom-middle (RRHO2, square)
       plot_layout(design = layout,
                   widths  = c(4.85, 4.85, 6.5, 4.4, 4.4),
                   heights = c(1, 1))

# -- Panel tags via cowplot --------------------------------------------------
# x positions from widths: cols 1-2 = 10/25, col 3 starts at 10/25 = 0.400, cols 4-5 start at 16/25 = 0.640

composite_final <- ggdraw(fig) +
  draw_label("A", x = 0.005, y = 0.995 + 1/COMP_H,
             size = TAG_SZ, fontface = "bold", hjust = 0, vjust = 1) +
  draw_label("B", x = 0.388, y = 0.995 + 1/COMP_H,
             size = TAG_SZ, fontface = "bold", hjust = 0, vjust = 1) +
  draw_label("C", x = 0.6492, y = 0.995 + 1/COMP_H,
             size = TAG_SZ, fontface = "bold", hjust = 0, vjust = 1) +
  draw_label("D", x = 0.005, y = 0.495 + 1/COMP_H,
             size = TAG_SZ, fontface = "bold", hjust = 0, vjust = 1) +
  draw_label("E", x = 0.388, y = 0.495 + 1/COMP_H,
             size = TAG_SZ, fontface = "bold", hjust = 0, vjust = 1)

# -- Save -------------------------------------------------------------------

ggsave(file.path(RPT_PDF, "MAIN_F04_composite.pdf"), composite_final,
       width = COMP_W, height = COMP_H, units = "mm", device = pdf_device)
ggsave(file.path(RPT_PNG, "MAIN_F04_composite.png"), composite_final,
       width = COMP_W, height = COMP_H, units = "mm", dpi = 300)

message("F04 composite (5-panel, 3-column layout) saved")

# --- Supplementary Excel: one workbook, sheets keyed to figure panels ---
source("04_Figures/shared/figure_supplement_helpers.R")

cat("=== F04 supplementary workbook ===\n")
f04_specs <- list(
  list(name = "panel_A_ora_quadrant",       path = "04_Figures/F04/c_data/panel_A/ora_quadrant.csv"),
  list(name = "panel_B_nes_scatter",        path = "04_Figures/F04/c_data/panel_B/nes_scatter.csv"),
  list(name = "panel_C_pattern_class",      path = "04_Figures/F04/c_data/panel_C_heatmap/pattern_classification.csv"),
  list(name = "panel_C_sankey",             path = "04_Figures/F04/c_data/panel_C_heatmap/sankey_links.csv"),
  list(name = "panel_C_bar",                path = "04_Figures/F04/c_data/panel_C_heatmap/bar_data.csv"),
  list(name = "panel_C_pie",                path = "04_Figures/F04/c_data/panel_C_heatmap/pie_data.csv"),
  list(name = "panel_C_circos",             path = "04_Figures/F04/c_data/panel_C_heatmap/circos_links.csv"),
  list(name = "panel_D_fry_results",        path = "04_Figures/F04/c_data/panel_D_fry/fry_results_all.csv"),
  list(name = "panel_D_fry_driving",        path = "04_Figures/F04/c_data/panel_D_fry/driving_proteins.csv"),
  list(name = "panel_D_pi_p_overlap",       path = "04_Figures/F04/c_data/panel_D_fry/pi_p_overlap.csv"),
  list(name = "panel_E_rrho2_summary",      path = "04_Figures/F04/c_data/panel_E/rrho2_summary.csv"),
  list(name = "panel_E_rrho2_hotspot",      path = "04_Figures/F04/c_data/panel_E/rrho2_hotspot_genes.csv"),
  list(name = "panel_E_rrho2_ora_conc",     path = "04_Figures/F04/c_data/panel_E/rrho2_ora_concordant.csv"),
  list(name = "panel_E_rrho2_ora_disc",     path = "04_Figures/F04/c_data/panel_E/rrho2_ora_discordant.csv"),
  list(name = "panel_E_chord_combined",     path = "04_Figures/F04/c_data/panel_E/chord_training_combined.csv"),
  list(name = "panel_E_chord_pathways",     path = "04_Figures/F04/c_data/panel_E/chord_training_pathways.csv"),
  list(name = "SUPP_dir_asym_ora",          path = "04_Figures/F04/c_data/supp/h_directional_asymmetry_ora.csv"),
  list(name = "SUPP_dir_asym_tests",        path = "04_Figures/F04/c_data/supp/h_directional_asymmetry_tests.csv"),
  list(name = "SUPP_nes_scatter_collapsed", path = "04_Figures/F04/c_data/supp/nes_scatter_collapsed.csv"),
  list(name = "SUPP_enrichment_blunting",   path = "04_Figures/F04/c_data/panel_supp/enrichment_blunting.csv"),
  list(name = "data_dictionary",            path = "04_Figures/F04/c_data/00_data_dictionary.csv")
)
build_workbook(
  "04_Figures/F04/c_data/F04_supplementary.xlsx",
  title = "F04 \u2014 Figure 4 source data",
  description = "Training-concordance diagnostics: quadrant ORA, pathway NES scatter, per-protein pattern classification, fry rotation test, RRHO2.",
  overview_df = data.frame(
    Sheet = c(
      "panel_A_ora_quadrant",
      "panel_B_nes_scatter",
      "panel_C_pattern_class", "panel_C_sankey", "panel_C_bar", "panel_C_pie", "panel_C_circos",
      "panel_D_fry_results", "panel_D_fry_driving", "panel_D_pi_p_overlap",
      "panel_E_rrho2_summary", "panel_E_rrho2_hotspot",
      "panel_E_rrho2_ora_conc", "panel_E_rrho2_ora_disc",
      "panel_E_chord_combined", "panel_E_chord_pathways",
      "SUPP_dir_asym_ora", "SUPP_dir_asym_tests",
      "SUPP_nes_scatter_collapsed",
      "SUPP_enrichment_blunting",
      "data_dictionary"),
    Description = c(
      "Panel A: ORA by training-concordance scatter quadrant",
      "Panel B: NES scatter (Training_Young vs Training_Old) per pathway, Spearman + Fisher Z",
      "Panel C: per-protein concordance pattern classification (heatmap source)",
      "Panel C: pathway-protein sankey links (heatmap source)",
      "Panel C: per-pattern bar chart counts",
      "Panel C: per-pattern pie chart fractions",
      "Panel C: circos-style protein-pathway links",
      "Panel D: fry rotation test p-values + directions per pathway",
      "Panel D: concordant driving proteins (gene, t-stats, logFC)",
      "Panel D: Pi-score x raw-p overlap (per-protein significance audit)",
      "Panel E: RRHO2 quadrant summary (max -log10p + hotspot counts)",
      "Panel E: RRHO2 hotspot gene lists per quadrant",
      "Panel E: ORA on concordant RRHO2 quadrants",
      "Panel E: ORA on discordant RRHO2 quadrants",
      "Panel E: chord protein-pathway combined edges",
      "Panel E: chord pathway-level edges (NES-weighted)",
      "SUPP: directional asymmetry ORA (Young vs Old training response sign-flip)",
      "SUPP: directional asymmetry tests (binomial + Wilcoxon, BH-corrected)",
      "SUPP: NES scatter collapsed (representative pathway per database, deduplicated)",
      "SUPP: blunting enrichment (Pi-score-weighted ORA on Training_Old DEPs)",
      "Column definitions for each sheet"),
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
