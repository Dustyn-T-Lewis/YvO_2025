# F03 — Composite: 4 volcano rings (Aging, Training_Young, Training_Old, Interaction)
# Sources 4 per-contrast panels and composes into a 2x2 grid + NES gradient legend.
#
# Outputs:
#   b_reports/main/{pdf,png}/MAIN_F03_composite.{pdf,png}
#   c_data/F03_supplementary.xlsx

setwd(rprojroot::find_rstudio_root_file())
source("04_Figures/shared/style.R")
source("04_Figures/shared/volcano_ring.R")

suppressPackageStartupMessages({
  library(patchwork)
  library(cowplot)
  library(ggplot2)
})

source("04_Figures/F03/a_script/main/panels/panel_A_aging.R")           # -> pA
source("04_Figures/F03/a_script/main/panels/panel_B_training_young.R")  # -> pB
source("04_Figures/F03/a_script/main/panels/panel_C_training_old.R")    # -> pC
source("04_Figures/F03/a_script/main/panels/panel_D_interaction.R")     # -> pD

RPT_PDF <- "04_Figures/F03/b_reports/main/pdf"
RPT_PNG <- "04_Figures/F03/b_reports/main/png"
pdf_device <- get_pdf_device()

# --- Per-contrast stats for subtitles ---
dbs_used <- c("Hallmark", "GO Slim", "GO:BP", "KEGG", "Reactome")
contrast_stats <- function(ctr) {
  pi_col    <- paste0("pi_score_", ctr)
  n_dep     <- if (pi_col %in% names(dep_df)) sum(dep_df[[pi_col]] < 0.05, na.rm = TRUE) else 0
  ctr_rows  <- fgsea_all[fgsea_all$contrast == ctr & fgsea_all$database %in% dbs_used, ]
  n_total   <- sum(!is.na(ctr_rows$padj))
  n_sig     <- sum(ctr_rows$padj < 0.05, na.rm = TRUE)
  sprintf("%d DEPs (\u03a0 < 0.05)  |  %d / %d pathways (FDR < 0.05)", n_dep, n_sig, n_total)
}
sub_A <- contrast_stats("Aging")
sub_B <- contrast_stats("Training_Young")
sub_C <- contrast_stats("Training_Old")
sub_D <- contrast_stats("Interaction")

# --- NES gradient bar legend (overlaid at bottom, not in patchwork) ---
nes_legend <- build_nes_legend_bar(text_size = 5, title_size = 5,
                                   bar_margin = margin(0, 0, 0, 0, "mm"))

# --- 2x2 composite (legend placed as overlay below) ---
# Panels arrive pre-stripped (no title, subtitle, tag, legend) via strip_for_composite().
# Tight vertical margins: top-row pushes down, bottom-row pushes up → close the gap.
top_row <- (pA | pB) & theme(plot.margin = margin(8, 0, 0, 0, "mm"))
bot_row <- (pC | pD) & theme(plot.margin = margin(0, 0, 2, 0, "mm"))
composite <- (top_row / bot_row) +
  plot_layout(heights = c(1, 1))

COMP_W <- 178          # J Physiol double-column
COMP_H <- 186          # 2×89 rings + tag strips + legend

TAG_SZ     <- composite_text_sizes(COMP_H)$tag + 4       # 12pt (visually match F02 on taller canvas)
TTL_SZ     <- composite_text_sizes(COMP_H)$title + 2    # ~9.3pt
SUB_SZ     <- composite_text_sizes(COMP_H)$subtitle + 2 # ~6.9pt
X_LEFT     <- 0.020
X_RIGHT    <- 0.510
Y_TOP      <- 0.960
Y_BOT      <- 0.505
X_TTL      <- 0.040     # title starts right of tag (F02 convention)
TAG_DY     <- -0.002    # raise tag to baseline-align with smaller title (F02 convention)
SUB_OFFSET <- 0.022     # standard title-to-subtitle gap (F01/F02 convention)

composite <- ggdraw(composite) +
  # Tags (baseline-raised via TAG_DY to align with title)
  draw_label("A", x = X_LEFT,  y = Y_TOP - TAG_DY, size = TAG_SZ, fontface = "bold", hjust = 0, vjust = 1) +
  draw_label("B", x = X_RIGHT, y = Y_TOP - TAG_DY, size = TAG_SZ, fontface = "bold", hjust = 0, vjust = 1) +
  draw_label("C", x = X_LEFT,  y = Y_BOT - TAG_DY, size = TAG_SZ, fontface = "bold", hjust = 0, vjust = 1) +
  draw_label("D", x = X_RIGHT, y = Y_BOT - TAG_DY, size = TAG_SZ, fontface = "bold", hjust = 0, vjust = 1) +
  # Titles
  draw_label("Aging Effect",                    x = X_LEFT  + X_TTL, y = Y_TOP, size = TTL_SZ, fontface = "bold", hjust = 0, vjust = 1) +
  draw_label("Training Response (Young)",       x = X_RIGHT + X_TTL, y = Y_TOP, size = TTL_SZ, fontface = "bold", hjust = 0, vjust = 1) +
  draw_label("Training Response (Old)",         x = X_LEFT  + X_TTL, y = Y_BOT, size = TTL_SZ, fontface = "bold", hjust = 0, vjust = 1) +
  draw_label("Age \u00d7 Training Interaction", x = X_RIGHT + X_TTL, y = Y_BOT, size = TTL_SZ, fontface = "bold", hjust = 0, vjust = 1) +
  # Subtitles (stats)
  draw_label(sub_A, x = X_LEFT  + X_TTL, y = Y_TOP - SUB_OFFSET, size = SUB_SZ, fontface = "bold.italic", colour = "grey40", hjust = 0, vjust = 1) +
  draw_label(sub_B, x = X_RIGHT + X_TTL, y = Y_TOP - SUB_OFFSET, size = SUB_SZ, fontface = "bold.italic", colour = "grey40", hjust = 0, vjust = 1) +
  draw_label(sub_C, x = X_LEFT  + X_TTL, y = Y_BOT - SUB_OFFSET, size = SUB_SZ, fontface = "bold.italic", colour = "grey40", hjust = 0, vjust = 1) +
  draw_label(sub_D, x = X_RIGHT + X_TTL, y = Y_BOT - SUB_OFFSET, size = SUB_SZ, fontface = "bold.italic", colour = "grey40", hjust = 0, vjust = 1) +
  # NES legend bar — thin overlay tucked under panels c/d
  draw_plot(nes_legend, x = 0.35, y = 0.025, width = 0.30, height = 0.025)

dir.create(RPT_PDF, recursive = TRUE, showWarnings = FALSE)
dir.create(RPT_PNG, recursive = TRUE, showWarnings = FALSE)
ggsave(file.path(RPT_PDF, "MAIN_F03_composite.pdf"), composite,
       width = COMP_W, height = COMP_H, units = "mm", device = pdf_device,
       limitsize = FALSE)
ggsave(file.path(RPT_PNG, "MAIN_F03_composite.png"), composite,
       width = COMP_W, height = COMP_H, units = "mm", dpi = 300,
       limitsize = FALSE)
message("F03 composite (2x2 volcano rings) saved")

# --- Supplementary Excel: per-contrast DEP tables + ring terms + supp diagnostics ---
source("04_Figures/shared/figure_supplement_helpers.R")

f03_contrasts <- c("Aging", "Training_Young", "Training_Old", "Interaction")
f03_specs <- lapply(f03_contrasts, function(ctr) {
  list(name = paste0("MAIN_", ctr),
       path = sprintf("03_DEP/c_data/04_per_contrast_results/%s.csv", ctr))
})
f03_specs <- c(f03_specs, lapply(c("A", "B", "C", "D"), function(tag) {
  list(name = paste0("RING_panel_", tag),
       path = sprintf("04_Figures/F03/c_data/panel_%s/ring_terms.csv", tag))
}))
f03_specs <- c(f03_specs, list(
  list(name = "SUPP_panel_B_phist",       path = "04_Figures/F03/c_data/supp/panel_B_phist.csv"),
  list(name = "SUPP_panel_C_ma",          path = "04_Figures/F03/c_data/supp/panel_C_ma.csv"),
  list(name = "SUPP_panel_D_sensitivity", path = "04_Figures/F03/c_data/supp/panel_D_sensitivity.csv"),
  list(name = "SUPP_panel_E_outlier",     path = "04_Figures/F03/c_data/supp/panel_E_outlier.csv"),
  list(name = "fGSEA_cache_used",         path = "04_Figures/shared/fgsea_tstat_all_v2.csv")
))

cat("=== F03 supplementary workbook ===\n")
build_workbook(
  "04_Figures/F03/c_data/F03_supplementary.xlsx",
  title = "F03 \u2014 Figure 3 source data",
  description = "Per-contrast volcano ring source data and supplementary diagnostics (p-value histograms, MA plots, imputation and outlier sensitivity).",
  overview_df = data.frame(
    Sheet = c(paste0("MAIN_", f03_contrasts),
              paste0("RING_panel_", c("A","B","C","D")),
              "SUPP_panel_B_phist", "SUPP_panel_C_ma",
              "SUPP_panel_D_sensitivity", "SUPP_panel_E_outlier",
              "fGSEA_cache_used"),
    Description = c(
      "Panel A volcano: Aging per-protein DEP table (logFC, CI, t, p, FDR, Pi)",
      "Panel B volcano: Training_Young per-protein DEP table",
      "Panel C volcano: Training_Old per-protein DEP table",
      "Panel D volcano: Age x Training Interaction per-protein DEP table",
      "Panel A ring terms: fGSEA top pathways shown on Aging volcano ring",
      "Panel B ring terms: fGSEA top pathways shown on Training_Young volcano ring",
      "Panel C ring terms: fGSEA top pathways shown on Training_Old volcano ring",
      "Panel D ring terms: fGSEA top pathways shown on Interaction volcano ring",
      "SUPP B: per-contrast raw p-values (source data for p-histogram)",
      "SUPP C: per-contrast MA-plot source data (logFC vs mean intensity, sig flag)",
      "SUPP D: imputation sensitivity Spearman rho per contrast",
      "SUPP E: outlier sensitivity DEP counts (full vs reduced cohort)",
      "fGSEA cache used by panel rings (4 contrasts x 4 databases, pathway-level)"),
    stringsAsFactors = FALSE),
  sheet_specs = f03_specs
)
cleanup_after_workbook(f03_specs,
  extra_subdirs = c("04_Figures/F03/c_data/panel_A",
                     "04_Figures/F03/c_data/panel_B",
                     "04_Figures/F03/c_data/panel_C",
                     "04_Figures/F03/c_data/panel_D",
                     "04_Figures/F03/c_data/supp"))
