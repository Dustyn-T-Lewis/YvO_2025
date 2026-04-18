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

VW <- 190
RPT_PDF <- "04_Figures/F03/b_reports/main/pdf"
RPT_PNG <- "04_Figures/F03/b_reports/main/png"
pdf_device <- get_pdf_device()

# --- NES gradient bar legend (centered below bottom row, narrow band) ---
nes_legend <- build_nes_legend_bar(bar_margin = margin(-4, 160, 4, 160, "mm"))

# --- 2x2 composite + NES legend ---
# Strip per-subplot tags so the outer canvas can place uniform tags.
volcano_plots_nt <- list(pA, pB, pC, pD) |> lapply(function(p) p + labs(tag = NULL))

composite <- ((volcano_plots_nt[[1]] | volcano_plots_nt[[2]]) /
              (volcano_plots_nt[[3]] | volcano_plots_nt[[4]]) /
              wrap_elements(full = nes_legend)) +
  plot_layout(heights = c(0.4970, 0.4970, 0.0060)) &
  theme(plot.margin = margin(12, 0, 0, 0, "mm"))

COMP_W <- VW * 2       # 380mm
COMP_H <- 380          # +20mm headroom for top tag band

TAG_SZ  <- composite_text_sizes(380)$tag
X_LEFT  <- 0.050
X_RIGHT <- 0.505
Y_TOP   <- 0.965 - 4/380
Y_BOT   <- 0.5016 - 4/380

composite <- ggdraw(composite) +
  draw_label("A", x = X_LEFT,  y = Y_TOP, size = TAG_SZ, fontface = "bold", hjust = 0, vjust = 1) +
  draw_label("B", x = X_RIGHT, y = Y_TOP, size = TAG_SZ, fontface = "bold", hjust = 0, vjust = 1) +
  draw_label("C", x = X_LEFT,  y = Y_BOT, size = TAG_SZ, fontface = "bold", hjust = 0, vjust = 1) +
  draw_label("D", x = X_RIGHT, y = Y_BOT, size = TAG_SZ, fontface = "bold", hjust = 0, vjust = 1)

dir.create(RPT_PDF, recursive = TRUE, showWarnings = FALSE)
dir.create(RPT_PNG, recursive = TRUE, showWarnings = FALSE)
ggsave(file.path(RPT_PDF, "MAIN_F03_composite.pdf"), composite,
       width = COMP_W, height = COMP_H, units = "mm", device = pdf_device)
ggsave(file.path(RPT_PNG, "MAIN_F03_composite.png"), composite,
       width = COMP_W, height = COMP_H, units = "mm", dpi = 300)
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
