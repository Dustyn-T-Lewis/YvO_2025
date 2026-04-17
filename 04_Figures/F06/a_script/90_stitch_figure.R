# Figure 6 — Composite: WGCNA Module-Trait Associations
#
# Layout:
#   Top row:    Panel A — combined 18-column heatmap (contrasts + phenotype)
#   Bottom row: Panel B (concordance scatter) | Panel C (reversal scatter)
#
# Generates: MAIN_F06_composite.pdf/png

setwd(rprojroot::find_rstudio_root_file())
source("04_Figures/shared/style.R")

library(patchwork)
library(cowplot)
library(png)
library(grid)

RPT      <- "04_Figures/F06/b_reports"
RPT_PNLS <- "04_Figures/F06/b_reports/panels"
WRITING_DIR <- "/Users/dtl0018/Library/CloudStorage/OneDrive-AuburnUniversity/YvO_writing/Figures/F06"
BOX_DIR     <- "/Users/dtl0018/Library/CloudStorage/Box-Box/YvO_proteomics_manuscript_2026-04-16/02_Figures/F06_WGCNA"
dir.create(WRITING_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(BOX_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(BOX_DIR, "pdf"), recursive = TRUE, showWarnings = FALSE)

read_panel <- function(file, dir = RPT_PNLS) {
  path <- file.path(dir, file)
  if (!file.exists(path)) stop("Missing: ", path)
  rasterGrob(readPNG(path), interpolate = TRUE)
}

pA       <- read_panel("MAIN_panel_A_heatmap.png")
pB       <- read_panel("MAIN_panel_B_concordance.png")
pC       <- read_panel("MAIN_panel_C_reversal.png")
p_BC_leg <- read_panel("MAIN_panel_BC_size_legend.png")

# Layout: A on top (full width), B+C side-by-side below.
# Bottom row reduced internal gutter so scatters sit visually paired.
composite <- (wrap_elements(pA) + theme(plot.margin = margin(0, 0, -12, 0, "mm"))) /
  (wrap_elements(pB) | wrap_elements(pC)) +
  plot_layout(heights = c(0.47, 0.53)) &
  theme(plot.margin = margin(0, 0, 0, 0))

COMP_W <- 520  # mm (matches heatmap width)
COMP_H <- 535  # +15mm headroom so the A tag can sit above the title without clipping
TAG_SZ <- composite_text_sizes(COMP_H)$tag   # F01-proportional

# Panel letter tags + shared size legend via cowplot.
# Heatmap row ends near y = 1 - 0.45 = 0.55 (scatter row starts at y = 0.55).
# B/C letters are pulled left into the new 14 mm gutter on the scatter panels.
composite_final <- ggdraw(composite) +
  # Shared "Proteins" size legend centered between B and C, nudged up above
  # the x-axis title band so it does not collide with panel axes.
  draw_grob(p_BC_leg,
            x = 0.40, y = 0.010, width = 0.20, height = 0.045) +
  draw_label("A", x = 0.012, y = 1.000,
             size = TAG_SZ, fontface = "bold", hjust = 0, vjust = 1) +
  draw_label("B", x = 0.010, y = 0.520,
             size = TAG_SZ, fontface = "bold", hjust = 0, vjust = 1) +
  draw_label("C", x = 0.504, y = 0.520,
             size = TAG_SZ, fontface = "bold", hjust = 0, vjust = 1)

pdf_device <- get_pdf_device()

ggsave(file.path(RPT, "MAIN_F06_composite.pdf"), composite_final,
       width = COMP_W, height = COMP_H, units = "mm",
       device = pdf_device, limitsize = FALSE)
ggsave(file.path(RPT, "MAIN_F06_composite.png"), composite_final,
       width = COMP_W, height = COMP_H, units = "mm",
       dpi = 300, limitsize = FALSE)
file.copy(file.path(RPT, "MAIN_F06_composite.pdf"), file.path(WRITING_DIR, "01_main_composite_f06.pdf"), overwrite = TRUE)
file.copy(file.path(RPT, "MAIN_F06_composite.png"), file.path(WRITING_DIR, "01_main_composite_f06.png"), overwrite = TRUE)
file.copy(file.path(RPT, "MAIN_F06_composite.png"), file.path(BOX_DIR, "F06_main.png"), overwrite = TRUE)
file.copy(file.path(RPT, "MAIN_F06_composite.pdf"), file.path(BOX_DIR, "pdf", "F06_main.pdf"), overwrite = TRUE)

message("F06 composite saved: ", file.path(RPT, "MAIN_F06_composite.{pdf,png}"))

# --- Supplementary Excel: WGCNA module-trait analysis + per-panel data ---
source("04_Figures/shared/figure_supplement_helpers.R")

f06 <- function(p) file.path("04_Figures/F06/c_data", p)

cat("=== F06 supplementary workbook ===\n")

# Convert F07-consumed matrix RDS files to Excel-embeddable data frames (2026-04-15).
# Keeps these objects in the workbook so F07 can read them via readxl,
# and the RDS files can be deleted from F06/c_data.
.shared       <- readRDS(f06("shared_objects.rds"))
.MEs          <- matrix_to_df(as.matrix(readRDS(f06("MEs.rds"))),      "sample_id")
.me_pre       <- matrix_to_df(as.matrix(readRDS(f06("me_pre.rds"))),   "subject_key")
.me_post      <- matrix_to_df(as.matrix(readRDS(f06("me_post.rds"))),  "subject_key")
.delta_me     <- matrix_to_df(as.matrix(readRDS(f06("delta_me.rds"))), "subject_key")
.common_subj  <- data.frame(subject_key = .shared$common_subj, stringsAsFactors = FALSE)

f06_specs <- list(
  # --- Panel/supp outputs written by active panel + supp scripts ---
  list(name="panel_A_heatmap",             path=f06("01_panel_A_heatmap_data.csv")),
  list(name="panel_A_correlation_CIs",     path=f06("01_panel_A_correlation_CIs.csv")),
  list(name="panel_B_module_fgsea",        path=f06("panel_B_module_fgsea.csv")),
  list(name="panel_B_triptych_zscores",    path=f06("03_panel_B_heatmap_zscores.csv")),
  list(name="panel_B_triptych_eigengene",  path=f06("03_panel_B_eigengene_data.csv")),
  list(name="panel_B_triptych_enrichment", path=f06("03_panel_B_triptych_enrichment.csv")),
  list(name="panel_D_hub_network",         path=f06("04_panel_D_hub_network.csv")),
  list(name="panel_D_hub_CIs",             path=f06("04_panel_D_hub_CIs.csv")),
  list(name="panel_E_preservation",        path=f06("05_panel_E_preservation.csv")),
  list(name="SUPP_a01_dendrogram",         path=f06("asupp_B_QC_dendrogram_SUPP_data.csv")),
  list(name="SUPP_a02_bicor_sensitivity",  path=f06("supp/a02_bicor_sensitivity.csv")),
  list(name="SUPP_a03_compartment_enrichment", path=f06("supp/a03_compartment_enrichment.csv")),
  list(name="SUPP_a05_sft_fit_indices",    path=f06("supp/a05_sft_fit_indices.csv")),
  list(name="SUPP_e02_dep_module_counts",  path=f06("supp/e02_dep_module_counts.csv")),
  list(name="SUPP_e02_dep_module_proportions", path=f06("supp/e02_dep_module_proportions.csv")),
  list(name="SUPP_lmm_stratified",         path=f06("supp_F_lmm_stratified_audit.csv")),
  # --- WGCNA reference outputs from YvO_WGCNA_run.R (wgcna/ subdir) ---
  list(name="WGCNA_module_assignments",    path=f06("wgcna/wgcna_module_assignments.csv")),
  list(name="WGCNA_mod_bio_labels",        path=f06("mod_bio_labels.csv")),
  list(name="WGCNA_hub_proteins",          path=f06("wgcna/wgcna_hub_proteins.csv")),
  list(name="WGCNA_module_enrichment",     path=f06("wgcna/wgcna_module_enrichment.csv")),
  list(name="WGCNA_module_GO_enrichment",  path=f06("wgcna/wgcna_module_GO_enrichment.csv")),
  list(name="WGCNA_module_trait_cor",      path=f06("wgcna/wgcna_module_trait_correlations.csv")),
  list(name="WGCNA_module_trait_pval_bh",  path=f06("wgcna/wgcna_module_trait_pvalues_bh.csv")),
  list(name="WGCNA_baseline_trait_cor",         path=f06("wgcna/wgcna_baseline_trait_correlations.csv")),
  list(name="WGCNA_baseline_trait_cor_young",   path=f06("wgcna/wgcna_baseline_trait_correlations_young.csv")),
  list(name="WGCNA_baseline_trait_cor_old",     path=f06("wgcna/wgcna_baseline_trait_correlations_old.csv")),
  list(name="WGCNA_baseline_pval_bh",         path=f06("wgcna/wgcna_baseline_trait_pvalues_bh.csv")),
  list(name="WGCNA_baseline_pval_bh_young",   path=f06("wgcna/wgcna_baseline_trait_pvalues_bh_young.csv")),
  list(name="WGCNA_baseline_pval_bh_old",     path=f06("wgcna/wgcna_baseline_trait_pvalues_bh_old.csv")),
  list(name="WGCNA_baseline_pval_raw_young",  path=f06("wgcna/wgcna_baseline_trait_pvalues_raw_young.csv")),
  list(name="WGCNA_baseline_pval_raw_old",    path=f06("wgcna/wgcna_baseline_trait_pvalues_raw_old.csv")),
  list(name="WGCNA_change_trait_cor",         path=f06("wgcna/wgcna_change_trait_correlations.csv")),
  list(name="WGCNA_change_trait_cor_young",   path=f06("wgcna/wgcna_change_trait_correlations_young.csv")),
  list(name="WGCNA_change_trait_cor_old",     path=f06("wgcna/wgcna_change_trait_correlations_old.csv")),
  list(name="WGCNA_change_pval_bh",           path=f06("wgcna/wgcna_change_trait_pvalues_bh.csv")),
  list(name="WGCNA_change_pval_bh_young",     path=f06("wgcna/wgcna_change_trait_pvalues_bh_young.csv")),
  list(name="WGCNA_change_pval_bh_old",       path=f06("wgcna/wgcna_change_trait_pvalues_bh_old.csv")),
  list(name="WGCNA_change_pval_raw_young",    path=f06("wgcna/wgcna_change_trait_pvalues_raw_young.csv")),
  list(name="WGCNA_change_pval_raw_old",      path=f06("wgcna/wgcna_change_trait_pvalues_raw_old.csv")),
  list(name="WGCNA_sft_summary",           path=f06("wgcna/wgcna_sft_summary.csv")),
  list(name="WGCNA_lmm_contrast_audit",    path=f06("wgcna/wgcna_lmm_contrast_audit.csv")),
  list(name="WGCNA_lmm_stratified_audit",  path=f06("wgcna/wgcna_lmm_stratified_audit.csv")),
  list(name="WGCNA_gs_phenotype_choices",  path=f06("wgcna/gs_phenotype_choices.csv")),
  # --- Metadata (module_df omitted — duplicates WGCNA_module_assignments) ---
  list(name="metadata_samples",            path=f06("meta.csv")),
  list(name="metadata_subj_age",           path=f06("subj_age.csv")),
  list(name="metadata_pheno_wide",         path=f06("pheno_wide.csv")),
  list(name="metadata_imp_annotations",    path=f06("imp_annotations.csv")),
  # --- ME matrices (from RDS → data.frame for Excel) ---
  list(name="MEs",         df=.MEs),
  list(name="me_pre",      df=.me_pre),
  list(name="me_post",     df=.me_post),
  list(name="delta_me",    df=.delta_me),
  list(name="common_subj", df=.common_subj)
)

build_workbook(
  f06("F06_supplementary.xlsx"),
  title = "F06 \u2014 Figure 6 source data",
  description = "WGCNA module assignments, hub proteins, pathway enrichment, trait correlations, eigengene metadata (MEs, me_pre, me_post, \u0394me).",
  overview_df = data.frame(
    Sheet = c(
      # Panel + supp outputs (active pipeline)
      "panel_A_heatmap", "panel_A_correlation_CIs",
      "panel_B_module_fgsea",
      "panel_B_triptych_zscores", "panel_B_triptych_eigengene",
      "panel_B_triptych_enrichment",
      "panel_D_hub_network", "panel_D_hub_CIs",
      "panel_E_preservation",
      "SUPP_a01_dendrogram", "SUPP_a02_bicor_sensitivity",
      "SUPP_a03_compartment_enrichment", "SUPP_a05_sft_fit_indices",
      "SUPP_e02_dep_module_counts", "SUPP_e02_dep_module_proportions",
      "SUPP_lmm_stratified",
      # WGCNA reference outputs (from YvO_WGCNA_run.R)
      "WGCNA_module_assignments", "WGCNA_mod_bio_labels", "WGCNA_hub_proteins",
      "WGCNA_module_enrichment", "WGCNA_module_GO_enrichment",
      "WGCNA_module_trait_cor", "WGCNA_module_trait_pval_bh",
      "WGCNA_baseline_trait_cor", "WGCNA_baseline_trait_cor_young", "WGCNA_baseline_trait_cor_old",
      "WGCNA_baseline_pval_bh", "WGCNA_baseline_pval_bh_young", "WGCNA_baseline_pval_bh_old",
      "WGCNA_baseline_pval_raw_young", "WGCNA_baseline_pval_raw_old",
      "WGCNA_change_trait_cor", "WGCNA_change_trait_cor_young", "WGCNA_change_trait_cor_old",
      "WGCNA_change_pval_bh", "WGCNA_change_pval_bh_young", "WGCNA_change_pval_bh_old",
      "WGCNA_change_pval_raw_young", "WGCNA_change_pval_raw_old",
      "WGCNA_sft_summary", "WGCNA_lmm_contrast_audit", "WGCNA_lmm_stratified_audit",
      "WGCNA_gs_phenotype_choices",
      # Metadata
      "metadata_samples", "metadata_subj_age",
      "metadata_pheno_wide", "metadata_imp_annotations",
      # ME matrices
      "MEs", "me_pre", "me_post", "delta_me", "common_subj"),
    Description = c(
      "Panel A: 18-column heatmap data (modules x traits)",
      "Panel A: bootstrap CIs for heatmap correlations",
      "Panel B: per-module fGSEA collapsed pathway list (NES scatter input)",
      "SUPP triptych: heatmap z-scores for key modules",
      "SUPP triptych: per-module ME eigengene long-format (Pre vs Post by age)",
      "SUPP triptych: per-module enrichment results",
      "SUPP D: hub protein network edges per key module",
      "SUPP D: bootstrap CIs for hub kME",
      "SUPP E: module preservation Zsummary scores (a04 output)",
      "SUPP a01: sample/module dendrogram for QC composite",
      "SUPP a02: bicor vs Pearson sensitivity heatmap data",
      "SUPP a03: subcellular compartment enrichment per module",
      "SUPP a05: SFT fit indices per power",
      "SUPP e02: DEP-module counts per contrast",
      "SUPP e02: DEP-module proportions",
      "SUPP F: LMM stratified audit (legacy preservation)",
      "WGCNA: protein -> module color + gene symbol",
      "WGCNA: module ID, color, biological label, n_proteins",
      "WGCNA: top-10 hub proteins per module by kME",
      "WGCNA: multi-database ORA per module (H+KEGG+Reactome+GO:BP)",
      "WGCNA: GO:BP enrichment per module",
      "WGCNA: combined Pearson r (modules x traits)",
      "WGCNA: BH-corrected p-values (modules x traits)",
      "WGCNA: baseline (Pre) trait correlations - combined cohort",
      "WGCNA: baseline (Pre) trait correlations - Young only",
      "WGCNA: baseline (Pre) trait correlations - Old only",
      "WGCNA: baseline BH-corrected p-values - combined",
      "WGCNA: baseline BH-corrected p-values - Young only",
      "WGCNA: baseline BH-corrected p-values - Old only",
      "WGCNA: baseline raw p-values - Young only",
      "WGCNA: baseline raw p-values - Old only",
      "WGCNA: change (Post-Pre) trait correlations - combined cohort",
      "WGCNA: change trait correlations - Young only",
      "WGCNA: change trait correlations - Old only",
      "WGCNA: change BH-corrected p-values - combined",
      "WGCNA: change BH-corrected p-values - Young only",
      "WGCNA: change BH-corrected p-values - Old only",
      "WGCNA: change raw p-values - Young only",
      "WGCNA: change raw p-values - Old only",
      "WGCNA: soft-threshold power selection (R^2, connectivity)",
      "WGCNA: LMM contrast audit (lmer + emmeans + KR df, BH)",
      "WGCNA: LMM age-stratified audit",
      "WGCNA: gene-significance phenotype choices (panel_D_hub config)",
      "Metadata: sample-level (Col_ID, Group, Timepoint, Subject_ID)",
      "Metadata: subject -> age group mapping",
      "Metadata: per-subject phenotype data (baseline + change for all traits)",
      "Metadata: imputed-matrix protein annotations (uniprot -> gene + description)",
      "Module eigengenes per sample (rows = sample_id, cols = MEturquoise ... MEgrey)",
      "Pre-timepoint module eigengenes per subject (for F07 readers)",
      "Post-timepoint module eigengenes per subject (for F07 readers)",
      "Change in module eigengenes (Post - Pre) per subject (for F07 readers)",
      "Common subject keys across Pre/Post for subject-paired analyses"),
    stringsAsFactors = FALSE),
  sheet_specs = f06_specs
)
file.copy(f06("F06_supplementary.xlsx"), file.path(WRITING_DIR, "04_supp_data_f06.xlsx"), overwrite = TRUE)
file.copy(f06("F06_supplementary.xlsx"), file.path(BOX_DIR, "F06_source_data.xlsx"), overwrite = TRUE)
# F06 keeps its wgcna/ subdir + a few root CSVs as an internal intra-run pipeline
# cache consumed by ~20 F06 panel/supp scripts (a01_dendrogram, a06_module_trait,
# panel_D_hub, etc.). These read from CSVs during a single pipeline pass; the
# stitcher absorbs them into F06_supplementary.xlsx for cross-figure readers
# (F07, manuscript) at the end of each run. The CSVs persist between runs for
# developer ergonomics (partial stitcher reruns stay idempotent without
# re-fitting the ~15-min WGCNA network). F06 is structurally different from
# F01-F05/F07 because it has a heavy upstream producer (YvO_WGCNA_run.R) —
# see 04_Figures/shared/README.md for the full rationale.
cleanup_after_workbook(f06_specs,
  extra_subdirs = character(0),
  preserve_patterns = c(
    "^00_input/", "^01_normalization/", "^02_Imputation/", "^03_DEP/",
    "^04_Figures/shared/",
    "^04_Figures/F06/c_data/wgcna/",
    "^04_Figures/F06/c_data/mod_bio_labels\\.csv$",
    "^04_Figures/F06/c_data/meta\\.csv$",
    "^04_Figures/F06/c_data/imp_annotations\\.csv$",
    "^04_Figures/F06/c_data/subj_age\\.csv$",
    "^04_Figures/F06/c_data/pheno_wide\\.csv$"))
