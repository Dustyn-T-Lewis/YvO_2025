#!/usr/bin/env Rscript
# Rebuild supplementary XLSX workbooks for F03-F06 without re-rendering figures.

setwd(rprojroot::find_rstudio_root_file())
source("04_Figures/shared/figure_supplement_helpers.R")

BOX <- "/Users/dtl0018/Library/CloudStorage/Box-Box/YvO_proteomics_manuscript"

copy_to_box <- function(local_xlsx, fig_tag, supp_name) {
  f_dir <- switch(fig_tag,
    F00 = "F00_pipeline_QC", F01 = "F01_phenotype",
    F02 = "F02_proteome_DEP", F03 = "F03_volcanoes",
    F04 = "F04_training_concordance", F05 = "F05_aging_reversal",
    F06 = "F06_WGCNA", F07 = "F07_phenotype_prediction")
  sd <- file.path(BOX, "02_Figures", f_dir, paste0(fig_tag, "_source_data.xlsx"))
  file.copy(local_xlsx, sd, overwrite = TRUE)
  if (!is.null(supp_name)) {
    st <- file.path(BOX, "03_Supplementary_Tables", supp_name)
    file.copy(local_xlsx, st, overwrite = TRUE)
  }
  cat(sprintf("  Copied %s -> Box (%s + %s)\n", fig_tag, f_dir, supp_name))
}

# F03
cat("=== Rebuilding F03 workbook ===\n")
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
  list(name = "SUPP_panel_E_outlier",     path = "04_Figures/F03/c_data/supp/panel_E_outlier.csv")
))
build_workbook(
  "04_Figures/F03/c_data/F03_supplementary.xlsx",
  title = "F03 \u2014 Figure 3 source data",
  description = "Per-contrast volcano ring source data and supplementary diagnostics (p-value histograms, MA plots, imputation and outlier sensitivity).",
  overview_df = data.frame(
    Sheet = vapply(f03_specs, `[[`, "", "name"),
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
      "SUPP E: outlier sensitivity DEP counts (full vs reduced cohort)"),
    stringsAsFactors = FALSE),
  sheet_specs = f03_specs
)
copy_to_box("04_Figures/F03/c_data/F03_supplementary.xlsx", "F03", "S09_F03_volcanoes.xlsx")

# F04
cat("=== Rebuilding F04 workbook ===\n")
f04_specs <- list(
  list(name = "panel_A_ora_quadrant",       path = "04_Figures/F04/c_data/panel_A/ora_quadrant.csv"),
  list(name = "panel_B_pattern_class",      path = "04_Figures/F04/c_data/panel_B_heatmap/pattern_classification.csv"),
  list(name = "panel_B_sankey",             path = "04_Figures/F04/c_data/panel_B_heatmap/sankey_links.csv"),
  list(name = "panel_B_bar",                path = "04_Figures/F04/c_data/panel_B_heatmap/bar_data.csv"),
  list(name = "panel_C_fry_results",        path = "04_Figures/F04/c_data/panel_C_fry/fry_results_all.csv"),
  list(name = "panel_C_fry_driving",        path = "04_Figures/F04/c_data/panel_C_fry/driving_proteins.csv"),
  list(name = "panel_D_nes_scatter",        path = "04_Figures/F04/c_data/panel_D/nes_scatter.csv"),
  list(name = "panel_E_rrho2_summary",      path = "04_Figures/F04/c_data/panel_E/rrho2_summary.csv"),
  list(name = "SUPP_enrichment_blunting",   path = "04_Figures/F04/c_data/panel_supp/enrichment_blunting.csv"),
  list(name = "SUPP_ora_dedup_sensitivity", path = "04_Figures/F04/c_data/panel_supp/SUPP_ora_dedup_sensitivity.csv"),
  list(name = "SUPP_rho_bootstrap",         path = "04_Figures/F04/c_data/panel_supp/SUPP_rho_bootstrap.csv"),
  list(name = "SUPP_threshold_sensitivity", path = "04_Figures/F04/c_data/panel_supp/SUPP_threshold_sensitivity.csv"),
  list(name = "SUPP_goslim_distribution",   path = "04_Figures/F04/c_data/panel_supp/SUPP_goslim_distribution.csv"),
  list(name = "SUPP_fry_leading_edge",      path = "04_Figures/F04/c_data/panel_supp/SUPP_fry_leading_edge.csv")
)
build_workbook(
  "04_Figures/F04/c_data/F04_supplementary.xlsx",
  title = "F04 \u2014 Figure 4 source data",
  description = "Training-concordance diagnostics: quadrant ORA, pathway NES scatter, per-protein pattern classification, fry rotation test, RRHO2, and supplementary sensitivity analyses.",
  overview_df = data.frame(
    Sheet = vapply(f04_specs, `[[`, "", "name"),
    Description = c(
      "Panel A: ORA by training-concordance scatter quadrant",
      "Panel B: per-protein concordance pattern classification (heatmap source)",
      "Panel B: pathway-protein sankey links (heatmap source)",
      "Panel B: per-pattern bar chart counts",
      "Panel C: fry rotation test p-values + directions per pathway",
      "Panel C: concordant driving proteins (gene, t-stats, logFC)",
      "Panel D: NES scatter (Training_Young vs Training_Old) per pathway, Spearman + Fisher Z",
      "Panel E: RRHO2 quadrant summary (max -log10p per quadrant)",
      "SUPP: blunting enrichment (Pi-score-weighted ORA on Training_Old DEPs)",
      "SUPP: ORA dedup Jaccard cutoff sensitivity (0.3, 0.5, 0.7, 1.0)",
      "SUPP: Spearman rho bootstrap (1000 resamples with 95% CI)",
      "SUPP: DEP threshold sensitivity (Pi, FDR, nominal p)",
      "SUPP: GO Slim category distribution by concordance quadrant",
      "SUPP: fry leading-edge driving proteins (top 20 by |t-stat|)"),
    stringsAsFactors = FALSE),
  sheet_specs = f04_specs
)
copy_to_box("04_Figures/F04/c_data/F04_supplementary.xlsx", "F04", "S10_F04_training_concord.xlsx")

# F05
cat("=== Rebuilding F05 workbook ===\n")
f05_specs <- list(
  list(name = "panel_A_ora_quadrant",       path = "04_Figures/F05/c_data/panel_A/ora_quadrant.csv"),
  list(name = "panel_B_pattern_class",      path = "04_Figures/F05/c_data/panel_B_heatmap/pattern_classification.csv"),
  list(name = "panel_B_sankey",             path = "04_Figures/F05/c_data/panel_B_heatmap/sankey_links.csv"),
  list(name = "panel_B_bar",                path = "04_Figures/F05/c_data/panel_B_heatmap/bar_data.csv"),
  list(name = "panel_C_fry_results",        path = "04_Figures/F05/c_data/panel_C_fry/fry_results_all.csv"),
  list(name = "panel_C_fry_driving",        path = "04_Figures/F05/c_data/panel_C_fry/driving_proteins.csv"),
  list(name = "panel_D_nes_scatter",        path = "04_Figures/F05/c_data/panel_D/nes_scatter.csv"),
  list(name = "panel_E_rrho2_summary",      path = "04_Figures/F05/c_data/panel_E/rrho2_summary.csv"),
  list(name = "SUPP_enrichment_reversal",   path = "04_Figures/F05/c_data/panel_supp/enrichment_reversal.csv"),
  list(name = "SUPP_reversal_pathway_stats", path = "04_Figures/F05/c_data/panel_supp/reversal_pathway_stats.csv"),
  list(name = "SUPP_ora_dedup_sensitivity", path = "04_Figures/F05/c_data/panel_supp/SUPP_ora_dedup_sensitivity.csv"),
  list(name = "SUPP_r_bootstrap",           path = "04_Figures/F05/c_data/panel_supp/SUPP_r_bootstrap.csv"),
  list(name = "SUPP_fry_circularity",       path = "04_Figures/F05/c_data/panel_supp/SUPP_fry_circularity.csv"),
  list(name = "SUPP_reversal_threshold",    path = "04_Figures/F05/c_data/panel_supp/SUPP_reversal_threshold.csv"),
  list(name = "SUPP_goslim_distribution",   path = "04_Figures/F05/c_data/panel_supp/SUPP_goslim_distribution.csv"),
  list(name = "SUPP_fry_leading_edge",      path = "04_Figures/F05/c_data/panel_supp/SUPP_fry_leading_edge.csv")
)
build_workbook(
  "04_Figures/F05/c_data/F05_supplementary.xlsx",
  title = "F05 \u2014 Figure 5 source data",
  description = "Aging-reversal diagnostics: quadrant ORA, pathway NES scatter, per-protein pattern classification, fry rotation test, RRHO2, and supplementary sensitivity analyses.",
  overview_df = data.frame(
    Sheet = vapply(f05_specs, `[[`, "", "name"),
    Description = c(
      "Panel A: ORA by reversal-quadrant scatter (Reversed Up/Down, Exacerbated Up/Down)",
      "Panel B: per-protein reversal pattern classification (heatmap source)",
      "Panel B: pathway-protein sankey links (heatmap source)",
      "Panel B: per-pattern bar chart counts",
      "Panel C: fry rotation test for reversal contrasts",
      "Panel C: reversal driving proteins (concordant sign-flip)",
      "Panel D: NES scatter (Aging vs Training_Old) per pathway with Spearman + Fisher Z",
      "Panel E: RRHO2 Aging vs Training_Old quadrant summary (max -log10p per quadrant)",
      "SUPP: pathway-level reversal enrichment",
      "SUPP: per-pathway reversal stats (NES Aging, NES TO, sign agreement)",
      "SUPP: ORA dedup Jaccard cutoff sensitivity (0.3, 0.5, 0.7, 1.0)",
      "SUPP: Pearson r bootstrap (1000 resamples with 95% CI)",
      "SUPP: fry circularity diagnostic (protein-permuted null distribution)",
      "SUPP: reversal classification threshold sensitivity (logFC 0.05-0.3)",
      "SUPP: GO Slim category distribution by reversal quadrant",
      "SUPP: fry leading-edge driving proteins (top 25 by |t-stat|)"),
    stringsAsFactors = FALSE),
  sheet_specs = f05_specs
)
copy_to_box("04_Figures/F05/c_data/F05_supplementary.xlsx", "F05", "S11_F05_aging_reversal.xlsx")

# F06
cat("=== Rebuilding F06 workbook ===\n")
# F06 has complex specs with RDS → df conversions; source the main stitch
# which already has the corrected Overview text, but skip the figure rendering.
# We run just the workbook section by loading the needed objects.
f06 <- function(p) file.path("04_Figures/F06/c_data", p)
.shared       <- readRDS(f06("shared_objects.rds"))
.MEs          <- matrix_to_df(as.matrix(readRDS(f06("MEs.rds"))),      "sample_id")
.me_pre       <- matrix_to_df(as.matrix(readRDS(f06("me_pre.rds"))),   "subject_key")
.me_post      <- matrix_to_df(as.matrix(readRDS(f06("me_post.rds"))),  "subject_key")
.delta_me     <- matrix_to_df(as.matrix(readRDS(f06("delta_me.rds"))), "subject_key")
.common_subj  <- data.frame(subject_key = .shared$common_subj, stringsAsFactors = FALSE)

f06_specs <- list(
  list(name="panel_A_heatmap",             path=f06("01_panel_A_heatmap_data.csv")),
  list(name="panel_B_module_fgsea",        path=f06("panel_B_module_fgsea.csv")),
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
  list(name="WGCNA_lmm_contrast_check",    path=f06("wgcna/wgcna_lmm_contrast_check.csv")),
  list(name="WGCNA_lmm_stratified_check",  path=f06("wgcna/wgcna_lmm_stratified_check.csv")),
  list(name="WGCNA_gs_phenotype_choices",  path=f06("wgcna/gs_phenotype_choices.csv")),
  list(name="metadata_samples",            path=f06("meta.csv")),
  list(name="metadata_subj_age",           path=f06("subj_age.csv")),
  list(name="metadata_pheno_wide",         path=f06("pheno_wide.csv")),
  list(name="metadata_imp_annotations",    path=f06("imp_annotations.csv")),
  list(name="MEs",         df=.MEs),
  list(name="me_pre",      df=.me_pre),
  list(name="me_post",     df=.me_post),
  list(name="delta_me",    df=.delta_me),
  list(name="common_subj", df=.common_subj)
)

overview_sheets <- vapply(f06_specs, `[[`, "", "name")
overview_descs <- c(
  "Panel A: 18-column heatmap data (modules x traits)",
  "Panels B & C: per-module fGSEA pathway NES for Training-Young, Training-Old, and Aging contrasts (source data for both NES scatter panels)",
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
  "Module eigengenes per sample",
  "Pre-timepoint module eigengenes per subject",
  "Post-timepoint module eigengenes per subject",
  "Change in module eigengenes (Post - Pre) per subject",
  "Common subject keys across Pre/Post for subject-paired analyses"
)

build_workbook(
  f06("F06_supplementary.xlsx"),
  title = "F06 \u2014 Figure 6 source data",
  description = "WGCNA module assignments, hub proteins, pathway enrichment, trait correlations, eigengene metadata.",
  overview_df = data.frame(Sheet = overview_sheets, Description = overview_descs,
                           stringsAsFactors = FALSE),
  sheet_specs = f06_specs
)
copy_to_box(f06("F06_supplementary.xlsx"), "F06", "S12_F06_WGCNA.xlsx")

cat("\n=== All workbooks rebuilt ===\n")
