# YvO WGCNA Runner — Generates all upstream data for Figures 5 & 6
# Inputs:  02_Imputation/c_data/01_imputed.csv, 01_DAList_imputed.rds
# Outputs: 04_Figures/F06/c_data/wgcna/ (network, assignments, hubs, correlations, GO)
# Reports: 04_Figures/F06/b_reports/ (soft threshold, dendrogram, heatmap as _SUPP)

library(WGCNA)
library(tidyverse)
library(lme4)
library(emmeans)
source("04_Figures/shared/style.R")
source("04_Figures/shared/pathway_utils.R")

allowWGCNAThreads()
set.seed(42)
setwd(rprojroot::find_rstudio_root_file())

DATA_FILE  <- "02_Imputation/c_data/01_imputed.csv"
DALIST_RDS <- "02_Imputation/c_data/01_DAList_imputed.rds"
REPORT_DIR      <- "04_Figures/F06/b_reports"
REPORT_SUPP_DIR <- "04_Figures/F06/b_reports/supp/01_QC"
DATA_DIR        <- "04_Figures/F06/c_data/wgcna"

dir.create(REPORT_DIR,      recursive = TRUE, showWarnings = FALSE)
dir.create(REPORT_SUPP_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(DATA_DIR,        recursive = TRUE, showWarnings = FALSE)

stopifnot(file.exists(DATA_FILE), file.exists(DALIST_RDS))

df <- read_csv(DATA_FILE, show_col_types = FALSE)
ann_cols   <- c("uniprot_id", "protein", "gene", "description")
ann        <- df[, ann_cols]
samp_names <- setdiff(names(df), ann_cols)
mat        <- as.matrix(df[, samp_names])
rownames(mat) <- ann$uniprot_id

# Transpose: samples as rows, proteins as columns (WGCNA convention)
datExpr <- t(mat)

message(sprintf("Data: %d samples x %d proteins", nrow(datExpr), ncol(datExpr)))
dal      <- readRDS(DALIST_RDS)
dal_meta <- as.data.frame(dal$metadata)

meta <- tibble(
  sample_id = dal_meta$Col_ID,
  subject   = sub("_(Pre|Post)$", "", dal_meta$Col_ID),
  age       = dal_meta$Group,
  time      = dal_meta$Timepoint,
  group     = dal_meta$Group_Time
)

# goodSamplesGenes check
gsg <- goodSamplesGenes(datExpr, verbose = 3)
if (!gsg$allOK) {
  datExpr <- datExpr[gsg$goodSamples, gsg$goodGenes]
  ann <- ann %>% filter(uniprot_id %in% colnames(datExpr))
  message(sprintf("After goodSamplesGenes: %d samples x %d proteins",
                  nrow(datExpr), ncol(datExpr)))
}

# WGCNA::cor / stats::cor dispatch conflict workaround (WGCNA 1.74 / R 4.5+)
cor <- WGCNA::cor

powers <- c(1:20)
sft <- pickSoftThreshold(datExpr, powerVector = powers,
                          networkType = "signed", verbose = 2)

# Find first power with R^2 > 0.85
# NOTE: Conventional threshold is R^2 > 0.90 (Zhang & Horvath 2005), relaxed to
# 0.85 for small-n proteomics per Langfelder & Horvath (2008) guidance.
# Selected power yields R^2 ~ 0.875 — acceptable for n=31 subjects, 2138 proteins.
r2_values <- -sign(sft$fitIndices$slope) * sft$fitIndices$SFT.R.sq
power_idx <- which(r2_values > 0.85)[1]
soft_power <- if (!is.na(power_idx)) powers[power_idx] else 6
# Log-log slope diagnostic: scale-free topology expects slope ~ -1 to -2
sft_slope <- sft$fitIndices$slope[soft_power]
r2_note <- paste0(
  ifelse(sft_slope > -1, " [NOTE: slope > -1, weak scale-free fit]", ""),
  ifelse(r2_values[soft_power] < 0.90,
         " [below 0.90 convention; acceptable for small-n, see Langfelder & Horvath 2008]", ""))
cat(sprintf("  Soft power: %d (R^2 = %.3f, slope = %.2f%s)\n",
            soft_power, r2_values[soft_power], sft_slope, r2_note))

plot_sft <- function() {
  par(mfrow = c(1, 2))
  plot(sft$fitIndices$Power, r2_values,
       xlab = "Soft Threshold (power)", ylab = "Scale Free Topology Model Fit (R^2)",
       main = "Scale independence", type = "n")
  text(sft$fitIndices$Power, r2_values, labels = powers, cex = 0.9, col = "red")
  abline(h = 0.85, col = "red", lty = 2)

  plot(sft$fitIndices$Power, sft$fitIndices$mean.k.,
       xlab = "Soft Threshold (power)", ylab = "Mean Connectivity",
       main = "Mean connectivity", type = "n")
  text(sft$fitIndices$Power, sft$fitIndices$mean.k., labels = powers, cex = 0.9, col = "red")
}
png(file.path(REPORT_SUPP_DIR, "SUPP_soft_threshold.png"), width = 3000, height = 1500, res = 300)
plot_sft()
dev.off()

# METHODS NOTE — Correlation choice: Pearson (not bicor)
# Rationale: DIA-MS with DIA-NN library-free search yields complete, approximately
# normally-distributed abundance matrices with fewer extreme outliers than label-based
# (TMT/iTRAQ) workflows, reducing the advantage of robust estimators like biweight
# midcorrelation. Pearson is the WGCNA default and is consistent with:
#   - O'Leary et al. 2024 (PMID 39663727): Pearson, muscle proteomics WGCNA
#   - Willis et al. 2020 (PMID 31910159): Pearson, muscle transcriptomics WGCNA
# Johnson et al. 2020 (Nat Med) / 2022 (Nat Neurosci) use bicor for TMT data across
# heterogeneous multi-cohort brain samples — different data characteristics.
# Empirical confirmation: supp/a02_bicor_sensitivity.R shows high module overlap (Jaccard)
# between Pearson and bicor networks for this dataset.
net <- blockwiseModules(
  datExpr,
  power             = soft_power,
  networkType       = "signed",
  TOMType           = "signed",
  minModuleSize     = 30,
  mergeCutHeight    = 0.25,
  numericLabels     = TRUE,
  pamRespectsDendro = FALSE,
  saveTOMs          = FALSE,
  verbose           = 3
)

module_colors <- labels2colors(net$colors)
n_modules <- length(unique(net$colors)) - (0 %in% net$colors)

cat(sprintf("  Modules detected: %d (+ grey/unassigned)\n", n_modules))
print(table(module_colors))

# Dendrogram plot removed — superseded by supp/a01_dendrogram.R (panel_D_dendrogram_SUPP)

# Enhanced trait matrix (3 design + 6 phenotype)
traits <- meta %>%
  mutate(
    age_num     = if_else(age == "Old", 1, 0),
    time_num    = if_else(time == "Post", 1, 0),
    interaction = age_num * time_num
  )

# Continuous phenotypes from DAList metadata
pheno_cols <- c("VL_thick_cm", "DXA_LBM_kg", "BMI",
                "Type_I_fCSA", "Type_II_fCSA", "deadlift_1rm_kg")

for (pc in pheno_cols) {
  if (pc %in% names(dal_meta)) {
    vals <- dal_meta[[pc]]
    if (!is.numeric(vals)) vals <- as.numeric(as.character(vals))
    traits[[pc]] <- vals[match(meta$sample_id, dal_meta$Col_ID)]
  }
}

trait_cols <- c("age_num", "time_num", "interaction", pheno_cols)
trait_cols <- intersect(trait_cols, names(traits))
traits_mat <- as.data.frame(traits[, trait_cols])
rownames(traits_mat) <- meta$sample_id

traits_mat <- traits_mat[rownames(datExpr), ]

MEs <- moduleEigengenes(datExpr, colors = module_colors)$eigengenes
MEs <- orderMEs(MEs)

# Correlations and raw p-values
module_trait_cor  <- cor(MEs, traits_mat, use = "pairwise.complete.obs")
module_trait_pval <- corPvalueStudent(module_trait_cor, nrow(datExpr))

# STAT AUDIT: corPvalueStudent() uses n=62 for ALL traits — slightly liberal
# for phenotype traits with missing data. Panel B recomputes cor.test() CIs.
# BH correction: GLOBAL across all modules x all traits.
pval_vec <- as.vector(module_trait_pval)
pval_bh_vec <- p.adjust(pval_vec, method = "BH")
module_trait_pval_bh <- matrix(pval_bh_vec, nrow = nrow(module_trait_pval),
                                ncol = ncol(module_trait_pval))
rownames(module_trait_pval_bh) <- rownames(module_trait_pval)
colnames(module_trait_pval_bh) <- colnames(module_trait_pval)

cat(sprintf("  BH correction: %d tests (%d modules x %d traits)\n",
            length(pval_vec), nrow(module_trait_pval), ncol(module_trait_pval)))
star_matrix <- ifelse(module_trait_pval_bh < 0.001, "***",
               ifelse(module_trait_pval_bh < 0.01, "**",
               ifelse(module_trait_pval_bh < 0.05, "*", "")))
text_matrix <- paste(signif(module_trait_cor, 2), star_matrix, sep = "\n")
dim(text_matrix) <- dim(module_trait_cor)

plot_trait_heatmap <- function() {
  par(mar = c(6, 10, 3, 3))
  labeledHeatmap(
    Matrix     = module_trait_cor,
    xLabels    = colnames(traits_mat),
    yLabels    = colnames(MEs),
    ySymbols   = colnames(MEs),
    colorLabels = FALSE,
    colors     = blueWhiteRed(50),
    textMatrix = text_matrix,
    setStdMargins = FALSE,
    cex.text   = 0.5,
    zlim       = c(-1, 1),
    main       = "Module-trait correlations (* BH < 0.05)"
  )
}
png(file.path(REPORT_SUPP_DIR, "SUPP_module_trait_heatmap.png"), width = 3000, height = 3000, res = 300)
plot_trait_heatmap()
dev.off()

trait_cor_df <- as.data.frame(module_trait_cor) %>%
  rownames_to_column("module")
write_csv(trait_cor_df, file.path(DATA_DIR, "wgcna_module_trait_correlations.csv"))

pval_bh_df <- as.data.frame(module_trait_pval_bh) %>%
  rownames_to_column("module")
write_csv(pval_bh_df, file.path(DATA_DIR, "wgcna_module_trait_pvalues_bh.csv"))

kME <- signedKME(datExpr, MEs)

module_df <- tibble(
  uniprot_id   = colnames(datExpr),
  module_color = module_colors,
  module_num   = net$colors
) %>% left_join(ann %>% dplyr::select(uniprot_id, gene), by = "uniprot_id")

hub_rows <- list()
unique_modules <- setdiff(unique(module_colors), "grey")

for (mod in unique_modules) {
  mod_proteins <- module_df$uniprot_id[module_df$module_color == mod]
  kme_col <- paste0("kME", mod)
  if (!(kme_col %in% colnames(kME))) next

  mod_kme <- tibble(
    uniprot_id = rownames(kME),
    kME = kME[, kme_col]
  ) %>%
    filter(uniprot_id %in% mod_proteins) %>%
    arrange(desc(abs(kME))) %>%
    head(10) %>%
    mutate(module = mod) %>%
    left_join(ann %>% dplyr::select(uniprot_id, gene), by = "uniprot_id")

  hub_rows <- c(hub_rows, list(mod_kme))
}

hub_df <- bind_rows(hub_rows)
cor <- stats::cor  # restore after WGCNA computations

# Multi-database ORA enrichment (H + KEGG + Reactome + WP + GO:BP + GO Slim)
# Replaces enrichGO (GO:BP only) with unified run_ora_deduplicated()
bg_genes <- ann$gene[ann$uniprot_id %in% colnames(datExpr)]
bg_genes <- unique(bg_genes[!is.na(bg_genes) & bg_genes != ""])

pw_collection <- build_pathway_collection(min_size = 15, max_size = 500,
                                          include_goslim = FALSE)

ora_results_list <- list()

for (mod in unique_modules) {
  mod_genes <- module_df$gene[module_df$module_color == mod]
  mod_genes <- unique(mod_genes[!is.na(mod_genes) & mod_genes != ""])

  if (length(mod_genes) < 5) next

  ora_res <- tryCatch(
    run_ora_deduplicated(
      genes          = mod_genes,
      universe       = bg_genes,
      pathways       = pw_collection,
      jaccard_cutoff = 0.5,
      min_size       = 15,
      max_size       = 500,
      padj_cutoff    = 0.05
    ),
    error = function(e) { warning(sprintf("ORA failed for '%s': %s", mod, e$message)); NULL }
  )

  if (!is.null(ora_res) && nrow(ora_res) > 0) {
    ora_res$module <- mod
    # Backward-compat columns for downstream panels
    ora_res$Description <- clean_pathway_name(ora_res$pathway)
    ora_res$geneID <- vapply(ora_res$overlapGenes, function(g) {
      paste(g, collapse = "/")
    }, character(1))
    ora_res$Count     <- ora_res$overlap
    ora_res$p.adjust  <- ora_res$padj
    ora_res$ID        <- ora_res$pathway
    ora_results_list  <- c(ora_results_list, list(ora_res))
  }
}

enrich_df <- bind_rows(ora_results_list)
# Drop list column (overlapGenes) for CSV serialization
enrich_df$overlapGenes <- NULL

write_csv(module_df,  file.path(DATA_DIR, "wgcna_module_assignments.csv"))
write_csv(hub_df,     file.path(DATA_DIR, "wgcna_hub_proteins.csv"))
saveRDS(net,          file.path(DATA_DIR, "wgcna_network.rds"))
write_csv(enrich_df,  file.path(DATA_DIR, "wgcna_module_enrichment.csv"))
# Legacy alias for backward compatibility
write_csv(enrich_df,  file.path(DATA_DIR, "wgcna_module_GO_enrichment.csv"))

# Soft-threshold summary for Methods (Langfelder & Horvath 2017)
sft_summary <- tibble(
  selected_power    = soft_power,
  R_squared         = r2_values[soft_power],
  mean_connectivity = sft$fitIndices$mean.k.[soft_power],
  n_proteins        = ncol(datExpr),
  n_samples         = nrow(datExpr)
)
write_csv(sft_summary, file.path(DATA_DIR, "wgcna_sft_summary.csv"))

# ── Panel-ready intermediates (consumed directly by self-contained panels) ──
# These replace the old prepare_data.R — every panel reads from c_data/.

PANEL_DIR <- "04_Figures/F06/c_data"
dir.create(PANEL_DIR, recursive = TRUE, showWarnings = FALSE)

# Group factor for consistent ordering
meta$group <- factor(meta$group,
                     levels = c("Young_Pre", "Young_Post", "Old_Pre", "Old_Post"))

# Phenotype columns from DAList metadata
pheno_cols_panel <- c("VL_thick_cm", "DXA_LBM_kg", "BMI",
                      "Type_I_fCSA", "Type_II_fCSA", "deadlift_1rm_kg")
for (pc in pheno_cols_panel) {
  if (pc %in% names(dal_meta) && !(pc %in% names(meta))) {
    vals <- dal_meta[[pc]]
    if (!is.numeric(vals)) vals <- as.numeric(as.character(vals))
    meta[[pc]] <- vals[match(meta$sample_id, dal_meta$Col_ID)]
  }
}

# Biological module labels — M1-M9 ordered by module size (largest first)
mod_sizes <- sort(table(module_colors[module_colors != "grey"]), decreasing = TRUE)
mod_bio_lookup <- c(
  blue      = "Cell Cycle/Proteostasis",
  brown     = "Redox/Glycolysis",
  turquoise = "Lipid Catabolism",
  green     = "Oxidative Phosphorylation",
  black     = "Translation Machinery",
  pink      = "Mitochondrial Biogenesis",
  yellow    = "Muscle Contraction",
  red       = "Immunity/Heme",
  magenta   = "Translation Initiation"
)
mod_bio_labels <- tibble(
  module_color = names(mod_sizes),
  module_id    = paste0("M", seq_along(mod_sizes)),
  bio_label    = mod_bio_lookup[names(mod_sizes)],
  n_proteins   = as.integer(mod_sizes),
  display_label = paste0(module_id, ": ", bio_label)
)

# Core objects (panels A-E). imp_mat.rds was formerly written here but had no
# active consumers (superseded by datExpr.rds); removed 2026-04-16 to keep
# c_data/ lean. See .archive_orphan/ for any legacy copies.
saveRDS(MEs,            file.path(PANEL_DIR, "MEs.rds"))
saveRDS(kME,            file.path(PANEL_DIR, "kME_all.rds"))
saveRDS(datExpr,        file.path(PANEL_DIR, "datExpr.rds"))
saveRDS(module_colors,  file.path(PANEL_DIR, "module_colors.rds"))
write_csv(meta,         file.path(PANEL_DIR, "meta.csv"))
write_csv(mod_bio_labels, file.path(PANEL_DIR, "mod_bio_labels.csv"))
write_csv(ann,          file.path(PANEL_DIR, "imp_annotations.csv"))

# Pre/Post eigengene split + phenotype pairing (panels F, G, supp)
pre_meta  <- meta %>% filter(time == "Pre")
post_meta <- meta %>% filter(time == "Post")

me_pre_raw  <- MEs[pre_meta$sample_id, , drop = FALSE]
me_post_raw <- MEs[post_meta$sample_id, , drop = FALSE]
pre_subjects  <- sub("_(Pre|Post)$", "", pre_meta$sample_id)
post_subjects <- sub("_(Pre|Post)$", "", post_meta$sample_id)
rownames(me_pre_raw)  <- pre_subjects
rownames(me_post_raw) <- post_subjects

common_subj <- intersect(rownames(me_pre_raw), rownames(me_post_raw))
me_pre  <- me_pre_raw[common_subj, , drop = FALSE]
me_post <- me_post_raw[common_subj, , drop = FALSE]
delta_me <- me_post - me_pre

# Expanded phenotype pairing (all available phenotypes)
pheno_pre <- meta %>%
  filter(time == "Pre") %>%
  mutate(subject_key = sub("_(Pre|Post)$", "", sample_id)) %>%
  dplyr::select(subject_key, VL_thick_cm, DXA_LBM_kg, BMI,
                deadlift_1rm_kg, Type_I_fCSA, Type_II_fCSA) %>%
  rename(VL_Pre = VL_thick_cm, LBM_Pre = DXA_LBM_kg, BMI_Pre = BMI,
         DL_Pre = deadlift_1rm_kg, T1_Pre = Type_I_fCSA, T2_Pre = Type_II_fCSA)

pheno_post <- meta %>%
  filter(time == "Post") %>%
  mutate(subject_key = sub("_(Pre|Post)$", "", sample_id)) %>%
  dplyr::select(subject_key, VL_thick_cm, DXA_LBM_kg,
                deadlift_1rm_kg, Type_I_fCSA, Type_II_fCSA) %>%
  rename(VL_Post = VL_thick_cm, LBM_Post = DXA_LBM_kg,
         DL_Post = deadlift_1rm_kg, T1_Post = Type_I_fCSA, T2_Post = Type_II_fCSA)

pheno_wide <- inner_join(pheno_pre, pheno_post, by = "subject_key") %>%
  mutate(delta_VL  = VL_Post  - VL_Pre,
         delta_LBM = LBM_Post - LBM_Pre,
         delta_DL  = DL_Post  - DL_Pre,
         delta_T1  = T1_Post  - T1_Pre,
         delta_T2  = T2_Post  - T2_Pre) %>%
  filter(subject_key %in% common_subj)

subj_age <- meta %>%
  filter(time == "Pre") %>%
  mutate(subject_key = sub("_(Pre|Post)$", "", sample_id)) %>%
  dplyr::select(subject_key, age) %>%
  distinct()

# Baseline expression (Pre only, keyed by subject)
pre_expr <- datExpr[pre_meta$sample_id, , drop = FALSE]
rownames(pre_expr) <- pre_subjects
pre_expr <- pre_expr[common_subj, , drop = FALSE]

# Gene significance: baseline expression vs delta phenotype (all available)
delta_vl_vec  <- pheno_wide$delta_VL[match(common_subj, pheno_wide$subject_key)]
delta_lbm_vec <- pheno_wide$delta_LBM[match(common_subj, pheno_wide$subject_key)]
delta_dl_vec  <- pheno_wide$delta_DL[match(common_subj, pheno_wide$subject_key)]
delta_t1_vec  <- pheno_wide$delta_T1[match(common_subj, pheno_wide$subject_key)]
delta_t2_vec  <- pheno_wide$delta_T2[match(common_subj, pheno_wide$subject_key)]

gs_vl <- cor(pre_expr, delta_vl_vec, use = "pairwise.complete.obs")
colnames(gs_vl) <- "GS_deltaVL"
gs_lbm <- cor(pre_expr, delta_lbm_vec, use = "pairwise.complete.obs")
colnames(gs_lbm) <- "GS_deltaLBM"
gs_dl <- cor(pre_expr, delta_dl_vec, use = "pairwise.complete.obs")
colnames(gs_dl) <- "GS_deltaDL"
gs_t1 <- cor(pre_expr, delta_t1_vec, use = "pairwise.complete.obs")
colnames(gs_t1) <- "GS_deltaT1"
gs_t2 <- cor(pre_expr, delta_t2_vec, use = "pairwise.complete.obs")
colnames(gs_t2) <- "GS_deltaT2"

# Per-module prediction correlations (baseline ME vs all delta phenotypes)
pred_cor <- tibble(module = colnames(me_pre)) %>%
  rowwise() %>%
  mutate(
    r_vl  = cor(me_pre[common_subj, module],
                pheno_wide$delta_VL[match(common_subj, pheno_wide$subject_key)],
                use = "complete.obs"),
    r_lbm = cor(me_pre[common_subj, module],
                pheno_wide$delta_LBM[match(common_subj, pheno_wide$subject_key)],
                use = "complete.obs"),
    r_dl  = cor(me_pre[common_subj, module],
                pheno_wide$delta_DL[match(common_subj, pheno_wide$subject_key)],
                use = "complete.obs"),
    r_t1  = cor(me_pre[common_subj, module],
                pheno_wide$delta_T1[match(common_subj, pheno_wide$subject_key)],
                use = "complete.obs"),
    r_t2  = cor(me_pre[common_subj, module],
                pheno_wide$delta_T2[match(common_subj, pheno_wide$subject_key)],
                use = "complete.obs"),
    max_r = max(abs(r_vl), abs(r_lbm), na.rm = TRUE)
  ) %>% ungroup() %>%
  mutate(across(c(r_vl, r_lbm, r_dl, r_t1, r_t2, max_r), as.numeric))

top3 <- pred_cor %>% arrange(desc(max_r)) %>% head(3) %>% pull(module)

# Focus modules for triptych (Panel C) and hub networks (Panel D)
# Key modules for triptych/hub supplementary panels.
# Canonical source: c_data/wgcna/key_modules.txt (written during network build).
# Read from wgcna/ rather than hardcoding to keep a single source of truth.
KEY_MODULES <- readLines(file.path(DATA_DIR, "key_modules.txt"))
KEY_MODULES <- KEY_MODULES[nzchar(trimws(KEY_MODULES))]
writeLines(KEY_MODULES, file.path(PANEL_DIR, "key_modules.txt"))

# ── Baseline & change trait correlation matrices (for restructured Panel A) ──
# Separates cross-sectional (design) from baseline and longitudinal signals.

all_mods <- colnames(me_pre_raw)  # all modules including grey

# Baseline correlations: ME_Pre (all Pre subjects) vs baseline phenotype
bl_meta <- meta %>% filter(time == "Pre") %>%
  mutate(subject_key = sub("_(Pre|Post)$", "", sample_id))

baseline_traits <- data.frame(
  BMI_Pre = bl_meta$BMI[match(pre_subjects, bl_meta$subject_key)],
  VL_Pre  = bl_meta$VL_thick_cm[match(pre_subjects, bl_meta$subject_key)],
  LBM_Pre = bl_meta$DXA_LBM_kg[match(pre_subjects, bl_meta$subject_key)],
  row.names = pre_subjects
)

bl_cor_mat <- cor(me_pre_raw[pre_subjects, all_mods, drop = FALSE],
                  baseline_traits, use = "pairwise.complete.obs")
bl_pval_mat <- matrix(NA_real_, nrow = length(all_mods), ncol = ncol(baseline_traits),
                       dimnames = list(all_mods, colnames(baseline_traits)))
for (trait in colnames(baseline_traits)) {
  n_ok <- sum(!is.na(baseline_traits[pre_subjects, trait]))
  bl_pval_mat[, trait] <- corPvalueStudent(bl_cor_mat[, trait], n_ok)
}

# Change correlations: delta_ME (paired subjects) vs delta phenotype
change_traits <- data.frame(
  delta_VL  = pheno_wide$delta_VL[match(common_subj, pheno_wide$subject_key)],
  delta_LBM = pheno_wide$delta_LBM[match(common_subj, pheno_wide$subject_key)],
  row.names = common_subj
)

ch_cor_mat <- cor(delta_me[common_subj, all_mods, drop = FALSE],
                  change_traits, use = "pairwise.complete.obs")
ch_pval_mat <- matrix(NA_real_, nrow = length(all_mods), ncol = ncol(change_traits),
                       dimnames = list(all_mods, colnames(change_traits)))
for (trait in colnames(change_traits)) {
  n_ok <- sum(!is.na(change_traits[common_subj, trait]))
  ch_pval_mat[, trait] <- corPvalueStudent(ch_cor_mat[, trait], n_ok)
}

# ── Stratified baseline & change correlations (Young vs Old separately) ──
# Unmasks age-specific signals hidden when pooling Young+Old together.
young_pre_subj <- subj_age$subject_key[subj_age$age == "Young"]
old_pre_subj   <- subj_age$subject_key[subj_age$age == "Old"]

# Young baseline (n≈16 Pre)
bl_subj_y   <- intersect(pre_subjects, young_pre_subj)
bl_traits_y <- baseline_traits[bl_subj_y, , drop = FALSE]
me_pre_y    <- me_pre_raw[bl_subj_y, all_mods, drop = FALSE]
bl_cor_young <- cor(me_pre_y, bl_traits_y, use = "pairwise.complete.obs")
bl_pval_young <- matrix(NA_real_, nrow = length(all_mods), ncol = ncol(bl_traits_y),
                         dimnames = list(all_mods, colnames(bl_traits_y)))
for (trait in colnames(bl_traits_y)) {
  n_ok <- sum(!is.na(bl_traits_y[bl_subj_y, trait]))
  bl_pval_young[, trait] <- corPvalueStudent(bl_cor_young[, trait], n_ok)
}

# Old baseline (n≈15 Pre)
bl_subj_o   <- intersect(pre_subjects, old_pre_subj)
bl_traits_o <- baseline_traits[bl_subj_o, , drop = FALSE]
me_pre_o    <- me_pre_raw[bl_subj_o, all_mods, drop = FALSE]
bl_cor_old  <- cor(me_pre_o, bl_traits_o, use = "pairwise.complete.obs")
bl_pval_old <- matrix(NA_real_, nrow = length(all_mods), ncol = ncol(bl_traits_o),
                       dimnames = list(all_mods, colnames(bl_traits_o)))
for (trait in colnames(bl_traits_o)) {
  n_ok <- sum(!is.na(bl_traits_o[bl_subj_o, trait]))
  bl_pval_old[, trait] <- corPvalueStudent(bl_cor_old[, trait], n_ok)
}

# Young change (n≈15 paired)
common_young <- intersect(common_subj, young_pre_subj)
ch_traits_y  <- change_traits[common_young, , drop = FALSE]
dme_y        <- delta_me[common_young, all_mods, drop = FALSE]
ch_cor_young <- cor(dme_y, ch_traits_y, use = "pairwise.complete.obs")
ch_pval_young <- matrix(NA_real_, nrow = length(all_mods), ncol = ncol(ch_traits_y),
                         dimnames = list(all_mods, colnames(ch_traits_y)))
for (trait in colnames(ch_traits_y)) {
  n_ok <- sum(!is.na(ch_traits_y[common_young, trait]))
  ch_pval_young[, trait] <- corPvalueStudent(ch_cor_young[, trait], n_ok)
}

# Old change (n≈15 paired)
common_old  <- intersect(common_subj, old_pre_subj)
ch_traits_o <- change_traits[common_old, , drop = FALSE]
dme_o       <- delta_me[common_old, all_mods, drop = FALSE]
ch_cor_old  <- cor(dme_o, ch_traits_o, use = "pairwise.complete.obs")
ch_pval_old <- matrix(NA_real_, nrow = length(all_mods), ncol = ncol(ch_traits_o),
                       dimnames = list(all_mods, colnames(ch_traits_o)))
for (trait in colnames(ch_traits_o)) {
  n_ok <- sum(!is.na(ch_traits_o[common_old, trait]))
  ch_pval_old[, trait] <- corPvalueStudent(ch_cor_old[, trait], n_ok)
}

cat(sprintf("  Stratified: Young baseline n=%d, Old baseline n=%d, Young change n=%d, Old change n=%d\n",
            length(bl_subj_y), length(bl_subj_o), length(common_young), length(common_old)))

# ── LMM contrasts (replaces Pearson design section for Panel A) ──
# corPvalueStudent() treats 62 paired samples as independent — anti-conservative.
# LMM properly models repeated measures: eigengene ~ group + (1|subject)
# Contrasts match DEP definitions: Aging, Training_Young, Training_Old, Interaction.
# Reference: Li et al. 2018 (WGCNA + mixed models for repeated-measures designs).

lmm_data <- meta %>%
  mutate(group = factor(group, levels = c("Young_Pre", "Young_Post", "Old_Pre", "Old_Post")))

lmm_contrast_list <- list(
  Aging          = c(-1,  0,  1,  0),
  Training_Young = c(-1,  1,  0,  0),
  Training_Old   = c( 0,  0, -1,  1),
  Interaction    = c( 1, -1, -1,  1)
)

lmm_rows <- list()
for (mod in all_mods) {
  lmm_data[[mod]] <- MEs[lmm_data$sample_id, mod]

  fit <- tryCatch(
    suppressWarnings(
      lmer(as.formula(paste0("`", mod, "` ~ group + (1 | subject)")),
           data = lmm_data, REML = TRUE)
    ),
    error = function(e) { warning(sprintf("LMM failed for %s: %s", mod, e$message)); NULL }
  )
  if (is.null(fit)) next

  singular <- isSingular(fit)
  if (singular) message(sprintf("  Note: singular fit for %s (near-zero random variance)", mod))

  emm <- emmeans(fit, ~ group)

  for (cname in names(lmm_contrast_list)) {
    ctr <- contrast(emm, list(ctr = lmm_contrast_list[[cname]]))
    s <- summary(ctr, ddf = "Kenward-Roger")
    est  <- s$estimate
    se   <- s$SE
    df_k <- s$df
    tval <- s$t.ratio
    praw <- s$p.value
    r_equiv <- sign(est) * sqrt(tval^2 / (tval^2 + df_k))

    lmm_rows <- c(lmm_rows, list(tibble(
      module   = mod,
      contrast = cname,
      estimate = round(est, 5),
      SE       = round(se, 5),
      df       = round(df_k, 2),
      t_ratio  = round(tval, 4),
      p_raw    = praw,
      r_equiv  = round(r_equiv, 4),
      singular = singular
    )))
  }
}

lmm_df <- bind_rows(lmm_rows)
cat(sprintf("  LMM contrasts: %d tests (%d modules x %d contrasts)\n",
            nrow(lmm_df), length(all_mods), length(lmm_contrast_list)))

# BH correction strategy:
#   LMM: per-section BH across 40 tests (confirmatory, pre-specified contrasts)
#   Stratified: per-column BH across modules within each trait (exploratory)
# Rationale: pooling confirmatory LMM (n=62) with exploratory stratified (n~15)
# into one global family is overly conservative — the LMM p-values inflate the
# BH threshold, penalizing the underpowered stratified sections.
# Per-column BH matches WGCNA vignette convention ("which module for this trait?").
lmm_df$p_bh <- p.adjust(lmm_df$p_raw, method = "BH")

per_col_bh <- function(pmat) {
  bh_mat <- pmat
  for (j in seq_len(ncol(pmat))) {
    bh_mat[, j] <- p.adjust(pmat[, j], method = "BH")
  }
  bh_mat
}
bl_pval_bh_young <- per_col_bh(bl_pval_young)
bl_pval_bh_old   <- per_col_bh(bl_pval_old)
ch_pval_bh_young <- per_col_bh(ch_pval_young)
ch_pval_bh_old   <- per_col_bh(ch_pval_old)

cat(sprintf("  BH correction: LMM per-section (%d tests) | stratified per-column (%d modules)\n",
            nrow(lmm_df), nrow(bl_pval_young)))

# Backward-compat BH for old global CSVs (same 90-test pool as before)
n_l <- nrow(lmm_df)
compat_raw <- c(lmm_df$p_raw, as.vector(bl_pval_mat), as.vector(ch_pval_mat))
compat_bh  <- p.adjust(compat_raw, method = "BH")
n_b <- length(as.vector(bl_pval_mat))
n_c <- length(as.vector(ch_pval_mat))
bl_pval_bh <- matrix(compat_bh[(n_l + 1):(n_l + n_b)], nrow = nrow(bl_pval_mat),
                      dimnames = dimnames(bl_pval_mat))
ch_pval_bh <- matrix(compat_bh[(n_l + n_b + 1):(n_l + n_b + n_c)], nrow = nrow(ch_pval_mat),
                      dimnames = dimnames(ch_pval_mat))

# Save LMM contrast audit (full emmeans output for reproducibility)
write_csv(lmm_df, file.path(DATA_DIR, "wgcna_lmm_contrast_audit.csv"))

# ── Stratified LMM: within-age training contrasts ──
# Fits ME ~ time + (1|subject) within each age group separately.
# Produces wgcna_lmm_stratified_audit.csv consumed by panel_A.R and a08.
strat_rows <- list()
for (age_grp in c("Young", "Old")) {
  strat_data <- lmm_data %>%
    filter(grepl(age_grp, group)) %>%
    mutate(time = factor(ifelse(grepl("Post", group), "Post", "Pre"),
                         levels = c("Pre", "Post")))

  for (mod in all_mods) {
    strat_data[["me_val"]] <- MEs[strat_data$sample_id, mod]

    fit_s <- tryCatch(
      suppressWarnings(
        lmer(me_val ~ time + (1 | subject), data = strat_data, REML = TRUE)
      ),
      error = function(e) NULL
    )
    if (is.null(fit_s)) next

    sing_s <- isSingular(fit_s)
    emm_s  <- emmeans(fit_s, ~ time)
    ctr_s  <- contrast(emm_s, list(training = c(-1, 1)))
    s_s    <- summary(ctr_s, ddf = "Kenward-Roger")
    t_s    <- s_s$t.ratio
    df_s   <- s_s$df
    r_eq_s <- sign(s_s$estimate) * sqrt(t_s^2 / (t_s^2 + df_s))

    strat_rows <- c(strat_rows, list(tibble(
      age_group = age_grp,
      module    = mod,
      estimate  = round(s_s$estimate, 5),
      SE        = round(s_s$SE, 5),
      df        = round(df_s, 2),
      t_ratio   = round(t_s, 4),
      p_raw     = s_s$p.value,
      r_equiv   = round(r_eq_s, 4),
      singular  = sing_s
    )))
  }
}
strat_df <- bind_rows(strat_rows)
strat_df$p_bh <- p.adjust(strat_df$p_raw, method = "BH")
write_csv(strat_df, file.path(DATA_DIR, "wgcna_lmm_stratified_audit.csv"))
cat(sprintf("  Stratified LMM: %d tests (%d modules x %d age groups)\n",
            nrow(strat_df), length(all_mods), 2))

# Save baseline correlation + global BH p-values
bl_cor_df  <- as.data.frame(bl_cor_mat) %>% rownames_to_column("module")
bl_pval_df <- as.data.frame(bl_pval_bh) %>% rownames_to_column("module")
write_csv(bl_cor_df,  file.path(DATA_DIR, "wgcna_baseline_trait_correlations.csv"))
write_csv(bl_pval_df, file.path(DATA_DIR, "wgcna_baseline_trait_pvalues_bh.csv"))

# Save change correlation + global BH p-values
ch_cor_df  <- as.data.frame(ch_cor_mat) %>% rownames_to_column("module")
ch_pval_df <- as.data.frame(ch_pval_bh) %>% rownames_to_column("module")
write_csv(ch_cor_df,  file.path(DATA_DIR, "wgcna_change_trait_correlations.csv"))
write_csv(ch_pval_df, file.path(DATA_DIR, "wgcna_change_trait_pvalues_bh.csv"))

# Stratified CSVs (per-column BH + raw p-values for nominal tier display)
strat_outputs <- list(
  list(prefix = "baseline_trait_correlations", young = bl_cor_young,     old = bl_cor_old),
  list(prefix = "baseline_trait_pvalues_bh",   young = bl_pval_bh_young, old = bl_pval_bh_old),
  list(prefix = "baseline_trait_pvalues_raw",  young = bl_pval_young,   old = bl_pval_old),
  list(prefix = "change_trait_correlations",   young = ch_cor_young,    old = ch_cor_old),
  list(prefix = "change_trait_pvalues_bh",     young = ch_pval_bh_young, old = ch_pval_bh_old),
  list(prefix = "change_trait_pvalues_raw",    young = ch_pval_young,   old = ch_pval_old)
)
for (out in strat_outputs) {
  for (age in c("young", "old")) {
    write_csv(as.data.frame(out[[age]]) %>% rownames_to_column("module"),
              file.path(DATA_DIR, paste0("wgcna_", out$prefix, "_", age, ".csv")))
  }
}

# Save panel-ready intermediates
# Panel-ready intermediates. Removed 2026-04-16 (no active consumers):
#   pre_expr.rds, gs_{vl,lbm,dl,t1,t2}.rds, pred_cor.csv
# Legacy copies swept to c_data/.archive_orphan/ during the same sweep.
saveRDS(me_pre,        file.path(PANEL_DIR, "me_pre.rds"))
saveRDS(me_post,       file.path(PANEL_DIR, "me_post.rds"))
saveRDS(delta_me,      file.path(PANEL_DIR, "delta_me.rds"))
write_csv(pheno_wide,  file.path(PANEL_DIR, "pheno_wide.csv"))
write_csv(subj_age,    file.path(PANEL_DIR, "subj_age.csv"))
# module_df.csv (root) removed 2026-04-16 — duplicated wgcna/wgcna_module_assignments.csv.
# Consumers read the wgcna/ version or the WGCNA_module_assignments xlsx sheet.

saveRDS(list(
  top3           = top3,
  common_subj    = common_subj,
  pre_subjects   = pre_subjects,
  mod_bio_labels = setNames(mod_bio_labels$bio_label, mod_bio_labels$module_color),
  outcome_labels = c(delta_VL = "Delta VL (cm)", delta_LBM = "Delta LBM (kg)"),
  delta_vl_vec   = delta_vl_vec,
  delta_lbm_vec  = delta_lbm_vec,
  delta_dl_vec   = delta_dl_vec,
  delta_t1_vec   = delta_t1_vec,
  delta_t2_vec   = delta_t2_vec
), file.path(PANEL_DIR, "shared_objects.rds"))

cat(sprintf("Done: %d modules, %d hub proteins, %d enriched pathways, %d paired subjects\n",
            n_modules, nrow(hub_df), nrow(enrich_df), length(common_subj)))
