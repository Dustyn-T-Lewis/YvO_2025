################################################################################
#
#   YvO WGCNA Runner — Generates all upstream data for Figures 5 & 6
#
#   Workflow:
#     1.  Setup & data loading (imputed matrix + DAList metadata)
#     2.  Soft-thresholding power selection (R^2 > 0.85)
#     3.  Network construction (blockwiseModules, signed)
#     4.  Module dendrogram plot
#     5.  Enhanced trait matrix (9 traits: 3 design + 6 phenotype)
#     6.  Module-trait correlations + BH-corrected p-values
#     7.  Hub proteins (top 10 per module by kME)
#     8.  GO:BP enrichment per module (clusterProfiler)
#     9.  Export all 6 files
#
#   Inputs:
#     02_Imputation/c_data/01_imputed.csv       (2124 proteins x 62 samples)
#     02_Imputation/c_data/01_DAList_imputed.rds (metadata with phenotypes)
#
#   Outputs → 04_Figures/F5/c_data/wgcna/
#     wgcna_network.rds                    blockwiseModules object
#     wgcna_module_assignments.csv         uniprot_id, module_color, module_num, gene
#     wgcna_hub_proteins.csv               uniprot_id, kME, module, gene (top 10/module)
#     wgcna_module_trait_correlations.csv   module x trait correlation matrix
#     wgcna_module_trait_pvalues_bh.csv     BH-corrected p-values
#     wgcna_module_GO_enrichment.csv        GO:BP enrichment per module
#
#   Reports → 04_Figures/F5/b_reports/wgcna/
#     01_soft_threshold.pdf
#     02_module_dendrogram.pdf
#     03_module_trait_heatmap.pdf
#
################################################################################

# ==============================================================================
# 1.  SETUP & DATA LOADING
# ==============================================================================

cat("=== YvO WGCNA Runner ===\n\n")
cat(">> 1 -- Setup & Data Loading\n")

suppressPackageStartupMessages({
  library(WGCNA)
  library(tidyverse)
  library(clusterProfiler)
  library(org.Hs.eg.db)
})

allowWGCNAThreads()
set.seed(42)
setwd(rprojroot::find_rstudio_root_file())

DATA_FILE  <- "02_Imputation/c_data/01_imputed.csv"
DALIST_RDS <- "02_Imputation/c_data/01_DAList_imputed.rds"
REPORT_DIR <- "04_Figures/F5/b_reports/wgcna"
DATA_DIR   <- "04_Figures/F5/c_data/wgcna"

dir.create(REPORT_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(DATA_DIR,   recursive = TRUE, showWarnings = FALSE)

stopifnot(file.exists(DATA_FILE), file.exists(DALIST_RDS))

df <- read_csv(DATA_FILE, show_col_types = FALSE)
ann_cols   <- c("uniprot_id", "protein", "gene", "description")
ann        <- df[, ann_cols]
samp_names <- setdiff(names(df), ann_cols)
mat        <- as.matrix(df[, samp_names])
rownames(mat) <- ann$uniprot_id

# Transpose: samples as rows, proteins as columns (WGCNA convention)
datExpr <- t(mat)

cat(sprintf("   Data: %d samples x %d proteins\n", nrow(datExpr), ncol(datExpr)))

# Build metadata from DAList (authoritative source for phenotypes)
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
  cat(sprintf("   After goodSamplesGenes: %d samples x %d proteins\n",
              nrow(datExpr), ncol(datExpr)))
} else {
  cat("   goodSamplesGenes: all OK\n")
}

# ==============================================================================
# 2.  SOFT-THRESHOLDING POWER
# ==============================================================================

cat("\n>> 2 -- Soft-Thresholding Power\n")

# WGCNA 1.74 / R 4.5+ compatibility: WGCNA::cor has extra args (weights.x,
# cosine) that conflict with stats::cor dispatch inside blockwiseModules.
# Standard workaround: temporarily override cor in the global environment.
cor <- WGCNA::cor

powers <- c(1:20)
sft <- pickSoftThreshold(datExpr, powerVector = powers,
                          networkType = "signed", verbose = 2)

# Find first power with R^2 > 0.85
r2_values <- -sign(sft$fitIndices$slope) * sft$fitIndices$SFT.R.sq
power_idx <- which(r2_values > 0.85)[1]
soft_power <- if (!is.na(power_idx)) powers[power_idx] else 6
cat(sprintf("   Selected soft-thresholding power: %d (R^2 = %.3f)\n",
            soft_power, r2_values[soft_power]))

# Plot
pdf(file.path(REPORT_DIR, "01_soft_threshold.pdf"), width = 10, height = 5)
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
dev.off()

cat("   Saved: 01_soft_threshold.pdf\n")

# ==============================================================================
# 3.  NETWORK CONSTRUCTION
# ==============================================================================

cat("\n>> 3 -- Network Construction\n")

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

cat(sprintf("   Modules detected: %d (+ grey/unassigned)\n", n_modules))
cat("   Module sizes:\n")
print(table(module_colors))

# ==============================================================================
# 4.  MODULE DENDROGRAM
# ==============================================================================

cat("\n>> 4 -- Module Dendrogram\n")

pdf(file.path(REPORT_DIR, "02_module_dendrogram.pdf"), width = 12, height = 6)
plotDendroAndColors(net$dendrograms[[1]], module_colors[net$blockGenes[[1]]],
                    "Module colors",
                    main = "Protein dendrogram and module colors",
                    dendroLabels = FALSE, hang = 0.03, addGuide = TRUE,
                    guideHang = 0.05)
dev.off()

cat("   Saved: 02_module_dendrogram.pdf\n")

# ==============================================================================
# 5.  ENHANCED TRAIT MATRIX (9 traits)
# ==============================================================================

cat("\n>> 5 -- Enhanced Trait Matrix\n")

# Binary/categorical design factors
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

# Match sample order to datExpr
traits_mat <- traits_mat[rownames(datExpr), ]

cat(sprintf("   Trait matrix: %d samples x %d traits\n",
            nrow(traits_mat), ncol(traits_mat)))
cat(sprintf("   Traits: %s\n", paste(trait_cols, collapse = ", ")))

# ==============================================================================
# 6.  MODULE-TRAIT CORRELATIONS + BH P-VALUES
# ==============================================================================

cat("\n>> 6 -- Module-Trait Correlations\n")

# Module eigengenes
MEs <- moduleEigengenes(datExpr, colors = module_colors)$eigengenes
MEs <- orderMEs(MEs)

# Correlations and raw p-values
module_trait_cor  <- cor(MEs, traits_mat, use = "pairwise.complete.obs")
module_trait_pval <- corPvalueStudent(module_trait_cor, nrow(datExpr))

# BH correction across all tests
pval_vec <- as.vector(module_trait_pval)
pval_bh_vec <- p.adjust(pval_vec, method = "BH")
module_trait_pval_bh <- matrix(pval_bh_vec, nrow = nrow(module_trait_pval),
                                ncol = ncol(module_trait_pval))
rownames(module_trait_pval_bh) <- rownames(module_trait_pval)
colnames(module_trait_pval_bh) <- colnames(module_trait_pval)

# Heatmap with stars for BH-significant
star_matrix <- ifelse(module_trait_pval_bh < 0.001, "***",
               ifelse(module_trait_pval_bh < 0.01, "**",
               ifelse(module_trait_pval_bh < 0.05, "*", "")))
text_matrix <- paste(signif(module_trait_cor, 2), star_matrix, sep = "\n")
dim(text_matrix) <- dim(module_trait_cor)

pdf(file.path(REPORT_DIR, "03_module_trait_heatmap.pdf"), width = 10, height = 10)
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
dev.off()

# Export correlation matrix
trait_cor_df <- as.data.frame(module_trait_cor) %>%
  rownames_to_column("module")
write_csv(trait_cor_df, file.path(DATA_DIR, "wgcna_module_trait_correlations.csv"))

# Export BH p-values
pval_bh_df <- as.data.frame(module_trait_pval_bh) %>%
  rownames_to_column("module")
write_csv(pval_bh_df, file.path(DATA_DIR, "wgcna_module_trait_pvalues_bh.csv"))

cat("   Saved: 03_module_trait_heatmap.pdf\n")
cat("   Saved: wgcna_module_trait_correlations.csv, wgcna_module_trait_pvalues_bh.csv\n")

# ==============================================================================
# 7.  HUB PROTEINS
# ==============================================================================

cat("\n>> 7 -- Hub Proteins\n")

kME <- signedKME(datExpr, MEs)

# Module assignments
module_df <- tibble(
  uniprot_id   = colnames(datExpr),
  module_color = module_colors,
  module_num   = net$colors
) %>% left_join(ann %>% dplyr::select(uniprot_id, gene), by = "uniprot_id")

# Top 10 hub proteins per module (by |kME|)
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
cat(sprintf("   Hub proteins: %d across %d modules\n", nrow(hub_df), length(unique_modules)))

# Restore stats::cor now that all WGCNA-specific computations are done
cor <- stats::cor

# ==============================================================================
# 8.  GO:BP ENRICHMENT PER MODULE
# ==============================================================================

cat("\n>> 8 -- GO:BP Enrichment per Module\n")

# Background gene universe
bg_genes <- ann$gene[ann$uniprot_id %in% colnames(datExpr)]
bg_genes <- unique(bg_genes[!is.na(bg_genes) & bg_genes != ""])

go_results_list <- list()

for (mod in unique_modules) {
  mod_genes <- module_df$gene[module_df$module_color == mod]
  mod_genes <- unique(mod_genes[!is.na(mod_genes) & mod_genes != ""])

  if (length(mod_genes) < 5) {
    cat(sprintf("   Skipping module '%s' (only %d genes)\n", mod, length(mod_genes)))
    next
  }

  ego <- tryCatch(
    enrichGO(
      gene         = mod_genes,
      universe     = bg_genes,
      OrgDb        = org.Hs.eg.db,
      keyType      = "SYMBOL",
      ont          = "BP",
      pAdjustMethod = "BH",
      pvalueCutoff = 0.05,
      qvalueCutoff = 0.2,
      minGSSize    = 10,
      maxGSSize    = 500,
      readable     = FALSE
    ),
    error = function(e) {
      cat(sprintf("   Warning: GO enrichment failed for module '%s': %s\n", mod, e$message))
      NULL
    }
  )

  if (!is.null(ego) && nrow(as.data.frame(ego)) > 0) {
    ego_df <- as.data.frame(ego) %>%
      mutate(ONTOLOGY = "BP", module = mod)
    go_results_list <- c(go_results_list, list(ego_df))
    cat(sprintf("   Module '%s': %d significant GO:BP terms\n", mod, nrow(ego_df)))
  } else {
    cat(sprintf("   Module '%s': no significant GO:BP terms\n", mod))
  }
}

go_df <- bind_rows(go_results_list)
cat(sprintf("   Total GO:BP terms across all modules: %d\n", nrow(go_df)))

# ==============================================================================
# 9.  EXPORT
# ==============================================================================

cat("\n>> 9 -- Export\n")

write_csv(module_df, file.path(DATA_DIR, "wgcna_module_assignments.csv"))
write_csv(hub_df,    file.path(DATA_DIR, "wgcna_hub_proteins.csv"))
saveRDS(net,         file.path(DATA_DIR, "wgcna_network.rds"))
write_csv(go_df,     file.path(DATA_DIR, "wgcna_module_GO_enrichment.csv"))

cat("   Saved: wgcna_module_assignments.csv\n")
cat("   Saved: wgcna_hub_proteins.csv\n")
cat("   Saved: wgcna_network.rds\n")
cat("   Saved: wgcna_module_trait_correlations.csv\n")
cat("   Saved: wgcna_module_trait_pvalues_bh.csv\n")
cat("   Saved: wgcna_module_GO_enrichment.csv\n")

cat(sprintf("\nSummary: %d modules, %d hub proteins, %d GO terms\n",
            n_modules, nrow(hub_df), nrow(go_df)))

cat("\n")
sessionInfo()

cat("\n=== YvO WGCNA Runner Complete ===\n")
