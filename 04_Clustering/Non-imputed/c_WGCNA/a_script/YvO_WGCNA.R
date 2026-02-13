################################################################################
#
#   YvO WGCNA — Non-imputed Data
#
#   Workflow:
#     1.  Setup & data loading, transpose, goodSamplesGenes check
#     2.  Soft-thresholding power (R^2 > 0.85)
#     3.  Network construction (blockwiseModules, signed)
#     4.  Module dendrogram
#     5.  Module-trait correlations
#     6.  Module eigengene boxplots by group
#     7.  Hub proteins (top 10 per module by kME)
#     8.  Export
#
#   Input: 01_normalization/c_data/01_normalized.csv (2124 proteins, NAs present)
#   Note: goodSamplesGenes filters proteins with >50% NA across samples.
#
################################################################################

# ==============================================================================
# 1.  SETUP & DATA LOADING
# ==============================================================================

cat("=== YvO WGCNA (Non-imputed) ===\n\n")
cat(">> 1 -- Setup & Data Loading\n")

suppressPackageStartupMessages({
  library(WGCNA)
  library(tidyverse)
  library(pheatmap)
  library(patchwork)
  library(ggrepel)
})

allowWGCNAThreads()

base_dir   <- normalizePath(file.path(dirname(getwd()), "..", "..", ".."), mustWork = TRUE)
DATA_FILE  <- file.path(base_dir, "01_normalization", "c_data", "01_normalized.csv")
REPORT_DIR <- file.path(base_dir, "04_Clustering", "Non-imputed", "c_WGCNA", "b_reports")
DATA_DIR   <- file.path(base_dir, "04_Clustering", "Non-imputed", "c_WGCNA", "c_data")

dir.create(REPORT_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(DATA_DIR,   recursive = TRUE, showWarnings = FALSE)

df <- read_csv(DATA_FILE, show_col_types = FALSE)
ann_cols   <- c("uniprot_id", "protein", "gene", "description")
ann        <- df[, ann_cols]
samp_names <- setdiff(names(df), ann_cols)
mat        <- as.matrix(df[, samp_names])
rownames(mat) <- ann$uniprot_id

n_missing <- sum(is.na(mat))
pct_missing <- round(100 * n_missing / length(mat), 1)
cat(sprintf("   Data: %d proteins x %d samples (%d NAs, %.1f%%)\n",
            nrow(mat), ncol(mat), n_missing, pct_missing))

# Transpose: samples as rows (WGCNA convention)
datExpr <- t(mat)

# Build metadata
meta <- tibble(sample_id = samp_names) %>%
  mutate(
    prefix   = str_extract(sample_id, "^[A-Z]+"),
    subj_num = str_extract(sample_id, "S\\d+"),
    time     = str_extract(sample_id, "(Pre|Post)$"),
    age      = if_else(str_detect(prefix, "^O"), "Old", "Young"),
    subject  = paste0(prefix, "_", subj_num),
    group    = paste(age, time, sep = "_")
  )

# goodSamplesGenes check — filters proteins with >50% NA
gsg <- goodSamplesGenes(datExpr, verbose = 3)
n_before <- ncol(datExpr)
if (!gsg$allOK) {
  datExpr <- datExpr[gsg$goodSamples, gsg$goodGenes]
  n_excluded <- n_before - ncol(datExpr)
  cat(sprintf("   goodSamplesGenes: excluded %d proteins (>50%% NA)\n", n_excluded))
  cat(sprintf("   After filtering: %d samples x %d proteins\n",
              nrow(datExpr), ncol(datExpr)))
} else {
  cat("   goodSamplesGenes: all OK\n")
}

# Update annotation to match filtered proteins
ann_filtered <- ann %>% filter(uniprot_id %in% colnames(datExpr))

# ==============================================================================
# 2.  SOFT-THRESHOLDING POWER
# ==============================================================================

cat("\n>> 2 -- Soft-Thresholding Power\n")

powers <- c(1:20)
sft <- pickSoftThreshold(datExpr, powerVector = powers,
                          networkType = "signed", verbose = 2)

r2_values <- -sign(sft$fitIndices$slope) * sft$fitIndices$SFT.R.sq
power_idx <- which(r2_values > 0.85)[1]
soft_power <- if (!is.na(power_idx)) powers[power_idx] else 6
cat(sprintf("   Selected soft-thresholding power: %d (R^2 = %.3f)\n",
            soft_power, r2_values[soft_power]))

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
  power            = soft_power,
  networkType      = "signed",
  TOMType          = "signed",
  minModuleSize    = 30,
  mergeCutHeight   = 0.25,
  numericLabels    = TRUE,
  pamRespectsDendro = FALSE,
  saveTOMs         = FALSE,
  verbose          = 3
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
                    main = "Protein dendrogram and module colors (Non-imputed)",
                    dendroLabels = FALSE, hang = 0.03, addGuide = TRUE,
                    guideHang = 0.05)
dev.off()

cat("   Saved: 02_module_dendrogram.pdf\n")

# ==============================================================================
# 5.  MODULE-TRAIT CORRELATIONS
# ==============================================================================

cat("\n>> 5 -- Module-Trait Correlations\n")

traits <- meta %>%
  mutate(
    age_num  = if_else(age == "Old", 1, 0),
    time_num = if_else(time == "Post", 1, 0),
    interaction = age_num * time_num
  ) %>%
  dplyr::select(age_num, time_num, interaction)
traits <- as.data.frame(traits)
rownames(traits) <- meta$sample_id

# Match sample order
traits <- traits[rownames(datExpr), ]

MEs <- moduleEigengenes(datExpr, colors = module_colors)$eigengenes
MEs <- orderMEs(MEs)

module_trait_cor <- cor(MEs, traits, use = "p")
module_trait_pval <- corPvalueStudent(module_trait_cor, nrow(datExpr))

text_matrix <- paste(signif(module_trait_cor, 2), "\n(",
                      signif(module_trait_pval, 1), ")", sep = "")
dim(text_matrix) <- dim(module_trait_cor)

pdf(file.path(REPORT_DIR, "03_module_trait_heatmap.pdf"), width = 8, height = 10)
labeledHeatmap(
  Matrix = module_trait_cor,
  xLabels = colnames(traits),
  yLabels = colnames(MEs),
  ySymbols = colnames(MEs),
  colorLabels = FALSE,
  colors = blueWhiteRed(50),
  textMatrix = text_matrix,
  setStdMargins = FALSE,
  cex.text = 0.6,
  zlim = c(-1, 1),
  main = "Module-trait correlations (Non-imputed)"
)
dev.off()

trait_cor_df <- as.data.frame(module_trait_cor) %>%
  rownames_to_column("module")
write_csv(trait_cor_df, file.path(DATA_DIR, "wgcna_module_trait_correlations.csv"))

cat("   Saved: 03_module_trait_heatmap.pdf, wgcna_module_trait_correlations.csv\n")

# ==============================================================================
# 6.  MODULE EIGENGENE BOXPLOTS
# ==============================================================================

cat("\n>> 6 -- Module Eigengene Boxplots\n")

ME_df <- MEs %>%
  rownames_to_column("sample_id") %>%
  left_join(meta %>% dplyr::select(sample_id, group), by = "sample_id") %>%
  pivot_longer(cols = starts_with("ME"), names_to = "module", values_to = "eigengene")

ME_df$group <- factor(ME_df$group,
                       levels = c("Young_Pre", "Young_Post", "Old_Pre", "Old_Post"))

p_me <- ggplot(ME_df, aes(x = group, y = eigengene, fill = group)) +
  geom_boxplot(outlier.size = 1) +
  facet_wrap(~ module, scales = "free_y") +
  scale_fill_manual(values = c("Young_Pre" = "#4A90D9", "Young_Post" = "#2E5FA1",
                                "Old_Pre" = "#D94A4A", "Old_Post" = "#A12E2E")) +
  labs(x = NULL, y = "Module Eigengene", title = "Module eigengenes by group (Non-imputed)") +
  theme_bw(base_size = 9) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "none")

ggsave(file.path(REPORT_DIR, "04_eigengene_boxplots.pdf"), p_me,
       width = 14, height = 10)

cat("   Saved: 04_eigengene_boxplots.pdf\n")

# ==============================================================================
# 7.  HUB PROTEINS
# ==============================================================================

cat("\n>> 7 -- Hub Proteins\n")

kME <- signedKME(datExpr, MEs)

module_df <- tibble(
  uniprot_id = colnames(datExpr),
  module_color = module_colors,
  module_num = net$colors
) %>% left_join(ann_filtered %>% dplyr::select(uniprot_id, gene), by = "uniprot_id")

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
    left_join(ann_filtered %>% dplyr::select(uniprot_id, gene), by = "uniprot_id")

  hub_rows <- c(hub_rows, list(mod_kme))
}

hub_df <- bind_rows(hub_rows)
cat(sprintf("   Hub proteins: %d across %d modules\n", nrow(hub_df), length(unique_modules)))

# ==============================================================================
# 8.  EXPORT
# ==============================================================================

cat("\n>> 8 -- Export\n")

write_csv(module_df, file.path(DATA_DIR, "wgcna_module_assignments.csv"))
write_csv(hub_df, file.path(DATA_DIR, "wgcna_hub_proteins.csv"))
saveRDS(net, file.path(DATA_DIR, "wgcna_network.rds"))

cat("   Saved: wgcna_module_assignments.csv, wgcna_hub_proteins.csv\n")
cat("   Saved: wgcna_module_trait_correlations.csv, wgcna_network.rds\n")
if (exists("n_excluded")) {
  cat(sprintf("   NOTE: %d proteins excluded by goodSamplesGenes (>50%% NA)\n", n_excluded))
}

cat("\n")
sessionInfo()

cat("\n=== YvO WGCNA (Non-imputed) Complete ===\n")
