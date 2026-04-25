# WGCNA F06 Supplementary — Bicor Sensitivity Analysis
#
# Compares Pearson (main analysis) and biweight midcorrelation (bicor) module
# assignments to confirm robustness of correlation choice.
#
# Rationale: Johnson et al. 2020 (PMID 32284590, Nat Med) and Johnson et al.
# 2022 (PMID 35115731, Nat Neurosci) use bicor as standard for large-scale
# proteomics WGCNA. Our main analysis uses Pearson because DIA-MS data has
# fewer outliers than TMT (see YvO_WGCNA_run.R methods note). This script
# provides empirical validation by showing module overlap between the two.
#
# Inputs:  02_Imputation/c_data/01_imputed.csv
#          04_Figures/F06/c_data/wgcna/wgcna_module_assignments.csv (Pearson)
# Outputs: c_data/supp/a02_bicor_sensitivity.csv
#          b_reports/supp/01_QC/SUPP_bicor_sensitivity.png

setwd(rprojroot::find_rstudio_root_file())
source("04_Figures/shared/style.R")

library(WGCNA)
library(tidyverse)

allowWGCNAThreads()
set.seed(42)

# --- Paths
DATA_FILE    <- "02_Imputation/c_data/01_imputed.csv"
PEARSON_MODS <- "04_Figures/F06/c_data/wgcna/wgcna_module_assignments.csv"
RPT_PNG      <- "04_Figures/F06/b_reports/supp/png/panels"
RPT_PDF      <- "04_Figures/F06/b_reports/supp/pdf/panels"
DAT_OUT      <- "04_Figures/F06/c_data/supp"

dir.create(RPT_PNG, recursive = TRUE, showWarnings = FALSE)
dir.create(RPT_PDF, recursive = TRUE, showWarnings = FALSE)
dir.create(DAT_OUT, recursive = TRUE, showWarnings = FALSE)

stopifnot(file.exists(DATA_FILE), file.exists(PEARSON_MODS))

# --- Load & preprocess (identical to YvO_WGCNA_run.R lines 29-61)
df <- read_csv(DATA_FILE, show_col_types = FALSE)
ann_cols   <- c("uniprot_id", "protein", "gene", "description")
samp_names <- setdiff(names(df), ann_cols)
mat        <- as.matrix(df[, samp_names])
rownames(mat) <- df$uniprot_id

datExpr <- t(mat)

gsg <- goodSamplesGenes(datExpr, verbose = 0)
if (!gsg$allOK) datExpr <- datExpr[gsg$goodSamples, gsg$goodGenes]

# WGCNA::cor dispatch workaround
cor <- WGCNA::cor

message(sprintf("Bicor sensitivity: %d samples x %d proteins", nrow(datExpr), ncol(datExpr)))

# --- Bicor soft threshold selection
powers <- 1:20
sft_bicor <- pickSoftThreshold(datExpr, powerVector = powers,
                                networkType = "signed",
                                corFnc = "bicor",
                                corOptions = list(maxPOutliers = 0.1),
                                verbose = 2)

r2_bicor <- -sign(sft_bicor$fitIndices$slope) * sft_bicor$fitIndices$SFT.R.sq
power_idx <- which(r2_bicor > 0.85)[1]
bicor_power <- if (!is.na(power_idx)) powers[power_idx] else 6

cat(sprintf("  Bicor soft power: %d (R^2 = %.3f)\n", bicor_power, r2_bicor[bicor_power]))

# --- Bicor network construction
net_bicor <- blockwiseModules(
  datExpr,
  power             = bicor_power,
  networkType       = "signed",
  TOMType           = "signed",
  corType           = "bicor",
  maxPOutliers      = 0.1,
  minModuleSize     = 30,
  mergeCutHeight    = 0.25,
  numericLabels     = TRUE,
  pamRespectsDendro = FALSE,
  saveTOMs          = FALSE,
  verbose           = 3
)

bicor_colors <- labels2colors(net_bicor$colors)
n_bicor <- length(unique(net_bicor$colors)) - (0 %in% net_bicor$colors)
cat(sprintf("  Bicor modules: %d (+ grey)\n", n_bicor))

# --- Load Pearson assignments
pearson_df <- read_csv(PEARSON_MODS, show_col_types = FALSE)

bicor_df <- tibble(
  uniprot_id   = colnames(datExpr),
  bicor_module = bicor_colors
)

# Merge
compare <- pearson_df %>%
  select(uniprot_id, pearson_module = module_color) %>%
  inner_join(bicor_df, by = "uniprot_id")

# --- Jaccard overlap matrix
pearson_mods <- setdiff(unique(compare$pearson_module), "grey")
bicor_mods   <- setdiff(unique(compare$bicor_module),   "grey")

jaccard_mat <- matrix(0, nrow = length(pearson_mods), ncol = length(bicor_mods),
                      dimnames = list(pearson_mods, bicor_mods))

for (pm in pearson_mods) {
  p_set <- compare$uniprot_id[compare$pearson_module == pm]
  for (bm in bicor_mods) {
    b_set <- compare$uniprot_id[compare$bicor_module == bm]
    intersection <- length(intersect(p_set, b_set))
    union_size   <- length(union(p_set, b_set))
    jaccard_mat[pm, bm] <- if (union_size > 0) intersection / union_size else 0
  }
}

# --- Greedy max-Jaccard matching
matched <- tibble(
  pearson_module = character(),
  bicor_module   = character(),
  jaccard        = numeric(),
  n_shared       = integer(),
  n_pearson      = integer(),
  n_bicor        = integer()
)

jac_temp <- jaccard_mat
used_p <- character()
used_b <- character()

while (TRUE) {
  remaining <- jac_temp[!rownames(jac_temp) %in% used_p,
                        !colnames(jac_temp) %in% used_b, drop = FALSE]
  if (nrow(remaining) == 0 || ncol(remaining) == 0) break
  best_idx <- which(remaining == max(remaining), arr.ind = TRUE)[1, ]
  best_pm  <- rownames(remaining)[best_idx[1]]
  best_bm  <- colnames(remaining)[best_idx[2]]

  p_set <- compare$uniprot_id[compare$pearson_module == best_pm]
  b_set <- compare$uniprot_id[compare$bicor_module   == best_bm]

  matched <- bind_rows(matched, tibble(
    pearson_module = best_pm,
    bicor_module   = best_bm,
    jaccard        = remaining[best_idx[1], best_idx[2]],
    n_shared       = length(intersect(p_set, b_set)),
    n_pearson      = length(p_set),
    n_bicor        = length(b_set)
  ))
  used_p <- c(used_p, best_pm)
  used_b <- c(used_b, best_bm)
}

matched <- matched %>% arrange(desc(jaccard))

write_csv(matched, file.path(DAT_OUT, "a02_bicor_sensitivity.csv"))
cat(sprintf("\n  Matched modules: %d\n", nrow(matched)))
cat(sprintf("  Mean Jaccard (matched): %.3f\n", mean(matched$jaccard)))
cat(sprintf("  Modules with Jaccard > 0.5: %d / %d\n",
            sum(matched$jaccard > 0.5), nrow(matched)))

# --- Heatmap visualization
jac_long <- as_tibble(jaccard_mat, rownames = "pearson_module") %>%
  pivot_longer(-pearson_module, names_to = "bicor_module", values_to = "jaccard")

# Order by matched pairs
pm_order <- c(matched$pearson_module,
              setdiff(pearson_mods, matched$pearson_module))
bm_order <- c(matched$bicor_module,
              setdiff(bicor_mods, matched$bicor_module))

jac_long <- jac_long %>%
  mutate(pearson_module = factor(pearson_module, levels = rev(pm_order)),
         bicor_module   = factor(bicor_module,   levels = bm_order))

p <- ggplot(jac_long, aes(bicor_module, pearson_module, fill = jaccard)) +
  geom_tile(colour = "white", linewidth = 0.3) +
  geom_text(aes(label = sprintf("%.2f", jaccard)),
            size = 2.8, colour = ifelse(jac_long$jaccard > 0.5, "white", "grey30")) +
  scale_fill_gradient2(low = "white", mid = "#B2DFDB", high = "#00695C",
                       midpoint = 0.35, limits = c(0, 1),
                       name = "Jaccard\nIndex") +
  labs(x = "Bicor module", y = "Pearson module",
       title = "Pearson vs Bicor Module Overlap",
       subtitle = sprintf("Mean matched Jaccard = %.3f  |  %d/%d modules > 0.5",
                          mean(matched$jaccard),
                          sum(matched$jaccard > 0.5), nrow(matched))) +
  FIG_THEME +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
        axis.text.y = element_text(size = 8),
        plot.title    = element_text(size = 13, face = "bold"),
        plot.subtitle = element_text(size = 10))

pdf_device <- get_pdf_device()

ggsave(file.path(RPT_PNG, "SUPP_bicor_sensitivity.png"), p,
       width = 180, height = 140, units = "mm", dpi = 300)
ggsave(file.path(RPT_PDF, "SUPP_bicor_sensitivity.pdf"), p,
       width = 180, height = 140, units = "mm", device = pdf_device)

cat("  Bicor sensitivity analysis complete.\n")
