# YvO 2026 — install all R packages used across the pipeline.
# Run once after cloning: Rscript setup.R

cran <- c(
  "here", "rprojroot",
  "readxl", "readr", "openxlsx",
  "dplyr", "tidyr", "tibble", "stringr", "purrr", "tidyverse",
  "ggplot2", "ggrepel", "patchwork", "cowplot", "gridExtra", "scales", "shadowtext",
  "missForest", "missMDA", "boot", "pwr",
  "ggsignif", "rstatix", "ggbeeswarm", "vegan",
  "WGCNA", "lme4", "emmeans", "msigdbr",
  "pROC", "ggtext", "ppcor",
  "ggnewscale", "ggforce", "concaveman", "graphlayouts",
  "circlize", "igraph", "ggraph", "tidygraph",
  "colorspace", "png",
  # 02_imputation benchmark methods (optional; methods skip if unavailable)
  "mice", "sn", "imp4p", "magrittr", "abind",
  "glmnet", "foreach", "doParallel"
)

bioc <- c(
  "limma", "fgsea", "GO.db", "AnnotationDbi", "impute",
  "ComplexHeatmap", "MsCoreUtils", "msImpute",
  # 02_imputation benchmark methods
  "pcaMethods", "imputeLCMD"
)

github <- c(
  proteoDA = "ByrumLab/proteoDA",
  RRHO2    = "RRHO2/RRHO2",
  # Optional: DreamAI ensemble imputation (benchmark falls back if missing)
  DreamAI  = "WangLab-MSSM/DreamAI"
)

installed <- rownames(installed.packages())

need_cran <- setdiff(cran, installed)
if (length(need_cran)) install.packages(need_cran)

need_bioc <- setdiff(bioc, installed)
if (length(need_bioc)) {
  if (!"BiocManager" %in% installed) install.packages("BiocManager")
  BiocManager::install(need_bioc, ask = FALSE, update = FALSE)
}

for (pkg in names(github)) {
  if (!pkg %in% installed) {
    if (!"remotes" %in% installed) install.packages("remotes")
    remotes::install_github(github[[pkg]])
  }
}
