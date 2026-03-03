# YvO_F1_setup.R — Shared setup for Figure 1
# Sources: 05_Figures/shared/style.R
# Provides: norm_df, imp_df, dep_df, meta, samp_names, imp_mat,
#           RPT_DIR, DAT_DIR, BEST_IMP_METHOD, CONTRASTS

setwd(rprojroot::find_rstudio_root_file())

suppressPackageStartupMessages({
  library(tidyverse)
  library(patchwork)
  library(ComplexHeatmap)
  library(fgsea)
  library(msigdbr)
  library(rrvgo)
  library(GOSemSim)
  library(org.Hs.eg.db)
  library(grid)
  library(ggsignif)
  library(vegan)
})

source("05_Figures/shared/style.R")

# ── Paths ──────────────────────────────────────────────────────────────────────
NORM_FILE <- "01_normalization/c_data/02_normalized.csv"
IMP_FILE  <- "02_Imputation/c_data/01_imputed.csv"
DEP_FILE  <- "03_DEP/c_data/03_combined_results.csv"

RPT_DIR <- "05_Figures/F1/b_reports"
DAT_DIR <- "05_Figures/F1/c_data"
dir.create(RPT_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(DAT_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(RPT_DIR, "supplementary"), showWarnings = FALSE)

CONTRASTS <- c("Aging", "Training_Young", "Training_Old", "Interaction")

# ── Load data ──────────────────────────────────────────────────────────────────
norm_df <- read_csv(NORM_FILE, show_col_types = FALSE)
imp_df  <- read_csv(IMP_FILE,  show_col_types = FALSE)
dep_df  <- read_csv(DEP_FILE,  show_col_types = FALSE)

ann_cols   <- c("uniprot_id", "protein", "gene", "description")
samp_names <- setdiff(names(norm_df), ann_cols)

meta <- tibble(sample_id = samp_names) |>
  mutate(
    prefix   = str_extract(sample_id, "^[A-Z]+"),
    subj_num = str_extract(sample_id, "S\\d+"),
    time     = str_extract(sample_id, "(Pre|Post)$"),
    age      = if_else(str_detect(prefix, "^O"), "Old", "Young"),
    subject  = paste0(prefix, "_", subj_num),
    group    = paste(age, time, sep = "_")
  )
meta$age   <- factor(meta$age,  levels = c("Young", "Old"))
meta$time  <- factor(meta$time, levels = c("Pre", "Post"))
meta$group <- factor(meta$group,
                     levels = c("Young_Pre", "Young_Post", "Old_Pre", "Old_Post"))

cat(sprintf("Loaded: %d proteins, %d samples, %d DEP rows\n",
            nrow(norm_df), length(samp_names), nrow(dep_df)))

# Winning imputation method
imp_summary <- readLines("02_Imputation/c_data/09_imputation_summary.txt")
BEST_IMP_METHOD <- toupper(trimws(sub(".*=\\s*", "", grep("^best_method", imp_summary, value = TRUE))))
if (length(BEST_IMP_METHOD) == 0) BEST_IMP_METHOD <- "IMPUTED"

# Imputed matrix (used by Panel C for PCA)
imp_mat <- as.matrix(imp_df[, samp_names])
rownames(imp_mat) <- imp_df$gene
