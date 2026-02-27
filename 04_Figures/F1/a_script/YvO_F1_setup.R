# YvO_F1_setup.R — Shared setup for Figure 1
# Sources: palettes.R
# Provides: norm_df, imp_df, dep_df, meta, samp_names, imp_mat,
#           GROUP_FILL, CONTRAST_COLORS, DIR_COLORS, DB_COLORS,
#           THEME_PUB, RPT_DIR, DAT_DIR, BEST_IMP_METHOD

# 0. SETUP ====================================================================

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

setwd(rprojroot::find_rstudio_root_file())

NORM_FILE <- "01_normalization/c_data/02_normalized.csv"
IMP_FILE  <- "02_Imputation/c_data/01_imputed.csv"
DEP_FILE  <- "03_DEP/c_data/03_combined_results.csv"

RPT_DIR <- "04_Figures/F1/b_reports"
DAT_DIR <- "04_Figures/F1/c_data"
dir.create(RPT_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(DAT_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(RPT_DIR, "supplementary"), showWarnings = FALSE)

CONTRASTS <- c("Aging", "Training_Young", "Training_Old", "Interaction")
CTR_AXIS  <- c(Aging = "Aging", Training_Young = "Training\n(Young)",
               Training_Old = "Training\n(Old)", Interaction = "Interaction")
CTR_FACET <- c(Aging = "Aging", Training_Young = "Training (Young)",
               Training_Old = "Training (Old)", Interaction = "Interaction")
CTR_SHORT <- c(Aging = "Aging", Training_Young = "Tr. (Y)",
               Training_Old = "Tr. (O)", Interaction = "Inter.")

# Palettes
GROUP_FILL <- c(Young_Pre = alpha("#4393C3", 0.5), Young_Post = "#4393C3",
                Old_Pre   = alpha("#D6604D", 0.5), Old_Post   = "#D6604D")
CONTRAST_COLORS <- c(Aging = "#4CAF50", Training_Young = "#E05A4E",
                     Training_Old = "#5DA5DA", Interaction = "#9B7FBF")
DIR_COLORS <- c(Up = "#D6604D", Down = "#4393C3")
DB_COLORS  <- c(Hallmark = "#AA336A", "GO:BP" = "#00796B",
                "GO:CC" = "#26A69A", "GO:MF" = "#CD5C5C")

# PCA: 4 distinct group colors + shapes (circle=Pre, triangle=Post)
PCA_COLORS <- c(Young_Pre = "#93C4DE", Young_Post = "#2166AC",
                Old_Pre   = "#F4A582", Old_Post   = "#B2182B")
PCA_SHAPES <- c(Young_Pre = 16, Young_Post = 17, Old_Pre = 16, Old_Post = 17)

# Key constants centralized in palettes.R
source("04_Figures/shared/palettes.R")

THEME_PUB <- theme_bw(base_size = 8) +
  theme(plot.title       = element_text(face = "bold", size = 9),
        plot.subtitle    = element_text(size = 6.5, color = "grey30", face = "italic"),
        strip.background = element_blank(),
        strip.text       = element_text(face = "bold", size = 6.5),
        legend.key.size  = unit(3, "mm"))

# 1. LOAD DATA ================================================================

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

# Read winning imputation method from summary
imp_summary <- readLines("02_Imputation/c_data/09_imputation_summary.txt")
BEST_IMP_METHOD <- toupper(trimws(sub(".*=\\s*", "", grep("^best_method", imp_summary, value = TRUE))))
if (length(BEST_IMP_METHOD) == 0) BEST_IMP_METHOD <- "IMPUTED"

# Imputed matrix (used by Panel C for PCA, and potentially other panels)
imp_mat <- as.matrix(imp_df[, samp_names])
rownames(imp_mat) <- imp_df$gene
