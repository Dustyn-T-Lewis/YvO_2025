#!/usr/bin/env Rscript
# F02 Supp — CV scatter (A) + CV violin (B) + Imputed variability (C)
# Updated to read from DAList .rds and xlsx instead of deleted CSVs

library(dplyr)
library(tidyr)
library(stringr)
library(readr)
library(readxl)
library(ggplot2)
library(ggrepel)
library(ggbeeswarm)
library(patchwork)
library(cowplot)

source(here::here("04_Figures", "shared", "style.R"))

# F02-specific overrides
HEATMAP_LO <- "#2166AC"; HEATMAP_HI <- "#B2182B"
BASE_COUNT <- BASE_COUNT + 1.0
BASE_GENE  <- BASE_GENE  + 0.8
BASE_STAT  <- BASE_STAT  + 0.5

BASE    <- here::here("04_Figures", "F02")
RPT_PNG <- file.path(BASE, "b_reports", "supp", "png", "panels")
RPT_PDF <- file.path(BASE, "b_reports", "supp", "pdf", "panels")
DAT     <- file.path(BASE, "c_data")
for (d in c(RPT_PNG, RPT_PDF, DAT)) dir.create(d, recursive = TRUE, showWarnings = FALSE)

pdf_dev <- get_pdf_device()

# Load normalized data from DAList (replaces deleted 02_normalized.csv)
dal_norm <- readRDS(here::here("01_normalization", "c_data", "03_DAList_normalized.rds"))
norm_mat <- as.matrix(dal_norm$data)
norm_meta <- as_tibble(dal_norm$metadata) |>
  mutate(age     = factor(Group, levels = c("Young", "Old")),
         time    = factor(Timepoint, levels = c("Pre", "Post")),
         subject = sub("_(Pre|Post)$", "", Col_ID),
         group   = factor(Group_Time,
           levels = c("Young_Pre", "Young_Post", "Old_Pre", "Old_Post"))) |>
  rename(sample_id = Col_ID)

samp_names <- norm_meta$sample_id

# Annotation from DAList
ann_df <- as_tibble(dal_norm$annotation) |>
  select(uniprot_id, gene, protein, description)
norm_df <- bind_cols(ann_df, as_tibble(norm_mat))

# ══════════════════════════════════════════════════════════════════════════════
# Panel A: CV% scatter Pre vs Post by age
# ══════════════════════════════════════════════════════════════════════════════

source(here::here("04_Figures", "F02", "a_script", "_supp_A_cv_scatter.R"))

# ══════════════════════════════════════════════════════════════════════════════
# Panel B: CV% violin
# ══════════════════════════════════════════════════════════════════════════════

source(here::here("04_Figures", "F02", "a_script", "_supp_B_cv_violin.R"))

# ══════════════════════════════════════════════════════════════════════════════
# Panel C: Imputed variability
# ══════════════════════════════════════════════════════════════════════════════

source(here::here("04_Figures", "F02", "a_script", "_supp_C_imputed.R"))

# ══════════════════════════════════════════════════════════════════════════════
# Supp composite (A triptych full-width; B+C side by side)
# ══════════════════════════════════════════════════════════════════════════════

SUPP_PNG <- file.path(BASE, "b_reports", "supp", "png")
SUPP_PDF <- file.path(BASE, "b_reports", "supp", "pdf")
dir.create(SUPP_PNG, recursive = TRUE, showWarnings = FALSE)
dir.create(SUPP_PDF, recursive = TRUE, showWarnings = FALSE)

COMP_W <- 178; COMP_H <- 115
txt <- composite_text_sizes(COMP_H)

# pC12/pC3 from supp_C, pD from supp_D, pE from supp_E
pSA_full <- (pC12 | pC3) + plot_layout(widths = c(2, 1))
composite <- (wrap_elements(pSA_full) / (pD | pE)) +
  plot_layout(heights = c(1, 0.7))

Y_TOP <- 0.985; Y_BOT <- 0.500
composite <- ggdraw(composite) +
  draw_label("A", x = 0.01, y = Y_TOP, size = txt$tag, fontface = "bold", hjust = 0, vjust = 1) +
  draw_label(pSA_title, x = 0.04, y = Y_TOP, size = txt$title, fontface = "bold", hjust = 0, vjust = 1) +
  draw_label("B", x = 0.01, y = Y_BOT, size = txt$tag, fontface = "bold", hjust = 0, vjust = 1) +
  draw_label(pSB_title, x = 0.04, y = Y_BOT, size = txt$title, fontface = "bold", hjust = 0, vjust = 1) +
  draw_label("C", x = 0.52, y = Y_BOT, size = txt$tag, fontface = "bold", hjust = 0, vjust = 1) +
  draw_label(pSC_title, x = 0.55, y = Y_BOT, size = txt$title, fontface = "bold", hjust = 0, vjust = 1)

ggsave(file.path(SUPP_PDF, "SUPP_F02_composite.pdf"), composite,
       width = COMP_W, height = COMP_H, units = "mm", device = pdf_dev)
ggsave(file.path(SUPP_PNG, "SUPP_F02_composite.png"), composite,
       width = COMP_W, height = COMP_H, units = "mm", dpi = 300)

message("F02 supp composite done")
