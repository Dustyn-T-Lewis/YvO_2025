# Supplementary Enrichment Gallery — UpSet Intersection Plot (F04: Blunting)
# Shows sig pathway overlaps across training contrasts using UpSetR.
setwd(rprojroot::find_rstudio_root_file())
source("04_Figures/shared/style.R")

suppressPackageStartupMessages({
  library(tidyverse)
  library(UpSetR)
})

RPT <- "04_Figures/F04/b_reports/supp"
DAT <- "04_Figures/F04/c_data/supp"
dir.create(RPT, recursive = TRUE, showWarnings = FALSE)

blunt <- readRDS(file.path(DAT, "prep_blunting.rds"))
blunt_upset <- blunt %>%
  transmute(
    pathway,
    `Training (Young)` = as.integer(sig_TY),
    `Training (Old)`   = as.integer(sig_TO),
    Interaction        = as.integer(sig_Int)
  ) %>%
  as.data.frame()

pdf(file.path(RPT, "a_upset.pdf"), width = 200/25.4, height = 140/25.4)
upset(
  blunt_upset,
  sets          = c("Training (Young)", "Training (Old)", "Interaction"),
  order.by      = "freq",
  keep.order    = TRUE,
  mb.ratio      = c(0.6, 0.4),
  sets.bar.color = unname(CONTRAST_COLORS[c("Training_Young", "Training_Old", "Interaction")]),
  main.bar.color = "grey30",
  text.scale    = c(1.3, 1, 1, 1, 1.2, 1),
  mainbar.y.label = "Intersection size\n(Training Blunting)",
  sets.x.label    = "Set size"
)
dev.off()

png(file.path(RPT, "a_upset.png"), width = 200, height = 140,
    units = "mm", res = 300)
upset(
  blunt_upset,
  sets          = c("Training (Young)", "Training (Old)", "Interaction"),
  order.by      = "freq",
  keep.order    = TRUE,
  mb.ratio      = c(0.6, 0.4),
  sets.bar.color = unname(CONTRAST_COLORS[c("Training_Young", "Training_Old", "Interaction")]),
  main.bar.color = "grey30",
  text.scale    = c(1.3, 1, 1, 1, 1.2, 1),
  mainbar.y.label = "Intersection size\n(Training Blunting)",
  sets.x.label    = "Set size"
)
dev.off()

cat("UpSet plot saved.\n")
