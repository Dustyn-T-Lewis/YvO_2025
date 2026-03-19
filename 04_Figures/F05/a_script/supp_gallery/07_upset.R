# Supplementary Enrichment Gallery — UpSet Intersection Plot (F05: Reversal)
# Shows sig pathway overlaps across aging/training contrasts using UpSetR.
setwd(rprojroot::find_rstudio_root_file())
source("04_Figures/shared/style.R")

suppressPackageStartupMessages({
  library(tidyverse)
  library(UpSetR)
})

RPT <- "04_Figures/F05/b_reports/supp"
DAT <- "04_Figures/F05/c_data/supp"
dir.create(RPT, recursive = TRUE, showWarnings = FALSE)

rev <- readRDS(file.path(DAT, "prep_reversal.rds"))
rev_upset <- rev %>%
  transmute(
    pathway,
    Aging              = as.integer(sig_Aging),
    `Training (Old)`   = as.integer(sig_TO)
  ) %>%
  as.data.frame()

pdf(file.path(RPT, "a_upset.pdf"), width = 200/25.4, height = 140/25.4)
upset(
  rev_upset,
  sets          = c("Aging", "Training (Old)"),
  order.by      = "freq",
  keep.order    = TRUE,
  mb.ratio      = c(0.6, 0.4),
  sets.bar.color = unname(CONTRAST_COLORS[c("Aging", "Training_Old")]),
  main.bar.color = "grey30",
  text.scale    = c(1.3, 1, 1, 1, 1.2, 1),
  mainbar.y.label = "Intersection size\n(Aging Reversal)",
  sets.x.label    = "Set size"
)
dev.off()

png(file.path(RPT, "a_upset.png"), width = 200, height = 140,
    units = "mm", res = 300)
upset(
  rev_upset,
  sets          = c("Aging", "Training (Old)"),
  order.by      = "freq",
  keep.order    = TRUE,
  mb.ratio      = c(0.6, 0.4),
  sets.bar.color = unname(CONTRAST_COLORS[c("Aging", "Training_Old")]),
  main.bar.color = "grey30",
  text.scale    = c(1.3, 1, 1, 1, 1.2, 1),
  mainbar.y.label = "Intersection size\n(Aging Reversal)",
  sets.x.label    = "Set size"
)
dev.off()

cat("UpSet plot saved.\n")
