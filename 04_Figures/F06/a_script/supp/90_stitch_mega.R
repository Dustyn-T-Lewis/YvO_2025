#!/usr/bin/env Rscript
# F06 SUPP — Mega orchestrator for supplementary outputs
#
# Produces:
#   supp/{pdf,png}/SUPP_F06_composite.{pdf,png}   (QC composite)
#
# Per-module triptychs and hub networks in supp/modules/ are generated
# by their own scripts (panel_B_triptych.R, panel_D_hub.R) and are NOT
# assembled into a composite — they serve as pathway-naming validation
# reference material and workbook-ready panels.

setwd(rprojroot::find_rstudio_root_file())

cat("=== F06 SUPP: Running QC composite ===\n")
source("04_Figures/F06/a_script/supp/01_QC/90_stitch_QC.R")

cat("=== F06 SUPP mega complete ===\n")
