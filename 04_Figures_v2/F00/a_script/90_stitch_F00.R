#!/usr/bin/env Rscript
# F00 — Pipeline QC: Orchestrator
# Sources the supp composite script
# F00 has no main figure — produces 2 supp composites:
#   SUPP_F00_normalization (panels A–G)
#   SUPP_F00_imputation    (panels H–N)

source(here::here("04_Figures_v2", "F00", "a_script", "01_supp_panels.R"))

message("F00 complete: composites in 04_Figures_v2/F00/b_reports, table in 04_Figures_v2/F00/c_data")
