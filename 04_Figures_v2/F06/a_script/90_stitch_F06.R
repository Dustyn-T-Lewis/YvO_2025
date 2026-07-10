# F06 master orchestrator: run YvO_WGCNA_run.R separately first.



# supp panels first: 01 assembles the workbook from their 03_/04_/05_ CSVs
source(here::here("04_Figures_v2", "F06", "a_script", "02_supp_panels.R"))
source(here::here("04_Figures_v2", "F06", "a_script", "01_main_panels.R"))

message("F06 complete: composites in 04_Figures_v2/F06/b_reports, table in 04_Figures_v2/F06/c_data")
