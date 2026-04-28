# F06 master orchestrator — run YvO_WGCNA_run.R separately first

library(here)

source(here::here("04_Figures", "F06", "a_script", "02_supp_panels.R"))
source(here::here("04_Figures", "F06", "a_script", "01_main_panels.R"))

BASE <- here::here("04_Figures", "F06")
BOX <- "/Users/dtl0018/Library/CloudStorage/Box-Box/YvO_proteomics_manuscript"
RPT <- file.path(BASE, "b_reports")

file.copy(file.path(RPT, "main/pdf/MAIN_F06_composite.pdf"),
          file.path(BOX, "02_Figures/F06_WGCNA/main/MAIN_F06_composite.pdf"), overwrite = TRUE)
file.copy(file.path(RPT, "main/png/MAIN_F06_composite.png"),
          file.path(BOX, "02_Figures/F06_WGCNA/main/MAIN_F06_composite.png"), overwrite = TRUE)
file.copy(file.path(RPT, "supp/pdf/SUPP_F06_composite.pdf"),
          file.path(BOX, "02_Figures/F06_WGCNA/supp/SUPP_F06_composite.pdf"), overwrite = TRUE)
file.copy(file.path(RPT, "supp/png/SUPP_F06_composite.png"),
          file.path(BOX, "02_Figures/F06_WGCNA/supp/SUPP_F06_composite.png"), overwrite = TRUE)
# Copy individual module panels (hub networks + triptychs)
box_modules <- file.path(BOX, "02_Figures/F06_WGCNA/supp/modules")
dir.create(box_modules, recursive = TRUE, showWarnings = FALSE)
mod_pdfs <- list.files(file.path(RPT, "supp/pdf/modules"),
                       pattern = "^SUPP_(hub|triptych|networks)_.*\\.(pdf)$", full.names = TRUE)
mod_pngs <- list.files(file.path(RPT, "supp/png/modules"),
                       pattern = "^SUPP_(hub|triptych|networks)_.*\\.(png)$", full.names = TRUE)
for (f in c(mod_pdfs, mod_pngs)) file.copy(f, file.path(box_modules, basename(f)), overwrite = TRUE)
file.copy(file.path(BASE, "c_data/F06_supplementary.xlsx"),
          file.path(BOX, "02_Figures/F06_WGCNA/F06_source_data.xlsx"), overwrite = TRUE)
file.copy(file.path(BASE, "c_data/F06_supplementary.xlsx"),
          file.path(BOX, "03_Supplementary_Tables/S12_F06_WGCNA.xlsx"), overwrite = TRUE)

message("F06 complete")
