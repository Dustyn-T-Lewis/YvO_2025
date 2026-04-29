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
file.copy(file.path(RPT, "supp/pdf/SUPP_F06_composite.pdf"),
          file.path(BOX, "02_Figures/F06_WGCNA/supp/S5_Figure.pdf"), overwrite = TRUE)
file.copy(file.path(RPT, "supp/png/SUPP_F06_composite.png"),
          file.path(BOX, "02_Figures/F06_WGCNA/supp/S5_Figure.png"), overwrite = TRUE)
# Copy individual module panels (hub networks + triptychs) + journal-named copies
box_modules <- file.path(BOX, "02_Figures/F06_WGCNA/supp/modules")
dir.create(box_modules, recursive = TRUE, showWarnings = FALSE)
mod_pdfs <- list.files(file.path(RPT, "supp/pdf/modules"),
                       pattern = "^SUPP_(hub|triptych|networks)_.*\\.(pdf)$", full.names = TRUE)
mod_pngs <- list.files(file.path(RPT, "supp/png/modules"),
                       pattern = "^SUPP_(hub|triptych|networks)_.*\\.(png)$", full.names = TRUE)
# triptych_<color>_<name> -> S6_Figure_<color>_<name>; hub_<color>_<name> -> S7_Figure_hub_<color>_<name>;
# networks_composite -> S7_Figure
journal_name <- function(orig) {
  base <- basename(orig)
  if (startsWith(base, "SUPP_triptych_"))   return(sub("^SUPP_triptych_", "S6_Figure_",   base))
  if (startsWith(base, "SUPP_hub_"))        return(sub("^SUPP_hub_",      "S7_Figure_hub_", base))
  if (startsWith(base, "SUPP_networks_composite.")) return(sub("^SUPP_networks_composite\\.", "S7_Figure.", base))
  base
}
for (f in c(mod_pdfs, mod_pngs)) {
  file.copy(f, file.path(box_modules, basename(f)), overwrite = TRUE)
  file.copy(f, file.path(box_modules, journal_name(f)), overwrite = TRUE)
}
file.copy(file.path(BASE, "c_data/F06_supplementary.xlsx"),
          file.path(BOX, "02_Figures/F06_WGCNA/F06_source_data.xlsx"), overwrite = TRUE)
file.copy(file.path(BASE, "c_data/F06_supplementary.xlsx"),
          file.path(BOX, "03_Supplementary/S10_F06_WGCNA.xlsx"), overwrite = TRUE)

message("F06 complete")
