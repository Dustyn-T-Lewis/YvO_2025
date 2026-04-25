# F03 SUPP — Composite stitch
# Sources 6 supplementary panels and composes a 3-row x 2-col grid:
#   Row 1: A (p-value histograms)    | B (Pi-score distributions)
#   Row 2: C (FDR distributions)     | D (MA plots)
#   Row 3: E (imputation sensitivity)| F (outlier sensitivity)
# No supra-title (consistent with other supplementary figures).
# Output: b_reports/supp/SUPP_F03_composite.{pdf,png} (178 x 225 mm)

setwd(rprojroot::find_rstudio_root_file())
source("04_Figures/shared/style.R")

suppressPackageStartupMessages({
  library(patchwork)
  library(ggplot2)
  library(cowplot)
})

source("04_Figures/F03/a_script/supp/panels/panel_B_phist.R")       # -> pB  (p-value histograms)
source("04_Figures/F03/a_script/supp/panels/panel_B2_pi_dist.R")    # -> pPi (Pi-score distributions)
source("04_Figures/F03/a_script/supp/panels/panel_C_fdr_dist.R")    # -> pFDR (FDR distributions)
source("04_Figures/F03/a_script/supp/panels/panel_C_ma.R")          # -> pC  (MA plots)
source("04_Figures/F03/a_script/supp/panels/panel_D_sensitivity.R") # -> pD  (imputation sensitivity)
source("04_Figures/F03/a_script/supp/panels/panel_E_outlier.R")     # -> pE  (outlier sensitivity)

RPT_PDF <- "04_Figures/F03/b_reports/supp/pdf"
RPT_PNG <- "04_Figures/F03/b_reports/supp/png"
dir.create(RPT_PDF, recursive = TRUE, showWarnings = FALSE)
dir.create(RPT_PNG, recursive = TRUE, showWarnings = FALSE)
pdf_device <- get_pdf_device()

COMP_W <- 178
COMP_H <- 225
CTS    <- composite_text_sizes(COMP_H)

# 3-row x 2-col grid — top margin on each panel for title group clearance
grid <- (pB | pPi) / (pFDR | pC) / (pD | pE) &
  theme(plot.margin = margin(9, 4, 4, 4))

# --- Tag + title + subtitle placement via cowplot ---
# Shifted right 0.5mm from prior positions
X_LEFT  <- 0.012     # right 0.5mm (was 0.009)
X_RIGHT <- 0.512     # right 0.5mm (was 0.509)
X_TTL   <- 0.029
SUB_OFFSET <- 0.014  # tight title-subtitle gap

# Row Y positions — R1 up 0.5mm (was 0.994)
Y_R1 <- 0.996
Y_R2 <- 0.673
Y_R3 <- 0.343

composite <- ggdraw(grid) +
  # --- Panel A (top-left): p-value histograms ---
  draw_label("A", x = X_LEFT, y = Y_R1,
             fontface = "bold", size = CTS$tag, hjust = 0, vjust = 1) +
  draw_label(pB_title, x = X_LEFT + X_TTL, y = Y_R1,
             fontface = "bold", size = CTS$title, hjust = 0, vjust = 1) +
  draw_label(pB_subtitle, x = X_LEFT + X_TTL, y = Y_R1 - SUB_OFFSET,
             fontface = "bold.italic", size = CTS$subtitle, colour = "grey30",
             hjust = 0, vjust = 1) +
  # --- Panel B (top-right): Pi-score distributions ---
  draw_label("B", x = X_RIGHT, y = Y_R1,
             fontface = "bold", size = CTS$tag, hjust = 0, vjust = 1) +
  draw_label(pPi_title, x = X_RIGHT + X_TTL, y = Y_R1,
             fontface = "bold", size = CTS$title, hjust = 0, vjust = 1) +
  draw_label(pPi_subtitle, x = X_RIGHT + X_TTL, y = Y_R1 - SUB_OFFSET,
             fontface = "bold.italic", size = CTS$subtitle, colour = "grey30",
             hjust = 0, vjust = 1) +
  # --- Panel C (mid-left): FDR distributions ---
  draw_label("C", x = X_LEFT, y = Y_R2,
             fontface = "bold", size = CTS$tag, hjust = 0, vjust = 1) +
  draw_label(pFDR_title, x = X_LEFT + X_TTL, y = Y_R2,
             fontface = "bold", size = CTS$title, hjust = 0, vjust = 1) +
  draw_label(pFDR_subtitle, x = X_LEFT + X_TTL, y = Y_R2 - SUB_OFFSET,
             fontface = "bold.italic", size = CTS$subtitle, colour = "grey30",
             hjust = 0, vjust = 1) +
  # --- Panel D (mid-right): MA plots ---
  draw_label("D", x = X_RIGHT, y = Y_R2,
             fontface = "bold", size = CTS$tag, hjust = 0, vjust = 1) +
  draw_label(pC_title, x = X_RIGHT + X_TTL, y = Y_R2,
             fontface = "bold", size = CTS$title, hjust = 0, vjust = 1) +
  draw_label(pC_subtitle, x = X_RIGHT + X_TTL, y = Y_R2 - SUB_OFFSET,
             fontface = "bold.italic", size = CTS$subtitle, colour = "grey30",
             hjust = 0, vjust = 1) +
  # --- Panel E (bottom-left): imputation sensitivity ---
  draw_label("E", x = X_LEFT, y = Y_R3,
             fontface = "bold", size = CTS$tag, hjust = 0, vjust = 1) +
  draw_label(pD_title, x = X_LEFT + X_TTL, y = Y_R3,
             fontface = "bold", size = CTS$title, hjust = 0, vjust = 1) +
  draw_label(pD_subtitle, x = X_LEFT + X_TTL, y = Y_R3 - SUB_OFFSET,
             fontface = "bold.italic", size = CTS$subtitle, colour = "grey30",
             hjust = 0, vjust = 1) +
  # --- Panel F (bottom-right): outlier sensitivity ---
  draw_label("F", x = X_RIGHT, y = Y_R3,
             fontface = "bold", size = CTS$tag, hjust = 0, vjust = 1) +
  draw_label(pE_title, x = X_RIGHT + X_TTL, y = Y_R3,
             fontface = "bold", size = CTS$title, hjust = 0, vjust = 1) +
  draw_label(pE_subtitle, x = X_RIGHT + X_TTL, y = Y_R3 - SUB_OFFSET,
             fontface = "bold.italic", size = CTS$subtitle, colour = "grey30",
             hjust = 0, vjust = 1)

ggsave(file.path(RPT_PDF, "SUPP_F03_composite.pdf"), composite,
       width = COMP_W, height = COMP_H, units = "mm", device = pdf_device)
ggsave(file.path(RPT_PNG, "SUPP_F03_composite.png"), composite,
       width = COMP_W, height = COMP_H, units = "mm", dpi = 300)

message(sprintf("F03 SUPP composite saved -> %s / %s", RPT_PDF, RPT_PNG))
