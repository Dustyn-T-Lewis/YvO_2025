#!/usr/bin/env Rscript
# F01 Supp — Deadlift 1RM (A) + Type II fCSA (B) + Type I fCSA (C)

setwd(rprojroot::find_root(rprojroot::has_file("setup.R")))

library(patchwork)
library(cowplot)

source("04_Figures/shared/style.R")

BASE    <- "04_Figures/F01"
RPT_PNG <- file.path(BASE, "b_reports", "supp", "png")
RPT_PDF <- file.path(BASE, "b_reports", "supp", "pdf")
PNL_PNG <- file.path(RPT_PNG, "panels")
PNL_PDF <- file.path(RPT_PDF, "panels")
DAT     <- file.path(BASE, "c_data", "supp")
for (d in c(PNL_PNG, PNL_PDF, DAT)) dir.create(d, recursive = TRUE, showWarnings = FALSE)

TMPL <- "04_Figures/F01/a_script/_prepost_template.R"

cfg <- list(
  dv_col = "deadlift_1rm_kg", y_label = "Deadlift 1RM (kg)",
  delta_label = expression(bold(Delta ~ "1RM (kg)")),
  title = "Deadlift 1RM", tag = "a", output_prefix = "pSA",
  file_tag = "panel_A_deadlift_1rm", audit_file = "panel_A_deadlift_1rm.csv",
  file_prefix = "SUPP", rpt_png = PNL_PNG, rpt_pdf = PNL_PDF, dat = DAT,
  coerce_cols = TRUE, filter_complete = TRUE)
source(TMPL)

cfg <- list(
  dv_col = "Type_II_fCSA",
  y_label = expression(bold("Type II fCSA (" * mu * m^2 * ")")),
  delta_label = expression(bold(Delta ~ "Type II fCSA (" * mu * m^2 * ")")),
  title = "Type II Fiber CSA", tag = "b", output_prefix = "pSB",
  file_tag = "panel_B_type_II_fcsa", audit_file = "panel_B_type_II_fcsa.csv",
  file_prefix = "SUPP", rpt_png = PNL_PNG, rpt_pdf = PNL_PDF, dat = DAT,
  coerce_cols = TRUE, filter_complete = TRUE)
source(TMPL)

cfg <- list(
  dv_col = "Type_I_fCSA",
  y_label = expression(bold("Type I fCSA (" * mu * m^2 * ")")),
  delta_label = expression(bold(Delta ~ "Type I fCSA (" * mu * m^2 * ")")),
  title = "Type I Fiber CSA", tag = "c", output_prefix = "pSC",
  file_tag = "panel_C_type_I_fcsa", audit_file = "panel_C_type_I_fcsa.csv",
  file_prefix = "SUPP", rpt_png = PNL_PNG, rpt_pdf = PNL_PDF, dat = DAT,
  coerce_cols = TRUE, filter_complete = TRUE)
source(TMPL)

composite <- (pSA_left + pSA_right +
              pSB_left + pSB_right +
              pSC_left + pSC_right) +
  plot_layout(ncol = 2, widths = c(0.65, 0.35), heights = c(1, 1, 1)) &
  theme(plot.margin = margin(6, 2, 2, 2))

COMP_W <- 85; COMP_H <- 88
txt <- composite_text_sizes(COMP_H)
TAG_SZ <- txt$tag - 3; TTL_SZ <- txt$title - 2; SUB_SZ <- txt$subtitle - 1
Y_A <- 0.983; Y_B <- 0.660; Y_C <- 0.336
X_TAG <- 0.056; X_TTL <- 0.116; SUB_OFF <- 0.016

composite <- ggdraw(composite) +
  draw_label("A", x = X_TAG, y = Y_A, size = TAG_SZ, fontface = "bold", hjust = 0, vjust = 1) +
  draw_label(pSA_title, x = X_TTL, y = Y_A, size = TTL_SZ, fontface = "bold", hjust = 0, vjust = 1) +
  draw_label(pSA_subtitle, x = X_TTL, y = Y_A - SUB_OFF, size = SUB_SZ, fontface = "bold.italic", hjust = 0, vjust = 1, colour = "grey30") +
  draw_label("B", x = X_TAG, y = Y_B, size = TAG_SZ, fontface = "bold", hjust = 0, vjust = 1) +
  draw_label(pSB_title, x = X_TTL, y = Y_B, size = TTL_SZ, fontface = "bold", hjust = 0, vjust = 1) +
  draw_label(pSB_subtitle, x = X_TTL, y = Y_B - SUB_OFF, size = SUB_SZ, fontface = "bold.italic", hjust = 0, vjust = 1, colour = "grey30") +
  draw_label("C", x = X_TAG, y = Y_C, size = TAG_SZ, fontface = "bold", hjust = 0, vjust = 1) +
  draw_label(pSC_title, x = X_TTL, y = Y_C, size = TTL_SZ, fontface = "bold", hjust = 0, vjust = 1) +
  draw_label(pSC_subtitle, x = X_TTL, y = Y_C - SUB_OFF, size = SUB_SZ, fontface = "bold.italic", hjust = 0, vjust = 1, colour = "grey30")

ggsave(file.path(RPT_PDF, "SUPP_F01_composite.pdf"), composite,
       width = COMP_W, height = COMP_H, units = "mm", device = get_pdf_device())
ggsave(file.path(RPT_PNG, "SUPP_F01_composite.png"), composite,
       width = COMP_W, height = COMP_H, units = "mm", dpi = 300)

message("F01 supp composite done")
