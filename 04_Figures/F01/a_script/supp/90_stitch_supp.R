# F01 — Supplementary Figure Composite (Panels A + B + C)
# Stacks three phenotype panels vertically:
#   A = Deadlift 1RM (was Panel D)
#   B = Type II Fiber CSA (was Panel E)
#   C = Type I Fiber CSA (was Panel F)

setwd(rprojroot::find_rstudio_root_file())

suppressPackageStartupMessages(library(patchwork))

source("04_Figures/F01/a_script/supp/panels/panel_D.R")
source("04_Figures/F01/a_script/supp/panels/panel_E.R")
source("04_Figures/F01/a_script/supp/panels/panel_F.R")

RPT_DIR <- "04_Figures/F01/b_reports/supp"
dir.create(RPT_DIR, recursive = TRUE, showWarnings = FALSE)

# Override tags: D/E/F -> A/B/C for supplementary figure
pD_left_a <- pD_left + labs(tag = "A")
pE_left_b <- pE_left + labs(tag = "B")
pF_left_c <- pF_left + labs(tag = "C")

# Trim margins to reduce inter-row gap
pD_left_a <- pD_left_a + theme(plot.margin = margin(5, 5, 0, 5))
pD_right  <- pD_right  + theme(plot.margin = margin(5, 5, 0, 5))
pE_left_b <- pE_left_b + theme(plot.margin = margin(0, 5, 0, 5))
pE_right  <- pE_right  + theme(plot.margin = margin(0, 5, 0, 5))
pF_left_c <- pF_left_c + theme(plot.margin = margin(0, 5, 5, 5))
pF_right  <- pF_right  + theme(plot.margin = margin(0, 5, 5, 5))

# Re-compose each row with updated tags (same widths as standalone panels)
row_A <- (pD_left_a | pD_right) + plot_layout(widths = c(0.65, 0.35))
row_B <- (pE_left_b | pE_right) + plot_layout(widths = c(0.65, 0.35))
row_C <- (pF_left_c | pF_right) + plot_layout(widths = c(0.65, 0.35))

# Stack vertically — wrap_elements needed for nested patchwork objects
composite <- wrap_elements(row_A) / wrap_elements(row_B) / wrap_elements(row_C) +
  plot_layout(heights = c(1, 1, 1))

COMP_W <- 170   # matches individual panel width
COMP_H <- 195   # tighter vertical spacing

ggsave(file.path(RPT_DIR, "SUPP_F01_composite.pdf"), composite,
       width = COMP_W, height = COMP_H, units = "mm", device = get_pdf_device())
ggsave(file.path(RPT_DIR, "SUPP_F01_composite.png"), composite,
       width = COMP_W, height = COMP_H, units = "mm", dpi = 300)

cat("F01 supp composite done\n")
