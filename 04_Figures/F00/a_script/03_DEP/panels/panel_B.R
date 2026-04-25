# F00 Stage 03 Panel B: Effect-size bootstrap forest
setwd(rprojroot::find_rstudio_root_file())
source("04_Figures/shared/style.R")

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(ggplot2)
})

PW <- 100; PH <- 80
RPT_PNG <- "04_Figures/F00/b_reports/03_DEP/png/panels"
RPT_PDF <- "04_Figures/F00/b_reports/03_DEP/pdf/panels"
DAT     <- "04_Figures/F00/c_data/03_DEP"
dir.create(RPT_PNG, recursive = TRUE, showWarnings = FALSE)
dir.create(RPT_PDF, recursive = TRUE, showWarnings = FALSE)
dir.create(DAT, recursive = TRUE, showWarnings = FALSE)

# CONTRAST_ORDER defined in shared/style.R

if (!exists("eff"))
  eff <- read_csv("03_DEP/c_data/07_effect_size_bootstrap.csv", show_col_types = FALSE)

eff2 <- eff %>%
  mutate(contrast = factor(contrast,
                           levels = rev(intersect(CONTRAST_ORDER, unique(contrast)))))

write.csv(eff2, file.path(DAT, "panel_B_effect_size.csv"), row.names = FALSE)

pB_s03 <- ggplot(eff2, aes(median_absLFC, contrast, color = contrast)) +
  geom_errorbar(aes(xmin = ci_lower, xmax = ci_upper),
                width = 0.18, linewidth = 0.6, orientation = "y") +
  geom_point(size = 3) +
  geom_text(aes(label = sprintf("%.3f", median_absLFC)),
            vjust = -1.0, size = 2.6, color = "grey20") +
  scale_color_manual(values = CONTRAST_COLORS, guide = "none") +
  scale_x_continuous(expand = expansion(mult = c(0.05, 0.10))) +
  labs(title = "Effect-size bootstrap",
       subtitle = sprintf("median |logFC| with 95%% CI (n=%s)",
                          format(as.integer(median(eff$n_proteins, na.rm = TRUE)),
                                 big.mark = ",")),
       x = "median |logFC| (log2)", y = NULL, tag = "B") +
  FIG_THEME +
  theme(axis.text.y = element_text(size = 7))

ggsave(file.path(RPT_PNG, "panel_B_effect_size.png"), pB_s03,
       width = PW, height = PH, units = "mm", dpi = 300)
ggsave(file.path(RPT_PDF, "panel_B_effect_size.pdf"), pB_s03,
       width = PW, height = PH, units = "mm", device = get_pdf_device())
cat("F00 Stage03 Panel B done\n")
