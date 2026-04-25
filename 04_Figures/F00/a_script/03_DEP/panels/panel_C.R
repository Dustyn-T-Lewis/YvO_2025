# F00 Stage 03 Panel C: Power analysis per contrast
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

if (!exists("pwr_dat"))
  pwr_dat <- read_csv("03_DEP/c_data/08_power_analysis.csv", show_col_types = FALSE)

pwr2 <- pwr_dat %>%
  mutate(contrast = factor(contrast,
                           levels = rev(intersect(CONTRAST_ORDER, unique(contrast)))))

write.csv(pwr2, file.path(DAT, "panel_C_power_analysis.csv"), row.names = FALSE)

pC_s03 <- ggplot(pwr2, aes(min_detectable_d, contrast, fill = contrast)) +
  geom_col(width = 0.65) +
  geom_text(aes(label = sprintf("d=%.2f, power=%.2f, n=%d",
                                min_detectable_d, power, n_subjects)),
            hjust = -0.05, size = 2.5, color = "grey20") +
  scale_fill_manual(values = CONTRAST_COLORS, guide = "none") +
  scale_x_continuous(expand = expansion(mult = c(0, 0.38))) +
  labs(title = "Power analysis",
       subtitle = sprintf("min detectable Cohen d | alpha=%.2f | power=%.2f",
                          median(pwr_dat$alpha, na.rm = TRUE),
                          median(pwr_dat$power, na.rm = TRUE)),
       x = "min detectable d", y = NULL, tag = "C") +
  FIG_THEME +
  theme(axis.text.y = element_text(size = 7))

ggsave(file.path(RPT_PNG, "panel_C_power_analysis.png"), pC_s03,
       width = PW, height = PH, units = "mm", dpi = 300)
ggsave(file.path(RPT_PDF, "panel_C_power_analysis.pdf"), pC_s03,
       width = PW, height = PH, units = "mm", device = get_pdf_device())
cat("F00 Stage03 Panel C done\n")
