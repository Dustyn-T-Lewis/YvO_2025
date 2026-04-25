# F00 Stage 02 Panel A: Per-protein %missing x MAR/MNAR classification
setwd(rprojroot::find_rstudio_root_file())
source("04_Figures/shared/style.R")

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
})

PW <- 100; PH <- 80
RPT_PNG <- "04_Figures/F00/b_reports/02_imputation/png/panels"
RPT_PDF <- "04_Figures/F00/b_reports/02_imputation/pdf/panels"
DAT     <- "04_Figures/F00/c_data/02_imputation"
dir.create(RPT_PNG, recursive = TRUE, showWarnings = FALSE)
dir.create(RPT_PDF, recursive = TRUE, showWarnings = FALSE)
dir.create(DAT, recursive = TRUE, showWarnings = FALSE)

if (!exists("int02"))
  int02 <- readRDS("02_Imputation/c_data/00_report_intermediates.rds")

mc <- int02$miss_class %>%
  mutate(classification = factor(classification,
                                 levels = c("Complete", "MAR", "MNAR")))

n_total <- nrow(mc)
class_counts <- mc %>%
  count(classification) %>%
  mutate(lab = sprintf("%s: %d (%.1f%%)", classification, n, 100*n/n_total))

write.csv(mc, file.path(DAT, "panel_A_miss_class_scatter.csv"), row.names = FALSE)

pA_s02 <- ggplot(mc, aes(mean_intensity, pct_miss, color = classification)) +
  geom_point(alpha = 0.45, size = 0.9) +
  scale_color_manual(values = PAL_CLASS, name = NULL,
                     labels = class_counts$lab) +
  scale_y_continuous(labels = function(x) paste0(x, "%")) +
  labs(title = "Per-protein missingness x class",
       subtitle = sprintf("%s proteins | k-means on mean intensity",
                          format(n_total, big.mark=",")),
       x = "Mean log2 intensity", y = "% missing", tag = "A") +
  FIG_THEME +
  theme(legend.position = "top",
        legend.key.size = unit(3, "mm"),       # F00: larger for QC readability
        legend.text = element_text(size = 7))  # F00: larger for QC readability

ggsave(file.path(RPT_PNG, "panel_A_miss_class.png"), pA_s02,
       width = PW, height = PH, units = "mm", dpi = 300)
ggsave(file.path(RPT_PDF, "panel_A_miss_class.pdf"), pA_s02,
       width = PW, height = PH, units = "mm", device = get_pdf_device())
cat("F00 Stage02 Panel A done\n")
