# F00 Stage 01 Panel D: Outlier consensus diagnostic
setwd(rprojroot::find_rstudio_root_file())
source("04_Figures/shared/style.R")

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
})

PW <- 100; PH <- 80
RPT_PNG <- "04_Figures/F00/b_reports/01_normalization/png/panels"
RPT_PDF <- "04_Figures/F00/b_reports/01_normalization/pdf/panels"
DAT     <- "04_Figures/F00/c_data/01_normalization"
dir.create(RPT_PNG, recursive = TRUE, showWarnings = FALSE)
dir.create(RPT_PDF, recursive = TRUE, showWarnings = FALSE)
dir.create(DAT, recursive = TRUE, showWarnings = FALSE)

if (!exists("int01"))
  int01 <- readRDS("01_normalization/c_data/00_report_intermediates.rds")

od <- int01$outlier_diag %>%
  mutate(age = ifelse(grepl("^O", prefix), "Old", "Young"),
         age = factor(age, levels = c("Young", "Old")),
         consensus_lab = ifelse(consensus_outlier, "Consensus outlier",
                                ifelse(n_flags > 0, "Single flag", "Clean")))

write.csv(od, file.path(DAT, "panel_D_outlier_consensus.csv"), row.names = FALSE)

pD_s01 <- ggplot(od, aes(mahal_dist, n_flags,
                     color = consensus_lab, shape = Timepoint)) +
  geom_jitter(width = 0.04, height = 0.12, size = 2.2, alpha = 0.85) +
  scale_color_manual(values = c("Clean" = "grey55",
                                "Single flag" = "#E6A100",
                                "Consensus outlier" = "#B2182B"),
                     name = NULL) +
  scale_shape_manual(values = SHAPE_TP, guide = "none") +
  scale_y_continuous(breaks = 0:4) +
  labs(title = "Outlier consensus diagnostic",
       subtitle = sprintf("%d samples | %d consensus outliers",
                          nrow(od), sum(od$consensus_outlier)),
       x = "Mahalanobis distance (PC1-3)",
       y = "Number of QC flags",
       tag = "D") +
  FIG_THEME +
  theme(legend.position = "top",
        legend.key.size = unit(3, "mm"),
        legend.text = element_text(size = 7))

ggsave(file.path(RPT_PNG, "panel_D_outlier_consensus.png"), pD_s01,
       width = PW, height = PH, units = "mm", dpi = 300)
ggsave(file.path(RPT_PDF, "panel_D_outlier_consensus.pdf"), pD_s01,
       width = PW, height = PH, units = "mm", device = get_pdf_device())
cat("F00 Stage01 Panel D done\n")
