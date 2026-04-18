# F00 Stage 02 Panel C: MNAR audit shift distribution
setwd(rprojroot::find_rstudio_root_file())
source("04_Figures/shared/style.R")

suppressPackageStartupMessages({
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

mn_audit <- int02$mnar_audit

write.csv(mn_audit, file.path(DAT, "panel_C_mnar_audit.csv"), row.names = FALSE)

pC_s02 <- ggplot(mn_audit, aes(shift)) +
  geom_histogram(binwidth = 0.025, fill = "#E41A1C", alpha = 0.70,
                 color = "white", linewidth = 0.2) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey30",
             linewidth = 0.4) +
  geom_vline(xintercept = mean(mn_audit$shift, na.rm = TRUE),
             linetype = "solid", color = "#B2182B", linewidth = 0.4) +
  annotate("text", x = mean(mn_audit$shift, na.rm=TRUE), y = Inf,
           label = sprintf("mean = %+.3f", mean(mn_audit$shift, na.rm=TRUE)),
           vjust = 1.4, hjust = -0.1, size = 2.7, color = "#B2182B") +
  labs(title = "MNAR audit: imputation shift",
       subtitle = sprintf("%d MNAR proteins | pre-mean vs post-mean (log2)",
                          nrow(mn_audit)),
       x = "Shift (log2 intensity)", y = "MNAR proteins", tag = "C") +
  FIG_THEME

ggsave(file.path(RPT_PNG, "panel_C_mnar_audit.png"), pC_s02,
       width = PW, height = PH, units = "mm", dpi = 300)
ggsave(file.path(RPT_PDF, "panel_C_mnar_audit.pdf"), pC_s02,
       width = PW, height = PH, units = "mm", device = get_pdf_device())
cat("F00 Stage02 Panel C done\n")
