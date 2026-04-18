# F00 Stage 01 Panel B: Per-protein missingness distribution
setwd(rprojroot::find_rstudio_root_file())
source("04_Figures/shared/style.R")

suppressPackageStartupMessages({
  library(tibble)
  library(ggplot2)
  library(scales)
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

prot_miss <- rowSums(is.na(int01$data_pre_outlier)) / ncol(int01$data_pre_outlier) * 100
miss_df <- tibble(pct_miss = prot_miss)
thresh_pct <- quantile(prot_miss, 0.99, na.rm = TRUE)

write.csv(miss_df, file.path(DAT, "panel_B_missingness.csv"), row.names = FALSE)

pB_s01 <- ggplot(miss_df, aes(pct_miss)) +
  geom_histogram(binwidth = 2.5, fill = "#4393C3", color = "white", alpha = 0.85) +
  geom_vline(xintercept = thresh_pct, linetype = "dashed", color = "#B2182B", linewidth = 0.4) +
  annotate("text", x = thresh_pct, y = Inf,
           label = sprintf("99th pct = %.1f%%", thresh_pct),
           vjust = 1.6, hjust = 1.05, color = "#B2182B", size = 2.7) +
  scale_y_continuous(labels = comma) +
  labs(title = "Per-protein missingness (post-filter)",
       subtitle = sprintf("%s proteins | median = %.1f%%",
                          format(nrow(miss_df), big.mark=","),
                          median(prot_miss, na.rm=TRUE)),
       x = "% missing across samples", y = "Proteins", tag = "B") +
  FIG_THEME

ggsave(file.path(RPT_PNG, "panel_B_missingness.png"), pB_s01,
       width = PW, height = PH, units = "mm", dpi = 300)
ggsave(file.path(RPT_PDF, "panel_B_missingness.pdf"), pB_s01,
       width = PW, height = PH, units = "mm", device = get_pdf_device())
cat("F00 Stage01 Panel B done\n")
