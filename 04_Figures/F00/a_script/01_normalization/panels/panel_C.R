# F00 Stage 01 Panel C: Eta-squared distribution
setwd(rprojroot::find_rstudio_root_file())
source("04_Figures/shared/style.R")

suppressPackageStartupMessages({
  library(tibble)
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

eta_df <- tibble(eta2 = as.numeric(int01$eta2_vals))
eta_med <- median(eta_df$eta2, na.rm = TRUE)
eta_90  <- quantile(eta_df$eta2, 0.90, na.rm = TRUE)

write.csv(eta_df, file.path(DAT, "panel_C_eta_squared.csv"), row.names = FALSE)

pC_s01 <- ggplot(eta_df, aes(eta2)) +
  geom_density(fill = "#4CAF50", alpha = 0.55, linewidth = 0.4) +
  geom_vline(xintercept = eta_med, linetype = "solid", color = "grey30",
             linewidth = 0.4) +
  geom_vline(xintercept = eta_90, linetype = "dashed", color = "#E05A4E",
             linewidth = 0.4) +
  annotate("text", x = eta_med, y = Inf, label = sprintf("median = %.2f", eta_med),
           vjust = 1.5, hjust = -0.1, size = 2.7, color = "grey20") +
  annotate("text", x = eta_90, y = Inf, label = sprintf("90th = %.2f", eta_90),
           vjust = 3.0, hjust = -0.1, size = 2.7, color = "#E05A4E") +
  labs(title = expression("Biological signal ("*eta^2*")"),
       subtitle = sprintf("%s proteins | primary factor explained variance",
                          format(nrow(eta_df), big.mark=",")),
       x = expression(eta^2), y = "Density", tag = "C") +
  FIG_THEME

ggsave(file.path(RPT_PNG, "panel_C_eta_squared.png"), pC_s01,
       width = PW, height = PH, units = "mm", dpi = 300)
ggsave(file.path(RPT_PDF, "panel_C_eta_squared.pdf"), pC_s01,
       width = PW, height = PH, units = "mm", device = get_pdf_device())
cat("F00 Stage01 Panel C done\n")
