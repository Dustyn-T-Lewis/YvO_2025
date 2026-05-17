#!/usr/bin/env Rscript
# Sourced by 02_supp_panels.R — expects style.R already loaded.
# F05 Supplementary Panel B: Pearson r Bootstrap
# Diagnostic for main Panels A & D — bootstrap confidence interval for Pearson r
# between logFC_Aging and logFC_Training_Old.

suppressPackageStartupMessages({
  library(tidyverse)
  library(boot)
})

BASE    <- "04_Figures/F05"
RPT_PNG <- file.path(BASE, "b_reports", "supp", "png", "panels")
RPT_PDF <- file.path(BASE, "b_reports", "supp", "pdf", "panels")
DAT     <- file.path(BASE, "c_data", "panel_supp")
dir.create(RPT_PNG, recursive = TRUE, showWarnings = FALSE)
dir.create(RPT_PDF, recursive = TRUE, showWarnings = FALSE)
dir.create(DAT,     recursive = TRUE, showWarnings = FALSE)
pdf_device <- get_pdf_device()

# Data
dep <- read_csv("03_DEP/c_data/03_combined_results.csv",
                show_col_types = FALSE)
fc_df <- dep |>
  transmute(gene, logFC_Aging, logFC_TO = logFC_Training_Old) |>
  filter(!is.na(logFC_Aging), !is.na(logFC_TO))

# Bootstrap
set.seed(42)
boot_r <- function(data, indices) {
  d <- data[indices, ]
  cor(d$logFC_Aging, d$logFC_TO, use = "complete.obs")
}

b <- boot(fc_df, statistic = boot_r, R = 1000)
obs_r  <- b$t0
ci     <- boot.ci(b, type = "perc", conf = 0.95)
ci_lo  <- ci$percent[4]
ci_hi  <- ci$percent[5]

boot_df <- tibble(replicate = seq_len(1000), r = as.numeric(b$t))
write_csv(boot_df, file.path(DAT, "SUPP_r_bootstrap.csv"))

# Plot
pS_r_boot <- ggplot(boot_df, aes(x = r)) +
  geom_histogram(bins = 40, fill = "#5DA5DA", color = "white", linewidth = 0.3) +
  geom_vline(xintercept = obs_r, color = "#D6604D", linewidth = 0.7, linetype = "solid") +
  geom_vline(xintercept = c(ci_lo, ci_hi), color = "#D6604D",
             linewidth = 0.5, linetype = "dashed") +
  annotate("label", x = obs_r, y = Inf, vjust = 1.5, hjust = 1.1,
           label = sprintf("r = %.3f\n95%% CI [%.3f, %.3f]", obs_r, ci_lo, ci_hi),
           size = BASE_STAT, fontface = "bold", fill = alpha("white", 0.9),
           label.padding = unit(1.5, "pt")) +
  labs(title    = "Reversal Pearson r Bootstrap",
       subtitle = sprintf("boot::boot() R = 1000 | n = %d proteins | percentile CI", nrow(fc_df)),
       x = "Pearson r (logFC Aging vs logFC Training Old)",
       y = "Count") +
  FIG_THEME

PW <- 89; PH <- 70
ggsave(file.path(RPT_PNG, "SUPP_r_bootstrap.png"), pS_r_boot,
       width = PW, height = PH, units = "mm", dpi = 300)
ggsave(file.path(RPT_PDF, "SUPP_r_bootstrap.pdf"), pS_r_boot,
       width = PW, height = PH, units = "mm", device = pdf_device)

message("F05 SUPP Panel B (Pearson r bootstrap) saved")

# Expose for composite
pS_rboot_title    <- "Reversal Pearson r Bootstrap"
pS_rboot_subtitle <- sprintf("boot::boot() R = 1000 | n = %d | r = %.3f [%.3f, %.3f]",
                              nrow(fc_df), obs_r, ci_lo, ci_hi)
pS_r_boot         <- strip_for_composite(pS_r_boot)
