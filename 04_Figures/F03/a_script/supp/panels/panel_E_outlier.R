# F03 SUPP — Panel E: DEP count retention after outlier removal
# Paired bars per contrast: counts under FULL vs REDUCED (leave-outliers-out)
# cohort for both FDR < 0.05 and Pi-score < 0.05 criteria.
#
# Reads:  03_DEP/c_data/11_outlier_sensitivity.csv
# Writes: b_reports/supp/SUPP_panel_E_outlier.png
#         c_data/supp/panel_E_outlier.csv

setwd(rprojroot::find_rstudio_root_file())
source("04_Figures/shared/style.R")

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(readr)
  library(ggplot2)
})

RPT_DIR <- "04_Figures/F03/b_reports/supp/panels"
DAT_DIR <- "04_Figures/F03/c_data/supp"
dir.create(RPT_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(DAT_DIR, recursive = TRUE, showWarnings = FALSE)

out_sens <- read_csv("03_DEP/c_data/11_outlier_sensitivity.csv",
                     show_col_types = FALSE) %>%
  mutate(Contrast = factor(Contrast,
                           levels = c("Aging", "Training_Young",
                                      "Training_Old", "Interaction",
                                      "Reversal")))

write_csv(out_sens, file.path(DAT_DIR, "panel_E_outlier.csv"))

long_df <- out_sens %>%
  pivot_longer(c(FDR_full, FDR_reduced, Pi_full, Pi_reduced),
               names_to  = "metric_cohort", values_to = "n") %>%
  mutate(metric  = sub("_.*", "", metric_cohort),
         cohort  = sub(".*_", "", metric_cohort),
         cohort  = factor(cohort, levels = c("full", "reduced")),
         metric  = factor(metric, levels = c("FDR", "Pi"),
                          labels = c("FDR < 0.05", "Pi < 0.05")))

rho_df <- out_sens %>%
  transmute(Contrast,
            label = sprintf("rho = %.4f", Spearman_rho))

pE <- ggplot(long_df, aes(cohort, n, fill = cohort)) +
  geom_col(width = 0.65) +
  geom_text(aes(label = n), vjust = -0.3, size = 2.8, fontface = "bold") +
  facet_grid(metric ~ Contrast, scales = "free_y", switch = "y",
             labeller = labeller(Contrast = c(Aging = "Aging",
                                              Training_Young = "Tr.(Y)",
                                              Training_Old = "Tr.(O)",
                                              Interaction = "Tr.(O)-Tr.(Y)",
                                              Reversal = "Reversal"))) +
  scale_fill_manual(values = c(full = "#2166AC", reduced = "#B2182B"),
                    labels = c(full = "Full cohort", reduced = "Outlier-removed"),
                    name = NULL) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.22))) +
  labs(title = "DEP count retention after outlier removal",
       subtitle = "Bars show DEP counts under FULL vs outlier-REDUCED cohort; per-contrast Spearman rho printed in logs",
       x = NULL, y = "DEPs",
       tag = "D") +
  FIG_THEME +
  theme(strip.background = element_blank(),
        strip.text = element_text(face = "bold", size = FIG_STRIP_SIZE - 1),
        legend.position = "top",
        legend.key.size = unit(3, "mm"),
        axis.text.x = element_text(size = 7))

message(sprintf("Outlier sensitivity rho summary:\n%s",
                paste(sprintf("  %s: %s", rho_df$Contrast, rho_df$label),
                      collapse = "\n")))

PW <- 180; PH <- 120
ggsave(file.path(RPT_DIR, "SUPP_panel_E_outlier.png"), pE,
       width = PW, height = PH, units = "mm", dpi = 300)

message("F03 SUPP panel E (outlier sensitivity) saved")
