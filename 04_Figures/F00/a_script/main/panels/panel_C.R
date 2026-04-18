# F00 Panel C: Pre-normalization PCA (Mahalanobis outlier flags)
setwd(rprojroot::find_rstudio_root_file())
source("04_Figures/F00/a_script/style.R")

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
})

PW <- 100; PH <- 100
RPT <- "04_Figures/F00/b_reports/panels"
DAT <- "04_Figures/F00/c_data"
dir.create(RPT, recursive = TRUE, showWarnings = FALSE)
dir.create(DAT, recursive = TRUE, showWarnings = FALSE)

int_norm <- readRDS("01_normalization/c_data/00_report_intermediates.rds")

pca_pre_df <- int_norm$pca_pre$scores %>%
  left_join(int_norm$outlier_diag %>% select(Col_ID, pca_flag, prefix), by = "Col_ID") %>%
  mutate(age = ifelse(grepl("^O", prefix), "Old", "Young"),
         age = factor(age, levels = c("Young", "Old")))

write.csv(pca_pre_df, file.path(DAT, "panel_C_pca_pre.csv"), row.names = FALSE)

pC <- ggplot(pca_pre_df, aes(PC1, PC2, color = age, shape = Timepoint)) +
  geom_point(size = 2.2, alpha = 0.85) +
  scale_color_manual(values = AGE_COLORS, name = NULL) +
  scale_shape_manual(values = SHAPE_TP, name = NULL) +
  labs(title = "Pre-normalization PCA",
       subtitle = sprintf("PC1 %.1f%% | PC2 %.1f%%",
                          int_norm$pca_pre$var_exp[1],
                          int_norm$pca_pre$var_exp[2]),
       x = sprintf("PC1 (%.1f%%)", int_norm$pca_pre$var_exp[1]),
       y = sprintf("PC2 (%.1f%%)", int_norm$pca_pre$var_exp[2]),
       tag = "C") +
  coord_fixed() +
  FIG_THEME + theme(legend.position = "top",
                    legend.key.size = unit(3, "mm"))

ggsave(file.path(RPT, "MAIN_panel_C_pca_pre.png"), pC,
       width = PW, height = PH, units = "mm", dpi = 300)
cat("F00 Panel C done\n")
