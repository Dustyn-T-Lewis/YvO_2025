# F00 Panel D: Post-normalization PCA with 80% ellipses
setwd(rprojroot::find_rstudio_root_file())
source("04_Figures/shared/style.R")

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
})

PW <- 100; PH <- 100
RPT_PNG <- "04_Figures/F00/b_reports/main/png/panels"
RPT_PDF <- "04_Figures/F00/b_reports/main/pdf/panels"
DAT     <- "04_Figures/F00/c_data"
dir.create(RPT_PNG, recursive = TRUE, showWarnings = FALSE)
dir.create(RPT_PDF, recursive = TRUE, showWarnings = FALSE)
dir.create(DAT, recursive = TRUE, showWarnings = FALSE)

if (!exists("int_norm"))
  int_norm <- readRDS("01_normalization/c_data/00_report_intermediates.rds")

pca_post_df <- int_norm$pca_post$scores

write.csv(pca_post_df, file.path(DAT, "panel_D_pca_post.csv"), row.names = FALSE)

pD <- ggplot(pca_post_df, aes(PC1, PC2, color = Group_Time,
                              fill = Group_Time, shape = Timepoint)) +
  stat_ellipse(aes(group = Group_Time), geom = "polygon",
               alpha = 0.12, level = 0.80, linewidth = 0.3,
               show.legend = FALSE) +
  geom_point(size = 2.2, alpha = 0.85) +
  scale_color_manual(values = PCA_COLORS,
                     labels = c("Young Pre", "Young Post", "Old Pre", "Old Post"),
                     name = NULL) +
  scale_fill_manual(values = PCA_COLORS, guide = "none") +
  scale_shape_manual(values = SHAPE_TP, guide = "none") +
  labs(title = "Post-normalization PCA",
       subtitle = sprintf("PC1 %.1f%% | PC2 %.1f%% | cyclic loess",
                          int_norm$pca_post$var_exp[1],
                          int_norm$pca_post$var_exp[2]),
       x = sprintf("PC1 (%.1f%%)", int_norm$pca_post$var_exp[1]),
       y = sprintf("PC2 (%.1f%%)", int_norm$pca_post$var_exp[2]),
       tag = "D") +
  coord_fixed() +
  FIG_THEME + theme(legend.position = "top",
                    legend.key.size = unit(3, "mm"),
                    legend.text = element_text(size = 7))

ggsave(file.path(RPT_PNG, "MAIN_panel_D_pca_post.png"), pD,
       width = PW, height = PH, units = "mm", dpi = 300)
ggsave(file.path(RPT_PDF, "MAIN_panel_D_pca_post.pdf"), pD,
       width = PW, height = PH, units = "mm", device = get_pdf_device())
cat("F00 Panel D done\n")
