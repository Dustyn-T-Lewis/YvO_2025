# F00 Panel B: Per-sample missingness
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

miss_bar <- int_norm$miss_bar_data %>%
  filter(status == "Missing") %>%
  mutate(Group_Time = factor(Group_Time,
                             levels = c("Young_Pre", "Young_Post",
                                        "Old_Pre", "Old_Post")))

write.csv(miss_bar, file.path(DAT, "panel_B_sample_missingness.csv"), row.names = FALSE)

pB <- ggplot(miss_bar, aes(reorder(Col_ID, -n), n, fill = Group_Time)) +
  geom_col(aes(alpha = is_outlier), width = 0.8) +
  scale_fill_manual(values = GROUP_FILL, name = NULL) +
  scale_alpha_manual(values = c("FALSE" = 1, "TRUE" = 0.3), guide = "none") +
  labs(title = "Per-sample missingness",
       subtitle = sprintf("%d samples | consensus outliers faded",
                          length(unique(miss_bar$Col_ID))),
       x = NULL, y = "Missing proteins", tag = "B") +
  FIG_THEME +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        legend.position = "top",
        legend.key.size = unit(3, "mm"))

ggsave(file.path(RPT, "MAIN_panel_B_sample_missingness.png"), pB,
       width = PW, height = PH, units = "mm", dpi = 300)
cat("F00 Panel B done\n")
