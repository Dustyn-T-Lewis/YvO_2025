# F00 Panel A: Protein filter cascade
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

filter_df <- int_norm$filter_bar_data
write.csv(filter_df, file.path(DAT, "panel_A_filter_cascade.csv"), row.names = FALSE)

pA <- ggplot(filter_df, aes(step, n, fill = status)) +
  geom_col(width = 0.7) +
  geom_text(aes(label = format(n, big.mark = ",")),
            position = position_stack(vjust = 0.5),
            size = 3.0, color = "white", fontface = "bold") +
  scale_fill_manual(values = c(Retained = "#2166AC", Removed = "#B2182B"), name = NULL) +
  labs(title = "Protein filter cascade",
       subtitle = sprintf("%s raw -> %s retained",
                          format(int_norm$n_raw, big.mark = ","),
                          format(int_norm$dal_nrow, big.mark = ",")),
       x = NULL, y = "Proteins", tag = "A") +
  FIG_THEME +
  theme(axis.text.x = element_text(angle = 25, hjust = 1),
        legend.position = "top",
        legend.key.size = unit(3, "mm"))

ggsave(file.path(RPT, "MAIN_panel_A_filter_cascade.png"), pA,
       width = PW, height = PH, units = "mm", dpi = 300)
cat("F00 Panel A done\n")
