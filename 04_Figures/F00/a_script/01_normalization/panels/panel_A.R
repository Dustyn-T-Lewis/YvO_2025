# F00 Stage 01 Panel A: Filter cascade waterfall
setwd(rprojroot::find_rstudio_root_file())
source("04_Figures/shared/style.R")

suppressPackageStartupMessages({
  library(dplyr)
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

fcasc <- int01$filter_log %>%
  mutate(step   = factor(step, levels = step),
         retained = n_after,
         removed  = ifelse(is.na(n_removed), 0, n_removed))

write.csv(fcasc, file.path(DAT, "panel_A_filter_cascade.csv"), row.names = FALSE)

pA_s01 <- ggplot(fcasc, aes(step, retained)) +
  geom_col(fill = "#2166AC", width = 0.7) +
  geom_text(aes(label = format(retained, big.mark = ",")),
            vjust = -1.8, size = 2.8, fontface = "bold") +
  geom_text(aes(label = ifelse(removed > 0,
                               sprintf("-%s", format(removed, big.mark=",")),
                               "")),
            vjust = -0.4, size = 2.4, color = "#B2182B") +
  scale_y_continuous(expand = expansion(mult = c(0, 0.18)),
                     labels = comma) +
  labs(title = "Filter cascade",
       subtitle = sprintf("%s raw -> %s retained (%.1f%%)",
                          format(fcasc$retained[1], big.mark=","),
                          format(tail(fcasc$retained,1), big.mark=","),
                          tail(fcasc$pct_of_raw,1)),
       x = NULL, y = "Proteins retained", tag = "A") +
  FIG_THEME +
  theme(axis.text.x = element_text(angle = 22, hjust = 1, size = 7))

ggsave(file.path(RPT_PNG, "panel_A_filter_cascade.png"), pA_s01,
       width = PW, height = PH, units = "mm", dpi = 300)
ggsave(file.path(RPT_PDF, "panel_A_filter_cascade.pdf"), pA_s01,
       width = PW, height = PH, units = "mm", device = get_pdf_device())
cat("F00 Stage01 Panel A done\n")
