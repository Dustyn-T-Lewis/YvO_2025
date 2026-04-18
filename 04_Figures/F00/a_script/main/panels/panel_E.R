# F00 Panel E: MAR/MNAR classification summary
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

if (!exists("int_imp"))
  int_imp <- readRDS("02_Imputation/c_data/00_report_intermediates.rds")

class_counts <- int_imp$miss_class %>%
  count(classification) %>%
  mutate(classification = factor(classification,
                                 levels = c("Complete", "MAR", "MNAR")),
         pct = 100 * n / sum(n))

write.csv(class_counts, file.path(DAT, "panel_E_miss_classification.csv"), row.names = FALSE)

pE <- ggplot(class_counts, aes(classification, n, fill = classification)) +
  geom_col(width = 0.6) +
  geom_text(aes(label = sprintf("%d\n(%.1f%%)", n, pct)),
            vjust = -0.2, size = 3.2, fontface = "bold") +
  scale_fill_manual(values = PAL_CLASS, guide = "none") +
  scale_y_continuous(expand = expansion(mult = c(0, 0.22))) +
  labs(title = "Missingness classification",
       subtitle = sprintf("%s proteins | %s classification",
                          format(sum(class_counts$n), big.mark = ","),
                          int_imp$classification_method %||% "k-means"),
       x = NULL, y = "Proteins", tag = "E") +
  FIG_THEME

ggsave(file.path(RPT_PNG, "MAIN_panel_E_miss_classification.png"), pE,
       width = PW, height = PH, units = "mm", dpi = 300)
ggsave(file.path(RPT_PDF, "MAIN_panel_E_miss_classification.pdf"), pE,
       width = PW, height = PH, units = "mm", device = get_pdf_device())
cat("F00 Panel E done\n")
