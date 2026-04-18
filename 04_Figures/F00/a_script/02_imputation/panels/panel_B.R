# F00 Stage 02 Panel B: Benchmark top-5 methods
setwd(rprojroot::find_rstudio_root_file())
source("04_Figures/shared/style.R")

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(ggplot2)
})

PW <- 100; PH <- 80
RPT_PNG <- "04_Figures/F00/b_reports/02_imputation/png/panels"
RPT_PDF <- "04_Figures/F00/b_reports/02_imputation/pdf/panels"
DAT     <- "04_Figures/F00/c_data/02_imputation"
dir.create(RPT_PNG, recursive = TRUE, showWarnings = FALSE)
dir.create(RPT_PDF, recursive = TRUE, showWarnings = FALSE)
dir.create(DAT, recursive = TRUE, showWarnings = FALSE)

if (!exists("int02"))
  int02 <- readRDS("02_Imputation/c_data/00_report_intermediates.rds")

bench <- read_csv("02_Imputation/c_data/benchmark/04_composite_ranking.csv",
                   show_col_types = FALSE)

topN <- bench %>%
  arrange(rank) %>%
  slice_head(n = 5) %>%
  mutate(method = factor(method, levels = rev(method)),
         is_winner = method == levels(method)[length(levels(method))])

write.csv(topN, file.path(DAT, "panel_B_benchmark_top5.csv"), row.names = FALSE)

pB_s02 <- ggplot(topN, aes(composite, method,
                      fill = ifelse(is_winner, "Winner", "Top 5"))) +
  geom_col(width = 0.7) +
  geom_text(aes(label = sprintf("%.3f", composite)),
            hjust = -0.1, size = 2.7, fontface = "bold") +
  scale_fill_manual(values = c(Winner = "#E41A1C", "Top 5" = "#377EB8"),
                    name = NULL) +
  scale_x_continuous(expand = expansion(mult = c(0, 0.15)),
                     limits = c(0, max(topN$composite) * 1.15)) +
  labs(title = "23-method benchmark (top 5)",
       subtitle = sprintf("Composite score | rank 1: %s (%.3f)",
                          as.character(topN$method[which.max(topN$composite)]),
                          max(topN$composite)),
       x = "Composite score", y = NULL, tag = "B") +
  FIG_THEME +
  theme(legend.position = "top",
        legend.key.size = unit(3, "mm"),
        legend.text = element_text(size = 7),
        axis.text.y = element_text(size = 7))

ggsave(file.path(RPT_PNG, "panel_B_benchmark.png"), pB_s02,
       width = PW, height = PH, units = "mm", dpi = 300)
ggsave(file.path(RPT_PDF, "panel_B_benchmark.pdf"), pB_s02,
       width = PW, height = PH, units = "mm", device = get_pdf_device())
cat("F00 Stage02 Panel B done\n")
