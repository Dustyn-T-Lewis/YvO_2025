# F00 Stage 03 Panel A: DEP count matrix heatmap (5 contrasts x 4 thresholds)
setwd(rprojroot::find_rstudio_root_file())
source("04_Figures/shared/style.R")

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(readr)
  library(ggplot2)
})

PW <- 100; PH <- 80
RPT_PNG <- "04_Figures/F00/b_reports/03_DEP/png/panels"
RPT_PDF <- "04_Figures/F00/b_reports/03_DEP/pdf/panels"
DAT     <- "04_Figures/F00/c_data/03_DEP"
dir.create(RPT_PNG, recursive = TRUE, showWarnings = FALSE)
dir.create(RPT_PDF, recursive = TRUE, showWarnings = FALSE)
dir.create(DAT, recursive = TRUE, showWarnings = FALSE)

CONTRAST_ORDER <- c("Aging", "Training_Young", "Training_Old",
                    "Interaction", "Reversal")

if (!exists("da_summ"))
  da_summ <- read_csv("03_DEP/c_data/02_DA_summary.csv", show_col_types = FALSE)

counts <- da_summ %>%
  filter(type %in% c("up", "down")) %>%
  group_by(contrast) %>%
  summarise(`p<0.05`    = sum(sig.PVal),
            `FDR<0.10`  = sum(sig.FDR.10),
            `FDR<0.05`  = sum(sig.FDR.05),
            `Pi<0.05`   = sum(sig.Pi),
            .groups = "drop") %>%
  pivot_longer(-contrast, names_to = "threshold", values_to = "n") %>%
  mutate(contrast  = factor(contrast,
                            levels = intersect(CONTRAST_ORDER, unique(contrast))),
         threshold = factor(threshold,
                            levels = c("p<0.05", "FDR<0.10", "FDR<0.05", "Pi<0.05")))

write.csv(counts, file.path(DAT, "panel_A_dep_counts.csv"), row.names = FALSE)

pA_s03 <- ggplot(counts, aes(threshold, contrast, fill = n)) +
  geom_tile(color = "white", linewidth = 0.8) +
  geom_text(aes(label = format(n, big.mark = ","),
                color = n > max(n) * 0.5),
            size = 2.8, fontface = "bold") +
  scale_color_manual(values = c("TRUE" = "white", "FALSE" = "grey20"),
                     guide = "none") +
  scale_fill_gradient(low = "#DEEBF7", high = "#08519C", name = "DEPs") +
  labs(title = "DEP counts: contrast x threshold",
       subtitle = "up + down (nonsig excluded)",
       x = NULL, y = NULL, tag = "A") +
  FIG_THEME +
  theme(axis.text.x = element_text(angle = 15, hjust = 1, size = 7),
        axis.text.y = element_text(size = 7),
        legend.position = "right",
        legend.key.size = unit(3, "mm"),
        legend.text = element_text(size = 7))

ggsave(file.path(RPT_PNG, "panel_A_dep_counts.png"), pA_s03,
       width = PW, height = PH, units = "mm", dpi = 300)
ggsave(file.path(RPT_PDF, "panel_A_dep_counts.pdf"), pA_s03,
       width = PW, height = PH, units = "mm", device = get_pdf_device())
cat("F00 Stage03 Panel A done\n")
