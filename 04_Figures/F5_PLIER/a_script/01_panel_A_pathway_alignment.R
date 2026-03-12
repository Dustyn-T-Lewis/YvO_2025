# Panel A — Pathway alignment heatmap.
# Which LVs align cleanly to known biology?

setwd(rprojroot::find_rstudio_root_file())
source("04_Figures/shared/style.R")

library(readr)
library(dplyr)
library(tidyr)
library(ggplot2)

DAT <- "04_Figures/F5_PLIER/c_data"
RPT <- "04_Figures/F5_PLIER/b_results"
dir.create(RPT, recursive = TRUE, showWarnings = FALSE)

# --- Load PLIER outputs ---
plier_res <- readRDS(file.path(DAT, "01_plier_fit.rds"))
align_df <- read_csv(file.path(DAT, "05_lv_pathway_alignment_summary.csv"),
                     show_col_types = FALSE)

# --- Identify pathway-aligned LVs ---
# Use FDR < 0.05 and AUC > 0.60 as quality thresholds
sig_align <- align_df |>

  filter(FDR < 0.05, AUC > 0.60) |>
  arrange(LV, desc(AUC))

# Top pathway per LV
top_per_lv <- sig_align |>
  group_by(LV) |>
  slice_max(AUC, n = 1, with_ties = FALSE) |>
  ungroup() |>
  mutate(pathway_clean = clean_pathway_name(pathway))

# For heatmap: show top 3 pathways per aligned LV
top_heatmap <- sig_align |>
  group_by(LV) |>
  slice_max(AUC, n = 3, with_ties = FALSE) |>
  ungroup() |>
  mutate(pathway_clean = clean_pathway_name(pathway),
         lv_label = paste0("LV", LV))

if (nrow(top_heatmap) == 0) {
  cat("No significant pathway-LV associations found. Skipping Panel A.\n")
  quit(save = "no")
}

# --- Heatmap of AUC values ---
pA <- ggplot(top_heatmap, aes(x = lv_label, y = reorder(pathway_clean, AUC),
                               fill = AUC)) +
  geom_tile(color = "white", linewidth = 0.3) +
  scale_fill_gradient(low = "grey95", high = "#D6604D", limits = c(0.5, 1),
                      oob = scales::squish, name = "AUC") +
  labs(title = "A  Pathway-Aligned Latent Variables",
       subtitle = "PLIER U-matrix, FDR < 0.05 & AUC > 0.60",
       x = NULL, y = NULL) +
  FIG_THEME +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 7),
        axis.text.y = element_text(size = 6),
        legend.position = "right")

ggsave(file.path(RPT, "11_panel_A_pathway_alignment.pdf"), pA,
       width = 200, height = 160, units = "mm")
ggsave(file.path(RPT, "11_panel_A_pathway_alignment.png"), pA,
       width = 200, height = 160, units = "mm", dpi = 300)

write_csv(top_per_lv, file.path(DAT, "11_panel_A_top_pathway_alignment.csv"))

cat(sprintf("Panel A: %d pathway-aligned LVs, %d significant associations\n",
            n_distinct(top_per_lv$LV), nrow(sig_align)))
