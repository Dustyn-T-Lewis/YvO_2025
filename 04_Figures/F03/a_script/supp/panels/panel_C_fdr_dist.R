# F03 SUPP — Panel C: FDR (adj.P.Val) distributions per contrast
# Same layout as p-value and Pi-score histogram panels.
# Contrast titles and stats rendered inside the plot.
#
# Reads:  03_DEP/c_data/04_per_contrast_results/{Aging,...}.csv
# Writes: b_reports/supp/{png,pdf}/panels/SUPP_panel_C_fdr_dist.{png,pdf}
#         c_data/supp/panel_C_fdr_dist.csv

setwd(rprojroot::find_rstudio_root_file())
source("04_Figures/shared/style.R")

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(ggplot2)
  library(tidyr)
})

CTRS    <- c("Aging", "Training_Young", "Training_Old", "Interaction")
RPT_PNG <- "04_Figures/F03/b_reports/supp/png/panels"
RPT_PDF <- "04_Figures/F03/b_reports/supp/pdf/panels"
DAT_DIR <- "04_Figures/F03/c_data/supp"
dir.create(RPT_PNG, recursive = TRUE, showWarnings = FALSE)
dir.create(RPT_PDF, recursive = TRUE, showWarnings = FALSE)
dir.create(DAT_DIR, recursive = TRUE, showWarnings = FALSE)

stopifnot(
  "03_DEP per-contrast results missing" =
    all(file.exists(sprintf("03_DEP/c_data/04_per_contrast_results/%s.csv", CTRS)))
)

fdr_df <- lapply(CTRS, function(ctr) {
  read_csv(sprintf("03_DEP/c_data/04_per_contrast_results/%s.csv", ctr),
           show_col_types = FALSE) %>%
    transmute(contrast = ctr, gene, FDR = adj.P.Val)
}) %>% bind_rows() %>%
  filter(!is.na(FDR)) %>%
  mutate(contrast = factor(contrast, levels = CTRS))

write_csv(fdr_df, file.path(DAT_DIR, "panel_C_fdr_dist.csv"))

# Count of FDR-significant proteins per contrast
n_sig_fdr <- fdr_df %>%
  group_by(contrast) %>%
  summarize(n_sig = sum(FDR < 0.05), .groups = "drop")

n_bins <- 20

pFDR <- ggplot(fdr_df, aes(FDR)) +
  geom_histogram(breaks = seq(0, 1, length.out = n_bins + 1),
                 fill = "#9B7FBF", color = "white", linewidth = 0.3) +
  geom_vline(xintercept = 0.05, linetype = "dashed", color = "grey30",
             linewidth = 0.4) +
  # Contrast title inside plot
  geom_text(data = n_sig_fdr,
            aes(x = 0.5, y = Inf, label = CTR_SHORT[as.character(contrast)]),
            inherit.aes = FALSE, hjust = 0.5, vjust = 1.2,
            size = 2.0, fontface = "bold", color = "grey20") +
  # Stats below contrast title (larger, bold, shifted down)
  geom_text(data = n_sig_fdr,
            aes(x = 0.5, y = Inf, label = sprintf("FDR < 0.05: %d", n_sig)),
            inherit.aes = FALSE, hjust = 0.5, vjust = 3.8,
            size = 2.8, fontface = "bold", color = "grey40") +
  facet_wrap(~ contrast, ncol = 2, scales = "free_y",
             labeller = labeller(contrast = CTR_SHORT)) +
  labs(title = "FDR distribution by contrast",
       subtitle = sprintf("%s proteins per contrast | 20 bins | dashed line = FDR = 0.05",
                          format(round(nrow(fdr_df) / length(CTRS)), big.mark = ",")),
       x = "FDR (BH-corrected)",
       y = "Proteins",
       tag = "c") +
  FIG_THEME +
  theme(strip.background = element_blank(),
        strip.text = element_blank())

PW <- 89; PH <- 75
ggsave(file.path(RPT_PNG, "SUPP_panel_C_fdr_dist.png"), pFDR,
       width = PW, height = PH, units = "mm", dpi = 300)
ggsave(file.path(RPT_PDF, "SUPP_panel_C_fdr_dist.pdf"), pFDR,
       width = PW, height = PH, units = "mm", device = get_pdf_device())

message("F03 SUPP panel C (FDR distributions) saved")

# --- Export for composite ---
pFDR_title    <- "FDR distribution by contrast"
pFDR_subtitle <- sprintf("%s proteins per contrast | 20 bins | dashed line = FDR = 0.05",
                         format(round(nrow(fdr_df) / length(CTRS)), big.mark = ","))
pFDR_legend   <- NULL
pFDR          <- strip_for_composite(pFDR)
