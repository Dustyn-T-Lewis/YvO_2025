# F03 SUPP — Panel A: p-value histograms per contrast
# One facet per contrast (Aging, Training_Young, Training_Old, Interaction).
# Contrast titles and stats (p < 0.05 count) rendered inside the plot.
#
# Reads:  03_DEP/c_data/04_per_contrast_results/{Aging,...}.csv
# Writes: b_reports/supp/{png,pdf}/panels/SUPP_panel_B_phist.{png,pdf}
#         c_data/supp/panel_B_phist.csv

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

phist_df <- lapply(CTRS, function(ctr) {
  read_csv(sprintf("03_DEP/c_data/04_per_contrast_results/%s.csv", ctr),
           show_col_types = FALSE) %>%
    transmute(contrast = ctr, gene, p_value = P.Value, FDR = adj.P.Val, pi_score, sig_pi)
}) %>% bind_rows() %>%
  mutate(contrast = factor(contrast, levels = CTRS))

write_csv(phist_df, file.path(DAT_DIR, "panel_B_phist.csv"))

# Uniform reference density = n_proteins / n_bins
n_bins <- 20
uniform_ref <- phist_df %>%
  group_by(contrast) %>%
  summarize(h = n() / n_bins, .groups = "drop")

# Count of p < 0.05 proteins per contrast (stats inside plot)
n_sig_p <- phist_df %>%
  group_by(contrast) %>%
  summarize(n_sig = sum(p_value < 0.05, na.rm = TRUE), .groups = "drop")

pB <- ggplot(phist_df, aes(p_value)) +
  geom_histogram(breaks = seq(0, 1, length.out = n_bins + 1),
                 fill = "#5DA5DA", color = "white", linewidth = 0.3) +
  geom_hline(data = uniform_ref, aes(yintercept = h),
             linetype = "dashed", color = "grey40", linewidth = 0.4) +
  # Contrast title inside plot
  geom_text(data = n_sig_p,
            aes(x = 0.5, y = Inf, label = CTR_SHORT[as.character(contrast)]),
            inherit.aes = FALSE, hjust = 0.5, vjust = 1.2,
            size = 2.0, fontface = "bold", color = "grey20") +
  # Stats below contrast title (larger, bold, shifted down)
  geom_text(data = n_sig_p,
            aes(x = 0.5, y = Inf, label = sprintf("p < 0.05: %d", n_sig)),
            inherit.aes = FALSE, hjust = 0.5, vjust = 3.8,
            size = 2.8, fontface = "bold", color = "grey40") +
  facet_wrap(~ contrast, ncol = 2, scales = "free_y",
             labeller = labeller(contrast = CTR_SHORT)) +
  labs(title = "Raw p-value distribution by contrast",
       subtitle = sprintf("%s proteins per contrast | 20 bins | dashed line = uniform null",
                          format(nrow(phist_df) / length(CTRS), big.mark = ",")),
       x = "p_value",
       y = "Proteins",
       tag = "a") +
  FIG_THEME +
  theme(strip.background = element_blank(),
        strip.text = element_blank())

PW <- 89; PH <- 75
ggsave(file.path(RPT_PNG, "SUPP_panel_B_phist.png"), pB,
       width = PW, height = PH, units = "mm", dpi = 300)
ggsave(file.path(RPT_PDF, "SUPP_panel_B_phist.pdf"), pB,
       width = PW, height = PH, units = "mm", device = get_pdf_device())

message("F03 SUPP panel A (p-histograms) saved")

# --- Export for composite ---
pB_title    <- "Raw p-value distribution by contrast"
pB_subtitle <- sprintf("%s proteins per contrast | 20 bins | dashed line = uniform null",
                       format(nrow(phist_df) / length(CTRS), big.mark = ","))
pB_legend   <- NULL
pB          <- strip_for_composite(pB)
