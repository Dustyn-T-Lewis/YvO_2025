# F03 SUPP — Panel B: p-value histograms per contrast
# One facet per contrast (Aging, Training_Young, Training_Old, Interaction).
# A well-calibrated test produces a near-uniform p-value distribution on the
# null; a spike near 0 indicates many true signals. Facets show uniform + Pi
# reference fraction.
#
# Reads:  03_DEP/c_data/04_per_contrast_results/{Aging,Training_Young,Training_Old,Interaction}.csv
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

phist_df <- lapply(CTRS, function(ctr) {
  read_csv(sprintf("03_DEP/c_data/04_per_contrast_results/%s.csv", ctr),
           show_col_types = FALSE) %>%
    transmute(contrast = ctr,
              gene,
              P.Value,
              adj.P.Val,
              pi_score,
              sig_pi)
}) %>% bind_rows() %>%
  mutate(contrast = factor(contrast, levels = CTRS))

write_csv(phist_df, file.path(DAT_DIR, "panel_B_phist.csv"))

# Uniform reference density = n_proteins / n_bins
n_bins <- 20
uniform_ref <- phist_df %>%
  group_by(contrast) %>%
  summarize(h = n() / n_bins, .groups = "drop")

pB <- ggplot(phist_df, aes(P.Value)) +
  geom_histogram(breaks = seq(0, 1, length.out = n_bins + 1),
                 fill = "#5DA5DA", color = "white", linewidth = 0.3) +
  geom_hline(data = uniform_ref, aes(yintercept = h),
             linetype = "dashed", color = "grey40", linewidth = 0.4) +
  facet_wrap(~ contrast, ncol = 2, scales = "free_y",
             labeller = labeller(contrast = CTR_SHORT)) +
  labs(title = "Raw p-value distribution by contrast",
       subtitle = sprintf("%s proteins per contrast | 20 bins | dashed line = uniform null",
                          format(nrow(phist_df) / length(CTRS), big.mark = ",")),
       x = "P.Value",
       y = "Proteins",
       tag = "A") +
  FIG_THEME +
  theme(strip.background = element_blank(),
        strip.text = element_text(face = "bold", size = FIG_STRIP_SIZE))

PW <- 160; PH <- 120
ggsave(file.path(RPT_PNG, "SUPP_panel_B_phist.png"), pB,
       width = PW, height = PH, units = "mm", dpi = 300)
ggsave(file.path(RPT_PDF, "SUPP_panel_B_phist.pdf"), pB,
       width = PW, height = PH, units = "mm", device = get_pdf_device())

message("F03 SUPP panel B (p-histograms) saved")
