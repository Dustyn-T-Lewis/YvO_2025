#!/usr/bin/env Rscript
# Sourced by 02_supp_panels.R — expects style.R already loaded.
# F05 Supplementary Panel C: Circularity Diagnostic
# Diagnostic for main Panel C (fry) — addresses the shared Old_Pre contrast term.
#
# The Aging (Old_Pre - Young_Pre) and Training_Old (Old_Post - Old_Pre) contrasts
# share Old_Pre, creating structural negative correlation. This panel tests
# whether the observed correlation is more extreme than expected by chance
# from the shared contrast structure alone.
#
# Method: protein-label permutation (shuffle row labels of logFC_Aging, compute
# Pearson r with logFC_Training_Old, repeat 1000x). Directly tests whether the
# observed correlation is explained by shared contrast structure.

suppressPackageStartupMessages({
  library(tidyverse)
})

BASE    <- "04_Figures/F05"
RPT_PNG <- file.path(BASE, "b_reports", "supp", "png", "panels")
RPT_PDF <- file.path(BASE, "b_reports", "supp", "pdf", "panels")
DAT     <- file.path(BASE, "c_data", "panel_supp")
dir.create(RPT_PNG, recursive = TRUE, showWarnings = FALSE)
dir.create(RPT_PDF, recursive = TRUE, showWarnings = FALSE)
dir.create(DAT,     recursive = TRUE, showWarnings = FALSE)
pdf_device <- get_pdf_device()

# Data
dep <- read_csv("03_DEP/c_data/03_combined_results.csv",
                show_col_types = FALSE)
fc_df <- dep |>
  transmute(gene, logFC_Aging, logFC_TO = logFC_Training_Old) |>
  filter(!is.na(logFC_Aging), !is.na(logFC_TO))

obs_r <- as.numeric(cor(fc_df$logFC_Aging, fc_df$logFC_TO, use = "complete.obs"))

# Protein-label permutation
set.seed(42)
n_perm   <- 1000
perm_r   <- numeric(n_perm)
to_vec   <- fc_df$logFC_TO
for (i in seq_len(n_perm)) {
  shuffled_aging <- sample(fc_df$logFC_Aging)
  perm_r[i]      <- cor(shuffled_aging, to_vec, use = "complete.obs")
}

null_mean <- mean(perm_r)
null_sd   <- sd(perm_r)
p_perm    <- mean(abs(perm_r) >= abs(obs_r))

circ_df <- tibble(replicate = seq_len(n_perm), perm_r = perm_r)
write_csv(circ_df, file.path(DAT, "SUPP_fry_circularity.csv"))

# Plot
sub_text <- sprintf("Observed r = %.3f | Null mean = %.4f | p_perm = %.3f",
                    obs_r, null_mean, p_perm)

pS_circ <- ggplot(circ_df, aes(x = perm_r)) +
  geom_histogram(bins = 50, fill = "grey70", color = "white", linewidth = 0.3) +
  geom_vline(xintercept = obs_r, color = "#D6604D", linewidth = 0.8) +
  geom_vline(xintercept = null_mean, color = "#4393C3", linewidth = 0.5,
             linetype = "dashed") +
  coord_cartesian(clip = "off") +
  annotate("label", x = obs_r, y = Inf, vjust = 1.3, hjust = 1.1,
           label = sprintf("Observed\nr = %.3f", obs_r),
           size = 2.5, fontface = "bold", color = "#D6604D",
           fill = alpha("white", 0.9), label.padding = unit(2, "pt")) +
  annotate("label", x = null_mean, y = Inf, vjust = 1.3, hjust = -0.1,
           label = sprintf("Null mean\n= %.4f", null_mean),
           size = 2.2, fontface = "bold", color = "#4393C3",
           fill = alpha("white", 0.9), label.padding = unit(2, "pt")) +
  labs(title    = "Circularity Diagnostic: Protein-Permuted Null",
       subtitle = sub_text,
       x = "Permuted Pearson r",
       y = "Count") +
  FIG_THEME

PW <- 89; PH <- 70
ggsave(file.path(RPT_PNG, "SUPP_fry_circularity.png"), pS_circ,
       width = PW, height = PH, units = "mm", dpi = 300)
ggsave(file.path(RPT_PDF, "SUPP_fry_circularity.pdf"), pS_circ,
       width = PW, height = PH, units = "mm", device = pdf_device)

message("F05 SUPP Panel C (circularity diagnostic) saved")

# Expose for composite
pS_circ_title    <- "Circularity Diagnostic: Protein-Permuted Null"
pS_circ_subtitle <- sub_text
pS_circ          <- strip_for_composite(pS_circ)
