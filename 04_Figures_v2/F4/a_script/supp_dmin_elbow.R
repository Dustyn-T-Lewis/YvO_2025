################################################################################
#   Figure 4 — Supplementary: Dmin Elbow Plot
#   Self-contained: loads pre-computed data from c_data/
#   Generates: supp_dmin_elbow.pdf/png
################################################################################

setwd(rprojroot::find_rstudio_root_file())
source("04_Figures_v2/shared/style.R")

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(tibble)
})

# ── Paths & constants ────────────────────────────────────────────────────────
RPT <- "04_Figures_v2/F4/b_reports"
DAT <- "04_Figures_v2/F4/c_data"
optimal_k <- 4

# ── Load data ────────────────────────────────────────────────────────────────
dmin_vals <- read_csv(file.path(DAT, "dmin_values.csv"), show_col_types = FALSE) %>%
  pull(dmin)

boot_ari <- read_csv(file.path(DAT, "07_cluster_stability.csv"), show_col_types = FALSE) %>%
  pull(ari)

pdf_device <- get_pdf_device()

message("Supplementary: Dmin elbow plot...")

PS_W <- 120
PS_H <- 90

boot_ari_ci <- quantile(boot_ari, probs = c(0.025, 0.975))
ari_label <- sprintf("k = %d: ARI = %.3f\n95%% CI [%.3f, %.3f]",
                     optimal_k, mean(boot_ari),
                     boot_ari_ci[1], boot_ari_ci[2])

txt_annot <- scale_text(BASE_STAT, PS_W)

p_dmin <- ggplot(tibble(k = 2:6, dmin = dmin_vals), aes(k, dmin)) +
  geom_line(linewidth = 0.8) + geom_point(size = 2) +
  geom_vline(xintercept = optimal_k, linetype = "dashed", color = "red") +
  annotate("text", x = 5.5, y = max(dmin_vals) * 0.9,
           label = ari_label, size = txt_annot, hjust = 1, color = "red") +
  labs(x = "k", y = "Dmin", title = "Cluster Selection: Dmin Elbow") +
  FIG_THEME +
  theme(legend.position = "none")

ggsave(file.path(RPT, "supp_dmin_elbow.pdf"), p_dmin,
       width = PS_W, height = PS_H, units = "mm", device = pdf_device)
ggsave(file.path(RPT, "supp_dmin_elbow.png"), p_dmin,
       width = PS_W, height = PS_H, units = "mm", dpi = 300)

message("  Supplementary Dmin elbow saved")
