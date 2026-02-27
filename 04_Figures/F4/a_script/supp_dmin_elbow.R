################################################################################
#   Figure 4 — Supplementary: Dmin Elbow Plot
#
#   Requires from setup: dmin_vals, optimal_k, THEME_PUB, RPT_DIR
#   Outputs: p_dmin (ggplot object), saved as supp_dmin_elbow.pdf
################################################################################

if (!exists("core_proteins")) source("04_Figures/F4/a_script/YvO_F4_setup.R")

cat("=== Building supplementary Dmin elbow plot ===\n")

p_dmin <- ggplot(tibble(k = 2:6, dmin = dmin_vals), aes(k, dmin)) +
  geom_line(linewidth = 0.8) + geom_point(size = 2) +
  geom_vline(xintercept = optimal_k, linetype = "dashed", color = "red") +
  labs(x = "k", y = "Dmin", title = "Cluster Selection: Dmin Elbow") + THEME_PUB

ggsave(file.path(RPT_DIR, "supp_dmin_elbow.pdf"), p_dmin,
       width = 120, height = 90, units = "mm")

cat("  Supplementary Dmin elbow plot saved\n")
