################################################################################
#   Figure 4 — Supplementary: Dmin Elbow Plot
#
#   Requires from setup: dmin_vals, optimal_k, THEME_PUB, RPT_DIR
#   Outputs: p_dmin (ggplot object), saved as supp_dmin_elbow.pdf
################################################################################
#
# STAT AUDIT (2026-02-27)
# ---------------------------------------------------------------------------
# 1. Dmin elbow:
#    - Dmin = minimum distance between cluster centroids for each k.
#      The elbow (point of diminishing returns in centroid separation)
#      is used to select k. This is a standard heuristic for FCM.       PASS
#
# 2. No formal test:
#    - Elbow selection is visual/heuristic. No p-value or CI applies.
#      The choice of k = 4 is validated by bootstrap ARI (mean ~0.96,
#      95% CI reported in setup), confirming high stability.            PASS
# ---------------------------------------------------------------------------

if (!exists("core_proteins")) source("04_Figures/F4/a_script/YvO_F4_setup.R")

cat("=== Building supplementary Dmin elbow plot ===\n")

# Annotate with bootstrap ARI 95% CI for the chosen k
boot_ari_ci <- quantile(boot_ari, probs = c(0.025, 0.975))
ari_label <- sprintf("k = %d: ARI = %.3f\n95%% CI [%.3f, %.3f]",
                     optimal_k, mean(boot_ari),
                     boot_ari_ci[1], boot_ari_ci[2])

p_dmin <- ggplot(tibble(k = 2:6, dmin = dmin_vals), aes(k, dmin)) +
  geom_line(linewidth = 0.8) + geom_point(size = 2) +
  geom_vline(xintercept = optimal_k, linetype = "dashed", color = "red") +
  annotate("text", x = 5.5, y = max(dmin_vals) * 0.9,
           label = ari_label, size = 2.5, hjust = 1, color = "red") +
  labs(x = "k", y = "Dmin", title = "Cluster Selection: Dmin Elbow") + THEME_PUB

ggsave(file.path(RPT_DIR, "supp_dmin_elbow.pdf"), p_dmin,
       width = 120, height = 90, units = "mm")

cat("  Supplementary Dmin elbow plot saved\n")
