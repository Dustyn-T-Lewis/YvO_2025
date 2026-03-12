################################################################################
#   Figure 1 — Supplementary S1.7: Pi-score Distributions
#
#   Requires from setup: dep_df, CTR_AXIS, THEME_PUB, RPT_DIR
#   Outputs: s17 (patchwork composite), p_hist, p_rank (ggplot objects)
#   Ref: Xiao et al. 2014 — pi-value
################################################################################

if (!exists("meta")) source("04_Figures/F1/a_script/YvO_F1_setup.R")

pi_long <- dep_df |>
  dplyr::select(gene, starts_with("pi_score_")) |>
  pivot_longer(starts_with("pi_score_"), names_to = "contrast", values_to = "pi_score") |>
  mutate(contrast = str_remove(contrast, "pi_score_")) |>
  filter(!is.na(pi_score))

p_hist <- ggplot(pi_long, aes(x = pi_score)) +
  geom_histogram(bins = 50, fill = "grey60", color = "white", linewidth = 0.2) +
  geom_vline(xintercept = 0.05, linetype = "dashed", color = "red", linewidth = 0.4) +
  facet_wrap(~ contrast, scales = "free_y", ncol = 2,
             labeller = labeller(contrast = CTR_AXIS)) +
  labs(x = expression(pi~"-score"), y = "Count") +
  THEME_PUB

pi_ranked <- pi_long |>
  group_by(contrast) |> arrange(pi_score) |>
  mutate(rank = row_number()) |> ungroup()

n_sig <- pi_long |>
  group_by(contrast) |>
  summarise(n = sum(pi_score < 0.05), .groups = "drop")

p_rank <- ggplot(pi_ranked, aes(x = rank, y = pi_score)) +
  geom_point(size = 0.2, alpha = 0.4, color = "grey40") +
  geom_hline(yintercept = 0.05, linetype = "dashed", color = "red", linewidth = 0.4) +
  geom_text(data = n_sig, aes(label = sprintf("n = %d", n)),
            x = Inf, y = 0.10, hjust = 1.1, vjust = 0, size = 2.5, color = "red") +
  facet_wrap(~ contrast, scales = "free_x", ncol = 2,
             labeller = labeller(contrast = CTR_AXIS)) +
  labs(x = "Protein rank", y = expression(pi~"-score")) +
  THEME_PUB

s17 <- p_hist / p_rank +
  plot_annotation(title = "S1.7  Pi-score distributions",
                  theme = theme(plot.title = element_text(face = "bold", size = 10)))

cat("Saved S1_7_pi_score_distributions.pdf\n")

ggsave(file.path(RPT_DIR, "supplementary", "S1_7_pi_score_distributions.pdf"), s17,
       width = 180, height = 200, units = "mm", device = pdf)
ggsave(file.path(RPT_DIR, "supplementary", "S1_7_pi_score_distributions.png"), s17,
       width = 180, height = 200, units = "mm", dpi = 300)
