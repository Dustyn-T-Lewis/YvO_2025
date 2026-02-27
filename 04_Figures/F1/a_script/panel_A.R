################################################################################
#   Figure 1 — Panel A: CV% Violins (Inter-Individual Variability)
#
#   Requires from setup: norm_df, meta, samp_names, GROUP_FILL, THEME_PUB,
#                         RPT_DIR, KEY_TEXT
#   Outputs: pA (ggplot object)
#   Ref: Brenes 2024 — CV on linear scale
################################################################################

if (!exists("meta")) source("04_Figures/F1/a_script/YvO_F1_setup.R")

# CV on linear (not log) scale per Brenes 2024
lin_mat <- 2^as.matrix(norm_df[, samp_names])

cv_list <- lapply(levels(meta$group), function(g) {
  idx <- meta$sample_id[meta$group == g]
  sub <- lin_mat[, idx, drop = FALSE]
  cv_pct <- apply(sub, 1, function(x) {
    x <- x[!is.na(x)]
    if (length(x) < 2) return(NA_real_)
    sd(x) / mean(x) * 100
  })
  tibble(group = g, cv = cv_pct)
})
cv_df <- bind_rows(cv_list) |> filter(!is.na(cv))
cv_df$group <- factor(cv_df$group,
                      levels = c("Young_Pre", "Young_Post", "Old_Pre", "Old_Post"))

cv_med <- cv_df |> group_by(group) |> summarise(med = median(cv), .groups = "drop")

# Wilcoxon brackets — only significant (p < 0.05)
bracket_comps <- list(c("Young_Pre", "Young_Post"), c("Old_Pre", "Old_Post"),
                      c("Young_Pre", "Old_Pre"),    c("Young_Post", "Old_Post"))
bracket_pvals <- sapply(bracket_comps, function(pair)
  wilcox.test(cv_df$cv[cv_df$group == pair[1]],
              cv_df$cv[cv_df$group == pair[2]])$p.value)
sig_idx <- which(bracket_pvals < 0.05)

cv_ymax      <- quantile(cv_df$cv, 0.99)
bracket_base <- cv_ymax * 0.85
bracket_step <- cv_ymax * 0.10

pA <- ggplot(cv_df, aes(x = group, y = cv, fill = group)) +
  geom_violin(alpha = 0.7, linewidth = 0.3, scale = "width") +
  geom_boxplot(width = 0.15, outlier.size = 0.3, linewidth = 0.3, fill = "white") +
  geom_hline(yintercept = 25, linetype = "dashed", color = "grey40", linewidth = 0.3) +
  geom_label(data = cv_med, aes(x = group, y = med, label = sprintf("%.0f%%", med)),
             vjust = -0.3, size = 2.5, fontface = "bold", fill = alpha("white", 0.8),
             linewidth = 0.2, label.padding = unit(1, "pt")) +
  scale_fill_manual(values = GROUP_FILL) +
  scale_x_discrete(labels = c("Young Pre", "Young Post", "Old Pre", "Old Post")) +
  labs(title = "A  Inter-Individual Variability (CV%)",
       subtitle = "Protein-level CV on cycloess-normalized intensities (linear scale)",
       x = NULL, y = "CV (%)") +
  THEME_PUB + theme(legend.position = "none")

if (length(sig_idx) > 0) {
  sig_ypos <- bracket_base + (seq_along(sig_idx) - 1) * bracket_step
  for (i in seq_along(sig_idx)) {
    pA <- pA + geom_signif(
      comparisons = list(bracket_comps[[sig_idx[i]]]),
      y_position = sig_ypos[i], tip_length = 0.01,
      test = "wilcox.test", textsize = KEY_TEXT, size = 0.3,
      map_signif_level = c("***" = 0.001, "**" = 0.01, "*" = 0.05))
  }
  pA <- pA + coord_cartesian(ylim = c(0, max(sig_ypos) + bracket_step * 1.5))
} else {
  pA <- pA + coord_cartesian(ylim = c(0, bracket_base + bracket_step))
}

cat("Panel A done\n")

ggsave(file.path(RPT_DIR, "panel_A_cv.pdf"), pA,
       width = 130, height = 60, units = "mm", device = pdf)
ggsave(file.path(RPT_DIR, "panel_A_cv.png"), pA,
       width = 130, height = 60, units = "mm", dpi = 300)
