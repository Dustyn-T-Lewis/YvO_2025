################################################################################
#   Figure 0 — Panel A: Training Volume
#   Bar chart comparing total training volume (kg) between Young and Old groups.
#   Training volume recorded at Post timepoint only (one value per subject).
################################################################################

if (!exists("tv_df")) source("04_Figures_v1/F0/a_script/YvO_F0_setup.R")

# Helper: format p-values for display
fmt_p <- function(p) {
  if (p < 0.001) return("p < 0.001")
  if (p < 0.01)  return(sprintf("p = %.3f", p))
  sprintf("p = %.2f", p)
}

# Bar colors: use Post-timepoint colors from GROUP_FILL
bar_colors <- c(Young = unname(GROUP_FILL["Young_Post"]),
                Old   = unname(GROUP_FILL["Old_Post"]))

pA <- ggplot(tv_df, aes(x = Group, y = tv, fill = Group)) +
  geom_bar(stat = "summary", fun = mean, width = 0.6, color = "black",
           linewidth = 0.3) +
  geom_errorbar(stat = "summary", fun.data = mean_se, width = 0.2,
                linewidth = 0.4) +
  geom_jitter(width = 0.15, size = 1.5, alpha = 0.5, shape = 16,
              color = "grey30") +
  geom_signif(
    comparisons = list(c("Young", "Old")),
    annotations = fmt_p(stats_A$p.value),
    textsize    = KEY_TEXT, tip_length = 0.02,
    y_position  = bracket_pos(tv_df$tv)
  ) +
  scale_fill_manual(values = bar_colors) +
  scale_x_discrete(labels = c(
    Young = sprintf("Younger (n = %d)", sum(tv_df$Group == "Young")),
    Old   = sprintf("Older (n = %d)",   sum(tv_df$Group == "Old")))) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.15)),
                     labels = label_comma()) +
  labs(title = "A  Total Training Volume",
       y     = "Total training volume (kg)",
       x     = NULL) +
  THEME_PUB +
  theme(legend.position = "none")

ggsave(file.path(RPT_DIR, "panel_A_training_volume.pdf"), pA,
       width = 70, height = 80, units = "mm", device = pdf)
ggsave(file.path(RPT_DIR, "panel_A_training_volume.png"), pA,
       width = 70, height = 80, units = "mm", dpi = 300)

cat("Panel A done\n")
