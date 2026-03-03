################################################################################
#   Figure 0 — Panel B: DXA Lean Body Mass
#   Left:  4 bars (Young_Pre, Young_Post, Old_Pre, Old_Post) with paired brackets
#   Right: 2 delta bars (Young, Old) with between-group t-test bracket
################################################################################

if (!exists("meta")) source("04_Figures/F0/a_script/YvO_F0_setup.R")

# Helper: format p-values for display
fmt_p <- function(p) {
  if (p < 0.001) return("p < 0.001")
  if (p < 0.01)  return(sprintf("p = %.3f", p))
  sprintf("p = %.2f", p)
}

# ---- ANOVA subtitle (Age + Time + Interaction) ----
anova_tbl <- as.data.frame(stats_B_anova)
anova_sub <- sprintf("Age %s   Time %s   Interaction %s",
                     fmt_p(anova_tbl$p[anova_tbl$Effect == "Group"]),
                     fmt_p(anova_tbl$p[anova_tbl$Effect == "Timepoint"]),
                     fmt_p(anova_tbl$p[anova_tbl$Effect == "Group:Timepoint"]))

# ---- Left panel: Pre/Post bars per group ----
# Summary subset for mean labels
group_summary_B <- group_summary %>%
  select(Group, Timepoint, Group_Time, dxa_mean, dxa_sem)

# Y-axis range for bracket positioning
y_max_left <- max(group_summary_B$dxa_mean + group_summary_B$dxa_sem) * 1.02

pB_left <- ggplot(meta, aes(x = Group_Time, y = DXA_LBM_kg, fill = Group_Time)) +
  geom_bar(stat = "summary", fun = mean,
           width = 0.65, color = "grey30", linewidth = 0.3) +
  geom_errorbar(stat = "summary", fun.data = mean_se,
                width = 0.2, linewidth = 0.4) +
  geom_jitter(width = 0.12, size = 1.2, alpha = 0.5,
              shape = 21, color = "black", stroke = 0.3) +
  geom_text(data = group_summary_B,
            aes(x = Group_Time, y = 0, label = sprintf("%.1f", dxa_mean)),
            vjust = 1.5, size = 2.5, color = "grey30") +
  # Paired brackets within Young
  geom_signif(
    comparisons  = list(c("Young_Pre", "Young_Post")),
    annotations  = fmt_p(stats_B_paired_young$p.value),
    y_position   = y_max_left + 2,
    textsize     = 2.5, tip_length = 0.01
  ) +
  # Paired brackets within Old
  geom_signif(
    comparisons  = list(c("Old_Pre", "Old_Post")),
    annotations  = fmt_p(stats_B_paired_old$p.value),
    y_position   = y_max_left + 2,
    textsize     = 2.5, tip_length = 0.01
  ) +
  scale_fill_manual(values = GROUP_FILL) +
  scale_x_discrete(labels = c("Young_Pre" = "Pre", "Young_Post" = "Post",
                               "Old_Pre"   = "Pre", "Old_Post"   = "Post")) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.18))) +
  labs(title    = "b",
       subtitle = anova_sub,
       y        = "DXA LSTM (kg)",
       x        = NULL) +
  THEME_PUB +
  theme(plot.title    = element_text(face = "bold", size = 14),
        plot.subtitle = element_text(size = 7, color = "grey40"),
        legend.position = "none")

# ---- Right panel: Delta bars ----
delta_bar_colors <- c(Young = unname(GROUP_FILL["Young_Post"]),
                      Old   = unname(GROUP_FILL["Old_Post"]))

y_max_right <- pheno_wide %>%
  group_by(Group) %>%
  summarise(m = mean(delta_DXA, na.rm = TRUE),
            s = sd(delta_DXA, na.rm = TRUE) / sqrt(sum(!is.na(delta_DXA))),
            .groups = "drop") %>%
  summarise(ymax = max(m + s)) %>%
  pull(ymax)

pB_right <- ggplot(pheno_wide, aes(x = Group, y = delta_DXA, fill = Group)) +
  geom_bar(stat = "summary", fun = mean,
           width = 0.55, color = "grey30", linewidth = 0.3) +
  geom_errorbar(stat = "summary", fun.data = mean_se,
                width = 0.15, linewidth = 0.4) +
  geom_jitter(width = 0.12, size = 1.2, alpha = 0.5,
              shape = 21, color = "black", stroke = 0.3) +
  geom_signif(
    comparisons = list(c("Young", "Old")),
    annotations = fmt_p(stats_B_delta$p.value),
    textsize    = 2.5, tip_length = 0.02,
    y_position  = y_max_right * 1.20
  ) +
  scale_fill_manual(values = delta_bar_colors) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.25))) +
  labs(y = "change in DXA LSTM (kg)",
       x = NULL) +
  THEME_PUB +
  theme(legend.position = "none")

# ---- Combine ----
pB <- (pB_left | pB_right) + plot_layout(widths = c(0.65, 0.35))

ggsave(file.path(RPT_DIR, "panel_B_dxa_lbm.pdf"), pB,
       width = 170, height = 80, units = "mm", device = pdf)
ggsave(file.path(RPT_DIR, "panel_B_dxa_lbm.png"), pB,
       width = 170, height = 80, units = "mm", dpi = 300)

cat("Panel B done\n")
