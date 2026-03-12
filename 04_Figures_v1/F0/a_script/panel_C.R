################################################################################
#   Figure 0 — Panel C: Vastus Lateralis Thickness
#   Left:  4 bars (Young_Pre, Young_Post, Old_Pre, Old_Post) with paired brackets
#   Right: 2 delta bars (Young, Old) with between-group t-test bracket
################################################################################

if (!exists("meta")) source("04_Figures_v1/F0/a_script/YvO_F0_setup.R")

# Helper: format p-values for display
fmt_p <- function(p) {
  if (p < 0.001) return("p < 0.001")
  if (p < 0.01)  return(sprintf("p = %.3f", p))
  sprintf("p = %.2f", p)
}

# ---- ANOVA subtitle (Age + Time + Interaction) ----
anova_tbl <- as.data.frame(stats_C_anova)
anova_sub <- sprintf("Age %s   Time %s   Interaction %s",
                     fmt_p(anova_tbl$p[anova_tbl$Effect == "Group"]),
                     fmt_p(anova_tbl$p[anova_tbl$Effect == "Timepoint"]),
                     fmt_p(anova_tbl$p[anova_tbl$Effect == "Group:Timepoint"]))

# ---- Left panel: Pre/Post bars per group ----
sig_y_C <- bracket_pos(meta$VL_thick_cm)

pC_left <- ggplot(meta, aes(x = Group_Time, y = VL_thick_cm, fill = Group_Time)) +
  geom_bar(stat = "summary", fun = mean,
           width = 0.65, color = "black", linewidth = 0.3) +
  geom_errorbar(stat = "summary", fun.data = mean_se,
                width = 0.2, linewidth = 0.4) +
  geom_jitter(width = 0.12, size = 1.2, alpha = 0.5,
              shape = 21, color = "black", stroke = 0.3) +
  geom_signif(
    comparisons  = list(c("Young_Pre", "Young_Post")),
    annotations  = fmt_p(stats_C_paired_young$p.value),
    y_position   = sig_y_C,
    textsize     = KEY_TEXT, tip_length = 0.01
  ) +
  geom_signif(
    comparisons  = list(c("Old_Pre", "Old_Post")),
    annotations  = fmt_p(stats_C_paired_old$p.value),
    y_position   = sig_y_C,
    textsize     = KEY_TEXT, tip_length = 0.01
  ) +
  # Group annotations below x-axis
  annotate("text", x = 1.5, y = -Inf, label = "Younger",
           vjust = 4.2, fontface = "bold", size = 3.2, color = "grey25") +
  annotate("text", x = 3.5, y = -Inf, label = "Older",
           vjust = 4.2, fontface = "bold", size = 3.2, color = "grey25") +
  scale_fill_manual(values = GROUP_FILL) +
  scale_x_discrete(labels = c("Young_Pre" = "Pre", "Young_Post" = "Post",
                               "Old_Pre"   = "Pre", "Old_Post"   = "Post")) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.22))) +
  coord_cartesian(clip = "off") +
  labs(title    = "C  Vastus Lateralis Thickness",
       subtitle = anova_sub,
       y        = "VL thickness (cm)",
       x        = NULL) +
  THEME_PUB +
  theme(plot.margin      = margin(5, 5, 20, 5),
        legend.position  = "none")

# ---- Right panel: Delta bars ----
delta_bar_colors <- c(Young = unname(GROUP_FILL["Young_Post"]),
                      Old   = unname(GROUP_FILL["Old_Post"]))

pC_right <- ggplot(pheno_wide, aes(x = Group, y = delta_VL, fill = Group)) +
  geom_bar(stat = "summary", fun = mean,
           width = 0.55, color = "black", linewidth = 0.3) +
  geom_errorbar(stat = "summary", fun.data = mean_se,
                width = 0.15, linewidth = 0.4) +
  geom_jitter(width = 0.12, size = 1.2, alpha = 0.5,
              shape = 21, color = "black", stroke = 0.3) +
  geom_signif(
    comparisons = list(c("Young", "Old")),
    annotations = fmt_p(stats_C_delta$p.value),
    textsize    = KEY_TEXT, tip_length = 0.02,
    y_position  = bracket_pos(pheno_wide$delta_VL)
  ) +
  scale_fill_manual(values = delta_bar_colors) +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.25))) +
  labs(y = "change in VL thickness (cm)",
       x = NULL) +
  THEME_PUB +
  theme(legend.position = "none")

# ---- Combine ----
pC <- (pC_left | pC_right) + plot_layout(widths = c(0.65, 0.35))

ggsave(file.path(RPT_DIR, "panel_C_vl_thickness.pdf"), pC,
       width = 170, height = 80, units = "mm", device = pdf)
ggsave(file.path(RPT_DIR, "panel_C_vl_thickness.png"), pC,
       width = 170, height = 80, units = "mm", dpi = 300)

cat("Panel C done\n")
