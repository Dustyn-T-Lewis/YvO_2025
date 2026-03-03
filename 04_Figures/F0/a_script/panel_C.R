################################################################################
#   Figure 0 — Panel C: Vastus Lateralis Thickness
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

# ---- ANOVA subtitle ----
anova_tbl <- as.data.frame(stats_C_anova)
time_row  <- anova_tbl[anova_tbl$Effect == "Timepoint", ]
int_row   <- anova_tbl[anova_tbl$Effect == "Group:Timepoint", ]
anova_sub <- sprintf("Time: F(%g,%g) = %.2f, %s | Interaction: F(%g,%g) = %.2f, %s",
                     time_row$DFn, time_row$DFd, time_row$F, fmt_p(time_row$p),
                     int_row$DFn, int_row$DFd, int_row$F, fmt_p(int_row$p))

# ---- Left panel: Pre/Post bars per group ----
vl_bar_data <- group_summary %>%
  select(Group, Timepoint, Group_Time, vl_mean, vl_sem)

# Y-axis range for bracket positioning
y_max_left <- max(vl_bar_data$vl_mean + vl_bar_data$vl_sem) * 1.02

pC_left <- ggplot(vl_bar_data, aes(x = Group_Time, y = vl_mean, fill = Group_Time)) +
  geom_col(width = 0.65, color = "grey30", linewidth = 0.3) +
  geom_errorbar(aes(ymin = vl_mean - vl_sem, ymax = vl_mean + vl_sem),
                width = 0.2, linewidth = 0.4) +
  geom_text(aes(label = sprintf("%.2f", vl_mean)),
            vjust = -1.8, size = 2.5, color = "grey30") +
  # Paired brackets within Young
  geom_signif(
    comparisons  = list(c("Young_Pre", "Young_Post")),
    annotations  = fmt_p(stats_C_paired_young$p.value),
    y_position   = y_max_left + 0.15,
    textsize     = 2.5, tip_length = 0.01
  ) +
  # Paired brackets within Old
  geom_signif(
    comparisons  = list(c("Old_Pre", "Old_Post")),
    annotations  = fmt_p(stats_C_paired_old$p.value),
    y_position   = y_max_left + 0.15,
    textsize     = 2.5, tip_length = 0.01
  ) +
  scale_fill_manual(values = GROUP_FILL) +
  scale_x_discrete(labels = c("Young_Pre" = "Pre", "Young_Post" = "Post",
                               "Old_Pre"   = "Pre", "Old_Post"   = "Post")) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.18))) +
  labs(title    = "c",
       subtitle = anova_sub,
       y        = "VL thickness (cm)",
       x        = NULL) +
  THEME_PUB +
  theme(plot.title    = element_text(face = "bold", size = 14),
        plot.subtitle = element_text(size = 7, color = "grey40"),
        legend.position = "none")

# ---- Right panel: Delta bars ----
delta_vl <- pheno_wide %>%
  group_by(Group) %>%
  summarise(
    mean_delta = mean(delta_VL, na.rm = TRUE),
    sem_delta  = sd(delta_VL, na.rm = TRUE) / sqrt(sum(!is.na(delta_VL))),
    .groups    = "drop"
  )

delta_bar_colors <- c(Young = unname(GROUP_FILL["Young_Post"]),
                      Old   = unname(GROUP_FILL["Old_Post"]))

y_max_right <- max(delta_vl$mean_delta + delta_vl$sem_delta)

pC_right <- ggplot(delta_vl, aes(x = Group, y = mean_delta, fill = Group)) +
  geom_col(width = 0.55, color = "grey30", linewidth = 0.3) +
  geom_errorbar(aes(ymin = mean_delta - sem_delta,
                    ymax = mean_delta + sem_delta),
                width = 0.15, linewidth = 0.4) +
  geom_signif(
    comparisons = list(c("Young", "Old")),
    annotations = fmt_p(stats_C_delta$p.value),
    textsize    = 2.5, tip_length = 0.02,
    y_position  = y_max_right * 1.20
  ) +
  scale_fill_manual(values = delta_bar_colors) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.25))) +
  labs(y = expression(Delta * " VL (cm)"),
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
