# F02 Supp Panel B: CV% Violins (Inter-Individual Variability)
# Faceted by Age (Young | Old), Pre/Post on x-axis. Median labels w/ bootstrap CIs.
# Outputs: pB (ggplot object), SUPP_panel_B_cv.{pdf,png}

# Assumes style.R sourced, packages loaded, norm_df/norm_meta/samp_names set by parent

PB_W <- 110
PB_H <- 80
RPT_PNG <- here::here("04_Figures_v2", "F02", "b_reports", "supp", "png", "panels")
RPT_PDF <- here::here("04_Figures_v2", "F02", "b_reports", "supp", "pdf", "panels")
DAT_DIR <- here::here("04_Figures_v2", "F02", "c_data")

meta <- norm_meta
meta$group <- factor(meta$group,
  levels = c("Young_Pre", "Young_Post", "Old_Pre", "Old_Post")
)

pdf_device <- get_pdf_device()

# CV on linear scale (Brenes 2024)
lin_mat <- 2^as.matrix(norm_df[, samp_names])

cv_list <- lapply(levels(meta$group), function(g) {
  idx <- meta$sample_id[meta$group == g]
  sub <- lin_mat[, idx, drop = FALSE]
  cv_pct <- apply(sub, 1, function(x) {
    x <- x[!is.na(x)]
    if (length(x) < 2) {
      return(NA_real_)
    }
    sd(x) / mean(x) * 100
  })
  tibble(group = g, cv = cv_pct)
})
cv_df <- bind_rows(cv_list) |> filter(!is.na(cv))
cv_df$group <- factor(cv_df$group,
  levels = c("Young_Pre", "Young_Post", "Old_Pre", "Old_Post")
)

cv_df$age <- factor(ifelse(grepl("Young", cv_df$group), "Young", "Old"),
  levels = c("Young", "Old")
)
cv_df$time <- factor(ifelse(grepl("Pre", cv_df$group), "Pre", "Post"),
  levels = c("Pre", "Post")
)

# Bootstrap 95% CI on median CV per group
set.seed(42)
# boot_median_ci() defined in shared/style.R

# Pairwise Wilcoxon tests, BH corrected (audit only — not on figure)
bracket_comps <- list(
  c("Young_Pre", "Young_Post"), c("Old_Pre", "Old_Post"),
  c("Young_Pre", "Old_Pre"), c("Young_Post", "Old_Post")
)
bracket_pvals_raw <- sapply(bracket_comps, function(pair) {
  wilcox.test(
    cv_df$cv[cv_df$group == pair[1]],
    cv_df$cv[cv_df$group == pair[2]]
  )$p.value
})
bracket_pvals <- p.adjust(bracket_pvals_raw, method = "BH")

cliffs_delta <- function(x, y) {
  nx <- length(x)
  ny <- length(y)
  d <- outer(x, y, function(a, b) sign(a - b))
  sum(d) / (nx * ny)
}
cliff_results <- data.frame(
  comparison = sapply(bracket_comps, paste, collapse = " vs "),
  p_raw = bracket_pvals_raw,
  p_bh = bracket_pvals,
  cliffs_d = sapply(bracket_comps, function(pair) {
    cliffs_delta(
      cv_df$cv[cv_df$group == pair[1]],
      cv_df$cv[cv_df$group == pair[2]]
    )
  })
)

cv_ci <- cv_df |>
  group_by(age, time, group) |>
  summarise(
    med = median(cv),
    ci_lo = boot_median_ci(cv)[["lower"]],
    ci_hi = boot_median_ci(cv)[["upper"]],
    cv_max = max(cv),
    .groups = "drop"
  )

n_prot <- nrow(norm_df)
grand_med <- median(cv_df$cv)
grand_ci <- boot_median_ci(cv_df$cv)

# Delta median CV per age group (Post - Pre)
delta_cv <- cv_ci |>
  select(age, time, med) |>
  pivot_wider(names_from = time, values_from = med) |>
  mutate(
    delta = Post - Pre,
    arrow_label = sprintf("%+.1f%%", delta)
  )

# Arrow annotation data: one arrow per age facet from Pre median to Post median
arrow_df <- delta_cv |>
  mutate(
    x = 1, xend = 2,
    y_mid = (Pre + Post) / 2
  )

sub_txt <- sprintf(
  paste0(
    "Group CV%% per Age \u00d7 Time | %s proteins\n",
    "median %.0f%% [%.0f\u2013%.0f CI] | Y %+.1f%%, O %+.1f%%"
  ),
  format(n_prot, big.mark = ","), grand_med, grand_ci[1], grand_ci[2],
  delta_cv$delta[delta_cv$age == "Young"],
  delta_cv$delta[delta_cv$age == "Old"]
)

# Consolidated single plot: group on x-axis, age shading in background
GROUP_LABELS <- c(
  Young_Pre = "Pre", Young_Post = "Post",
  Old_Pre = "Pre", Old_Post = "Post"
)

# Arrow annotations spanning Pre→Post within each age group
arrow_df_single <- delta_cv |>
  mutate(
    x = ifelse(age == "Young", 1, 3),
    xend = ifelse(age == "Young", 2, 4),
    y_mid = (Pre + Post) / 2
  )

pB <- ggplot(cv_df, aes(x = group, y = cv, fill = group)) +
  # Age group labels at top
  annotate("text",
    x = 1.5, y = Inf, label = "Young", vjust = 1.3,
    size = 2.0, fontface = "bold", color = "grey25"
  ) +
  annotate("text",
    x = 3.5, y = Inf, label = "Old", vjust = 1.3,
    size = 2.0, fontface = "bold", color = "grey25"
  ) +
  geom_violin(alpha = 0.5, linewidth = 0.3, color = "black", scale = "width") +
  geom_quasirandom(aes(color = group),
    alpha = 0.15, size = 0.5,
    width = 0.25, show.legend = FALSE
  ) +
  geom_boxplot(
    width = 0.15, outlier.shape = NA, linewidth = 0.3,
    color = "black", fill = "white", coef = 0
  ) +
  geom_hline(
    yintercept = 25, linetype = "dashed", color = "grey50",
    linewidth = 0.4
  ) +
  geom_label(
    data = cv_ci,
    aes(
      x = group, y = cv_max + 3,
      label = sprintf("%.0f%% [%.0f\u2013%.0f]", med, ci_lo, ci_hi)
    ),
    size = scale_text(BASE_COUNT - 1.5, PB_W),
    fontface = "bold", fill = alpha("white", 0.8),
    linewidth = 0.2, label.padding = unit(1.0, "pt"),
    hjust = 0.5
  ) +
  geom_segment(
    data = arrow_df_single,
    aes(x = x, xend = xend, y = Pre, yend = Post),
    inherit.aes = FALSE, color = "grey30",
    arrow = arrow(length = unit(1.5, "mm"), type = "closed"),
    linewidth = 0.6
  ) +
  geom_label(
    data = arrow_df_single,
    aes(x = (x + xend) / 2, y = y_mid, label = arrow_label),
    inherit.aes = FALSE, size = scale_text(BASE_COUNT - 1.5, PB_W),
    fontface = "bold.italic", fill = alpha("white", 0.85),
    label.padding = unit(1.5, "pt"), linewidth = 0.2,
    color = "grey30"
  ) +
  scale_fill_manual(values = GROUP_FILL) +
  scale_color_manual(values = GROUP_FILL) +
  scale_x_discrete(labels = GROUP_LABELS) +
  coord_cartesian(ylim = c(0, max(cv_ci$cv_max) + 15)) +
  labs(
    title = "Inter-Individual Variability (CV%)",
    subtitle = sub_txt,
    x = NULL, y = "CV (%)",
    tag = "b"
  ) +
  FIG_THEME +
  theme(
    legend.position = "none",
    plot.title = element_text(margin = margin(b = 0)),
    plot.subtitle = element_text(
      size = FIG_SUBTITLE_SIZE - 1.0,
      face = "bold.italic", color = "grey40",
      margin = margin(t = 0, b = 1)
    ),
    plot.margin = margin(t = 0, r = 5.5, b = 5.5, l = 5.5)
  )

write.csv(as.data.frame(cv_ci),
  file.path(DAT_DIR, "SUPP_panel_B_median_cv_ci.csv"),
  row.names = FALSE
)
write.csv(cliff_results,
  file.path(DAT_DIR, "SUPP_panel_B_wilcoxon.csv"),
  row.names = FALSE
)

ggsave(file.path(RPT_PNG, "SUPP_panel_B_cv.png"), pB,
  width = PB_W, height = PB_H, units = "mm", dpi = 300
)
ggsave(file.path(RPT_PDF, "SUPP_panel_B_cv.pdf"), pB,
  width = PB_W, height = PB_H, units = "mm", device = pdf_device
)

# Export for composite
pSB_title <- "Inter-Individual Variability (CV%)"
pSB_subtitle <- sub_txt
pSB_legend <- NULL
pB <- strip_for_composite(pB)

message("F02 Supp Panel B done")
