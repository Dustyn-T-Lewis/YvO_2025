################################################################################
#   Figure 1 — Panel A: CV% Violins (Inter-Individual Variability)
#
#   Requires from setup: norm_df, meta, samp_names, GROUP_FILL, THEME_PUB,
#                         RPT_DIR, KEY_TEXT
#   Outputs: pA (ggplot object)
#   Ref: Brenes 2024 — CV on linear scale
################################################################################
#
# STAT AUDIT (2026-02-27)
# ---------------------------------------------------------------------------
# 1. Test appropriateness:
#    - Wilcoxon rank-sum (unpaired) used for all four comparisons.    ISSUE
#      Pre vs Post comparisons (Young_Pre vs Young_Post, Old_Pre vs
#      Old_Post) involve paired samples (same subjects measured twice).
#      However, CV% is computed across subjects within a group — each CV
#      value is a per-protein statistic, not per-subject. Proteins are
#      independent observations. Unpaired Wilcoxon is therefore CORRECT for
#      comparing CV distributions across groups.                       PASS
#    - Wilcoxon is appropriate for non-normal, skewed CV distributions. PASS
#
# 2. Assumption checking:
#    - No distributional assumptions needed for Wilcoxon.              PASS
#    - Independence: each protein contributes one CV per group; proteins
#      may share biological regulation, but this is standard practice. PASS
#
# 3. Multiple comparison correction:
#    - Four pairwise Wilcoxon tests performed; no BH/Bonferroni applied
#      to the p-values used for bracket display.                       ISSUE
#      FIX: Apply BH correction across the 4 comparisons.
#
# 4. Effect sizes:
#    - Median CV reported on violins, but no formal effect size.       ISSUE
#      FIX: Add Cliff's delta (non-parametric effect size) for each pair.
#
# 5. Sample size adequacy:
#    - Each group has ~2100 protein CV values — ample power.           PASS
#
# 6. Confidence intervals:
#    - No CI on median CV per group.                                   ISSUE
#      FIX: Add bootstrap 95% CI on median CV for each group.
#
# 7. Reproducibility:
#    - No randomization in this panel; deterministic.                  PASS
# ---------------------------------------------------------------------------

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

# --- AUDIT FIX: Bootstrap 95% CI on median CV per group ---
set.seed(42)
boot_median_ci <- function(x, R = 2000, conf = 0.95) {
  meds <- replicate(R, median(sample(x, replace = TRUE)))
  qs   <- quantile(meds, c((1 - conf) / 2, (1 + conf) / 2))
  c(lower = unname(qs[1]), upper = unname(qs[2]))
}
cv_ci <- cv_df |>
  group_by(group) |>
  summarise(
    med   = median(cv),
    ci_lo = boot_median_ci(cv)[["lower"]],
    ci_hi = boot_median_ci(cv)[["upper"]],
    .groups = "drop"
  )
cat("Median CV% with 95% bootstrap CI:\n"); print(as.data.frame(cv_ci))

# Wilcoxon brackets — with BH correction across 4 comparisons
bracket_comps <- list(c("Young_Pre", "Young_Post"), c("Old_Pre", "Old_Post"),
                      c("Young_Pre", "Old_Pre"),    c("Young_Post", "Old_Post"))
bracket_pvals_raw <- sapply(bracket_comps, function(pair)
  wilcox.test(cv_df$cv[cv_df$group == pair[1]],
              cv_df$cv[cv_df$group == pair[2]])$p.value)
bracket_pvals <- p.adjust(bracket_pvals_raw, method = "BH")  # AUDIT FIX: BH correction
sig_idx <- which(bracket_pvals < 0.05)

# --- AUDIT FIX: Cliff's delta effect sizes ---
cliffs_delta <- function(x, y) {
  # Cliff's delta: proportion of (x_i > y_j) minus proportion of (x_i < y_j)
  nx <- length(x); ny <- length(y)
  d <- outer(x, y, function(a, b) sign(a - b))
  sum(d) / (nx * ny)
}
cliff_results <- data.frame(
  comparison = sapply(bracket_comps, paste, collapse = " vs "),
  p_raw      = bracket_pvals_raw,
  p_bh       = bracket_pvals,
  cliffs_d   = sapply(bracket_comps, function(pair)
    cliffs_delta(cv_df$cv[cv_df$group == pair[1]],
                 cv_df$cv[cv_df$group == pair[2]]))
)
cat("Wilcoxon tests with BH correction and Cliff's delta:\n")
print(cliff_results)

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

# AUDIT FIX: Use pre-computed BH-adjusted p-values instead of re-running
# unadjusted tests inside geom_signif
.signif_label <- function(p) {
  if (p < 0.001) "***" else if (p < 0.01) "**" else if (p < 0.05) "*" else "ns"
}

if (length(sig_idx) > 0) {
  sig_ypos <- bracket_base + (seq_along(sig_idx) - 1) * bracket_step
  for (i in seq_along(sig_idx)) {
    pA <- pA + geom_signif(
      comparisons = list(bracket_comps[[sig_idx[i]]]),
      y_position = sig_ypos[i], tip_length = 0.01,
      annotations = .signif_label(bracket_pvals[sig_idx[i]]),
      textsize = KEY_TEXT, size = 0.3)
  }
  pA <- pA + coord_cartesian(ylim = c(0, max(sig_ypos) + bracket_step * 1.5))
} else {
  pA <- pA + coord_cartesian(ylim = c(0, bracket_base + bracket_step))
}

# --- AUDIT: Export CI and effect-size tables ---
write.csv(as.data.frame(cv_ci),
          file.path(DAT_DIR, "audit_panel_A_median_cv_ci.csv"), row.names = FALSE)
write.csv(cliff_results,
          file.path(DAT_DIR, "audit_panel_A_wilcoxon_effects.csv"), row.names = FALSE)

cat("Panel A done\n")

ggsave(file.path(RPT_DIR, "panel_A_cv.pdf"), pA,
       width = 130, height = 60, units = "mm", device = pdf)
ggsave(file.path(RPT_DIR, "panel_A_cv.png"), pA,
       width = 130, height = 60, units = "mm", dpi = 300)
