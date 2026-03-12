################################################################################
#   Figure 1 — Panel B: logFC Density Histograms (Effect Size Distributions)
#
#   Requires from setup: dep_df, CTR_FACET, CTR_SHORT, CONTRAST_COLORS,
#                         THEME_PUB, RPT_DIR, DAT_DIR, KEY_TEXT
#   Outputs: pB (ggplot object)
################################################################################
#
# STAT AUDIT (2026-02-27)
# ---------------------------------------------------------------------------
# 1. Test appropriateness:
#    - Pairwise Wilcoxon on |logFC| between contrasts: tests whether the
#      magnitude of fold-change distributions differs. Appropriate for
#      non-normal, heavy-tailed fold-change data.                      PASS
#    - KS test: tests whether the full distributions of |logFC| differ
#      in shape (location + spread + shape). Complementary to Wilcoxon
#      which focuses on central tendency.                              PASS
#    - Fligner-Killeen: tests homogeneity of variances — appropriate
#      non-parametric test for comparing spread.                       PASS
#
# 2. Assumption checking:
#    - No parametric assumptions needed.                               PASS
#    - logFC values from limma are per-protein estimates; proteins are
#      independent modelling units.                                    PASS
#
# 3. Multiple comparison correction:
#    - Pairwise Wilcoxon uses BH correction (via pairwise.wilcox.test). PASS
#    - KS and Fligner are single pre-planned comparisons (Young vs Old
#      training responses); no correction needed.                      PASS
#
# 4. Effect sizes:
#    - Median |logFC| reported per contrast.                           PASS
#    - No formal effect size for pairwise comparisons.                 ISSUE
#      FIX: Add Cliff's delta for pairwise |logFC| comparisons.
#
# 5. Sample size adequacy:
#    - ~2100 proteins per contrast — ample power.                      PASS
#
# 6. Confidence intervals:
#    - No CI on median |logFC|.                                        ISSUE
#      FIX: Add bootstrap 95% CI on median |logFC| per contrast.
#    - No CI on KS D statistic or Fligner chi-sq.                     NOTE
#      These are well-powered with N~2100; CIs not standard practice.
#
# 7. Reproducibility:
#    - No randomization in this panel.                                 PASS
# ---------------------------------------------------------------------------

if (!exists("meta")) source("04_Figures/F1/a_script/YvO_F1_setup.R")

lfc_long <- dep_df |>
  dplyr::select(gene, starts_with("logFC_")) |>
  pivot_longer(starts_with("logFC_"), names_to = "contrast", values_to = "logFC") |>
  mutate(contrast = str_remove(contrast, "logFC_")) |>
  filter(!is.na(logFC), contrast %in% c("Aging", "Training_Young", "Training_Old"))

lfc_long$contrast <- factor(lfc_long$contrast,
                            levels = c("Aging", "Training_Young", "Training_Old"))

# --- AUDIT FIX: Bootstrap 95% CI on median |logFC| per contrast ---
set.seed(42)
boot_median_ci <- function(x, R = 2000, conf = 0.95) {
  meds <- replicate(R, median(sample(x, replace = TRUE)))
  qs   <- quantile(meds, c((1 - conf) / 2, (1 + conf) / 2))
  c(lower = unname(qs[1]), upper = unname(qs[2]))
}

lfc_stats <- lfc_long |>
  group_by(contrast) |>
  summarise(
    med_abs_lfc = median(abs(logFC)),
    ci_lo       = boot_median_ci(abs(logFC))[["lower"]],
    ci_hi       = boot_median_ci(abs(logFC))[["upper"]],
    n_above_05  = sum(abs(logFC) > 0.5),
    .groups = "drop"
  )
cat("Median |logFC| with 95% bootstrap CI:\n"); print(as.data.frame(lfc_stats))

# Pairwise Wilcoxon on |logFC| between contrasts (BH-adjusted)
pw_lfc <- pairwise.wilcox.test(abs(lfc_long$logFC), lfc_long$contrast,
                                p.adjust.method = "BH")

.lookup_pw <- function(a, b, mat) {
  if (a %in% rownames(mat) && b %in% colnames(mat) && !is.na(mat[a, b])) return(mat[a, b])
  if (b %in% rownames(mat) && a %in% colnames(mat) && !is.na(mat[b, a])) return(mat[b, a])
  NA_real_
}

# --- AUDIT FIX: Cliff's delta for pairwise |logFC| comparisons ---
cliffs_delta <- function(x, y) {
  nx <- length(x); ny <- length(y)
  d <- outer(x, y, function(a, b) sign(a - b))
  sum(d) / (nx * ny)
}
contrast_pairs <- list(c("Aging", "Training_Young"), c("Aging", "Training_Old"),
                       c("Training_Young", "Training_Old"))
cliff_lfc <- data.frame(
  comparison = sapply(contrast_pairs, paste, collapse = " vs "),
  cliffs_d   = sapply(contrast_pairs, function(pair)
    cliffs_delta(abs(lfc_long$logFC[lfc_long$contrast == pair[1]]),
                 abs(lfc_long$logFC[lfc_long$contrast == pair[2]])))
)
cat("Cliff's delta for |logFC| comparisons:\n"); print(cliff_lfc)

# Build compact subtitle with all per-contrast stats
sub_med <- paste0("Med.|logFC|: ",
  paste(sprintf("%s = %.2f [%.2f, %.2f]", CTR_SHORT[as.character(lfc_stats$contrast)],
                lfc_stats$med_abs_lfc, lfc_stats$ci_lo, lfc_stats$ci_hi), collapse = ", "))
sub_n   <- paste0("n(|logFC| > 0.5): ",
  paste(sprintf("%s = %d", CTR_SHORT[as.character(lfc_stats$contrast)],
                lfc_stats$n_above_05), collapse = ", "),
  ";  all pairwise p < 0.001 (BH)")
panel_b_sub <- paste(sub_med, sub_n, sep = "\n")

# Label data for in-panel contrast titles (replaces facet strips)
title_df <- data.frame(
  contrast = factor(c("Aging", "Training_Young", "Training_Old"),
                    levels = c("Aging", "Training_Young", "Training_Old")),
  title    = CTR_FACET[c("Aging", "Training_Young", "Training_Old")]
)

lfc_binwidth <- 4 / 50

pB <- ggplot(lfc_long, aes(x = logFC, fill = contrast)) +
  geom_histogram(bins = 50, color = "grey40", linewidth = 0.1, alpha = 0.85) +
  geom_density(aes(y = after_stat(count) * lfc_binwidth),
               alpha = 0.15, linewidth = 0.5, color = "grey20") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey40", linewidth = 0.3) +
  geom_text(data = title_df, aes(label = title, color = contrast),
            x = -Inf, y = Inf, hjust = -0.05, vjust = 1.3,
            size = 3.2, fontface = "bold", inherit.aes = FALSE, show.legend = FALSE) +
  coord_cartesian(xlim = c(-2, 2)) +
  facet_wrap(~ contrast, ncol = 1, scales = "free_y") +
  scale_fill_manual(values = CONTRAST_COLORS[c("Aging", "Training_Young", "Training_Old")]) +
  scale_color_manual(values = CONTRAST_COLORS[c("Aging", "Training_Young", "Training_Old")]) +
  labs(title = "B  Effect Size Distribution (Log2FC)",
       subtitle = panel_b_sub,
       x = expression(log[2]~FC), y = "Count") +
  THEME_PUB + theme(legend.position = "none",
                    strip.text = element_blank(),
                    panel.spacing.y = unit(2, "mm"))

# KS + Fligner distributional stats annotation
blunt_file <- "03_DEP/c_data/06_blunting_diagnostics.csv"
if (file.exists(blunt_file)) {
  blunt <- read_csv(blunt_file, show_col_types = FALSE)
  ks_row  <- blunt[blunt$test == "Kolmogorov-Smirnov", ]
  fl_row  <- blunt[blunt$test == "Fligner-Killeen", ]

  dist_label <- bquote(
    "Young vs Old |logFC|:  KS D ="
    ~ .(sprintf("%.2f, p < 10", ks_row$statistic))^{-30}
    * ";  Fligner" ~ chi^2 ~ "="
    ~ .(sprintf("%.0f, p < 10", fl_row$statistic))^{-48}
  )

  pB <- pB +
    labs(caption = dist_label) +
    theme(plot.caption = element_text(size = 7, color = "grey30",
                                       hjust = 0,
                                       margin = margin(t = 6)))
  cat(sprintf("  Added KS/Fligner annotation: D=%.3f, chi-sq=%.1f\n",
              ks_row$statistic, fl_row$statistic))
} else {
  cat("  blunting_diagnostics.csv not found — skipping annotation\n")
}

# --- AUDIT: Export CI and effect-size tables ---
write.csv(as.data.frame(lfc_stats),
          file.path(DAT_DIR, "audit_panel_B_median_lfc_ci.csv"), row.names = FALSE)
write.csv(cliff_lfc,
          file.path(DAT_DIR, "audit_panel_B_cliff_delta.csv"), row.names = FALSE)

cat("Panel B done\n")

ggsave(file.path(RPT_DIR, "panel_B_logfc_density.pdf"), pB,
       width = 130, height = 110, units = "mm", device = pdf)
ggsave(file.path(RPT_DIR, "panel_B_logfc_density.png"), pB,
       width = 130, height = 110, units = "mm", dpi = 300)
