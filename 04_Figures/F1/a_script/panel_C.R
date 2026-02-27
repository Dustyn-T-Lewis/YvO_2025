################################################################################
#   Figure 1 — Panel C: PCA Biplot + PERMANOVA + Legend Key
#
#   Requires from setup: imp_df, imp_mat, meta, samp_names, PCA_COLORS,
#                         PCA_SHAPES, BEST_IMP_METHOD, THEME_PUB, RPT_DIR,
#                         DAT_DIR, KEY_TEXT, KEY_TITLE
#   Outputs: pC (ggplot object), pC_key (ggplot legend key)
################################################################################
#
# STAT AUDIT (2026-02-27)
# ---------------------------------------------------------------------------
# 1. Test appropriateness:
#    - PCA on centered + scaled imputed matrix: appropriate for
#      exploratory dimensionality reduction.                           PASS
#    - PERMANOVA (adonis2) with age * time: appropriate for testing
#      multivariate group differences on distance matrix.              PASS
#    - `by = "terms"` gives sequential (Type I) sums of squares.
#      For balanced designs this is fine; for unbalanced, Type III
#      (by = "margin") is more robust. Design may be slightly
#      unbalanced (N differs between Young/Old). Worth checking but
#      unlikely to change conclusions with this effect magnitude.      NOTE
#
# 2. Assumption checking:
#    - PERMANOVA assumes homogeneity of multivariate dispersions.      ISSUE
#      FIX: Add betadisper test (Anderson 2006) to verify.
#
# 3. Multiple comparison correction:
#    - Three PERMANOVA terms (age, time, age:time) tested
#      simultaneously as part of one model — not a multiple-testing
#      issue (analogous to ANOVA terms).                               PASS
#
# 4. Effect sizes:
#    - R-squared reported for each term.                               PASS
#
# 5. Sample size adequacy:
#    - N = 31 samples for 3 PERMANOVA terms: adequate. Anderson (2017)
#      recommends >=10 per group; we have 15-16 per age group.         PASS
#
# 6. Confidence intervals:
#    - No CI on % variance explained by PC1, PC2.                     ISSUE
#      FIX: Add bootstrap CI on variance proportions.
#    - No CI on PERMANOVA R-squared.                                  NOTE
#      Not standard; permutation p-value suffices.
#
# 7. Reproducibility:
#    - set.seed(42) before PERMANOVA.                                  PASS
# ---------------------------------------------------------------------------

if (!exists("meta")) source("04_Figures/F1/a_script/YvO_F1_setup.R")

pca <- prcomp(t(imp_mat), center = TRUE, scale. = TRUE)
var_pct <- round(100 * summary(pca)$importance[2, 1:2], 1)

# --- AUDIT FIX: Bootstrap 95% CI on % variance explained (PC1, PC2) ---
set.seed(42)
n_boot <- 1000
n_prot <- nrow(imp_mat)
boot_var <- matrix(NA_real_, nrow = n_boot, ncol = 2)
for (b in seq_len(n_boot)) {
  idx <- sample(n_prot, replace = TRUE)
  bp  <- prcomp(t(imp_mat[idx, ]), center = TRUE, scale. = TRUE)
  boot_var[b, ] <- 100 * summary(bp)$importance[2, 1:2]
}
var_ci <- data.frame(
  PC      = c("PC1", "PC2"),
  var_pct = var_pct,
  ci_lo   = apply(boot_var, 2, quantile, 0.025),
  ci_hi   = apply(boot_var, 2, quantile, 0.975)
)
cat("PCA variance explained with 95% bootstrap CI:\n"); print(var_ci)

pca_df <- as.data.frame(pca$x[, 1:2]) |>
  mutate(sample_id = rownames(pca$x)) |>
  left_join(meta, by = "sample_id")

# PERMANOVA
dist_mat <- dist(t(imp_mat))
set.seed(42)
perm_res <- adonis2(dist_mat ~ age * time, data = meta,
                    permutations = 999, by = "terms")

.fmt_p <- function(p) {
  if (is.na(p)) return("NA")
  if (p < 0.001) return("< 0.001")
  sprintf("%.3f", p)
}
perm_terms <- c("age", "time", "age:time")
perm_r2 <- perm_res[perm_terms, "R2"]
perm_pv <- perm_res[perm_terms, "Pr(>F)"]
perm_label <- sprintf(
  "PERMANOVA (999 perm.)\nAge  R\u00b2 = %.3f,  p = %s\nTime  R\u00b2 = %.3f,  p = %s\nAge\u00d7Time  R\u00b2 = %.3f,  p = %s",
  perm_r2[1], .fmt_p(perm_pv[1]),
  perm_r2[2], .fmt_p(perm_pv[2]),
  perm_r2[3], .fmt_p(perm_pv[3]))
# --- AUDIT FIX: betadisper — homogeneity of multivariate dispersions ---
# Anderson (2006): PERMANOVA is sensitive to dispersion differences;
# betadisper tests the assumption that groups have similar spread.
bd_age  <- betadisper(dist_mat, meta$age)
bd_time <- betadisper(dist_mat, meta$time)
bd_grp  <- betadisper(dist_mat, meta$group)
bd_age_p  <- permutest(bd_age,  pairwise = FALSE, permutations = 999)$tab$`Pr(>F)`[1]
bd_time_p <- permutest(bd_time, pairwise = FALSE, permutations = 999)$tab$`Pr(>F)`[1]
bd_grp_p  <- permutest(bd_grp,  pairwise = FALSE, permutations = 999)$tab$`Pr(>F)`[1]
cat(sprintf("betadisper p-values: Age=%.3f, Time=%.3f, Group=%.3f\n",
            bd_age_p, bd_time_p, bd_grp_p))
if (bd_age_p < 0.05 || bd_time_p < 0.05)
  warning("Heterogeneous dispersions detected — interpret PERMANOVA with caution")

cat("PERMANOVA results:\n"); print(perm_res)

pC <- ggplot(pca_df, aes(x = PC1, y = PC2, color = group, shape = group)) +
  stat_ellipse(aes(fill = group), geom = "polygon",
               alpha = 0.10, level = 0.80, show.legend = FALSE) +
  stat_ellipse(aes(group = group), level = 0.80, linewidth = 0.4,
               linetype = "dashed", show.legend = FALSE) +
  geom_point(size = 2.0, alpha = 0.85) +
  annotate("text", x = -Inf, y = Inf, label = perm_label,
           hjust = -0.05, vjust = 1.2,
           size = KEY_TEXT, fontface = "bold", color = "grey25") +
  scale_color_manual(values = PCA_COLORS) +
  scale_fill_manual(values = PCA_COLORS) +
  scale_shape_manual(values = PCA_SHAPES) +
  labs(title = "C  Principal Component Analysis (PCA)",
       subtitle = sprintf("%s proteins (HPA-filtered, cycloess, %s-imputed); n = %d samples",
                          format(nrow(imp_df), big.mark = ","), BEST_IMP_METHOD, nrow(meta)),
       x = sprintf("PC1 (%.1f%% [%.1f, %.1f])", var_pct[1], var_ci$ci_lo[1], var_ci$ci_hi[1]),
       y = sprintf("PC2 (%.1f%% [%.1f, %.1f])", var_pct[2], var_ci$ci_lo[2], var_ci$ci_hi[2])) +
  THEME_PUB + theme(legend.position = "none")

# PCA key — 2x2 grid with group-colored backgrounds
pca_leg <- tibble(
  col   = c(0, 0.50, 0, 0.50),
  row   = c(1, 1, 0, 0),
  label = c("Young Pre", "Young Post", "Old Pre", "Old Post"),
  group = c("Young_Pre", "Young_Post", "Old_Pre", "Old_Post"),
  shape = c(16, 17, 16, 17)
) |> mutate(fill = PCA_COLORS[group])

pC_key <- ggplot(pca_leg) +
  geom_rect(aes(xmin = col, xmax = col + 0.38, ymin = row - 0.38, ymax = row + 0.38),
            fill = alpha(pca_leg$fill, 0.30), color = "grey70", linewidth = 0.2) +
  geom_point(aes(x = col + 0.05, y = row),
             shape = pca_leg$shape, color = pca_leg$fill, size = 2) +
  geom_text(aes(x = col + 0.10, y = row, label = label),
            hjust = 0, size = KEY_TEXT, fontface = "bold", color = "grey15") +
  annotate("text", x = 0, y = 1.75, label = "Group \u00d7 Time:",
           hjust = 0, size = KEY_TITLE, fontface = "bold", color = "grey25") +
  scale_x_continuous(limits = c(-0.02, 0.92)) +
  scale_y_continuous(limits = c(-0.45, 1.92)) +
  theme_void() +
  theme(plot.margin = margin(0, 5.5, 0, 5.5))

# --- AUDIT: Export CI and betadisper results ---
write.csv(var_ci,
          file.path(DAT_DIR, "audit_panel_C_pca_variance_ci.csv"), row.names = FALSE)
betadisper_results <- data.frame(
  factor  = c("age", "time", "group"),
  p_value = c(bd_age_p, bd_time_p, bd_grp_p),
  significant = c(bd_age_p < 0.05, bd_time_p < 0.05, bd_grp_p < 0.05)
)
write.csv(betadisper_results,
          file.path(DAT_DIR, "audit_panel_C_betadisper.csv"), row.names = FALSE)

cat("Panel C done\n")

pC_combined <- pC / pC_key + plot_layout(heights = c(0.83, 0.17))

ggsave(file.path(RPT_DIR, "panel_C_pca.pdf"), pC_combined,
       width = 130, height = 120, units = "mm", device = pdf)
ggsave(file.path(RPT_DIR, "panel_C_pca.png"), pC_combined,
       width = 130, height = 120, units = "mm", dpi = 300)
