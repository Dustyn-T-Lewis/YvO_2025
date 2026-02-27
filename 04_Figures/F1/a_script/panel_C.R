################################################################################
#   Figure 1 — Panel C: PCA Biplot + PERMANOVA + Legend Key
#
#   Requires from setup: imp_df, imp_mat, meta, samp_names, PCA_COLORS,
#                         PCA_SHAPES, BEST_IMP_METHOD, THEME_PUB, RPT_DIR,
#                         KEY_TEXT, KEY_TITLE
#   Outputs: pC (ggplot object), pC_key (ggplot legend key)
################################################################################

if (!exists("meta")) source("04_Figures/F1/a_script/YvO_F1_setup.R")

pca <- prcomp(t(imp_mat), center = TRUE, scale. = TRUE)
var_pct <- round(100 * summary(pca)$importance[2, 1:2], 1)

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
       x = sprintf("PC1 (%.1f%%)", var_pct[1]),
       y = sprintf("PC2 (%.1f%%)", var_pct[2])) +
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

cat("Panel C done\n")

pC_combined <- pC / pC_key + plot_layout(heights = c(0.83, 0.17))

ggsave(file.path(RPT_DIR, "panel_C_pca.pdf"), pC_combined,
       width = 130, height = 120, units = "mm", device = pdf)
ggsave(file.path(RPT_DIR, "panel_C_pca.png"), pC_combined,
       width = 130, height = 120, units = "mm", dpi = 300)
