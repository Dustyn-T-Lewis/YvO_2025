################################################################################
#   Figure 1 — Panel B: logFC Density Histograms (Effect Size Distributions)
#
#   Requires from setup: dep_df, CTR_FACET, CTR_SHORT, CONTRAST_COLORS,
#                         THEME_PUB, RPT_DIR, KEY_TEXT
#   Outputs: pB (ggplot object)
################################################################################

if (!exists("meta")) source("04_Figures/F1/a_script/YvO_F1_setup.R")

lfc_long <- dep_df |>
  dplyr::select(gene, starts_with("logFC_")) |>
  pivot_longer(starts_with("logFC_"), names_to = "contrast", values_to = "logFC") |>
  mutate(contrast = str_remove(contrast, "logFC_")) |>
  filter(!is.na(logFC), contrast %in% c("Aging", "Training_Young", "Training_Old"))

lfc_long$contrast <- factor(lfc_long$contrast,
                            levels = c("Aging", "Training_Young", "Training_Old"))

lfc_stats <- lfc_long |>
  group_by(contrast) |>
  summarise(
    med_abs_lfc = median(abs(logFC)),
    n_above_05  = sum(abs(logFC) > 0.5),
    .groups = "drop"
  )

# Pairwise Wilcoxon on |logFC| between contrasts (BH-adjusted)
pw_lfc <- pairwise.wilcox.test(abs(lfc_long$logFC), lfc_long$contrast,
                                p.adjust.method = "BH")

.lookup_pw <- function(a, b, mat) {
  if (a %in% rownames(mat) && b %in% colnames(mat) && !is.na(mat[a, b])) return(mat[a, b])
  if (b %in% rownames(mat) && a %in% colnames(mat) && !is.na(mat[b, a])) return(mat[b, a])
  NA_real_
}

lfc_stats$label <- sapply(seq_len(nrow(lfc_stats)), function(i) {
  ctr <- as.character(lfc_stats$contrast[i])
  others <- setdiff(c("Aging", "Training_Young", "Training_Old"), ctr)
  pw_lines <- sapply(others, function(o) {
    p <- .lookup_pw(ctr, o, pw_lfc$p.value)
    p_str <- if (is.na(p)) "p = NA" else if (p < 0.001) "p < 0.001" else sprintf("p = %.3f", p)
    sprintf("vs %s: %s", CTR_SHORT[o], p_str)
  })
  sprintf("Med.|logFC| = %.2f, n(>0.5) = %d\n%s",
          lfc_stats$med_abs_lfc[i], lfc_stats$n_above_05[i],
          paste(pw_lines, collapse = "\n"))
})

lfc_binwidth <- 4 / 50

pB <- ggplot(lfc_long, aes(x = logFC, fill = contrast)) +
  geom_histogram(bins = 50, color = "grey40", linewidth = 0.1, alpha = 0.85) +
  geom_density(aes(y = after_stat(count) * lfc_binwidth),
               alpha = 0.15, linewidth = 0.5, color = "grey20") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey40", linewidth = 0.3) +
  geom_text(data = lfc_stats, aes(label = label),
            x = Inf, y = Inf, hjust = 1.05, vjust = 1.2,
            size = KEY_TEXT, fontface = "bold", inherit.aes = FALSE, color = "grey25") +
  coord_cartesian(xlim = c(-2, 2)) +
  facet_wrap(~ contrast, ncol = 1, scales = "free_y",
             labeller = labeller(contrast = CTR_FACET)) +
  scale_fill_manual(values = CONTRAST_COLORS[c("Aging", "Training_Young", "Training_Old")]) +
  labs(title = "B  Effect Size Distribution (Log2FC)",
       subtitle = "Log2 fold-changes from limma (Interaction excluded)",
       x = expression(log[2]~FC), y = "Count") +
  THEME_PUB + theme(legend.position = "none")

# KS + Fligner distributional stats annotation
blunt_file <- "03_DEP/c_data/06_blunting_diagnostics.csv"
if (file.exists(blunt_file)) {
  blunt <- read_csv(blunt_file, show_col_types = FALSE)
  ks_row  <- blunt[blunt$test == "Kolmogorov-Smirnov", ]
  fl_row  <- blunt[blunt$test == "Fligner-Killeen", ]

  dist_label <- sprintf(
    "Young vs Old |logFC|:  KS D = %.2f, p < 10^{-30};  Fligner \u03C7\u00B2 = %.0f, p < 10^{-48}\n\u2192 Young shows larger, more variable training response",
    ks_row$statistic, fl_row$statistic)

  pB <- pB +
    labs(caption = dist_label) +
    theme(plot.caption = element_text(size = 5.5, color = "grey30",
                                       hjust = 0, face = "italic",
                                       margin = margin(t = 4)))
  cat(sprintf("  Added KS/Fligner annotation: D=%.3f, chi-sq=%.1f\n",
              ks_row$statistic, fl_row$statistic))
} else {
  cat("  blunting_diagnostics.csv not found — skipping annotation\n")
}

cat("Panel B done\n")

ggsave(file.path(RPT_DIR, "panel_B_logfc_density.pdf"), pB,
       width = 130, height = 110, units = "mm", device = pdf)
ggsave(file.path(RPT_DIR, "panel_B_logfc_density.png"), pB,
       width = 130, height = 110, units = "mm", dpi = 300)
