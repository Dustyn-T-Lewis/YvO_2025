################################################################################
#   Figure 3 — Panel C: Reversal Scatter
#   logFC_Aging (x) vs logFC_Training_Old (y)
#   Generates: panel_C_reversal_scatter.pdf/png
#              + c_data/panel_C/reversal_scatter.csv
################################################################################
#
# STAT AUDIT (2026-02-27)
# ---------------------------------------------------------------------------
# 1. Pearson r (logFC_Aging vs logFC_Training_Old):
#    - cor.test() used, which provides 95% CI via Fisher z-transform.    PASS
#    - CI was computed but NOT displayed in subtitle.                     ISSUE
#      FIX: Add 95% CI for Pearson r to subtitle annotation.
#
# 2. Spearman rho:
#    - cor.test(method = "spearman") used.                               PASS
#    - Spearman rho CI not available from cor.test(); needs bootstrap.   ISSUE
#      FIX: Add bootstrap 95% CI for Spearman rho.
#
# 3. Reversal %:
#    - Proportion of proteins with opposite-sign logFC (all proteins).   PASS
#    - Bootstrap 95% CI on reversal proportion.                          PASS
#
# 4. Multiple testing:
#    - No correction needed: these are descriptive global statistics
#      (whole-proteome correlation), not per-protein tests.              PASS
#
# 5. Independence assumption:
#    - Proteins are not fully independent (co-regulated), which may
#      inflate correlation significance. This is standard practice in
#      proteomics; noted but not correctable without pathway-level
#      blocking.                                                         NOTE
# ---------------------------------------------------------------------------

if (!exists("dep_df")) source("04_Figures/F3/a_script/YvO_F3_setup.R")

# Capture Melov magnitude-reversal values before local reversal_pct overwrites
melov_rev_pct <- reversal_pct   # from setup: (d_pre - d_post) / d_pre * 100
melov_p       <- perm_pvalue    # from setup: permutation p-value

message("Panel C: reversal scatter...")

scatter_df <- dep_df %>%
  transmute(
    gene,
    logFC_Aging        = logFC_Aging,
    logFC_Training_Old = logFC_Training_Old,
    pi_Aging       = pi_score_Aging,
    pi_Training_Old = pi_score_Training_Old,
    pi_Reversal    = pi_score_Reversal
  ) %>%
  filter(!is.na(logFC_Aging), !is.na(logFC_Training_Old)) %>%
  left_join(imputation_df, by = "gene") %>%
  mutate(
    imputed = replace_na(imputed, FALSE),
    significance = classify_proteins_f3(pi_Aging, pi_Training_Old, pi_Reversal),
    quadrant = case_when(
      logFC_Aging > 0 & logFC_Training_Old < 0 ~ "Reversed (Aging Up / Training Down)",
      logFC_Aging < 0 & logFC_Training_Old > 0 ~ "Reversed (Aging Down / Training Up)",
      logFC_Aging > 0 & logFC_Training_Old > 0 ~ "Exacerbated Up",
      TRUE ~ "Exacerbated Down"
    ),
    border_col = ifelse(imputed, "black", "grey75"),
    point_size = case_when(
      significance == "NS" ~ 1.8,
      TRUE                 ~ 2.3
    ),
    point_stroke = case_when(
      significance == "NS" ~ 0.6,
      TRUE                 ~ 0.9
    ),
    bubble_alpha = case_when(
      significance == "NS"                ~ 0.30,
      significance == "Reversal"          ~ 0.55,
      significance == "Sig Both"          ~ 0.75,
      significance == "Sig Aging only"    ~ 0.85,
      significance == "Sig Training only" ~ 0.85,
      TRUE ~ 0.30
    )
  )

# Correlation
cor_r   <- cor.test(scatter_df$logFC_Aging, scatter_df$logFC_Training_Old, method = "pearson")
cor_rho <- cor.test(scatter_df$logFC_Aging, scatter_df$logFC_Training_Old, method = "spearman")

# AUDIT FIX: Bootstrap 95% CI for Spearman rho (not available from cor.test)
set.seed(42)
n_boot_rho <- 2000
boot_rho <- replicate(n_boot_rho, {
  idx <- sample(nrow(scatter_df), replace = TRUE)
  cor(scatter_df$logFC_Aging[idx], scatter_df$logFC_Training_Old[idx],
      method = "spearman")
})
rho_ci <- quantile(boot_rho, c(0.025, 0.975))

# Reversal %: proportion of proteins with opposite-sign logFC (all proteins)
reversal_pct <- mean(sign(scatter_df$logFC_Aging) !=
                     sign(scatter_df$logFC_Training_Old)) * 100

# Bootstrap 95% CI on reversal %
boot_rev_pct_c <- replicate(n_boot_rho, {
  idx <- sample(nrow(scatter_df), replace = TRUE)
  mean(sign(scatter_df$logFC_Aging[idx]) !=
       sign(scatter_df$logFC_Training_Old[idx])) * 100
})
rev_pct_ci_c <- quantile(boot_rev_pct_c, c(0.025, 0.975))

message(sprintf("  Pearson r = %.3f [%.3f, %.3f], p = %.2g",
                cor_r$estimate, cor_r$conf.int[1], cor_r$conf.int[2], cor_r$p.value))
message(sprintf("  Spearman rho = %.3f [%.3f, %.3f]",
                cor_rho$estimate, rho_ci[1], rho_ci[2]))
message(sprintf("  Reversal %% = %.1f%% [%.1f, %.1f]",
                reversal_pct, rev_pct_ci_c[1], rev_pct_ci_c[2]))

# Axis ranges
xlim_range <- range(scatter_df$logFC_Aging, na.rm = TRUE) * 1.05
ylim_range <- range(scatter_df$logFC_Training_Old, na.rm = TRUE) * 1.05

# Quadrant counts (total and sig)
q_df <- scatter_df %>%
  mutate(q = case_when(
    logFC_Aging > 0 & logFC_Training_Old < 0 ~ "BR",  # bottom-right = reversed
    logFC_Aging < 0 & logFC_Training_Old > 0 ~ "TL",  # top-left = reversed
    logFC_Aging > 0 & logFC_Training_Old > 0 ~ "TR",  # exacerbated up
    TRUE ~ "BL"                                         # exacerbated down
  ))
q_counts <- q_df %>% count(q) %>% deframe()
q_sig    <- q_df %>% filter(significance != "NS") %>% count(q) %>% deframe()
for (qq in c("BR","TL","TR","BL")) { if (is.na(q_sig[qq])) q_sig[qq] <- 0 }

# Labels for significant proteins
label_df <- scatter_df %>%
  filter(significance != "NS") %>%
  group_by(significance) %>%
  arrange(desc(abs(logFC_Aging) + abs(logFC_Training_Old))) %>%
  slice_head(n = 5) %>%
  ungroup() %>%
  mutate(
    label_fill     = SIG_LABEL_FILL[as.character(significance)],
    label_text_col = SIG_LABEL_TEXT[as.character(significance)]
  )

# Split NS vs significant for separate geom_point layers
ns_df  <- scatter_df %>% filter(significance == "NS")
sig_df <- scatter_df %>% filter(significance != "NS")

# Symmetric axis limits for coord_fixed
ax_max <- max(abs(c(xlim_range, ylim_range)))
sym_lim <- c(-ax_max, ax_max)

pC <- ggplot(mapping = aes(x = logFC_Aging, y = logFC_Training_Old)) +
  # Quadrant shading: blue = reversal (TL + BR), red = exacerbation (TR + BL)
  annotate("rect", xmin = 0, xmax = Inf,  ymin = -Inf, ymax = 0,
           fill = "#DCEEFF", alpha = 0.55) +   # bottom-right: reversed
  annotate("rect", xmin = -Inf, xmax = 0, ymin = 0, ymax = Inf,
           fill = "#DCEEFF", alpha = 0.55) +   # top-left: reversed
  annotate("rect", xmin = 0, xmax = Inf,  ymin = 0, ymax = Inf,
           fill = "#FFE0E0", alpha = 0.55) +   # top-right: exacerbated
  annotate("rect", xmin = -Inf, xmax = 0, ymin = -Inf, ymax = 0,
           fill = "#FFE0E0", alpha = 0.55) +   # bottom-left: exacerbated
  geom_hline(yintercept = 0, color = "grey60", linewidth = 0.2) +
  geom_vline(xintercept = 0, color = "grey60", linewidth = 0.2) +
  # Anti-diagonal reference (y = -x = perfect reversal)
  geom_abline(slope = -1, intercept = 0, linetype = "dashed",
              color = "black", linewidth = 0.3) +
  # NS proteins: subdued background layer
  geom_point(data = ns_df, shape = 21,
             color = "grey80", fill = "grey85", alpha = 0.3,
             size = 1.0, stroke = 0.2) +
  # Significant proteins: coloured foreground layer
  geom_point(data = sig_df, aes(fill = significance), shape = 21,
             size = sig_df$point_size, color = sig_df$border_col,
             alpha = sig_df$bubble_alpha, stroke = sig_df$point_stroke) +
  scale_fill_manual(values = SIG_COLORS, name = "Significance",
                    guide = guide_legend(
                      order = 1,
                      override.aes = list(size = 3.5, alpha = 0.85,
                                          stroke = 0.6, color = "black"))) +
  geom_label_repel(data = label_df, aes(label = gene),
                   fill = label_df$label_fill, color = label_df$label_text_col,
                   size = TXT_GENE, fontface = "italic",
                   max.overlaps = 40,
                   segment.size = 0.2, segment.color = "grey50",
                   min.segment.length = 0, show.legend = FALSE,
                   box.padding = 0.6, force = 3, force_pull = 0.5,
                   label.padding = unit(1.5, "pt"),
                   label.r = unit(1, "pt"),
                   label.size = 0.15, seed = 42) +
  # Quadrant labels
  annotate("label", x = Inf, y = -Inf,
           label = sprintf("Reversed  n = %s/%s", q_sig["BR"], q_counts["BR"]),
           hjust = 1, vjust = 0, size = TXT_QUADRANT, fontface = "bold",
           color = "#2563EB", fill = alpha("white", 0.9),
           label.padding = unit(2.5, "pt")) +
  annotate("label", x = -Inf, y = Inf,
           label = sprintf("Reversed  n = %s/%s", q_sig["TL"], q_counts["TL"]),
           hjust = 0, vjust = 1, size = TXT_QUADRANT, fontface = "bold",
           color = "#2563EB", fill = alpha("white", 0.9),
           label.padding = unit(2.5, "pt")) +
  annotate("label", x = Inf, y = Inf,
           label = sprintf("Exacerbated  n = %s/%s", q_sig["TR"], q_counts["TR"]),
           hjust = 1, vjust = 1, size = TXT_QUADRANT, fontface = "bold",
           color = "#DC2626", fill = alpha("white", 0.9),
           label.padding = unit(2.5, "pt")) +
  annotate("label", x = -Inf, y = -Inf,
           label = sprintf("Exacerbated  n = %s/%s", q_sig["BL"], q_counts["BL"]),
           hjust = 0, vjust = 0, size = TXT_QUADRANT, fontface = "bold",
           color = "#DC2626", fill = alpha("white", 0.9),
           label.padding = unit(2.5, "pt")) +
  # Melov magnitude-reversal annotation (along anti-diagonal)
  annotate("label",
           x = xlim_range[2] * 0.45, y = -xlim_range[2] * 0.45 - diff(ylim_range) * 0.06,
           label = sprintf("Magnitude reversal: %.1f%%, p = %.2f (n.s.)\n(Melov permutation, %d aging-sig. proteins)",
                           melov_rev_pct, melov_p, n_aging),
           hjust = 0.5, vjust = 1, size = TXT_STATS, fontface = "italic",
           color = "grey35", fill = alpha("white", 0.85),
           label.padding = unit(2, "pt")) +
  coord_fixed(ratio = 1, xlim = sym_lim, ylim = sym_lim, expand = FALSE) +
  labs(
    tag = "C",
    title = "Protein-Level Reversal of Aging by Training",
    subtitle = sprintf("r = %.2f [%.2f, %.2f], rho = %.2f [%.2f, %.2f]\n%s proteins | reversal = %.0f%% [%.0f, %.0f]",
                       cor_r$estimate, cor_r$conf.int[1], cor_r$conf.int[2],
                       cor_rho$estimate, rho_ci[1], rho_ci[2],
                       format(nrow(scatter_df), big.mark = ","),
                       reversal_pct, rev_pct_ci_c[1], rev_pct_ci_c[2]),
    x = expression(log[2]*FC ~ "(Aging)"),
    y = expression(log[2]*FC ~ "(Training Old)")
  ) +
  THEME_FIG +
  theme(
    legend.position = "none"
  )

ggsave(file.path(RPT_DIR, "panel_C_reversal_scatter.pdf"), pC,
       width = 200, height = 200, units = "mm", device = cairo_pdf)
ggsave(file.path(RPT_DIR, "panel_C_reversal_scatter.png"), pC,
       width = 200, height = 200, units = "mm", dpi = 300)

scatter_df %>%
  transmute(
    gene,
    logFC_Aging        = round(logFC_Aging, 4),
    logFC_Training_Old = round(logFC_Training_Old, 4),
    pi_score_Aging       = round(pi_Aging, 6),
    pi_score_Training_Old = round(pi_Training_Old, 6),
    pi_score_Reversal    = round(pi_Reversal, 6),
    significance         = as.character(significance),
    quadrant,
    imputed
  ) %>%
  arrange(significance, desc(abs(logFC_Aging) + abs(logFC_Training_Old))) %>%
  write_csv(file.path(DAT_DIR, "panel_C", "reversal_scatter.csv"))

message("  Panel C saved")
