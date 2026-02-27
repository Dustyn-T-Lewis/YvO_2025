################################################################################
#   Figure 3 — Panel C: Reversal Scatter
#   logFC_Aging (x) vs logFC_Training_Old (y)
#   Generates: panel_C_reversal_scatter.pdf/png
#              + c_data/panel_C/reversal_scatter.csv
################################################################################

if (!exists("dep_df")) source("04_Figures/F3/a_script/YvO_F3_setup.R")

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

# Reversal %: proteins with opposite signs, |logFC| > 0.2 in at least one contrast
reversal_set <- scatter_df %>%
  filter(abs(logFC_Aging) > 0.2 | abs(logFC_Training_Old) > 0.2)
reversal_pct <- mean(sign(reversal_set$logFC_Aging) !=
                     sign(reversal_set$logFC_Training_Old)) * 100

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

# Sort: NS first (bottom layer)
plot_order <- scatter_df %>% arrange(desc(as.integer(significance)))

pC <- ggplot(plot_order, aes(x = logFC_Aging, y = logFC_Training_Old)) +
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
  geom_point(aes(fill = significance),
             shape = 21, size = plot_order$point_size,
             color = plot_order$border_col,
             alpha = plot_order$bubble_alpha,
             stroke = plot_order$point_stroke) +
  scale_fill_manual(values = SIG_COLORS, name = "Significance",
                    guide = guide_legend(
                      order = 1,
                      override.aes = list(size = 3.5, alpha = 0.85,
                                          stroke = 0.6, color = "black"))) +
  geom_label_repel(data = label_df, aes(label = gene),
                   fill = label_df$label_fill, color = label_df$label_text_col,
                   size = 2.2, fontface = "bold",
                   max.overlaps = 30,
                   segment.size = 0.2, segment.color = "grey50",
                   min.segment.length = 0, show.legend = FALSE,
                   box.padding = 0.5, force = 3, force_pull = 0.5,
                   label.padding = unit(1.5, "pt"),
                   label.r = unit(1, "pt"),
                   label.size = 0.15, seed = 42) +
  # Quadrant labels
  annotate("label", x = Inf, y = -Inf,
           label = sprintf("Reversed\u2002n = %s/%s", q_sig["BR"], q_counts["BR"]),
           hjust = 1, vjust = 0, size = 2.5, fontface = "bold",
           color = "#2563EB", fill = alpha("white", 0.9),
           label.padding = unit(2.5, "pt")) +
  annotate("label", x = -Inf, y = Inf,
           label = sprintf("Reversed\u2002n = %s/%s", q_sig["TL"], q_counts["TL"]),
           hjust = 0, vjust = 1, size = 2.5, fontface = "bold",
           color = "#2563EB", fill = alpha("white", 0.9),
           label.padding = unit(2.5, "pt")) +
  annotate("label", x = Inf, y = Inf,
           label = sprintf("Exacerbated\u2002n = %s/%s", q_sig["TR"], q_counts["TR"]),
           hjust = 1, vjust = 1, size = 2.5, fontface = "bold",
           color = "#DC2626", fill = alpha("white", 0.9),
           label.padding = unit(2.5, "pt")) +
  annotate("label", x = -Inf, y = -Inf,
           label = sprintf("Exacerbated\u2002n = %s/%s", q_sig["BL"], q_counts["BL"]),
           hjust = 0, vjust = 0, size = 2.5, fontface = "bold",
           color = "#DC2626", fill = alpha("white", 0.9),
           label.padding = unit(2.5, "pt")) +
  coord_cartesian(xlim = xlim_range, ylim = ylim_range, expand = FALSE) +
  labs(
    title = "Protein-Level Reversal of Aging by Training",
    subtitle = sprintf("logFC Aging vs Training Old | %s proteins | r = %.2f, \u03c1 = %.2f, reversal = %.0f%%",
                       format(nrow(scatter_df), big.mark = ","),
                       cor_r$estimate, cor_rho$estimate, reversal_pct),
    x = expression(log[2]*FC ~ "(Aging)"),
    y = expression(log[2]*FC ~ "(Training Old)")
  ) +
  THEME_PUB +
  theme(legend.position = "bottom",
        legend.direction = "horizontal",
        legend.box = "horizontal",
        legend.key.size = unit(3, "mm"),
        legend.text = element_text(size = 5.5),
        legend.title = element_text(size = 6, face = "bold"),
        legend.spacing.x = unit(4, "mm"))

# Imputation key strip
pC_imp_key <- ggplot(tibble(x = c(1, 3), y = c(0, 0),
                             label = c("Imputed", "Non-imputed")),
                      aes(x = x, y = y)) +
  annotate("text", x = 0.3, y = 0, label = "Border:",
           hjust = 0, size = 2.0, fontface = "bold", color = "grey30") +
  geom_point(shape = 21, size = 3.5, fill = "grey70",
             color = c("black", "grey75"), stroke = c(0.8, 1.2)) +
  geom_text(aes(label = label), hjust = -0.3, size = 1.8, color = "grey30") +
  scale_x_continuous(limits = c(-0.5, 5)) +
  theme_void() +
  theme(plot.margin = margin(0, 0, 0, 0))

pC_combined <- pC / pC_imp_key + plot_layout(heights = c(0.96, 0.04))

ggsave(file.path(RPT_DIR, "panel_C_reversal_scatter.pdf"), pC_combined,
       width = 200, height = 200, units = "mm", device = pdf)
ggsave(file.path(RPT_DIR, "panel_C_reversal_scatter.png"), pC_combined,
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
