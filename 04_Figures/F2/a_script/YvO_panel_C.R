################################################################################
#   Figure 2 — Panel C: Concordance Scatter
#   Generates: panel_C_concordance.pdf, panel_C_concordance.png
#              + c_data/panel_C/concordance.csv
################################################################################

if (!exists("dep_df")) source("04_Figures/F2/a_script/YvO_figure2_shared.R")

# ==============================================================================
# PANEL C — Concordance Scatter
# ==============================================================================

message("Panel C: concordance scatter...")

scatter_df <- dep_df %>%
  transmute(
    gene,
    logFC_Training_Young = logFC_Training_Young,
    logFC_Training_Old   = logFC_Training_Old,
    pi_Young     = pi_score_Training_Young,
    pi_Old       = pi_score_Training_Old,
    pi_Interaction = pi_score_Interaction
  ) %>%
  filter(!is.na(logFC_Training_Young), !is.na(logFC_Training_Old)) %>%
  left_join(imputation_df, by = "gene") %>%
  mutate(
    imputed = replace_na(imputed, FALSE),
    significance = classify_proteins(pi_Young, pi_Old, pi_Interaction),
    quadrant = case_when(
      logFC_Training_Young > 0 & logFC_Training_Old > 0 ~ "Concordant Up",
      logFC_Training_Young < 0 & logFC_Training_Old < 0 ~ "Concordant Down",
      logFC_Training_Young > 0 & logFC_Training_Old < 0 ~ "Discordant (Young Up / Old Down)",
      TRUE ~ "Discordant (Young Down / Old Up)"
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
      significance == "NS"             ~ 0.30,
      significance == "Interaction"    ~ 0.55,
      significance == "Sig Both"       ~ 0.75,
      significance == "Sig Young only" ~ 0.85,
      significance == "Sig Old only"   ~ 0.85,
      TRUE ~ 0.30
    )
  )

# Correlation
cor_r   <- cor.test(scatter_df$logFC_Training_Young, scatter_df$logFC_Training_Old, method = "pearson")
cor_rho <- cor.test(scatter_df$logFC_Training_Young, scatter_df$logFC_Training_Old, method = "spearman")

# Sign concordance among proteins with |logFC| > 0.2 in at least one contrast
concordance_set <- scatter_df %>%
  filter(abs(logFC_Training_Young) > 0.2 | abs(logFC_Training_Old) > 0.2)
sign_concordance <- mean(sign(concordance_set$logFC_Training_Young) ==
                         sign(concordance_set$logFC_Training_Old)) * 100

xlim_range <- c(-2.5, 2)
ylim_range <- c(-1, 2)

# Quadrant counts (total and sig)
q_df <- scatter_df %>%
  mutate(q = case_when(
    logFC_Training_Young > 0 & logFC_Training_Old > 0 ~ "Q1",
    logFC_Training_Young < 0 & logFC_Training_Old < 0 ~ "Q3",
    logFC_Training_Young > 0 & logFC_Training_Old < 0 ~ "Q4",
    TRUE ~ "Q2"
  ))
q_counts <- q_df %>% count(q) %>% deframe()
q_sig    <- q_df %>% filter(significance != "NS") %>% count(q) %>% deframe()
# Fill missing quadrants with 0
for (qq in c("Q1","Q2","Q3","Q4")) { if (is.na(q_sig[qq])) q_sig[qq] <- 0 }

# Labels for significant proteins
label_df <- scatter_df %>%
  filter(significance != "NS") %>%
  group_by(significance) %>%
  arrange(desc(abs(logFC_Training_Young) + abs(logFC_Training_Old))) %>%
  slice_head(n = 5) %>%
  ungroup() %>%
  mutate(
    label_fill     = SIG_LABEL_FILL[as.character(significance)],
    label_text_col = SIG_LABEL_TEXT[as.character(significance)],
    nudge_y = case_when(
      significance == "Interaction"    ~ -0.03,
      significance == "Sig Young only" ~  0.03,
      significance == "Sig Both"       ~  0.04,
      significance == "Sig Old only"   ~ -0.04,
      TRUE ~ 0
    )
  )

# --- Sort: NS first (bottom layer), significant on top ---
plot_order <- scatter_df %>% arrange(desc(as.integer(significance)))

pC <- ggplot(plot_order, aes(x = logFC_Training_Young, y = logFC_Training_Old)) +
  # Quadrant shading
  annotate("rect", xmin = 0, xmax = Inf,  ymin = 0, ymax = Inf,
           fill = "#FFE0E0", alpha = 0.55) +
  annotate("rect", xmin = -Inf, xmax = 0, ymin = -Inf, ymax = 0,
           fill = "#FFE0E0", alpha = 0.55) +
  annotate("rect", xmin = 0, xmax = Inf,  ymin = -Inf, ymax = 0,
           fill = "#DCEEFF", alpha = 0.55) +
  annotate("rect", xmin = -Inf, xmax = 0, ymin = 0, ymax = Inf,
           fill = "#DCEEFF", alpha = 0.55) +
  geom_hline(yintercept = 0, color = "grey60", linewidth = 0.2) +
  geom_vline(xintercept = 0, color = "grey60", linewidth = 0.2) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed",
              color = "black", linewidth = 0.3) +
  # Bubble layer: fill = significance, border = imputation
  geom_point(aes(fill = significance),
             shape = 21,
             size = plot_order$point_size,
             color = plot_order$border_col,
             alpha = plot_order$bubble_alpha,
             stroke = plot_order$point_stroke) +
  scale_fill_manual(values = SIG_COLORS, name = "Significance",
                    guide = guide_legend(
                      order = 1,
                      override.aes = list(size = 3.5, alpha = 0.85,
                                          stroke = 0.6, color = "black"))) +
  # Gene labels (colored label boxes, matching Panel D style)
  geom_label_repel(data = label_df, aes(label = gene),
                   fill = label_df$label_fill, color = label_df$label_text_col,
                   nudge_y = label_df$nudge_y,
                   size = 2.2, fontface = "bold",
                   max.overlaps = 30,
                   segment.size = 0.2, segment.color = "grey50",
                   min.segment.length = 0, show.legend = FALSE,
                   box.padding = 0.5, force = 3, force_pull = 0.5,
                   label.padding = unit(1.5, "pt"),
                   label.r = unit(1, "pt"),
                   label.size = 0.15, seed = 42,
                   xlim = c(xlim_range[1] * 0.9, xlim_range[2] * 0.9),
                   ylim = c(ylim_range[1] * 0.9, ylim_range[2] * 0.9)) +
  # Quadrant labels (flush to panel corners)
  annotate("label", x = Inf, y = Inf,
           label = sprintf("Concordant Up\u2002n = %s/%s", q_sig["Q1"], q_counts["Q1"]),
           hjust = 1, vjust = 1, size = 2.5, fontface = "bold",
           color = "#DC2626", fill = alpha("white", 0.9),
           label.padding = unit(2.5, "pt")) +
  annotate("label", x = -Inf, y = -Inf,
           label = sprintf("Concordant Down\u2002n = %s/%s", q_sig["Q3"], q_counts["Q3"]),
           hjust = 0, vjust = 0, size = 2.5, fontface = "bold",
           color = "#DC2626", fill = alpha("white", 0.9),
           label.padding = unit(2.5, "pt")) +
  annotate("label", x = -Inf, y = Inf,
           label = sprintf("Discordant\u2002n = %s/%s", q_sig["Q2"], q_counts["Q2"]),
           hjust = 0, vjust = 1, size = 2.5, fontface = "bold",
           color = "#2563EB", fill = alpha("white", 0.9),
           label.padding = unit(2.5, "pt")) +
  annotate("label", x = Inf, y = -Inf,
           label = sprintf("Discordant\u2002n = %s/%s", q_sig["Q4"], q_counts["Q4"]),
           hjust = 1, vjust = 0, size = 2.5, fontface = "bold",
           color = "#2563EB", fill = alpha("white", 0.9),
           label.padding = unit(2.5, "pt")) +
  coord_cartesian(xlim = xlim_range, ylim = ylim_range, expand = FALSE) +
  labs(
    title = "Protein-Level Concordance of Training Response",
    subtitle = sprintf("logFC Training Young vs Old | %s proteins | r = %.2f, rho = %.2f, concordance = %.0f%% | border = imputation status",
                       format(nrow(scatter_df), big.mark = ","),
                       cor_r$estimate, cor_rho$estimate, sign_concordance),
    x = expression(log[2]*FC ~ "(Training Young)"),
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

# Imputation key strip (mirrors Panel D database key strip)
pC_imp_key_strip <- ggplot(tibble(x = c(1, 3), y = c(0, 0),
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

pC_combined <- pC / pC_imp_key_strip + plot_layout(heights = c(0.96, 0.04))

ggsave(file.path(RPT_DIR, "panel_C_concordance.pdf"), pC_combined,
       width = 200, height = 200, units = "mm", device = pdf)
ggsave(file.path(RPT_DIR, "panel_C_concordance.png"), pC_combined,
       width = 200, height = 200, units = "mm", dpi = 300)

# Clean CSV (now includes imputed column)
scatter_df %>%
  transmute(
    gene,
    logFC_Training_Young = round(logFC_Training_Young, 4),
    logFC_Training_Old   = round(logFC_Training_Old, 4),
    pi_score_Young       = round(pi_Young, 6),
    pi_score_Old         = round(pi_Old, 6),
    pi_score_Interaction = round(pi_Interaction, 6),
    significance         = as.character(significance),
    quadrant,
    imputed
  ) %>%
  arrange(significance, desc(abs(logFC_Training_Young) + abs(logFC_Training_Old))) %>%
  write_csv(file.path(DAT_DIR, "panel_C", "concordance.csv"))

message("  Panel C saved")
