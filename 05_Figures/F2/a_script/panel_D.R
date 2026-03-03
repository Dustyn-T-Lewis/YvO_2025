################################################################################
#   Figure 2 — Panel D: Concordance Scatter (Hexbin + Category Overlays)
#   Generates: panel_D_concordance.pdf/png,
#              c_data/panel_D/concordance.csv, concordance_stats.csv
################################################################################

if (!exists("dep_df")) source("05_Figures/F2/a_script/YvO_F2_setup.R")

message("Panel D: concordance scatter...")

# === 1. PREPARE DATA =========================================================

scatter_df <- dep_df %>%
  transmute(gene,
            logFC_Training_Young, logFC_Training_Old,
            pi_Young       = pi_score_Training_Young,
            pi_Old         = pi_score_Training_Old,
            pi_Interaction = pi_score_Interaction) %>%
  filter(!is.na(logFC_Training_Young), !is.na(logFC_Training_Old)) %>%
  left_join(imputation_df, by = "gene") %>%
  mutate(
    imputed      = replace_na(imputed, FALSE),
    significance = classify_proteins(pi_Young, pi_Old, pi_Interaction),
    quadrant = case_when(
      logFC_Training_Young > 0 & logFC_Training_Old > 0 ~ "Concordant Up",
      logFC_Training_Young < 0 & logFC_Training_Old < 0 ~ "Concordant Down",
      logFC_Training_Young > 0 & logFC_Training_Old < 0 ~ "Discordant (Young Up / Old Down)",
      TRUE ~ "Discordant (Young Down / Old Up)"),
    border_col   = ifelse(imputed, "black", "grey75"),
    point_size   = ifelse(significance == "NS", 1.8, 2.3),
    point_stroke = ifelse(significance == "NS", 0.6, 0.9),
    bubble_alpha = case_when(
      significance == "NS"             ~ 0.30,
      significance == "Interaction"    ~ 0.55,
      significance == "Sig Both"       ~ 0.75,
      TRUE                             ~ 0.85)
  )

# === 2. STATISTICS ============================================================

# Pearson r with 95% CI
cor_r   <- cor.test(scatter_df$logFC_Training_Young, scatter_df$logFC_Training_Old,
                    method = "pearson", conf.level = 0.95)
# Spearman rho with Fisher z-transformation CI
cor_rho <- cor.test(scatter_df$logFC_Training_Young, scatter_df$logFC_Training_Old,
                    method = "spearman", conf.level = 0.95)
n_obs   <- nrow(scatter_df)
rho_ci  <- tanh(atanh(cor_rho$estimate) + c(-1, 1) * qnorm(0.975) / sqrt(n_obs - 3))

# Sign concordance: |logFC| > 0.2 gate, BCa bootstrap 95% CI
concordance_set <- scatter_df %>%
  filter(abs(logFC_Training_Young) > 0.2 | abs(logFC_Training_Old) > 0.2)
sign_concordance <- mean(sign(concordance_set$logFC_Training_Young) ==
                         sign(concordance_set$logFC_Training_Old)) * 100

set.seed(42)
boot_sign_conc <- boot::boot(
  data = concordance_set,
  statistic = function(d, i)
    mean(sign(d$logFC_Training_Young[i]) == sign(d$logFC_Training_Old[i])) * 100,
  R = 10000)
boot_ci <- tryCatch(
  boot::boot.ci(boot_sign_conc, type = "bca", conf = 0.95)$bca[4:5],
  error = function(e) quantile(boot_sign_conc$t, c(0.025, 0.975)))

# === 3. SUBTITLE STRING =======================================================

sub_txt <- sprintf("n = %s | r = %.2f [%.2f, %.2f] | concordance = %.0f%%",
                   format(n_obs, big.mark = ","),
                   cor_r$estimate, cor_r$conf.int[1], cor_r$conf.int[2],
                   sign_concordance)

# === 4. LABELS & QUADRANT COUNTS =============================================

label_df <- scatter_df %>%
  filter(significance != "NS") %>%
  group_by(significance) %>%
  arrange(desc(abs(logFC_Training_Young) + abs(logFC_Training_Old))) %>%
  slice_head(n = 5) %>%
  ungroup() %>%
  mutate(label_fill     = SIG_LABEL_FILL[as.character(significance)],
         label_text_col = SIG_LABEL_TEXT[as.character(significance)],
         nudge_y = case_when(
           significance == "Interaction"    ~ -0.03,
           significance == "Sig Young only" ~  0.03,
           significance == "Sig Both"       ~  0.04,
           significance == "Sig Old only"   ~ -0.04,
           TRUE ~ 0))

q_df     <- scatter_df %>%
  mutate(q = case_when(
    logFC_Training_Young > 0 & logFC_Training_Old > 0 ~ "Q1",
    logFC_Training_Young < 0 & logFC_Training_Old < 0 ~ "Q3",
    logFC_Training_Young > 0 & logFC_Training_Old < 0 ~ "Q4",
    TRUE ~ "Q2"))
q_counts <- q_df %>% count(q) %>% deframe()
q_sig    <- q_df %>% filter(significance != "NS") %>% count(q) %>% deframe()
for (qq in c("Q1","Q2","Q3","Q4")) if (is.na(q_sig[qq])) q_sig[qq] <- 0

# === 5. BUILD PLOT ============================================================

ns_df  <- scatter_df %>% filter(significance == "NS")
sig_df <- scatter_df %>% filter(significance != "NS")

axis_max   <- max(abs(c(scatter_df$logFC_Training_Young,
                        scatter_df$logFC_Training_Old)), na.rm = TRUE) * 1.05
xlim_range <- c(-axis_max, axis_max)
ylim_range <- c(-axis_max, axis_max)

pD <- ggplot(mapping = aes(x = logFC_Training_Young, y = logFC_Training_Old)) +
  # Quadrant shading
  annotate("rect", xmin = 0, xmax = Inf,  ymin = 0,    ymax = Inf,
           fill = "#FFE0E0", alpha = 0.55) +
  annotate("rect", xmin = -Inf, xmax = 0, ymin = -Inf, ymax = 0,
           fill = "#FFE0E0", alpha = 0.55) +
  annotate("rect", xmin = 0, xmax = Inf,  ymin = -Inf, ymax = 0,
           fill = "#DCEEFF", alpha = 0.55) +
  annotate("rect", xmin = -Inf, xmax = 0, ymin = 0,    ymax = Inf,
           fill = "#DCEEFF", alpha = 0.55) +
  geom_hline(yintercept = 0, color = "grey60", linewidth = 0.2) +
  geom_vline(xintercept = 0, color = "grey60", linewidth = 0.2) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed",
              color = "black", linewidth = 0.3) +
  # Hexbin background: NS proteins
  geom_hex(data = ns_df, bins = 40, alpha = 0.7) +
  scale_fill_gradient(low = "grey95", high = "grey50", guide = "none") +
  ggnewscale::new_scale_fill() +
  # Significant category points
  geom_point(data = sig_df, aes(fill = significance), shape = 21,
             size = sig_df$point_size, color = sig_df$border_col,
             alpha = sig_df$bubble_alpha, stroke = sig_df$point_stroke) +
  scale_fill_manual(values = SIG_COLORS, name = "Significance",
                    guide = guide_legend(
                      order = 1,
                      override.aes = list(size = 3.5, alpha = 0.85,
                                          stroke = 0.6, color = "black"))) +
  # Gene labels
  geom_label_repel(data = label_df, aes(label = gene),
                   fill = label_df$label_fill, color = label_df$label_text_col,
                   nudge_y = label_df$nudge_y,
                   size = 2.8, fontface = "italic", max.overlaps = 40,
                   segment.size = 0.2, segment.color = "grey50",
                   min.segment.length = 0, show.legend = FALSE,
                   box.padding = 0.6, force = 3, force_pull = 0.5,
                   label.padding = unit(1.5, "pt"), label.r = unit(1, "pt"),
                   label.size = 0.15, seed = 42,
                   xlim = xlim_range * 0.9, ylim = ylim_range * 0.9) +
  # Quadrant corner labels
  annotate("label", x = Inf, y = Inf,
           label = sprintf("Concordant Up\u2002n = %s/%s", q_sig["Q1"], q_counts["Q1"]),
           hjust = 1, vjust = 1, size = 3.5, fontface = "bold",
           color = "#DC2626", fill = alpha("white", 0.9),
           label.padding = unit(2.5, "pt")) +
  annotate("label", x = -Inf, y = -Inf,
           label = sprintf("Concordant Down\u2002n = %s/%s", q_sig["Q3"], q_counts["Q3"]),
           hjust = 0, vjust = 0, size = 3.5, fontface = "bold",
           color = "#DC2626", fill = alpha("white", 0.9),
           label.padding = unit(2.5, "pt")) +
  annotate("label", x = -Inf, y = Inf,
           label = sprintf("Discordant\u2002n = %s/%s", q_sig["Q2"], q_counts["Q2"]),
           hjust = 0, vjust = 1, size = 3.5, fontface = "bold",
           color = "#2563EB", fill = alpha("white", 0.9),
           label.padding = unit(2.5, "pt")) +
  annotate("label", x = Inf, y = -Inf,
           label = sprintf("Discordant\u2002n = %s/%s", q_sig["Q4"], q_counts["Q4"]),
           hjust = 1, vjust = 0, size = 3.5, fontface = "bold",
           color = "#2563EB", fill = alpha("white", 0.9),
           label.padding = unit(2.5, "pt")) +
  coord_fixed(ratio = 1, xlim = xlim_range, ylim = ylim_range, expand = FALSE) +
  labs(title    = "Protein-Level Concordance of Training Response",
       subtitle = sub_txt,
       x = expression(log[2]*FC ~ "(Training Young)"),
       y = expression(log[2]*FC ~ "(Training Old)")) +
  FIG_THEME +
  theme(legend.position    = "bottom",
        legend.direction   = "horizontal",
        legend.box         = "horizontal",
        legend.key.size    = unit(3, "mm"),
        legend.text        = element_text(size = 9),
        legend.title       = element_text(size = 10, face = "bold"),
        legend.spacing.x   = unit(4, "mm"))

# Imputation border key strip
pD_imp_key <- ggplot(tibble(x = c(1, 3), y = c(0, 0),
                            label = c("Imputed", "Non-imputed")),
                     aes(x = x, y = y)) +
  annotate("text", x = 0.3, y = 0, label = "Border:",
           hjust = 0, size = 2.8, fontface = "bold", color = "grey30") +
  geom_point(shape = 21, size = 3.5, fill = "grey70",
             color = c("black", "grey75"), stroke = c(0.8, 1.2)) +
  geom_text(aes(label = label), hjust = -0.3, size = 2.8, color = "grey30") +
  scale_x_continuous(limits = c(-0.5, 5)) +
  theme_void() +
  theme(plot.margin = margin(0, 0, 0, 0))

pD_combined <- pD / pD_imp_key + plot_layout(heights = c(0.96, 0.04))

# === 6. SAVE OUTPUTS ==========================================================

ggsave(file.path(RPT_DIR, "panel_D_concordance.pdf"), pD_combined,
       width = 200, height = 200, units = "mm", device = cairo_pdf)
ggsave(file.path(RPT_DIR, "panel_D_concordance.png"), pD_combined,
       width = 200, height = 200, units = "mm", dpi = 300)

# Concordance data CSV
scatter_df %>%
  transmute(gene,
            logFC_Training_Young = round(logFC_Training_Young, 4),
            logFC_Training_Old   = round(logFC_Training_Old, 4),
            pi_score_Young       = round(pi_Young, 6),
            pi_score_Old         = round(pi_Old, 6),
            pi_score_Interaction = round(pi_Interaction, 6),
            significance         = as.character(significance),
            quadrant, imputed) %>%
  arrange(significance, desc(abs(logFC_Training_Young) + abs(logFC_Training_Old))) %>%
  write_csv(file.path(DAT_DIR, "panel_D", "concordance.csv"))

# Concordance statistics CSV
tibble(
  metric   = c("Pearson_r", "Spearman_rho", "Sign_concordance_pct"),
  estimate = c(cor_r$estimate, cor_rho$estimate, sign_concordance),
  ci_lower = c(cor_r$conf.int[1], rho_ci[1], boot_ci[1]),
  ci_upper = c(cor_r$conf.int[2], rho_ci[2], boot_ci[2]),
  p_value  = c(cor_r$p.value, cor_rho$p.value, NA_real_),
  n        = c(n_obs, n_obs, nrow(concordance_set)),
  note     = c("95% CI from cor.test()",
               "95% CI via Fisher z-transformation",
               "95% BCa bootstrap CI (10000 replicates, |logFC|>0.2 filter)")
) %>%
  write_csv(file.path(DAT_DIR, "panel_D", "concordance_stats.csv"))

message("  Panel D saved")
