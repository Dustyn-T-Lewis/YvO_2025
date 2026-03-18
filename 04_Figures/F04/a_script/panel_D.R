# F2 Panel D: Concordance Scatter (Point + Category Overlays)
setwd(rprojroot::find_rstudio_root_file())
source("04_Figures/shared/style.R")
suppressPackageStartupMessages({
  library(tidyverse)
  library(ggrepel)
  library(boot)
})

PD_W <- 200

RPT <- "04_Figures/F04/b_reports"
DAT <- "04_Figures/F04/c_data"
dir.create(RPT, recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(DAT, "panel_D"), recursive = TRUE, showWarnings = FALSE)

pdf_device <- get_pdf_device()

dep_df <- read_csv("03_DEP/c_data/03_combined_results.csv", show_col_types = FALSE)
imputation_df <- read_csv("02_Imputation/c_data/02_mar_mnar_classification.csv",
                           show_col_types = FALSE) %>%
  transmute(gene, imputed = classification != "Complete")

scatter_df <- dep_df %>%
  transmute(
    gene,
    logFC_Training_Young, logFC_Training_Old,
    pi_Young       = pi_score_Training_Young,
    pi_Old         = pi_score_Training_Old,
    pi_Interaction = pi_score_Interaction
  ) %>%
  filter(!is.na(logFC_Training_Young), !is.na(logFC_Training_Old)) %>%
  left_join(imputation_df, by = "gene") %>%
  mutate(
    imputed = replace_na(imputed, FALSE),
    significance = classify_proteins_f2(pi_Young, pi_Old, pi_Interaction),
    quadrant = case_when(
      logFC_Training_Young > 0 & logFC_Training_Old > 0 ~ "Concordant Up",
      logFC_Training_Young < 0 & logFC_Training_Old < 0 ~ "Concordant Down",
      logFC_Training_Young > 0 & logFC_Training_Old < 0 ~ "Discordant (Young Up / Old Down)",
      TRUE ~ "Discordant (Young Down / Old Up)"
    ),
    border_col = "grey75",
    point_size = ifelse(significance == "NS", 1.8, 2.3),
    point_stroke = ifelse(significance == "NS", 0.6, 0.9),
    bubble_alpha = case_when(
      significance == "NS"             ~ 0.30,
      significance == "Interaction"    ~ 0.55,
      significance == "Sig Both"       ~ 0.75,
      TRUE ~ 0.85
    )
  )

# Correlation with 95% CIs
cor_r   <- cor.test(scatter_df$logFC_Training_Young, scatter_df$logFC_Training_Old,
                    method = "pearson",  conf.level = 0.95)
cor_rho <- cor.test(scatter_df$logFC_Training_Young, scatter_df$logFC_Training_Old,
                    method = "spearman", conf.level = 0.95)
n_obs <- nrow(scatter_df)
rho_z <- atanh(cor_rho$estimate)
rho_se <- 1 / sqrt(n_obs - 3)
rho_ci <- tanh(rho_z + c(-1, 1) * qnorm(0.975) * rho_se)

# Sign concordance on ALL proteins (no |logFC| filter — see decision-log)
sign_concordance <- mean(sign(scatter_df$logFC_Training_Young) ==
                         sign(scatter_df$logFC_Training_Old)) * 100

# Bootstrap 95% CI for sign concordance (BCa, 10000 replicates)
set.seed(42)
boot_sign_conc <- boot::boot(
  data = scatter_df,
  statistic = function(d, i) {
    mean(sign(d$logFC_Training_Young[i]) == sign(d$logFC_Training_Old[i])) * 100
  },
  R = 10000
)
boot_ci <- tryCatch(
  boot::boot.ci(boot_sign_conc, type = "bca", conf = 0.95)$bca[4:5],
  error = function(e) quantile(boot_sign_conc$t, c(0.025, 0.975))
)

sig_mask <- scatter_df$significance != "NS"
conc_sig <- mean(sign(scatter_df$logFC_Training_Young[sig_mask]) ==
                 sign(scatter_df$logFC_Training_Old[sig_mask])) * 100
n_sig <- sum(sig_mask)

concordance_stats <- tibble(
  metric = c("Pearson_r", "Spearman_rho", "Sign_concordance_pct", "Sign_concordance_sig_pct"),
  estimate = c(cor_r$estimate, cor_rho$estimate, sign_concordance, conc_sig),
  ci_lower = c(cor_r$conf.int[1], rho_ci[1], boot_ci[1], NA_real_),
  ci_upper = c(cor_r$conf.int[2], rho_ci[2], boot_ci[2], NA_real_),
  p_value  = c(cor_r$p.value, cor_rho$p.value, NA_real_, NA_real_),
  n        = c(n_obs, n_obs, n_obs, n_sig),
  note     = c("95% CI from cor.test()",
               "95% CI via Fisher z-transformation",
               "95% BCa bootstrap CI (10000 replicates, all proteins)",
               "Sign concordance among significant proteins only")
)
write_csv(concordance_stats, file.path(DAT, "panel_D", "concordance_stats.csv"))

axis_max <- max(abs(c(scatter_df$logFC_Training_Young, scatter_df$logFC_Training_Old)),
                na.rm = TRUE) * 1.05
xlim_range <- c(-axis_max, axis_max)
ylim_range <- c(-axis_max, axis_max)

q_df <- scatter_df %>%
  mutate(q = case_when(
    logFC_Training_Young > 0 & logFC_Training_Old > 0 ~ "Q1",
    logFC_Training_Young < 0 & logFC_Training_Old < 0 ~ "Q3",
    logFC_Training_Young > 0 & logFC_Training_Old < 0 ~ "Q4",
    TRUE ~ "Q2"
  ))
q_counts <- q_df %>% count(q) %>% deframe()
q_sig    <- q_df %>% filter(significance != "NS") %>% count(q) %>% deframe()
for (qq in c("Q1","Q2","Q3","Q4")) { if (is.na(q_sig[qq])) q_sig[qq] <- 0 }

label_df <- scatter_df %>%
  filter(significance != "NS") %>%
  group_by(significance) %>%
  arrange(desc(abs(logFC_Training_Young) + abs(logFC_Training_Old))) %>%
  slice_head(n = 5) %>%
  ungroup() %>%
  mutate(
    label_fill     = SIG_LABEL_FILL_F2[as.character(significance)],
    label_text_col = SIG_LABEL_TEXT_F2[as.character(significance)],
    nudge_y = case_when(
      significance == "Interaction"    ~ -0.03,
      significance == "Sig Young only" ~  0.03,
      significance == "Sig Both"       ~  0.04,
      significance == "Sig Old only"   ~ -0.04,
      TRUE ~ 0
    )
  )

ns_df  <- scatter_df %>% filter(significance == "NS")
sig_df <- scatter_df %>% filter(significance != "NS")

txt_gene <- scale_text(BASE_GENE, PD_W)
txt_quad <- scale_text(BASE_QUADRANT, PD_W)

pD <- ggplot(mapping = aes(x = logFC_Training_Young, y = logFC_Training_Old)) +
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
  geom_point(data = ns_df,
             aes(x = logFC_Training_Young, y = logFC_Training_Old),
             color = "grey80", fill = "grey85", shape = 21,
             size = 1.0, alpha = 0.3, stroke = 0.2) +
  geom_point(data = sig_df,
             aes(fill = significance),
             shape = 21,
             size = sig_df$point_size,
             color = sig_df$border_col,
             alpha = sig_df$bubble_alpha,
             stroke = sig_df$point_stroke) +
  scale_fill_manual(values = SIG_COLORS_F2, name = "Significance") +
  geom_label_repel(data = label_df, aes(label = gene),
                   fill = label_df$label_fill, color = label_df$label_text_col,
                   nudge_y = label_df$nudge_y,
                   size = txt_gene, fontface = "italic",
                   max.overlaps = 40,
                   segment.size = 0.2, segment.color = "grey50",
                   min.segment.length = 0, show.legend = FALSE,
                   box.padding = 0.6, force = 3, force_pull = 0.5,
                   label.padding = unit(1.5, "pt"),
                   label.r = unit(1, "pt"),
                   label.size = 0.15, seed = 42,
                   xlim = xlim_range * 0.9, ylim = ylim_range * 0.9) +
  annotate("label", x = Inf, y = Inf,
           label = sprintf("Concordant Up\u2002n = %s/%s", q_sig["Q1"], q_counts["Q1"]),
           hjust = 1, vjust = 1, size = txt_quad, fontface = "bold",
           color = "#DC2626", fill = alpha("white", 0.9),
           label.padding = unit(2.5, "pt")) +
  annotate("label", x = -Inf, y = -Inf,
           label = sprintf("Concordant Down\u2002n = %s/%s", q_sig["Q3"], q_counts["Q3"]),
           hjust = 0, vjust = 0, size = txt_quad, fontface = "bold",
           color = "#DC2626", fill = alpha("white", 0.9),
           label.padding = unit(2.5, "pt")) +
  annotate("label", x = -Inf, y = Inf,
           label = sprintf("Discordant\u2002n = %s/%s", q_sig["Q2"], q_counts["Q2"]),
           hjust = 0, vjust = 1, size = txt_quad, fontface = "bold",
           color = "#2563EB", fill = alpha("white", 0.9),
           label.padding = unit(2.5, "pt")) +
  annotate("label", x = Inf, y = -Inf,
           label = sprintf("Discordant\u2002n = %s/%s", q_sig["Q4"], q_counts["Q4"]),
           hjust = 1, vjust = 0, size = txt_quad, fontface = "bold",
           color = "#2563EB", fill = alpha("white", 0.9),
           label.padding = unit(2.5, "pt")) +
  coord_fixed(ratio = 1, xlim = xlim_range, ylim = ylim_range, expand = FALSE) +
  labs(
    title = "Protein-Level Concordance of Training Response",
    subtitle = sprintf("n = %s | \u03c1 = %.3f [%.3f, %.3f] | concordance = %.0f%% (all), %.0f%% (%d sig.)",
                       format(n_obs, big.mark = ","),
                       cor_rho$estimate, rho_ci[1], rho_ci[2],
                       sign_concordance, conc_sig, n_sig),
    x = expression(log[2]*FC ~ "(Training Young)"),
    y = expression(log[2]*FC ~ "(Training Old)")
  ) +
  FIG_THEME +
  theme(legend.position = "none")

ggsave(file.path(RPT, "panel_D_concordance.pdf"), pD,
       width = PD_W, height = PD_W, units = "mm", device = pdf_device)
ggsave(file.path(RPT, "panel_D_concordance.png"), pD,
       width = PD_W, height = PD_W, units = "mm", dpi = 300)

scatter_df %>%
  transmute(
    gene,
    logFC_Training_Young = round(logFC_Training_Young, 4),
    logFC_Training_Old   = round(logFC_Training_Old, 4),
    pi_score_Young       = round(pi_Young, 6),
    pi_score_Old         = round(pi_Old, 6),
    pi_score_Interaction = round(pi_Interaction, 6),
    significance         = as.character(significance),
    quadrant, imputed
  ) %>%
  arrange(significance, desc(abs(logFC_Training_Young) + abs(logFC_Training_Old))) %>%
  write_csv(file.path(DAT, "panel_D", "concordance.csv"))

message("F2 Panel D done")
