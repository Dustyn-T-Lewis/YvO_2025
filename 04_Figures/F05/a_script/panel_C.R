# F3 Panel C: Reversal Concordance Scatter
setwd(rprojroot::find_rstudio_root_file())
source("04_Figures/shared/style.R")

suppressPackageStartupMessages({
  library(tidyverse)
  library(ggrepel)
  library(boot)
})

PW <- 200; PH <- 200
RPT <- "04_Figures/F05/b_reports"
DAT <- "04_Figures/F05/c_data"
dir.create(RPT, recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(DAT, "panel_C"), recursive = TRUE, showWarnings = FALSE)

pdf_device <- get_pdf_device()

dep_df <- read_csv("03_DEP/c_data/03_combined_results.csv", show_col_types = FALSE)

imputation_df <- read_csv("02_Imputation/c_data/02_mar_mnar_classification.csv",
                           show_col_types = FALSE) %>%
  transmute(gene, imputed = classification != "Complete")

melov_df <- read_csv("04_Figures/F05/c_data/reversal_tests/melov_permutation.csv",
                      show_col_types = FALSE)

scatter_df <- dep_df %>%
  transmute(gene,
            logFC_Aging, logFC_Training_Old,
            pi_Aging        = pi_score_Aging,
            pi_Training_Old = pi_score_Training_Old,
            pi_Reversal     = pi_score_Reversal) %>%
  filter(!is.na(logFC_Aging), !is.na(logFC_Training_Old)) %>%
  left_join(imputation_df, by = "gene") %>%
  mutate(
    imputed      = replace_na(imputed, FALSE),
    significance = classify_proteins_f3(pi_Aging, pi_Training_Old, pi_Reversal),
    quadrant = case_when(
      logFC_Aging > 0 & logFC_Training_Old < 0 ~ "Reversed (Aging Up / Training Down)",
      logFC_Aging < 0 & logFC_Training_Old > 0 ~ "Reversed (Aging Down / Training Up)",
      logFC_Aging > 0 & logFC_Training_Old > 0 ~ "Exacerbated Up",
      TRUE ~ "Exacerbated Down"),
    border_col   = ifelse(imputed, "black", "grey75"),
    point_size   = ifelse(significance == "NS", 1.8, 2.3),
    point_stroke = ifelse(significance == "NS", 0.6, 0.9),
    bubble_alpha = case_when(
      significance == "NS"                ~ 0.30,
      significance == "Reversal"          ~ 0.55,
      significance == "Sig Both"          ~ 0.75,
      TRUE                                ~ 0.85)
  )

cor_r   <- cor.test(scatter_df$logFC_Aging, scatter_df$logFC_Training_Old,
                    method = "pearson", conf.level = 0.95)
cor_rho <- cor.test(scatter_df$logFC_Aging, scatter_df$logFC_Training_Old,
                    method = "spearman", conf.level = 0.95)
n_obs   <- nrow(scatter_df)
rho_ci  <- tanh(atanh(cor_rho$estimate) + c(-1, 1) * qnorm(0.975) / sqrt(n_obs - 3))

reversal_pct <- mean(sign(scatter_df$logFC_Aging) !=
                     sign(scatter_df$logFC_Training_Old)) * 100

set.seed(42)
boot_rev <- boot::boot(
  data = scatter_df,
  statistic = function(d, i)
    mean(sign(d$logFC_Aging[i]) != sign(d$logFC_Training_Old[i])) * 100,
  R = 10000)
rev_ci <- tryCatch(
  boot::boot.ci(boot_rev, type = "bca", conf = 0.95)$bca[4:5],
  error = function(e) quantile(boot_rev$t, c(0.025, 0.975)))

# Correlations on significant proteins only
sig_mask <- scatter_df$significance != "NS"
n_sig    <- sum(sig_mask)
cor_r_sig   <- cor.test(scatter_df$logFC_Aging[sig_mask],
                        scatter_df$logFC_Training_Old[sig_mask],
                        method = "pearson", conf.level = 0.95)
cor_rho_sig <- cor.test(scatter_df$logFC_Aging[sig_mask],
                        scatter_df$logFC_Training_Old[sig_mask],
                        method = "spearman", conf.level = 0.95)
# Bootstrap CI for Spearman rho (more precise at small n than Fisher z)
set.seed(43)
boot_rho_sig <- boot::boot(
  data = scatter_df[sig_mask, ],
  statistic = function(d, i)
    cor(d$logFC_Aging[i], d$logFC_Training_Old[i], method = "spearman"),
  R = 10000
)
rho_sig_ci <- tryCatch(
  boot::boot.ci(boot_rho_sig, type = "bca", conf = 0.95)$bca[4:5],
  error = function(e) quantile(boot_rho_sig$t, c(0.025, 0.975))
)
rev_sig    <- mean(sign(scatter_df$logFC_Aging[sig_mask]) !=
                   sign(scatter_df$logFC_Training_Old[sig_mask])) * 100

message(sprintf("  Pearson r = %.3f [%.3f, %.3f], p = %.2g",
                cor_r$estimate, cor_r$conf.int[1], cor_r$conf.int[2], cor_r$p.value))
message(sprintf("  Spearman rho = %.3f [%.3f, %.3f]",
                cor_rho$estimate, rho_ci[1], rho_ci[2]))
message(sprintf("  Reversal %% = %.1f%% [%.1f, %.1f]",
                reversal_pct, rev_ci[1], rev_ci[2]))
message(sprintf("  Sig-only (n=%d): r = %.3f, rho = %.3f, reversal = %.1f%%",
                n_sig, cor_r_sig$estimate, cor_rho_sig$estimate, rev_sig))

sub_txt <- sprintf(
  "All (n = %s): r = %.2f [%.2f, %.2f], \u03c1 = %.2f [%.2f, %.2f] | reversal = %.0f%%\nSig. (n = %d): r = %.2f [%.2f, %.2f], \u03c1 = %.2f [%.2f, %.2f] | reversal = %.0f%%",
  format(n_obs, big.mark = ","),
  cor_r$estimate, cor_r$conf.int[1], cor_r$conf.int[2],
  cor_rho$estimate, rho_ci[1], rho_ci[2], reversal_pct,
  n_sig,
  cor_r_sig$estimate, cor_r_sig$conf.int[1], cor_r_sig$conf.int[2],
  cor_rho_sig$estimate, rho_sig_ci[1], rho_sig_ci[2], rev_sig)

txt_gene <- scale_text(BASE_GENE, PW)
txt_quad <- scale_text(BASE_QUADRANT, PW)
txt_stat <- scale_text(BASE_STAT, PW)

label_df <- scatter_df %>%
  filter(significance != "NS") %>%
  group_by(significance) %>%
  arrange(desc(abs(logFC_Aging) + abs(logFC_Training_Old))) %>%
  slice_head(n = 5) %>%
  ungroup() %>%
  mutate(label_fill     = SIG_LABEL_FILL_F3[as.character(significance)],
         label_text_col = SIG_LABEL_TEXT_F3[as.character(significance)])

q_df <- scatter_df %>%
  mutate(q = case_when(
    logFC_Aging > 0 & logFC_Training_Old < 0 ~ "BR",
    logFC_Aging < 0 & logFC_Training_Old > 0 ~ "TL",
    logFC_Aging > 0 & logFC_Training_Old > 0 ~ "TR",
    TRUE ~ "BL"))
q_counts <- q_df %>% count(q) %>% deframe()
q_sig    <- q_df %>% filter(significance != "NS") %>% count(q) %>% deframe()
for (qq in c("BR", "TL", "TR", "BL")) if (is.na(q_sig[qq])) q_sig[qq] <- 0


ns_df  <- scatter_df %>% filter(significance == "NS")
sig_df <- scatter_df %>% filter(significance != "NS")

axis_max <- max(abs(c(scatter_df$logFC_Aging, scatter_df$logFC_Training_Old)),
                na.rm = TRUE) * 1.05
xlim_range <- c(-axis_max, axis_max)
ylim_range <- c(-axis_max, axis_max)

melov_rev_pct <- melov_df$reversal_pct
melov_p       <- melov_df$p_value
melov_n       <- melov_df$n_aging_proteins

pC <- ggplot(mapping = aes(x = logFC_Aging, y = logFC_Training_Old)) +
  annotate("rect", xmin = 0, xmax = Inf,  ymin = -Inf, ymax = 0,
           fill = "#DCEEFF", alpha = 0.55) +
  annotate("rect", xmin = -Inf, xmax = 0, ymin = 0, ymax = Inf,
           fill = "#DCEEFF", alpha = 0.55) +
  annotate("rect", xmin = 0, xmax = Inf,  ymin = 0, ymax = Inf,
           fill = "#FFE0E0", alpha = 0.55) +
  annotate("rect", xmin = -Inf, xmax = 0, ymin = -Inf, ymax = 0,
           fill = "#FFE0E0", alpha = 0.55) +
  geom_hline(yintercept = 0, color = "grey60", linewidth = 0.2) +
  geom_vline(xintercept = 0, color = "grey60", linewidth = 0.2) +
  geom_abline(slope = -1, intercept = 0, linetype = "dashed",
              color = "black", linewidth = 0.3) +
  geom_point(data = ns_df, color = "grey70", size = 0.8, alpha = 0.25) +
  geom_point(data = sig_df, aes(fill = significance), shape = 21,
             size = sig_df$point_size, color = sig_df$border_col,
             alpha = sig_df$bubble_alpha, stroke = sig_df$point_stroke) +
  scale_fill_manual(values = SIG_COLORS_F3, name = "Significance") +
  geom_label_repel(data = label_df, aes(label = gene),
                   fill = label_df$label_fill, color = label_df$label_text_col,
                   size = txt_gene, fontface = "italic", max.overlaps = 40,
                   segment.size = 0.2, segment.color = "grey50",
                   min.segment.length = 0, show.legend = FALSE,
                   box.padding = 0.6, force = 3, force_pull = 0.5,
                   label.padding = unit(1.5, "pt"), label.r = unit(1, "pt"),
                   label.size = 0.15, seed = 42) +
  annotate("label", x = xlim_range[2] * 0.95, y = ylim_range[1] * 0.95,
           label = sprintf("Reversed\nn = %s/%s", q_sig["BR"], q_counts["BR"]),
           hjust = 1, vjust = 0, size = txt_quad, fontface = "bold",
           color = "#2563EB", fill = alpha("white", 0.92),
           label.padding = unit(3, "pt")) +
  annotate("label", x = xlim_range[1] * 0.95, y = ylim_range[2] * 0.95,
           label = sprintf("Reversed\nn = %s/%s", q_sig["TL"], q_counts["TL"]),
           hjust = 0, vjust = 1, size = txt_quad, fontface = "bold",
           color = "#2563EB", fill = alpha("white", 0.92),
           label.padding = unit(3, "pt")) +
  annotate("label", x = xlim_range[2] * 0.95, y = ylim_range[2] * 0.95,
           label = sprintf("Exacerbated\nn = %s/%s", q_sig["TR"], q_counts["TR"]),
           hjust = 1, vjust = 1, size = txt_quad, fontface = "bold",
           color = "#DC2626", fill = alpha("white", 0.92),
           label.padding = unit(3, "pt")) +
  annotate("label", x = xlim_range[1] * 0.95, y = ylim_range[1] * 0.95,
           label = sprintf("Exacerbated\nn = %s/%s", q_sig["BL"], q_counts["BL"]),
           hjust = 0, vjust = 0, size = txt_quad, fontface = "bold",
           color = "#DC2626", fill = alpha("white", 0.92),
           label.padding = unit(3, "pt")) +
  annotate("label",
           x = 0,
           y = ylim_range[1] * 0.85,
           label = sprintf("Melov magnitude reversal: %.1f%%, p = %.2f\n(%d aging-sig. proteins, n = 15 subjects)",
                           melov_rev_pct, melov_p, melov_n),
           hjust = 0.5, vjust = 1, size = txt_stat, fontface = "italic",
           color = "grey35", fill = alpha("white", 0.85),
           label.padding = unit(2, "pt")) +
  coord_fixed(ratio = 1, xlim = xlim_range, ylim = ylim_range, expand = FALSE) +
  labs(title = "Protein-Level Reversal of Aging by Training",
       subtitle = sub_txt,
       x = expression(log[2]*FC ~ "(Aging)"),
       y = expression(log[2]*FC ~ "(Training Old)")) +
  FIG_THEME +
  theme(legend.position = "none")

ggsave(file.path(RPT, "panel_C_reversal_scatter.pdf"), pC,
       width = PW, height = PH, units = "mm", device = pdf_device)
ggsave(file.path(RPT, "panel_C_reversal_scatter.png"), pC,
       width = PW, height = PH, units = "mm", dpi = 300)

scatter_df %>%
  transmute(gene,
            logFC_Aging        = round(logFC_Aging, 4),
            logFC_Training_Old = round(logFC_Training_Old, 4),
            pi_score_Aging        = round(pi_Aging, 6),
            pi_score_Training_Old = round(pi_Training_Old, 6),
            pi_score_Reversal     = round(pi_Reversal, 6),
            significance          = as.character(significance),
            quadrant, imputed) %>%
  arrange(significance, desc(abs(logFC_Aging) + abs(logFC_Training_Old))) %>%
  write_csv(file.path(DAT, "panel_C", "reversal_scatter.csv"))

tibble(
  metric   = c("Pearson_r", "Spearman_rho", "Reversal_pct",
               "Pearson_r_sig", "Spearman_rho_sig", "Reversal_pct_sig",
               "Melov_magnitude_reversal_pct"),
  estimate = c(cor_r$estimate, cor_rho$estimate, reversal_pct,
               cor_r_sig$estimate, cor_rho_sig$estimate, rev_sig,
               melov_rev_pct),
  ci_lower = c(cor_r$conf.int[1], rho_ci[1], rev_ci[1],
               cor_r_sig$conf.int[1], rho_sig_ci[1], NA_real_,
               melov_df$reversal_pct_ci_lower),
  ci_upper = c(cor_r$conf.int[2], rho_ci[2], rev_ci[2],
               cor_r_sig$conf.int[2], rho_sig_ci[2], NA_real_,
               melov_df$reversal_pct_ci_upper),
  p_value  = c(cor_r$p.value, cor_rho$p.value, NA_real_,
               cor_r_sig$p.value, cor_rho_sig$p.value, NA_real_,
               melov_p),
  n        = c(n_obs, n_obs, n_obs, n_sig, n_sig, n_sig, melov_n),
  note     = c("95% CI from cor.test()",
               "95% CI via Fisher z-transformation",
               "95% BCa bootstrap CI (10000 replicates, all proteins)",
               "Sig proteins only — 95% CI from cor.test()",
               "Sig proteins only — 95% BCa bootstrap CI (10000 replicates)",
               "Sig proteins only — no CI",
               sprintf("Melov permutation test (%d perms)", melov_df$n_permutations))
) %>%
  write_csv(file.path(DAT, "panel_C", "reversal_scatter_stats.csv"))

cat("F3 Panel C done\n")
