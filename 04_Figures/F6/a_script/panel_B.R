################################################################################
#   Figure 6 — Panel B: Baseline Eigengene Associates with Training Response
#   Scatter + facet of baseline module eigengene vs change in phenotype
################################################################################

source("04_Figures/F6/a_script/YvO_F6_setup.R")

# ---- Build long-format data for top 3 predictive modules ----
pred_long <- me_pre[common_subj, , drop = FALSE] %>%
  as.data.frame() %>%
  rownames_to_column("subject_key") %>%
  pivot_longer(cols = all_of(top3), names_to = "module",
               values_to = "baseline_ME") %>%
  left_join(pheno_wide %>% dplyr::select(subject_key, delta_VL, delta_LBM),
            by = "subject_key") %>%
  left_join(subj_age, by = "subject_key") %>%
  pivot_longer(cols = c(delta_VL, delta_LBM), names_to = "outcome",
               values_to = "delta_pheno") %>%
  mutate(
    outcome_label = outcome_labels[outcome],
    mod_color     = gsub("^ME", "", module)
  )

# Strip ME prefix for readable facet labels
pred_long$module_label <- mod_bio_labels[gsub("^ME", "", pred_long$module)]
# If any NA (unmapped), use the stripped name with title case
pred_long$module_label[is.na(pred_long$module_label)] <-
  str_to_title(gsub("^ME", "", pred_long$module[is.na(pred_long$module_label)]))

# ---- Compute per-facet correlation stats (raw + age-adjusted partial) ----
# Reference: Kim 2015, Commun Korean Stat Soc (ppcor package)
suppressPackageStartupMessages(library(ppcor))

pred_stats <- pred_long %>%
  group_by(module, module_label, outcome_label) %>%
  summarize(
    r = cor(baseline_ME, delta_pheno, use = "complete.obs"),
    p = tryCatch(cor.test(baseline_ME, delta_pheno)$p.value,
                 error = function(e) NA_real_),
    r_partial = tryCatch({
      age_numeric <- ifelse(age == "Old", 1, 0)
      pcor.test(baseline_ME, delta_pheno, age_numeric, method = "pearson")$estimate
    }, error = function(e) NA_real_),
    p_partial = tryCatch({
      age_numeric <- ifelse(age == "Old", 1, 0)
      pcor.test(baseline_ME, delta_pheno, age_numeric, method = "pearson")$p.value
    }, error = function(e) NA_real_),
    .groups = "drop"
  ) %>%
  mutate(label = sprintf("r = %.2f (p = %s)\nr_partial = %.2f (p = %s)",
                         r,
                         ifelse(p < 0.001, formatC(p, format = "e", digits = 1),
                                sprintf("%.3f", p)),
                         r_partial,
                         ifelse(is.na(p_partial), "NA",
                                ifelse(p_partial < 0.001,
                                       formatC(p_partial, format = "e", digits = 1),
                                       sprintf("%.3f", p_partial)))))

# ---- Plot ----
pB <- ggplot(pred_long, aes(x = baseline_ME, y = delta_pheno)) +
  geom_point(aes(color = age), size = 1.5, alpha = 0.7) +
  geom_smooth(method = "lm", se = TRUE, linewidth = 0.5, color = "grey30") +
  geom_text(data = pred_stats, aes(label = label),
            x = -Inf, y = Inf, hjust = -0.05, vjust = 1.3,
            size = 2.8, fontface = "bold", color = "grey25",
            inherit.aes = FALSE) +
  facet_grid(outcome_label ~ module_label, scales = "free") +
  scale_color_manual(values = AGE_COLORS) +
  labs(title    = "B  Baseline Eigengene Associates with Training Response",
       subtitle = "Baseline module eigengene (Pre) vs change in phenotype",
       x = "Baseline Module Eigengene", y = NULL,
       color = "Age") +
  THEME_PUB +
  LEGEND_THEME

# ---- Save ----
ggsave(file.path(RPT_DIR, "panel_B_baseline_prediction.pdf"), pB,
       width = 300, height = 250, units = "mm")
ggsave(file.path(RPT_DIR, "panel_B_baseline_prediction.png"), pB,
       width = 300, height = 250, units = "mm", dpi = 300)

write_csv(pred_long, file.path(DAT_DIR, "fig6_panel_B_baseline_prediction.csv"))

cat("Panel B done\n")
