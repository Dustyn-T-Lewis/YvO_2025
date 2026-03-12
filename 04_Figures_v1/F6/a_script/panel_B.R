################################################################################
#   Figure 6 — Panel B: Baseline Eigengene Associates with Training Response
#   Scatter + facet of baseline module eigengene vs change in phenotype
#
#   INTEGRITY NOTE: This panel is a transparent null finding. After BH
#   correction across all 6 facets (2 outcomes x 3 modules), 0/6 reach
#   p_adj < 0.05. Retained for completeness: the absence of significant
#   baseline-predictive associations is itself informative (modules track
#   concurrent change [Panel C] but do not predict it a priori).
################################################################################
#
# STAT AUDIT — Panel B (2026-02-27)
# ---------------------------------------------------------------------------
# 1. Pearson r: ADDED 95% CI via cor.test()$conf.int.
# 2. Partial r (ppcor): pcor.test() does not return CIs natively.
#    ADDED Fisher z-transform CI (fisher_z_ci from setup, k=1 covariate).
#    Ref: Bonett & Wright 2000, Psychometrika 65:23-28.
# 3. BH correction: applied across all 6 facets (2 outcomes x 3 modules)
#    for both raw and partial p-values.  PASS.
# 4. POWER NOTE: With N ~ length(common_subj) subjects, power to detect
#    r = 0.5 at alpha = 0.05 (two-sided) is approximately 0.46 (N=15) to
#    0.86 (N=31).  Null results should be interpreted cautiously.
# 5. REDUNDANCY WITH PANEL C: Panel B tests whether BASELINE module
#    eigengene predicts future phenotype change (predictive association).
#    Panel C tests whether CHANGE in module eigengene co-occurs with
#    CHANGE in phenotype (concurrent plasticity). These are complementary
#    hypotheses: B = "Can we predict response?", C = "Do modules and
#    phenotypes change together?" They may share variance but test
#    mechanistically distinct questions.
# ---------------------------------------------------------------------------
#

source("04_Figures/F6/a_script/YvO_F6_setup.R")

# ---- Build long-format data for top 3 associative modules ----
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
    n_obs = sum(!is.na(baseline_ME) & !is.na(delta_pheno)),
    r = cor(baseline_ME, delta_pheno, use = "complete.obs"),
    p = tryCatch(cor.test(baseline_ME, delta_pheno)$p.value,
                 error = function(e) NA_real_),
    # 95% CI for Pearson r via cor.test()
    r_ci_lo = tryCatch(cor.test(baseline_ME, delta_pheno)$conf.int[1],
                       error = function(e) NA_real_),
    r_ci_hi = tryCatch(cor.test(baseline_ME, delta_pheno)$conf.int[2],
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
  # 95% CI for partial r via Fisher z-transform (k=1 covariate for age)
  rowwise() %>%
  mutate(rp_ci_lo = fisher_z_ci(r_partial, n_obs, k = 1)[["lo"]],
         rp_ci_hi = fisher_z_ci(r_partial, n_obs, k = 1)[["hi"]]) %>%
  ungroup() %>%
  # BH correction across all facets for both raw and partial p-values
  mutate(p_adj         = p.adjust(p, method = "BH"),
         p_partial_adj = p.adjust(p_partial, method = "BH")) %>%
  mutate(
    partial_stars = case_when(
      is.na(p_partial_adj) ~ "",
      p_partial_adj < 0.001 ~ "***",
      p_partial_adj < 0.01  ~ "**",
      p_partial_adj < 0.05  ~ "*",
      TRUE ~ ""
    ),
    label = sprintf("r = %.2f, r_partial = %.2f%s",
                    r, r_partial, partial_stars)
  )

# ---- Plot ----
pB <- ggplot(pred_long, aes(x = baseline_ME, y = delta_pheno)) +
  geom_point(aes(color = age), size = 1.5, alpha = 0.7) +
  geom_smooth(method = "lm", se = TRUE, linewidth = 0.5, color = "grey30") +
  geom_text(data = pred_stats, aes(label = label),
            x = -Inf, y = Inf, hjust = -0.05, vjust = 1.3,
            size = 2.2, fontface = "bold", color = "grey25",
            inherit.aes = FALSE) +
  facet_grid(outcome_label ~ module_label, scales = "free") +
  scale_color_manual(values = AGE_COLORS) +
  labs(title    = "B  Baseline Eigengene vs Training Response",
       subtitle = "Baseline module eigengene (Pre) vs change in phenotype | No facet reaches p_adj < 0.05",
       x = "Baseline Module Eigengene", y = NULL,
       color = "Age") +
  THEME_PUB +
  LEGEND_THEME

# ---- Save ----
ggsave(file.path(RPT_DIR, "panel_B_baseline_association.pdf"), pB,
       width = 300, height = 250, units = "mm")
ggsave(file.path(RPT_DIR, "panel_B_baseline_association.png"), pB,
       width = 300, height = 250, units = "mm", dpi = 300)

write_csv(pred_long, file.path(DAT_DIR, "02_panel_B_baseline_association.csv"))

cat("Panel B done\n")
