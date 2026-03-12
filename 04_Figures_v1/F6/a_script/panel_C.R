################################################################################
#   Figure 6 — Panel C: Training-Induced Module Changes Track Phenotype
#   Scatter + facet of delta eigengene (Post-Pre) vs delta phenotype
################################################################################
#
# STAT AUDIT — Panel C (2026-02-27)
# ---------------------------------------------------------------------------
# 1. Pearson r: ADDED 95% CI via cor.test()$conf.int.
# 2. Partial r (ppcor): ADDED Fisher z-transform CI (k=1 covariate).
# 3. BH correction: applied across all 6 facets (2 outcomes x 3 modules)
#    for both raw and partial p-values.  PASS.
# 4. REDUNDANCY NOTE: Panel C tests CONCURRENT plasticity (delta ME vs
#    delta phenotype), while Panel B tests PREDICTIVE association (baseline
#    ME vs delta phenotype). These are complementary, not redundant:
#    - B addresses "Can baseline proteome predict who responds?"
#    - C addresses "Do module-level changes co-occur with phenotype gains?"
#    The top-3 modules may differ between B and C (selected independently
#    by max |r| in each context).
# 5. Same power caveats as Panel B (N ~ 15-31).
# ---------------------------------------------------------------------------
#

source("04_Figures/F6/a_script/YvO_F6_setup.R")

# ---- Plasticity correlations: top 3 modules by max |r| ----
plast_cor <- tibble(module = colnames(delta_me)) %>%
  rowwise() %>%
  mutate(
    r_vl  = cor(delta_me[common_subj, module],
                pheno_wide$delta_VL[match(common_subj, pheno_wide$subject_key)],
                use = "complete.obs"),
    r_lbm = cor(delta_me[common_subj, module],
                pheno_wide$delta_LBM[match(common_subj, pheno_wide$subject_key)],
                use = "complete.obs"),
    max_r = max(abs(r_vl), abs(r_lbm), na.rm = TRUE)
  ) %>% ungroup()

top3_plast <- plast_cor %>% arrange(desc(max_r)) %>% head(3) %>% pull(module)

# ---- Build long-format data ----
plast_long <- delta_me[common_subj, , drop = FALSE] %>%
  as.data.frame() %>%
  rownames_to_column("subject_key") %>%
  pivot_longer(cols = all_of(top3_plast), names_to = "module",
               values_to = "delta_ME") %>%
  left_join(pheno_wide %>% dplyr::select(subject_key, delta_VL, delta_LBM),
            by = "subject_key") %>%
  left_join(subj_age, by = "subject_key") %>%
  pivot_longer(cols = c(delta_VL, delta_LBM), names_to = "outcome",
               values_to = "delta_pheno") %>%
  mutate(outcome_label = outcome_labels[outcome])

# Strip ME prefix for readable facet labels
plast_long$module_label <- mod_bio_labels[gsub("^ME", "", plast_long$module)]
# If any NA (unmapped), use the stripped name with title case
plast_long$module_label[is.na(plast_long$module_label)] <-
  str_to_title(gsub("^ME", "", plast_long$module[is.na(plast_long$module_label)]))

# ---- Compute per-facet correlation stats (raw + age-adjusted partial) ----
# Reference: Kim 2015, Commun Korean Stat Soc (ppcor package)
suppressPackageStartupMessages(library(ppcor))

plast_stats <- plast_long %>%
  group_by(module, module_label, outcome_label) %>%
  summarize(
    n_obs = sum(!is.na(delta_ME) & !is.na(delta_pheno)),
    r = cor(delta_ME, delta_pheno, use = "complete.obs"),
    p = tryCatch(cor.test(delta_ME, delta_pheno)$p.value,
                 error = function(e) NA_real_),
    # 95% CI for Pearson r via cor.test()
    r_ci_lo = tryCatch(cor.test(delta_ME, delta_pheno)$conf.int[1],
                       error = function(e) NA_real_),
    r_ci_hi = tryCatch(cor.test(delta_ME, delta_pheno)$conf.int[2],
                       error = function(e) NA_real_),
    r_partial = tryCatch({
      age_numeric <- ifelse(age == "Old", 1, 0)
      pcor.test(delta_ME, delta_pheno, age_numeric, method = "pearson")$estimate
    }, error = function(e) NA_real_),
    p_partial = tryCatch({
      age_numeric <- ifelse(age == "Old", 1, 0)
      pcor.test(delta_ME, delta_pheno, age_numeric, method = "pearson")$p.value
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
pC <- ggplot(plast_long, aes(x = delta_ME, y = delta_pheno)) +
  geom_point(aes(color = age), size = 1.5, alpha = 0.7) +
  geom_smooth(method = "lm", se = TRUE, linewidth = 0.5, color = "grey30") +
  geom_text(data = plast_stats, aes(label = label),
            x = -Inf, y = Inf, hjust = -0.05, vjust = 1.3,
            size = 2.2, fontface = "bold", color = "grey25",
            inherit.aes = FALSE) +
  facet_grid(outcome_label ~ module_label, scales = "free") +
  scale_color_manual(values = AGE_COLORS) +
  labs(title    = "C  Training-Induced Module Changes Track Phenotype",
       subtitle = "Delta eigengene (Post\u2013Pre) vs delta phenotype",
       x = "Delta Module Eigengene", y = NULL,
       color = "Age") +
  THEME_PUB +
  LEGEND_THEME

# ---- Save ----
ggsave(file.path(RPT_DIR, "panel_C_plasticity.pdf"), pC,
       width = 300, height = 250, units = "mm")
ggsave(file.path(RPT_DIR, "panel_C_plasticity.png"), pC,
       width = 300, height = 250, units = "mm", dpi = 300)

write_csv(plast_long, file.path(DAT_DIR, "03_panel_C_plasticity.csv"))

cat("Panel C done\n")
