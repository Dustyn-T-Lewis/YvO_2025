# Figure 6 — Panel B: Baseline Eigengene vs Training Response
#
# STAT AUDIT (2026-02-27):
# 1. Pearson r + 95% CI via cor.test().
# 2. Partial r (ppcor) + Fisher z CI (Bonett & Wright 2000, k=1).
# 3. BH correction across all 6 facets (2 outcomes x 3 modules). PASS.
# 4. POWER: N ~15-31, power ~0.46-0.86 for r=0.5. Nulls need caution.
# 5. Complementary to Panel C: B = predictive (baseline ME -> delta pheno),
#    C = concurrent plasticity (delta ME ~ delta pheno).

setwd(rprojroot::find_rstudio_root_file())
source("04_Figures/shared/style.R")

suppressPackageStartupMessages({
  library(tidyverse)
  library(ppcor)
})

RPT <- "04_Figures/F6/b_reports"
DAT <- "04_Figures/F6/c_data"

me_pre    <- readRDS(file.path(DAT, "me_pre.rds"))
pheno_wide <- read_csv(file.path(DAT, "pheno_wide.csv"), show_col_types = FALSE)
subj_age  <- read_csv(file.path(DAT, "subj_age.csv"), show_col_types = FALSE)
shared    <- readRDS(file.path(DAT, "shared_objects.rds"))
top3           <- shared$top3
common_subj    <- shared$common_subj
outcome_labels <- shared$outcome_labels
mod_bio_labels <- shared$mod_bio_labels

pdf_device <- get_pdf_device()

message("Panel B: baseline eigengene vs training response...")

PB_W <- 300  # panel width mm
PB_H <- 250  # panel height mm

txt_stat <- scale_text(BASE_STAT, PB_W)
txt_axis <- scale_text(BASE_STAT, PB_W)
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

pred_long$module_label <- mod_bio_labels[gsub("^ME", "", pred_long$module)]
pred_long$module_label[is.na(pred_long$module_label)] <-
  str_to_title(gsub("^ME", "", pred_long$module[is.na(pred_long$module_label)]))

# Per-facet correlations: raw + age-adjusted partial (Kim 2015, ppcor)
pred_stats <- pred_long %>%
  group_by(module, module_label, outcome_label) %>%
  summarize(
    n_obs = sum(!is.na(baseline_ME) & !is.na(delta_pheno)),
    r = cor(baseline_ME, delta_pheno, use = "complete.obs"),
    p = tryCatch(cor.test(baseline_ME, delta_pheno)$p.value,
                 error = function(e) NA_real_),
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
  rowwise() %>%
  mutate(rp_ci_lo = fisher_z_ci(r_partial, n_obs, k = 1)[["lo"]],
         rp_ci_hi = fisher_z_ci(r_partial, n_obs, k = 1)[["hi"]]) %>%
  ungroup() %>%
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

pB <- ggplot(pred_long, aes(x = baseline_ME, y = delta_pheno)) +
  geom_point(aes(color = age), size = 1.5, alpha = 0.7) +
  geom_smooth(method = "lm", se = TRUE, linewidth = 0.5, color = "grey30") +
  geom_text(data = pred_stats, aes(label = label),
            x = -Inf, y = Inf, hjust = -0.05, vjust = 1.3,
            size = txt_stat, fontface = "bold", color = "grey25",
            inherit.aes = FALSE) +
  facet_grid(outcome_label ~ module_label, scales = "free") +
  scale_color_manual(values = AGE_COLORS) +
  labs(title    = "Baseline Eigengene vs Training Response",
       subtitle = "Baseline module eigengene (Pre) vs change in phenotype | No facet reaches p_adj < 0.05",
       x = "Baseline Module Eigengene", y = NULL,
       color = "Age") +
  FIG_THEME +
  theme(legend.position = "none")

ggsave(file.path(RPT, "panel_B_baseline_association.pdf"), pB,
       width = PB_W, height = PB_H, units = "mm",
       device = pdf_device, limitsize = FALSE)
ggsave(file.path(RPT, "panel_B_baseline_association.png"), pB,
       width = PB_W, height = PB_H, units = "mm", dpi = 300)

write_csv(pred_long, file.path(DAT, "02_panel_B_baseline_association.csv"))

message("  Panel B saved")
