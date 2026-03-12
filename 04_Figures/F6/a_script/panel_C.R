# Figure 6 — Panel C: Training-Induced Module Changes Track Phenotype
# Delta eigengene (Post-Pre) vs delta phenotype
#
# STAT AUDIT (2026-02-27):
# 1. Pearson r + 95% CI via cor.test().
# 2. Partial r (ppcor) + Fisher z CI (k=1).
# 3. BH correction across all 6 facets. PASS.
# 4. Complementary to Panel B: C = concurrent plasticity (delta ME ~ delta
#    pheno), B = predictive (baseline ME -> delta pheno). Top-3 may differ.
# 5. Same power caveats as Panel B (N ~15-31).

setwd(rprojroot::find_rstudio_root_file())
source("04_Figures/shared/style.R")

suppressPackageStartupMessages({
  library(tidyverse)
  library(ppcor)
})

RPT <- "04_Figures/F6/b_reports"
DAT <- "04_Figures/F6/c_data"

delta_me   <- readRDS(file.path(DAT, "delta_me.rds"))
pheno_wide <- read_csv(file.path(DAT, "pheno_wide.csv"), show_col_types = FALSE)
subj_age   <- read_csv(file.path(DAT, "subj_age.csv"), show_col_types = FALSE)
shared     <- readRDS(file.path(DAT, "shared_objects.rds"))
common_subj    <- shared$common_subj
outcome_labels <- shared$outcome_labels
mod_bio_labels <- shared$mod_bio_labels

pdf_device <- get_pdf_device()

message("Panel C: training-induced module changes track phenotype...")

PC_W <- 300  # panel width mm
PC_H <- 250  # panel height mm

txt_stat <- scale_text(BASE_STAT, PC_W)
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

plast_long$module_label <- mod_bio_labels[gsub("^ME", "", plast_long$module)]
plast_long$module_label[is.na(plast_long$module_label)] <-
  str_to_title(gsub("^ME", "", plast_long$module[is.na(plast_long$module_label)]))

# Per-facet correlations: raw + age-adjusted partial (Kim 2015, ppcor)
plast_stats <- plast_long %>%
  group_by(module, module_label, outcome_label) %>%
  summarize(
    n_obs = sum(!is.na(delta_ME) & !is.na(delta_pheno)),
    r = cor(delta_ME, delta_pheno, use = "complete.obs"),
    p = tryCatch(cor.test(delta_ME, delta_pheno)$p.value,
                 error = function(e) NA_real_),
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

pC <- ggplot(plast_long, aes(x = delta_ME, y = delta_pheno)) +
  geom_point(aes(color = age), size = 1.5, alpha = 0.7) +
  geom_smooth(method = "lm", se = TRUE, linewidth = 0.5, color = "grey30") +
  geom_text(data = plast_stats, aes(label = label),
            x = -Inf, y = Inf, hjust = -0.05, vjust = 1.3,
            size = txt_stat, fontface = "bold", color = "grey25",
            inherit.aes = FALSE) +
  facet_grid(outcome_label ~ module_label, scales = "free") +
  scale_color_manual(values = AGE_COLORS) +
  labs(title    = "Training-Induced Module Changes Track Phenotype",
       subtitle = "Delta eigengene (Post\u2013Pre) vs delta phenotype",
       x = "Delta Module Eigengene", y = NULL,
       color = "Age") +
  FIG_THEME +
  theme(legend.position = "none")

ggsave(file.path(RPT, "panel_C_plasticity.pdf"), pC,
       width = PC_W, height = PC_H, units = "mm",
       device = pdf_device, limitsize = FALSE)
ggsave(file.path(RPT, "panel_C_plasticity.png"), pC,
       width = PC_W, height = PC_H, units = "mm", dpi = 300)

write_csv(plast_long, file.path(DAT, "03_panel_C_plasticity.csv"))

message("  Panel C saved")
