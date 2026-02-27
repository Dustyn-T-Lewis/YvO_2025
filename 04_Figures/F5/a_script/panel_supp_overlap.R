################################################################################
#   YvO Figure 5 — Supplementary: Cross-Method Integration
#   4-axis alluvial/Sankey diagram linking:
#     Concordance (F2) → FCM Cluster (F4) → WGCNA Module (F5) → Reversal (F3)
#   Fisher's exact test for FCM × WGCNA pairs.
#
#   Generates: supp_fcm_wgcna_overlap.pdf
#              + c_data/fig5_fcm_wgcna_overlap.csv (Fisher's test results)
#              + c_data/fig5_cross_method_integration.csv (4-way merged data)
################################################################################
#
# ── STAT AUDIT (2026-02-27) ──────────────────────────────────────────────────
#
# FISHER'S EXACT TEST (FCM x WGCNA)
#   Pairwise 2x2 Fisher's exact test for each FCM cluster × WGCNA module
#   pair.  BH correction applied globally across all pairs (line 104).
#   Fisher's exact is appropriate (some cells may be sparse).
#   Odds ratios and CIs are available from fisher.test() but only the
#   point estimate is stored.  AUDIT ADDITION: store 95% CI on OR.
#
# CHI-SQUARED TESTS
#   Concordance × FCM and WGCNA × Reversal use chi-squared tests.
#   These are global omnibus tests.  Some expected cell counts may be
#   low (< 5) for smaller modules.  AUDIT ADDITION: add note when
#   expected counts < 5 (chi-squared approximation may be unreliable).
#
# AUDIT VERDICT: Fisher's test with BH correction is appropriate.
# CIs on OR added below.
# ─────────────────────────────────────────────────────────────────────────────

source("04_Figures/F5/a_script/YvO_F5_setup.R")

suppressPackageStartupMessages({
  library(ggalluvial)
})

# === 1. LOAD DATA =============================================================

cat("Panel Supp: Cross-method integration...\n")

# FCM assignments from Figure 4 (cluster column already has "C1"-"C4" labels)
fcm_path <- "04_Figures/F4/c_data/mfuzz_assignments.csv"
if (!file.exists(fcm_path)) stop("FCM assignments not found at ", fcm_path,
                                  ". Run Figure 4 pipeline first.")
fcm_df <- read_csv(fcm_path, show_col_types = FALSE) %>%
  dplyr::select(gene, fcm_cluster = cluster, fcm_membership = membership)

# WGCNA module assignments (already loaded via F5 setup as module_df)
wgcna_df <- module_df %>%
  dplyr::select(gene, wgcna_module = module_color) %>%
  filter(wgcna_module != "grey") %>%
  mutate(wgcna_module = str_to_title(wgcna_module))

# Concordance from Figure 2 (ASCII labels for PDF compatibility)
conc_path <- "04_Figures/F2/c_data/panel_C/concordance.csv"
if (!file.exists(conc_path)) stop("Concordance data not found at ", conc_path)
conc_df <- read_csv(conc_path, show_col_types = FALSE) %>%
  dplyr::select(gene, concordance_quadrant = quadrant) %>%
  mutate(concordance_quadrant = case_when(
    concordance_quadrant == "Concordant Up"                      ~ "Conc. Up",
    concordance_quadrant == "Concordant Down"                    ~ "Conc. Down",
    concordance_quadrant == "Discordant (Young Up / Old Down)"   ~ "Disc. YUp/ODn",
    concordance_quadrant == "Discordant (Young Down / Old Up)"   ~ "Disc. YDn/OUp",
    TRUE ~ concordance_quadrant
  ))

# Reversal/rejuvenation from Figure 3
rev_path <- "04_Figures/F3/c_data/fig3_protein_rejuvenation.csv"
if (!file.exists(rev_path)) stop("Rejuvenation data not found at ", rev_path)
rev_df <- read_csv(rev_path, show_col_types = FALSE) %>%
  dplyr::select(gene, reversal_quadrant = quadrant)

# === 2. MERGE ALL FOUR LAYERS ================================================

merged_df <- fcm_df %>%
  inner_join(wgcna_df, by = "gene") %>%
  inner_join(conc_df, by = "gene") %>%
  inner_join(rev_df, by = "gene")

cat(sprintf("  Merged: %d proteins across all 4 analyses\n", nrow(merged_df)))
cat(sprintf("    FCM: %d | WGCNA (non-grey): %d | Concordance: %d | Reversal: %d\n",
            nrow(fcm_df), nrow(wgcna_df), nrow(conc_df), nrow(rev_df)))

# Export merged data
write_csv(merged_df, file.path(DAT_DIR, "10_cross_method_integration.csv"))

# === 3. FISHER'S EXACT TEST: FCM × WGCNA =====================================

fcm_wgcna_df <- merged_df %>%
  dplyr::select(gene, fcm_cluster, wgcna_module) %>%
  distinct()

fcm_clusters  <- sort(unique(fcm_wgcna_df$fcm_cluster))
wgcna_modules <- sort(unique(fcm_wgcna_df$wgcna_module))
n_total <- nrow(fcm_wgcna_df)

fisher_results <- list()
for (cl in fcm_clusters) {
  for (mod in wgcna_modules) {
    in_cl  <- fcm_wgcna_df$fcm_cluster == cl
    in_mod <- fcm_wgcna_df$wgcna_module == mod

    ct <- matrix(c(
      sum(in_cl & in_mod),
      sum(in_cl & !in_mod),
      sum(!in_cl & in_mod),
      sum(!in_cl & !in_mod)
    ), nrow = 2)

    ft <- fisher.test(ct)

    fisher_results <- c(fisher_results, list(tibble(
      fcm_cluster  = cl,
      wgcna_module = mod,
      n_overlap    = ct[1, 1],
      odds_ratio   = round(ft$estimate, 3),
      or_ci_lo     = round(ft$conf.int[1], 3),
      or_ci_hi     = round(ft$conf.int[2], 3),
      p_value      = ft$p.value
    )))
  }
}

fisher_df <- bind_rows(fisher_results)
fisher_df$p_adj_bh <- p.adjust(fisher_df$p_value, method = "BH")
fisher_df$significant <- fisher_df$p_adj_bh < 0.05

write_csv(fisher_df, file.path(DAT_DIR, "09_fcm_wgcna_overlap.csv"))

n_sig <- sum(fisher_df$significant)
cat(sprintf("  Fisher's exact: %d / %d pairs significant (BH < 0.05)\n",
            n_sig, nrow(fisher_df)))

# === 4. ALLUVIAL / SANKEY PLOT ================================================

# Count proteins per unique combination of all 4 axes
alluvial_data <- merged_df %>%
  count(concordance_quadrant, fcm_cluster, wgcna_module, reversal_quadrant,
        name = "n")

# Order factors for visual clarity
alluvial_data$concordance_quadrant <- factor(
  alluvial_data$concordance_quadrant,
  levels = c("Conc. Up", "Conc. Down", "Disc. YUp/ODn", "Disc. YDn/OUp")
)
alluvial_data$fcm_cluster <- factor(alluvial_data$fcm_cluster,
                                     levels = paste0("C", 1:4))
alluvial_data$wgcna_module <- factor(
  alluvial_data$wgcna_module,
  levels = c("Turquoise", "Blue", "Brown", "Yellow", "Green",
             "Red", "Black", "Pink", "Magenta")
)
alluvial_data$reversal_quadrant <- factor(
  alluvial_data$reversal_quadrant,
  levels = c("Rejuvenation", "Reversal",
             "Non-restored (+/+)", "Non-restored (-/-)")
)

# Color palette: FCM cluster (matches Figure 4)
fcm_pal <- c("C1" = "#F8766D", "C2" = "#7CAE00", "C3" = "#00BFC4", "C4" = "#C77CFF")

p_sankey <- ggplot(alluvial_data,
       aes(axis1 = concordance_quadrant,
           axis2 = fcm_cluster,
           axis3 = wgcna_module,
           axis4 = reversal_quadrant,
           y = n)) +
  geom_alluvium(aes(fill = fcm_cluster), alpha = 0.5, width = 1/6,
                curve_type = "sigmoid") +
  geom_stratum(width = 1/6, fill = "grey90", color = "grey40", linewidth = 0.3) +
  geom_text(stat = "stratum", aes(label = after_stat(stratum)),
            size = 2.2, fontface = "bold") +
  scale_x_discrete(limits = c("Training\nConcordance",
                               "FCM\nCluster",
                               "WGCNA\nModule",
                               "Aging\nReversal"),
                   expand = c(0.12, 0.05)) +
  scale_fill_manual(values = fcm_pal, name = "FCM Cluster") +
  labs(
    title = "Cross-Method Integration: Concordance > FCM > WGCNA > Reversal",
    subtitle = sprintf(
      "%d proteins | %d sig. FCM x WGCNA pairs (Fisher BH < 0.05) | ribbons colored by FCM cluster",
      nrow(merged_df), n_sig
    ),
    y = "Number of Proteins"
  ) +
  THEME_PUB +
  theme(legend.position = "bottom",
        plot.title = element_text(size = 9, face = "bold"),
        plot.subtitle = element_text(size = 6, color = "grey30", face = "italic"),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.title.y = element_text(size = 7),
        plot.margin = margin(5, 15, 5, 15))

# Save using quartz PDF device (macOS native; supports colors + Unicode)
outpath <- file.path(RPT_DIR, "supp_fcm_wgcna_overlap.pdf")
quartz(type = "pdf", file = outpath, width = 280 / 25.4, height = 300 / 25.4)
print(p_sankey)
dev.off()

# === 5. SUMMARY STATISTICS ====================================================

# Cross-tabulation summaries for manuscript
cat("\n  === Cross-tabulation summaries ===\n")

# Concordance x FCM
cat("\n  Concordance > FCM:\n")
ct_conc_fcm <- table(merged_df$concordance_quadrant, merged_df$fcm_cluster)
print(ct_conc_fcm)
chi_cf <- chisq.test(ct_conc_fcm)
cat(sprintf("  Chi-sq = %.1f, df = %d, p = %.2e\n",
            chi_cf$statistic, chi_cf$parameter, chi_cf$p.value))
# STAT AUDIT: check expected counts
if (any(chi_cf$expected < 5))
  cat(sprintf("  WARNING: %d cells with expected count < 5 (chi-sq approximation may be unreliable)\n",
              sum(chi_cf$expected < 5)))

# FCM x WGCNA (already done via Fisher's)
cat(sprintf("\n  FCM x WGCNA: %d significant pairs\n", n_sig))

# WGCNA x Reversal
cat("\n  WGCNA > Reversal:\n")
ct_wgcna_rev <- table(merged_df$wgcna_module, merged_df$reversal_quadrant)
print(ct_wgcna_rev)
chi_wr <- chisq.test(ct_wgcna_rev)
cat(sprintf("  Chi-sq = %.1f, df = %d, p = %.2e\n",
            chi_wr$statistic, chi_wr$parameter, chi_wr$p.value))
# STAT AUDIT: check expected counts
if (any(chi_wr$expected < 5))
  cat(sprintf("  WARNING: %d cells with expected count < 5 (chi-sq approximation may be unreliable)\n",
              sum(chi_wr$expected < 5)))

# Full 4-way
cat(sprintf("\n  Total unique combinations: %d\n", nrow(alluvial_data)))
cat(sprintf("  Combinations with n >= 10: %d\n", sum(alluvial_data$n >= 10)))

cat("\n  Saved: supp_fcm_wgcna_overlap.pdf\n")
cat("         09_fcm_wgcna_overlap.csv\n")
cat("         10_cross_method_integration.csv\n")
cat("Panel Supp (cross-method integration) complete\n")
