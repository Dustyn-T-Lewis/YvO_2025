################################################################################
#   YvO Figure 5 — Supplementary: FCM-WGCNA Module Overlap
#   Sankey/alluvial diagram showing correspondence between FCM temporal
#   clusters (from Figure 4) and WGCNA co-expression modules (Figure 5).
#   Fisher's exact test for each cluster x module pair.
#
#   Reference: Broadbent et al. 2017 (PMID 28104750)
################################################################################

source("04_Figures/F5/a_script/YvO_F5_setup.R")

suppressPackageStartupMessages({
  library(ggalluvial)
})

# === 1. LOAD DATA ============================================================

cat("Panel Supp: FCM-WGCNA overlap...\n")

# FCM assignments from Figure 4
fcm_path <- "04_Figures/F4/c_data/mfuzz_assignments.csv"
if (!file.exists(fcm_path)) stop("FCM assignments not found at ", fcm_path,
                                  ". Run Figure 4 pipeline first.")
fcm_df <- read_csv(fcm_path, show_col_types = FALSE) %>%
  dplyr::select(gene, fcm_cluster = cluster, fcm_membership = membership) %>%
  mutate(fcm_cluster = as.character(fcm_cluster))

# WGCNA module assignments (already loaded via F5 setup as module_df)
wgcna_df <- module_df %>%
  dplyr::select(gene, wgcna_module = module_color) %>%
  filter(wgcna_module != "grey")

# Merge by gene
overlap_df <- inner_join(fcm_df, wgcna_df, by = "gene")
cat(sprintf("  Overlapping proteins: %d (FCM: %d, WGCNA: %d)\n",
            nrow(overlap_df), nrow(fcm_df), nrow(wgcna_df)))

# === 2. FISHER'S EXACT TEST FOR EACH PAIR ====================================

fcm_clusters  <- sort(unique(overlap_df$fcm_cluster))
wgcna_modules <- sort(unique(overlap_df$wgcna_module))
n_total <- nrow(overlap_df)

fisher_results <- list()
for (cl in fcm_clusters) {
  for (mod in wgcna_modules) {
    in_cl  <- overlap_df$fcm_cluster == cl
    in_mod <- overlap_df$wgcna_module == mod

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
      p_value      = ft$p.value
    )))
  }
}

fisher_df <- bind_rows(fisher_results)
fisher_df$p_adj_bh <- p.adjust(fisher_df$p_value, method = "BH")
fisher_df$significant <- fisher_df$p_adj_bh < 0.05

write_csv(fisher_df, file.path(DAT_DIR, "fig5_fcm_wgcna_overlap.csv"))

n_sig <- sum(fisher_df$significant)
cat(sprintf("  Significant overlaps (BH < 0.05): %d / %d pairs\n",
            n_sig, nrow(fisher_df)))

# === 3. ALLUVIAL / SANKEY PLOT ================================================

# Count proteins per FCM cluster x WGCNA module pair
alluvial_data <- overlap_df %>%
  count(fcm_cluster, wgcna_module, name = "n") %>%
  left_join(fisher_df %>% dplyr::select(fcm_cluster, wgcna_module, p_adj_bh, significant),
            by = c("fcm_cluster", "wgcna_module")) %>%
  filter(n >= 2)  # only show pairs with at least 2 proteins

# Color by FCM cluster
fcm_pal <- scales::hue_pal()(length(fcm_clusters))
names(fcm_pal) <- fcm_clusters

p_sankey <- ggplot(alluvial_data,
       aes(axis1 = fcm_cluster, axis2 = wgcna_module, y = n)) +
  geom_alluvium(aes(fill = fcm_cluster), alpha = 0.6, width = 1/8) +
  geom_stratum(width = 1/8, fill = "grey85", color = "grey40", linewidth = 0.3) +
  geom_text(stat = "stratum", aes(label = after_stat(stratum)),
            size = 2.5, fontface = "bold") +
  scale_x_discrete(limits = c("FCM Cluster", "WGCNA Module"),
                   expand = c(0.1, 0.05)) +
  scale_fill_manual(values = fcm_pal, name = "FCM Cluster") +
  labs(title = "FCM Temporal Cluster — WGCNA Module Overlap",
       subtitle = sprintf("%d proteins | ribbons show shared membership | %d significant pairs (BH < 0.05)",
                          nrow(overlap_df), n_sig),
       y = "Number of Proteins") +
  THEME_PUB +
  theme(legend.position = "bottom",
        plot.title = element_text(size = 10, face = "bold"),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank())

ggsave(file.path(RPT_DIR, "supp_fcm_wgcna_overlap.pdf"), p_sankey,
       width = 220, height = 260, units = "mm", device = pdf)

cat("  Saved: supp_fcm_wgcna_overlap.pdf, fig5_fcm_wgcna_overlap.csv\n")
cat("Panel Supp (FCM-WGCNA overlap) complete\n")
