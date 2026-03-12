################################################################################
#   Figure 5 — Supplementary: Cross-Method Integration
#   Self-contained panel script for 04_Figures
#   4-axis alluvial/Sankey: Concordance > FCM > WGCNA > Reversal
#   Generates: supp_fcm_wgcna_overlap.pdf, c_data/09_*.csv, c_data/10_*.csv
################################################################################

setwd(rprojroot::find_rstudio_root_file())
source("04_Figures/shared/style.R")

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(stringr)
  library(ggalluvial)
})

RPT <- "04_Figures/F5/b_reports"
DAT <- "04_Figures/F5/c_data"

pdf_device <- get_pdf_device()

message("Panel Supp: cross-method integration...")

PS_W <- 280  # panel width mm
PS_H <- 300  # panel height mm

# --- Text sizes ---
txt_stratum <- scale_text(BASE_GENE, PS_W)
txt_axis    <- scale_text(BASE_STAT, PS_W)

module_df <- read_csv(file.path(DAT, "wgcna/wgcna_module_assignments.csv"),
                      show_col_types = FALSE)

fcm_path <- "04_Figures/F4/c_data/06_mfuzz_assignments.csv"
if (!file.exists(fcm_path)) {
  fcm_path <- "04_Figures_v1/F4/c_data/06_mfuzz_assignments.csv"
}
if (!file.exists(fcm_path)) stop("FCM assignments not found. Run Figure 4 first.")
fcm_df <- read_csv(fcm_path, show_col_types = FALSE) %>%
  dplyr::select(gene, fcm_cluster = cluster, fcm_membership = membership)

wgcna_df <- module_df %>%
  dplyr::select(gene, wgcna_module = module_color) %>%
  filter(wgcna_module != "grey") %>%
  mutate(wgcna_module = str_to_title(wgcna_module))

conc_path <- "04_Figures/F2/c_data/panel_C/concordance.csv"
if (!file.exists(conc_path)) {
  conc_path <- "04_Figures_v1/F2/c_data/panel_C/concordance.csv"
}
if (!file.exists(conc_path)) stop("Concordance data not found.")
conc_df <- read_csv(conc_path, show_col_types = FALSE) %>%
  dplyr::select(gene, concordance_quadrant = quadrant) %>%
  mutate(concordance_quadrant = case_when(
    concordance_quadrant == "Concordant Up"                      ~ "Conc. Up",
    concordance_quadrant == "Concordant Down"                    ~ "Conc. Down",
    concordance_quadrant == "Discordant (Young Up / Old Down)"   ~ "Disc. YUp/ODn",
    concordance_quadrant == "Discordant (Young Down / Old Up)"   ~ "Disc. YDn/OUp",
    TRUE ~ concordance_quadrant
  ))

rev_path <- "04_Figures/F3/c_data/panel_C/reversal_scatter.csv"
if (!file.exists(rev_path)) {
  rev_path <- "04_Figures_v1/F3/c_data/panel_C/reversal_scatter.csv"
}
if (!file.exists(rev_path)) stop("Reversal data not found.")
rev_df <- read_csv(rev_path, show_col_types = FALSE) %>%
  dplyr::select(gene, reversal_quadrant = quadrant)

merged_df <- fcm_df %>%
  inner_join(wgcna_df, by = "gene") %>%
  inner_join(conc_df, by = "gene") %>%
  inner_join(rev_df, by = "gene")

cat(sprintf("  Merged: %d proteins across all 4 analyses\n", nrow(merged_df)))

write_csv(merged_df, file.path(DAT, "10_cross_method_integration.csv"))

# Fisher's exact test: FCM x WGCNA overlap
fcm_wgcna_df <- merged_df %>%
  dplyr::select(gene, fcm_cluster, wgcna_module) %>%
  distinct()

fcm_clusters  <- sort(unique(fcm_wgcna_df$fcm_cluster))
wgcna_modules <- sort(unique(fcm_wgcna_df$wgcna_module))

fisher_results <- list()
for (cl in fcm_clusters) {
  for (mod in wgcna_modules) {
    in_cl  <- fcm_wgcna_df$fcm_cluster == cl
    in_mod <- fcm_wgcna_df$wgcna_module == mod

    ct <- matrix(c(
      sum(in_cl & in_mod), sum(in_cl & !in_mod),
      sum(!in_cl & in_mod), sum(!in_cl & !in_mod)
    ), nrow = 2)

    ft <- fisher.test(ct)

    fisher_results <- c(fisher_results, list(tibble(
      fcm_cluster = cl, wgcna_module = mod,
      n_overlap = ct[1, 1],
      odds_ratio = round(ft$estimate, 3),
      or_ci_lo = round(ft$conf.int[1], 3),
      or_ci_hi = round(ft$conf.int[2], 3),
      p_value = ft$p.value
    )))
  }
}

fisher_df <- bind_rows(fisher_results)
fisher_df$p_adj_bh <- p.adjust(fisher_df$p_value, method = "BH")
fisher_df$significant <- fisher_df$p_adj_bh < 0.05

write_csv(fisher_df, file.path(DAT, "09_fcm_wgcna_overlap.csv"))

n_sig <- sum(fisher_df$significant)

alluvial_data <- merged_df %>%
  count(concordance_quadrant, fcm_cluster, wgcna_module, reversal_quadrant,
        name = "n")

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
            size = txt_stratum, fontface = "bold") +
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
  FIG_THEME +
  theme(legend.position = "none",
        axis.text.y    = element_blank(),
        axis.ticks.y   = element_blank(),
        axis.title.y   = element_text(size = txt_axis),
        plot.margin    = margin(5, 15, 5, 15))

ggsave(file.path(RPT, "supp_fcm_wgcna_overlap.pdf"), p_sankey,
       width = PS_W, height = PS_H, units = "mm",
       device = pdf_device, limitsize = FALSE)
ggsave(file.path(RPT, "supp_fcm_wgcna_overlap.png"), p_sankey,
       width = PS_W, height = PS_H, units = "mm",
       dpi = 300, limitsize = FALSE)

message("  Panel Supp (cross-method integration) complete")
