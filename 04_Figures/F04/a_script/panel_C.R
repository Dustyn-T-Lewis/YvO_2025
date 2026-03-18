# F2 Panel C: Age x Training Interaction Volcano Ring
setwd(rprojroot::find_rstudio_root_file())
source("04_Figures/shared/style.R")
source("04_Figures/shared/volcano_ring.R")

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
})

PC_W <- 190
PH   <- 180

RPT <- "04_Figures/F04/b_reports"
DAT <- "04_Figures/F04/c_data"
dir.create(RPT, recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(DAT, "panel_C"), recursive = TRUE, showWarnings = FALSE)

dep_df <- read_csv("03_DEP/c_data/03_combined_results.csv", show_col_types = FALSE)

fgsea_cache <- file.path(DAT, "shared", "fgsea_tstat_all_v2.csv")
if (!file.exists(fgsea_cache)) {
  f1_cache <- "04_Figures/shared/fgsea_f1_panel_F.csv"
  if (file.exists(f1_cache)) {
    dir.create(file.path(DAT, "shared"), recursive = TRUE, showWarnings = FALSE)
    file.copy(f1_cache, fgsea_cache)
  } else stop("fGSEA cache not found")
}
fgsea_all <- read_csv(fgsea_cache, show_col_types = FALSE)

top_terms_C <- select_ring_terms(fgsea_all, "Interaction", n_each = 6)
ring_C      <- build_ring_with_gaps(top_terms_C, "Interaction", fgsea_all, n_each = 6)

pC <- make_volcano_ring(
  de_df = dep_df, go_df = fgsea_all, contrast = "Interaction",
  title = NULL,
  contrast_title    = "Age \u00d7 Training Interaction",
  contrast_subtitle = "Training_Old \u2212 Training_Young",
  ring_data_override = ring_C,
  label_size = scale_text(BASE_PATHWAY, PC_W),
  title_size = scale_text(BASE_TAG, PC_W),
  point_size = 1.2, point_alpha = 0.55,
  count_label_size = scale_text(BASE_COUNT, PC_W)
)

ggsave(file.path(RPT, "panel_C_volcano.pdf"), pC,
       width = PC_W, height = PH, units = "mm", device = get_pdf_device())
ggsave(file.path(RPT, "panel_C_volcano.png"), pC,
       width = PC_W, height = PH, units = "mm", dpi = 300)

ring_data_C <- attr(pC, "ring_data")
if (!is.null(ring_data_C) && nrow(ring_data_C) > 0) {
  write_csv(ring_data_C %>% dplyr::select(-gene_list),
            file.path(DAT, "panel_C", "ring_terms.csv"))
}

dep_df %>%
  transmute(
    gene,
    log2_fold_change = round(logFC_Interaction, 4),
    neg_log10_pvalue = round(-log10(P.Value_Interaction), 4),
    pi_score         = round(pi_score_Interaction, 6),
    adjusted_pvalue  = round(adj.P.Val_Interaction, 6),
    direction = case_when(
      pi_score_Interaction < 0.05 & logFC_Interaction > 0 ~ "Up",
      pi_score_Interaction < 0.05 & logFC_Interaction < 0 ~ "Down",
      TRUE ~ "NS"
    )
  ) %>%
  filter(!is.na(log2_fold_change), !is.na(neg_log10_pvalue)) %>%
  arrange(pi_score) %>%
  write_csv(file.path(DAT, "panel_C", "volcano_interaction.csv"))

message("F2 Panel C done")
