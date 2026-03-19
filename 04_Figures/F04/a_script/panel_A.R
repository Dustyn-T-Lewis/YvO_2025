# F4 Panel A: Training Response (Young) Volcano Ring
setwd(rprojroot::find_rstudio_root_file())
source("04_Figures/shared/style.R")
source("04_Figures/shared/volcano_ring.R")

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
})

PA_W <- 190
PH   <- 180

RPT <- "04_Figures/F04/b_reports"
DAT <- "04_Figures/F04/c_data"
dir.create(RPT, recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(DAT, "panel_A"), recursive = TRUE, showWarnings = FALSE)

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

top_terms_A <- select_ring_terms(fgsea_all, "Training_Young")
ring_A      <- build_ring_with_gaps(top_terms_A, "Training_Young", fgsea_all)

pA <- make_volcano_ring(
  de_df = dep_df, go_df = fgsea_all, contrast = "Training_Young",
  title = NULL,
  contrast_title    = "Training Response (Young)",
  contrast_subtitle = "Young_Post \u2212 Young_Pre",
  ring_data_override = ring_A,
  label_size = scale_text(BASE_PATHWAY, PA_W),
  title_size = scale_text(BASE_TAG, PA_W),
  point_size = 1.2, point_alpha = 0.55,
  count_label_size = scale_text(BASE_COUNT, PA_W)
)

ggsave(file.path(RPT, "panel_A_volcano.pdf"), pA,
       width = PA_W, height = PH, units = "mm", device = get_pdf_device())
ggsave(file.path(RPT, "panel_A_volcano.png"), pA,
       width = PA_W, height = PH, units = "mm", dpi = 300)

ring_data_A <- attr(pA, "ring_data")
if (!is.null(ring_data_A) && nrow(ring_data_A) > 0) {
  write_csv(ring_data_A %>% dplyr::select(-gene_list),
            file.path(DAT, "panel_A", "ring_terms.csv"))
}

dep_df %>%
  transmute(
    gene,
    log2_fold_change = round(logFC_Training_Young, 4),
    neg_log10_pvalue = round(-log10(P.Value_Training_Young), 4),
    pi_score         = round(pi_score_Training_Young, 6),
    adjusted_pvalue  = round(adj.P.Val_Training_Young, 6),
    direction = case_when(
      pi_score_Training_Young < 0.05 & logFC_Training_Young > 0 ~ "Up",
      pi_score_Training_Young < 0.05 & logFC_Training_Young < 0 ~ "Down",
      TRUE ~ "NS"
    )
  ) %>%
  filter(!is.na(log2_fold_change), !is.na(neg_log10_pvalue)) %>%
  arrange(pi_score) %>%
  write_csv(file.path(DAT, "panel_A", "volcano_young.csv"))

message("F4 Panel A done")
