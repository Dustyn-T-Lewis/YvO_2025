# F2 Panel B: Training Response (Old) Volcano Ring
setwd(rprojroot::find_rstudio_root_file())
source("04_Figures_v2/shared/style.R")
source("04_Figures_v2/shared/volcano_ring.R")

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
})

PB_W <- 190
PH   <- 180

RPT <- "04_Figures_v2/F2/b_reports"
DAT <- "04_Figures_v2/F2/c_data"
dir.create(RPT, recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(DAT, "panel_B"), recursive = TRUE, showWarnings = FALSE)

dep_df <- read_csv("03_DEP/c_data/03_combined_results.csv", show_col_types = FALSE)

fgsea_cache <- file.path(DAT, "shared", "fgsea_tstat_all_v2.csv")
if (!file.exists(fgsea_cache)) {
  f1_cache <- "04_Figures_v2/F1/c_data/06_panel_F_fgsea_results.csv"
  if (file.exists(f1_cache)) {
    dir.create(file.path(DAT, "shared"), recursive = TRUE, showWarnings = FALSE)
    file.copy(f1_cache, fgsea_cache)
  } else stop("fGSEA cache not found")
}
fgsea_all <- read_csv(fgsea_cache, show_col_types = FALSE)

top_terms_B <- select_ring_terms(fgsea_all, "Training_Old", n_each = 6)
ring_B      <- build_ring_with_gaps(top_terms_B, "Training_Old", fgsea_all, n_each = 6)

pB <- make_volcano_ring(
  de_df = dep_df, go_df = fgsea_all, contrast = "Training_Old",
  title = NULL,
  contrast_title    = "Training Response (Old)",
  contrast_subtitle = "Old_Post \u2212 Old_Pre",
  ring_data_override = ring_B,
  label_size = scale_text(BASE_PATHWAY, PB_W),
  title_size = scale_text(BASE_TAG, PB_W),
  point_size = 1.2, point_alpha = 0.55,
  count_label_size = scale_text(BASE_COUNT, PB_W)
)

ggsave(file.path(RPT, "panel_B_volcano.pdf"), pB,
       width = PB_W, height = PH, units = "mm", device = get_pdf_device())
ggsave(file.path(RPT, "panel_B_volcano.png"), pB,
       width = PB_W, height = PH, units = "mm", dpi = 300)

ring_data_B <- attr(pB, "ring_data")
if (!is.null(ring_data_B) && nrow(ring_data_B) > 0) {
  write_csv(ring_data_B %>% dplyr::select(-gene_list),
            file.path(DAT, "panel_B", "ring_terms.csv"))
}

dep_df %>%
  transmute(
    gene,
    log2_fold_change = round(logFC_Training_Old, 4),
    neg_log10_pvalue = round(-log10(P.Value_Training_Old), 4),
    pi_score         = round(pi_score_Training_Old, 6),
    adjusted_pvalue  = round(adj.P.Val_Training_Old, 6),
    direction = case_when(
      pi_score_Training_Old < 0.05 & logFC_Training_Old > 0 ~ "Up",
      pi_score_Training_Old < 0.05 & logFC_Training_Old < 0 ~ "Down",
      TRUE ~ "NS"
    )
  ) %>%
  filter(!is.na(log2_fold_change), !is.na(neg_log10_pvalue)) %>%
  arrange(pi_score) %>%
  write_csv(file.path(DAT, "panel_B", "volcano_old.csv"))

message("F2 Panel B done")
