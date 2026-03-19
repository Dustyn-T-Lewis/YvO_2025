# F5 Panel A: Aging Volcano Ring
setwd(rprojroot::find_rstudio_root_file())
source("04_Figures/shared/style.R")
source("04_Figures/shared/volcano_ring.R")

suppressPackageStartupMessages({
  library(tidyverse)
  library(fgsea)
})

PW <- 190; PH <- 180
RPT <- "04_Figures/F05/b_reports"
DAT <- "04_Figures/F05/c_data"
dir.create(RPT, recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(DAT, "panel_A"), recursive = TRUE, showWarnings = FALSE)

pdf_device <- get_pdf_device()

dep_df <- read_csv("03_DEP/c_data/03_combined_results.csv", show_col_types = FALSE)

fgsea_cache <- "04_Figures/F05/c_data/shared/fgsea_tstat_all_v2.csv"
if (!file.exists(fgsea_cache)) {
  f2_cache <- "04_Figures/F04/c_data/shared/fgsea_tstat_all_v2.csv"
  f1_cache <- "04_Figures/shared/fgsea_f1_panel_F.csv"
  if (file.exists(f2_cache)) {
    dir.create(dirname(fgsea_cache), recursive = TRUE, showWarnings = FALSE)
    file.copy(f2_cache, fgsea_cache)
  } else if (file.exists(f1_cache)) {
    dir.create(dirname(fgsea_cache), recursive = TRUE, showWarnings = FALSE)
    file.copy(f1_cache, fgsea_cache)
  } else {
    stop("fGSEA cache not found")
  }
}
fgsea_all <- read_csv(fgsea_cache, show_col_types = FALSE)

top_terms_aging <- select_ring_terms(fgsea_all, "Aging")
ring_aging      <- build_ring_with_gaps(top_terms_aging, "Aging", fgsea_all)

txt_label <- scale_text(BASE_PATHWAY, PW)

pA <- make_volcano_ring(
  de_df              = dep_df,
  go_df              = fgsea_all,
  contrast           = "Aging",
  title              = NULL,
  contrast_title     = "Aging Effect",
  contrast_subtitle  = "Old_Pre \u2212 Young_Pre",
  title_size         = scale_text(BASE_TAG, PW),
  label_size         = txt_label,
  point_size         = 1.2,
  point_alpha        = 0.55,
  count_label_size   = scale_text(BASE_COUNT, PW),
  ring_data_override = ring_aging
)

ggsave(file.path(RPT, "panel_A_volcano.pdf"), pA,
       width = PW, height = PH, units = "mm", device = pdf_device)
ggsave(file.path(RPT, "panel_A_volcano.png"), pA,
       width = PW, height = PH, units = "mm", dpi = 300)

ring_data_A <- attr(pA, "ring_data")
if (!is.null(ring_data_A) && nrow(ring_data_A) > 0) {
  write_csv(ring_data_A %>% dplyr::select(-gene_list),
            file.path(DAT, "panel_A", "ring_terms.csv"))
}

dep_df %>%
  transmute(
    gene,
    log2_fold_change = round(logFC_Aging, 4),
    neg_log10_pvalue = round(-log10(P.Value_Aging), 4),
    pi_score         = round(pi_score_Aging, 6),
    adjusted_pvalue  = round(adj.P.Val_Aging, 6),
    direction = case_when(
      pi_score_Aging < 0.05 & logFC_Aging > 0 ~ "Up",
      pi_score_Aging < 0.05 & logFC_Aging < 0 ~ "Down",
      TRUE ~ "NS"
    )
  ) %>%
  filter(!is.na(log2_fold_change), !is.na(neg_log10_pvalue)) %>%
  arrange(pi_score) %>%
  write_csv(file.path(DAT, "panel_A", "volcano_aging.csv"))

cat("F5 Panel A done\n")
