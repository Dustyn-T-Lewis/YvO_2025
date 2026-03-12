# F3 Panel B: Reversal Volcano Ring
setwd(rprojroot::find_rstudio_root_file())
source("04_Figures/shared/style.R")
source("04_Figures/shared/volcano_ring.R")

suppressPackageStartupMessages({
  library(tidyverse)
  library(fgsea)
})

PW <- 190; PH <- 180
RPT <- "04_Figures/F3/b_reports"
DAT <- "04_Figures/F3/c_data"
dir.create(RPT, recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(DAT, "panel_B"), recursive = TRUE, showWarnings = FALSE)

pdf_device <- get_pdf_device()

dep_df <- read_csv("03_DEP/c_data/03_combined_results.csv", show_col_types = FALSE)

fgsea_cache <- "04_Figures/F3/c_data/shared/fgsea_tstat_all_v2.csv"
if (!file.exists(fgsea_cache)) {
  f2_cache <- "04_Figures/F2/c_data/shared/fgsea_tstat_all_v2.csv"
  f1_cache <- "04_Figures/F1/c_data/06_panel_F_fgsea_results.csv"
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

if (!"Reversal" %in% unique(fgsea_all$contrast)) {
  stop("Reversal contrast not found in fGSEA cache. Run F3 setup or panel_A first.")
}

top_terms_rev <- select_ring_terms(fgsea_all, "Reversal", n_each = 6)
ring_rev      <- build_ring_with_gaps(top_terms_rev, "Reversal", fgsea_all, n_each = 6)

txt_label <- scale_text(BASE_PATHWAY, PW)

pB <- make_volcano_ring(
  de_df              = dep_df,
  go_df              = fgsea_all,
  contrast           = "Reversal",
  title              = NULL,
  contrast_title     = "Reversal (Training \u2212 Aging)",
  contrast_subtitle  = "Training_Old \u2212 Aging",
  title_size         = scale_text(BASE_TAG, PW),
  label_size         = txt_label,
  point_size         = 1.2,
  point_alpha        = 0.55,
  count_label_size   = scale_text(BASE_COUNT, PW),
  ring_data_override = ring_rev
)

ggsave(file.path(RPT, "panel_B_volcano.pdf"), pB,
       width = PW, height = PH, units = "mm", device = pdf_device)
ggsave(file.path(RPT, "panel_B_volcano.png"), pB,
       width = PW, height = PH, units = "mm", dpi = 300)

ring_data_B <- attr(pB, "ring_data")
if (!is.null(ring_data_B) && nrow(ring_data_B) > 0) {
  write_csv(ring_data_B %>% dplyr::select(-gene_list),
            file.path(DAT, "panel_B", "ring_terms.csv"))
}

dep_df %>%
  transmute(
    gene,
    log2_fold_change = round(logFC_Reversal, 4),
    neg_log10_pvalue = round(-log10(P.Value_Reversal), 4),
    pi_score         = round(pi_score_Reversal, 6),
    adjusted_pvalue  = round(adj.P.Val_Reversal, 6),
    direction = case_when(
      pi_score_Reversal < 0.05 & logFC_Reversal > 0 ~ "Up",
      pi_score_Reversal < 0.05 & logFC_Reversal < 0 ~ "Down",
      TRUE ~ "NS"
    )
  ) %>%
  filter(!is.na(log2_fold_change), !is.na(neg_log10_pvalue)) %>%
  arrange(pi_score) %>%
  write_csv(file.path(DAT, "panel_B", "volcano_reversal.csv"))

cat("F3 Panel B done\n")
