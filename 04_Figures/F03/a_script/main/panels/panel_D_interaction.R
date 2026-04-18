# F03 Panel D: Age x Training interaction volcano ring (Training_Old - Training_Young)
# Self-contained — exposes pD for the composite stitcher.
# Outputs: b_reports/main/{png,pdf}/panels/MAIN_panel_D_age_x_training_interaction.{png,pdf} + c_data/panel_D/ring_terms.csv

setwd(rprojroot::find_rstudio_root_file())
source("04_Figures/shared/style.R")
source("04_Figures/shared/volcano_ring.R")

suppressPackageStartupMessages({
  library(readr); library(dplyr); library(ggplot2)
})

VW <- 190; VH <- 180
RPT_PNL_PNG <- "04_Figures/F03/b_reports/main/png/panels"
RPT_PNL_PDF <- "04_Figures/F03/b_reports/main/pdf/panels"
DAT     <- "04_Figures/F03/c_data"
dir.create(RPT_PNL_PNG, recursive = TRUE, showWarnings = FALSE)
dir.create(RPT_PNL_PDF, recursive = TRUE, showWarnings = FALSE)
dir.create(DAT,     recursive = TRUE, showWarnings = FALSE)

if (!exists("dep_df")) {
  dep_df <- read_csv("03_DEP/c_data/03_combined_results.csv", show_col_types = FALSE)
}

if (!exists("fgsea_all")) {
  fgsea_cache <- "04_Figures/shared/fgsea_tstat_all_v2.csv"
  stopifnot("fGSEA cache missing: frozen cache — see 04_Figures/shared/README.md" =
    file.exists(fgsea_cache))
  fgsea_all <- read_csv(fgsea_cache, show_col_types = FALSE)
}

spec <- list(contrast = "Interaction", title = "Age \u00d7 Training Interaction",
             subtitle = "Training_Old \u2212 Training_Young", tag = "D")

pi_col <- paste0("pi_score_", spec$contrast)
n_dep  <- if (pi_col %in% names(dep_df)) sum(dep_df[[pi_col]] < 0.05, na.rm = TRUE) else 0
n_path <- sum(!is.na(fgsea_all$padj) & fgsea_all$padj < 0.05 &
              fgsea_all$contrast == spec$contrast &
              fgsea_all$database %in% c("Hallmark", "GO Slim", "GO:BP", "KEGG", "Reactome"))
enriched_sub <- sprintf("%s | %d DEPs, %d pathways", spec$subtitle, n_dep, n_path)

top_terms <- select_ring_terms(fgsea_all, spec$contrast)
ring_data <- build_ring_with_gaps(top_terms, spec$contrast, fgsea_all)
max_arc      <- max(ring_data$arc_r1_var, na.rm = TRUE)
adaptive_gap <- 0.7 + 0.3 * (max_arc - 4.8) / 1.6

pD <- make_volcano_ring(
  de_df              = dep_df,
  go_df              = fgsea_all,
  contrast           = spec$contrast,
  title              = NULL,
  contrast_title     = spec$title,
  contrast_subtitle  = enriched_sub,
  ring_data_override = ring_data,
  label_size         = 4.5,
  label_gap          = adaptive_gap,
  title_size         = composite_text_sizes(380)$title,
  subtitle_size      = composite_text_sizes(380)$subtitle,
  point_size         = 1.2,
  point_alpha        = 0.55,
  count_label_size   = scale_text(BASE_COUNT + 1, VW),
  bg_color           = unname(CONTRAST_COLORS[spec$contrast]),
  bg_alpha           = 0.20,
  show_legend        = FALSE
) + labs(tag = spec$tag)

fname <- tolower(gsub(" ", "_", spec$title))
ggsave(file.path(RPT_PNL_PNG, sprintf("MAIN_panel_%s_%s.png", spec$tag, fname)),
       pD, width = VW, height = VH, units = "mm", dpi = 300)
ggsave(file.path(RPT_PNL_PDF, sprintf("MAIN_panel_%s_%s.pdf", spec$tag, fname)),
       pD, width = VW, height = VH, units = "mm", device = get_pdf_device())

ring_out <- attr(pD, "ring_data")
if (!is.null(ring_out) && nrow(ring_out) > 0) {
  dir.create(file.path(DAT, paste0("panel_", spec$tag)), showWarnings = FALSE)
  write_csv(ring_out %>% dplyr::select(-gene_list),
            file.path(DAT, paste0("panel_", spec$tag), "ring_terms.csv"))
}

message("F03 panel D (Interaction) done")
