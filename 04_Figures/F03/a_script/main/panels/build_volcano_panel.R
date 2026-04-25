# F03 — Shared volcano-ring builder for panels A-D
# Expects: `spec` list with fields: contrast, title, subtitle, tag
# Produces: plot object + standalone PNG/PDF + ring_terms CSV
# Assigns p<tag> into the caller's environment for composite stitching.

# --- Setup (idempotent when sourced by stitcher) -------------------------
setwd(rprojroot::find_rstudio_root_file())
source("04_Figures/shared/style.R")
source("04_Figures/shared/volcano_ring.R")

suppressPackageStartupMessages({
  library(readr); library(dplyr); library(ggplot2)
})

VW <- 89; VH <- 89   # half of 178mm double-col -> no journal scaling
RPT_PNL_PNG <- "04_Figures/F03/b_reports/main/png/panels"
RPT_PNL_PDF <- "04_Figures/F03/b_reports/main/pdf/panels"
DAT     <- "04_Figures/F03/c_data"
dir.create(RPT_PNL_PNG, recursive = TRUE, showWarnings = FALSE)
dir.create(RPT_PNL_PDF, recursive = TRUE, showWarnings = FALSE)
dir.create(DAT,     recursive = TRUE, showWarnings = FALSE)

# Load DEP results once per session (shared by panels A-D when sourced together)
if (!exists("dep_df")) {
  dep_df <- read_csv("03_DEP/c_data/03_combined_results.csv", show_col_types = FALSE)
}

# Load shared fGSEA cache once per session
if (!exists("fgsea_all")) {
  fgsea_cache <- "04_Figures/shared/fgsea_tstat_all_v2.csv"
  stopifnot("fGSEA cache missing: frozen cache -- see 04_Figures/shared/README.md" =
    file.exists(fgsea_cache))
  fgsea_all <- read_csv(fgsea_cache, show_col_types = FALSE)
}

# --- Build volcano ring ---------------------------------------------------
stopifnot("spec must be defined before sourcing build_volcano_panel.R" = exists("spec"))

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

.panel_plot <- make_volcano_ring(
  de_df              = dep_df,
  go_df              = fgsea_all,
  contrast           = spec$contrast,
  title              = NULL,
  contrast_title     = spec$title,
  contrast_subtitle  = enriched_sub,
  ring_data_override = ring_data,
  label_size         = 2.7,
  label_gap          = adaptive_gap,
  title_size         = 5,
  subtitle_size      = 3.5,
  point_size         = 0.5,
  point_alpha        = 0.55,
  count_label_size   = scale_text(BASE_COUNT, VW) + 0.4,
  count_y_mult       = 0.75,
  count_x_mult       = 0.85,
  bg_color           = unname(CONTRAST_COLORS[spec$contrast]),
  bg_alpha           = 0.20,
  show_legend        = FALSE
) + labs(tag = spec$tag)

# --- Save standalone panel ------------------------------------------------
fname <- tolower(gsub(" ", "_", spec$title))
ggsave(file.path(RPT_PNL_PNG, sprintf("MAIN_panel_%s_%s.png", spec$tag, fname)),
       .panel_plot, width = VW, height = VH, units = "mm", dpi = 300)
ggsave(file.path(RPT_PNL_PDF, sprintf("MAIN_panel_%s_%s.pdf", spec$tag, fname)),
       .panel_plot, width = VW, height = VH, units = "mm", device = get_pdf_device())

# --- Write ring terms CSV -------------------------------------------------
ring_out <- attr(.panel_plot, "ring_data")
if (!is.null(ring_out) && nrow(ring_out) > 0) {
  dir.create(file.path(DAT, paste0("panel_", spec$tag)), showWarnings = FALSE)
  write_csv(ring_out %>% dplyr::select(-gene_list,
              -any_of(c("term_idx", "start_deg", "end_deg", "mid_deg",
                         "start_rad", "end_rad", "mid_rad", "arc_r1_var"))),
            file.path(DAT, paste0("panel_", spec$tag), "ring_terms.csv"))
}

# --- Export for composite (strip titles/tags/legend) ----------------------
.panel_plot <- strip_for_composite(.panel_plot)

# Assign p<tag> into the parent environment so the stitcher can find it
.var_name <- paste0("p", spec$tag)
assign(.var_name, .panel_plot, envir = .GlobalEnv)

message(sprintf("F03 panel %s (%s) done", spec$tag, spec$contrast))

# Clean up builder temporaries
rm(.panel_plot, .var_name)
