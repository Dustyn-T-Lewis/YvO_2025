# F03 — Single-Column Composite (4 volcano rings in 2x2 at column width)
# Re-calls make_volcano_ring() with scaled-down text sizes for ~84mm width.
# Does NOT modify panel scripts or shared volcano_ring.R.

setwd(rprojroot::find_rstudio_root_file())
source("04_Figures/shared/style.R")
source("04_Figures/shared/volcano_ring.R")

suppressPackageStartupMessages({
  library(patchwork)
  library(cowplot)
  library(ggplot2)
  library(readr)
  library(dplyr)
})

RPT_PDF     <- "04_Figures/F03/b_reports/main/pdf"
RPT_PNG     <- "04_Figures/F03/b_reports/main/png"
pdf_device <- get_pdf_device()

# --- Load data (same as panel scripts) ---
dep_df    <- read_csv("03_DEP/c_data/03_combined_results.csv", show_col_types = FALSE)
fgsea_all <- read_csv("04_Figures/shared/fgsea_tstat_all_v2.csv", show_col_types = FALSE)

# --- Build each panel at compact scale ---
specs <- list(
  list(contrast = "Aging",          title = "Aging Effect",
       subtitle = "Old_Pre \u2212 Young_Pre", tag = "A"),
  list(contrast = "Training_Young", title = "Training Response (Young)",
       subtitle = "Young_Post \u2212 Young_Pre", tag = "B"),
  list(contrast = "Training_Old",   title = "Training Response (Old)",
       subtitle = "Old_Post \u2212 Old_Pre", tag = "C"),
  list(contrast = "Interaction",    title = "Age \u00d7 Training Interaction",
       subtitle = "Training_Old \u2212 Training_Young", tag = "D")
)

panels_sc <- lapply(specs, function(sp) {
  # Count DEPs and pathways for subtitle
  pi_col <- paste0("pi_score_", sp$contrast)
  n_dep  <- if (pi_col %in% names(dep_df)) sum(dep_df[[pi_col]] < 0.05, na.rm = TRUE) else 0
  n_path <- sum(!is.na(fgsea_all$padj) & fgsea_all$padj < 0.05 &
                fgsea_all$contrast == sp$contrast &
                fgsea_all$database %in% c("Hallmark", "GO Slim", "GO:BP", "KEGG", "Reactome"))
  enriched_sub <- sprintf("%s | %d DEPs, %d pathways", sp$subtitle, n_dep, n_path)

  # Build ring data — limit to 4 up + 4 down terms for clean spacing at compact scale
  top_terms <- select_ring_terms(fgsea_all, sp$contrast, n_each = 4)
  ring_data <- build_ring_with_gaps(top_terms, sp$contrast, fgsea_all)
  max_arc   <- max(ring_data$arc_r1_var, na.rm = TRUE)
  # Larger gap than full-width version to give labels room at compact scale
  adaptive_gap <- 1.0 + 0.4 * (max_arc - 4.8) / 1.6

  # Re-call make_volcano_ring with scaled-down sizes
  p <- make_volcano_ring(
    de_df              = dep_df,
    go_df              = fgsea_all,
    contrast           = sp$contrast,
    title              = NULL,
    contrast_title     = sp$title,
    contrast_subtitle  = enriched_sub,
    ring_data_override = ring_data,
    label_size          = 1.6,            # readable at 42mm panel width without edge clipping
    label_gap           = adaptive_gap,
    label_padding       = 1,             # tight box around pathway labels
    min_angle_gap       = 32,            # wider gap threshold to prevent overlap at compact scale
    nudge_outward       = 2.4,           # push overlapping labels further out
    title_size          = 5,
    subtitle_size       = 3,
    point_size          = 0.3,
    point_alpha         = 0.45,
    count_label_size    = 1.0,
    count_label_padding = 1,             # tight box around DEP counts
    count_border_width  = 0.1,           # hairline border
    count_y_mult        = 0.80,          # within upper scatter region
    count_x_mult        = 0.78,          # near ring edges
    bg_color            = unname(CONTRAST_COLORS[sp$contrast]),
    bg_alpha            = 0.20,
    show_legend         = FALSE
  ) + labs(tag = NULL) +
    theme(
      plot.title    = element_text(margin = margin(b = 0, unit = "mm")),
      plot.subtitle = element_text(margin = margin(b = 0, unit = "mm")),
      plot.margin   = margin(-1, 0.5, -1, 0.5, "mm")
    )

  p
})

# --- NES gradient bar (taller bar, narrower width) ---
nes_legend_sc <- build_nes_legend_bar(
  text_size  = 5,
  title_size = 5.5,
  bar_margin = margin(-1, 22, 0, 22, "mm")
)

# --- 2x2 composite + NES legend ---
composite_sc <- ((panels_sc[[1]] | panels_sc[[2]]) /
                 (panels_sc[[3]] | panels_sc[[4]]) /
                 wrap_elements(full = nes_legend_sc)) +
  plot_layout(heights = c(0.473, 0.473, 0.054)) &
  theme(plot.margin = margin(0, 0, 0, 0, "mm"))

# Dimensions: column width ~84mm, compact height
SC_W <- 85
SC_H <- 90

# Add tags
composite_sc <- ggdraw(composite_sc) +
  draw_label("A", x = 0.02, y = 0.97, size = 6, fontface = "bold", hjust = 0, vjust = 1) +
  draw_label("B", x = 0.51, y = 0.97, size = 6, fontface = "bold", hjust = 0, vjust = 1) +
  draw_label("C", x = 0.02, y = 0.50, size = 6, fontface = "bold", hjust = 0, vjust = 1) +
  draw_label("D", x = 0.51, y = 0.50, size = 6, fontface = "bold", hjust = 0, vjust = 1)

dir.create(RPT_PDF, recursive = TRUE, showWarnings = FALSE)
dir.create(RPT_PNG, recursive = TRUE, showWarnings = FALSE)
ggsave(file.path(RPT_PDF, "MAIN_F03_composite_single_col.pdf"), composite_sc,
       width = SC_W, height = SC_H, units = "mm", device = pdf_device)
ggsave(file.path(RPT_PNG, "MAIN_F03_composite_single_col.png"), composite_sc,
       width = SC_W, height = SC_H, units = "mm", dpi = 300)

cat("F03 single-column composite done\n")
