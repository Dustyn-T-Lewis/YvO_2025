# F03: All Volcano Rings (2x2: Aging, Training Young, Training Old, Interaction)
# Self-contained script — generates all 4 volcano rings and saves individually + composite
# Outputs: panel_A-D individual PDFs/PNGs + F03_volcanoes_MAIN composite

setwd(rprojroot::find_rstudio_root_file())
source("04_Figures/shared/style.R")
source("04_Figures/shared/volcano_ring.R")

library(readr)
library(dplyr)
library(patchwork)
library(cowplot)

VW <- 190; VH <- 180
RPT     <- "04_Figures/F03/b_reports"
RPT_PNL <- "04_Figures/F03/b_reports/panels"
DAT     <- "04_Figures/F03/c_data"
dir.create(RPT_PNL, recursive = TRUE, showWarnings = FALSE)
dir.create(DAT, recursive = TRUE, showWarnings = FALSE)

pdf_device <- get_pdf_device()

dep_df <- read_csv("03_DEP/c_data/03_combined_results.csv", show_col_types = FALSE)

# --- Load fGSEA cache ---
# Prefer F03's local shared/ copy; fall back to F04's live copy.
fgsea_cache <- file.path(DAT, "shared", "fgsea_tstat_all_v2.csv")
if (!file.exists(fgsea_cache)) {
  message("WARNING: F03 fgsea cache missing; copying from F04")
  stopifnot("F04 fgsea cache also missing" =
    file.exists("04_Figures/F04/c_data/shared/fgsea_tstat_all_v2.csv"))
  dir.create(file.path(DAT, "shared"), recursive = TRUE, showWarnings = FALSE)
  file.copy("04_Figures/F04/c_data/shared/fgsea_tstat_all_v2.csv", fgsea_cache)
}
fgsea_all <- read_csv(fgsea_cache, show_col_types = FALSE)

# --- Contrast definitions ---
volcano_specs <- list(
  list(contrast = "Aging",
       title = "Aging Effect",
       subtitle = "Old_Pre \u2212 Young_Pre",
       tag = "A"),
  list(contrast = "Training_Young",
       title = "Training Response (Young)",
       subtitle = "Young_Post \u2212 Young_Pre",
       tag = "B"),
  list(contrast = "Training_Old",
       title = "Training Response (Old)",
       subtitle = "Old_Post \u2212 Old_Pre",
       tag = "C"),
  list(contrast = "Interaction",
       title = "Age \u00d7 Training Interaction",
       subtitle = "Training_Old \u2212 Training_Young",
       tag = "D")
)

# --- Generate each volcano ring ---
volcano_plots <- lapply(volcano_specs, function(spec) {
  message(sprintf("  Building volcano: %s", spec$contrast))

  # Enrich subtitle with DEP + pathway counts
  pi_col <- paste0("pi_score_", spec$contrast)
  n_dep  <- if (pi_col %in% names(dep_df)) sum(dep_df[[pi_col]] < 0.05, na.rm = TRUE) else 0
  n_path <- sum(!is.na(fgsea_all$padj) & fgsea_all$padj < 0.05 &
                fgsea_all$contrast == spec$contrast)
  enriched_sub <- sprintf("%s | %d DEPs, %d pathways", spec$subtitle, n_dep, n_path)

  top_terms <- select_ring_terms(fgsea_all, spec$contrast)
  ring_data <- build_ring_with_gaps(top_terms, spec$contrast, fgsea_all)

  # Adaptive label gap: shorter arcs → labels closer to ring
  max_arc <- max(ring_data$arc_r1_var, na.rm = TRUE)
  adaptive_gap <- 0.7 + 0.3 * (max_arc - 4.8) / 1.6

  p <- make_volcano_ring(
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
    bg_alpha           = 0.12,
    show_legend        = FALSE
  ) + labs(tag = spec$tag)

  # Save individual panel
  fname <- tolower(gsub(" ", "_", spec$title))
  ggsave(file.path(RPT_PNL, sprintf("MAIN_panel_%s_%s.png", spec$tag, fname)),
         p, width = VW, height = VH, units = "mm", dpi = 300)

  # Save ring data
  ring_out <- attr(p, "ring_data")
  if (!is.null(ring_out) && nrow(ring_out) > 0) {
    dir.create(file.path(DAT, paste0("panel_", spec$tag)), showWarnings = FALSE)
    write_csv(ring_out |> dplyr::select(-gene_list),
              file.path(DAT, paste0("panel_", spec$tag), "ring_terms.csv"))
  }

  p
})

# --- NES gradient bar legend (compact, centered below bottom row) ---
# Narrower bar: widen L/R margins so the gradient spans only the center ~60mm
nes_legend <- build_nes_legend_bar(bar_margin = margin(-4, 160, 4, 160, "mm"))

# --- 2x2 composite + NES legend ---
# Strip per-subplot tags so we can place them in the outer canvas margin
# (avoids clipping seen with plot.tag.position inside each ring). Individual
# panel PNGs (saved above) keep their in-plot tags for standalone use.
volcano_plots_nt <- lapply(volcano_plots, function(p) p + labs(tag = NULL))

composite <- ((volcano_plots_nt[[1]] | volcano_plots_nt[[2]]) /
              (volcano_plots_nt[[3]] | volcano_plots_nt[[4]]) /
              wrap_elements(full = nes_legend)) +
  plot_layout(heights = c(0.4970, 0.4970, 0.0060)) &
  # Add top padding to each subplot so the panel title drops ~12mm into
  # the cell, leaving the top canvas band clear for the bold panel letter.
  theme(plot.margin = margin(12, 0, 0, 0, "mm"))

COMP_W <- VW * 2       # 380mm
COMP_H <- 380           # +20mm (10mm per row) to absorb new top padding

# Composite-level tag placement (matches F02/F04/F05/F06/F07 convention)
# Letters sit in the outer top band ABOVE each subplot's title.
#   X_LEFT  ~3.8mm from left edge; X_RIGHT puts tag just past column break.
#   Y_TOP / Y_BOT placed near the top edge of each row's padding band.
TAG_SZ  <- composite_text_sizes(380)$tag
X_LEFT  <- 0.050
X_RIGHT <- 0.505
Y_TOP   <- 0.965 - 4/380                     # -4 mm tag shift (user-tuned)
# Y_BOT offset from row2 top (0.5066) matches Y_TOP offset from row1 top (1.0)
# so C/D sit the same distance above their titles as A/B.
Y_BOT   <- 0.5016 - 4/380

composite <- ggdraw(composite) +
  draw_label("A", x = X_LEFT,  y = Y_TOP, size = TAG_SZ, fontface = "bold", hjust = 0, vjust = 1) +
  draw_label("B", x = X_RIGHT, y = Y_TOP, size = TAG_SZ, fontface = "bold", hjust = 0, vjust = 1) +
  draw_label("C", x = X_LEFT,  y = Y_BOT, size = TAG_SZ, fontface = "bold", hjust = 0, vjust = 1) +
  draw_label("D", x = X_RIGHT, y = Y_BOT, size = TAG_SZ, fontface = "bold", hjust = 0, vjust = 1)

ggsave(file.path(RPT, "MAIN_F03_volcanoes.pdf"), composite,
       width = COMP_W, height = COMP_H, units = "mm", device = pdf_device)
ggsave(file.path(RPT, "MAIN_F03_volcanoes.png"), composite,
       width = COMP_W, height = COMP_H, units = "mm", dpi = 300)

message("F03 all volcanoes done")
