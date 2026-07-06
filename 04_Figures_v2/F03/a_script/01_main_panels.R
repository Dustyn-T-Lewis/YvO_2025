#!/usr/bin/env Rscript
# F03 Main — 4 volcano rings (Aging, Training_Young, Training_Old, Interaction)
# 2×2 composite + NES gradient legend

withr::local_dir(rprojroot::find_root(rprojroot::has_file("setup.R")))

library(readr)
library(dplyr)
library(ggplot2)
library(patchwork)
library(cowplot)

source("04_Figures_v2/shared/style.R")
source("04_Figures_v2/shared/volcano_ring.R")
source("04_Figures_v2/shared/build_fgsea_cache.R")

BUILDER <- "04_Figures_v2/F03/a_script/_build_volcano_panel.R"

PI_THRESH <- 0.05 # pi-score significance threshold
FDR_THRESH <- 0.05 # BH-adjusted p threshold for fGSEA
DBS_USED <- c("Hallmark", "GO Slim", "GO:BP", "KEGG", "Reactome")

spec <- list(
  contrast = "Aging", title = "Aging Effect",
  subtitle = "Old_Pre \u2212 Young_Pre", tag = "A"
)
source(BUILDER)

spec <- list(
  contrast = "Training_Young", title = "Training Response (Young)",
  subtitle = "Young_Post \u2212 Young_Pre", tag = "B"
)
source(BUILDER)

spec <- list(
  contrast = "Training_Old", title = "Training Response (Old)",
  subtitle = "Old_Post \u2212 Old_Pre", tag = "C"
)
source(BUILDER)

spec <- list(
  contrast = "Interaction", title = "Age \u00d7 Training Interaction",
  subtitle = "Training_Old \u2212 Training_Young", tag = "D"
)
source(BUILDER)

RPT_PDF <- "04_Figures_v2/F03/b_reports/main/pdf"
RPT_PNG <- "04_Figures_v2/F03/b_reports/main/png"
dir.create(RPT_PDF, recursive = TRUE, showWarnings = FALSE)
dir.create(RPT_PNG, recursive = TRUE, showWarnings = FALSE)

contrast_stats <- function(ctr) {
  pi_col <- paste0("pi_score_", ctr)
  n_dep <- if (pi_col %in% names(dep_df)) sum(dep_df[[pi_col]] < PI_THRESH, na.rm = TRUE) else 0
  ctr_rows <- fgsea_all[fgsea_all$contrast == ctr & fgsea_all$database %in% DBS_USED, ]
  sprintf(
    "%d DEPs (\u03a0 < 0.05)  |  %d / %d pathways (FDR < 0.05)",
    n_dep, sum(ctr_rows$padj < FDR_THRESH, na.rm = TRUE), sum(!is.na(ctr_rows$padj))
  )
}

nes_legend <- build_nes_legend_bar(
  text_size = 5, title_size = 5,
  bar_margin = margin(0, 0, 0, 0, "mm")
)

pA <- pA + theme(plot.margin = margin(4, -9, 0, 9, "mm"))
pB <- pB + theme(plot.margin = margin(4, 0, 0, 0, "mm"))
pC <- pC + theme(plot.margin = margin(0, -9, 0, 9, "mm"))
pD <- pD + theme(plot.margin = margin(0, 0, 0, 0, "mm"))

composite <- ((pA | pB) / (pC | pD)) + plot_layout(heights = c(1, 1))

layout_cfg <- list(
  # Canvas dimensions (mm)
  w = 178, h = 180,
  # Text size adjustments relative to composite_text_sizes() base
  tag_bump = 4, # added to txt$tag
  ttl_bump = 2, # added to txt$title
  sub_bump = 2, # added to txt$subtitle
  # Label x positions: left and right columns, plus title offset from tag
  x_l = 0.070, # left-column tag x (panels A, C)
  x_r = 0.510, # right-column tag x (panels B, D)
  x_ttl = 0.040, # title x offset from tag
  # Label y positions: top and bottom rows; subtitle drops by sub_off
  y_top = 0.960,
  y_bot = 0.505,
  sub_off = 0.022,
  # NES legend position (draw_plot args)
  nes_x = 0.35, nes_y = 0.058, nes_w = 0.30, nes_h = 0.028,
  # Methods footnote (bottom gutter, below the NES legend)
  note_x = 0.5, note_y = 0.038, note_lh = 0.9
)

COMP_W <- layout_cfg$w
COMP_H <- layout_cfg$h
txt <- composite_text_sizes(COMP_H)
TAG_SZ <- txt$tag + layout_cfg$tag_bump
TTL_SZ <- txt$title + layout_cfg$ttl_bump
SUB_SZ <- txt$subtitle + layout_cfg$sub_bump
X_L <- layout_cfg$x_l
X_R <- layout_cfg$x_r
X_TTL <- layout_cfg$x_ttl
Y_TOP <- layout_cfg$y_top
Y_BOT <- layout_cfg$y_bot
SUB_OFF <- layout_cfg$sub_off

titles <- c(
  "Aging Effect", "Training Response (Young)",
  "Training Response (Old)", "Age \u00d7 Training Interaction"
)
subs <- c(
  contrast_stats("Aging"), contrast_stats("Training_Young"),
  contrast_stats("Training_Old"), contrast_stats("Interaction")
)
tags <- LETTERS[1:4]
xs <- c(X_L, X_R, X_L, X_R)
ys <- c(Y_TOP, Y_TOP, Y_BOT, Y_BOT)

composite <- ggdraw(composite)
for (i in seq_along(tags)) {
  composite <- composite +
    draw_label(tags[i],
      x = xs[i], y = ys[i] + 0.002, size = TAG_SZ,
      fontface = "bold", hjust = 0, vjust = 1
    ) +
    draw_label(titles[i],
      x = xs[i] + X_TTL, y = ys[i], size = TTL_SZ,
      fontface = "bold", hjust = 0, vjust = 1
    ) +
    draw_label(subs[i],
      x = xs[i] + X_TTL, y = ys[i] - SUB_OFF, size = SUB_SZ,
      fontface = "bold.italic", colour = "grey40", hjust = 0, vjust = 1
    )
}
composite <- composite +
  draw_plot(nes_legend,
    x = layout_cfg$nes_x, y = layout_cfg$nes_y,
    width = layout_cfg$nes_w, height = layout_cfg$nes_h
  )

fig_note <- paste0(
  "Point colour marks the π < 0.05 frontier (π = P^|log2FC|, a curved ",
  "magnitude–significance boundary); volcano axes are per-contrast relative, not absolute.\n",
  "Ring arc height encodes term rank, not padj. Panel A is a cross-sectional baseline age ",
  "difference (Old_Pre − Young_Pre); B–D are longitudinal training and interaction effects."
)
composite <- composite +
  draw_label(fig_note,
    x = layout_cfg$note_x, y = layout_cfg$note_y, size = txt$subtitle,
    colour = "grey40", hjust = 0.5, vjust = 1, lineheight = layout_cfg$note_lh
  )

ggsave(file.path(RPT_PDF, "MAIN_F03_composite.pdf"), composite,
  width = COMP_W, height = COMP_H, units = "mm", device = get_pdf_device(),
  limitsize = FALSE
)
ggsave(file.path(RPT_PNG, "MAIN_F03_composite.png"), composite,
  width = COMP_W, height = COMP_H, units = "mm", dpi = 300, limitsize = FALSE
)

message("F03 main composite done")
