# F04 Supplementary Composite: Diagnostic Panels (A-E)
# 3+2 grid layout:
#   Top row:    A (ORA dedup) | B (rho bootstrap) | C (threshold sens.)
#   Bottom row: D (GO Slim bars)                   | E (fry leading edge)
setwd(rprojroot::find_rstudio_root_file())
source("04_Figures/shared/style.R")

suppressPackageStartupMessages({
  library(patchwork)
  library(cowplot)
})

# ── Source panels (before output dirs to avoid RPT_PDF/RPT_PNG clobbering) ──
message("=== F04 SUPP Composite: sourcing panels ===")
source("04_Figures/F04/a_script/supp/panels/SUPP_ora_dedup.R")       # -> pS_ora_dedup
source("04_Figures/F04/a_script/supp/panels/SUPP_rho_bootstrap.R")   # -> pS_rho_boot
source("04_Figures/F04/a_script/supp/panels/SUPP_threshold_sens.R")  # -> pS_thresh
source("04_Figures/F04/a_script/supp/panels/SUPP_goslim_bars.R")     # -> pS_goslim
source("04_Figures/F04/a_script/supp/panels/SUPP_fry_leading.R")     # -> pS_fry_lead

# ── Output directories (set after panels so their RPT_PDF/RPT_PNG don't persist)
RPT_PDF <- "04_Figures/F04/b_reports/supp/pdf"
RPT_PNG <- "04_Figures/F04/b_reports/supp/png"
BOX_DIR <- "02_Figures/F04_training_concordance/supp"
dir.create(RPT_PDF, recursive = TRUE, showWarnings = FALSE)
dir.create(RPT_PNG, recursive = TRUE, showWarnings = FALSE)
dir.create(BOX_DIR, recursive = TRUE, showWarnings = FALSE)

pdf_device <- get_pdf_device()

# ── Composite layout ─────────────────────────────────────────────────────────
COMP_W <- 280    # wide — 3 panels across need room
COMP_H <- 275    # tall enough for 2 rows with titles + bottom legend

cts <- composite_text_sizes(COMP_H)
TAG_SZ <- cts$tag
TTL_SZ <- cts$title + 1    # slightly larger for readability at this width
SUB_SZ <- cts$subtitle + 0.5

# Fix individual panels before combining
# Panel C: angle x-axis labels so threshold names don't overlap
pS_thresh <- pS_thresh +
  theme(axis.text.x = element_text(angle = 40, hjust = 1, vjust = 1, size = 7))

# Tighten y-axis title spacing + suppress legends
axis_fix <- theme(legend.position = "none",
                  axis.title.y = element_text(margin = margin(0, 2, 0, 0)),
                  axis.title.x = element_text(margin = margin(2, 0, 0, 0)))
pS_ora_dedup <- pS_ora_dedup + axis_fix +
  theme(plot.margin = margin(18, 10, 2, 5))
pS_rho_boot  <- pS_rho_boot + axis_fix +
  theme(plot.margin = margin(18, 10, 2, 5))
pS_thresh    <- pS_thresh + axis_fix +
  theme(plot.margin = margin(18, 14, 2, 5))
pS_goslim    <- pS_goslim + axis_fix +
  theme(plot.margin = margin(18, 10, 2, 5))
pS_fry_lead  <- pS_fry_lead + axis_fix +
  theme(plot.margin = margin(18, 10, 2, 5))

# Build grid — guides = "keep" prevents patchwork from collecting/displaying legends
top_row <- (pS_ora_dedup | pS_rho_boot | pS_thresh) +
  plot_layout(widths = c(1, 1, 1.15), guides = "keep")
bot_row <- (pS_goslim | pS_fry_lead) +
  plot_layout(guides = "keep")

fig <- (top_row / bot_row) +
  plot_layout(heights = c(1, 1.05), guides = "keep")

# ── Panel tags, titles, subtitles via cowplot ────────────────────────────────
# Top row: A/B/C left edges aligned with patchwork column left borders
X_A <- 0.020;  X_B <- 0.350;  X_C <- 0.675
# Bottom row: D/E left edges aligned with patchwork column left borders
X_D <- 0.020;  X_E <- 0.520
X_TTL      <- 0.022        # title offset right of tag
SUB_OFFSET <- 0.016        # subtitle gap below title
Y_TOP <- 0.993;  Y_BOT <- 0.513

composite_final <- ggdraw(fig) +
  # Panel A
  draw_label("A", x = X_A,         y = Y_TOP,               size = TAG_SZ, fontface = "bold", hjust = 0, vjust = 1) +
  draw_label("ORA Dedup Sensitivity", x = X_A + X_TTL, y = Y_TOP, size = TTL_SZ, fontface = "bold", hjust = 0, vjust = 1) +
  draw_label("Enriched pathways by Jaccard cutoff", x = X_A + X_TTL, y = Y_TOP - SUB_OFFSET, size = SUB_SZ, fontface = "bold.italic", hjust = 0, vjust = 1, colour = "grey40") +
  # Panel B
  draw_label("B", x = X_B,         y = Y_TOP,               size = TAG_SZ, fontface = "bold", hjust = 0, vjust = 1) +
  draw_label("Spearman \u03c1 Bootstrap", x = X_B + X_TTL, y = Y_TOP, size = TTL_SZ, fontface = "bold", hjust = 0, vjust = 1) +
  draw_label("1000 replicates, 95% CI", x = X_B + X_TTL, y = Y_TOP - SUB_OFFSET, size = SUB_SZ, fontface = "bold.italic", hjust = 0, vjust = 1, colour = "grey40") +
  # Panel C
  draw_label("C", x = X_C,         y = Y_TOP,               size = TAG_SZ, fontface = "bold", hjust = 0, vjust = 1) +
  draw_label("Threshold Sensitivity", x = X_C + X_TTL, y = Y_TOP, size = TTL_SZ, fontface = "bold", hjust = 0, vjust = 1) +
  draw_label("Quadrant counts by significance criterion", x = X_C + X_TTL, y = Y_TOP - SUB_OFFSET, size = SUB_SZ, fontface = "bold.italic", hjust = 0, vjust = 1, colour = "grey40") +
  # Panel D
  draw_label("D", x = X_D,         y = Y_BOT,               size = TAG_SZ, fontface = "bold", hjust = 0, vjust = 1) +
  draw_label("GO Slim Distribution", x = X_D + X_TTL, y = Y_BOT, size = TTL_SZ, fontface = "bold", hjust = 0, vjust = 1) +
  draw_label("Proteins per GO Slim category by quadrant", x = X_D + X_TTL, y = Y_BOT - SUB_OFFSET, size = SUB_SZ, fontface = "bold.italic", hjust = 0, vjust = 1, colour = "grey40") +
  # Panel E
  draw_label("E", x = X_E,         y = Y_BOT,               size = TAG_SZ, fontface = "bold", hjust = 0, vjust = 1) +
  draw_label("fry Leading-Edge Proteins", x = X_E + X_TTL, y = Y_BOT, size = TTL_SZ, fontface = "bold", hjust = 0, vjust = 1) +
  draw_label("Top 20 by |t-stat| in Training Old", x = X_E + X_TTL, y = Y_BOT - SUB_OFFSET, size = SUB_SZ, fontface = "bold.italic", hjust = 0, vjust = 1, colour = "grey40")

# ── Consolidated legend bar at bottom ────────────────────────────────────────
# Build a tiny ggplot with the shared legend, extract, and overlay
legend_df <- data.frame(
  x = 1:3,
  fill_lab = factor(c("Concordant Up", "Concordant Down", "Discordant"),
                    levels = c("Concordant Up", "Concordant Down", "Discordant")))
legend_plot <- ggplot(legend_df, aes(x = x, y = 1, fill = fill_lab)) +
  geom_col() +
  scale_fill_manual(values = c("Concordant Up" = "#E57373",
                               "Concordant Down" = "#64B5F6",
                               "Discordant" = "#FFB74D"),
                    name = "Quadrant") +
  theme_void() +
  theme(legend.position = "bottom",
        legend.direction = "horizontal",
        legend.title = element_text(face = "bold", size = 8),
        legend.text = element_text(size = 7),
        legend.key.size = unit(4, "mm"),
        legend.spacing.x = unit(3, "mm"))

# Extract just the legend grob
legend_grob <- cowplot::get_legend(legend_plot)
composite_final <- composite_final +
  draw_plot(legend_grob, x = 0.08, y = -0.003, width = 0.45, height = 0.035)

# ── Save ─────────────────────────────────────────────────────────────────────
ggsave(file.path(RPT_PDF, "SUPP_F04_diagnostics.pdf"), composite_final,
       width = COMP_W, height = COMP_H, units = "mm", device = pdf_device)
ggsave(file.path(RPT_PNG, "SUPP_F04_diagnostics.png"), composite_final,
       width = COMP_W, height = COMP_H, units = "mm", dpi = 300)

# Copy to Box
BOX <- "/Users/dtl0018/Library/CloudStorage/Box-Box/YvO_proteomics_manuscript"
box_pdf <- file.path(BOX, BOX_DIR, "pdf")
box_png <- file.path(BOX, BOX_DIR, "png")
dir.create(box_pdf, recursive = TRUE, showWarnings = FALSE)
dir.create(box_png, recursive = TRUE, showWarnings = FALSE)
file.copy(file.path(RPT_PDF, "SUPP_F04_diagnostics.pdf"),
          file.path(box_pdf, "SUPP_F04_diagnostics.pdf"), overwrite = TRUE)
file.copy(file.path(RPT_PNG, "SUPP_F04_diagnostics.png"),
          file.path(box_png, "SUPP_F04_diagnostics.png"), overwrite = TRUE)

message(sprintf("F04 SUPP composite (5-panel) saved: %s x %s mm", COMP_W, COMP_H))
message(sprintf("  PDF: %s/SUPP_F04_diagnostics.pdf", RPT_PDF))
message(sprintf("  PNG: %s/SUPP_F04_diagnostics.png", RPT_PNG))
message(sprintf("  Box: %s", BOX_DIR))
