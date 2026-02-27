################################################################################
#   YvO Figure 5 — Panel A: Protein Dendrogram & Module Colors
#   Shows: dendrogram, Dynamic Tree Cut colors, Merged module colors
#   Rendered via base-R plotDendroAndColors, embedded as raster in ggplot
################################################################################
#
# ── STAT AUDIT (2026-02-27) ──────────────────────────────────────────────────
#
# Panel A is a visualization of the WGCNA dendrogram — no inferential
# statistics are displayed.  Key parameters (signed network, power = 14,
# minModuleSize = 30, mergeCutHeight = 0.25) are set in YvO_WGCNA_run.R.
#
# No changes required.
# ─────────────────────────────────────────────────────────────────────────────

source("04_Figures/F5/a_script/YvO_F5_setup.R")

# --- Retrieve soft-thresholding power from the WGCNA run ---
# (power is not stored in net; re-derive from the script's default logic)
soft_power <- 14   # signed network, R^2 > 0.85 (see YvO_WGCNA_run.R)

# --- Prepare color bars: Dynamic Tree Cut (unmerged) + Merged modules ---
block_genes  <- net$blockGenes[[1]]
merged_cols  <- module_colors[block_genes]          # final merged assignments
unmerged_cols <- net$unmergedColors[block_genes]     # pre-merge assignments

color_matrix <- cbind(unmerged_cols, merged_cols)
color_labels <- c("Dynamic Tree Cut", "Merged Modules")

# --- Module summary for subtitle ---
n_mods   <- length(unique(merged_cols[merged_cols != "grey"]))
n_genes  <- length(merged_cols)
n_grey   <- sum(merged_cols == "grey")

# --- Render base-R dendrogram to temp PNG ---
dendro_tmp <- tempfile(fileext = ".png")
tryCatch({
  png(dendro_tmp, width = 2800, height = 1600, res = 300)
  par(mar = c(1, 5, 1, 1))
  plotDendroAndColors(net$dendrograms[[1]],
                      color_matrix,
                      color_labels,
                      main = "",
                      dendroLabels = FALSE, hang = 0.03,
                      addGuide = TRUE, guideHang = 0.05,
                      cex.colorLabels = 0.7,
                      cex.axis = 0.8)
  dev.off()
}, error = function(e) {
  try(dev.off(), silent = TRUE)
  message("Panel A dendrogram render failed: ", e$message)
})

if (file.exists(dendro_tmp) && file.size(dendro_tmp) > 0) {
  dendro_img <- readPNG(dendro_tmp)
} else {
  dendro_img <- NULL
}

# --- Wrap in ggplot for consistent figure pipeline ---
subtitle_text <- sprintf(
  "Signed network | power = %d | %d modules | %s proteins (%d unassigned)",
  soft_power, n_mods, format(n_genes, big.mark = ","), n_grey
)

pA <- ggplot() +
  { if (!is.null(dendro_img))
      annotation_raster(dendro_img, xmin = 0, xmax = 1, ymin = 0, ymax = 1)
    else
      annotate("text", x = 0.5, y = 0.5, label = "Dendrogram not available",
               size = 4, color = "grey50")
  } +
  labs(title    = "A  Protein Dendrogram & Module Colors",
       subtitle = subtitle_text) +
  coord_cartesian(xlim = c(0, 1), ylim = c(0, 1), expand = FALSE) +
  theme_void() +
  theme(plot.title    = element_text(face = "bold", size = 9),
        plot.subtitle = element_text(size = 6.5, color = "grey30", face = "italic"),
        plot.margin   = margin(2, 2, 2, 2))

# --- Export panel_A data ---
dendro_data <- tibble(
  uniprot_id    = names(module_colors)[block_genes],
  unmerged_color = unmerged_cols,
  merged_color   = merged_cols
)
write_csv(dendro_data, file.path(DAT_DIR, "01_panel_A_dendrogram_data.csv"))

# --- Save (2:1 aspect ratio) ---
ggsave(file.path(RPT_DIR, "panel_A_dendrogram.pdf"), pA,
       width = 240, height = 120, units = "mm",
       device = pdf, limitsize = FALSE)
ggsave(file.path(RPT_DIR, "panel_A_dendrogram.png"), pA,
       width = 240, height = 120, units = "mm",
       dpi = 300, limitsize = FALSE)

cat("Panel A saved\n")
