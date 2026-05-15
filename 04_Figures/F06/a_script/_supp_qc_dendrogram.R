# Sourced by 02_supp_panels.R — expects style.R already loaded.

library(readr)
library(dplyr)
library(tibble)
library(WGCNA)
library(png)

BASE <- "04_Figures/F06"

RPT_PNG <- file.path(BASE, "b_reports", "supp", "png", "panels")
RPT_PDF <- file.path(BASE, "b_reports", "supp", "pdf", "panels")
DAT     <- file.path(BASE, "c_data")
dir.create(RPT_PNG, recursive = TRUE, showWarnings = FALSE)
dir.create(RPT_PDF, recursive = TRUE, showWarnings = FALSE)
dir.create(DAT, recursive = TRUE, showWarnings = FALSE)

pdf_device <- get_pdf_device()

net           <- readRDS(file.path(DAT, "wgcna/wgcna_network.rds"))
module_df     <- read_csv(file.path(DAT, "wgcna/wgcna_module_assignments.csv"))
module_colors <- readRDS(file.path(DAT, "module_colors.rds"))

message("Supplementary: protein dendrogram & module colors...")

PA_W <- 240
PA_H <- 120

txt_title <- scale_text(BASE_STAT, PA_W)
txt_sub   <- scale_text(BASE_GENE, PA_W)

sft_csv <- read.csv(file.path(DAT, "wgcna/wgcna_sft_summary.csv"))
soft_power <- sft_csv$selected_power[1]

block_genes  <- net$blockGenes[[1]]
merged_cols  <- module_colors[block_genes]
unmerged_cols <- net$unmergedColors[block_genes]

color_matrix <- cbind(unmerged_cols, merged_cols)
color_labels <- c("Dynamic Tree Cut", "Merged Modules")

n_mods   <- length(unique(merged_cols[merged_cols != "grey"]))
n_genes  <- length(merged_cols)
n_grey   <- sum(merged_cols == "grey")

dendro_tmp <- tempfile(fileext = ".png")
tryCatch({
  png(dendro_tmp, width = 3200, height = 1600, res = 300)
  par(mar = c(1, 4, 1, 0.5))
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
  message("Supplementary dendrogram render failed: ", e$message)
})

if (file.exists(dendro_tmp) && file.size(dendro_tmp) > 0) {
  dendro_img <- readPNG(dendro_tmp)
} else {
  dendro_img <- NULL
}

subtitle_text <- sprintf(
  "Signed network | power = %d | %d modules | %s proteins (%d unassigned)",
  soft_power, n_mods, format(n_genes, big.mark = ","), n_grey
)

pA <- ggplot() +
  { if (!is.null(dendro_img))
      annotation_raster(dendro_img, xmin = 0, xmax = 1, ymin = 0, ymax = 1)
    else
      annotate("text", x = 0.5, y = 0.5, label = "Dendrogram not available",
               size = txt_title, color = "grey50")
  } +
  labs(title    = "Protein Dendrogram & Module Colors",
       subtitle = subtitle_text) +
  coord_cartesian(xlim = c(0, 1), ylim = c(0, 1), expand = FALSE) +
  theme_void() +
  theme(plot.title    = element_text(face = "bold", size = 13),
        plot.subtitle = element_text(size = 10, color = "grey30",
                                     face = "italic"),
        plot.margin   = margin(2, 2, 2, 2),
        legend.position = "none")

dendro_data <- tibble(
  uniprot_id    = names(module_colors)[block_genes],
  unmerged_color = unmerged_cols,
  merged_color   = merged_cols
)
write_csv(dendro_data, file.path(DAT, "asupp_B_QC_dendrogram_SUPP_data.csv"))

ggsave(file.path(RPT_PNG, "SUPP_dendrogram.png"), pA,
       width = PA_W, height = PA_H, units = "mm", dpi = 300)
ggsave(file.path(RPT_PDF, "SUPP_dendrogram.pdf"), pA,
       width = PA_W, height = PA_H, units = "mm", device = pdf_device)

message("  Supplementary dendrogram saved")
