# F06 SUPP — Hub-protein chord diagram across KEY_MODULES
#
# Alternate representation to the per-module hub networks (panel_D_hub.R):
# a single circlize chord with one sector per key module and chords linking
# top-kME hub proteins to their module sector. Width-encodes kME.
#
# Output: 04_Figures/F06/b_reports/supp/03_module/SUPP_hub_chord.{pdf,png}

setwd(rprojroot::find_rstudio_root_file())
source("04_Figures/shared/style.R")

library(readr)
library(dplyr)
library(circlize)
library(stringr)

DAT <- "04_Figures/F06/c_data"
RPT <- "04_Figures/F06/b_reports/supp/03_module"
dir.create(RPT, recursive = TRUE, showWarnings = FALSE)

key_mods <- readLines(file.path(DAT, "key_modules.txt")) |>
  trimws() |> (\(x) x[nzchar(x)])()

hubs <- read_csv(file.path(DAT, "wgcna/wgcna_hub_proteins.csv"),
                 show_col_types = FALSE) %>%
  filter(module %in% key_mods) %>%
  group_by(module) %>%
  slice_max(kME, n = 8, with_ties = FALSE) %>%
  ungroup() %>%
  mutate(gene_label = ifelse(is.na(gene) | gene == "", uniprot_id, gene))

# Build adjacency: module sector → hub gene sector, weight = kME
link_df <- hubs %>%
  transmute(from = str_to_title(module),
            to   = gene_label,
            value = kME)

mod_colors <- setNames(unique(hubs$module), unique(hubs$module))
sectors    <- c(str_to_title(names(mod_colors)), link_df$to)
grid_cols  <- c(setNames(unname(mod_colors), str_to_title(names(mod_colors))),
                setNames(rep("grey80", nrow(link_df)), link_df$to))

pdf_path <- file.path(RPT, "SUPP_hub_chord.pdf")
png_path <- file.path(RPT, "SUPP_hub_chord.png")

draw_chord <- function() {
  circos.clear()
  circos.par(start.degree = 90, gap.degree = 2, track.margin = c(0.005, 0.005))
  chordDiagram(
    link_df,
    order            = sectors,
    grid.col         = grid_cols,
    transparency     = 0.3,
    directional      = 1,
    direction.type   = c("diffHeight"),
    diffHeight       = -0.04,
    annotationTrack  = "grid",
    preAllocateTracks = list(track.height = 0.08)
  )
  circos.trackPlotRegion(
    track.index = 1, panel.fun = function(x, y) {
      sector <- get.cell.meta.data("sector.index")
      xlim   <- get.cell.meta.data("xlim")
      ylim   <- get.cell.meta.data("ylim")
      is_mod <- sector %in% str_to_title(names(mod_colors))
      circos.text(mean(xlim), ylim[1] + (if (is_mod) 0.8 else 0.4),
                  sector,
                  facing    = "clockwise",
                  niceFacing = TRUE,
                  adj       = c(0, 0.5),
                  cex       = if (is_mod) 0.9 else 0.55,
                  font      = if (is_mod) 2 else 1,
                  col       = "grey15")
    }, bg.border = NA)
  title("Key-Module Hub Proteins (top 8 by kME per module)",
        cex.main = 0.95, line = -1)
}

pdf(pdf_path, width = 7.5, height = 7.5)
draw_chord(); dev.off()

png(png_path, width = 7.5, height = 7.5, units = "in", res = 300)
draw_chord(); dev.off()

message(sprintf("  SUPP hub chord saved (%d modules, %d hub links)",
                length(unique(link_df$from)), nrow(link_df)))
