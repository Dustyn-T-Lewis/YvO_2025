# Sourced by 02_supp_panels.R — expects style.R already loaded.

library(readr)
library(dplyr)
library(circlize)
library(stringr)

BASE <- "04_Figures/F06"

DAT     <- file.path(BASE, "c_data")
RPT_PNG <- file.path(BASE, "b_reports", "supp", "03_module", "png", "panels")
RPT_PDF <- file.path(BASE, "b_reports", "supp", "03_module", "pdf", "panels")
dir.create(RPT_PNG, recursive = TRUE, showWarnings = FALSE)
dir.create(RPT_PDF, recursive = TRUE, showWarnings = FALSE)

key_mods <- readLines(file.path(DAT, "key_modules.txt")) |>
  trimws() |> (\(x) x[nzchar(x)])()

hubs <- read_csv(file.path(DAT, "wgcna/wgcna_hub_proteins.csv")) |>
  filter(module %in% key_mods) |>
  group_by(module) |>
  slice_max(kME, n = 8, with_ties = FALSE) |>
  ungroup() |>
  mutate(gene_label = ifelse(is.na(gene) | gene == "", uniprot_id, gene))

link_df <- hubs |>
  transmute(from = str_to_title(module),
            to   = gene_label,
            value = kME)

mod_colors <- setNames(unique(hubs$module), unique(hubs$module))
sectors    <- c(str_to_title(names(mod_colors)), link_df$to)
grid_cols  <- c(setNames(unname(mod_colors), str_to_title(names(mod_colors))),
                setNames(rep("grey80", nrow(link_df)), link_df$to))

pdf_path <- file.path(RPT_PDF, "SUPP_hub_chord.pdf")
png_path <- file.path(RPT_PNG, "SUPP_hub_chord.png")

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
