# F06 Supplementary — Module Triptych Composite
# Dynamic: reads module list from key_modules.txt + mod_bio_labels.csv
# Layout: vertical stack with protein-count-weighted heights
#
# Generates: 03_module/SUPP_triptych_composite.{pdf,png}

setwd(rprojroot::find_rstudio_root_file())
source("04_Figures/shared/style.R")

suppressPackageStartupMessages({
  library(patchwork)
  library(cowplot)
  library(png)
  library(grid)
  library(readr)
  library(dplyr)
  library(stringr)
})

RPT_SRC <- "04_Figures/F06/b_reports/supp/03_module/png/panels"
RPT_PDF <- "04_Figures/F06/b_reports/supp/03_module/pdf"
RPT_PNG <- "04_Figures/F06/b_reports/supp/03_module/png"
DAT     <- "04_Figures/F06/c_data"
dir.create(RPT_PDF, recursive = TRUE, showWarnings = FALSE)
dir.create(RPT_PNG, recursive = TRUE, showWarnings = FALSE)

# --- Dynamic module list ---
km_file <- file.path(DAT, "key_modules.txt")
KEY_MODULES <- if (file.exists(km_file)) {
  km <- readLines(km_file)
  km[nzchar(trimws(km))]
} else c("turquoise", "brown", "magenta", "blue", "green")

mod_bio <- read_csv(file.path(DAT, "mod_bio_labels.csv"), show_col_types = FALSE)
pathway_slug <- setNames(
  gsub("[/ ]+", "_", tolower(mod_bio$bio_label)),
  mod_bio$module_color
)

# Module protein counts for weighted heights
mod_assign <- read_csv(file.path(DAT, "wgcna/wgcna_module_assignments.csv"), show_col_types = FALSE)
mod_counts <- mod_assign %>%
  filter(module_color %in% KEY_MODULES) %>%
  count(module_color) %>%
  mutate(module_color = factor(module_color, levels = KEY_MODULES)) %>%
  arrange(module_color)

# --- Read panels ---
read_panel <- function(file, dir = RPT_SRC) {
  path <- file.path(dir, file)
  if (!file.exists(path)) stop("Missing: ", path)
  rasterGrob(readPNG(path), interpolate = TRUE)
}

supp_letters <- setNames(LETTERS[seq_along(KEY_MODULES)], KEY_MODULES)
panels <- list()
for (mod in KEY_MODULES) {
  letter <- supp_letters[mod]
  slug <- pathway_slug[mod]
  if (is.na(slug)) slug <- mod
  fname <- sprintf("SUPP_triptych_%s_%s.png", mod, slug)
  panels[[mod]] <- read_panel(fname)
}

n_panels <- length(KEY_MODULES)

# --- Weighted heights based on sqrt(protein count) ---
# Floor at 0.5 so even small modules get readable space
counts <- mod_counts$n[match(KEY_MODULES, mod_counts$module_color)]
weights <- pmax(sqrt(counts), 0.5)
weights <- weights / sum(weights)

# Build heights vector: panel, spacer, panel, spacer, ...
SPACER <- 0.008
n_spacers <- n_panels - 1
total_spacer <- n_spacers * SPACER
panel_heights <- weights * (1 - total_spacer)

heights_vec <- numeric(0)
for (i in seq_along(panel_heights)) {
  heights_vec <- c(heights_vec, panel_heights[i])
  if (i < n_panels) heights_vec <- c(heights_vec, SPACER)
}

# Build layout string: A, #, B, #, C, ...
layout_chars <- LETTERS[1:n_panels]
layout_parts <- character(0)
for (i in seq_along(layout_chars)) {
  layout_parts <- c(layout_parts, layout_chars[i])
  if (i < n_panels) layout_parts <- c(layout_parts, "#")
}
layout_str <- paste(layout_parts, collapse = "\n")

# --- Compose ---
composite <- Reduce(`+`, lapply(panels, function(p) wrap_elements(full = p)))
composite <- composite +
  plot_layout(
    design  = layout_str,
    heights = heights_vec,
    widths  = c(1)
  ) +
  plot_annotation(
    title    = "Supplementary Figure S-Triptych: Per-Module Deep Dives",
    subtitle = "Heatmap | eigengene dynamics | pathway enrichment \u2014 key modules",
    theme = theme(
      plot.title    = element_text(face = "bold", size = 11, hjust = 0),
      plot.subtitle = element_text(size = 8, face = "bold.italic", color = "grey40", hjust = 0),
      plot.margin   = margin(2, 4, 2, 4)
    )
  )

# --- Tag placement via cowplot ---
TAG_SZ <- 12
X_LEFT <- 0.02

# Calculate Y positions from cumulative heights (top of each panel)
cum_h <- cumsum(heights_vec)
total_h <- sum(heights_vec)
# Y = 1 - (cumulative height before panel start) / total, adjusted for annotation header
header_frac <- 0.03  # approximate space for title/subtitle
y_positions <- numeric(n_panels)
y_positions[1] <- 1 - header_frac
for (i in 2:n_panels) {
  # Sum of all heights before panel i (all prior panels + spacers)
  idx_before <- 2 * (i - 1)  # heights_vec indices: 1, (2=spacer), 3, (4=spacer), ...
  frac_before <- sum(heights_vec[1:idx_before]) / total_h
  y_positions[i] <- 1 - header_frac - frac_before * (1 - header_frac)
}

composite <- ggdraw(composite)
for (i in seq_along(KEY_MODULES)) {
  composite <- composite +
    draw_label(supp_letters[KEY_MODULES[i]],
               x = X_LEFT, y = y_positions[i],
               size = TAG_SZ, fontface = "bold", hjust = 0, vjust = 1)
}

COMP_W <- 300
COMP_H <- 550

pdf_device <- get_pdf_device()

ggsave(file.path(RPT_PDF, "SUPP_triptych_composite.pdf"), composite,
       width = COMP_W, height = COMP_H, units = "mm",
       device = pdf_device, limitsize = FALSE)
ggsave(file.path(RPT_PNG, "SUPP_triptych_composite.png"), composite,
       width = COMP_W, height = COMP_H, units = "mm",
       dpi = 300, limitsize = FALSE)

message("F06 supp module triptych composite saved (",
        paste(KEY_MODULES, collapse = ", "), ")")
