# F06 Supplementary — Hub Network Composite
# Dynamic: reads module list from key_modules.txt + mod_bio_labels.csv
# Layout: 2-row grid (ceil(n/3) rows x 3 cols) + networks composite in last slot
#
# Generates: 03_module/SUPP_hub_composite.{pdf,png}

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

# --- Read panels ---
read_panel <- function(file, dir = RPT_SRC) {
  path <- file.path(dir, file)
  if (!file.exists(path)) stop("Missing: ", path)
  rasterGrob(readPNG(path), interpolate = TRUE)
}

# Hub letters start after triptych letters (A-E for triptych, F+ for hub)
n_triptych <- length(KEY_MODULES)
hub_letters <- setNames(letters[(n_triptych + 1):(2 * n_triptych)], KEY_MODULES)
net_letter  <- letters[2 * n_triptych + 1]  # k for 5 modules

panels <- list()
for (mod in KEY_MODULES) {
  letter <- hub_letters[mod]
  slug <- pathway_slug[mod]
  if (is.na(slug)) slug <- mod
  fname <- sprintf("SUPP_hub_%s_%s.png", mod, slug)
  panels[[letter]] <- read_panel(fname)
}
# Networks composite
net_fname <- "SUPP_networks_composite.png"
panels[[net_letter]] <- read_panel(net_fname)

all_letters <- c(hub_letters, setNames(net_letter, "networks"))
n_total <- length(panels)

# --- Layout: 2 rows x 3 cols ---
n_cols <- 3
n_rows <- ceiling(n_total / n_cols)

# Build layout string
layout_rows <- character(n_rows)
letter_idx <- 1
all_panel_letters <- names(panels)
for (r in seq_len(n_rows)) {
  row_chars <- character(n_cols)
  for (c in seq_len(n_cols)) {
    if (letter_idx <= n_total) {
      row_chars[c] <- all_panel_letters[letter_idx]
      letter_idx <- letter_idx + 1
    } else {
      row_chars[c] <- "#"  # empty cell
    }
  }
  layout_rows[r] <- paste(row_chars, collapse = "")
}
# Add spacer rows between content rows
layout_parts <- character(0)
for (r in seq_along(layout_rows)) {
  layout_parts <- c(layout_parts, layout_rows[r])
  if (r < length(layout_rows)) layout_parts <- c(layout_parts, paste(rep("#", n_cols), collapse = ""))
}
layout_str <- paste(layout_parts, collapse = "\n")

# Heights: equal rows with spacers
SPACER <- 0.015
row_h <- (1 - SPACER * (n_rows - 1)) / n_rows
heights_vec <- numeric(0)
for (r in seq_len(n_rows)) {
  heights_vec <- c(heights_vec, row_h)
  if (r < n_rows) heights_vec <- c(heights_vec, SPACER)
}

# --- Compose ---
composite <- Reduce(`+`, lapply(panels, function(p) wrap_elements(full = p)))
composite <- composite +
  plot_layout(
    design  = layout_str,
    heights = heights_vec,
    widths  = rep(1, n_cols)
  ) +
  plot_annotation(
    title    = "Supplementary Figure S-Hub: Hub Protein Networks",
    subtitle = paste0(length(KEY_MODULES), " key modules + composite \u2014 node color = gene significance; magenta (translation initiation) shows few hubs because TOM co-expression within the module is dominated by one connected component anchored at EIF2S3"),
    theme = theme(
      plot.title    = element_text(face = "bold", size = 7, hjust = 0),
      plot.subtitle = element_text(size = 4, face = "bold.italic", color = "grey40", hjust = 0),
      plot.margin   = margin(2, 4, 2, 4)
    )
  )

# --- Tag placement via cowplot ---
TAG_SZ <- 8
# X positions for 3 columns
x_pos <- (seq_len(n_cols) - 1) / n_cols + 0.02
# Y positions for each row top
header_frac <- 0.03
y_pos <- numeric(n_rows)
for (r in seq_len(n_rows)) {
  frac_above <- (r - 1) * (row_h + SPACER) / sum(heights_vec)
  y_pos[r] <- 1 - header_frac - frac_above * (1 - header_frac)
}

composite <- ggdraw(composite)
panel_idx <- 1
for (r in seq_len(n_rows)) {
  for (c in seq_len(n_cols)) {
    if (panel_idx <= n_total) {
      composite <- composite +
        draw_label(all_panel_letters[panel_idx],
                   x = x_pos[c], y = y_pos[r],
                   size = TAG_SZ, fontface = "bold", hjust = 0, vjust = 1)
      panel_idx <- panel_idx + 1
    }
  }
}

COMP_W <- 178
COMP_H <- 120

pdf_device <- get_pdf_device()

ggsave(file.path(RPT_PDF, "SUPP_hub_composite.pdf"), composite,
       width = COMP_W, height = COMP_H, units = "mm",
       device = pdf_device, limitsize = FALSE)
ggsave(file.path(RPT_PNG, "SUPP_hub_composite.png"), composite,
       width = COMP_W, height = COMP_H, units = "mm",
       dpi = 300, limitsize = FALSE)

message("F06 supp module hub composite saved (",
        paste(KEY_MODULES, collapse = ", "), " + networks)")
