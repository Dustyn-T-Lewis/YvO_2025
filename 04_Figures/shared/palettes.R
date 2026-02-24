################################################################################
#   YvO Figures — Shared Palettes
#   Canonical colour definitions for all figures and pipeline QC plots.
#   Source this file from any script that needs group/age/direction colours.
################################################################################

AGE_COLORS <- c(Young = "#4393C3", Old = "#D6604D")
DIR_COLORS <- c(Up = "#D6604D", Down = "#4393C3", NS = "grey70")

GROUP_FILL <- c(
  Young_Pre  = scales::alpha("#4393C3", 0.5),
  Young_Post = "#4393C3",
  Old_Pre    = scales::alpha("#D6604D", 0.5),
  Old_Post   = "#D6604D"
)

# Shape mappings
SHAPE_TP <- c(Pre = 16, Post = 17)

# Publication theme (shared across F1-F6)
THEME_PUB <- ggplot2::theme_bw(base_size = 8) +
  ggplot2::theme(
    plot.title       = ggplot2::element_text(face = "bold", size = 9),
    plot.subtitle    = ggplot2::element_text(size = 6.5, color = "grey30",
                                             face = "italic"),
    strip.background = ggplot2::element_blank(),
    strip.text       = ggplot2::element_text(face = "bold", size = 6.5),
    legend.key.size  = grid::unit(3, "mm")
  )

# Shared helpers
sig_stars <- function(padj) {
  dplyr::case_when(padj < 0.001 ~ "***",
                   padj < 0.01  ~ "**",
                   padj < 0.05  ~ "*",
                   TRUE         ~ "")
}
