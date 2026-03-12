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
THEME_PUB <- ggplot2::theme_bw(base_size = 10) +
  ggplot2::theme(
    plot.title       = ggplot2::element_text(face = "bold", size = 12),
    plot.subtitle    = ggplot2::element_text(size = 9, color = "grey30",
                                             face = "bold.italic"),
    strip.background = ggplot2::element_blank(),
    strip.text       = ggplot2::element_text(face = "bold", size = 10),
    legend.key.size  = grid::unit(3, "mm")
  )

# Figure-level theme for F2/F3 composite panels (larger base for dense figures)
THEME_FIG <- ggplot2::theme_bw(base_size = 14) +
  ggplot2::theme(
    plot.title       = ggplot2::element_text(face = "bold", size = 16),
    plot.subtitle    = ggplot2::element_text(size = 11, color = "grey30",
                                             face = "bold.italic"),
    axis.title       = ggplot2::element_text(face = "bold", size = 14),
    axis.text        = ggplot2::element_text(size = 12),
    strip.background = ggplot2::element_blank(),
    strip.text       = ggplot2::element_text(face = "bold", size = 12),
    legend.text      = ggplot2::element_text(size = 12),
    legend.title     = ggplot2::element_text(size = 14, face = "bold"),
    legend.key.size  = grid::unit(4, "mm"),
    plot.margin      = ggplot2::margin(2, 2, 2, 2, "mm")
  )

# Unified annotation text constants for F2/F3 panels (geom_text/annotate size in mm)
# Fontface conventions:
#   genes = italic | pathways = bold | quadrant labels = bold
#   statistics = italic | direction labels = plain | headers = bold
TXT_TAG   <- 18    # panel letter tags (A, B, C...)
TXT_LABEL <- 3.5   # all annotation/label text (genes, pathways, quadrants, stats, ORA counts)

# Panel F text: compensate for wrap_elements() shrinkage in composite
# Standalone ~280mm → composite ~215mm = 0.77x scale; moderate boost sufficient
TXT_PF       <- 5.0     # Panel F annotation text (pathways, counts, key labels)
TXT_PF_GENE  <- 16      # Panel F gene symbol axis text (theme element, in pt)

# Legacy aliases (backwards compatibility with any scripts still referencing old names)
TXT_PATHWAY   <- TXT_LABEL
TXT_GENE      <- TXT_LABEL
TXT_QUADRANT  <- TXT_LABEL
TXT_STATS     <- TXT_LABEL
TXT_ORA_BAR   <- TXT_LABEL
TXT_ORA_AXIS  <- TXT_LABEL
TXT_ORA_STRIP <- TXT_LABEL

# Key/legend styling constants (F1 reference design)
KEY_TEXT      <- 2.8   # annotate() size for tip text in custom keys
KEY_TITLE     <- 3.5   # annotate() size for section headers
KEY_ITEM      <- 3.5   # geom_text size for item labels in swatch keys
KEY_BOX_HALF  <- 0.18  # half-height of swatch boxes
KEY_HDR_COL   <- "grey25"
KEY_ITEM_COL  <- "grey15"
KEY_BORDER    <- "grey70"
KEY_LW        <- 0.2

# Built-in ggplot legend theme (for panels that use guide/legend, not custom keys)
LEGEND_THEME <- ggplot2::theme(
  legend.text      = ggplot2::element_text(size = 9, color = "grey15"),
  legend.title     = ggplot2::element_text(size = 10, face = "bold", color = "grey25"),
  legend.key.size  = grid::unit(3, "mm"),
  legend.position  = "bottom",
  legend.box       = "horizontal",
  legend.margin    = ggplot2::margin(0, 0, 0, 0)
)

COLORBAR_GUIDE <- ggplot2::guide_colorbar(
  barwidth  = grid::unit(15, "mm"),
  barheight = grid::unit(2, "mm"),
  title.position = "left",
  title.hjust = 1
)

# Shared helpers
sig_stars <- function(padj) {
  dplyr::case_when(padj < 0.001 ~ "***",
                   padj < 0.01  ~ "**",
                   padj < 0.05  ~ "*",
                   TRUE         ~ "")
}

clean_pathway_name <- function(name) {
  name |>
    stringr::str_remove("^HALLMARK_") |>
    stringr::str_remove("^GOBP_") |>
    stringr::str_remove("^GOCC_") |>
    stringr::str_remove("^GOMF_") |>
    stringr::str_remove("^REACTOME_") |>
    stringr::str_remove("^KEGG_MEDICUS_") |>
    stringr::str_remove("^KEGG_") |>
    stringr::str_remove("^WP_") |>
    stringr::str_replace_all("_", " ") |>
    stringr::str_to_title() |>
    stringr::str_replace("Mtorc1", "mTORC1") |>
    stringr::str_replace("Myc ", "MYC ") |>
    stringr::str_replace("E2f ", "E2F ") |>
    stringr::str_replace("Dna ", "DNA ") |>
    stringr::str_replace("Rna ", "RNA ") |>
    stringr::str_replace("Tnfa ", "TNFa ") |>
    stringr::str_replace("Uv ", "UV ") |>
    stringr::str_replace("G2m ", "G2M ") |>
    stringr::str_replace("Il6 ", "IL6 ") |>
    stringr::str_replace("Il2 ", "IL2 ") |>
    stringr::str_replace("Kras ", "KRAS ") |>
    stringr::str_replace("P53 ", "p53 ") |>
    stringr::str_replace("Tgf ", "TGF ") |>
    stringr::str_replace("Nf Kb", "NF-kB") |>
    stringr::str_replace("Atp ", "ATP ") |>
    stringr::str_replace("Nadh ", "NADH ") |>
    stringr::str_replace("Oxidative Phosphorylation", "OXPHOS") |>
    stringr::str_replace("External Encapsulating Structure Or.*",
                         "Extracellular Matrix Organization") |>
    stringr::str_replace("Enzyme Linked Receptor Protein Signaling.*",
                         "Receptor Protein Signaling") |>
    stringr::str_trunc(45, ellipsis = "...")
}

darken_color <- function(col, factor = 0.5) {
  rgb_vals <- grDevices::col2rgb(col) / 255
  sapply(seq_along(col), function(i)
    grDevices::rgb(rgb_vals[1, i] * factor, rgb_vals[2, i] * factor,
                   rgb_vals[3, i] * factor))
}
