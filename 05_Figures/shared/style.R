################################################################################
#   05_Figures — Shared Style
#   Self-contained colour palettes, sizing constants, theme, and helpers.
#   Source this file from any 05_Figures script that needs shared styling.
################################################################################

suppressPackageStartupMessages({
  library(ggplot2)
  library(scales)
  library(grid)
})

# ── Colour palettes ──────────────────────────────────────────────────────────

AGE_COLORS <- c(Young = "#4393C3", Old = "#D6604D")

DIR_COLORS <- c(Up = "#D6604D", Down = "#4393C3", NS = "grey70")

GROUP_FILL <- c(
  Young_Pre  = scales::alpha("#4393C3", 0.5),
  Young_Post = "#4393C3",
  Old_Pre    = scales::alpha("#D6604D", 0.5),
  Old_Post   = "#D6604D"
)

SHAPE_TP <- c(Pre = 16, Post = 17)

CONTRAST_COLORS <- c(
  Aging          = "#4CAF50",
  Training_Young = "#E05A4E",
  Training_Old   = "#5DA5DA",
  Interaction    = "#9B7FBF",
  Reversal       = "#FF8F00"
)

DB_COLORS <- c(Hallmark = "#AA336A", "GO:BP" = "#00796B")

PCA_COLORS <- c(
  Young_Pre  = "#93C4DE",
  Young_Post = "#2166AC",
  Old_Pre    = "#F4A582",
  Old_Post   = "#B2182B"
)

PCA_SHAPES <- c(Young_Pre = 16, Young_Post = 17, Old_Pre = 16, Old_Post = 17)

# ── Sizing constants ─────────────────────────────────────────────────────────

FIG_BASE_SIZE     <- 10
FIG_TITLE_SIZE    <- 12
FIG_SUBTITLE_SIZE <- 9
FIG_STRIP_SIZE    <- 10
FIG_AXIS_TITLE    <- 10
FIG_AXIS_TEXT     <- 8.5
FIG_LEGEND_TITLE  <- 9.5
FIG_LEGEND_TEXT   <- 8.5
FIG_TAG_SIZE      <- 15

# Key/legend styling constants (custom annotation-based legends)
KEY_TEXT      <- 2.8
KEY_TITLE     <- 3.2
KEY_ITEM      <- 0.40
KEY_BOX_HALF  <- 0.18
KEY_LW        <- 0.2

# ── Publication theme ────────────────────────────────────────────────────────

FIG_THEME <- theme_bw(base_size = FIG_BASE_SIZE) +
  theme(
    plot.title         = element_text(face = "bold", size = FIG_TITLE_SIZE),
    plot.subtitle      = element_text(face = "bold.italic", size = FIG_SUBTITLE_SIZE,
                                      color = "grey30"),
    plot.tag           = element_text(face = "bold", size = FIG_TAG_SIZE),
    strip.background   = element_blank(),
    strip.text         = element_text(face = "bold", size = FIG_STRIP_SIZE),
    axis.title         = element_text(face = "bold", size = FIG_AXIS_TITLE),
    axis.text          = element_text(size = FIG_AXIS_TEXT, color = "grey15"),
    legend.title       = element_text(face = "bold", size = FIG_LEGEND_TITLE,
                                      color = "grey20"),
    legend.text        = element_text(size = FIG_LEGEND_TEXT, color = "grey15"),
    legend.key.size    = unit(3, "mm"),
    panel.grid.minor   = element_blank()
  )

# ── Helper functions ─────────────────────────────────────────────────────────

sig_stars <- function(padj) {
  dplyr::case_when(
    padj < 0.001 ~ "***",
    padj < 0.01  ~ "**",
    padj < 0.05  ~ "*",
    TRUE         ~ ""
  )
}

clean_pathway_name <- function(name, max_chars = 45) {
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
    stringr::str_trunc(max_chars, ellipsis = "...")
}

darken_color <- function(col, factor = 0.7) {
  rgb_vals <- grDevices::col2rgb(col) / 255
  sapply(seq_along(col), function(i)
    grDevices::rgb(rgb_vals[1, i] * factor, rgb_vals[2, i] * factor,
                   rgb_vals[3, i] * factor))
}

strip_plot_meta <- function(p) {
  p + theme(plot.title = element_blank(), plot.subtitle = element_blank())
}

# ── Contrast label mappings ──────────────────────────────────────────────────

CTR_AXIS <- c(
  Aging          = "Aging",
  Training_Young = "Training\n(Young)",
  Training_Old   = "Training\n(Old)",
  Interaction    = "Interaction"
)

CTR_FACET <- c(
  Aging          = "Aging",
  Training_Young = "Training (Young)",
  Training_Old   = "Training (Old)",
  Interaction    = "Interaction"
)

CTR_SHORT <- c(
  Aging          = "Aging",
  Training_Young = "Tr. (Y)",
  Training_Old   = "Tr. (O)",
  Interaction    = "Inter."
)
