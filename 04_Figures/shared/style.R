# 04_Figures — Unified Style
# Single source of truth: palettes, themes, sizing constants, helpers.

suppressPackageStartupMessages({
  library(ggplot2)
  library(scales)
  library(grid)
})

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

DB_COLORS <- c(Hallmark = "#AA336A", KEGG = "#E65100", Reactome = "#1565C0",
               WikiPathways = "#6A1B9A", "GO:BP" = "#00796B",
               BioCarta = "#795548", PID = "#455A64")

PCA_COLORS <- c(
  Young_Pre  = "#93C4DE",
  Young_Post = "#2166AC",
  Old_Pre    = "#F4A582",
  Old_Post   = "#B2182B"
)

PCA_SHAPES <- c(Young_Pre = 16, Young_Post = 17, Old_Pre = 16, Old_Post = 17)

SIG_COLORS_F2 <- c(
  "Interaction"    = "#7B5EA7",
  "Sig Both"       = "#2E7D32",
  "Sig Young only" = "#E05A4E",
  "Sig Old only"   = "#5DA5DA",
  "NS"             = "grey70"
)

SIG_LABEL_FILL_F2 <- c(
  "Interaction"    = scales::alpha("#7B5EA7", 0.75),
  "Sig Both"       = scales::alpha("#2E7D32", 0.75),
  "Sig Young only" = scales::alpha("#E05A4E", 0.75),
  "Sig Old only"   = scales::alpha("#5DA5DA", 0.75),
  "NS"             = scales::alpha("grey70",  0.75)
)
SIG_LABEL_TEXT_F2 <- c(
  "Interaction"    = "white",
  "Sig Both"       = "white",
  "Sig Young only" = "white",
  "Sig Old only"   = "white",
  "NS"             = "white"
)

SIG_COLORS_F3 <- c(
  "Reversal"           = "#7B5EA7",
  "Sig Both"           = "#2E7D32",
  "Sig Aging only"     = "#E05A4E",
  "Sig Training only"  = "#5DA5DA",
  "NS"                 = "grey70"
)

SIG_LABEL_FILL_F3 <- c(
  "Reversal"           = scales::alpha("#7B5EA7", 0.75),
  "Sig Both"           = scales::alpha("#2E7D32", 0.75),
  "Sig Aging only"     = scales::alpha("#E05A4E", 0.75),
  "Sig Training only"  = scales::alpha("#5DA5DA", 0.75),
  "NS"                 = scales::alpha("grey70",  0.75)
)
SIG_LABEL_TEXT_F3 <- c(
  "Reversal"           = "white",
  "Sig Both"           = "white",
  "Sig Aging only"     = "white",
  "Sig Training only"  = "white",
  "NS"                 = "white"
)

EXPANDED_ORDER <- c("Interaction",
                    "Sig Both Up", "Sig Both Down",
                    "Sig Young Up", "Sig Young Down",
                    "Sig Old Up", "Sig Old Down")

EXPANDED_COLORS <- c(
  "Interaction"    = "#7B5EA7",
  "Sig Both Up"    = "#2E7D32",
  "Sig Both Down"  = "#1B5E20",
  "Sig Young Up"   = "#E05A4E",
  "Sig Young Down" = "#B71C1C",
  "Sig Old Up"     = "#5DA5DA",
  "Sig Old Down"   = "#1565C0"
)

ORA_QUAD_COLORS_F2 <- c(
  "Concordant Up"               = "#E57373",
  "Concordant Down"             = "#64B5F6",
  "Discordant (Y Up / O Down)"  = "#FFB74D",
  "Discordant (Y Down / O Up)"  = "#81C784"
)

ORA_QUAD_COLORS_F3 <- c(
  "Reversed (Aging Up / Training Down)"  = "#E57373",
  "Reversed (Aging Down / Training Up)"  = "#64B5F6",
  "Exacerbated Up"                       = "#FFB74D",
  "Exacerbated Down"                     = "#81C784"
)

CLUSTER_COLORS <- c(C1 = "#E74C3C", C2 = "#3498DB", C3 = "#2ECC71",
                     C4 = "#F39C12", C5 = "#9B59B6", C6 = "#1ABC9C",
                     C7 = "#E67E22", C8 = "#34495E", C9 = "#D35400",
                     C10 = "#7F8C8D")

THEME_COLORS <- c(
  "Mitochondrial & Energy Metabolism"    = "#E57373",
  "Muscle Structure & Myogenesis"        = "#64B5F6",
  "Proteostasis & Stress Response"       = "#81C784",
  "Cytoskeletal & Cell Division"         = "#CE93D8",
  "Immune & Complement"                  = "#FFB74D",
  "ECM & Tissue Remodeling"             = "#F06292",
  "Metabolic & Redox Regulation"         = "#FFD54F",
  "Intracellular Transport & Signaling"  = "#4DB6AC"
)


PANEL_SM  <- 120   # small panels (e.g., PCA, bar charts)
PANEL_MD  <- 180   # medium panels (e.g., volcano rings, scatter)
PANEL_LG  <- 280   # large panels (e.g., classification triptych, heatmaps)

# Base annotation sizes calibrated for PANEL_MD (180mm).
# Use scale_text() to adjust for panels of different widths.
BASE_PATHWAY  <- 4.0   # pathway/term labels inside figures
BASE_GENE     <- 3.2   # gene name labels
BASE_STAT     <- 3.5   # statistical annotations (r=, p=, etc.)
BASE_QUADRANT <- 4.0   # quadrant labels (RRHO, scatter)
BASE_COUNT    <- 3.5   # count labels on bars
BASE_TAG      <- 18    # panel letter tags (A, B, C...)

scale_text <- function(base_size, panel_width_mm, ref_width = PANEL_MD) {
  base_size * sqrt(panel_width_mm / ref_width)
}


FIG_BASE_SIZE     <- 10
FIG_TITLE_SIZE    <- 12
FIG_SUBTITLE_SIZE <- 9
FIG_STRIP_SIZE    <- 10
FIG_AXIS_TITLE    <- 10
FIG_AXIS_TEXT     <- 8.5
FIG_LEGEND_TITLE  <- 9.5
FIG_LEGEND_TEXT   <- 8.5
FIG_TAG_SIZE      <- 15

# Standard panels (base_size = 10)
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

# Dense/large panels (base_size = 14): RRHO, classification, heatmaps
THEME_FIG <- theme_bw(base_size = 14) +
  theme(
    plot.title       = element_text(face = "bold", size = 16),
    plot.subtitle    = element_text(size = 11, color = "grey30",
                                    face = "bold.italic"),
    axis.title       = element_text(face = "bold", size = 14),
    axis.text        = element_text(size = 12),
    strip.background = element_blank(),
    strip.text       = element_text(face = "bold", size = 12),
    legend.text      = element_text(size = 12),
    legend.title     = element_text(size = 14, face = "bold"),
    legend.key.size  = unit(4, "mm"),
    plot.margin      = margin(2, 2, 2, 2, "mm")
  )

# Legacy alias for scripts referencing old name
THEME_PUB <- FIG_THEME


LEGEND_THEME <- theme(
  legend.text      = element_text(size = 9, color = "grey15"),
  legend.title     = element_text(size = 10, face = "bold", color = "grey25"),
  legend.key.size  = unit(3, "mm"),
  legend.position  = "bottom",
  legend.box       = "horizontal",
  legend.margin    = margin(0, 0, 0, 0)
)

COLORBAR_GUIDE <- guide_colorbar(
  barwidth  = unit(15, "mm"),
  barheight = unit(2, "mm"),
  title.position = "left",
  title.hjust = 1
)


KEY_TEXT      <- 2.8
KEY_TITLE     <- 3.2
KEY_ITEM      <- 0.40
KEY_BOX_HALF  <- 0.18
KEY_LW        <- 0.2
KEY_HDR_COL   <- "grey25"
KEY_ITEM_COL  <- "grey15"
KEY_BORDER    <- "grey70"


fmt_p <- function(p) {
  if (p < 0.001) return("p < 0.001")
  if (p < 0.01)  return(sprintf("p = %.3f", p))
  sprintf("p = %.2f", p)
}

classify_proteins_f2 <- function(pi_Y, pi_O, pi_int, threshold = 0.05) {
  dplyr::case_when(
    pi_int < threshold                    ~ "Interaction",
    pi_Y < threshold & pi_O < threshold   ~ "Sig Both",
    pi_Y < threshold                      ~ "Sig Young only",
    pi_O < threshold                      ~ "Sig Old only",
    TRUE                                  ~ "NS"
  ) |>
    factor(levels = c("Interaction", "Sig Both",
                       "Sig Young only", "Sig Old only", "NS"))
}

classify_expanded_f2 <- function(pi_Y, pi_O, pi_int, logFC_Y, logFC_O,
                                  threshold = 0.05) {
  base <- classify_proteins_f2(pi_Y, pi_O, pi_int, threshold)
  dplyr::case_when(
    base == "Interaction"      ~ "Interaction",
    base == "Sig Both"   & logFC_Y > 0  ~ "Sig Both Up",
    base == "Sig Both"   & logFC_Y <= 0 ~ "Sig Both Down",
    base == "Sig Young only" & logFC_Y > 0  ~ "Sig Young Up",
    base == "Sig Young only" & logFC_Y <= 0 ~ "Sig Young Down",
    base == "Sig Old only"   & logFC_O > 0  ~ "Sig Old Up",
    base == "Sig Old only"   & logFC_O <= 0 ~ "Sig Old Down",
    TRUE ~ NA_character_
  ) |>
    factor(levels = EXPANDED_ORDER)
}

classify_proteins_f3 <- function(pi_aging, pi_training_old, pi_reversal,
                                  threshold = 0.05) {
  dplyr::case_when(
    pi_reversal < threshold                            ~ "Reversal",
    pi_aging < threshold & pi_training_old < threshold ~ "Sig Both",
    pi_aging < threshold                               ~ "Sig Aging only",
    pi_training_old < threshold                        ~ "Sig Training only",
    TRUE                                               ~ "NS"
  ) |>
    factor(levels = c("Reversal", "Sig Both",
                       "Sig Aging only", "Sig Training only", "NS"))
}

# Fisher z-transform CI for Pearson or partial r
# For partial r, effective n = n - k (k = number of covariates)
# Reference: Bonett & Wright 2000, Psychometrika 65:23-28
fisher_z_ci <- function(r, n, k = 0, level = 0.95) {
  n_eff <- n - k
  if (n_eff < 4 || is.na(r)) return(c(lo = NA_real_, hi = NA_real_))
  z <- atanh(r)
  se <- 1 / sqrt(n_eff - 3)
  crit <- qnorm(1 - (1 - level) / 2)
  c(lo = tanh(z - crit * se), hi = tanh(z + crit * se))
}

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

make_sigmoid_ribbon <- function(x0, x1, y0_top, y0_bot, y1_top, y1_bot,
                                n_pts = 50, ribbon_id) {
  t <- seq(0, 1, length.out = n_pts)
  blend <- (1 - cos(pi * t)) / 2
  tibble::tibble(
    x = c(x0 + (x1 - x0) * t, rev(x0 + (x1 - x0) * t)),
    y = c(y0_top + (y1_top - y0_top) * blend,
          rev(y0_bot + (y1_bot - y0_bot) * blend)),
    ribbon_id = ribbon_id
  )
}

reorder_within <- function(x, by, within, fun = mean, sep = "___", ...) {
  new_x <- paste(x, within, sep = sep)
  stats::reorder(new_x, by, FUN = fun, ...)
}

scale_y_reordered <- function(..., sep = "___") {
  reg <- paste0(sep, ".+$")
  ggplot2::scale_y_discrete(labels = function(x) gsub(reg, "", x), ...)
}

assign_theme <- function(pathway_name) {
  pw <- tolower(pathway_name)
  dplyr::case_when(
    stringr::str_detect(pw, "mitochon|oxidative.phosph|respiratory|electron.transport|tca|citrate|nadh|atp|fatty.acid|lipid|adipogen|acetyl.coa|amide.metabol") ~
      "Mitochondrial & Energy Metabolism",
    stringr::str_detect(pw, "myogen|muscle|contract|myofib|sarco|neuromuscul") ~
      "Muscle Structure & Myogenesis",
    stringr::str_detect(pw, "mtorc|unfold|chaper|heat.shock|protein.stabili|proteasom|ubiquitin|protein.fold|apoptosis|programmed.cell.death|cell.death") ~
      "Proteostasis & Stress Response",
    stringr::str_detect(pw, "microtub|spindle|mitotic|cell.divis|cell.cycle|cytoskelet|tubulin|actin") ~
      "Cytoskeletal & Cell Division",
    stringr::str_detect(pw, "immun|inflam|complement|cytokine|interferon|heme|blood|coagulat") ~
      "Immune & Complement",
    stringr::str_detect(pw, "extracellular|matrix|collagen|adhesion|integrin|mesenchym|epithelial") ~
      "ECM & Tissue Remodeling",
    stringr::str_detect(pw, "glycol|metabol|xenobiot|aldehyde|oxidant|detox|pyridine|reactive.oxygen|peroxide") ~
      "Metabolic & Redox Regulation",
    stringr::str_detect(pw, "vesicle|transport|endosom|golgi|lysosom|signal|kinase|androgen") ~
      "Intracellular Transport & Signaling",
    TRUE ~ "Other"
  )
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


get_pdf_device <- function() {
  tryCatch(
    { cairo_pdf(tempfile()); dev.off(); cairo_pdf },
    error = function(e) "pdf"
  )
}
