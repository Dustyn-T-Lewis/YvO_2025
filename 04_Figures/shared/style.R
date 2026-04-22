# 04_Figures — Unified Style
# Single source of truth: palettes, themes, sizing constants, helpers.

library(ggplot2)
library(scales)
library(grid)

# Headless-script guard: redirect R's default graphics device to a null pdf
# so `Rplots.pdf` is never auto-created in the working dir on implicit plot
# calls (circlize, ComplexHeatmap, cowplot, and any cairo_pdf probe fallback).
# Works even if a downstream dev.off() chain closes whatever device is open:
# the NEXT implicit open also routes to nullfile(). ggsave() and explicit
# pdf("file.pdf") / cairo_pdf("file.pdf") are unaffected — they manage their
# own devices. Safe in interactive sessions too (RStudio uses its own device).
options(device = function(...) grDevices::pdf(file = nullfile(), ...))

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
  Interaction    = "#9B7FBF"
)

DB_COLORS <- c(Hallmark = "#AA336A", KEGG = "#E65100", Reactome = "#1565C0",
               "GO:BP" = "#00796B",
               "GO Slim" = "#8D6E63", BioCarta = "#795548", PID = "#455A64")

PAL_CLASS <- c(Complete = "#4DAF4A", MAR = "#377EB8", MNAR = "#E41A1C")

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
  "Sig Both"           = "#2E7D32",
  "Sig Aging only"     = "#E05A4E",
  "Sig Training only"  = "#5DA5DA",
  "NS"                 = "grey70"
)

SIG_LABEL_FILL_F3 <- c(
  "Sig Both"           = scales::alpha("#2E7D32", 0.75),
  "Sig Aging only"     = scales::alpha("#E05A4E", 0.75),
  "Sig Training only"  = scales::alpha("#5DA5DA", 0.75),
  "NS"                 = scales::alpha("grey70",  0.75)
)
SIG_LABEL_TEXT_F3 <- c(
  "Sig Both"           = "white",
  "Sig Aging only"     = "white",
  "Sig Training only"  = "white",
  "NS"                 = "white"
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


PANEL_MD  <- 178   # J Physiol double-column width (reference for scale_text)

# Base annotation sizes (mm) calibrated for PANEL_MD (178mm) at T4 (7pt).
# Use scale_text() to adjust for panels of different widths.
# Conversion: 7pt / 2.835 ≈ 2.5mm; 8pt / 2.835 ≈ 2.8mm
BASE_PATHWAY  <- 2.8   # pathway/term labels inside figures (~8pt)
BASE_GENE     <- 2.5   # gene name labels (~7pt)
BASE_STAT     <- 2.5   # statistical annotations (r=, p=, etc.) (~7pt)
BASE_QUADRANT <- 2.8   # quadrant labels (RRHO, scatter) (~8pt)
BASE_COUNT    <- 2.5   # count labels on bars (~7pt)
BASE_TAG      <- 8     # panel letter tags (a, b, c...)

scale_text <- function(base_size, panel_width_mm, ref_width = PANEL_MD) {
  base_size * sqrt(panel_width_mm / ref_width)
}

# Strip titles, subtitles, tags, and legends for composite stitching.
# Call AFTER standalone ggsave so standalone PNGs keep decorations.
strip_for_composite <- function(p) {
  p + labs(title = NULL, subtitle = NULL, tag = NULL) +
    theme(legend.position = "none")
}

# Luminance-based check for light module colors (text contrast)
is_light_color <- function(color_name) {
  rgb_val <- col2rgb(color_name)
  (0.299 * rgb_val[1] + 0.587 * rgb_val[2] + 0.114 * rgb_val[3]) / 255 > 0.6
}


# Journal-normalised text hierarchy (J Physiol / Wiley)
# T1 = 7pt   panel tags, panel titles
# T2 = 5pt   axis titles, axis text, strip text, legend titles
# T3 = 4pt   subtitles, legend text, in-plot annotations
FIG_TITLE_SIZE    <- 7      # T1
FIG_SUBTITLE_SIZE <- 4      # T3
FIG_STRIP_SIZE    <- 5      # T2
FIG_AXIS_TEXT     <- 5      # T2
FIG_LEGEND_TITLE  <- 5      # T2
FIG_LEGEND_TEXT   <- 4      # T3

# Journal-calibrated composite text sizes.
# J Physiol spec: 10–12 pt sans-serif at reproduction size; 12 pt for tags.
# Composites range ~60–240 mm tall at journal widths (85/178 mm).
# Gentle height scaling keeps text proportional; clamped to spec range.
composite_text_sizes <- function(comp_h_mm) {
  list(
    title    = pmax(6, pmin(8, round(5 + comp_h_mm / 80))),
    subtitle = pmax(4, pmin(6, round(3 + comp_h_mm / 100))),
    tag      = 8
  )
}

# Standard panels — Helvetica sans-serif per J Physiol spec
FIG_THEME <- theme_bw(base_size = 6, base_family = "Helvetica") +
  theme(
    plot.title         = element_text(face = "bold", size = FIG_TITLE_SIZE,
                                      margin = margin(b = 1)),
    plot.subtitle      = element_text(face = "bold.italic", size = FIG_SUBTITLE_SIZE,
                                      color = "grey30", margin = margin(t = 0, b = 2)),
    plot.tag           = element_text(face = "bold", size = BASE_TAG),
    strip.background   = element_blank(),
    strip.text         = element_text(face = "bold", size = FIG_STRIP_SIZE),
    axis.title.x       = element_text(face = "bold", size = 5,
                                      margin = margin(t = 0)),
    axis.title.y       = element_text(face = "bold", size = 5,
                                      margin = margin(r = -1)),
    axis.text          = element_text(size = FIG_AXIS_TEXT, color = "grey15"),
    legend.title       = element_text(face = "bold", size = FIG_LEGEND_TITLE,
                                      color = "grey20"),
    legend.text        = element_text(size = FIG_LEGEND_TEXT, color = "grey15"),
    legend.key.size    = unit(2.5, "mm"),
    panel.grid.minor   = element_blank()
  )


fmt_p <- function(p) {
  if (p < 0.001) return("p < 0.001")
  if (p < 0.01)  return(sprintf("p = %.3f", p))
  sprintf("p = %.2f", p)
}

# Plotmath-parseable p-value: bold if significant, plain otherwise.
# Use with parse = TRUE in geom_signif or annotate("text").
fmt_p_plot <- function(p, threshold = 0.05) {
  label <- fmt_p(p)
  if (p < threshold) paste0('bold("', label, '")') else paste0('"', label, '"')
}

# Plotmath-parseable RM-ANOVA subtitle: each effect bold.italic if sig, italic if NS.
fmt_anova_sub <- function(age_p, time_p, int_p, threshold = 0.05) {
  wrap <- function(label, p) {
    txt <- sprintf("%s %s", label, fmt_p(p))
    if (p < threshold) sprintf('bold(italic("%s"))', txt) else sprintf('italic("%s")', txt)
  }
  paste('italic("RM-ANOVA:")', wrap("Age", age_p), wrap("Time", time_p),
        wrap("Int.", int_p), sep = " ~~ ")
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

classify_proteins_f3 <- function(pi_aging, pi_training_old,
                                  threshold = 0.05) {
  dplyr::case_when(
    pi_aging < threshold & pi_training_old < threshold ~ "Sig Both",
    pi_aging < threshold                               ~ "Sig Aging only",
    pi_training_old < threshold                        ~ "Sig Training only",
    TRUE                                               ~ "NS"
  ) |>
    factor(levels = c("Sig Both",
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

clean_pathway_name <- function(name, max_chars = NULL) {
  name |>
    stringr::str_remove("^HALLMARK_") |>
    stringr::str_remove("^GOSLIM_") |>
    stringr::str_remove("^GOBP_") |>
    stringr::str_remove("^GOCC_") |>
    stringr::str_remove("^GOMF_") |>
    stringr::str_remove("^REACTOME_") |>
    stringr::str_remove("^KEGG_MEDICUS_") |>
    stringr::str_remove("^KEGG_") |>
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
                         "Receptor Protein Signaling")
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

# Contrast short-labels used in bars/facets/labeller. Single definition —
# CTR_FACET is kept as an alias because F02/panel_F.R refers to it by that
# name in its labeller() call.
CTR_SHORT <- c(
  Aging          = "Aging",
  Training_Young = "Tr.(Y)",
  Training_Old   = "Tr.(O)",
  Interaction    = "Inter."
)
CTR_FACET <- CTR_SHORT


get_pdf_device <- function() {
  # Prefer cairo_pdf (Linux/most setups). On macOS without XQuartz, fall back
  # to quartz(type = "pdf"), which embeds Unicode glyphs (Pi, rho, arrows)
  # via CoreText. Last resort: base pdf (no Greek/arrow glyph support).
  tryCatch(
    { cairo_pdf(tempfile()); dev.off(); cairo_pdf },
    error = function(e) {
      tryCatch(
        {
          fp <- tempfile(fileext = ".pdf")
          quartz(type = "pdf", file = fp); dev.off()
          function(filename, width, height, ...) {
            quartz(file = filename, type = "pdf",
                   width = width, height = height)
          }
        },
        error = function(e) "pdf"
      )
    }
  )
}
