# style.R — palettes, themes, sizing, helpers for all figure scripts

library(ggplot2)
library(scales)
library(grid)

# suppress stray Rplots.pdf from implicit device opens
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
SIG_LABEL_TEXT_F2 <- setNames(rep("white", 5), names(SIG_LABEL_FILL_F2))

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
SIG_LABEL_TEXT_F3 <- setNames(rep("white", 4), names(SIG_LABEL_FILL_F3))

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


PANEL_MD  <- 178   # J Physiol double-column width

# base annotation sizes (mm) at PANEL_MD; scale_text() adjusts for other widths
BASE_PATHWAY  <- 2.8   # ~8pt pathway labels
BASE_GENE     <- 2.5   # ~7pt gene labels
BASE_STAT     <- 2.5   # ~7pt stat annotations
BASE_QUADRANT <- 2.8   # ~8pt quadrant labels
BASE_COUNT    <- 2.5   # ~7pt bar counts
BASE_TAG      <- 8     # panel tags

scale_text <- function(base_size, panel_width_mm, ref_width = PANEL_MD) {
  base_size * sqrt(panel_width_mm / ref_width)
}

# strip decorations before embedding in composite (call after standalone save)
strip_for_composite <- function(p) {
  p + labs(title = NULL, subtitle = NULL, tag = NULL) +
    theme(legend.position = "none")
}

is_light_color <- function(color_name) {
  rgb_val <- col2rgb(color_name)
  (0.299 * rgb_val[1] + 0.587 * rgb_val[2] + 0.114 * rgb_val[3]) / 255 > 0.6
}


# text hierarchy — J Physiol spec
FIG_TITLE_SIZE    <- 7
FIG_SUBTITLE_SIZE <- 4
FIG_STRIP_SIZE    <- 5
FIG_AXIS_TEXT     <- 5
FIG_LEGEND_TITLE  <- 5
FIG_LEGEND_TEXT   <- 4

composite_text_sizes <- function(comp_h_mm) {
  list(
    title    = pmax(6, pmin(8, round(5 + comp_h_mm / 80))),
    subtitle = pmax(4, pmin(6, round(3 + comp_h_mm / 100))),
    tag      = 8
  )
}

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

fmt_p_plot <- function(p, threshold = 0.05) {
  label <- fmt_p(p)
  if (p < threshold) paste0('bold("', label, '")') else paste0('"', label, '"')
}

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

# Bonett & Wright 2000 — Fisher z CI for r (k = number of covariates)
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

boot_median_ci <- function(x, R = 2000, conf = 0.95) {
  meds <- replicate(R, median(sample(x, replace = TRUE)))
  qs   <- quantile(meds, c((1 - conf) / 2, (1 + conf) / 2))
  c(lower = unname(qs[1]), upper = unname(qs[2]))
}

.DB_PREFIXES <- c("^HALLMARK_", "^GOSLIM_", "^GOBP_", "^GOCC_", "^GOMF_",
                   "^REACTOME_", "^KEGG_MEDICUS_", "^KEGG_")

.SCI_CAPS <- c(
  "Mtorc1"  = "mTORC1", "Myc "    = "MYC ",   "E2f "    = "E2F ",
  "Dna "    = "DNA ",   "Rna "    = "RNA ",   "Tnfa "   = "TNFa ",
  "Uv "     = "UV ",    "G2m "    = "G2M ",   "Il6 "    = "IL6 ",
  "Il2 "    = "IL2 ",   "Kras "   = "KRAS ",  "P53 "    = "p53 ",
  "Tgf "    = "TGF ",   "Nfkb"    = "NF-kB",  "Atp "    = "ATP ",
  "Nadh "   = "NADH ",
  "Oxidative Phosphorylation"                        = "OXPHOS",
  "External Encapsulating Structure Or.*"            = "Extracellular Matrix Organization",
  "Enzyme Linked Receptor Protein Signaling.*"       = "Receptor Protein Signaling"
)

clean_pathway_name <- function(name) {
  out <- name
  for (pfx in .DB_PREFIXES) out <- stringr::str_remove(out, pfx)
  out <- stringr::str_replace_all(out, "_", " ")
  out <- stringr::str_to_title(out)
  for (i in seq_along(.SCI_CAPS)) {
    out <- stringr::str_replace(out, names(.SCI_CAPS)[i], .SCI_CAPS[i])
  }
  out
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

CTR_SHORT <- c(
  Aging          = "Aging",
  Training_Young = "Tr.(Y)",
  Training_Old   = "Tr.(O)",
  Interaction    = "Inter."
)
CTR_FACET <- CTR_SHORT
CONTRAST_ORDER <- c("Aging", "Training_Young", "Training_Old", "Interaction")


get_pdf_device <- function() {
  # cairo_pdf > quartz > base pdf
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
