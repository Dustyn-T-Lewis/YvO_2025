#!/usr/bin/env Rscript
# F07 Supplementary — Panel B full-sweep grid
#
# Shows the full 2-source x 10-module x 5-outcome scatter landscape
# underlying the 6 hero cells in MAIN Panel B.
#
# Sourced by 01_main_panels.R — expects style.R + figure_supplement_helpers.R
# already loaded.
# Reads panel_B_full_screen_bh.csv from c_data (or falls back to xlsx).

suppressPackageStartupMessages({
  library(tidyverse); library(patchwork)
})

BASE     <- here::here("04_Figures", "F07")
DAT_F07  <- file.path(BASE, "c_data")
RPT_PNG  <- file.path(BASE, "b_reports", "supp", "png", "panels")
RPT_PDF  <- file.path(BASE, "b_reports", "supp", "pdf", "panels")
dir.create(RPT_PNG, recursive = TRUE, showWarnings = FALSE)
dir.create(RPT_PDF, recursive = TRUE, showWarnings = FALSE)

F06_SUPP <- here::here("04_Figures", "F06", "c_data", "F06_supplementary.xlsx")
stopifnot("F06 stitcher must run first: missing F06_supplementary.xlsx" =
  file.exists(F06_SUPP))

me_pre     <- read_matrix_sheet(F06_SUPP, "me_pre",   "subject_key")
delta_me   <- read_matrix_sheet(F06_SUPP, "delta_me", "subject_key")
pheno_wide <- read_sheet_df(F06_SUPP, "metadata_pheno_wide")
subj_age   <- read_sheet_df(F06_SUPP, "metadata_subj_age")
common_subj    <- read_vector_sheet(F06_SUPP, "common_subj")
mod_bio_df     <- read_sheet_df(F06_SUPP, "WGCNA_mod_bio_labels")
mod_bio_labels <- setNames(mod_bio_df$bio_label, mod_bio_df$module_color)

F07_SUPP <- file.path(DAT_F07, "F07_supplementary.xlsx")
screen_src <- file.path(DAT_F07, "panel_B_full_screen_bh.csv")
screen <- if (file.exists(screen_src)) {
  read_csv(screen_src, show_col_types = FALSE)
} else {
  read_sheet_df(F07_SUPP, "panel_B_full_screen_bh")
}

OUTCOMES <- c("delta_VL", "delta_LBM", "delta_DL", "delta_T1", "delta_T2")
outcome_nice <- c(delta_VL = "\u0394VL", delta_LBM = "\u0394LBM",
                  delta_DL = "\u0394DL", delta_T1 = "\u0394T1",
                  delta_T2 = "\u0394T2")

MODULES <- grep("^ME", colnames(delta_me), value = TRUE)
MODULES <- setdiff(MODULES, "MEgrey")
# Stable order by WGCNA convention
mod_order <- c("MEturquoise","MEblue","MEbrown","MEyellow","MEgreen",
               "MEred","MEblack","MEpink","MEmagenta","MEpurple")
MODULES <- mod_order[mod_order %in% MODULES]

pretty_mod <- function(m) {
  short <- gsub("^ME", "", m)
  bio   <- mod_bio_labels[short]
  if (is.na(bio)) stringr::str_to_title(short) else bio
}

age_vec <- subj_age$age[match(common_subj, subj_age$subject_key)]

# Build one cell: scatter + border from Young-stratum p_raw / p_bh
build_cell <- function(source, module, outcome) {
  src_mat <- if (source == "baseline") me_pre else delta_me
  x <- as.numeric(src_mat[common_subj, module])
  y <- pheno_wide[[outcome]][match(common_subj, pheno_wide$subject_key)]
  d <- tibble(x = x, y = y, age = factor(age_vec, levels = c("Young","Old"))) |>
    filter(!is.na(x), !is.na(y))

  rY <- screen |> filter(source == !!source, module == !!module,
                          outcome == !!outcome, stratum == "Young")
  rO <- screen |> filter(source == !!source, module == !!module,
                          outcome == !!outcome, stratum == "Old")

  y_r <- if (nrow(rY) == 1) sprintf("Y %+.2f", rY$r) else "Y n/a"
  o_r <- if (nrow(rO) == 1) sprintf("O %+.2f", rO$r) else "O n/a"

  p_bh_y  <- if (nrow(rY) == 1) rY$p_bh  else NA_real_
  p_raw_y <- if (nrow(rY) == 1) rY$p_raw else NA_real_
  sig <- case_when(
    is.na(p_bh_y)   ~ "ns",
    p_bh_y  < 0.05  ~ "q<.05",
    p_raw_y < 0.05  ~ "p<.05",
    TRUE            ~ "ns")
  border_color <- ifelse(sig == "ns", "grey85", "black")
  border_lty   <- ifelse(sig == "p<.05", "dashed", "solid")
  border_lw    <- case_when(sig == "q<.05" ~ 1.6,
                            sig == "p<.05" ~ 1.6,
                            TRUE           ~ 0.25)

  ggplot(d, aes(x = x, y = y, color = age)) +
    geom_smooth(data = filter(d, age == "Young"),
                method = "lm", se = FALSE, linewidth = 0.6) +
    geom_smooth(data = filter(d, age == "Old"),
                method = "lm", se = FALSE, linewidth = 0.6, linetype = "dashed") +
    geom_point(size = 1.1, alpha = 0.9) +
    scale_color_manual(values = AGE_COLORS, guide = "none") +
    annotate("label", x = -Inf, y = Inf,
             label = paste(y_r, o_r, sep = "\n"),
             hjust = 0, vjust = 1, size = 3.6, fontface = "bold",
             color = "grey15",
             fill = alpha("white", 0.88), linewidth = 0,
             label.r = unit(0.6, "mm"),
             label.padding = unit(1.0, "mm"),
             lineheight = 1.1) +
    theme_void() +
    theme(panel.border = element_rect(fill = NA,
                                      color = border_color,
                                      linewidth = border_lw,
                                      linetype = border_lty),
          plot.margin = margin(2, 2, 2, 2))
}

mod_header <- function(mod) {
  short <- gsub("^ME", "", mod)
  ggplot() + theme_void() +
    annotate("rect", xmin = 0, xmax = 1, ymin = 0.55, ymax = 1, fill = short) +
    annotate("text", x = 0.5, y = 0.26, label = stringr::str_to_title(short),
             size = 3.0, fontface = "bold", color = "grey10") +
    coord_cartesian(xlim = c(0, 1), ylim = c(0, 1), expand = FALSE) +
    theme(plot.margin = margin(1, 2, 1, 2))
}

row_label <- function(txt) {
  ggplot() + theme_void() +
    annotate("text", x = 0.5, y = 0.5, label = txt, angle = 90,
             fontface = "bold", size = 4, color = "grey10") +
    coord_cartesian(xlim = c(0, 1), ylim = c(0, 1), expand = FALSE)
}

block_label <- function(txt) {
  ggplot() + theme_void() +
    annotate("text", x = 0.5, y = 0.5, label = txt,
             fontface = "bold", size = 5, color = "grey10") +
    coord_cartesian(xlim = c(0, 1), ylim = c(0, 1), expand = FALSE)
}

build_block <- function(source_name) {
  header_cells <- c(list(patchwork::plot_spacer()),
                    lapply(MODULES, mod_header))
  rows <- lapply(OUTCOMES, function(out) {
    c(list(row_label(outcome_nice[out])),
      lapply(MODULES, function(m) build_cell(source_name, m, out)))
  })
  all_cells <- c(header_cells, unlist(rows, recursive = FALSE))
  wrap_plots(all_cells,
             nrow = 1 + length(OUTCOMES), ncol = 1 + length(MODULES),
             widths  = c(0.45, rep(1, length(MODULES))),
             heights = c(0.30, rep(1, length(OUTCOMES))))
}

message("Building block A (Delta ME) ...")
block_A <- build_block("delta_ME")
message("Building block B (Baseline ME) ...")
block_B <- build_block("baseline")

lbl_A <- block_label("Predictor: \u0394 ME")
lbl_B <- block_label("Predictor: Baseline ME")

composite <- (lbl_A / block_A / lbl_B / block_B) +
  plot_layout(heights = c(0.04, 1, 0.04, 1)) +
  plot_annotation(
    title = "F07 Supplementary \u2014 Panel B full sweep: module eigengenes vs \u0394phenotype, age-stratified",
    subtitle = paste0(
      "Each cell = per-age Pearson r (Young solid / Old dashed OLS).\n",
      "Border: solid black = BH q<0.05; dashed = nominal p<0.05 in Young; grey = ns.  Screen: 200 tests; 6 raw-sig, 0 BH-sig."),
    caption = paste0(
      "Inline stats = Pearson r per stratum (Y, O). ",
      "Data: 04_Figures/F07/c_data/panel_B_full_screen_bh.csv. ",
      "Hero picks in MAIN Panel B (MEgreen\u2192\u0394VL, MEturquoise\u2192\u0394VL base, ",
      "MEbrown\u2192\u0394LBM, MEturquoise\u2192\u0394T1, MEyellow\u2192\u0394T2 base, MEyellow\u2192\u0394T2)."),
    theme = theme(
      plot.title    = element_text(face = "bold", size = FIG_TITLE_SIZE + 4, color = "grey10",
                                   margin = margin(t = 14, l = 52, b = 4, unit = "pt")),
      plot.subtitle = element_text(size = FIG_SUBTITLE_SIZE + 4, face = "italic", color = "grey30",
                                   margin = margin(l = 52, b = 12, unit = "pt"),
                                   lineheight = 1.25),
      plot.caption  = element_text(size = FIG_AXIS_TEXT, color = "grey45", hjust = 0,
                                   lineheight = 1.3)))

W_in <- 14; H_in <- 16
ggsave(file.path(RPT_PDF, "SUPP_F07_panel_B_grid.pdf"), composite,
       width = W_in, height = H_in, units = "in",
       device = get_pdf_device())
ggsave(file.path(RPT_PNG, "SUPP_F07_panel_B_grid.png"), composite,
       width = W_in, height = H_in, units = "in", dpi = 300)

message("Wrote: SUPP_F07_panel_B_grid.{pdf,png}")
