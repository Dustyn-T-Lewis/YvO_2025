#!/usr/bin/env Rscript
# Figure 7 — Panel B: Module x Phenotype coupling is age-dependent (2x3)
#
# 6 hero scatters — one per strong age-flipped (or dual-significant)
# delta-ME -> delta-pheno relationship.
#
# Sourced by 01_main_panels.R — expects style.R + figure_supplement_helpers.R
# already loaded.
# Exports: pB, pB_title, pB_subtitle, pB_legend

library(tidyverse)
library(patchwork)
library(ppcor)
library(ggtext)

BASE     <- "04_Figures/F07"
RPT_PNG  <- file.path(BASE, "b_reports", "main", "png", "panels")
RPT_PDF  <- file.path(BASE, "b_reports", "main", "pdf", "panels")
DAT_OUT  <- file.path(BASE, "c_data")
dir.create(RPT_PNG, recursive = TRUE, showWarnings = FALSE)
dir.create(RPT_PDF, recursive = TRUE, showWarnings = FALSE)
dir.create(DAT_OUT, recursive = TRUE, showWarnings = FALSE)

F06_SUPP <- "04_Figures/F06/c_data/F06_supplementary.xlsx"
stopifnot("F06 stitcher must run first: missing F06_supplementary.xlsx" =
  file.exists(F06_SUPP))

me_pre     <- read_matrix_sheet(F06_SUPP, "me_pre",   "subject_key")
delta_me   <- read_matrix_sheet(F06_SUPP, "delta_me", "subject_key")
pheno_wide <- read_sheet_df(F06_SUPP, "metadata_pheno_wide")
subj_age   <- read_sheet_df(F06_SUPP, "metadata_subj_age")
common_subj    <- read_vector_sheet(F06_SUPP, "common_subj")
mod_bio_df     <- read_sheet_df(F06_SUPP, "WGCNA_mod_bio_labels")
mod_bio_labels <- setNames(mod_bio_df$bio_label, mod_bio_df$module_color)

pdf_device <- get_pdf_device()

message("Panel B: 2x3 hero scatter grid (\u0394ME x \u0394pheno)...")

# -- Outcome + predictor labels -----------------------------------------------
outcome_nice <- c(delta_VL = "\u0394VL (cm)",
                  delta_LBM = "\u0394LBM (kg)",
                  delta_DL  = "\u0394DL (kg)",
                  delta_T1  = "\u0394T1",
                  delta_T2  = "\u0394T2")

pretty_mod <- function(m) {
  m_short <- gsub("^ME", "", m)
  bio <- mod_bio_labels[m_short]
  if (is.na(bio)) str_to_title(m_short) else bio
}

# -- Hero picks (top-6 by min within-age p, sign-flipping preferred) -----------
HERO_PICKS <- tribble(
  ~x_source,   ~module,        ~outcome,
  "delta_ME",  "MEgreen",      "delta_VL",   # Y p=0.004 ***
  "baseline",  "MEturquoise",  "delta_VL",   # Y p=0.012 *  (was Panel C)
  "delta_ME",  "MEbrown",      "delta_LBM",  # Y p=0.014 *
  "delta_ME",  "MEturquoise",  "delta_T1",   # Y p=0.021 *
  "baseline",  "MEyellow",     "delta_T2",   # Y p=0.022 *
  "delta_ME",  "MEyellow",     "delta_T2"    # dual-sig: Y *, O dagger
)

get_x <- function(x_source, mod) {
  src <- if (x_source == "baseline") me_pre else delta_me
  as.numeric(src[common_subj, mod])
}

# -- Helpers -------------------------------------------------------------------
fmt_p <- function(p) {
  if (is.na(p)) return("n/a")
  if (p < 0.001) "p<0.001" else sprintf("p=%.3f", p)
}
sig_mark <- function(p) {
  if (is.na(p)) ""
  else if (p < 0.001) "***"
  else if (p < 0.01)  "**"
  else if (p < 0.05)  "*"
  else if (p < 0.10)  "\u2020"
  else ""
}

build_hero_mini <- function(x_source, module, outcome) {
  mod <- module
  x <- get_x(x_source, mod)
  y <- pheno_wide[[outcome]][match(common_subj, pheno_wide$subject_key)]
  age <- subj_age$age[match(common_subj, subj_age$subject_key)]

  d <- tibble(x = x, y = y, age = factor(age, levels = c("Young", "Old"))) |>
    filter(!is.na(x), !is.na(y))

  cor_by_age <- function(grp) {
    sub <- d |> filter(age == grp)
    if (nrow(sub) < 4) return(c(r = NA_real_, p = NA_real_, n = nrow(sub)))
    ct <- suppressWarnings(cor.test(sub$x, sub$y))
    c(r = unname(ct$estimate), p = ct$p.value, n = nrow(sub))
  }
  cY <- cor_by_age("Young"); cO <- cor_by_age("Old")

  y_stat <- sprintf("Y (n=%d): r=%+.2f, %s%s",
                    cY["n"], cY["r"], fmt_p(cY["p"]), sig_mark(cY["p"]))
  o_stat <- sprintf("O (n=%d): r=%+.2f, %s%s",
                    cO["n"], cO["r"], fmt_p(cO["p"]), sig_mark(cO["p"]))

  x_prefix <- if (x_source == "baseline") "Baseline" else "\u0394"
  x_lab <- "Eigengene value (a.u.)"
  y_lab <- unname(outcome_nice[outcome])
  title <- sprintf("%s %s \u2192 %s",
                   x_prefix, pretty_mod(mod),
                   str_replace(y_lab, " \\(.*\\)", ""))

  stat_label <- sprintf(
    "<span style='color:%s'><b>%s</b></span><br><span style='color:%s'><b>%s</b></span>",
    AGE_COLORS[["Young"]], y_stat, AGE_COLORS[["Old"]], o_stat)

  ggplot(d, aes(x = x, y = y, color = age, fill = age)) +
    geom_hline(yintercept = 0, linetype = "dotted",
               color = "grey70", linewidth = 0.2) +
    geom_smooth(data = filter(d, age == "Young"),
                method = "lm", se = TRUE, linewidth = 0.6, alpha = 0.13) +
    geom_smooth(data = filter(d, age == "Old"),
                method = "lm", se = TRUE, linewidth = 0.6,
                alpha = 0.11, linetype = "dashed") +
    geom_point(shape = 21, size = 2.4, stroke = 0.4,
               color = "white", alpha = 0.95) +
    scale_color_manual(values = AGE_COLORS) +
    scale_fill_manual(values = AGE_COLORS) +
    scale_y_continuous(expand = expansion(mult = c(0.25, 0.12))) +
    labs(x = x_lab, y = y_lab) +
    # Title inside plot area at top centre
    annotate("text", x = mean(range(d$x, na.rm = TRUE)), y = Inf,
             label = title, hjust = 0.5, vjust = 1.3,
             size = 2.8, fontface = "bold", color = "grey10") +
    ggtext::geom_richtext(
      inherit.aes = FALSE,
      data = tibble(x = mean(range(d$x, na.rm = TRUE)), y = -Inf, lab = stat_label),
      aes(x = x, y = y, label = lab),
      hjust = 0.5, vjust = -0.1,
      size = 2.2, lineheight = 1.0,
      label.color = "grey70", label.r = unit(1.2, "pt"),
      label.padding = unit(c(1.5, 2.5, 1.5, 2.5), "pt"),
      label.margin = unit(c(2, 2, 2, 2), "pt"),
      fill = scales::alpha("white", 0.88)) +
    FIG_THEME +
    coord_cartesian(clip = "off") +
    theme(axis.title.x = element_text(size = 5, face = "bold",
                                      margin = margin(t = 3)),
          axis.title.y = element_text(size = 5, face = "bold",
                                      margin = margin(r = -16, unit = "pt")),
          axis.text  = element_text(size = FIG_AXIS_TEXT),
          legend.position = "none",
          plot.margin = margin(4, 2, 4, 2))
}

# -- Full screening + BH correction audit -------------------------------------
MODULES_SCREEN <- grep("^ME", colnames(delta_me), value = TRUE)
MODULES_SCREEN <- setdiff(MODULES_SCREEN, "MEgrey")
SOURCES <- c("delta_ME", "baseline")
OUTCOMES <- names(outcome_nice)

screen_rows <- list()
for (src in SOURCES) {
  src_mat <- if (src == "baseline") me_pre else delta_me
  for (mod in MODULES_SCREEN) {
    x <- as.numeric(src_mat[common_subj, mod])
    for (out in OUTCOMES) {
      y <- pheno_wide[[out]][match(common_subj, pheno_wide$subject_key)]
      for (grp in c("Young", "Old")) {
        age <- subj_age$age[match(common_subj, subj_age$subject_key)]
        keep <- !is.na(x) & !is.na(y) & age == grp
        if (sum(keep) < 4) next
        ct <- suppressWarnings(cor.test(x[keep], y[keep]))
        screen_rows[[length(screen_rows) + 1]] <- tibble(
          source = src, module = mod, outcome = out, stratum = grp,
          n = sum(keep), r = unname(ct$estimate), p_raw = ct$p.value)
      }
    }
  }
}

screen_df <- bind_rows(screen_rows) |>
  mutate(p_bh = p.adjust(p_raw, method = "BH")) |>
  arrange(p_raw)

n_screen   <- nrow(screen_df)
n_raw_sig  <- sum(screen_df$p_raw < 0.05)
n_bh_sig   <- sum(screen_df$p_bh < 0.05)
message(sprintf("  Full screen: %d tests, %d raw p<0.05, %d BH p<0.05",
                n_screen, n_raw_sig, n_bh_sig))

write_csv(screen_df, file.path(DAT_OUT, "panel_B_full_screen_bh.csv"))

# Validate hero picks against full screening — warn if any pick has no raw p<0.10
hero_check <- HERO_PICKS |>
  mutate(key = paste(x_source, module, outcome, sep = "|")) |>
  left_join(screen_df |> mutate(key = paste(source, module, outcome, sep = "|")),
            by = "key", multiple = "all") |>
  filter(p_raw < 0.10)
if (nrow(hero_check) < nrow(HERO_PICKS)) {
  warning("Some HERO_PICKS have no stratum with raw p<0.10 — consider updating")
}

hero_plots <- pmap(HERO_PICKS, build_hero_mini)

# Remove x-axis title from all but the bottom row (indices 4,5,6 in ncol=3 layout)
for (i in seq_along(hero_plots)) {
  if (i <= length(hero_plots) - 3) {
    hero_plots[[i]] <- hero_plots[[i]] +
      theme(axis.title.x = element_blank())
  }
}

# Shared legend at the bottom (one line)
legend_data <- tibble(age = factor(c("Young", "Old"),
                                    levels = c("Young", "Old")))
p_legend <- ggplot(legend_data, aes(x = 1, y = 1, color = age)) +
  geom_point(size = 3) +
  scale_color_manual(values = AGE_COLORS, name = "Age") +
  FIG_THEME +
  theme(legend.position = "bottom",
        legend.title = element_text(face = "bold", size = FIG_LEGEND_TITLE),
        legend.text  = element_text(size = FIG_LEGEND_TEXT))
legend_grob <- cowplot::get_legend(p_legend)

panel_B <- wrap_plots(hero_plots, ncol = 3) +
  plot_annotation(
    title    = "Module\u2013Phenotype Coupling (Age-Dependent)",
    subtitle = sprintf("Stratified Pearson r | %d/%d raw p<0.05 | 0/%d BH-sig",
                       n_raw_sig, n_screen, n_screen),
    theme = theme(
      plot.title    = element_text(face = "bold", size = FIG_TITLE_SIZE,
                                   lineheight = 1.15,
                                   margin = margin(t = 6, l = 34, b = 3,
                                                   unit = "pt")),
      plot.subtitle = element_text(size = FIG_SUBTITLE_SIZE,
                                   face = "bold.italic", color = "grey30",
                                   lineheight = 1.2,
                                   margin = margin(t = 2, l = 34, b = 6,
                                                   unit = "pt"))
    )
  )

PB_W <- 178
PB_H <- 90

ggsave(file.path(RPT_PNG, "MAIN_panel_B_hero_grid.png"),
       panel_B, width = PB_W, height = PB_H, units = "mm", dpi = 300)
ggsave(file.path(RPT_PDF, "MAIN_panel_B_hero_grid.pdf"),
       panel_B, width = PB_W, height = PB_H, units = "mm",
       device = get_pdf_device())

message("  MAIN_panel_B_hero_grid saved (2x3 hero grid)")

# Export for composite
pB_title    <- "Module\u2013Phenotype Coupling (Age-Dependent)"
pB_subtitle <- sprintf("Stratified Pearson r | %d/%d raw p<0.05 | 0/%d BH-sig",
                       n_raw_sig, n_screen, n_screen)
pB_legend   <- NULL
pB          <- panel_B +
  plot_annotation(title = NULL, subtitle = NULL)
