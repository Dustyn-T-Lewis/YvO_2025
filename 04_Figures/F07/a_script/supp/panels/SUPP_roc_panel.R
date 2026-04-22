# F07 Supplementary — ROC panel (standalone SUPP)
#
# Presents all pilot ROC classifiers as a single grouped-column figure.
# Columns = classifier family; each cell = ROC curve + shaded AUC ribbon with
# AUC / CI / permutation p printed centered inside the curve.
#
# Reads: 04_Figures/F07/c_data/roc_pilot_{summary,curves}.csv (produced by pilot_rocs.R)
# Writes: 04_Figures/F07/b_reports/supp/{png,pdf}/panels/SUPP_F07_roc_panel.{png,pdf}
#
# Regen pilot data: Rscript 04_Figures/F07/a_script/supp/panels/pilot_rocs.R
#
# Method framework: same LOOCV + top-k + permutation engine as F07 panel_A
# (Ambroise & McLachlan 2002 PNAS, Varma & Simon 2006 BMC Bioinformatics).

setwd(rprojroot::find_rstudio_root_file())
source("04_Figures/shared/style.R")

suppressPackageStartupMessages({
  library(tidyverse); library(patchwork)
})

DAT     <- "04_Figures/F07/c_data"
RPT_PNG <- "04_Figures/F07/b_reports/supp/png/panels"
RPT_PDF <- "04_Figures/F07/b_reports/supp/pdf/panels"
dir.create(RPT_PNG, recursive = TRUE, showWarnings = FALSE)
dir.create(RPT_PDF, recursive = TRUE, showWarnings = FALSE)

stopifnot("pilot ROC data missing: run F07/a_script/supp/panels/pilot_rocs.R first" =
  file.exists(file.path(DAT, "roc_pilot_summary.csv")) &&
  file.exists(file.path(DAT, "roc_pilot_curves.csv")))
summ   <- read_csv(file.path(DAT, "roc_pilot_summary.csv"), show_col_types = FALSE)
curves <- read_csv(file.path(DAT, "roc_pilot_curves.csv"),  show_col_types = FALSE)

# ---- Classifier → column group + display name + color -----------------------
clf_meta <- tibble::tribble(
  ~classifier,                 ~group,              ~display,                   ~color,
  "Age|proteins",              "Age (subject)",     "Proteome (2113)",          "#D6604D",
  "Age|ME+pheno",              "Age (subject)",     "MEs + Phenotype",          "#9B7FBF",
  "Age|pheno",                 "Age (subject)",     "Phenotype only",           "#8D6E63",
  "AgingDEP|prot_features",    "Protein architecture", "Aging-DEP status",      "#1565C0",
  "Reversed|prot_features",    "Protein architecture", "Reversed vs Exacerb.",  "#00796B",
  "ΔVL-resp|ME_Pre",           "Training response", "Δ VL responder",          "#E65100",
  "ΔLBM-resp|ME_Pre",          "Training response", "Δ LBM responder",         "#AA336A"
)

col_order <- c("Age (subject)", "Protein architecture", "Training response")

summ <- summ %>% left_join(clf_meta, by = "classifier")
curves <- curves %>% left_join(clf_meta, by = "classifier")

# Significance marker (pilot_rocs uses n_perm=200, Phipson-Smyth floor ~0.005;
# <0.001 tier dropped as unreachable under current permutation budget)
sig_sym <- function(p) ifelse(is.na(p), "",
                       ifelse(p < 0.01,  "**",
                       ifelse(p < 0.05,  "*",
                       ifelse(p < 0.10,  "\u2020", "ns"))))
summ <- summ %>% mutate(sig = sig_sym(perm_p))

# ---- Single-ROC builder ------------------------------------------------------
build_roc <- function(row, df) {
  col <- row$color
  lab <- sprintf("AUC = %.2f %s\n95%% CI  %.2f \u2013 %.2f\np = %.3f\nn = %d  (%d / %d)",
                 row$auc, row$sig, row$ci_lo, row$ci_hi,
                 row$perm_p,
                 row$n, row$n_pos, row$n_neg)
  ggplot(df, aes(fpr, tpr)) +
    geom_ribbon(aes(ymin = 0, ymax = tpr), fill = col, alpha = 0.32) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed",
                color = "grey70", linewidth = 0.35) +
    geom_line(color = col, linewidth = 1.2) +
    annotate("text", x = 0.62, y = 0.28, label = lab,
             hjust = 0.5, vjust = 0.5, size = 3.3, fontface = "bold",
             color = "grey15", lineheight = 1.1) +
    scale_x_continuous(limits = c(0,1), breaks = c(0, 0.5, 1),
                       expand = c(0,0)) +
    scale_y_continuous(limits = c(0,1), breaks = c(0, 0.5, 1),
                       expand = c(0,0)) +
    coord_fixed() +
    labs(x = "1 − Specificity", y = "Sensitivity",
         title = row$display) +
    theme_classic(base_size = 10) +
    theme(plot.title   = element_text(face = "bold", size = FIG_TITLE_SIZE,
                                      color = "grey10",
                                      margin = margin(b = 2)),
          axis.title   = element_text(size = 5, face = "bold",
                                      color = "grey25"),
          axis.text    = element_text(size = FIG_AXIS_TEXT, color = "grey30"),
          axis.line    = element_line(color = "grey30", linewidth = 0.4),
          axis.ticks   = element_line(color = "grey30", linewidth = 0.4),
          plot.margin  = margin(4, 6, 4, 6))
}

# ---- Build each plot, keyed by classifier -----------------------------------
plot_list <- setNames(vector("list", nrow(summ)), summ$classifier)
for (i in seq_len(nrow(summ))) {
  cl <- summ$classifier[i]
  df <- curves %>% filter(classifier == cl)
  plot_list[[cl]] <- build_roc(summ[i,], df)
}

# ---- Assemble column blocks (3 columns, each stacked vertically) -------------
# Max rows = 3 (age column has 3 ROCs)
blank <- patchwork::plot_spacer()

col_plots <- lapply(col_order, function(g) {
  ids <- summ$classifier[summ$group == g]
  # order within group by descending AUC
  ids <- ids[order(-summ$auc[match(ids, summ$classifier)])]
  ps <- plot_list[ids]
  # pad to 3 rows
  while (length(ps) < 3) ps <- c(ps, list(blank))
  # Add column header via patchwork::plot_annotation on stacked col
  stacked <- ps[[1]] / ps[[2]] / ps[[3]]
  wrap_elements(stacked) +
    plot_annotation(title = g,
                    theme = theme(
                      plot.title = element_text(face = "bold", size = FIG_TITLE_SIZE,
                                                color = "grey10", hjust = 0.5)))
})

# Compose: 3 columns side by side with group headers
col_with_header <- function(g) {
  ids <- summ$classifier[summ$group == g]
  ids <- ids[order(-summ$auc[match(ids, summ$classifier)])]
  ps <- plot_list[ids]
  while (length(ps) < 3) ps <- c(ps, list(blank))
  header <- ggplot() + theme_void() +
    annotate("text", x = 0.5, y = 0.5, label = g,
             fontface = "bold", size = 3.5, color = "grey10") +
    theme(plot.margin = margin(2, 2, 2, 2))
  header / ps[[1]] / ps[[2]] / ps[[3]] +
    plot_layout(heights = c(0.18, 1, 1, 1))
}

composite <- col_with_header("Age (subject)") |
             col_with_header("Protein architecture") |
             col_with_header("Training response")

composite <- composite +
  plot_annotation(
    title    = "F07 Supplementary — ROC menu: what is predictable, what is not",
    subtitle = paste0(
      "LOOCV top-k logistic (subject-level, n=30) | 10-fold CV (protein-level) | ",
      "permutation p (Phipson\u2013Smyth, 200 perms)  **<0.01  *<0.05  \u2020<0.10  ns\u22650.10"),
    caption = paste0(
      "Columns group classifiers by outcome family. ROC shading = area under curve; ",
      "dashed diagonal = chance. Pilot data: 04_Figures/F07/c_data/roc_pilot/."),
    theme = theme(
      plot.title    = element_text(face = "bold", size = FIG_TITLE_SIZE, color = "grey10",
                                   margin = margin(b = 2)),
      plot.subtitle = element_text(size = FIG_SUBTITLE_SIZE, face = "italic", color = "grey30",
                                   margin = margin(b = 6)),
      plot.caption  = element_text(size = FIG_AXIS_TEXT, color = "grey45", hjust = 0)))

W_in <- 11; H_in <- 11
ggsave(file.path(RPT_PDF, "SUPP_F07_roc_panel.pdf"), composite,
       width = W_in, height = H_in, units = "in",
       device = get_pdf_device())
ggsave(file.path(RPT_PNG, "SUPP_F07_roc_panel.png"), composite,
       width = W_in, height = H_in, units = "in", dpi = 300)

message("SUPP_F07_roc_panel written")
