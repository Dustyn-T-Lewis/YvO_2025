# Figure 6 — Panel B: Combined NES Scatter (stacked)
# Top:    Training Concordance (TY vs TO)
# Bottom: Aging Reversal (Aging vs TO)
# Output: b_reports/main/png/panels/MAIN_panel_B_scatters.png

setwd(rprojroot::find_rstudio_root_file())
source("04_Figures/shared/style.R")
source("04_Figures/shared/pathway_utils.R")

library(readr)
library(dplyr)
library(tidyr)
library(tibble)
library(stringr)
library(ggrepel)
library(patchwork)
library(cowplot)

RPT_PNG <- "04_Figures/F06/b_reports/main/png/panels"
RPT_PDF <- "04_Figures/F06/b_reports/main/pdf/panels"
DAT <- "04_Figures/F06/c_data"
dir.create(RPT_PNG, recursive = TRUE, showWarnings = FALSE)
dir.create(RPT_PDF, recursive = TRUE, showWarnings = FALSE)

# --- Input validation ---
stopifnot(
  "DEP results missing: 03_DEP/c_data/03_combined_results.csv" =
    file.exists("03_DEP/c_data/03_combined_results.csv"),
  "WGCNA module assignments missing — run YvO_WGCNA_run.R first" =
    file.exists(file.path(DAT, "wgcna/wgcna_module_assignments.csv")),
  "Module bio-labels missing — run YvO_WGCNA_run.R first" =
    file.exists(file.path(DAT, "mod_bio_labels.csv"))
)

# --- Load data (same as panel_B_module_nes_scatter.R) ---
module_df <- read_csv(file.path(DAT, "wgcna/wgcna_module_assignments.csv"), show_col_types = FALSE)
mod_bio   <- read_csv(file.path(DAT, "mod_bio_labels.csv"), show_col_types = FALSE)
combined  <- read_csv("03_DEP/c_data/03_combined_results.csv", show_col_types = FALSE)

module_sets <- module_df %>%
  filter(module_color != "grey", !is.na(gene), gene != "") %>%
  split(.$module_color) %>%
  lapply(function(x) x$gene)
mod_sizes <- sapply(module_sets, length)

build_ranks <- function(df, col) {
  vals <- df[[col]]; names(vals) <- df$gene
  vals <- vals[!is.na(vals)]; sort(vals, decreasing = TRUE)
}
ranks_TY    <- build_ranks(combined, "t_Training_Young")
ranks_TO    <- build_ranks(combined, "t_Training_Old")
ranks_Aging <- build_ranks(combined, "t_Aging")

run_module_fgsea <- function(ranks, module_sets) {
  res <- fgsea::fgseaMultilevel(pathways = module_sets, stats = ranks,
                                 minSize = 15, maxSize = 500, nPermSimple = 10000, eps = 0)
  # padj from fgseaMultilevel is BH-adjusted internally — no re-adjustment needed
  as.data.frame(res)
}
fgsea_TY    <- run_module_fgsea(ranks_TY, module_sets)
fgsea_TO    <- run_module_fgsea(ranks_TO, module_sets)
fgsea_Aging <- run_module_fgsea(ranks_Aging, module_sets)

merge_fgsea <- function(res, suffix) {
  res %>% as_tibble() %>%
    dplyr::select(pathway, NES, padj, size) %>%
    dplyr::rename_with(~ paste0(., "_", suffix), c(NES, padj, size))
}
fgsea_wide <- merge_fgsea(fgsea_TY, "TY") %>%
  left_join(merge_fgsea(fgsea_TO, "TO"), by = "pathway") %>%
  left_join(merge_fgsea(fgsea_Aging, "Aging"), by = "pathway") %>%
  mutate(module_color = pathway, n_proteins = mod_sizes[pathway])
fgsea_wide <- fgsea_wide %>%
  left_join(mod_bio %>% dplyr::select(module_color, bio_label), by = "module_color") %>%
  mutate(bio_label = ifelse(is.na(bio_label), stringr::str_to_title(module_color), bio_label),
         bio_label = gsub("Cell Cycle/Proteostasis", "Cell Cycle/\nProteostasis", bio_label),
         bio_label = gsub("Mitochondrial Biogenesis", "Mitochondrial\nBiogenesis", bio_label),
         bio_label = gsub("Oxidative Phosphorylation", "Oxidative\nPhosphorylation", bio_label),
         sig_conc = case_when(
           !is.na(padj_TY) & padj_TY < 0.05 & !is.na(padj_TO) & padj_TO < 0.05 ~ "Both",
           !is.na(padj_TY) & padj_TY < 0.05 ~ "Young only",
           !is.na(padj_TO) & padj_TO < 0.05 ~ "Old only", TRUE ~ "NS"),
         sig_rev = case_when(
           !is.na(padj_Aging) & padj_Aging < 0.05 & !is.na(padj_TO) & padj_TO < 0.05 ~ "Both",
           !is.na(padj_Aging) & padj_Aging < 0.05 ~ "Aging only",
           !is.na(padj_TO) & padj_TO < 0.05 ~ "Training only", TRUE ~ "NS"))

# --- Scatter plot helper (adapted for tighter pilot formatting) ---
build_scatter <- function(df, x_col, y_col, x_lab, y_lab, title, quad_labels, sig_col) {
  x_vals <- df[[x_col]]; y_vals <- df[[y_col]]
  nes_lim <- max(abs(c(x_vals, y_vals)), na.rm = TRUE) * 1.35

  sp <- cor.test(x_vals, y_vals, method = "spearman")
  subtitle_txt <- sprintf("n = %d modules | \u03c1 = %.2f%s",
    nrow(df), sp$estimate,
    ifelse(sp$p.value < 0.001, ", p < 0.001", sprintf(", p = %.3f", sp$p.value)))

  q_tr <- sum(x_vals > 0 & y_vals > 0, na.rm = TRUE)
  q_bl <- sum(x_vals < 0 & y_vals < 0, na.rm = TRUE)
  q_tl <- sum(x_vals < 0 & y_vals > 0, na.rm = TRUE)
  q_br <- sum(x_vals > 0 & y_vals < 0, na.rm = TRUE)

  df$label_col <- ifelse(sapply(df$module_color, is_light_color), "black", "white")
  # Force green and pink to black text — light enough for readability
  df$label_col[df$module_color == "green"] <- "black"
  df$label_col[df$module_color == "pink"] <- "black"
  df$label_col[df$module_color == "magenta"] <- "white"

  txt_lab  <- 5.0
  txt_quad <- txt_lab   # same size as pathway labels

  ggplot(df, aes(x = .data[[x_col]], y = .data[[y_col]])) +
    annotate("rect", xmin = 0, xmax = Inf, ymin = 0, ymax = Inf,
             fill = quad_labels$fill[1], alpha = 0.20) +
    annotate("rect", xmin = -Inf, xmax = 0, ymin = -Inf, ymax = 0,
             fill = quad_labels$fill[2], alpha = 0.20) +
    annotate("rect", xmin = 0, xmax = Inf, ymin = -Inf, ymax = 0,
             fill = quad_labels$fill[3], alpha = 0.20) +
    annotate("rect", xmin = -Inf, xmax = 0, ymin = 0, ymax = Inf,
             fill = quad_labels$fill[4], alpha = 0.20) +
    geom_hline(yintercept = 0, color = "grey60", linewidth = 0.2) +
    geom_vline(xintercept = 0, color = "grey60", linewidth = 0.2) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "black", linewidth = 0.3) +
    geom_point(aes(size = n_proteins), fill = df$module_color, color = "black",
               shape = 21, alpha = 0.85, stroke = 0.5) +
    geom_label_repel(aes(label = bio_label),
      fill = scales::alpha(df$module_color, 0.85), color = df$label_col,
      size = txt_lab, fontface = "bold", lineheight = 0.7, max.overlaps = Inf,
      segment.size = 0.3, segment.color = "grey40", min.segment.length = 0,
      box.padding = 1.4, point.padding = 1.0, force = 100, force_pull = 0.04,
      label.padding = unit(1.5, "pt"), label.r = unit(1.5, "pt"),
      linewidth = 0, seed = 7, show.legend = FALSE) +
    annotate("label", x = Inf, y = Inf,
             label = sprintf("%s  n=%d", quad_labels$label[1], q_tr),
             hjust = 1, vjust = 1, size = txt_quad, fontface = "bold",
             color = quad_labels$text_col[1], fill = scales::alpha("white", 0.92),
             label.padding = unit(2, "pt")) +
    annotate("label", x = -Inf, y = -Inf,
             label = sprintf("%s  n=%d", quad_labels$label[2], q_bl),
             hjust = 0, vjust = 0, size = txt_quad, fontface = "bold",
             color = quad_labels$text_col[2], fill = scales::alpha("white", 0.92),
             label.padding = unit(2, "pt")) +
    annotate("label", x = Inf, y = -Inf,
             label = sprintf("%s  n=%d", quad_labels$label[3], q_br),
             hjust = 1, vjust = 0, size = txt_quad, fontface = "bold",
             color = quad_labels$text_col[3], fill = scales::alpha("white", 0.92),
             label.padding = unit(2, "pt")) +
    annotate("label", x = -Inf, y = Inf,
             label = sprintf("%s  n=%d", quad_labels$label[4], q_tl),
             hjust = 0, vjust = 1, size = txt_quad, fontface = "bold",
             color = quad_labels$text_col[4], fill = scales::alpha("white", 0.92),
             label.padding = unit(2, "pt")) +
    scale_size_continuous(range = c(3, 9), name = "Proteins",
                          breaks = c(50, 100, 200, 300)) +
    scale_x_continuous(expand = expansion(mult = 0.02)) +
    scale_y_continuous(expand = expansion(mult = 0.02)) +
    coord_fixed(ratio = 1, xlim = c(-nes_lim, nes_lim), ylim = c(-nes_lim, nes_lim), clip = "off") +
    labs(title = NULL, subtitle = NULL, x = x_lab, y = y_lab) +
    FIG_THEME +
    theme(
      axis.text     = element_text(size = 11, face = "bold", color = "grey30"),
      axis.title.x  = element_text(size = 12, face = "bold"),
      axis.title.y  = element_text(size = 12, face = "bold", margin = margin(r = -4)),
      legend.position = "none",
      plot.margin   = margin(2, 4, 2, -8)
    )
}

# --- Build both scatters ---
quad_conc <- list(
  label = c("Concordant Up", "Concordant Down", "Discordant", "Discordant"),
  fill = c(AGE_COLORS["Old"], AGE_COLORS["Old"], AGE_COLORS["Young"], AGE_COLORS["Young"]),
  text_col = c("#D6604D", "#D6604D", "#4393C3", "#4393C3"))

p_top <- build_scatter(fgsea_wide, "NES_TY", "NES_TO",
  "NES (Training Young)", "NES (Training Old)",
  "Training Concordance", quad_conc, "sig_conc")

quad_rev <- list(
  label = c("Exacerbated", "Exacerbated", "Reversed", "Reversed"),
  fill = c(AGE_COLORS["Old"], AGE_COLORS["Old"], "#2E7D32", "#2E7D32"),
  text_col = c("#D6604D", "#D6604D", "#4393C3", "#4393C3"))

p_bottom <- build_scatter(fgsea_wide, "NES_Aging", "NES_TO",
  "NES (Aging)", "NES (Training Old)",
  "Aging Reversal", quad_rev, "sig_rev")

# --- Scatters only (NO title, NO legend) ---
scatters_panel <- (p_top / p_bottom) +
  plot_layout(heights = c(1, 1))

PB_W <- 220
PB_H <- 270

ggsave(file.path(RPT_PNG, "MAIN_panel_B_scatters.png"), scatters_panel,
       width = PB_W, height = PB_H, units = "mm", dpi = 300, limitsize = FALSE)
ggsave(file.path(RPT_PDF, "MAIN_panel_B_scatters.pdf"), scatters_panel,
       width = PB_W, height = PB_H, units = "mm", device = get_pdf_device(), limitsize = FALSE)

# --- Legend only (separate file for stitcher alignment) ---
p_legend_src <- p_top +
  scale_size_continuous(range = c(3, 9), name = "Proteins",
    breaks = c(50, 100, 200, 300),
    guide = guide_legend(nrow = 1,
      override.aes = list(alpha = 0.7, fill = "grey60", stroke = 0))) +
  theme(legend.position = "bottom",
        legend.title = element_text(size = 12, face = "bold"),
        legend.text = element_text(size = 11),
        legend.key = element_rect(fill = NA, color = NA),
        legend.key.size = unit(5, "mm"),
        legend.background = element_rect(fill = NA, color = NA))
legend_grob <- cowplot::get_plot_component(p_legend_src, "guide-box-bottom", return_all = FALSE)
p_legend <- cowplot::ggdraw(legend_grob)

ggsave(file.path(RPT_PNG, "MAIN_panel_B_legend.png"), p_legend,
       width = 90, height = 16, units = "mm", dpi = 300, limitsize = FALSE)

# Write fGSEA data for supplementary workbook
write_csv(fgsea_wide, file.path(DAT, "panel_B_module_fgsea.csv"))

message("Panel B: scatters + legend saved")
