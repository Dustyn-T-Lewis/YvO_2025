# Sourced by 01_main_panels.R — expects style.R already loaded.

source(here::here("04_Figures_v2", "shared_functions", "F04-F06_pathway_utils.R"))

library(readr)
library(dplyr)
library(tidyr)
library(tibble)
library(stringr)
library(ggrepel)
library(patchwork)
library(cowplot)

BASE <- here::here("04_Figures_v2", "F06")

RPT_PNG <- file.path(BASE, "b_reports", "main", "png", "panels")
RPT_PDF <- file.path(BASE, "b_reports", "main", "pdf", "panels")
DAT <- file.path(BASE, "c_data")
dir.create(RPT_PNG, recursive = TRUE, showWarnings = FALSE)
dir.create(RPT_PDF, recursive = TRUE, showWarnings = FALSE)

stopifnot(
  "DEP results missing: 03_DEP/c_data/03_combined_results.csv" =
    file.exists(here::here("03_DEP", "c_data", "03_combined_results.csv")),
  "WGCNA module assignments missing — run YvO_WGCNA_run.R first" =
    file.exists(file.path(DAT, "wgcna/wgcna_module_assignments.csv")),
  "Module bio-labels missing — run YvO_WGCNA_run.R first" =
    file.exists(file.path(DAT, "mod_bio_labels.csv"))
)

module_df <- read_csv(file.path(DAT, "wgcna/wgcna_module_assignments.csv"), show_col_types = FALSE)
mod_bio <- read_csv(file.path(DAT, "mod_bio_labels.csv"), show_col_types = FALSE)
combined <- read_csv(here::here("03_DEP", "c_data", "03_combined_results.csv"))

module_df_filt <- module_df |>
  filter(module_color != "grey", !is.na(gene), gene != "")
module_sets <- split(module_df_filt, module_df_filt$module_color) |>
  lapply(function(x) x$gene)
mod_sizes <- sapply(module_sets, length)

build_ranks <- function(df, col) {
  vals <- df[[col]]
  names(vals) <- df$gene
  vals <- vals[!is.na(vals)]
  sort(vals, decreasing = TRUE)
}
ranks_TY <- build_ranks(combined, "t_Training_Young")
ranks_TO <- build_ranks(combined, "t_Training_Old")
ranks_Aging <- build_ranks(combined, "t_Aging")

run_module_fgsea <- function(ranks, module_sets) {
  set.seed(42) # fgseaMultilevel is permutation-based; seed per call for reproducible NES/p
  res <- fgsea::fgseaMultilevel(
    pathways = module_sets, stats = ranks,
    minSize = 15, maxSize = 500, nPermSimple = 10000, eps = 0
  )
  as.data.frame(res)
}
fgsea_TY <- run_module_fgsea(ranks_TY, module_sets)
fgsea_TO <- run_module_fgsea(ranks_TO, module_sets)
fgsea_Aging <- run_module_fgsea(ranks_Aging, module_sets)

merge_fgsea <- function(res, suffix) {
  res |>
    as_tibble() |>
    dplyr::select(pathway, NES, padj, size) |>
    dplyr::rename_with(~ paste0(., "_", suffix), c(NES, padj, size))
}
fgsea_wide <- merge_fgsea(fgsea_TY, "TY") |>
  left_join(merge_fgsea(fgsea_TO, "TO"), by = "pathway") |>
  left_join(merge_fgsea(fgsea_Aging, "Aging"), by = "pathway") |>
  mutate(module_color = pathway, n_proteins = mod_sizes[pathway])
fgsea_wide <- fgsea_wide |>
  left_join(mod_bio |> dplyr::select(module_color, module_id, bio_label), by = "module_color") |>
  mutate(
    bio_label = ifelse(is.na(bio_label), stringr::str_to_title(module_color), bio_label),
    module_id = ifelse(is.na(module_id), toupper(substr(module_color, 1, 1)), module_id),
    bio_label = stringr::str_wrap(bio_label, width = 12),
    sig_conc = case_when(
      !is.na(padj_TY) & padj_TY < 0.05 & !is.na(padj_TO) & padj_TO < 0.05 ~ "Both",
      !is.na(padj_TY) & padj_TY < 0.05 ~ "Young only",
      !is.na(padj_TO) & padj_TO < 0.05 ~ "Old only", TRUE ~ "NS"
    ),
    sig_rev = case_when(
      !is.na(padj_Aging) & padj_Aging < 0.05 & !is.na(padj_TO) & padj_TO < 0.05 ~ "Both",
      !is.na(padj_Aging) & padj_Aging < 0.05 ~ "Aging only",
      !is.na(padj_TO) & padj_TO < 0.05 ~ "Training only", TRUE ~ "NS"
    )
  )

build_scatter <- function(df, x_col, y_col, x_lab, y_lab, title, quad_labels, sig_col) {
  x_vals <- df[[x_col]]
  y_vals <- df[[y_col]]
  nes_lim <- max(abs(c(x_vals, y_vals)), na.rm = TRUE) * 1.35

  sp <- cor.test(x_vals, y_vals, method = "spearman")
  subtitle_txt <- sprintf(
    "n = %d modules | \u03c1 = %.2f%s",
    nrow(df), sp$estimate,
    ifelse(sp$p.value < 0.001, ", p < 0.001", sprintf(", p = %.3f", sp$p.value))
  )

  q_tr <- sum(x_vals > 0 & y_vals > 0, na.rm = TRUE)
  q_bl <- sum(x_vals < 0 & y_vals < 0, na.rm = TRUE)
  q_tl <- sum(x_vals < 0 & y_vals > 0, na.rm = TRUE)
  q_br <- sum(x_vals > 0 & y_vals < 0, na.rm = TRUE)

  df$label_col <- ifelse(vapply(df$module_color, is_light_color, logical(1)), "black", "white")
  # green and pink are visually light despite is_light_color returning FALSE
  df$label_col[df$module_color == "green"] <- "black"
  df$label_col[df$module_color == "pink"] <- "black"
  df$label_col[df$module_color == "magenta"] <- "white"

  txt_lab <- 5.0
  txt_quad <- txt_lab

  ggplot(df, aes(x = .data[[x_col]], y = .data[[y_col]])) +
    annotate("rect",
      xmin = 0, xmax = Inf, ymin = 0, ymax = Inf,
      fill = quad_labels$fill[1], alpha = 0.20
    ) +
    annotate("rect",
      xmin = -Inf, xmax = 0, ymin = -Inf, ymax = 0,
      fill = quad_labels$fill[2], alpha = 0.20
    ) +
    annotate("rect",
      xmin = 0, xmax = Inf, ymin = -Inf, ymax = 0,
      fill = quad_labels$fill[3], alpha = 0.20
    ) +
    annotate("rect",
      xmin = -Inf, xmax = 0, ymin = 0, ymax = Inf,
      fill = quad_labels$fill[4], alpha = 0.20
    ) +
    geom_hline(yintercept = 0, color = "grey60", linewidth = 0.2) +
    geom_vline(xintercept = 0, color = "grey60", linewidth = 0.2) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "black", linewidth = 0.3) +
    geom_point(aes(size = n_proteins),
      fill = df$module_color, color = "black",
      shape = 21, alpha = 0.85, stroke = 0.5
    ) +
    # Compact M# identifier on each bubble (drop-cluttered inline pathway names);
    # full module->pathway mapping lives in the color-block key below the panel.
    geom_text(aes(label = module_id),
      color = df$label_col,
      size = txt_lab * 0.65, fontface = "bold",
      show.legend = FALSE
    ) +
    annotate("label",
      x = Inf, y = Inf,
      label = sprintf("%s  n=%d", quad_labels$label[1], q_tr),
      hjust = 1, vjust = 1, size = txt_quad, fontface = "bold",
      color = quad_labels$text_col[1], fill = scales::alpha("white", 0.92),
      label.padding = unit(2, "pt")
    ) +
    annotate("label",
      x = -Inf, y = -Inf,
      label = sprintf("%s  n=%d", quad_labels$label[2], q_bl),
      hjust = 0, vjust = 0, size = txt_quad, fontface = "bold",
      color = quad_labels$text_col[2], fill = scales::alpha("white", 0.92),
      label.padding = unit(2, "pt")
    ) +
    annotate("label",
      x = Inf, y = -Inf,
      label = sprintf("%s  n=%d", quad_labels$label[3], q_br),
      hjust = 1, vjust = 0, size = txt_quad, fontface = "bold",
      color = quad_labels$text_col[3], fill = scales::alpha("white", 0.92),
      label.padding = unit(2, "pt")
    ) +
    annotate("label",
      x = -Inf, y = Inf,
      label = sprintf("%s  n=%d", quad_labels$label[4], q_tl),
      hjust = 0, vjust = 1, size = txt_quad, fontface = "bold",
      color = quad_labels$text_col[4], fill = scales::alpha("white", 0.92),
      label.padding = unit(2, "pt")
    ) +
    scale_size_continuous(
      range = c(3, 9), name = "Proteins",
      breaks = c(50, 100, 200, 300)
    ) +
    scale_x_continuous(expand = expansion(mult = 0.02)) +
    scale_y_continuous(expand = expansion(mult = 0.02)) +
    coord_fixed(ratio = 1, xlim = c(-nes_lim, nes_lim), ylim = c(-nes_lim, nes_lim), clip = "off") +
    labs(title = NULL, subtitle = NULL, x = x_lab, y = y_lab) +
    FIG_THEME +
    theme(
      axis.text = element_text(size = 11, face = "bold", color = "grey30"),
      axis.title.x = element_text(size = 12, face = "bold"),
      axis.title.y = element_text(size = 12, face = "bold", margin = margin(r = -4)),
      legend.position = "none",
      plot.margin = margin(2, 4, 2, -8)
    )
}

quad_conc <- list(
  label = c("Concordant Up", "Concordant Down", "Discordant", "Discordant"),
  fill = c(AGE_COLORS["Old"], AGE_COLORS["Old"], AGE_COLORS["Young"], AGE_COLORS["Young"]),
  text_col = c("#D6604D", "#D6604D", "#4393C3", "#4393C3")
)

p_top <- build_scatter(
  fgsea_wide, "NES_TY", "NES_TO",
  "NES (Training Young)", "NES (Training Old)",
  "Training Concordance", quad_conc, "sig_conc"
)

# Aging axis is cross-sectional (Old vs Young), so training cannot literally
# "reverse" it within subjects. Quadrants describe directional overlap only:
# same sign = tracks the aging signature, opposite sign = counters it.
quad_rev <- list(
  label = c("Tracks aging", "Tracks aging", "Counters aging", "Counters aging"),
  fill = c(AGE_COLORS["Old"], AGE_COLORS["Old"], "#2E7D32", "#2E7D32"),
  text_col = c("#D6604D", "#D6604D", "#4393C3", "#4393C3")
)

p_bottom <- build_scatter(
  fgsea_wide, "NES_Aging", "NES_TO",
  "NES (Aging, cross-sectional)", "NES (Training Old)",
  "Directional overlap with the aging signature", quad_rev, "sig_rev"
)

# Module color key: data-driven; lists exactly the 9 modules in the scatter
# in canonical M# order (drop-clutter equivalent of F07 Panel A's key strip).
mod_key_df <- mod_bio |>
  dplyr::filter(module_color %in% fgsea_wide$module_color) |>
  dplyr::arrange(match(module_id, sprintf("M%d", 1:9)))

build_key_cell <- function(mod, mid, lab) {
  ggplot() +
    annotate("rect",
      xmin = 0.02, xmax = 0.10, ymin = 0.25, ymax = 0.75,
      fill = mod, color = "grey20", linewidth = 0.4
    ) +
    annotate("text",
      x = 0.06, y = 0.50, label = mid,
      color = ifelse(is_light_color(mod) | mod %in% c("green", "pink"), "black", "white"),
      fontface = "bold", size = 3.2
    ) +
    annotate("text",
      x = 0.13, y = 0.50, label = lab,
      hjust = 0, vjust = 0.5, size = 3.0, color = "grey10"
    ) +
    scale_x_continuous(limits = c(0, 1), expand = c(0, 0)) +
    scale_y_continuous(limits = c(0, 1), expand = c(0, 0)) +
    coord_cartesian(clip = "off") +
    theme_void() +
    theme(plot.margin = margin_auto(4, 1))
}
key_cells <- purrr::pmap(
  list(mod_key_df$module_color, mod_key_df$module_id, mod_key_df$bio_label),
  build_key_cell
)
# 3 cols x ceil(n/3) rows -> ~73 mm per cell at 220 mm wide,
# matching F07 panel A key cell density (prevents label overlap).
key_ncol <- 3L
key_nrow <- ceiling(length(key_cells) / key_ncol)
module_key_strip <- wrap_plots(key_cells, nrow = key_nrow, ncol = key_ncol)

scatters_panel <- (p_top / p_bottom / module_key_strip) +
  plot_layout(heights = c(1, 1, 0.32))

PB_W <- 220
PB_H <- 320 # +50 mm for the 3-row key strip

ggsave(file.path(RPT_PNG, "MAIN_panel_B_scatters.png"), scatters_panel,
  width = PB_W, height = PB_H, units = "mm", dpi = 300
)
ggsave(file.path(RPT_PDF, "MAIN_panel_B_scatters.pdf"), scatters_panel,
  width = PB_W, height = PB_H, units = "mm", device = get_pdf_device()
)

p_legend_src <- p_top +
  scale_size_continuous(
    range = c(3, 9), name = "Proteins",
    breaks = c(50, 100, 200, 300),
    guide = guide_legend(
      nrow = 1,
      override.aes = list(alpha = 0.7, fill = "grey60", stroke = 0)
    )
  ) +
  theme_sub_legend(
    position = "bottom",
    title = element_text(size = 12, face = "bold"),
    text = element_text(size = 11),
    key = element_rect(fill = NA, color = NA),
    key.size = unit(5, "mm"),
    background = element_rect(fill = NA, color = NA)
  )
legend_grob <- cowplot::get_plot_component(p_legend_src, "guide-box-bottom", return_all = FALSE)
p_legend <- cowplot::ggdraw(legend_grob)

ggsave(file.path(RPT_PNG, "MAIN_panel_B_legend.png"), p_legend,
  width = 90, height = 16, units = "mm", dpi = 300
)

write_csv(fgsea_wide, file.path(DAT, "panel_B_module_fgsea.csv"))

message("Panel B: scatters + legend saved")
