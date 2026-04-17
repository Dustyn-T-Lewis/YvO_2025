# Figure 6 — Panel B: Module-Level NES Scatter (fGSEA, two-panel)
#
# Left:  Training_Young vs Training_Old  — concordance/blunting
# Right: Aging vs Training_Old           — reversal
#
# Uses WGCNA modules as custom gene sets for fGSEA, ranked by limma t-statistics.
# Points colored by module color, sized by protein count, labeled with bio_label.
#
# Generates:
#   b_reports/panels/MAIN_panel_B_concordance.png
#   b_reports/panels/MAIN_panel_C_reversal.png
#   c_data/panel_B_module_fgsea.csv

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

RPT <- "04_Figures/F06/b_reports/panels"
DAT <- "04_Figures/F06/c_data"
dir.create(RPT, recursive = TRUE, showWarnings = FALSE)

# =============================================================================
# 1. Load data
# =============================================================================

# Module assignments → gene sets
module_df <- read_csv(file.path(DAT, "wgcna/wgcna_module_assignments.csv"),
                      show_col_types = FALSE)
mod_bio   <- read_csv(file.path(DAT, "mod_bio_labels.csv"), show_col_types = FALSE)

# Combined DEP results (t-statistics)
combined <- read_csv("03_DEP/c_data/03_combined_results.csv", show_col_types = FALSE)

# =============================================================================
# 2. Build module gene sets (exclude grey)
# =============================================================================

module_sets <- module_df %>%
  filter(module_color != "grey", !is.na(gene), gene != "") %>%
  split(.$module_color) %>%
  lapply(function(x) x$gene)

mod_sizes <- sapply(module_sets, length)
message(sprintf("  %d modules, sizes: %s",
                length(module_sets),
                paste(sprintf("%s=%d", names(mod_sizes), mod_sizes), collapse = ", ")))

# =============================================================================
# 3. Build t-statistic rankings for 3 contrasts
# =============================================================================

build_ranks <- function(df, col) {
  vals <- df[[col]]
  names(vals) <- df$gene
  vals <- vals[!is.na(vals)]
  sort(vals, decreasing = TRUE)
}

ranks_TY    <- build_ranks(combined, "t_Training_Young")
ranks_TO    <- build_ranks(combined, "t_Training_Old")
ranks_Aging <- build_ranks(combined, "t_Aging")

message(sprintf("  Rank vectors: TY=%d, TO=%d, Aging=%d genes",
                length(ranks_TY), length(ranks_TO), length(ranks_Aging)))

# =============================================================================
# 4. Run fGSEA per contrast
# =============================================================================

run_module_fgsea <- function(ranks, module_sets) {
  res <- fgsea::fgseaMultilevel(
    pathways    = module_sets,
    stats       = ranks,
    minSize     = 15,
    maxSize     = 500,
    nPermSimple = 10000,
    eps         = 0
  )
  res <- as.data.frame(res)
  res$padj <- p.adjust(res$pval, method = "BH")
  res
}

fgsea_TY    <- run_module_fgsea(ranks_TY, module_sets)
fgsea_TO    <- run_module_fgsea(ranks_TO, module_sets)
fgsea_Aging <- run_module_fgsea(ranks_Aging, module_sets)

message(sprintf("  fGSEA: TY %d sig, TO %d sig, Aging %d sig (padj<0.05)",
                sum(fgsea_TY$padj < 0.05),
                sum(fgsea_TO$padj < 0.05),
                sum(fgsea_Aging$padj < 0.05)))

# =============================================================================
# 5. Merge to wide format
# =============================================================================

merge_fgsea <- function(res, suffix) {
  res %>%
    as_tibble() %>%
    dplyr::select(pathway, NES, padj, size) %>%
    dplyr::rename_with(~ paste0(., "_", suffix), c(NES, padj, size))
}

fgsea_wide <- merge_fgsea(fgsea_TY, "TY") %>%
  left_join(merge_fgsea(fgsea_TO, "TO"), by = "pathway") %>%
  left_join(merge_fgsea(fgsea_Aging, "Aging"), by = "pathway") %>%
  mutate(
    module_color = pathway,
    n_proteins   = mod_sizes[pathway]
  )

# Add bio labels
fgsea_wide <- fgsea_wide %>%
  left_join(mod_bio %>% dplyr::select(module_color, bio_label), by = "module_color") %>%
  mutate(bio_label = ifelse(is.na(bio_label),
                            stringr::str_to_title(module_color),
                            bio_label))

# Significance flags
fgsea_wide <- fgsea_wide %>%
  mutate(
    sig_TY    = !is.na(padj_TY) & padj_TY < 0.05,
    sig_TO    = !is.na(padj_TO) & padj_TO < 0.05,
    sig_Aging = !is.na(padj_Aging) & padj_Aging < 0.05,
    # Shape encoding for concordance panel
    sig_conc = case_when(
      sig_TY & sig_TO ~ "Both",
      sig_TY          ~ "Young only",
      sig_TO          ~ "Old only",
      TRUE            ~ "NS"
    ),
    # Shape encoding for reversal panel
    sig_rev = case_when(
      sig_Aging & sig_TO ~ "Both",
      sig_Aging          ~ "Aging only",
      sig_TO             ~ "Training only",
      TRUE               ~ "NS"
    )
  )

# =============================================================================
# 6. Export data
# =============================================================================

fgsea_wide %>%
  transmute(
    module = module_color, bio_label, n_proteins,
    NES_Training_Young = round(NES_TY, 3),
    NES_Training_Old   = round(NES_TO, 3),
    NES_Aging          = round(NES_Aging, 3),
    padj_Training_Young = signif(padj_TY, 4),
    padj_Training_Old   = signif(padj_TO, 4),
    padj_Aging          = signif(padj_Aging, 4)
  ) %>%
  arrange(module) %>%
  write_csv(file.path(DAT, "panel_B_module_fgsea.csv"))

# =============================================================================
# 7. Scatter plot helper
# =============================================================================

build_module_scatter <- function(df, x_col, y_col, x_lab, y_lab,
                                 title, quad_labels, sig_col) {
  x_vals <- df[[x_col]]
  y_vals <- df[[y_col]]

  # Symmetric axis limits
  nes_lim <- max(abs(c(x_vals, y_vals)), na.rm = TRUE) * 1.25

  # Spearman correlation
  sp <- cor.test(x_vals, y_vals, method = "spearman")
  sp_ci <- fisher_z_ci(sp$estimate, nrow(df))
  rho_label <- sprintf("\u03c1 = %.2f, %s",
                       sp$estimate,
                       ifelse(sp$p.value < 0.001, "p < 0.001",
                              sprintf("p = %.3f", sp$p.value)))

  # Panel subtitle (harmonized with F04/F05 panel B sibling subtitles)
  subtitle_txt <- sprintf(
    "WGCNA modules (n = %d) | \u03c1 = %.2f%s | fGSEA on module-member t-stat ranks\n(distinct from ME\u2013trait r shown in panel A)",
    nrow(df), sp$estimate,
    ifelse(sp$p.value < 0.001, ", p < 0.001", sprintf(", p = %.3f", sp$p.value))
  )

  # Quadrant counts
  q_tr <- sum(x_vals > 0 & y_vals > 0, na.rm = TRUE)
  q_bl <- sum(x_vals < 0 & y_vals < 0, na.rm = TRUE)
  q_tl <- sum(x_vals < 0 & y_vals > 0, na.rm = TRUE)
  q_br <- sum(x_vals > 0 & y_vals < 0, na.rm = TRUE)

  # Significance shape: filled circle = sig for axis contrast, open = NS
  sig_status <- df[[sig_col]]
  df$pt_shape <- ifelse(sig_status == "Both", 21,
                  ifelse(sig_status == "NS", 1, 21))
  df$pt_alpha <- ifelse(sig_status == "NS", 0.6, 0.9)
  df$pt_stroke <- ifelse(sig_status == "NS", 1.0, 0.6)

  # Black text for light modules, white for dark
  df$label_col <- ifelse(sapply(df$module_color, is_light_color), "black", "white")

  txt_quad <- scale_text(BASE_QUADRANT, 200) * 0.85
  txt_lab  <- scale_text(BASE_GENE, 200) * 1.5

  p <- ggplot(df, aes(x = .data[[x_col]], y = .data[[y_col]])) +
    # Quadrant shading
    annotate("rect", xmin = 0, xmax = Inf,  ymin = 0, ymax = Inf,
             fill = quad_labels$fill[1], alpha = 0.20) +
    annotate("rect", xmin = -Inf, xmax = 0, ymin = -Inf, ymax = 0,
             fill = quad_labels$fill[2], alpha = 0.20) +
    annotate("rect", xmin = 0, xmax = Inf,  ymin = -Inf, ymax = 0,
             fill = quad_labels$fill[3], alpha = 0.20) +
    annotate("rect", xmin = -Inf, xmax = 0, ymin = 0, ymax = Inf,
             fill = quad_labels$fill[4], alpha = 0.20) +
    # Reference lines
    geom_hline(yintercept = 0, color = "grey60", linewidth = 0.2) +
    geom_vline(xintercept = 0, color = "grey60", linewidth = 0.2) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed",
                color = "black", linewidth = 0.3) +
    # Module points
    geom_point(aes(size = n_proteins),
               fill = df$module_color, color = "black",
               shape = 21, alpha = df$pt_alpha, stroke = df$pt_stroke) +
    # Labels
    geom_label_repel(
      aes(label = bio_label),
      fill = scales::alpha(df$module_color, 0.85),
      color = df$label_col,
      size = txt_lab, fontface = "bold",
      max.overlaps = 20,
      segment.size = 0.4, segment.color = "grey40",
      min.segment.length = 0,
      box.padding = 0.9, point.padding = 0.6,
      force = 15, force_pull = 0.25,
      label.padding = unit(2, "pt"),
      label.r = unit(1.5, "pt"),
      label.size = 0.2, seed = 42,
      show.legend = FALSE
    ) +
    # Quadrant labels (F04/F05 style)
    annotate("label", x = nes_lim, y = nes_lim,
             label = sprintf("%s  n = %d", quad_labels$label[1], q_tr),
             hjust = 1, vjust = 1, size = txt_quad, fontface = "bold",
             color = quad_labels$text_col[1],
             fill = scales::alpha("white", 0.92),
             label.padding = unit(2.5, "pt")) +
    annotate("label", x = -nes_lim, y = -nes_lim,
             label = sprintf("%s  n = %d", quad_labels$label[2], q_bl),
             hjust = 0, vjust = 0, size = txt_quad, fontface = "bold",
             color = quad_labels$text_col[2],
             fill = scales::alpha("white", 0.92),
             label.padding = unit(2.5, "pt")) +
    annotate("label", x = nes_lim, y = -nes_lim,
             label = sprintf("%s  n = %d", quad_labels$label[3], q_br),
             hjust = 1, vjust = 0, size = txt_quad, fontface = "bold",
             color = quad_labels$text_col[3],
             fill = scales::alpha("white", 0.92),
             label.padding = unit(2.5, "pt")) +
    annotate("label", x = -nes_lim, y = nes_lim,
             label = sprintf("%s  n = %d", quad_labels$label[4], q_tl),
             hjust = 0, vjust = 1, size = txt_quad, fontface = "bold",
             color = quad_labels$text_col[4],
             fill = scales::alpha("white", 0.92),
             label.padding = unit(2.5, "pt")) +
    scale_size_continuous(range = c(4, 12), name = "Proteins",
                          breaks = c(50, 100, 200, 350)) +
    scale_x_continuous(expand = expansion(mult = 0.02)) +
    scale_y_continuous(expand = expansion(mult = 0.02)) +
    coord_fixed(ratio = 1,
                xlim = c(-nes_lim, nes_lim),
                ylim = c(-nes_lim, nes_lim)) +
    labs(title = title, subtitle = subtitle_txt, x = x_lab, y = y_lab) +
    FIG_THEME +
    theme(
      axis.text     = element_text(size = 16, face = "bold", color = "grey30"),
      axis.title    = element_text(size = 17, face = "bold"),
      legend.position = "none",
      # Title/subtitle line-spacing harmonized with panel A: bold 22pt title
      # with a ~3 mm gap to the italic 15pt subtitle, no extra bottom pad.
      plot.title    = element_text(size = composite_text_sizes(535)$title, face = "bold", hjust = 0,
                                   margin = margin(b = 3),
                                   lineheight = 1.05),
      plot.subtitle = element_text(size = composite_text_sizes(535)$subtitle, face = "bold.italic",
                                   color = "grey40", hjust = 0,
                                   margin = margin(b = 6),
                                   lineheight = 1.15),
      plot.margin   = margin(-8, 4, 2, 14)
    )

  p
}

# =============================================================================
# 8. Build two panels
# =============================================================================

message("Panel B: module NES scatter (two-panel)...")

# Left: Training concordance (TY vs TO)
quad_conc <- list(
  label    = c("Concordant Up", "Concordant Down", "Discordant", "Discordant"),
  fill     = c(AGE_COLORS["Old"], AGE_COLORS["Old"], AGE_COLORS["Young"], AGE_COLORS["Young"]),
  text_col = c("#D6604D", "#D6604D", "#4393C3", "#4393C3")
)

p_left <- build_module_scatter(
  fgsea_wide, "NES_TY", "NES_TO",
  x_lab = "NES (Training Young)",
  y_lab = "NES (Training Old)",
  title = "Training Concordance",
  quad_labels = quad_conc,
  sig_col = "sig_conc"
)

# Right: Reversal (Aging vs TO)
quad_rev <- list(
  label    = c("Exacerbated", "Exacerbated", "Reversed", "Reversed"),
  fill     = c(AGE_COLORS["Old"], AGE_COLORS["Old"], "#2E7D32", "#2E7D32"),
  text_col = c("#D6604D", "#D6604D", "#4393C3", "#4393C3")
)

p_right <- build_module_scatter(
  fgsea_wide, "NES_Aging", "NES_TO",
  x_lab = "NES (Aging)",
  y_lab = "NES (Training Old)",
  title = "Aging Reversal",
  quad_labels = quad_rev,
  sig_col = "sig_rev"
)

# =============================================================================
# 9. Save individual panels + composite
# =============================================================================

PG_W <- 260   # single scatter width (for composite: two side-by-side = 520mm)
PG_H <- 310   # enlarged from 250 so scatters dominate the lower row

# Individual panels (PNG only, for composite stitch)
ggsave(file.path(RPT, "MAIN_panel_B_concordance.png"), p_left,
       width = PG_W, height = PG_H, units = "mm",
       dpi = 300, limitsize = FALSE)

ggsave(file.path(RPT, "MAIN_panel_C_reversal.png"), p_right,
       width = PG_W, height = PG_H, units = "mm",
       dpi = 300, limitsize = FALSE)

# Shared size legend for B/C — built from p_left and stripped down.
p_legend_src <- p_left +
  scale_size_continuous(range = c(4, 12), name = "Proteins",
                        breaks = c(50, 100, 200, 350),
                        guide = guide_legend(nrow = 1,
                                             override.aes = list(alpha = 0.7,
                                                                 fill = "grey60"))) +
  theme(legend.position = "bottom",
        legend.title    = element_text(size = 13, face = "bold"),
        legend.text     = element_text(size = 12),
        legend.key.size = unit(6, "mm"))
legend_grob <- cowplot::get_plot_component(p_legend_src, "guide-box-bottom",
                                           return_all = FALSE)
p_BC_legend <- cowplot::ggdraw(legend_grob)

ggsave(file.path(RPT, "MAIN_panel_BC_size_legend.png"), p_BC_legend,
       width = 80, height = 18, units = "mm",
       dpi = 300, limitsize = FALSE)

message(sprintf("  Panel B (concordance) + Panel C (reversal) saved: %d x %d mm each",
                PG_W, PG_H))
message("  Shared size legend saved: MAIN_panel_BC_size_legend.png")
