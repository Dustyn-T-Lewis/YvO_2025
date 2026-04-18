# NOTE: Stratified LMM columns were absorbed into panel_A.R (2026-04-10).
# This script still generates the standalone SUPP visualization
# (SUPP_age_stratified_lmm) consumed by 90_stitch_analysis.R.
#
# Original description:
# Supplementary F: Age-Stratified LMM + Phenotype Heatmap (Panel-A-styled)
#
# Companion to Panel A. Four fixed LMM columns + dynamic phenotype columns,
# flanked by gene count bars (left) and preservation Zsummary bars (right).
#
# Fixed columns (4):
#   Training_Young / Training_Old  — from wgcna_lmm_stratified_audit.csv
#   Aging_Pre / Aging_Post         — unpaired t-tests (Young vs Old at each timepoint)
#
# Dynamic columns:
#   Baseline and change phenotype correlations, stratified by age group.
#   Included only if >= 1 cell has BH p < 0.05 across 10 modules.
#
# Reads:  F06/c_data/MEs.rds, meta.csv, mod_bio_labels.csv,
#         wgcna/wgcna_module_assignments.csv, wgcna/wgcna_lmm_stratified_audit.csv,
#         wgcna/wgcna_{baseline,change}_trait_{correlations,pvalues_bh}_{young,old}.csv,
#         05_panel_E_preservation.csv
# Writes: b_reports/supp/02_analysis/SUPP_age_stratified_lmm.png
#         c_data/supp_F_lmm_stratified_audit.csv

setwd(rprojroot::find_rstudio_root_file())
source("04_Figures/shared/style.R")

library(readr)
library(dplyr)
library(tidyr)
library(tibble)
library(stringr)
library(patchwork)

DAT     <- "04_Figures/F06/c_data"
RPT_PNG <- "04_Figures/F06/b_reports/supp/02_analysis/png/panels"
RPT_PDF <- "04_Figures/F06/b_reports/supp/02_analysis/pdf/panels"
dir.create(RPT_PNG, recursive = TRUE, showWarnings = FALSE)
dir.create(RPT_PDF, recursive = TRUE, showWarnings = FALSE)

pdf_device <- get_pdf_device()

# =============================================================================
# 1. Data loading
# =============================================================================

MEs       <- readRDS(file.path(DAT, "MEs.rds"))
meta      <- read_csv(file.path(DAT, "meta.csv"), show_col_types = FALSE)
mod_bio   <- read_csv(file.path(DAT, "mod_bio_labels.csv"), show_col_types = FALSE)
module_df <- read_csv(file.path(DAT, "wgcna/wgcna_module_assignments.csv"),
                      show_col_types = FALSE)

meta <- meta %>%
  mutate(time = factor(time, levels = c("Pre", "Post")),
         age  = factor(age,  levels = c("Young", "Old")))

# Display label map
if (!"display_label" %in% colnames(mod_bio)) {
  mod_bio <- mod_bio %>%
    mutate(display_label = paste0(bio_label, " (",
                                  stringr::str_to_title(module_color), ")"))
}
mod_display_vec <- setNames(mod_bio$display_label,
                            paste0("ME", mod_bio$module_color))

# Module ordering: ascending protein count (matches Panel A)
mod_size <- module_df %>%
  filter(module_color != "grey") %>%
  count(module_color, name = "n_proteins") %>%
  mutate(module = paste0("ME", module_color))
mod_order <- mod_size$module[order(mod_size$n_proteins)]

all_mods <- colnames(MEs)
non_grey <- mod_order[mod_order %in% all_mods]

message(sprintf("Supp F: age-stratified heatmap on %d modules...", length(non_grey)))

# =============================================================================
# 2. Fixed LMM columns: Training_Young, Training_Old (from audit CSV)
# =============================================================================

strat_audit <- read_csv(file.path(DAT, "wgcna/wgcna_lmm_stratified_audit.csv"),
                        show_col_types = FALSE)

trn_young <- strat_audit %>% filter(age_group == "Young") %>%
  dplyr::select(module, r_equiv, p_bh) %>%
  mutate(col_name = "Training_Young")
trn_old <- strat_audit %>% filter(age_group == "Old") %>%
  dplyr::select(module, r_equiv, p_bh) %>%
  mutate(col_name = "Training_Old")

# =============================================================================
# 3. Fixed LMM columns: Aging_Pre, Aging_Post (unpaired t-tests)
# =============================================================================

compute_aging_ttest <- function(me_data, meta_df, timepoint_label) {
  sub_meta <- meta_df %>% filter(time == timepoint_label)
  rows <- list()
  for (mod in non_grey) {
    me_vals <- me_data[sub_meta$sample_id, mod]
    young_vals <- me_vals[sub_meta$age == "Young"]
    old_vals   <- me_vals[sub_meta$age == "Old"]
    tt <- t.test(old_vals, young_vals)  # Old - Young direction (positive = higher in Old)
    df_val <- tt$parameter
    t_val  <- tt$statistic
    r_eq   <- sign(t_val) * sqrt(t_val^2 / (t_val^2 + df_val))
    rows[[mod]] <- tibble(
      module    = mod,
      estimate  = round(tt$estimate[1] - tt$estimate[2], 5),
      t_stat    = round(as.numeric(t_val), 4),
      df        = round(as.numeric(df_val), 2),
      p_raw     = tt$p.value,
      r_equiv   = round(as.numeric(r_eq), 4)
    )
  }
  bind_rows(rows)
}

aging_pre  <- compute_aging_ttest(MEs, meta, "Pre")
aging_post <- compute_aging_ttest(MEs, meta, "Post")

# BH correction across 10 modules per timepoint
aging_pre$p_bh  <- p.adjust(aging_pre$p_raw, method = "BH")
aging_post$p_bh <- p.adjust(aging_post$p_raw, method = "BH")

aging_pre$col_name  <- "Aging_Pre"
aging_post$col_name <- "Aging_Post"

message(sprintf("  Aging_Pre:  %d modules, %d with BH<0.05",
                nrow(aging_pre), sum(aging_pre$p_bh < 0.05)))
message(sprintf("  Aging_Post: %d modules, %d with BH<0.05",
                nrow(aging_post), sum(aging_post$p_bh < 0.05)))

# =============================================================================
# 4. Dynamic phenotype columns (filtered by >= 1 BH < 0.05 cell)
# =============================================================================

load_pheno_section <- function(cor_file, pval_file, suffix) {
  cor_mat  <- read_csv(file.path(DAT, "wgcna", cor_file),
                       show_col_types = FALSE) %>%
    column_to_rownames("module") %>% as.matrix()
  pval_mat <- read_csv(file.path(DAT, "wgcna", pval_file),
                       show_col_types = FALSE) %>%
    column_to_rownames("module") %>% as.matrix()

  # Filter to non-grey modules present in both
  shared <- intersect(non_grey, rownames(cor_mat))
  cor_mat  <- cor_mat[shared, , drop = FALSE]
  pval_mat <- pval_mat[shared, , drop = FALSE]

  # Keep only columns with >= 1 BH < 0.05
  keep <- apply(pval_mat, 2, function(p) any(p < 0.05, na.rm = TRUE))
  if (!any(keep)) return(NULL)

  cor_mat  <- cor_mat[, keep, drop = FALSE]
  pval_mat <- pval_mat[, keep, drop = FALSE]

  # Suffix column names
  colnames(cor_mat)  <- paste0(colnames(cor_mat), "_", suffix)
  colnames(pval_mat) <- paste0(colnames(pval_mat), "_", suffix)

  list(cor = cor_mat, pval = pval_mat)
}

bl_young <- load_pheno_section("wgcna_baseline_trait_correlations_young.csv",
                               "wgcna_baseline_trait_pvalues_bh_young.csv", "Y")
bl_old   <- load_pheno_section("wgcna_baseline_trait_correlations_old.csv",
                               "wgcna_baseline_trait_pvalues_bh_old.csv", "O")
ch_young <- load_pheno_section("wgcna_change_trait_correlations_young.csv",
                               "wgcna_change_trait_pvalues_bh_young.csv", "Y")
ch_old   <- load_pheno_section("wgcna_change_trait_correlations_old.csv",
                               "wgcna_change_trait_pvalues_bh_old.csv", "O")

# =============================================================================
# 5. Assemble combined r and p matrices
# =============================================================================

# Fixed columns: r and p matrices (rows = non_grey)
fixed_r <- matrix(NA, nrow = length(non_grey), ncol = 4,
                  dimnames = list(non_grey, c("Training_Young", "Training_Old",
                                              "Aging_Pre", "Aging_Post")))
fixed_p <- fixed_r

for (src in list(trn_young, trn_old, aging_pre, aging_post)) {
  cn <- src$col_name[1]
  for (i in seq_len(nrow(src))) {
    if (src$module[i] %in% non_grey) {
      fixed_r[src$module[i], cn] <- src$r_equiv[i]
      fixed_p[src$module[i], cn] <- src$p_bh[i]
    }
  }
}

# Concatenate with dynamic phenotype columns
pheno_sections <- list(bl_young, bl_old, ch_young, ch_old)
pheno_sections <- pheno_sections[!sapply(pheno_sections, is.null)]

if (length(pheno_sections) > 0) {
  pheno_r <- do.call(cbind, lapply(pheno_sections, function(s) s$cor[non_grey, , drop = FALSE]))
  pheno_p <- do.call(cbind, lapply(pheno_sections, function(s) s$pval[non_grey, , drop = FALSE]))
  cor_mat  <- cbind(fixed_r, pheno_r)
  pval_mat <- cbind(fixed_p, pheno_p)
} else {
  cor_mat  <- fixed_r
  pval_mat <- fixed_p
}

n_fixed <- 4
n_pheno <- ncol(cor_mat) - n_fixed
n_cols  <- ncol(cor_mat)

message(sprintf("  Total columns: %d fixed + %d phenotype = %d",
                n_fixed, n_pheno, n_cols))

# =============================================================================
# 6. Build heatmap data frame
# =============================================================================

trait_order <- colnames(cor_mat)

# Column labels: readable short names
col_labels <- c(
  Training_Young = "Trn_Y",
  Training_Old   = "Trn_O",
  Aging_Pre      = "Age_Pre",
  Aging_Post     = "Age_Post"
)
# Dynamic phenotype labels: strip suffix and add age marker
for (cn in trait_order[!trait_order %in% names(col_labels)]) {
  base <- sub("_(Y|O)$", "", cn)
  suffix <- sub(".*_(Y|O)$", "\\1", cn)
  # Pretty-print common phenotype names
  pretty <- gsub("_", " ", base)
  pretty <- sub("^delta ", "\u0394", pretty)
  pretty <- sub("^BMI Pre$", "BMI", pretty)
  pretty <- sub("^VL Pre$", "VL", pretty)
  pretty <- sub("^LBM Pre$", "LBM", pretty)
  col_labels[cn] <- paste0(pretty, "\n(", suffix, ")")
}

# X positions with gap between fixed LMM and phenotype sections
gap <- 0.4
xpos <- seq_len(n_fixed)
if (n_pheno > 0) {
  xpos <- c(xpos, n_fixed + gap + seq_len(n_pheno))
}
trait_xpos <- setNames(xpos, trait_order)

heat_df <- expand.grid(module = non_grey, trait = trait_order,
                       stringsAsFactors = FALSE) %>%
  mutate(cor      = as.vector(cor_mat[non_grey, ]),
         pval     = as.vector(pval_mat[non_grey, ]),
         stars    = sig_stars(pval),
         label    = sprintf("%.2f%s", cor, stars),
         xpos     = trait_xpos[trait],
         has_black = !is.na(pval) & pval < 0.05)
heat_df$module <- factor(heat_df$module, levels = mod_order)

# Y-axis labels: biological display labels
y_label_vec <- mod_display_vec[mod_order]
miss <- is.na(y_label_vec)
y_label_vec[miss] <- stringr::str_to_title(gsub("^ME", "", mod_order[miss]))

# Module color map for text color decisions
mod_color_raw <- setNames(gsub("^ME", "", mod_order), mod_order)
light_modules <- names(which(sapply(mod_color_raw, is_light_color)))

# =============================================================================
# 7. Save audit CSV
# =============================================================================

audit_out <- bind_rows(
  trn_young %>% dplyr::select(module, col_name, r_equiv, p_bh),
  trn_old   %>% dplyr::select(module, col_name, r_equiv, p_bh),
  aging_pre %>% dplyr::select(module, col_name, r_equiv, p_bh),
  aging_post %>% dplyr::select(module, col_name, r_equiv, p_bh)
)
write_csv(audit_out, file.path(DAT, "supp_F_lmm_stratified_audit.csv"))
message(sprintf("  Audit saved: %s", file.path(DAT, "supp_F_lmm_stratified_audit.csv")))

# =============================================================================
# 8. Figure dimensions (scaled by column count)
# =============================================================================

PA_W <- max(350, 200 + n_cols * 22)
PA_H <- 250

txt_cell  <- scale_text(BASE_GENE, PA_W) * 0.95
txt_count <- scale_text(BASE_COUNT, PA_W) * 0.85
txt_axis  <- scale_text(BASE_STAT, PA_W) * 1.5
txt_ylab  <- scale_text(BASE_GENE, PA_W) * 1.10
txt_brack <- scale_text(BASE_GENE, PA_W) * 1.15

# =============================================================================
# 9. Gene count bars (left flank) — same as Panel A
# =============================================================================

mod_counts <- module_df %>%
  filter(module_color != "grey") %>%
  count(module_color, name = "n_proteins") %>%
  mutate(module = paste0("ME", module_color)) %>%
  filter(module %in% mod_order) %>%
  mutate(module = factor(module, levels = mod_order),
         x_sqrt = sqrt(n_proteins))

# Pathway labels inside bars
pathway_label_map <- setNames(
  mod_bio$bio_label[match(gsub("^ME", "", mod_order), mod_bio$module_color)],
  mod_order
)
mod_counts$pathway_label <- pathway_label_map[as.character(mod_counts$module)]
mod_counts$pathway_label[is.na(mod_counts$pathway_label)] <- ""

# Wrap labels for smaller bars
wrap_mods <- c("magenta", "pink", "black", "green")
for (i in seq_len(nrow(mod_counts))) {
  if (mod_counts$module_color[i] %in% wrap_mods) {
    mod_counts$pathway_label[i] <- sub(" ", "\n", mod_counts$pathway_label[i])
  }
}

# Text color for in-bar pathway names
black_text_mods <- c("yellow", "green", "turquoise", "pink")
mod_counts$bar_text_col <- ifelse(mod_counts$module_color %in% black_text_mods,
                                   "black", "white")

# Color name labels above each bar
color_label_df <- mod_counts %>%
  mutate(color_name = stringr::str_to_title(module_color))

p_counts <- ggplot(mod_counts, aes(x = x_sqrt, y = module)) +
  geom_col(fill = mod_counts$module_color, color = "black",
           linewidth = 0.3, width = 0.65) +
  geom_text(aes(label = pathway_label, x = x_sqrt * 0.95),
            size = txt_count * 0.80, fontface = "bold",
            color = mod_counts$bar_text_col, hjust = 0, lineheight = 0.85) +
  geom_text(data = color_label_df,
            aes(label = color_name, x = x_sqrt / 2, y = module),
            nudge_y = 0.42, size = txt_count * 0.7, fontface = "bold",
            color = "black", hjust = 0.5) +
  scale_x_reverse(expand = c(0, 0),
                  limits = c(max(mod_counts$x_sqrt) * 1.30, 0),
                  breaks = sqrt(c(50, 100, 200, 300)),
                  labels = c(50, 100, 200, 300)) +
  scale_y_discrete(labels = NULL) +
  labs(y = NULL, x = "Protein Counts") +
  FIG_THEME +
  theme(axis.text.y        = element_blank(),
        axis.ticks.y       = element_blank(),
        axis.text.x        = element_text(size = txt_axis),
        axis.title.x       = element_text(size = txt_axis, margin = margin(t = -2)),
        panel.grid.major.y = element_blank(),
        panel.grid.minor   = element_blank(),
        panel.border       = element_blank(),
        axis.line.x        = element_line(color = "black", linewidth = 0.3),
        legend.position    = "none",
        plot.margin        = margin(2, 0, -5, 2))

# =============================================================================
# 10. Section brackets (above heatmap)
# =============================================================================

xmin_all <- min(xpos) - 0.55
xmax_all <- max(xpos) + 0.55

brackets <- tibble(
  label = c("Stratified LMM\n(within-age training)",
            "Unpaired t-test\n(cross-age at timepoint)"),
  start = c(xpos[1], xpos[3]),
  end   = c(xpos[2], xpos[4])
)

# Add phenotype bracket if dynamic columns exist
if (n_pheno > 0) {
  brackets <- bind_rows(brackets, tibble(
    label = "Phenotype Correlations\n(BH < 0.05 filter)",
    start = xpos[n_fixed + 1],
    end   = xpos[n_cols]
  ))
}
brackets <- brackets %>% mutate(mid = (start + end) / 2)

p_brackets <- ggplot() +
  geom_segment(data = brackets,
               aes(x = start - 0.4, xend = end + 0.4, y = 0.28, yend = 0.28),
               linewidth = 0.4, color = "grey30") +
  geom_segment(data = brackets,
               aes(x = start - 0.4, xend = start - 0.4, y = 0.28, yend = 0.14),
               linewidth = 0.4, color = "grey30") +
  geom_segment(data = brackets,
               aes(x = end + 0.4, xend = end + 0.4, y = 0.28, yend = 0.14),
               linewidth = 0.4, color = "grey30") +
  geom_text(data = brackets,
            aes(x = mid, y = 0.62, label = label),
            size = txt_brack, fontface = "bold", color = "grey25",
            lineheight = 0.85) +
  scale_x_continuous(limits = c(xmin_all, xmax_all), expand = c(0, 0)) +
  scale_y_continuous(limits = c(0, 1), expand = c(0, 0)) +
  theme_void() +
  theme(plot.margin = margin(2, 2, 0, 0))

# =============================================================================
# 11. Heatmap
# =============================================================================

p_heat <- ggplot(heat_df, aes(x = xpos, y = module, fill = cor)) +
  # Section background bands
  annotate("rect", xmin = xpos[1] - 0.5, xmax = xpos[2] + 0.5,
           ymin = -Inf, ymax = Inf, fill = "grey70", alpha = 0.06) +
  annotate("rect", xmin = xpos[3] - 0.5, xmax = xpos[4] + 0.5,
           ymin = -Inf, ymax = Inf, fill = "grey70", alpha = 0.06) +
  geom_tile(color = "black", linewidth = 0.3) +
  # Black border for BH FDR < 0.05
  geom_tile(data = heat_df %>% filter(has_black),
            color = "black", linewidth = 1.0, fill = NA) +
  # Non-significant cells: bold black text
  geom_text(data = heat_df %>% filter(!has_black),
            aes(label = label),
            size = txt_cell - 0.5, fontface = "bold", color = "black") +
  # Significant cells: bold white text (larger)
  geom_text(data = heat_df %>% filter(has_black),
            aes(label = label),
            size = (txt_cell - 0.5) * 1.15, fontface = "bold", color = "white") +
  scale_fill_gradient2(low = "#4393C3", mid = "white", high = "#D6604D",
                       midpoint = 0, limits = c(-0.8, 0.8),
                       oob = scales::squish,
                       name = "r / r-equiv.",
                       na.value = "white",
                       guide = guide_colorbar(
                         barwidth = unit(60, "mm"),
                         barheight = unit(4, "mm"),
                         title.position = "left",
                         title.vjust = 0.75
                       )) +
  scale_x_continuous(
    breaks = xpos,
    labels = col_labels[trait_order],
    limits = c(xmin_all, xmax_all),
    expand = c(0, 0)
  ) +
  scale_y_discrete(labels = NULL) +
  labs(x = NULL, y = NULL) +
  FIG_THEME +
  theme(axis.text.x        = element_text(angle = 45, hjust = 1,
                                           size = txt_cell * 1.8, face = "bold",
                                           lineheight = 0.9),
        axis.text.y        = element_blank(),
        axis.ticks.y       = element_blank(),
        axis.ticks.x       = element_blank(),
        panel.grid         = element_blank(),
        panel.border       = element_blank(),
        legend.position    = "bottom",
        legend.title       = element_text(size = txt_axis, face = "bold"),
        legend.text        = element_text(size = txt_cell),
        legend.margin      = margin(0, 0, 0, 0),
        legend.box.margin  = margin(t = -15, b = 0),
        legend.box.spacing = unit(0, "mm"),
        plot.margin        = margin(0, 2, 0, -3))

# =============================================================================
# 12. Preservation Zsummary bars (right flank)
# =============================================================================

pres_file <- file.path(DAT, "05_panel_E_preservation.csv")
has_preservation <- file.exists(pres_file)

if (has_preservation) {
  pres_raw <- read_csv(pres_file, show_col_types = FALSE) %>%
    mutate(module = paste0("ME", module)) %>%
    filter(module %in% non_grey)

  pres_df <- pres_raw %>%
    filter(module %in% mod_order) %>%
    mutate(module = factor(module, levels = mod_order),
           mod_color = gsub("^ME", "", as.character(module)))

  p_preservation <- ggplot(pres_df, aes(x = Zsummary, y = module)) +
    geom_vline(xintercept = c(2, 10), linetype = "dashed",
               color = "grey50", linewidth = 0.3) +
    geom_col(fill = pres_df$mod_color, color = "black",
             linewidth = 0.3, width = 0.65) +
    geom_text(aes(label = sprintf("%.1f", Zsummary), x = Zsummary / 2),
              size = txt_count, fontface = "bold",
              color = ifelse(pres_df$mod_color %in% light_modules,
                             "grey30", "white")) +
    scale_x_continuous(expand = c(0, 0),
                       limits = c(0, max(pres_df$Zsummary) * 1.05),
                       breaks = c(0, 2, 10, 20, 40)) +
    scale_y_discrete(labels = NULL) +
    labs(y = NULL,
         x = "Preservation (Zsummary)\n<2 none | 2\u201310 moderate | >10 strong\n(Langfelder et al. 2011)") +
    FIG_THEME +
    theme(axis.text.y        = element_blank(),
          axis.ticks.y       = element_blank(),
          axis.text.x        = element_text(size = txt_axis),
          axis.title.x       = element_text(size = txt_axis * 0.85, lineheight = 0.9,
                                            margin = margin(t = -2)),
          panel.grid.major.y = element_blank(),
          panel.grid.minor   = element_blank(),
          panel.border       = element_blank(),
          axis.line.x        = element_line(color = "black", linewidth = 0.3),
          legend.position    = "none",
          plot.margin        = margin(8, 2, -5, 0))
} else {
  warning("Preservation file not found: ", pres_file, " \u2014 skipping Zsummary column")
}

# =============================================================================
# 13. Patchwork layout assembly
# =============================================================================

# Proportional widths: gene counts (3) | heatmap (n_cols-scaled) | preservation (3)
heat_cols <- max(10, n_cols)

if (has_preservation) {
  design <- c(
    area(1, 4, 1, 3 + heat_cols),
    area(2, 1, 10, 3),
    area(2, 4, 10, 3 + heat_cols),
    area(2, 4 + heat_cols, 10, 6 + heat_cols)
  )
  p_composite <- wrap_elements(p_brackets) + p_counts + p_heat + p_preservation +
    plot_layout(design = design)
} else {
  design <- c(
    area(1, 4, 1, 3 + heat_cols),
    area(2, 1, 10, 3),
    area(2, 4, 10, 3 + heat_cols)
  )
  p_composite <- wrap_elements(p_brackets) + p_counts + p_heat +
    plot_layout(design = design)
}

p_composite <- p_composite +
  plot_annotation(
    title = "Age-Stratified Analysis \u2014 Training, Aging, & Phenotype Associations",
    subtitle = paste0(
      length(non_grey), " modules | ",
      "LMM: eigengene ~ time + (1|subject), BH per stratum | ",
      "Aging: unpaired t-test, BH per timepoint | ",
      "Phenotype: Pearson r, BH per column\n",
      "* p<.05  ** p<.01  *** p<.001  \u2020 p<.10 | ",
      "Black border = BH < 0.05"
    ),
    theme = theme(
      plot.title    = element_text(face = "bold", size = txt_axis * 2.5),
      plot.subtitle = element_text(size = txt_cell * 2, color = "grey30",
                                   face = "italic", lineheight = 0.9),
      plot.margin   = margin(2, 2, 2, 2)
    )
  )

# =============================================================================
# 14. Save outputs
# =============================================================================

out_base <- "SUPP_age_stratified_lmm"

ggsave(file.path(RPT_PNG, paste0(out_base, ".png")), p_composite,
       width = PA_W, height = PA_H, units = "mm",
       dpi = 300, limitsize = FALSE)
ggsave(file.path(RPT_PDF, paste0(out_base, ".pdf")), p_composite,
       width = PA_W, height = PA_H, units = "mm",
       device = pdf_device, limitsize = FALSE)

message(sprintf("\n  Saved: %s.{png,pdf}", out_base))
message(sprintf("  Dimensions: %d x %d mm (%d columns)", PA_W, PA_H, n_cols))
