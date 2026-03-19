#!/usr/bin/env Rscript
# YvO Normalization Reports — diagnostic plots + supplementary workbook
#
# Reads: c_data/00_report_intermediates.rds (from 01_run_normalization.R)
#
# Outputs:
#   b_reports/04_diagnostics.pdf    — 4-page custom QC report
#   c_data/05_normalization_supp.xlsx — supplementary workbook (3 sheets)

library(dplyr)
library(ggplot2)
library(ggrepel)
library(patchwork)
library(openxlsx)
library(cowplot)

setwd(rprojroot::find_rstudio_root_file())

# --- Configuration -----------------------------------------------------------

cfg <- list(
  report_dir = "01_normalization/b_reports",
  data_dir   = "01_normalization/c_data"
)

# --- Load intermediates ------------------------------------------------------

int <- readRDS(file.path(cfg$data_dir, "00_report_intermediates.rds"))

filter_bar_data  <- int$filter_bar_data
miss_bar_data    <- int$miss_bar_data
n_raw            <- int$n_raw
n_outliers       <- int$n_outliers
outlier_diag     <- int$outlier_diag
outlier_ids      <- int$outlier_ids
miss_thresh      <- int$miss_thresh
delta_thresh     <- int$delta_thresh
pca_pre          <- int$pca_pre
pca_post         <- int$pca_post
global_med       <- int$global_med
mad_val          <- int$mad_val
subj_var         <- int$subj_var
eta2_vals        <- int$eta2_vals
filter_log       <- int$filter_log
filtered_proteins <- int$filtered_proteins
dal_nrow         <- int$dal_nrow
dal_ncol         <- int$dal_ncol
cfg_full         <- int$cfg

# --- Palette & theme (self-contained) ----------------------------------------

pal_gt <- c(
  Young_Pre = scales::alpha("#4393C3", 0.5), Young_Post = "#4393C3",
  Old_Pre   = scales::alpha("#D6604D", 0.5), Old_Post   = "#D6604D"
)
shape_tp <- c(Pre = 16, Post = 17)
theme_qc <- theme_minimal(base_size = 12)

# --- PAGE 1: Filtering & missingness -----------------------------------------

p_filter <- ggplot(filter_bar_data, aes(step, n, fill = status)) +
  geom_col(width = 0.7) +
  geom_text(aes(label = n), position = position_stack(vjust = 0.5), size = 4) +
  scale_fill_manual(values = c(Retained = "#2166AC", Removed = "#B2182B")) +
  labs(x = NULL, y = "Proteins", fill = NULL, title = "Protein retention") +
  theme_qc + theme(axis.text.x = element_text(angle = 25, hjust = 1))

p_miss <- ggplot(miss_bar_data, aes(reorder(Col_ID, -n * (status == "Detected")),
                                     n, fill = status)) +
  geom_col(aes(alpha = is_outlier), width = 0.8) +
  scale_fill_manual(values = c(Detected = "#2166AC", Missing = "#D6604D")) +
  scale_alpha_manual(values = c("FALSE" = 1, "TRUE" = 0.4), guide = "none") +
  facet_grid(~ Group_Time, scales = "free_x", space = "free_x") +
  labs(x = NULL, y = "Proteins", fill = NULL,
       title = "Per-sample detection (all samples, outliers faded)") +
  theme_qc + theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 5),
                   strip.text = element_text(face = "bold"))

# --- PAGE 2: Outlier diagnostics ---------------------------------------------

# Shared aesthetics: blue = Young, red = Old, circle = Pre, triangle = Post
col_age <- c(Young = "#4393C3", Old = "#D6604D")
outlier_diag$age <- ifelse(grepl("^Y_", outlier_diag$prefix), "Young", "Old")

p_out_miss <- ggplot(outlier_diag, aes(pct_missing, delta_missing,
                                        color = age, shape = Timepoint)) +
  geom_point(size = 3) +
  geom_vline(xintercept = miss_thresh, linetype = "dashed", color = "red", alpha = 0.5) +
  geom_hline(yintercept = delta_thresh, linetype = "dashed", color = "red", alpha = 0.5) +
  geom_text_repel(data = \(d) filter(d, miss_flag),
                  aes(label = prefix), size = 2.5, show.legend = FALSE) +
  scale_color_manual(values = col_age) +
  scale_shape_manual(values = shape_tp) +
  labs(x = "Sample missingness (%)", y = "Delta missingness (|Post - Pre|)",
       title = "A: Missingness",
       subtitle = sprintf("IQR threshold: %.1f%% miss, %.1f%% delta | %d flagged",
                           miss_thresh, delta_thresh, sum(outlier_diag$miss_flag))) +
  coord_fixed(ratio = diff(range(outlier_diag$pct_missing)) /
                      max(1, diff(range(outlier_diag$delta_missing)))) +
  theme_qc

pca_outlier_df <- pca_pre$scores |>
  left_join(outlier_diag |> select(Col_ID, pca_flag, prefix, age), by = "Col_ID")

p_out_pca <- ggplot(pca_outlier_df, aes(PC1, PC2, color = age, shape = Timepoint)) +
  geom_point(size = 3.5, alpha = 0.85) +
  geom_text_repel(data = \(d) filter(d, pca_flag),
                  aes(label = prefix), size = 2.5, show.legend = FALSE) +
  scale_color_manual(values = col_age) +
  scale_shape_manual(values = shape_tp) +
  labs(x = sprintf("PC1 (%.1f%%)", pca_pre$var_exp[1]),
       y = sprintf("PC2 (%.1f%%)", pca_pre$var_exp[2]),
       title = "B: PCA Mahalanobis",
       subtitle = sprintf("chi-sq cutoff p < %.2f | %d flagged",
                           cfg_full$mahal_p, sum(pca_outlier_df$pca_flag))) +
  coord_fixed() +
  theme_qc

p_out_mad <- ggplot(outlier_diag, aes(reorder(prefix, sample_median),
                                       sample_median, color = age, shape = Timepoint)) +
  geom_point(size = 2.5) +
  geom_text_repel(data = \(d) filter(d, mad_flag),
                  aes(label = prefix), size = 2.5, show.legend = FALSE) +
  geom_hline(yintercept = global_med) +
  geom_hline(yintercept = global_med + c(-1, 1) * cfg_full$mad_k * mad_val,
             linetype = "dashed", color = "red", alpha = 0.5) +
  scale_color_manual(values = col_age) +
  scale_shape_manual(values = shape_tp) +
  labs(x = "Sample", y = "Median log2 intensity",
       title = "C: MAD median intensity",
       subtitle = sprintf("%dx MAD band | %d flagged",
                           cfg_full$mad_k, sum(outlier_diag$mad_flag))) +
  theme_qc + theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 4))

p_out_cor <- ggplot(outlier_diag, aes(reorder(prefix, median_cor),
                                       median_cor, color = age, shape = Timepoint)) +
  geom_point(size = 2.5) +
  geom_text_repel(data = \(d) filter(d, cor_flag),
                  aes(label = prefix), size = 2.5, show.legend = FALSE) +
  geom_hline(yintercept = median(outlier_diag$median_cor)) +
  geom_hline(yintercept = median(outlier_diag$median_cor) - cfg_full$mad_k * mad(outlier_diag$median_cor),
             linetype = "dashed", color = "red", alpha = 0.5) +
  scale_color_manual(values = col_age) +
  scale_shape_manual(values = shape_tp) +
  labs(x = "Sample", y = "Median pairwise correlation",
       title = "D: Inter-sample correlation",
       subtitle = sprintf("%dx MAD band | %d flagged",
                           cfg_full$mad_k, sum(outlier_diag$cor_flag))) +
  theme_qc + theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 4))

# --- PAGE 3: Post-normalization PCA ------------------------------------------

p_pca_post <- ggplot(pca_post$scores, aes(PC1, PC2,
                                            color = Group_Time, shape = Timepoint)) +
  geom_point(size = 3.5, alpha = 0.85) +
  stat_ellipse(aes(group = Group_Time), type = "norm", level = 0.68, linewidth = 0.7) +
  scale_color_manual(values = pal_gt) + scale_shape_manual(values = shape_tp) +
  labs(x = sprintf("PC1 (%.1f%%)", pca_post$var_exp[1]),
       y = sprintf("PC2 (%.1f%%)", pca_post$var_exp[2]),
       title = "Post-normalization PCA") +
  theme_qc + theme(legend.position = "bottom")

# --- PAGE 4: Variability -----------------------------------------------------

n_subjects <- length(unique(subj_var$Subject_ID))

p_cv <- ggplot(subj_var, aes(reorder(Subject_ID, iqr), iqr)) +
  geom_line(aes(group = Subject_ID), color = "gray60", linewidth = 0.4) +
  geom_point(aes(color = Group_Time, shape = Timepoint), size = 2.5) +
  scale_color_manual(values = pal_gt) +
  scale_shape_manual(values = shape_tp) +
  labs(x = "Subject", y = "IQR (log2 intensity)",
       title = "Per-subject variability",
       subtitle = sprintf("IQR of log2 intensities per sample | %d subjects | lines connect Pre-Post pairs",
                           n_subjects)) +
  theme_qc + theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 5),
                   legend.position = "bottom")

p_eta2 <- ggplot(data.frame(eta2 = eta2_vals[!is.na(eta2_vals)]), aes(eta2)) +
  geom_histogram(bins = 50, fill = "#2166AC", color = "white", alpha = 0.8) +
  geom_vline(xintercept = median(eta2_vals, na.rm = TRUE),
             linetype = "dashed", color = "red") +
  annotate("text", x = median(eta2_vals, na.rm = TRUE) + 0.02, y = Inf,
           vjust = 2, size = 3.5, color = "red",
           label = sprintf("median = %.2f", median(eta2_vals, na.rm = TRUE))) +
  labs(x = expression(eta^2 ~ "(between-group / total)"),
       y = "Proteins", title = "Variance partition by group",
       subtitle = sprintf("eta-sq = SS_between / SS_total | median = %.2f | higher = more group-driven",
                           median(eta2_vals, na.rm = TRUE))) +
  theme_qc

# --- ASSEMBLE PDF (4 pages) --------------------------------------------------

pdf(file.path(cfg$report_dir, "04_diagnostics.pdf"), width = 20, height = 10)

print(
  p_filter / p_miss +
    plot_layout(heights = c(1, 1.2)) +
    plot_annotation(
      title = "Protein Filtering & Detection",
      subtitle = sprintf("%d raw -> %d retained | %d samples",
                         n_raw, dal_nrow, dal_ncol),
      theme = theme(plot.title = element_text(size = 18, face = "bold"),
                    plot.subtitle = element_text(size = 14)))
)

# Extract shared legend from panel A, then suppress all individual legends
shared_legend <- cowplot::get_legend(
  p_out_miss + theme(legend.position = "bottom",
                     legend.justification = "center")
)

print(
  ((p_out_miss + theme(legend.position = "none")) |
   (p_out_pca  + theme(legend.position = "none"))) /
  ((p_out_mad  + theme(legend.position = "none")) |
   (p_out_cor  + theme(legend.position = "none"))) /
  wrap_elements(shared_legend) +
    plot_layout(heights = c(1, 1, 0.08)) +
    plot_annotation(
      title = "Outlier Diagnostics (4-method consensus)",
      subtitle = sprintf("Consensus rule: sample removed if \u22653 of 4 methods agree | %d removed",
                          n_outliers),
      theme = theme(plot.title = element_text(size = 18, face = "bold"),
                    plot.subtitle = element_text(size = 13)))
)

print(
  p_pca_post +
    plot_annotation(
      title = "Post-Normalization QC",
      theme = theme(plot.title = element_text(size = 18, face = "bold")))
)

print(
  (p_cv | p_eta2) +
    plot_annotation(
      title = "Variability Summary",
      theme = theme(plot.title = element_text(size = 18, face = "bold")))
)

dev.off()

# --- SUPPLEMENTARY WORKBOOK --------------------------------------------------

add_sheet <- function(wb, name, title, notes, df) {
  addWorksheet(wb, name)
  writeData(wb, name, x = title, startRow = 1)
  addStyle(wb, name, createStyle(fontSize = 13, textDecoration = "bold"), 1, 1)
  mergeCells(wb, name, cols = 1:ncol(df), rows = 1)
  writeData(wb, name, x = notes, startRow = 2)
  addStyle(wb, name, createStyle(fontSize = 10, fontColour = "#555555",
                                  wrapText = TRUE), 2, 1)
  mergeCells(wb, name, cols = 1:ncol(df), rows = 2)
  writeData(wb, name, x = df, startRow = 4,
            headerStyle = createStyle(textDecoration = "bold",
                                       border = "Bottom", fgFill = "#DCE6F1"))
  freezePane(wb, name, firstActiveRow = 5)
  setColWidths(wb, name, cols = 1:ncol(df), widths = "auto")
}

wb <- createWorkbook()

add_sheet(wb, "Pipeline_Summary",
  "Protein Filtering Pipeline",
  "step: filter stage | n_before/n_after: counts | pct_of_raw: cumulative retention",
  filter_log)

add_sheet(wb, "Outlier_Diagnostics",
  "Per-Sample Outlier Diagnostics",
  "miss_flag/pca_flag/mad_flag/cor_flag: per-method | consensus_outlier: TRUE if \u22653 of 4 agree",
  outlier_diag)

add_sheet(wb, "Filtered_Proteins",
  "Proteins Removed by Filtering",
  "removal_step: HPA or Missingness | identifiers: uniprot_id, gene, description",
  filtered_proteins)

saveWorkbook(wb, file.path(cfg$data_dir, "05_normalization_supp.xlsx"), overwrite = TRUE)

cat(sprintf("Done: 04_diagnostics.pdf + 05_normalization_supp.xlsx -> %s/, %s/\n",
            cfg$report_dir, cfg$data_dir))
