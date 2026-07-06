#!/usr/bin/env Rscript
# F00 — Pipeline QC Supplementary Figures
#
# Two supp composites (no main figure — F00 is QC-only):
#   SUPP_F00_normalization: A–G (filter → normalize → outlier → missingness)
#   SUPP_F00_imputation:    H–N (classify → benchmark → impute → DEP)
#
# Data: reads upstream 00_report_intermediates.rds + benchmark CSV + DEP xlsx
# Outputs: 2 SUPP composites (PDF+PNG), per-panel PNGs, F00_supplementary.xlsx

setwd(rprojroot::find_root(rprojroot::has_file("setup.R")))

library(dplyr)
library(tibble)
library(tidyr)
library(readr)
library(readxl)
library(ggplot2)
library(scales)
library(patchwork)
library(openxlsx)

source("04_Figures_v2/shared/style.R")
source("04_Figures_v2/shared/figure_supplement_helpers.R")

BASE <- "04_Figures_v2/F00"
RPT_PNG <- file.path(BASE, "b_reports", "supp", "png")
RPT_PDF <- file.path(BASE, "b_reports", "supp", "pdf")
PNL_PNG <- file.path(RPT_PNG, "panels")
DAT <- file.path(BASE, "c_data")
for (d in c(PNL_PNG, RPT_PDF, DAT)) dir.create(d, recursive = TRUE, showWarnings = FALSE)

PW <- 89
PH <- 65

int_norm <- readRDS("01_normalization/c_data/00_report_intermediates.rds")
int_imp <- readRDS("02_imputation/c_data/00_report_intermediates.rds")
bench <- read_csv("02_imputation/c_data/benchmark/04_composite_ranking.csv", show_col_types = FALSE)
DEP_XLSX <- "03_DEP/c_data/03_DEP_results.xlsx"
da_summ <- as.data.frame(read_excel(DEP_XLSX, sheet = "DA_summary"))

fcasc <- int_norm$filter_log |>
  mutate(
    step = factor(step, levels = step),
    removed = ifelse(is.na(n_removed), 0L, n_removed)
  )

pA <- ggplot(fcasc, aes(step, n_after)) +
  geom_col(fill = "#2166AC", width = 0.7) +
  geom_text(aes(label = comma(n_after)), vjust = -1.8, size = 2.2, fontface = "bold") +
  geom_text(aes(label = ifelse(removed > 0, sprintf("-%s", comma(removed)), "")),
    vjust = -0.4, size = 2.0, color = "#B2182B"
  ) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.18)), labels = comma) +
  labs(
    x = NULL, y = "Proteins retained", tag = "A",
    title = "Protein filter cascade",
    subtitle = sprintf(
      "%s \u2192 %s (%.1f%%)", comma(fcasc$n_after[1]),
      comma(tail(fcasc$n_after, 1)), tail(fcasc$pct_of_raw, 1)
    )
  ) +
  FIG_THEME +
  theme(axis.text.x = element_text(angle = 22, hjust = 1, size = 5.5))

ggsave(file.path(PNL_PNG, "SUPP_panel_A_filter_cascade.png"), pA,
  width = PW, height = PH, units = "mm", dpi = 300
)

prot_miss <- rowSums(is.na(int_norm$data_pre_outlier)) /
  ncol(int_norm$data_pre_outlier) * 100
miss_hist_df <- tibble(pct_miss = prot_miss)
q99 <- quantile(prot_miss, 0.99, na.rm = TRUE)

pB <- ggplot(miss_hist_df, aes(pct_miss)) +
  geom_histogram(binwidth = 2.5, fill = "#4393C3", color = "white", alpha = 0.85) +
  geom_vline(xintercept = q99, linetype = "dashed", color = "#B2182B", linewidth = 0.4) +
  annotate("text",
    x = q99, y = Inf, label = sprintf("99th = %.1f%%", q99),
    vjust = 1.6, hjust = 1.05, color = "#B2182B", size = 2.2
  ) +
  scale_y_continuous(labels = comma) +
  labs(
    x = "% missing", y = "Proteins", tag = "B",
    title = "Per-protein missingness distribution",
    subtitle = sprintf(
      "%s proteins | median = %.1f%%",
      comma(nrow(miss_hist_df)), median(prot_miss)
    )
  ) +
  FIG_THEME

ggsave(file.path(PNL_PNG, "SUPP_panel_B_protein_missingness.png"), pB,
  width = PW, height = PH, units = "mm", dpi = 300
)

pca_pre_df <- int_norm$pca_pre$scores |>
  left_join(int_norm$outlier_diag |> select(Col_ID, prefix),
    by = join_by(Col_ID)
  ) |>
  mutate(age = factor(ifelse(grepl("^O", prefix), "Old", "Young"),
    levels = c("Young", "Old")
  ))

pC <- ggplot(pca_pre_df, aes(PC1, PC2, color = age, fill = age)) +
  stat_ellipse(aes(group = age),
    geom = "polygon",
    alpha = 0.12, level = 0.80, linewidth = 0.3, show.legend = FALSE
  ) +
  geom_point(aes(shape = Timepoint), size = 1.6, alpha = 0.85) +
  scale_color_manual(values = AGE_COLORS, name = NULL) +
  scale_fill_manual(values = AGE_COLORS, guide = "none") +
  scale_shape_manual(values = SHAPE_TP, name = NULL) +
  labs(
    x = sprintf("PC1 (%.1f%%)", int_norm$pca_pre$var_exp[1]),
    y = sprintf("PC2 (%.1f%%)", int_norm$pca_pre$var_exp[2]),
    tag = "C", title = "PCA before normalization",
    caption = "Normalization QC; 80% ellipses, not a test of group separation (small n)"
  ) +
  FIG_THEME +
  theme_sub_legend(position = "top", key.size = unit(2.5, "mm"))

ggsave(file.path(PNL_PNG, "SUPP_panel_C_pca_pre.png"), pC,
  width = PW, height = PH, units = "mm", dpi = 300
)

pD <- ggplot(
  int_norm$pca_post$scores,
  aes(PC1, PC2, color = Group_Time, fill = Group_Time)
) +
  stat_ellipse(aes(group = Group_Time),
    geom = "polygon",
    alpha = 0.12, level = 0.80, linewidth = 0.3, show.legend = FALSE
  ) +
  geom_point(aes(shape = Timepoint), size = 1.6, alpha = 0.85) +
  scale_color_manual(
    values = PCA_COLORS, name = NULL,
    labels = c("Young Pre", "Young Post", "Old Pre", "Old Post")
  ) +
  scale_fill_manual(values = PCA_COLORS, guide = "none") +
  scale_shape_manual(values = SHAPE_TP, guide = "none") +
  labs(
    x = sprintf("PC1 (%.1f%%)", int_norm$pca_post$var_exp[1]),
    y = sprintf("PC2 (%.1f%%)", int_norm$pca_post$var_exp[2]),
    tag = "D", title = "PCA after cyclic loess",
    caption = "Normalization QC; 80% ellipses, not a test of group separation (small n)"
  ) +
  FIG_THEME +
  theme_sub_legend(
    position = "top", key.size = unit(2.5, "mm"),
    text = element_text(size = 5.5)
  )

ggsave(file.path(PNL_PNG, "SUPP_panel_D_pca_post.png"), pD,
  width = PW, height = PH, units = "mm", dpi = 300
)

eta_df <- tibble(eta2 = as.numeric(int_norm$eta2_vals))
eta_med <- median(eta_df$eta2, na.rm = TRUE)
eta_90 <- quantile(eta_df$eta2, 0.90, na.rm = TRUE)

pE <- ggplot(eta_df, aes(eta2)) +
  geom_density(fill = "#4CAF50", alpha = 0.55, linewidth = 0.4) +
  geom_vline(xintercept = eta_med, color = "grey30", linewidth = 0.4) +
  geom_vline(
    xintercept = eta_90, linetype = "dashed",
    color = "#E05A4E", linewidth = 0.4
  ) +
  annotate("text",
    x = eta_med, y = Inf, label = sprintf("median = %.2f", eta_med),
    vjust = 1.5, hjust = -0.1, size = 2.2
  ) +
  annotate("text",
    x = eta_90, y = Inf, label = sprintf("90th = %.2f", eta_90),
    vjust = 3.0, hjust = -0.1, size = 2.2, color = "#E05A4E"
  ) +
  labs(
    x = expression(eta^2), y = "Density", tag = "E",
    title = expression("Biological signal retention (" * eta^2 * ")"),
    subtitle = sprintf("Group effect | %s proteins", comma(nrow(eta_df)))
  ) +
  FIG_THEME

ggsave(file.path(PNL_PNG, "SUPP_panel_E_eta_squared.png"), pE,
  width = PW, height = PH, units = "mm", dpi = 300
)

od <- int_norm$outlier_diag |>
  mutate(
    age = factor(ifelse(grepl("^O", prefix), "Old", "Young"),
      levels = c("Young", "Old")
    ),
    status = case_when(
      consensus_outlier ~ "Outlier",
      n_flags > 0 ~ "Flagged",
      TRUE ~ "Clean"
    )
  )

pF <- ggplot(od, aes(mahal_dist, n_flags, color = status, shape = Timepoint)) +
  geom_jitter(width = 0.04, height = 0.12, size = 1.8, alpha = 0.85) +
  scale_color_manual(values = c(
    Clean = "grey55", Flagged = "#E6A100",
    Outlier = "#B2182B"
  ), name = NULL) +
  scale_shape_manual(values = SHAPE_TP, guide = "none") +
  scale_y_continuous(breaks = 0:4) +
  labs(
    x = "Mahalanobis distance", y = "QC flags", tag = "F",
    title = "Outlier detection consensus",
    subtitle = sprintf(
      "%d/%d removed: %s",
      sum(od$consensus_outlier), nrow(od),
      paste(int_norm$outlier_ids, collapse = ", ")
    )
  ) +
  FIG_THEME +
  theme_sub_legend(position = "top", key.size = unit(2.5, "mm"))

ggsave(file.path(PNL_PNG, "SUPP_panel_F_outlier_consensus.png"), pF,
  width = PW, height = PH, units = "mm", dpi = 300
)

miss_bar <- int_norm$miss_bar_data |>
  filter(status == "Missing") |>
  mutate(Group_Time = factor(Group_Time,
    levels = c("Young_Pre", "Young_Post", "Old_Pre", "Old_Post")
  ))

pG <- ggplot(miss_bar, aes(reorder(Col_ID, -n), n, fill = Group_Time)) +
  geom_col(aes(alpha = is_outlier), width = 0.8) +
  scale_fill_manual(
    values = GROUP_FILL, name = NULL,
    labels = c("Young Pre", "Young Post", "Old Pre", "Old Post")
  ) +
  scale_alpha_manual(values = c("FALSE" = 1, "TRUE" = 0.3), guide = "none") +
  labs(
    x = NULL, y = "Missing proteins", tag = "G",
    title = "Per-sample missing protein counts",
    subtitle = sprintf(
      "%d samples | %d outliers (faded) | %.1f%% overall",
      length(unique(miss_bar$Col_ID)),
      int_norm$n_outliers,
      int_imp$pct_miss
    )
  ) +
  FIG_THEME +
  theme(
    axis.text.x = element_blank(), axis.ticks.x = element_blank(),
    legend.position = "top", legend.key.size = unit(2.5, "mm"),
    legend.text = element_text(size = 5.5)
  )

ggsave(file.path(PNL_PNG, "SUPP_panel_G_sample_missingness.png"), pG,
  width = PW * 2, height = PH * 0.75, units = "mm", dpi = 300
)

mc <- int_imp$miss_class |>
  mutate(classification = factor(classification, levels = c("Complete", "MAR", "MNAR")))
n_total <- nrow(mc)
class_lab <- mc |>
  count(classification) |>
  mutate(lab = sprintf("%s: %d (%.0f%%)", classification, n, 100 * n / n_total))

pH <- ggplot(mc, aes(mean_intensity, pct_miss, color = classification)) +
  geom_point(alpha = 0.45, size = 0.7) +
  scale_color_manual(values = PAL_CLASS, name = NULL, labels = class_lab$lab) +
  scale_y_continuous(labels = \(x) paste0(x, "%")) +
  labs(
    x = "Mean log2 intensity", y = "% missing", tag = "H",
    title = "MAR vs. MNAR classification",
    subtitle = sprintf(
      "%s | %s proteins", int_imp$classification_method,
      comma(n_total)
    )
  ) +
  FIG_THEME +
  theme_sub_legend(
    position = "top", key.size = unit(2.5, "mm"),
    text = element_text(size = 5.5)
  )

ggsave(file.path(PNL_PNG, "SUPP_panel_H_miss_class_scatter.png"), pH,
  width = PW, height = PH, units = "mm", dpi = 300
)

class_counts <- mc |>
  count(classification) |>
  mutate(pct = 100 * n / sum(n))

pI <- ggplot(class_counts, aes(classification, n, fill = classification)) +
  geom_col(width = 0.6) +
  geom_text(aes(label = sprintf("%d\n(%.0f%%)", n, pct)),
    vjust = -0.2, size = 2.2, fontface = "bold"
  ) +
  scale_fill_manual(values = PAL_CLASS, name = NULL) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.22))) +
  labs(
    x = NULL, y = "Proteins", tag = "I",
    title = "Missingness class counts",
    subtitle = sprintf(
      "MAR: %s vals | MNAR: %s vals | %.1f%% total",
      comma(int_imp$mar_miss_vals), comma(int_imp$mnar_miss_vals),
      int_imp$pct_miss
    )
  ) +
  FIG_THEME

ggsave(file.path(PNL_PNG, "SUPP_panel_I_miss_class_bar.png"), pI,
  width = PW, height = PH, units = "mm", dpi = 300
)

bench_plot <- bench |>
  filter(method != "Non_imputed") |>
  arrange(rank) |>
  mutate(
    method = factor(method, levels = rev(method)),
    bar_col = case_when(
      method == "missForest" ~ "Selected",
      rank <= 5 ~ "Top 5",
      TRUE ~ "Other"
    )
  )

mf_rank <- bench$rank[bench$method == "missForest"]
mf_composite <- bench$composite[bench$method == "missForest"]

pJ <- ggplot(bench_plot, aes(composite, method, fill = bar_col)) +
  geom_col(width = 0.7) +
  geom_text(aes(label = sprintf("%.3f", composite)),
    hjust = -0.08, size = 1.9, fontface = "bold"
  ) +
  scale_fill_manual(values = c(
    Selected = "#E41A1C", `Top 5` = "#377EB8",
    Other = "grey70"
  ), name = NULL) +
  scale_x_continuous(
    expand = expansion(mult = c(0, 0.18)),
    limits = c(0, NA)
  ) +
  labs(
    x = "Composite score", y = NULL, tag = "J",
    title = sprintf("Imputation method benchmark (%d)", nrow(bench_plot)),
    subtitle = sprintf(
      "missForest selected (rank #%d, composite = %.3f)",
      mf_rank, mf_composite
    )
  ) +
  FIG_THEME +
  theme(
    legend.position = "top", legend.key.size = unit(2.5, "mm"),
    axis.text.y = element_text(size = 5)
  )

ggsave(file.path(PNL_PNG, "SUPP_panel_J_benchmark.png"), pJ,
  width = PW, height = PH * 1.3, units = "mm", dpi = 300
)

obs_vals <- as.numeric(int_imp$mat[!int_imp$was_na])
imp_vals <- as.numeric(int_imp$mat_imp[int_imp$was_na])
dens_df <- bind_rows(
  tibble(value = obs_vals, type = "Observed"),
  tibble(value = imp_vals, type = "Imputed")
) |>
  mutate(type = factor(type, levels = c("Observed", "Imputed")))

pK <- ggplot(dens_df, aes(value, fill = type, color = type)) +
  geom_density(alpha = 0.4, linewidth = 0.4) +
  scale_fill_manual(values = c(Observed = "#377EB8", Imputed = "#E41A1C"), name = NULL) +
  scale_color_manual(values = c(Observed = "#377EB8", Imputed = "#E41A1C"), name = NULL) +
  annotate("text",
    x = Inf, y = Inf,
    label = sprintf("OOB = %.3f", int_imp$oob_error),
    hjust = 1.1, vjust = 1.5, size = 2.2, fontface = "italic", color = "grey30"
  ) +
  labs(
    x = "log2 intensity", y = "Density", tag = "K",
    title = "Observed vs. imputed intensity",
    subtitle = sprintf(
      "%s observed | %s imputed",
      comma(length(obs_vals)), comma(length(imp_vals))
    )
  ) +
  FIG_THEME +
  theme_sub_legend(position = "top", key.size = unit(2.5, "mm"))

ggsave(file.path(PNL_PNG, "SUPP_panel_K_imputation_density.png"), pK,
  width = PW, height = PH, units = "mm", dpi = 300
)

audit <- int_imp$mnar_audit
shift_mean <- mean(audit$shift, na.rm = TRUE)

pL <- ggplot(audit, aes(shift)) +
  geom_histogram(
    binwidth = 0.025, fill = "#E41A1C", alpha = 0.70,
    color = "white", linewidth = 0.2
  ) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey30", linewidth = 0.4) +
  geom_vline(xintercept = shift_mean, color = "#B2182B", linewidth = 0.4) +
  annotate("text",
    x = shift_mean, y = Inf,
    label = sprintf("mean = %+.3f", shift_mean),
    vjust = 1.4, hjust = -0.1, size = 2.2, color = "#B2182B"
  ) +
  labs(
    x = "Shift (log2)", y = "MNAR proteins", tag = "L",
    title = "MNAR imputed value shift",
    subtitle = sprintf("%d MNAR proteins | missForest", nrow(audit))
  ) +
  FIG_THEME

ggsave(file.path(PNL_PNG, "SUPP_panel_L_mnar_shift.png"), pL,
  width = PW, height = PH, units = "mm", dpi = 300
)

samp_integrity <- tibble(
  Col_ID = colnames(int_imp$mat),
  pre = colMeans(int_imp$mat, na.rm = TRUE),
  post = colMeans(int_imp$mat_imp)
) |>
  left_join(int_imp$meta |> select(Col_ID, Group_Time) |> distinct(),
    by = join_by(Col_ID)
  )
r2 <- cor(samp_integrity$pre, samp_integrity$post)^2

pM <- ggplot(samp_integrity, aes(pre, post, color = Group_Time)) +
  geom_abline(
    slope = 1, intercept = 0, linetype = "dashed",
    color = "grey60", linewidth = 0.4
  ) +
  geom_point(size = 1.8, alpha = 0.85) +
  scale_color_manual(
    values = GROUP_FILL, name = NULL,
    labels = c("Young Pre", "Young Post", "Old Pre", "Old Post")
  ) +
  labs(
    x = "Mean log2 (observed)", y = "Mean log2 (imputed)", tag = "M",
    title = "Pre- vs. post-imputation means",
    subtitle = sprintf("R\u00b2 = %.4f | OOB = %.3f", r2, int_imp$oob_error)
  ) +
  FIG_THEME +
  theme_sub_legend(
    position = "top", key.size = unit(2.5, "mm"),
    text = element_text(size = 5.5)
  )

ggsave(file.path(PNL_PNG, "SUPP_panel_M_sample_integrity.png"), pM,
  width = PW, height = PH, units = "mm", dpi = 300
)

dep_counts <- da_summ |>
  filter(type %in% c("up", "down")) |>
  summarise(
    `p<0.05` = sum(sig.PVal),
    `FDR<0.10` = sum(sig.FDR),
    `FDR<0.05` = sum(sig.FDR.05),
    `Pi<0.05` = sum(sig.Pi),
    .by = contrast
  ) |>
  pivot_longer(-contrast, names_to = "threshold", values_to = "n") |>
  mutate(
    contrast = factor(contrast,
      levels = intersect(CONTRAST_ORDER, unique(contrast))
    ),
    threshold = factor(threshold,
      levels = c("p<0.05", "FDR<0.10", "FDR<0.05", "Pi<0.05")
    )
  )

pN <- ggplot(dep_counts, aes(threshold, contrast, fill = n)) +
  geom_tile(color = "white", linewidth = 0.8) +
  geom_text(aes(label = comma(n), color = n > max(n) * 0.5),
    size = 2.5, fontface = "bold"
  ) +
  scale_color_manual(values = c("TRUE" = "white", "FALSE" = "grey20"), guide = "none") +
  scale_fill_gradient(low = "#DEEBF7", high = "#08519C", name = "DEPs") +
  labs(
    x = NULL, y = NULL, tag = "N",
    title = "DEP counts by contrast and threshold",
    subtitle = sprintf(
      "limma + duplicateCorrelation | %s proteins",
      comma(int_norm$dal_nrow)
    )
  ) +
  FIG_THEME +
  theme(
    axis.text.x = element_text(angle = 15, hjust = 1, size = 6),
    legend.position = "right", legend.key.size = unit(3, "mm")
  )

ggsave(file.path(PNL_PNG, "SUPP_panel_N_dep_heatmap.png"), pN,
  width = PW * 2, height = PH * 0.75, units = "mm", dpi = 300
)

COMP_W <- 178
COMP_H <- 245
txt <- composite_text_sizes(COMP_H)

# Composite strip: keep titles + tags + compact legend, drop subtitles
# Title size 7pt matches FIG_TITLE_SIZE (journal T1 tier)
COMP_TITLE_SZ <- 7

sl <- function(p) {
  p + labs(subtitle = NULL) +
    theme(
      plot.title = element_text(
        face = "bold", size = COMP_TITLE_SZ,
        margin = margin(b = 1)
      ),
      legend.position = "bottom",
      legend.key.size = unit(2, "mm"),
      legend.text = element_text(size = 4.5),
      legend.title = element_text(size = 5),
      legend.margin = margin(t = -1, b = 0),
      legend.spacing.y = unit(0, "mm")
    )
}
# Variant with right-side legend (for full-width gradient panels)
sr <- function(p) {
  p + labs(subtitle = NULL) +
    theme(
      plot.title = element_text(
        face = "bold", size = COMP_TITLE_SZ,
        margin = margin(b = 1)
      ),
      legend.position = "right",
      legend.key.size = unit(3, "mm"),
      legend.key.width = unit(2.5, "mm"),
      legend.text = element_text(size = 4.5),
      legend.title = element_text(size = 5),
      legend.margin = margin(l = -1)
    )
}

# PCA panels: drop coord_fixed() for uniform alignment in composite
pC_comp <- pC + coord_cartesian()
pD_comp <- pD + coord_cartesian()

page1 <- (sl(pA) | sl(pB)) /
  (sl(pC_comp) | sl(pD_comp)) /
  (sl(pE) | sl(pF)) /
  sl(pG) +
  plot_layout(heights = c(1, 1, 1, 0.75)) +
  plot_annotation(
    title = "Pipeline QC \u2014 Pre-Processing",
    subtitle = sprintf(
      "DIA-NN \u2192 HPA filter \u2192 cyclic loess \u2192 outlier removal | %s \u2192 %s proteins \u00d7 %d samples",
      comma(int_norm$n_raw), comma(int_norm$dal_nrow), int_norm$dal_ncol
    ),
    theme = theme(
      plot.title = element_text(face = "bold", size = txt$title),
      plot.subtitle = element_text(
        face = "italic", size = txt$subtitle - 1,
        color = "grey30"
      )
    )
  ) &
  theme(plot.tag = element_text(face = "bold", size = txt$tag))

page2 <- (sl(pH) | sl(pI)) /
  (sl(pJ) | sl(pK)) /
  (sl(pL) | sl(pM)) /
  sr(pN) +
  plot_layout(heights = c(1, 1.2, 1, 0.75)) +
  plot_annotation(
    title = "Pipeline QC \u2014 Imputation & Differential Expression",
    subtitle = sprintf(
      "missForest (OOB = %.3f) | %d MAR + %d MNAR proteins | limma + duplicateCorrelation",
      int_imp$oob_error, int_imp$n_mar_prots, int_imp$n_mnar_prots
    ),
    theme = theme(
      plot.title = element_text(face = "bold", size = txt$title),
      plot.subtitle = element_text(
        face = "italic", size = txt$subtitle - 1,
        color = "grey30"
      )
    )
  ) &
  theme(plot.tag = element_text(face = "bold", size = txt$tag))

# Separate SUPP composites (PDF + PNG each)
pdf_dev <- get_pdf_device()

pdf_dev(file.path(RPT_PDF, "SUPP_F00_normalization.pdf"),
  width = COMP_W / 25.4, height = COMP_H / 25.4
)
print(page1)
dev.off()

pdf_dev(file.path(RPT_PDF, "SUPP_F00_imputation.pdf"),
  width = COMP_W / 25.4, height = COMP_H / 25.4
)
print(page2)
dev.off()

ggsave(file.path(RPT_PNG, "SUPP_F00_normalization.png"), page1,
  width = COMP_W, height = COMP_H, units = "mm", dpi = 300
)
ggsave(file.path(RPT_PNG, "SUPP_F00_imputation.png"), page2,
  width = COMP_W, height = COMP_H, units = "mm", dpi = 300
)

dens_summary <- dens_df |>
  summarise(
    n = n(), mean = mean(value), median = median(value),
    sd = sd(value), .by = type
  ) |>
  mutate(oob_error = ifelse(type == "Imputed", int_imp$oob_error, NA_real_))

build_workbook(
  file.path(DAT, "F00_supplementary.xlsx"),
  title = "F00 \u2014 Pipeline QC source data",
  description = "Panels A\u2013N source data for two supplementary QC composites: normalization (A\u2013G) and imputation (H\u2013N).",
  sheet_specs = list(
    list(name = "panel_A", df = as.data.frame(fcasc)),
    list(name = "panel_B", df = as.data.frame(miss_hist_df)),
    list(name = "panel_C", df = as.data.frame(pca_pre_df)),
    list(name = "panel_D", df = as.data.frame(int_norm$pca_post$scores)),
    list(name = "panel_E", df = as.data.frame(eta_df)),
    list(name = "panel_F", df = as.data.frame(od)),
    list(name = "panel_G", df = as.data.frame(miss_bar)),
    list(name = "panel_H", df = as.data.frame(mc)),
    list(name = "panel_I", df = as.data.frame(class_counts)),
    list(name = "panel_J", df = as.data.frame(bench)),
    list(name = "panel_K", df = as.data.frame(dens_summary)),
    list(name = "panel_L", df = as.data.frame(audit)),
    list(name = "panel_M", df = as.data.frame(samp_integrity)),
    list(name = "panel_N", df = as.data.frame(dep_counts))
  )
)

message("F00 supp composites done")
