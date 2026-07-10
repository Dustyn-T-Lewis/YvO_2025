#!/usr/bin/env Rscript
# F02 Main — Proteome + DEP overview (6-panel 3×2 grid)
# A: PCA biplot  B: logFC density  C: DEPs per contrast
# D: UpSet overlap  E: fGSEA pathways  F: Barcode rank

setwd(here::here())

library(dplyr)
library(tidyr)
library(tibble)
library(stringr)
library(readr)
library(readxl)
library(ggplot2)
library(patchwork)
library(cowplot)
library(vegan)
library(ComplexHeatmap)
library(purrr)

source("04_Figures/shared/style.R")
source("04_Figures/shared/build_fgsea_cache.R")

# F02-specific overrides (from the old F02/style.R)
HEATMAP_LO <- "#2166AC"; HEATMAP_HI <- "#B2182B"
BASE_COUNT <- BASE_COUNT + 1.0
BASE_GENE  <- BASE_GENE  + 0.8
BASE_STAT  <- BASE_STAT  + 0.5

BASE    <- "04_Figures/F02"
RPT_PNG <- file.path(BASE, "b_reports", "main", "png")
RPT_PDF <- file.path(BASE, "b_reports", "main", "pdf")
PNL_PNG <- file.path(RPT_PNG, "panels")
PNL_PDF <- file.path(RPT_PDF, "panels")
DAT     <- file.path(BASE, "c_data")
for (d in c(PNL_PNG, PNL_PDF, DAT)) dir.create(d, recursive = TRUE, showWarnings = FALSE)

pdf_dev <- get_pdf_device()

dep_df <- read_csv("03_DEP/c_data/03_combined_results.csv",
                   show_col_types = FALSE)

dal_imp <- readRDS("02_imputation/c_data/01_DAList_imputed.rds")
imp_mat <- as.matrix(dal_imp$data)
imp_meta <- as_tibble(dal_imp$metadata) |>
  mutate(age    = factor(Group, levels = c("Young", "Old")),
         time   = factor(Timepoint, levels = c("Pre", "Post")),
         group  = factor(Group_Time,
                   levels = c("Young_Pre", "Young_Post", "Old_Pre", "Old_Post")),
         subject = sub("_(Pre|Post)$", "", Col_ID)) |>
  rename(sample_id = Col_ID)

DEP_XLSX <- "03_DEP/c_data/03_DEP_results.xlsx"
CONTRASTS <- c("Aging", "Training_Young", "Training_Old", "Interaction")

SET_LABELS <- c(Aging = "Aging", Training_Young = "Tr.(Y)",
                Training_Old = "Tr.(O)", Interaction = "Inter.")

PC_W <- 67; PC_H <- 55

pca <- prcomp(t(imp_mat), center = TRUE, scale. = TRUE)
var_pct <- round(100 * summary(pca)$importance[2, 1:2], 1)

set.seed(42)
boot_var <- replicate(1000, {
  idx <- sample(nrow(imp_mat), replace = TRUE)
  100 * summary(prcomp(t(imp_mat[idx, ]), center = TRUE, scale. = TRUE))$importance[2, 1:2]
})
var_ci <- data.frame(PC = c("PC1", "PC2"), var_pct = var_pct,
                     ci_lo = apply(boot_var, 1, quantile, 0.025),
                     ci_hi = apply(boot_var, 1, quantile, 0.975))

pca_df <- as.data.frame(pca$x[, 1:2]) |>
  mutate(sample_id = rownames(pca$x)) |>
  left_join(imp_meta, by = join_by(sample_id))

dist_mat <- dist(scale(t(imp_mat)))
set.seed(42)
perm_res <- adonis2(dist_mat ~ age * time, data = imp_meta,
                    permutations = how(nperm = 999, blocks = imp_meta$subject),
                    by = "terms")
perm_terms <- c("age", "time", "age:time")
perm_r2 <- perm_res[perm_terms, "R2"]
perm_pv <- perm_res[perm_terms, "Pr(>F)"]
perm_label <- sprintf(
  " PERMANOVA\nAge        R\u00b2 = %.3f,  %s\nTime       R\u00b2 = %.3f,  %s\nAge\u00d7Time R\u00b2 = %.3f,  %s",
  perm_r2[1], fmt_p(perm_pv[1]), perm_r2[2], fmt_p(perm_pv[2]),
  perm_r2[3], fmt_p(perm_pv[3]))

bd_age_p  <- permutest(betadisper(dist_mat, imp_meta$age),  permutations = 999)$tab$`Pr(>F)`[1]
bd_time_p <- permutest(betadisper(dist_mat, imp_meta$time), permutations = 999)$tab$`Pr(>F)`[1]
if (bd_age_p < 0.05 || bd_time_p < 0.05)
  warning("Heterogeneous dispersions — interpret PERMANOVA with caution")

pA <- ggplot(pca_df, aes(PC1, PC2, color = group, shape = group)) +
  stat_ellipse(aes(fill = group), geom = "polygon",
               alpha = 0.10, level = 0.80, show.legend = FALSE) +
  stat_ellipse(aes(group = group), level = 0.80, linewidth = 0.4,
               linetype = "dashed", show.legend = FALSE) +
  geom_point(size = 1.8, alpha = 0.85) +
  annotate("label", x = -Inf, y = Inf, label = perm_label,
           hjust = -0.02, vjust = 1.05, lineheight = 0.9,
           size = scale_text(BASE_COUNT, 44) - 0.5, color = "grey20", fontface = "bold",
           fill = alpha("white", 0.85), linewidth = 0.2,
           label.padding = unit(0.12, "lines")) +
  scale_color_manual(values = PCA_COLORS,
                     labels = c("Young Pre", "Young Post", "Old Pre", "Old Post"),
                     guide = guide_legend(override.aes = list(size = 1.6))) +
  scale_fill_manual(values = PCA_COLORS, guide = "none") +
  scale_shape_manual(values = PCA_SHAPES,
                     labels = c("Young Pre", "Young Post", "Old Pre", "Old Post")) +
  labs(title = "Sample PCA",
       subtitle = sprintf("n = %d, %s proteins (imputed)",
                          nrow(imp_meta), format(nrow(imp_mat), big.mark = ",")),
       x = sprintf("PC1 (%.1f%% [%.1f, %.1f])", var_pct[1], var_ci$ci_lo[1], var_ci$ci_hi[1]),
       y = sprintf("PC2 (%.1f%% [%.1f, %.1f])", var_pct[2], var_ci$ci_lo[2], var_ci$ci_hi[2]),
       tag = "a") +
  FIG_THEME + theme(plot.subtitle = element_text(size = FIG_SUBTITLE_SIZE,
                                                  face = "bold.italic", color = "grey30"),
                    legend.position = c(0.88, 0.15), legend.background = element_blank(),
                    legend.key = element_blank(), legend.title = element_blank(),
                    legend.text = element_text(size = FIG_LEGEND_TEXT + 0.5),
                    legend.key.size = unit(3, "mm"),
                    plot.margin = margin(6, 6, 2, 8))

write.csv(var_ci, file.path(DAT, "panel_A_pca_variance_ci.csv"), row.names = FALSE)
ggsave(file.path(PNL_PNG, "MAIN_panel_A_pca.png"), pA,
       width = PC_W, height = PC_H, units = "mm", dpi = 300)

pA_title <- "Sample PCA"
pA_subtitle <- sprintf("n = %d, %s proteins (imputed)",
                        nrow(imp_meta), format(nrow(imp_mat), big.mark = ","))
pA <- pA + labs(title = NULL, subtitle = NULL, tag = NULL)

PD_W <- 48; PD_H <- 55

lfc_long_all <- dep_df |>
  select(any_of(c("uniprot_id", "gene")), starts_with("logFC_")) |>
  pivot_longer(starts_with("logFC_"), names_to = "contrast", values_to = "logFC") |>
  mutate(contrast = str_remove(contrast, "logFC_")) |>
  filter(!is.na(logFC))
write_csv(lfc_long_all, file.path(DAT, "panel_B_logfc_long.csv"))

lfc_long <- lfc_long_all |>
  filter(contrast %in% c("Aging", "Training_Young", "Training_Old")) |>
  mutate(contrast = factor(contrast, levels = c("Aging", "Training_Young", "Training_Old")))

set.seed(42)
lfc_stats <- lfc_long |>
  summarise(med_abs_lfc = median(abs(logFC)),
            ci_lo = boot_median_ci(abs(logFC))[["lower"]],
            ci_hi = boot_median_ci(abs(logFC))[["upper"]],
            n_above_05 = sum(abs(logFC) > 0.5), .by = contrast)
write_csv(lfc_stats, file.path(DAT, "panel_B_stats.csv"))

lfc_stats$annotation <- sprintf("Med.|logFC| = %.2f\nn(>0.5) = %d",
                                lfc_stats$med_abs_lfc, lfc_stats$n_above_05)

dist_subtitle <- NULL
blunt <- tryCatch(as.data.frame(read_excel(DEP_XLSX, sheet = "blunting")), error = \(e) NULL)
if (!is.null(blunt)) {
  ks_row <- blunt[grepl("KS", blunt$test), ]
  lfc_ty <- abs(lfc_long$logFC[lfc_long$contrast == "Training_Young"])
  lfc_to <- abs(lfc_long$logFC[lfc_long$contrast == "Training_Old"])
  obs_ratio <- median(lfc_to) / median(lfc_ty)
  set.seed(42)
  ratio_ci <- quantile(replicate(2000, median(sample(lfc_to, replace = TRUE)) /
                                       median(sample(lfc_ty, replace = TRUE))),
                       c(0.025, 0.975))
  dist_subtitle <- sprintf("n = %s | KS D = %.2f, %s | Blunting = %.2f [%.2f, %.2f]",
                           format(length(unique(lfc_long$gene)), big.mark = ","),
                           ks_row$statistic, fmt_p(ks_row$p_value),
                           obs_ratio, ratio_ci[1], ratio_ci[2])
}

lfc_binwidth <- 4 / 50
pB <- ggplot(lfc_long, aes(logFC, fill = contrast)) +
  geom_histogram(bins = 50, color = "black", linewidth = 0.2, alpha = 0.85) +
  geom_density(aes(y = after_stat(count) * lfc_binwidth),
               alpha = 0.15, linewidth = 0.5, color = "grey20") +
  geom_text(data = lfc_stats,
            aes(x = 0, y = 695, label = CTR_SHORT[as.character(contrast)]),
            inherit.aes = FALSE, hjust = 0.5, vjust = 1,
            size = scale_text(BASE_COUNT, PD_W) + 0.2, color = "grey30", fontface = "bold") +
  geom_label(data = lfc_stats,
             aes(x = -1.08, y = 695, label = annotation),
             inherit.aes = FALSE, hjust = 0, vjust = 1,
             size = scale_text(BASE_COUNT, PD_W) - 0.5, color = "grey20", fontface = "bold",
             lineheight = 0.9, fill = alpha("white", 0.85), linewidth = 0.2,
             label.padding = unit(0.12, "lines")) +
  facet_wrap(~contrast, ncol = 1, labeller = labeller(contrast = CTR_SHORT)) +
  coord_cartesian(xlim = c(-1, 1), ylim = c(0, 700)) +
  scale_fill_manual(values = CONTRAST_COLORS[c("Aging", "Training_Young", "Training_Old")]) +
  labs(title = "Effect Size Distribution", subtitle = dist_subtitle,
       x = expression(bold(log[2]~FC)), y = " ", tag = "b") +
  FIG_THEME + theme(legend.position = "none", strip.text = element_blank(),
                    strip.background = element_blank(), panel.spacing.y = unit(0, "pt"),
                    plot.subtitle = element_text(size = FIG_SUBTITLE_SIZE,
                                                face = "bold.italic", color = "grey40"),
                    axis.text.y = element_text(size = FIG_AXIS_TEXT - 1.5, color = "grey40"),
                    axis.ticks.y = element_blank(),
                    plot.margin = margin(6, 4, 0, 4))
ggsave(file.path(PNL_PNG, "MAIN_panel_B_logfc_density.png"), pB,
       width = PD_W, height = PD_H, units = "mm", dpi = 300)

pB_title <- "Effect Size Distribution"
pB_subtitle <- dist_subtitle %||% ""
pB <- strip_for_composite(pB)

PA_W <- 67; PA_H <- 55
all_genes <- unique(dep_df$gene[!is.na(dep_df$gene)])

sig_sets <- list(); dir_map <- list()
for (ctr in CONTRASTS) {
  pi_vals <- dep_df[[paste0("pi_score_", ctr)]]
  lfc_vals <- dep_df[[paste0("logFC_", ctr)]]
  is_sig <- !is.na(pi_vals) & pi_vals < 0.05
  sig_sets[[ctr]] <- dep_df$gene[is_sig]
  dir_map[[ctr]] <- setNames(ifelse(lfc_vals[is_sig] > 0, "Up", "Down"), dep_df$gene[is_sig])
}
n_total <- length(all_genes)

pi_total <- sum(sapply(sig_sets, length))
fdr_total <- sum(sapply(CONTRASTS, \(ctr) {
  fdr_col <- paste0("adj.P.Val_", ctr)
  if (fdr_col %in% names(dep_df)) sum(dep_df[[fdr_col]] < 0.05, na.rm = TRUE) else 0
}))
p_total <- sum(sapply(CONTRASTS, \(ctr) {
  p_col <- paste0("P.Value_", ctr)
  if (p_col %in% names(dep_df)) sum(dep_df[[p_col]] < 0.05, na.rm = TRUE) else 0
}))

SET_DISPLAY_COLORS <- c("Aging" = unname(CONTRAST_COLORS["Aging"]),
                        "Tr.(Y)" = unname(CONTRAST_COLORS["Training_Young"]),
                        "Tr.(O)" = unname(CONTRAST_COLORS["Training_Old"]),
                        "Inter." = unname(CONTRAST_COLORS["Interaction"]))

frac_df <- bind_rows(lapply(CONTRASTS, \(ctr) {
  fdr_col <- paste0("adj.P.Val_", ctr); p_col <- paste0("P.Value_", ctr)
  tibble(contrast = SET_LABELS[ctr], threshold = c("p < 0.05", "q < 0.05", "Pi < 0.05"),
         n = c(sum(!is.na(dep_df[[p_col]]) & dep_df[[p_col]] < 0.05),
               sum(!is.na(dep_df[[fdr_col]]) & dep_df[[fdr_col]] < 0.05),
               length(sig_sets[[ctr]])))
})) |>
  mutate(contrast = factor(contrast, levels = rev(c("Aging", "Tr.(Y)", "Tr.(O)", "Inter."))),
         threshold = factor(threshold, levels = c("p < 0.05", "q < 0.05", "Pi < 0.05")),
         pct = 100 * n / n_total, fill_key = paste(contrast, threshold, sep = "___")) |>
  filter(n > 1)

FRAC_FILL <- c()
for (cname in names(SET_DISPLAY_COLORS)) {
  col <- unname(SET_DISPLAY_COLORS[cname])
  FRAC_FILL[paste(cname, "p < 0.05",  sep = "___")] <- adjustcolor(col, alpha.f = 0.15)
  FRAC_FILL[paste(cname, "q < 0.05",  sep = "___")] <- adjustcolor(col, alpha.f = 0.40)
  FRAC_FILL[paste(cname, "Pi < 0.05", sep = "___")] <- col
}

THRESH_LABEL <- c("p < 0.05" = "p", "q < 0.05" = "FDR", "Pi < 0.05" = "\u03A0")
label_df <- frac_df |>
  arrange(contrast, threshold) |>
  mutate(next_pct = lead(pct, default = 0), seg_width = pct - next_pct,
         label_y = (next_pct + pct) / 2,
         label = THRESH_LABEL[as.character(threshold)],
         .by = contrast) |>
  filter(seg_width > 0.3)

pC <- ggplot(frac_df, aes(contrast, pct, fill = fill_key)) +
  annotate("rect", xmin = 3.5, xmax = 4.5, ymin = -Inf, ymax = Inf,
           fill = CONTRAST_COLORS["Aging"], alpha = 0.20, color = "grey70", linewidth = 0.2) +
  annotate("rect", xmin = 2.5, xmax = 3.5, ymin = -Inf, ymax = Inf,
           fill = CONTRAST_COLORS["Training_Young"], alpha = 0.20, color = "grey70", linewidth = 0.2) +
  annotate("rect", xmin = 1.5, xmax = 2.5, ymin = -Inf, ymax = Inf,
           fill = CONTRAST_COLORS["Training_Old"], alpha = 0.20, color = "grey70", linewidth = 0.2) +
  annotate("rect", xmin = 0.5, xmax = 1.5, ymin = -Inf, ymax = Inf,
           fill = CONTRAST_COLORS["Interaction"], alpha = 0.20, color = "grey70", linewidth = 0.2) +
  geom_col(position = "identity", width = 0.75, color = "black", linewidth = 0.3) +
  shadowtext::geom_shadowtext(data = label_df,
             aes(x = contrast, y = label_y, label = label),
             inherit.aes = FALSE, hjust = 0.5, size = 2.2, fontface = "bold",
             color = "white", bg.color = "grey30", bg.r = 0.12) +
  scale_fill_manual(values = FRAC_FILL) +
  scale_y_continuous(expand = expansion(mult = c(0, 0)), breaks = seq(0, 28, by = 7), limits = c(0, 28)) +
  coord_flip() +
  labs(title = "DEPs per Contrast",
       subtitle = sprintf("%s proteins | \u03A0 %d | FDR %d | p %d",
                          format(n_total, big.mark = ","), pi_total, fdr_total, p_total),
       x = NULL, y = "% of proteome", tag = "c") +
  FIG_THEME + theme(plot.subtitle = element_text(size = FIG_SUBTITLE_SIZE,
                                                  face = "bold.italic", color = "grey40"),
                    legend.position = "none",
                    axis.text.y = element_text(face = "bold", size = FIG_AXIS_TEXT - 0.5))
ggsave(file.path(PNL_PNG, "MAIN_panel_C_dep_counts.png"), pC,
       width = PA_W, height = PA_H, units = "mm", dpi = 300)

pC_title <- "DEPs per Contrast"
pC_subtitle <- sprintf("%s proteins | \u03A0 %d | FDR %d | p %d",
                        format(n_total, big.mark = ","), pi_total, fdr_total, p_total)
pC <- strip_for_composite(pC)

source("04_Figures/F02/a_script/_panel_D_upset.R")

source("04_Figures/F02/a_script/_panel_E_fgsea.R")

source("04_Figures/F02/a_script/_panel_F_barcode.R")

layout <- "ABC\n###\nDEF"
ROW_TOP <- 0.458
SPACER  <- 0.00

# Per-panel margins (top breathing room + per-panel width nudges).
pA <- pA + theme(plot.margin = margin(12,   2, 12, 2))
pB <- pB + theme(plot.margin = margin(12, -52,  0, 5))
pC <- pC + theme(plot.margin = margin(12,   2,  0, 0))
pF <- pF + theme(plot.margin = margin(-3,   2, -12, 0))

composite <- wrap_elements(full = pA) + pB + pC +
             wrap_elements(full = pD) +
             wrap_elements(full = pE) +
             pF +
  plot_layout(
    design  = layout,
    widths  = c(160, 127, 138),
    heights = c(ROW_TOP, SPACER, 1 - ROW_TOP - SPACER)
  )

# Manual tag + title + subtitle placement via cowplot.
COMP_W <- 178
COMP_H <- 115
TAG_SZ <- composite_text_sizes(COMP_H)$tag
TTL_SZ <- composite_text_sizes(COMP_H)$title
SUB_SZ <- composite_text_sizes(COMP_H)$subtitle
TOP_Y  <- 0.995 - 2/COMP_H + 0.020 + 0.002 - 0.005 - 0.009
BOT_Y  <- 1 - ROW_TOP - SPACER + 0.005 - 2/COMP_H + 0.028 + 0.015 - 0.001 - 0.005 - 0.005 - 0.015 - 0.0045 - 0.003 - 0.009 + 0.017 - 0.007
X_LEFT  <- 0.002
X_MID   <- 0.372
X_RIGHT <- 0.630
X_MID_BOT   <- X_MID
X_RIGHT_BOT <- X_RIGHT + 0.010
X_TTL      <- 0.04
TTL_NUDGE  <- -0.008
BE_NUDGE   <- 0.021
TAG_DY     <- -0.002
SUB_OFFSET <- 0.022

composite <- ggdraw(composite) +
  # Panel A
  draw_label("A",         x = X_LEFT,                            y = TOP_Y - TAG_DY,     size = TAG_SZ, fontface = "bold",        hjust = 0, vjust = 1) +
  draw_label(pA_title,    x = X_LEFT + X_TTL,                    y = TOP_Y,              size = TTL_SZ, fontface = "bold",        hjust = 0, vjust = 1) +
  draw_label(pA_subtitle, x = X_LEFT + X_TTL,                    y = TOP_Y - SUB_OFFSET, size = SUB_SZ, fontface = "bold.italic", hjust = 0, vjust = 1, colour = "grey30") +
  # Panel B
  draw_label("B",         x = X_MID,                             y = TOP_Y - TAG_DY,     size = TAG_SZ, fontface = "bold",        hjust = 0, vjust = 1) +
  draw_label(pB_title,    x = X_MID + X_TTL - BE_NUDGE,          y = TOP_Y,              size = TTL_SZ, fontface = "bold",        hjust = 0, vjust = 1) +
  draw_label(pB_subtitle, x = X_MID + X_TTL - BE_NUDGE,          y = TOP_Y - SUB_OFFSET, size = SUB_SZ, fontface = "bold.italic", hjust = 0, vjust = 1, colour = "grey30") +
  # Panel C
  draw_label("C",         x = X_RIGHT + 0.042,                   y = TOP_Y - TAG_DY,     size = TAG_SZ, fontface = "bold",        hjust = 0, vjust = 1) +
  draw_label(pC_title,    x = X_RIGHT + 0.042 + X_TTL + TTL_NUDGE, y = TOP_Y,            size = TTL_SZ, fontface = "bold",        hjust = 0, vjust = 1) +
  draw_label(pC_subtitle, x = X_RIGHT + 0.042 + X_TTL + TTL_NUDGE, y = TOP_Y - SUB_OFFSET, size = SUB_SZ, fontface = "bold.italic", hjust = 0, vjust = 1, colour = "grey30") +
  # Panel D
  draw_label("D",         x = X_LEFT,                            y = BOT_Y - TAG_DY,     size = TAG_SZ, fontface = "bold",        hjust = 0, vjust = 1) +
  draw_label(pD_title,    x = X_LEFT + X_TTL,                    y = BOT_Y,              size = TTL_SZ, fontface = "bold",        hjust = 0, vjust = 1) +
  draw_label(pD_subtitle, x = X_LEFT + X_TTL,                    y = BOT_Y - SUB_OFFSET, size = SUB_SZ, fontface = "bold.italic", hjust = 0, vjust = 1, colour = "grey30") +
  # Panel E
  draw_label("E",         x = X_MID_BOT,                         y = BOT_Y - TAG_DY,     size = TAG_SZ, fontface = "bold",        hjust = 0, vjust = 1) +
  draw_label(pE_title,    x = X_MID + X_TTL - BE_NUDGE,          y = BOT_Y,              size = TTL_SZ, fontface = "bold",        hjust = 0, vjust = 1) +
  draw_label(pE_subtitle, x = X_MID + X_TTL - BE_NUDGE,          y = BOT_Y - SUB_OFFSET, size = SUB_SZ, fontface = "bold.italic", hjust = 0, vjust = 1, colour = "grey30") +
  # Panel F
  draw_label("F",         x = X_RIGHT_BOT + 0.042,               y = BOT_Y - TAG_DY,     size = TAG_SZ, fontface = "bold",        hjust = 0, vjust = 1) +
  draw_label(pF_title,    x = X_RIGHT + 0.042 + X_TTL + TTL_NUDGE, y = BOT_Y,            size = TTL_SZ, fontface = "bold",        hjust = 0, vjust = 1) +
  draw_label(pF_subtitle, x = X_RIGHT + 0.042 + X_TTL + TTL_NUDGE, y = BOT_Y - SUB_OFFSET, size = SUB_SZ, fontface = "bold.italic", hjust = 0, vjust = 1, colour = "grey30")

# Reset to parent dir (helpers overwrite RPT_PDF/RPT_PNG to panels/)
RPT_PDF <- file.path(BASE, "b_reports", "main", "pdf")
RPT_PNG <- file.path(BASE, "b_reports", "main", "png")

ggsave(file.path(RPT_PDF, "MAIN_F02_composite.pdf"), composite,
       width = COMP_W, height = COMP_H, units = "mm", device = pdf_dev)
ggsave(file.path(RPT_PNG, "MAIN_F02_composite.png"), composite,
       width = COMP_W, height = COMP_H, units = "mm", dpi = 300)

message("F02 main composite done")
