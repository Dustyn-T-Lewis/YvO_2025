################################################################################
#   YvO Figure 1 — Proteomic Landscape of Aging & Resistance Training
#   Panels: A (CV%), B (logFC density), C (PCA), D (fraction altered),
#           E (UpSet), F (fgsea bar)
#   Supplemental: S1.7 (pi-score distributions)
#
#   References:
#     - Korotkevich et al. 2021, bioRxiv 060012 — fgsea
#     - Zyla et al. 2017, BMC Bioinform 18:304 — t-stat for GSEA ranking
#     - Sayols 2023, microPublication Biology PMC10155054 — rrvgo
#     - Liberzon et al. 2015, Cell Syst 1(6):417-425 — MSigDB Hallmark
#     - Gu & Hubschmann 2022, Bioinformatics 38(5):1474-1476 — ComplexHeatmap
#     - Brenes 2024, J Proteome Res PMC11629372 — CV on linear scale
#     - Xiao et al. 2014, Bioinformatics 30(6):801-807 — pi-value
################################################################################

# 0. SETUP ====================================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(patchwork)
  library(ComplexHeatmap)
  library(fgsea)
  library(msigdbr)
  library(rrvgo)
  library(GOSemSim)
  library(org.Hs.eg.db)
  library(grid)
  library(ggsignif)
  library(vegan)
})

setwd(rprojroot::find_rstudio_root_file())

NORM_FILE <- "01_normalization/c_data/01_normalized.csv"
IMP_FILE  <- "02_Imputation/c_data/01_imputed.csv"
DEP_FILE  <- "03_DEP/c_data/combined_results.csv"

RPT_DIR <- "04_Figures/F1/b_reports"
DAT_DIR <- "04_Figures/F1/c_data"
dir.create(RPT_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(DAT_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(RPT_DIR, "supplementary"), showWarnings = FALSE)

CONTRASTS <- c("Aging", "Training_Young", "Training_Old", "Interaction")
CTR_AXIS  <- c(Aging = "Aging", Training_Young = "Training\n(Young)",
               Training_Old = "Training\n(Old)", Interaction = "Interaction")
CTR_FACET <- c(Aging = "Aging", Training_Young = "Training (Young)",
               Training_Old = "Training (Old)", Interaction = "Interaction")
CTR_SHORT <- c(Aging = "Aging", Training_Young = "Tr. (Y)",
               Training_Old = "Tr. (O)", Interaction = "Inter.")

# Palettes
GROUP_FILL <- c(Young_Pre = alpha("#4393C3", 0.5), Young_Post = "#4393C3",
                Old_Pre   = alpha("#D6604D", 0.5), Old_Post   = "#D6604D")
CONTRAST_COLORS <- c(Aging = "#4CAF50", Training_Young = "#E05A4E",
                     Training_Old = "#5DA5DA", Interaction = "#9B7FBF")
DIR_COLORS <- c(Up = "#D6604D", Down = "#4393C3")
DB_COLORS  <- c(Hallmark = "#AA336A", "GO:BP" = "#00796B",
                "GO:CC" = "#26A69A", "GO:MF" = "#CD5C5C")

# PCA: 4 distinct group colors + shapes (circle=Pre, triangle=Post)
PCA_COLORS <- c(Young_Pre = "#93C4DE", Young_Post = "#2166AC",
                Old_Pre   = "#F4A582", Old_Post   = "#B2182B")
PCA_SHAPES <- c(Young_Pre = 16, Young_Post = 17, Old_Pre = 16, Old_Post = 17)

# Key constants centralized in palettes.R
source("04_Figures/shared/palettes.R")

THEME_PUB <- theme_bw(base_size = 8) +
  theme(plot.title       = element_text(face = "bold", size = 9),
        plot.subtitle    = element_text(size = 6.5, color = "grey30", face = "italic"),
        strip.background = element_blank(),
        strip.text       = element_text(face = "bold", size = 6.5),
        legend.key.size  = unit(3, "mm"))

# 1. LOAD DATA ================================================================

norm_df <- read_csv(NORM_FILE, show_col_types = FALSE)
imp_df  <- read_csv(IMP_FILE,  show_col_types = FALSE)
dep_df  <- read_csv(DEP_FILE,  show_col_types = FALSE)

ann_cols   <- c("uniprot_id", "protein", "gene", "description")
samp_names <- setdiff(names(norm_df), ann_cols)

meta <- tibble(sample_id = samp_names) |>
  mutate(
    prefix   = str_extract(sample_id, "^[A-Z]+"),
    subj_num = str_extract(sample_id, "S\\d+"),
    time     = str_extract(sample_id, "(Pre|Post)$"),
    age      = if_else(str_detect(prefix, "^O"), "Old", "Young"),
    subject  = paste0(prefix, "_", subj_num),
    group    = paste(age, time, sep = "_")
  )
meta$age   <- factor(meta$age,  levels = c("Young", "Old"))
meta$time  <- factor(meta$time, levels = c("Pre", "Post"))
meta$group <- factor(meta$group,
                     levels = c("Young_Pre", "Young_Post", "Old_Pre", "Old_Post"))

cat(sprintf("Loaded: %d proteins, %d samples, %d DEP rows\n",
            nrow(norm_df), length(samp_names), nrow(dep_df)))

# Read winning imputation method from summary
imp_summary <- readLines("02_Imputation/c_data/imputation_summary.txt")
BEST_IMP_METHOD <- toupper(trimws(sub(".*=\\s*", "", grep("^best_method", imp_summary, value = TRUE))))
if (length(BEST_IMP_METHOD) == 0) BEST_IMP_METHOD <- "IMPUTED"

# 2. PANEL A — CV% Violins ====================================================

# CV on linear (not log) scale per Brenes 2024
lin_mat <- 2^as.matrix(norm_df[, samp_names])

cv_list <- lapply(levels(meta$group), function(g) {
  idx <- meta$sample_id[meta$group == g]
  sub <- lin_mat[, idx, drop = FALSE]
  cv_pct <- apply(sub, 1, function(x) {
    x <- x[!is.na(x)]
    if (length(x) < 2) return(NA_real_)
    sd(x) / mean(x) * 100
  })
  tibble(group = g, cv = cv_pct)
})
cv_df <- bind_rows(cv_list) |> filter(!is.na(cv))
cv_df$group <- factor(cv_df$group,
                      levels = c("Young_Pre", "Young_Post", "Old_Pre", "Old_Post"))

cv_med <- cv_df |> group_by(group) |> summarise(med = median(cv), .groups = "drop")

# Wilcoxon brackets — only significant (p < 0.05)
bracket_comps <- list(c("Young_Pre", "Young_Post"), c("Old_Pre", "Old_Post"),
                      c("Young_Pre", "Old_Pre"),    c("Young_Post", "Old_Post"))
bracket_pvals <- sapply(bracket_comps, function(pair)
  wilcox.test(cv_df$cv[cv_df$group == pair[1]],
              cv_df$cv[cv_df$group == pair[2]])$p.value)
sig_idx <- which(bracket_pvals < 0.05)

cv_ymax      <- quantile(cv_df$cv, 0.99)
bracket_base <- cv_ymax * 0.85
bracket_step <- cv_ymax * 0.10

pA <- ggplot(cv_df, aes(x = group, y = cv, fill = group)) +
  geom_violin(alpha = 0.7, linewidth = 0.3, scale = "width") +
  geom_boxplot(width = 0.15, outlier.size = 0.3, linewidth = 0.3, fill = "white") +
  geom_hline(yintercept = 25, linetype = "dashed", color = "grey40", linewidth = 0.3) +
  geom_label(data = cv_med, aes(x = group, y = med, label = sprintf("%.0f%%", med)),
             vjust = -0.3, size = 2.5, fontface = "bold", fill = alpha("white", 0.8),
             linewidth = 0.2, label.padding = unit(1, "pt")) +
  scale_fill_manual(values = GROUP_FILL) +
  scale_x_discrete(labels = c("Young Pre", "Young Post", "Old Pre", "Old Post")) +
  labs(title = "A  Inter-Individual Variability (CV%)",
       subtitle = "Protein-level CV on cycloess-normalized intensities (linear scale)",
       x = NULL, y = "CV (%)") +
  THEME_PUB + theme(legend.position = "none")

if (length(sig_idx) > 0) {
  sig_ypos <- bracket_base + (seq_along(sig_idx) - 1) * bracket_step
  for (i in seq_along(sig_idx)) {
    pA <- pA + geom_signif(
      comparisons = list(bracket_comps[[sig_idx[i]]]),
      y_position = sig_ypos[i], tip_length = 0.01,
      test = "wilcox.test", textsize = KEY_TEXT, size = 0.3,
      map_signif_level = c("***" = 0.001, "**" = 0.01, "*" = 0.05))
  }
  pA <- pA + coord_cartesian(ylim = c(0, max(sig_ypos) + bracket_step * 1.5))
} else {
  pA <- pA + coord_cartesian(ylim = c(0, bracket_base + bracket_step))
}

cat("Panel A done\n")

# 3. PANEL B — logFC Density ===================================================

lfc_long <- dep_df |>
  dplyr::select(gene, starts_with("logFC_")) |>
  pivot_longer(starts_with("logFC_"), names_to = "contrast", values_to = "logFC") |>
  mutate(contrast = str_remove(contrast, "logFC_")) |>
  filter(!is.na(logFC), contrast %in% c("Aging", "Training_Young", "Training_Old"))

lfc_long$contrast <- factor(lfc_long$contrast,
                            levels = c("Aging", "Training_Young", "Training_Old"))

lfc_stats <- lfc_long |>
  group_by(contrast) |>
  summarise(
    med_abs_lfc = median(abs(logFC)),
    n_above_05  = sum(abs(logFC) > 0.5),
    .groups = "drop"
  )

# Pairwise Wilcoxon on |logFC| between contrasts (BH-adjusted)
pw_lfc <- pairwise.wilcox.test(abs(lfc_long$logFC), lfc_long$contrast,
                                p.adjust.method = "BH")

.lookup_pw <- function(a, b, mat) {
  if (a %in% rownames(mat) && b %in% colnames(mat) && !is.na(mat[a, b])) return(mat[a, b])
  if (b %in% rownames(mat) && a %in% colnames(mat) && !is.na(mat[b, a])) return(mat[b, a])
  NA_real_
}

lfc_stats$label <- sapply(seq_len(nrow(lfc_stats)), function(i) {
  ctr <- as.character(lfc_stats$contrast[i])
  others <- setdiff(c("Aging", "Training_Young", "Training_Old"), ctr)
  pw_lines <- sapply(others, function(o) {
    p <- .lookup_pw(ctr, o, pw_lfc$p.value)
    p_str <- if (is.na(p)) "p = NA" else if (p < 0.001) "p < 0.001" else sprintf("p = %.3f", p)
    sprintf("vs %s: %s", CTR_SHORT[o], p_str)
  })
  sprintf("Med.|logFC| = %.2f, n(>0.5) = %d\n%s",
          lfc_stats$med_abs_lfc[i], lfc_stats$n_above_05[i],
          paste(pw_lines, collapse = "\n"))
})

lfc_binwidth <- 4 / 50

pB <- ggplot(lfc_long, aes(x = logFC, fill = contrast)) +
  geom_histogram(bins = 50, color = "grey40", linewidth = 0.1, alpha = 0.85) +
  geom_density(aes(y = after_stat(count) * lfc_binwidth),
               alpha = 0.15, linewidth = 0.5, color = "grey20") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey40", linewidth = 0.3) +
  geom_text(data = lfc_stats, aes(label = label),
            x = Inf, y = Inf, hjust = 1.05, vjust = 1.2,
            size = KEY_TEXT, fontface = "bold", inherit.aes = FALSE, color = "grey25") +
  coord_cartesian(xlim = c(-2, 2)) +
  facet_wrap(~ contrast, ncol = 1, scales = "free_y",
             labeller = labeller(contrast = CTR_FACET)) +
  scale_fill_manual(values = CONTRAST_COLORS[c("Aging", "Training_Young", "Training_Old")]) +
  labs(title = "B  Effect Size Distribution (Log2FC)",
       subtitle = "Log2 fold-changes from limma (Interaction excluded)",
       x = expression(log[2]~FC), y = "Count") +
  THEME_PUB + theme(legend.position = "none")

# KS + Fligner distributional stats annotation
blunt_file <- "03_DEP/c_data/blunting_diagnostics.csv"
if (file.exists(blunt_file)) {
  blunt <- read_csv(blunt_file, show_col_types = FALSE)
  ks_row  <- blunt[blunt$test == "Kolmogorov-Smirnov", ]
  fl_row  <- blunt[blunt$test == "Fligner-Killeen", ]

  dist_label <- sprintf(
    "Young vs Old |logFC|:  KS D = %.2f, p < 10^{-30};  Fligner \u03C7\u00B2 = %.0f, p < 10^{-48}\n\u2192 Young shows larger, more variable training response",
    ks_row$statistic, fl_row$statistic)

  pB <- pB +
    labs(caption = dist_label) +
    theme(plot.caption = element_text(size = 5.5, color = "grey30",
                                       hjust = 0, face = "italic",
                                       margin = margin(t = 4)))
  cat(sprintf("  Added KS/Fligner annotation: D=%.3f, chi-sq=%.1f\n",
              ks_row$statistic, fl_row$statistic))
} else {
  cat("  blunting_diagnostics.csv not found — skipping annotation\n")
}

cat("Panel B done\n")

# 4. PANEL C — PCA =============================================================

imp_mat <- as.matrix(imp_df[, samp_names])
rownames(imp_mat) <- imp_df$gene

pca <- prcomp(t(imp_mat), center = TRUE, scale. = TRUE)
var_pct <- round(100 * summary(pca)$importance[2, 1:2], 1)

pca_df <- as.data.frame(pca$x[, 1:2]) |>
  mutate(sample_id = rownames(pca$x)) |>
  left_join(meta, by = "sample_id")

# PERMANOVA
dist_mat <- dist(t(imp_mat))
set.seed(42)
perm_res <- adonis2(dist_mat ~ age * time, data = meta,
                    permutations = 999, by = "terms")

.fmt_p <- function(p) {
  if (is.na(p)) return("NA")
  if (p < 0.001) return("< 0.001")
  sprintf("%.3f", p)
}
perm_terms <- c("age", "time", "age:time")
perm_r2 <- perm_res[perm_terms, "R2"]
perm_pv <- perm_res[perm_terms, "Pr(>F)"]
perm_label <- sprintf(
  "PERMANOVA (999 perm.)\nAge  R\u00b2 = %.3f,  p = %s\nTime  R\u00b2 = %.3f,  p = %s\nAge\u00d7Time  R\u00b2 = %.3f,  p = %s",
  perm_r2[1], .fmt_p(perm_pv[1]),
  perm_r2[2], .fmt_p(perm_pv[2]),
  perm_r2[3], .fmt_p(perm_pv[3]))
cat("PERMANOVA results:\n"); print(perm_res)

pC <- ggplot(pca_df, aes(x = PC1, y = PC2, color = group, shape = group)) +
  stat_ellipse(aes(fill = group), geom = "polygon",
               alpha = 0.10, level = 0.80, show.legend = FALSE) +
  stat_ellipse(aes(group = group), level = 0.80, linewidth = 0.4,
               linetype = "dashed", show.legend = FALSE) +
  geom_point(size = 2.0, alpha = 0.85) +
  annotate("text", x = -Inf, y = Inf, label = perm_label,
           hjust = -0.05, vjust = 1.2,
           size = KEY_TEXT, fontface = "bold", color = "grey25") +
  scale_color_manual(values = PCA_COLORS) +
  scale_fill_manual(values = PCA_COLORS) +
  scale_shape_manual(values = PCA_SHAPES) +
  labs(title = "C  Principal Component Analysis (PCA)",
       subtitle = sprintf("%s proteins (HPA-filtered, cycloess, %s-imputed); n = %d samples",
                          format(nrow(imp_df), big.mark = ","), BEST_IMP_METHOD, nrow(meta)),
       x = sprintf("PC1 (%.1f%%)", var_pct[1]),
       y = sprintf("PC2 (%.1f%%)", var_pct[2])) +
  THEME_PUB + theme(legend.position = "none")

# PCA key — 2x2 grid with group-colored backgrounds
pca_leg <- tibble(
  col   = c(0, 0.50, 0, 0.50),
  row   = c(1, 1, 0, 0),
  label = c("Young Pre", "Young Post", "Old Pre", "Old Post"),
  group = c("Young_Pre", "Young_Post", "Old_Pre", "Old_Post"),
  shape = c(16, 17, 16, 17)
) |> mutate(fill = PCA_COLORS[group])

pC_key <- ggplot(pca_leg) +
  geom_rect(aes(xmin = col, xmax = col + 0.38, ymin = row - 0.38, ymax = row + 0.38),
            fill = alpha(pca_leg$fill, 0.30), color = "grey70", linewidth = 0.2) +
  geom_point(aes(x = col + 0.05, y = row),
             shape = pca_leg$shape, color = pca_leg$fill, size = 2) +
  geom_text(aes(x = col + 0.10, y = row, label = label),
            hjust = 0, size = KEY_TEXT, fontface = "bold", color = "grey15") +
  annotate("text", x = 0, y = 1.75, label = "Group \u00d7 Time:",
           hjust = 0, size = KEY_TITLE, fontface = "bold", color = "grey25") +
  scale_x_continuous(limits = c(-0.02, 0.92)) +
  scale_y_continuous(limits = c(-0.45, 1.92)) +
  theme_void() +
  theme(plot.margin = margin(0, 5.5, 0, 5.5))

cat("Panel C done\n")

# 5. PANEL D — DEPs per Contrast ===============================================

SET_LABELS <- c(Aging = "Aging", Training_Young = "Training (Young)",
                Training_Old = "Training (Old)", Interaction = "Interaction")
all_genes <- unique(dep_df$gene[!is.na(dep_df$gene)])

# Pi-score significant sets (shared with Panel E)
sig_sets <- list()
dir_map  <- list()

for (ctr in CONTRASTS) {
  pi_vals  <- dep_df[[paste0("pi_score_", ctr)]]
  lfc_vals <- dep_df[[paste0("logFC_", ctr)]]
  is_sig   <- !is.na(pi_vals) & pi_vals < 0.05
  sig_sets[[ctr]] <- dep_df$gene[is_sig]
  dir_map[[ctr]]  <- setNames(ifelse(lfc_vals[is_sig] > 0, "Up", "Down"),
                               dep_df$gene[is_sig])
}

# Compute FDR and Pi totals for subtitle
pi_total  <- sum(sapply(sig_sets, length))
fdr_total <- sum(sapply(CONTRASTS, function(ctr) {
  fdr_col <- paste0("adj.P.Val_", ctr)
  if (fdr_col %in% names(dep_df)) sum(dep_df[[fdr_col]] < 0.05, na.rm = TRUE) else 0
}))

frac_list <- lapply(CONTRASTS, function(ctr) {
  fdr_col <- paste0("adj.P.Val_", ctr)
  p_col   <- paste0("P.Value_", ctr)
  tibble(
    contrast  = SET_LABELS[ctr],
    threshold = c("p < 0.05", "q < 0.05", "\u03A0 < 0.05"),
    n = c(sum(!is.na(dep_df[[p_col]])   & dep_df[[p_col]]   < 0.05),
          sum(!is.na(dep_df[[fdr_col]]) & dep_df[[fdr_col]] < 0.05),
          length(sig_sets[[ctr]]))
  )
})
frac_df <- bind_rows(frac_list) |>
  mutate(
    contrast  = factor(contrast,
                       levels = rev(c("Aging", "Training (Young)",
                                      "Training (Old)", "Interaction"))),
    threshold = factor(threshold, levels = c("p < 0.05", "q < 0.05", "\u03A0 < 0.05")),
    pct       = 100 * n / length(all_genes),
    fill_key  = paste(contrast, threshold, sep = "___")
  ) |>
  filter(n > 1)

SET_DISPLAY_COLORS <- c("Aging"            = unname(CONTRAST_COLORS["Aging"]),
                        "Training (Young)" = unname(CONTRAST_COLORS["Training_Young"]),
                        "Training (Old)"   = unname(CONTRAST_COLORS["Training_Old"]),
                        "Interaction"      = unname(CONTRAST_COLORS["Interaction"]))

FRAC_FILL <- c()
for (cname in names(SET_DISPLAY_COLORS)) {
  col <- unname(SET_DISPLAY_COLORS[cname])
  FRAC_FILL[paste(cname, "p < 0.05",      sep = "___")] <- adjustcolor(col, alpha.f = 0.25)
  FRAC_FILL[paste(cname, "q < 0.05",      sep = "___")] <- adjustcolor(col, alpha.f = 0.55)
  FRAC_FILL[paste(cname, "\u03A0 < 0.05", sep = "___")] <- col
}

THRESH_LABEL <- c("p < 0.05" = "p \u2264 0.05", "q < 0.05" = "FDR \u2264 0.05",
                  "\u03A0 < 0.05" = "\u03A0 \u2264 0.05")

label_df <- frac_df |>
  group_by(contrast) |> arrange(contrast, threshold) |>
  mutate(label    = THRESH_LABEL[as.character(threshold)],
         label_y  = (lead(pct, default = 0) + pct) / 2,
         text_col = if_else(threshold == "p < 0.05", "grey20", "white")) |>
  ungroup()

pD <- ggplot(frac_df, aes(x = contrast, y = pct, fill = fill_key)) +
  annotate("rect", xmin = 3.5, xmax = 4.5, ymin = -Inf, ymax = Inf,
           fill = CONTRAST_COLORS["Aging"], alpha = 0.20,
           color = "grey70", linewidth = 0.2) +
  annotate("rect", xmin = 2.5, xmax = 3.5, ymin = -Inf, ymax = Inf,
           fill = CONTRAST_COLORS["Training_Young"], alpha = 0.20,
           color = "grey70", linewidth = 0.2) +
  annotate("rect", xmin = 1.5, xmax = 2.5, ymin = -Inf, ymax = Inf,
           fill = CONTRAST_COLORS["Training_Old"], alpha = 0.20,
           color = "grey70", linewidth = 0.2) +
  annotate("rect", xmin = 0.5, xmax = 1.5, ymin = -Inf, ymax = Inf,
           fill = CONTRAST_COLORS["Interaction"], alpha = 0.20,
           color = "grey70", linewidth = 0.2) +
  geom_col(position = "identity", width = 0.75, color = "black", linewidth = 0.3) +
  geom_text(data = label_df,
            aes(x = contrast, y = label_y, label = label, color = I(text_col)),
            inherit.aes = FALSE, hjust = 0.5, size = 1.8, fontface = "bold") +
  scale_fill_manual(values = FRAC_FILL) +
  scale_y_continuous(trans = scales::trans_new("pow0.7",
                       transform = function(x) ifelse(x >= 0, x^0.7, 0),
                       inverse   = function(x) ifelse(x >= 0, x^(1/0.7), 0)),
                     expand = expansion(mult = c(0, 0.05)),
                     breaks = c(0, 5, 10, 15, 20, 25, 30)) +
  coord_flip() +
  labs(title = "D  DEPs per Contrast",
       subtitle = sprintf("Fraction of %s filtered proteins",
                          format(length(all_genes), big.mark = ",")),
       x = NULL, y = "% of proteome") +
  THEME_PUB + theme(legend.position = "none",
                    axis.text.y = element_text(face = "bold", size = 7))

cat("Panel D done\n")

# 6. KEYS =====================================================================

# Unified key — Contrast Definitions (left) + Direction (right)
# Single ggplot ensures perfect alignment and consistent box sizing.

ctr_leg <- tibble(
  y     = c(1.35, 0.90, 0.45, 0),
  label = c("Aging:  Old_Pre \u2212 Young_Pre",
            "Tr. (Young):  Young_Post \u2212 Young_Pre",
            "Tr. (Old):  Old_Post \u2212 Old_Pre",
            "Interaction:  (Old_Post\u2212Old_Pre) \u2212 (Young_Post\u2212Young_Pre)"),
  fill  = unname(CONTRAST_COLORS[c("Aging", "Training_Young",
                                    "Training_Old", "Interaction")])
)

dir_leg <- tibble(
  y     = c(0.90, 0.45),
  label = c("Up", "Down"),
  fill  = unname(DIR_COLORS[c("Up", "Down")])
)

# X layout: contrast boxes 0..5.2, direction boxes 5.6..6.5
CTR_XMAX <- 5.2    # wide enough for the Interaction text
DIR_XMIN <- 5.6    # small gap between sections
DIR_XMAX <- 6.5    # direction boxes narrower (only "Up"/"Down")
BOX_HALF <- 0.18   # half-height of each box (same for both key types)

pKeys <- ggplot() +
  # -- Contrast definition boxes (left) --
  geom_rect(data = ctr_leg,
            aes(xmin = 0, xmax = CTR_XMAX, ymin = y - BOX_HALF, ymax = y + BOX_HALF),
            fill = alpha(ctr_leg$fill, 0.25), color = "grey70", linewidth = 0.2) +
  geom_text(data = ctr_leg,
            aes(x = 0.06, y = y, label = label),
            hjust = 0, size = 2.8, fontface = "bold", color = "grey15") +

  # -- Direction boxes (right, same row height) --
  geom_rect(data = dir_leg,
            aes(xmin = DIR_XMIN, xmax = DIR_XMAX, ymin = y - BOX_HALF, ymax = y + BOX_HALF),
            fill = dir_leg$fill, color = "grey70", linewidth = 0.2) +
  geom_text(data = dir_leg,
            aes(x = (DIR_XMIN + DIR_XMAX) / 2, y = y, label = label),
            color = "white", fontface = "bold", size = 2.8) +

  # -- Section titles (same y, same font) --
  annotate("text", x = 0, y = 1.62,
           label = "Contrast Definitions:",
           hjust = 0, size = KEY_TITLE, fontface = "bold", color = "grey25") +
  annotate("text", x = DIR_XMIN, y = 1.18,
           label = "Direction:",
           hjust = 0, size = KEY_TITLE, fontface = "bold", color = "grey25") +

  scale_x_continuous(limits = c(-0.05, 6.8), expand = c(0, 0)) +
  scale_y_continuous(limits = c(-0.25, 1.80), expand = c(0, 0)) +
  theme_void() +
  theme(plot.margin = margin(0, 5.5, 0, 5.5))

cat("Keys done\n")

# 7. PANEL E — UpSet (pi-score overlap) =======================================

bin_mat <- sapply(sig_sets, function(s) as.integer(all_genes %in% s))
rownames(bin_mat) <- all_genes
colnames(bin_mat) <- SET_LABELS[colnames(bin_mat)]
label_to_contrast <- setNames(names(SET_LABELS), SET_LABELS)

cm <- make_comb_mat(bin_mat, mode = "intersect")
cs <- comb_size(cm)
keep <- cs > 0 & comb_degree(cm) > 0
cm_sub <- cm[keep]

comb_names_vec    <- comb_name(cm_sub)
n_comb            <- length(comb_names_vec)
set_names_ordered <- set_name(cm_sub)
up_counts   <- numeric(n_comb)
down_counts <- numeric(n_comb)

for (i in seq_len(n_comb)) {
  members <- extract_comb(cm_sub, comb_names_vec[i])
  if (length(members) == 0) next
  bits <- as.logical(as.integer(strsplit(comb_names_vec[i], "")[[1]]))
  active_contrasts <- label_to_contrast[set_names_ordered[bits]]
  gene_dirs <- sapply(members, function(g) {
    dirs <- sapply(active_contrasts, function(ctr) {
      if (g %in% names(dir_map[[ctr]])) dir_map[[ctr]][g] else NA
    })
    dirs <- dirs[!is.na(dirs)]
    if (all(dirs == "Up")) "Up" else if (all(dirs == "Down")) "Down" else "Mixed"
  })
  up_counts[i]   <- sum(gene_dirs == "Up")
  down_counts[i] <- sum(gene_dirs == "Down")
}

display_total <- up_counts + down_counts
keep_display  <- display_total > 0
comb_ord <- which(keep_display)[order(-display_total[keep_display])]
up_ord   <- up_counts[comb_ord]
down_ord <- down_counts[comb_ord]

set_display_order <- c("Aging", "Training (Young)", "Training (Old)", "Interaction")
n_int          <- length(comb_ord)
comb_names_ord <- comb_names_vec[comb_ord]
set_order_ch   <- set_name(cm_sub)
set_y_levels   <- rev(set_display_order)

bar_long <- tibble(
  x         = rep(seq_len(n_int), 2),
  direction = factor(rep(c("Up", "Down"), each = n_int), levels = c("Down", "Up")),
  count     = c(up_ord, down_ord)
)

dot_df <- expand_grid(x = seq_len(n_int), set = set_display_order)
dot_df$active <- vapply(seq_len(nrow(dot_df)), function(r) {
  bits <- as.integer(strsplit(comb_names_ord[dot_df$x[r]], "")[[1]])
  as.logical(bits[match(dot_df$set[r], set_order_ch)])
}, logical(1))
dot_df$set  <- factor(dot_df$set, levels = set_y_levels)
dot_df$ynum <- as.numeric(dot_df$set)

seg_list <- vector("list", n_int)
for (i in seq_len(n_int)) {
  bits <- as.logical(as.integer(strsplit(comb_names_ord[i], "")[[1]]))
  ypos <- match(set_order_ch[bits], set_y_levels)
  if (length(ypos) > 1)
    seg_list[[i]] <- tibble(x = i, ymin = min(ypos), ymax = max(ypos))
}
seg_df <- bind_rows(seg_list)

stripe_fills <- adjustcolor(unname(SET_DISPLAY_COLORS[set_y_levels]), alpha.f = 0.20)

bar_bg_list <- lapply(seq_len(n_int), function(i) {
  bits <- as.logical(as.integer(strsplit(comb_names_ord[i], "")[[1]]))
  if (sum(bits) != 1) return(NULL)
  tibble(xmin = i - 0.5, xmax = i + 0.5,
         fill = unname(SET_DISPLAY_COLORS[set_order_ch[bits]]))
})
bar_bg <- bind_rows(compact(bar_bg_list))

pE_bars <- ggplot(bar_long, aes(x, count, fill = direction)) +
  geom_rect(data = bar_bg,
            aes(xmin = xmin, xmax = xmax, ymin = -Inf, ymax = Inf),
            fill = bar_bg$fill, alpha = 0.20,
            color = "grey70", linewidth = 0.2, inherit.aes = FALSE) +
  geom_col(position = position_dodge(width = 0.7), width = 0.6,
           color = "black", linewidth = 0.3) +
  geom_text(aes(label = ifelse(count > 0, count, ""), y = count / 2),
            position = position_dodge(width = 0.7), vjust = 0.5,
            size = KEY_TEXT, color = "white", fontface = "bold") +
  scale_fill_manual(values = c(Up = unname(DIR_COLORS["Up"]),
                                Down = unname(DIR_COLORS["Down"]))) +
  scale_x_continuous(expand = expansion(add = 0.6)) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.05))) +
  labs(y = "Intersection\nsize",
       title = "E  Contrast Overlap (UpSet)",
       subtitle = sprintf("DEPs by pi < 0.05; bars split by direction | Pi: %d / FDR: %d total",
                          pi_total, fdr_total)) +
  THEME_PUB +
  theme(axis.text.x  = element_blank(), axis.ticks.x = element_blank(),
        axis.title.x = element_blank(),
        panel.grid.major.x = element_blank(), panel.grid.minor = element_blank(),
        legend.position = "none",
        plot.margin = margin(2, 5.5, 0, 5.5))

pE_dots <- ggplot() +
  annotate("rect",
           xmin = 0.4, xmax = n_int + 0.6,
           ymin = seq_along(set_y_levels) - 0.5,
           ymax = seq_along(set_y_levels) + 0.5,
           fill = stripe_fills) +
  {if (nrow(seg_df) > 0)
    geom_segment(data = seg_df,
                 aes(x = x, xend = x, y = ymin, yend = ymax),
                 linewidth = 0.6, color = "grey25")} +
  geom_point(data = dot_df |> filter(!active),
             aes(x = x, y = ynum), color = "grey78", size = 1.2) +
  geom_point(data = dot_df |> filter(active),
             aes(x = x, y = ynum), color = "grey15", size = 1.2) +
  scale_x_continuous(expand = expansion(add = 0.6)) +
  scale_y_continuous(breaks = seq_along(set_y_levels), labels = set_y_levels,
                     expand = expansion(add = 0.15)) +
  labs(x = NULL, y = NULL) +
  THEME_PUB +
  theme(axis.text.x  = element_blank(), axis.ticks = element_blank(),
        panel.grid   = element_blank(), panel.border = element_blank(),
        axis.text.y  = element_text(size = 6.5),
        plot.margin  = margin(0, 5.5, 0, 5.5))

cat("Panel E done\n")

# 8. PANEL F — fgsea Pathway Enrichment ========================================

go_bp_full <- msigdbr(species = "Homo sapiens", collection = "C5", subcollection = "GO:BP")
go_cc_full <- msigdbr(species = "Homo sapiens", collection = "C5", subcollection = "GO:CC")
go_mf_full <- msigdbr(species = "Homo sapiens", collection = "C5", subcollection = "GO:MF")

db_list <- list(
  Hallmark = msigdbr(species = "Homo sapiens", collection = "H") |>
    dplyr::select(gs_name, gene_symbol),
  "GO:BP" = go_bp_full |> dplyr::select(gs_name, gene_symbol),
  "GO:CC" = go_cc_full |> dplyr::select(gs_name, gene_symbol),
  "GO:MF" = go_mf_full |> dplyr::select(gs_name, gene_symbol)
)

go_id_lookup <- bind_rows(
  go_bp_full |> dplyr::select(gs_name, gs_exact_source),
  go_cc_full |> dplyr::select(gs_name, gs_exact_source),
  go_mf_full |> dplyr::select(gs_name, gs_exact_source)
) |> distinct(gs_name, .keep_all = TRUE) |> deframe()

orgdb <- org.Hs.eg.db
set.seed(42)
fgsea_all <- list()

for (ctr in CONTRASTS) {
  stats <- dep_df[[paste0("t_", ctr)]]
  names(stats) <- dep_df$gene
  stats <- sort(stats[!is.na(stats) & is.finite(stats)], decreasing = TRUE)

  for (db_name in names(db_list)) {
    pathway_list <- split(db_list[[db_name]]$gene_symbol, db_list[[db_name]]$gs_name)
    res <- fgseaMultilevel(pathways = pathway_list, stats = stats,
                           minSize = 15, maxSize = 200, nPermSimple = 10000, eps = 0)
    res$contrast <- ctr
    res$database <- db_name
    fgsea_all[[paste(ctr, db_name, sep = "_")]] <- as.data.frame(res)
  }
  cat(sprintf("  fgsea done: %s\n", ctr))
}

fgsea_combined <- bind_rows(fgsea_all)

fgsea_export <- fgsea_combined |>
  mutate(leadingEdge = sapply(leadingEdge, paste, collapse = ";")) |>
  arrange(database, contrast, padj)
write_csv(fgsea_export, file.path(DAT_DIR, "fgsea_results.csv"))
cat(sprintf("Saved fgsea_results.csv: %d rows\n", nrow(fgsea_export)))

# Reduce GO terms with rrvgo (threshold 0.7)
.reduce_go <- function(df, ont_short) {
  sig <- df |> filter(padj < 0.05)
  if (nrow(sig) == 0) return(sig)
  tryCatch({
    go_ids <- go_id_lookup[sig$pathway]
    valid  <- !is.na(go_ids)
    if (sum(valid) < 2) return(sig)
    scores <- setNames(-log10(sig$padj[valid]), go_ids[valid])
    simMat <- calculateSimMatrix(go_ids[valid], orgdb = orgdb,
                                 ont = ont_short, method = "Rel")
    keep_ids <- intersect(names(scores), rownames(simMat))
    if (length(keep_ids) < 2) return(sig)
    reduced <- reduceSimMatrix(simMat[keep_ids, keep_ids],
                               scores = scores[keep_ids],
                               threshold = 0.7, orgdb = orgdb)
    goid_to_name <- setNames(sig$pathway[valid], go_ids[valid])
    keep_pw <- goid_to_name[unique(reduced$parent)]
    if (length(keep_pw) > 0) sig |> filter(pathway %in% keep_pw) else sig
  }, error = function(e) { message("rrvgo skipped: ", e$message); sig })
}

ont_map <- c("GO:BP" = "BP", "GO:CC" = "CC", "GO:MF" = "MF")
count_df <- bind_rows(lapply(CONTRASTS, function(ctr) {
  bind_rows(lapply(names(db_list), function(db) {
    sub <- fgsea_combined |> filter(contrast == ctr, database == db)
    sig <- if (db %in% names(ont_map)) .reduce_go(sub, ont_map[db])
           else sub |> filter(padj < 0.05)
    tibble(contrast = ctr, database = db,
           direction = c("Up", "Down"),
           count = c(sum(sig$NES > 0), sum(sig$NES < 0)))
  }))
}))
count_df$contrast  <- factor(count_df$contrast, levels = CONTRASTS)
count_df$database  <- factor(count_df$database, levels = names(DB_COLORS))
count_df$direction <- factor(count_df$direction, levels = c("Up", "Down"))

pF <- ggplot(count_df, aes(x = contrast, y = count, fill = direction)) +
  annotate("rect", xmin = 0.5, xmax = 1.5, ymin = -Inf, ymax = Inf,
           fill = CONTRAST_COLORS["Aging"], alpha = 0.20,
           color = "grey70", linewidth = 0.2) +
  annotate("rect", xmin = 1.5, xmax = 2.5, ymin = -Inf, ymax = Inf,
           fill = CONTRAST_COLORS["Training_Young"], alpha = 0.20,
           color = "grey70", linewidth = 0.2) +
  annotate("rect", xmin = 2.5, xmax = 3.5, ymin = -Inf, ymax = Inf,
           fill = CONTRAST_COLORS["Training_Old"], alpha = 0.20,
           color = "grey70", linewidth = 0.2) +
  annotate("rect", xmin = 3.5, xmax = 4.5, ymin = -Inf, ymax = Inf,
           fill = CONTRAST_COLORS["Interaction"], alpha = 0.20,
           color = "grey70", linewidth = 0.2) +
  geom_col(position = position_dodge(width = 0.7), width = 0.6,
           color = "black", linewidth = 0.3) +
  geom_text(aes(label = ifelse(count > 0, count, ""), y = count / 2),
            position = position_dodge(width = 0.7), vjust = 0.5, size = KEY_TEXT,
            color = "white", fontface = "bold") +
  facet_grid(database ~ ., scales = "free_y") +
  scale_x_discrete(labels = CTR_SHORT) +
  scale_fill_manual(values = DIR_COLORS) +
  labs(title = "F  Pathway Enrichment",
       subtitle = "fGSEA (padj < 0.05); GO reduced\nvia rrvgo (threshold 0.7)",
       x = NULL, y = "Significant pathways") +
  THEME_PUB +
  theme(axis.text.x    = element_text(angle = 45, hjust = 1, size = 6.5),
        legend.position = "none",
        strip.text.y   = element_text(size = 6.5, angle = 0))

cat("Panel F done\n")

# 9. COMPOSITE =================================================================

left_col  <- (pA / pB / pC / pC_key) +
              plot_layout(heights = c(0.20, 0.38, 0.35, 0.07))
mid_col   <- (pD / pE_bars / pE_dots / pKeys) +
              plot_layout(heights = c(0.22, 0.46, 0.18, 0.14))
fig1 <- (left_col | mid_col | pF) +
  plot_layout(widths = c(0.32, 0.43, 0.25)) +
  plot_annotation(
    title = "Proteomic Landscape of Aging and Resistance Training Response",
    subtitle = "Data quality, effect magnitude, contrast overlap, and functional enrichment across the 2 \u00d7 2 factorial design.",
    theme = theme(
      plot.title    = element_text(face = "bold", size = 11, hjust = 0.5),
      plot.subtitle = element_text(size = 8, color = "grey30", hjust = 0.5)
    )
  )

ggsave(file.path(RPT_DIR, "Figure_1.pdf"), fig1,
       width = 380, height = 260, units = "mm", device = pdf)
ggsave(file.path(RPT_DIR, "Figure_1.png"), fig1,
       width = 380, height = 260, units = "mm", dpi = 300)

cat("Saved Figure_1.pdf and Figure_1.png\n")

# 10. S1.7 — pi-score Distributions ===========================================

pi_long <- dep_df |>
  dplyr::select(gene, starts_with("pi_score_")) |>
  pivot_longer(starts_with("pi_score_"), names_to = "contrast", values_to = "pi_score") |>
  mutate(contrast = str_remove(contrast, "pi_score_")) |>
  filter(!is.na(pi_score))

p_hist <- ggplot(pi_long, aes(x = pi_score)) +
  geom_histogram(bins = 50, fill = "grey60", color = "white", linewidth = 0.2) +
  geom_vline(xintercept = 0.05, linetype = "dashed", color = "red", linewidth = 0.4) +
  facet_wrap(~ contrast, scales = "free_y", ncol = 2,
             labeller = labeller(contrast = CTR_AXIS)) +
  labs(x = expression(pi~"-score"), y = "Count") +
  THEME_PUB

pi_ranked <- pi_long |>
  group_by(contrast) |> arrange(pi_score) |>
  mutate(rank = row_number()) |> ungroup()

n_sig <- pi_long |>
  group_by(contrast) |>
  summarise(n = sum(pi_score < 0.05), .groups = "drop")

p_rank <- ggplot(pi_ranked, aes(x = rank, y = pi_score)) +
  geom_point(size = 0.2, alpha = 0.4, color = "grey40") +
  geom_hline(yintercept = 0.05, linetype = "dashed", color = "red", linewidth = 0.4) +
  geom_text(data = n_sig, aes(label = sprintf("n = %d", n)),
            x = Inf, y = 0.10, hjust = 1.1, vjust = 0, size = 2.5, color = "red") +
  facet_wrap(~ contrast, scales = "free_x", ncol = 2,
             labeller = labeller(contrast = CTR_AXIS)) +
  labs(x = "Protein rank", y = expression(pi~"-score")) +
  THEME_PUB

s17 <- p_hist / p_rank +
  plot_annotation(title = "S1.7  Pi-score distributions",
                  theme = theme(plot.title = element_text(face = "bold", size = 10)))

ggsave(file.path(RPT_DIR, "supplementary", "S1_7_pi_score_distributions.pdf"), s17,
       width = 180, height = 200, units = "mm", device = pdf)

cat("Saved S1_7_pi_score_distributions.pdf\n")
cat("Figure 1 pipeline complete.\n")
