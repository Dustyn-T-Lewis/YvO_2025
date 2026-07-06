#!/usr/bin/env Rscript
# F07 Supplementary — Age Discrimination via Multi-Classifier Comparison
# Writes: b_reports/supp/{png,pdf}/panels/SUPP_F07_multivariate_classifier.{png,pdf}
#
# 4 classifiers compared head-to-head:
#   1. Pre-only:  baseline eigengenes
#   2. Post-only: post-training eigengenes
#   3. delta-ME:  training response (Post - Pre)
#   4. Combined:  mean(Pre, Post) per subject — AVERAGED not stacked
#
# Sourced by 01_main_panels.R — expects style.R + figure_supplement_helpers.R
# already loaded.

library(tidyverse)
library(patchwork)
library(pROC)

source("04_Figures_v2/F07/a_script/_f07_helpers.R")

BASE <- "04_Figures_v2/F07"
RPT_PNG <- file.path(BASE, "b_reports", "supp", "png", "panels")
RPT_PDF <- file.path(BASE, "b_reports", "supp", "pdf", "panels")
DAT <- "04_Figures_v2/F06/c_data"
DAT_OUT <- file.path(BASE, "c_data")
dir.create(RPT_PNG, recursive = TRUE, showWarnings = FALSE)
dir.create(RPT_PDF, recursive = TRUE, showWarnings = FALSE)
dir.create(DAT_OUT, recursive = TRUE, showWarnings = FALSE)

F06_SUPP <- "04_Figures_v2/F06/c_data/F06_supplementary.xlsx"
stopifnot(
  "F06 stitcher must run first: missing F06_supplementary.xlsx" =
    file.exists(F06_SUPP)
)

MEs <- read_matrix_sheet(F06_SUPP, "MEs", "sample_id")
me_pre <- read_matrix_sheet(F06_SUPP, "me_pre", "subject_key")
me_post <- read_matrix_sheet(F06_SUPP, "me_post", "subject_key")
subj_age <- read_sheet_df(F06_SUPP, "metadata_subj_age")
common_subj <- read_vector_sheet(F06_SUPP, "common_subj")

pdf_device <- get_pdf_device()

message("SUPP multivariate: multi-classifier age discrimination...")

me_delta <- as.matrix(me_post[common_subj, ]) - as.matrix(me_pre[common_subj, ])

true_labels <- ifelse(
  subj_age$age[match(common_subj, subj_age$subject_key)] == "Old", 1, 0
)
n_subj <- length(common_subj)
n_young <- sum(true_labels == 0)
n_old <- sum(true_labels == 1)

# run_topk_loocv lives in _f07_helpers.R.

# Run all 4 classifiers
set.seed(42)

clf_names <- c("Pre", "Post", "DeltaME", "Combined")
clf_colors <- c(
  Pre = "#D6604D", Post = "#5DA5DA",
  DeltaME = "#4CAF50", Combined = "#9B7FBF"
)

message("  Running Pre-only classifier...")
res_pre <- run_topk_loocv(true_labels, me_pre[common_subj, ])

message("  Running Post-only classifier...")
res_post <- run_topk_loocv(true_labels, me_post[common_subj, ])

message("  Running delta-ME classifier...")
res_delta <- run_topk_loocv(true_labels, me_delta)

message("  Running Combined (averaged Pre+Post) classifier...")
me_avg <- (as.matrix(me_pre[common_subj, ]) + as.matrix(me_post[common_subj, ])) / 2
res_combined <- run_topk_loocv(true_labels, me_avg)

all_res <- list(
  Pre = res_pre, Post = res_post,
  DeltaME = res_delta, Combined = res_combined
)

clf_summary <- map_dfr(clf_names, function(nm) {
  probs <- all_res[[nm]]$probs
  roc_obj <- tryCatch(roc(true_labels, probs, quiet = TRUE),
    error = function(e) NULL
  )
  if (is.null(roc_obj)) {
    return(tibble(
      classifier = nm, auc = 0.5,
      ci_lo = NA_real_, ci_hi = NA_real_
    ))
  }
  ci_obj <- ci.auc(roc_obj)
  acc <- mean((probs > 0.5) == true_labels)
  tibble(
    classifier = nm, auc = as.numeric(auc(roc_obj)),
    ci_lo = ci_obj[1], ci_hi = ci_obj[3],
    accuracy = acc, n = n_subj, n_young = n_young, n_old = n_old
  )
})

best_clf <- clf_summary |>
  arrange(desc(auc)) |>
  slice(1) |>
  pull(classifier)
best_auc <- clf_summary |>
  filter(classifier == best_clf) |>
  pull(auc)

message(sprintf("  AUC summary:"))
for (i in seq_len(nrow(clf_summary))) {
  row <- clf_summary[i, ]
  star <- if (row$classifier == best_clf) " *" else ""
  message(sprintf(
    "    %s: AUC=%.3f [%.3f-%.3f] Acc=%.1f%%%s",
    row$classifier, row$auc, row$ci_lo, row$ci_hi,
    row$accuracy * 100, star
  ))
}

feat_stab <- map_dfr(clf_names, function(nm) {
  freq <- all_res[[nm]]$feature_freq
  tibble(
    classifier = nm, feature = names(freq),
    n_selected = as.integer(freq),
    pct_selected = round(100 * freq / n_subj, 1)
  )
})
write_csv(feat_stab, file.path(DAT_OUT, "panel_A_feature_stability.csv"))

# Permutation test on best classifier (1000 perms); fast_loocv_auc from helpers
message(sprintf("  Permutation test on best classifier (%s)...", best_clf))

best_me <- switch(best_clf,
  Pre      = me_pre[common_subj, ],
  Post     = me_post[common_subj, ],
  DeltaME  = me_delta,
  Combined = me_avg
)
best_k_median <- as.integer(median(all_res[[best_clf]]$best_k_log))
message(sprintf("  Using fixed k=%d (median) for permutations", best_k_median))

n_perm <- 1000
null_aucs <- numeric(n_perm)

for (i in seq_len(n_perm)) {
  shuffled <- sample(true_labels)
  null_aucs[i] <- fast_loocv_auc(shuffled, best_me, best_k_median)
}

perm_pvalue <- (sum(null_aucs >= best_auc) + 1) / (n_perm + 1)
null_mean <- mean(null_aucs)
null_sd <- sd(null_aucs)

message(sprintf(
  "  %s AUC=%.3f, perm p=%.4f, null mean=%.3f +/-%.3f",
  best_clf, best_auc, perm_pvalue, null_mean, null_sd
))

# Also run permutation on all 4 classifiers (100 perms each for speed)
message("  Quick permutation test on all classifiers (100 perms)...")
n_perm_quick <- 100
all_perm_p <- setNames(numeric(length(clf_names)), clf_names)

for (nm in clf_names) {
  me_nm <- switch(nm,
    Pre = me_pre[common_subj, ],
    Post = me_post[common_subj, ],
    DeltaME = me_delta,
    Combined = me_avg
  )
  k_nm <- as.integer(median(all_res[[nm]]$best_k_log))
  obs_auc <- clf_summary |>
    filter(classifier == nm) |>
    pull(auc)
  quick_nulls <- numeric(n_perm_quick)
  for (j in seq_len(n_perm_quick)) {
    shuf <- sample(true_labels)
    quick_nulls[j] <- fast_loocv_auc(shuf, me_nm, k_nm)
  }
  all_perm_p[nm] <- (sum(quick_nulls >= obs_auc) + 1) / (n_perm_quick + 1)
}
all_perm_p[best_clf] <- perm_pvalue

clf_summary$perm_p <- all_perm_p[clf_summary$classifier]

write_csv(clf_summary, file.path(DAT_OUT, "panel_A_multi_classifier_auc.csv"))

perm_df <- tibble(
  best_classifier = best_clf,
  observed_auc    = best_auc,
  perm_pvalue     = perm_pvalue,
  null_auc_mean   = null_mean,
  null_auc_sd     = null_sd,
  null_auc_ci_lo  = null_mean - 1.96 * null_sd,
  null_auc_ci_hi  = null_mean + 1.96 * null_sd,
  n_permutations  = n_perm
)
write_csv(perm_df, file.path(DAT_OUT, "panel_A_permutation.csv"))


# Visualization: AUC bars + sparklines
PG_ROC_W <- 580
PG_ROC_H <- 553
txt_annot <- scale_text(BASE_STAT, PG_ROC_W)

mod_bio_df <- read_sheet_df(F06_SUPP, "WGCNA_mod_bio_labels")
mod_bio_lbl <- setNames(mod_bio_df$bio_label, mod_bio_df$module_color)
top2_features <- function(nm) {
  freq <- all_res[[nm]]$feature_freq
  top2 <- names(sort(freq, decreasing = TRUE))[seq_len(min(2, length(freq)))]
  bio <- mod_bio_lbl[gsub("^ME", "", top2)]
  bio[is.na(bio)] <- stringr::str_to_title(gsub(
    "^ME", "",
    top2[is.na(bio)]
  ))
  paste(bio, collapse = "\n")
}

clf_display <- c(
  Combined = "Pre + Post",
  Pre = "Pre",
  Post = "Post",
  DeltaME = "\u0394 Eigengene"
)

fmt_perm_p <- function(p) {
  if (is.na(p)) "" else if (p < 0.001) "p<0.001" else sprintf("p=%.3f", p)
}

bar_df <- clf_summary |>
  arrange(desc(auc)) |>
  mutate(
    classifier_name = clf_display[classifier],
    features = vapply(classifier, top2_features, character(1)),
    display = factor(classifier, levels = rev(classifier)),
    sig = case_when(
      is.na(perm_p) ~ "",
      perm_p < 0.05 ~ "*",
      perm_p < 0.10 ~ "\u2020",
      TRUE ~ ""
    ),
    auc_label = sprintf("%.2f%s", auc, sig),
    stats_line = sprintf(
      "CI %.2f\u2013%.2f \u00b7 %s",
      ci_lo, ci_hi, vapply(perm_p, fmt_perm_p, character(1))
    ),
    bar_text = paste(features, stats_line, sep = "\n")
  )

roc_all_df <- map_dfr(clf_names, function(nm) {
  probs <- all_res[[nm]]$probs
  roc_obj <- tryCatch(roc(true_labels, probs, quiet = TRUE),
    error = function(e) NULL
  )
  if (is.null(roc_obj)) {
    return(tibble())
  }
  tibble(
    classifier = nm,
    fpr = 1 - roc_obj$specificities,
    tpr = roc_obj$sensitivities
  )
})
write_csv(roc_all_df, file.path(DAT_OUT, "panel_A_roc_curves.csv"))

build_sparkline <- function(nm, auc_val = NA_real_, sig = "") {
  df <- roc_all_df |> filter(classifier == nm)
  col <- unname(clf_colors[nm])
  auc_lbl <- if (!is.na(auc_val)) sprintf("AUC %.2f%s", auc_val, sig) else ""
  ggplot(df, aes(x = fpr, y = tpr)) +
    geom_ribbon(aes(ymin = 0, ymax = tpr), fill = col, alpha = 0.35) +
    geom_abline(
      slope = 1, intercept = 0, linetype = "dashed",
      color = "grey70", linewidth = 0.25
    ) +
    geom_line(color = col, linewidth = 1.1) +
    annotate("text",
      x = 0.97, y = 0.06, label = auc_lbl,
      hjust = 1, vjust = 0, size = 8.5, fontface = "bold",
      color = "grey10"
    ) +
    scale_x_continuous(
      breaks = c(0, 1), labels = c("0", "1"),
      limits = c(0, 1), expand = c(0, 0)
    ) +
    scale_y_continuous(
      breaks = c(0, 1), labels = c("", "1"),
      limits = c(0, 1), expand = c(0, 0)
    ) +
    coord_fixed() +
    labs(x = "1\u2013Spec", y = "Sens") +
    theme_classic(base_size = 18) +
    theme(
      axis.title.x = element_text(
        size = 20, face = "bold",
        color = "grey15",
        margin = margin(t = -2)
      ),
      axis.title.y = element_text(
        size = 20, face = "bold",
        color = "grey15",
        margin = margin(r = -2)
      ),
      axis.text.x = element_text(
        size = 21, face = "bold",
        color = "grey25",
        margin = margin(t = 1)
      ),
      axis.text.y = element_text(
        size = 21, face = "bold",
        color = "grey25",
        margin = margin(r = 1)
      ),
      axis.ticks = element_line(color = "grey25", linewidth = 0.5),
      axis.ticks.length = unit(2, "pt"),
      axis.line = element_line(color = "grey25", linewidth = 0.6),
      plot.margin = margin(2, 2, 2, 2)
    )
}

disp_levels <- levels(bar_df$display)
bar_y <- setNames(seq_along(disp_levels), disp_levels)

spk_width <- 0.55
spk_half_h <- 0.48
spk_xmax <- 1.55

sparkline_layers <- lapply(seq_len(nrow(bar_df)), function(i) {
  nm <- as.character(bar_df$classifier[i])
  disp <- as.character(bar_df$display[i])
  y_center <- bar_y[[disp]]
  xmax_i <- spk_xmax
  xmin_i <- xmax_i - spk_width
  auc_i <- bar_df$auc[i]
  sig_i <- bar_df$sig[i]
  g <- ggplotGrob(build_sparkline(nm, auc_i, sig_i))
  annotation_custom(g,
    xmin = xmin_i, xmax = xmax_i,
    ymin = y_center - spk_half_h,
    ymax = y_center + spk_half_h
  )
})

y_label_map <- setNames(bar_df$classifier_name, as.character(bar_df$display))

pG_roc <- ggplot(bar_df, aes(x = auc, y = display)) +
  geom_col(aes(fill = classifier),
    color = "grey25",
    linewidth = 0.3, width = 0.6
  ) +
  geom_text(
    aes(
      x = 0.025,
      y = as.numeric(display),
      label = bar_text
    ),
    hjust = 0, vjust = 0.5,
    size = txt_annot * 1.55, fontface = "bold", color = "white",
    lineheight = 1.05
  ) +
  sparkline_layers +
  scale_fill_manual(values = clf_colors, guide = "none") +
  scale_y_discrete(labels = y_label_map, expand = expansion(add = 0.55)) +
  scale_x_continuous(
    limits = c(0, 1.55),
    breaks = c(0, 0.25, 0.5, 0.75, 1.0),
    expand = c(0, 0)
  ) +
  coord_cartesian(clip = "off") +
  labs(
    x = "AUC (LOOCV)", y = NULL,
    title = "Eigengenes Classify Age",
    subtitle = sprintf(
      "LOOCV logistic | n = %d | * p<0.05  \u2020 p<0.10",
      n_subj
    )
  ) +
  FIG_THEME +
  theme(
    panel.grid.major.y = element_blank(),
    axis.text.y = element_text(
      size = 26, face = "bold",
      color = "grey15",
      margin = margin(r = 1)
    ),
    axis.ticks.y = element_blank(),
    axis.text.x = element_text(
      size = 16, face = "bold",
      color = "grey25"
    ),
    axis.title.x = element_text(face = "bold", size = 18),
    plot.title = element_text(
      face = "bold",
      size = composite_text_sizes(640)$title,
      lineheight = 1.15,
      margin = margin(t = 10, b = 3, unit = "pt")
    ),
    plot.subtitle = element_text(
      size = composite_text_sizes(640)$subtitle,
      face = "bold.italic", color = "grey30",
      lineheight = 1.2,
      margin = margin(t = 2, b = 8, unit = "pt")
    ),
    plot.margin = margin(t = 8, r = 1, b = 1, l = 2)
  )

ggsave(file.path(RPT_PNG, "SUPP_F07_multivariate_classifier.png"), pG_roc,
  width = PG_ROC_W, height = PG_ROC_H, units = "mm", dpi = 300
)
ggsave(file.path(RPT_PDF, "SUPP_F07_multivariate_classifier.pdf"), pG_roc,
  width = PG_ROC_W, height = PG_ROC_H, units = "mm",
  device = get_pdf_device()
)

message("  SUPP_F07_multivariate_classifier saved")
