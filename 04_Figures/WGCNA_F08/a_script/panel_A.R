# Figure 6 — Panel A: Age Discrimination (PCA + ROC)
# PCA of baseline module eigengenes + LOOCV logistic-regression ROC
#
# STAT AUDIT (2026-02-27):
# 1. DATA LEAKAGE CHECK: PASS. Feature selection INSIDE each LOOCV fold
#    (Ambroise & McLachlan 2002, PNAS 99:6562).
# 2. AUC 95% CI via pROC::ci.auc() (DeLong 1988). PASS.
# 3. Permutation test: 1000 permutations + 95% CI on null distribution.
# 4. Pre vs Post AUC: bootstrap roc.test() (10k resamples) + CI on diff.
# 5. Small-sample: N ~15-31; LOOCV appropriate.

setwd(rprojroot::find_rstudio_root_file())
source("04_Figures/shared/style.R")

suppressPackageStartupMessages({
  library(tidyverse)
  library(patchwork)
  library(pROC)
})

RPT <- "04_Figures/WGCNA_F08/b_reports"
DAT <- "04_Figures/WGCNA_F08/c_data"

MEs       <- readRDS(file.path(DAT, "MEs.rds"))
me_pre    <- readRDS(file.path(DAT, "me_pre.rds"))
me_post   <- readRDS(file.path(DAT, "me_post.rds"))
subj_age  <- read_csv(file.path(DAT, "subj_age.csv"), show_col_types = FALSE)
shared    <- readRDS(file.path(DAT, "shared_objects.rds"))
common_subj <- shared$common_subj

pdf_device <- get_pdf_device()

message("Panel A: age discrimination (PCA + ROC)...")

PA_PCA_W <- 180   # PCA sub-panel width mm
PA_PCA_H <- 180
PA_ROC_W <- 180   # ROC sub-panel width mm
PA_ROC_H <- 180
PA_COMB_W <- 350  # combined width mm
PA_COMB_H <- 180
PA_PERM_W <- 180  # permutation histogram width mm
PA_PERM_H <- 120

txt_annot <- scale_text(BASE_STAT, PA_PCA_W)
txt_axis  <- scale_text(BASE_STAT, PA_PCA_W)
txt_perm  <- scale_text(BASE_STAT, PA_PERM_W)

me_pre_df <- as.data.frame(me_pre) %>%
  rownames_to_column("subject_key") %>%
  left_join(subj_age, by = "subject_key")

pca_me <- prcomp(me_pre, scale. = TRUE, center = TRUE)
pca_scores <- as.data.frame(pca_me$x[, 1:2]) %>%
  rownames_to_column("subject_key") %>%
  left_join(subj_age, by = "subject_key")
var_expl <- summary(pca_me)$importance[2, 1:2] * 100

pA_pca <- ggplot(pca_scores, aes(x = PC1, y = PC2, color = age)) +
  stat_ellipse(level = 0.80, alpha = 0.10, geom = "polygon",
               aes(fill = age), linewidth = 0.4, linetype = "dashed") +
  geom_point(size = 2.0, alpha = 0.85) +
  scale_color_manual(values = AGE_COLORS) +
  scale_fill_manual(values = AGE_COLORS) +
  labs(x = sprintf("PC1 (%.1f%%)", var_expl[1]),
       y = sprintf("PC2 (%.1f%%)", var_expl[2]),
       title = "Age Discrimination via Module Eigengenes",
       subtitle = sprintf(
         "PCA of %d baseline eigengenes | nested LOOCV (top %d per fold)",
         ncol(MEs), min(5, ncol(me_pre))),
       color = "Age", fill = "Age") +
  FIG_THEME +
  theme(legend.position = "none")

# Feature selection INSIDE each fold (Ambroise & McLachlan 2002)
n_top_eigengenes <- min(5, ncol(me_pre))
loocv_probs <- numeric(length(common_subj))

for (i in seq_along(common_subj)) {
  train_me  <- me_pre[-i, , drop = FALSE]
  test_me   <- me_pre[i, , drop = FALSE]
  train_age <- subj_age$age[match(rownames(train_me), subj_age$subject_key)]
  train_y   <- ifelse(train_age == "Old", 1, 0)

  train_cors <- abs(cor(train_me, train_y))
  top_mes_i  <- names(sort(train_cors[, 1], decreasing = TRUE))[
    seq_len(n_top_eigengenes)]

  fit_data <- cbind(y = train_y, as.data.frame(train_me[, top_mes_i, drop = FALSE]))
  fit <- tryCatch(
    glm(y ~ ., data = fit_data, family = binomial(link = "logit")),
    warning = function(w) suppressWarnings(
      glm(y ~ ., data = fit_data, family = binomial(link = "logit"))
    )
  )
  loocv_probs[i] <- predict(fit, newdata = as.data.frame(test_me[, top_mes_i, drop = FALSE]),
                              type = "response")
}

true_labels <- ifelse(
  subj_age$age[match(common_subj, subj_age$subject_key)] == "Old", 1, 0
)
roc_obj <- roc(true_labels, loocv_probs, quiet = TRUE)
auc_val <- auc(roc_obj)
ci_obj  <- ci.auc(roc_obj)

roc_df <- data.frame(fpr = 1 - roc_obj$specificities,
                     tpr = roc_obj$sensitivities)

pA_roc <- ggplot(roc_df, aes(x = fpr, y = tpr)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed",
              color = "grey40", linewidth = 0.3) +
  geom_line(linewidth = 0.8, color = "#D6604D") +
  labs(x = "1 \u2013 Specificity", y = "Sensitivity") +
  FIG_THEME +
  coord_equal() +
  theme(legend.position = "none")

write_csv(pca_scores, file.path(DAT, "01_panel_A_pca_scores.csv"))
write_csv(roc_df,     file.path(DAT, "01_panel_A_roc_curve.csv"))

# Permutation AUC test (Ojala & Garriga 2010, JMLR)
run_loocv_auc <- function(labels, me_data, n_top = n_top_eigengenes) {
  n <- length(labels)
  probs <- numeric(n)
  for (i in seq_len(n)) {
    train_me <- me_data[-i, , drop = FALSE]
    test_me  <- me_data[i, , drop = FALSE]
    train_y  <- labels[-i]

    train_cors <- abs(cor(train_me, train_y))
    top_mes_i  <- names(sort(train_cors[, 1], decreasing = TRUE))[seq_len(n_top)]

    fit_data <- cbind(y = train_y, as.data.frame(train_me[, top_mes_i, drop = FALSE]))
    fit <- tryCatch(
      glm(y ~ ., data = fit_data, family = binomial(link = "logit")),
      warning = function(w) suppressWarnings(
        glm(y ~ ., data = fit_data, family = binomial(link = "logit")))
    )
    probs[i] <- predict(fit, newdata = as.data.frame(test_me[, top_mes_i, drop = FALSE]),
                          type = "response")
  }
  tryCatch(as.numeric(auc(roc(labels, probs, quiet = TRUE))), error = function(e) 0.5)
}

observed_auc <- as.numeric(auc_val)

set.seed(42)
n_perm <- 1000
null_aucs <- numeric(n_perm)

for (i in seq_len(n_perm)) {
  shuffled_labels <- sample(true_labels)
  null_aucs[i] <- run_loocv_auc(shuffled_labels, me_pre)
}

perm_pvalue <- mean(null_aucs >= observed_auc)

null_ci_lo <- mean(null_aucs) - 1.96 * sd(null_aucs)
null_ci_hi <- mean(null_aucs) + 1.96 * sd(null_aucs)

message(sprintf("  AUC=%.3f, perm p=%.4f, null mean=%.3f [%.3f, %.3f]",
                observed_auc, perm_pvalue, mean(null_aucs), null_ci_lo, null_ci_hi))

perm_label <- if (perm_pvalue < 0.001) "< 0.001" else sprintf("= %.3f", perm_pvalue)
pA_roc <- pA_roc +
  annotate("text", x = 0.6, y = 0.2,
           label = sprintf("AUC = %.2f\n(95%% CI: %.2f\u2013%.2f)\nPermutation p %s\nN = %d",
                           auc_val, ci_obj[1], ci_obj[3],
                           perm_label, length(common_subj)),
           size = txt_annot, fontface = "bold", color = "grey25")

pA_combined <- pA_pca + pA_roc + plot_layout(widths = c(1, 1))
ggsave(file.path(RPT, "panel_A_pca.pdf"), pA_pca,
       width = PA_PCA_W, height = PA_PCA_H, units = "mm",
       device = pdf_device, limitsize = FALSE)
ggsave(file.path(RPT, "panel_A_roc.pdf"), pA_roc,
       width = PA_ROC_W, height = PA_ROC_H, units = "mm",
       device = pdf_device, limitsize = FALSE)
ggsave(file.path(RPT, "panel_A_age_discrimination.pdf"), pA_combined,
       width = PA_COMB_W, height = PA_COMB_H, units = "mm",
       device = pdf_device, limitsize = FALSE)
ggsave(file.path(RPT, "panel_A_age_discrimination.png"), pA_combined,
       width = PA_COMB_W, height = PA_COMB_H, units = "mm", dpi = 300)

perm_df <- tibble(
  observed_auc    = observed_auc,
  perm_pvalue     = perm_pvalue,
  null_auc_mean   = mean(null_aucs),
  null_auc_sd     = sd(null_aucs),
  null_auc_ci_lo  = null_ci_lo,
  null_auc_ci_hi  = null_ci_hi,
  n_permutations  = n_perm
)
write_csv(perm_df, file.path(DAT, "01_panel_A_permutation.csv"))

p_perm <- ggplot(tibble(auc = null_aucs), aes(x = auc)) +
  geom_histogram(bins = 40, fill = "grey70", color = "grey50", linewidth = 0.2) +
  geom_vline(xintercept = observed_auc, color = "#D6604D",
             linewidth = 1, linetype = "solid") +
  annotate("text", x = observed_auc, y = Inf, vjust = -0.5, hjust = -0.1,
           label = sprintf("Observed AUC = %.2f\nPermutation p = %.3f",
                           observed_auc, perm_pvalue),
           size = txt_perm, fontface = "bold", color = "#D6604D") +
  labs(title = "Permutation Test for Age Classification AUC",
       subtitle = sprintf("Null distribution from %d label permutations", n_perm),
       x = "AUC", y = "Count") +
  FIG_THEME +
  theme(legend.position = "none")

ggsave(file.path(RPT, "panel_A_permutation.pdf"), p_perm,
       width = PA_PERM_W, height = PA_PERM_H, units = "mm",
       device = pdf_device, limitsize = FALSE)
ggsave(file.path(RPT, "panel_A_permutation.png"), p_perm,
       width = PA_PERM_W, height = PA_PERM_H, units = "mm", dpi = 300)

# Pre vs Post AUC comparison (DeLong 1988; Robin 2011)
auc_pre <- observed_auc
roc_pre <- roc_obj

me_post_df <- as.data.frame(me_post) %>%
  rownames_to_column("subject_key") %>%
  left_join(subj_age, by = "subject_key")

post_labels <- ifelse(
  subj_age$age[match(common_subj, subj_age$subject_key)] == "Old", 1, 0)

loocv_probs_post <- numeric(length(common_subj))
for (i in seq_along(common_subj)) {
  train_me  <- me_post[-i, , drop = FALSE]
  test_me   <- me_post[i, , drop = FALSE]
  train_age <- subj_age$age[match(rownames(train_me), subj_age$subject_key)]
  train_y   <- ifelse(train_age == "Old", 1, 0)

  train_cors <- abs(cor(train_me, train_y))
  top_mes_i  <- names(sort(train_cors[, 1], decreasing = TRUE))[
    seq_len(n_top_eigengenes)]

  fit_data <- cbind(y = train_y, as.data.frame(train_me[, top_mes_i, drop = FALSE]))
  fit <- tryCatch(
    glm(y ~ ., data = fit_data, family = binomial(link = "logit")),
    warning = function(w) suppressWarnings(
      glm(y ~ ., data = fit_data, family = binomial(link = "logit")))
  )
  loocv_probs_post[i] <- predict(fit, newdata = as.data.frame(test_me[, top_mes_i, drop = FALSE]),
                                   type = "response")
}

roc_post <- roc(post_labels, loocv_probs_post, quiet = TRUE)
auc_post <- as.numeric(auc(roc_post))
ci_post  <- ci.auc(roc_post)

# Unpaired bootstrap test (different observations per timepoint)
roc_compare <- roc.test(roc_pre, roc_post, method = "bootstrap",
                         paired = FALSE, boot.n = 10000)

set.seed(42)
n_boot_diff <- 2000
boot_diffs  <- numeric(n_boot_diff)
for (b in seq_len(n_boot_diff)) {
  idx_pre  <- sample(seq_along(true_labels), replace = TRUE)
  idx_post <- sample(seq_along(post_labels), replace = TRUE)
  roc_b_pre  <- tryCatch(as.numeric(auc(roc(true_labels[idx_pre],
    loocv_probs[idx_pre], quiet = TRUE))), error = function(e) NA_real_)
  roc_b_post <- tryCatch(as.numeric(auc(roc(post_labels[idx_post],
    loocv_probs_post[idx_post], quiet = TRUE))), error = function(e) NA_real_)
  boot_diffs[b] <- roc_b_pre - roc_b_post
}
auc_diff_ci <- quantile(boot_diffs, c(0.025, 0.975), na.rm = TRUE)

message(sprintf("  AUC Pre=%.3f, Post=%.3f, diff=%.3f [%.3f, %.3f], p=%.4g",
                auc_pre, auc_post, auc_pre - auc_post,
                auc_diff_ci[1], auc_diff_ci[2], roc_compare$p.value))

pre_post_df <- tibble(
  auc_pre  = round(auc_pre, 3),
  auc_post = round(auc_post, 3),
  auc_diff = round(auc_pre - auc_post, 3),
  auc_diff_ci_lo = round(auc_diff_ci[1], 3),
  auc_diff_ci_hi = round(auc_diff_ci[2], 3),
  bootstrap_p = roc_compare$p.value,
  interpretation = ifelse(auc_post < auc_pre,
    "Post AUC lower: training narrows age gap (reversal evidence)",
    "Post AUC not lower: no evidence of age gap narrowing")
)
write_csv(pre_post_df, file.path(DAT, "01_panel_A_pre_vs_post.csv"))

roc_post_df <- data.frame(fpr = 1 - roc_post$specificities,
                           tpr = roc_post$sensitivities)

p_roc_compare <- ggplot() +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed",
              color = "grey40", linewidth = 0.3) +
  geom_line(data = roc_df, aes(x = fpr, y = tpr, color = "Pre"),
            linewidth = 0.8) +
  geom_line(data = roc_post_df, aes(x = fpr, y = tpr, color = "Post"),
            linewidth = 0.8, linetype = "dashed") +
  scale_color_manual(values = c(Pre = "#D6604D", Post = "#5DA5DA"),
                      name = "Timepoint") +
  annotate("text", x = 0.55, y = 0.15,
           label = sprintf("AUC Pre = %.2f\nAUC Post = %.2f\nBootstrap p = %.3g",
                           auc_pre, auc_post, roc_compare$p.value),
           size = txt_annot, fontface = "bold", color = "grey25") +
  labs(title = "Pre vs Post Age Classification",
       subtitle = "ROC curves from LOOCV logistic regression",
       x = "1 \u2013 Specificity", y = "Sensitivity") +
  FIG_THEME +
  coord_equal() +
  theme(legend.position = "none")

ggsave(file.path(RPT, "panel_A_roc_comparison.pdf"), p_roc_compare,
       width = PA_ROC_W, height = PA_ROC_H, units = "mm",
       device = pdf_device, limitsize = FALSE)
ggsave(file.path(RPT, "panel_A_roc_comparison.png"), p_roc_compare,
       width = PA_ROC_W, height = PA_ROC_H, units = "mm", dpi = 300)

message("  Panel A saved")
