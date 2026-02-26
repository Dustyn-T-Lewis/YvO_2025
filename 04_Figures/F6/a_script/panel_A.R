################################################################################
#   Figure 6 — Panel A: Age Discrimination (PCA + ROC)
#   PCA of baseline module eigengenes + LOOCV logistic-regression ROC
################################################################################

source("04_Figures/F6/a_script/YvO_F6_setup.R")

# ---- PCA on baseline eigengenes ----
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
       color = "Age", fill = "Age") +
  THEME_PUB +
  theme(legend.position  = c(0.02, 0.98),
        legend.justification = c(0, 1),
        legend.background = element_rect(
          fill = scales::alpha("white", 0.85), color = NA))

# ---- ROC from LOOCV logistic regression ----
# Feature selection happens INSIDE each fold to prevent data leakage
# (Ambroise & McLachlan 2002, PNAS 99:6562-6566)
n_top_eigengenes <- min(5, ncol(me_pre))
loocv_probs <- numeric(length(common_subj))

for (i in seq_along(common_subj)) {
  train_me  <- me_pre[-i, , drop = FALSE]
  test_me   <- me_pre[i, , drop = FALSE]
  train_age <- subj_age$age[match(rownames(train_me), subj_age$subject_key)]
  train_y   <- ifelse(train_age == "Old", 1, 0)

  # Select top eigengenes from TRAINING data only (no leakage)
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
  annotate("text", x = 0.6, y = 0.2,
           label = sprintf("AUC = %.2f\n(95%% CI: %.2f\u2013%.2f)",
                           auc_val, ci_obj[1], ci_obj[3]),
           size = 2.5, fontface = "bold", color = "grey25") +
  labs(x = "1 \u2013 Specificity", y = "Sensitivity") +
  THEME_PUB +
  coord_equal()

# ---- Add panel title to PCA sub-panel ----
pA_pca <- pA_pca +
  labs(title    = "A  Age Discrimination via Module Eigengenes",
       subtitle = sprintf(
         "PCA of %d baseline eigengenes | nested LOOCV (top %d per fold)",
         ncol(MEs), n_top_eigengenes))

# ---- Save individual sub-panels ----
ggsave(file.path(RPT_DIR, "panel_A_pca.pdf"), pA_pca,
       width = 180, height = 180, units = "mm")
ggsave(file.path(RPT_DIR, "panel_A_pca.png"), pA_pca,
       width = 180, height = 180, units = "mm", dpi = 300)

ggsave(file.path(RPT_DIR, "panel_A_roc.pdf"), pA_roc,
       width = 180, height = 180, units = "mm")
ggsave(file.path(RPT_DIR, "panel_A_roc.png"), pA_roc,
       width = 180, height = 180, units = "mm", dpi = 300)

# ---- Save combined side-by-side ----
pA_combined <- pA_pca + pA_roc + plot_layout(widths = c(1, 1))

ggsave(file.path(RPT_DIR, "panel_A_age_discrimination.pdf"), pA_combined,
       width = 350, height = 180, units = "mm")
ggsave(file.path(RPT_DIR, "panel_A_age_discrimination.png"), pA_combined,
       width = 350, height = 180, units = "mm", dpi = 300)

# ---- Export data ----
write_csv(pca_scores, file.path(DAT_DIR, "fig6_panel_A_pca_scores.csv"))
write_csv(roc_df,     file.path(DAT_DIR, "fig6_panel_A_roc_curve.csv"))

# ============================================================================
# PHASE 5b: PERMUTATION AUC TEST (1000 permutations)
# Tests whether the observed AUC is significantly better than chance.
# Reference: Ojala & Garriga 2010, JMLR
# ============================================================================

cat("\n--- Permutation AUC Test ---\n")

# Wrap existing LOOCV in a reusable function
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
  if (i %% 100 == 0) cat(sprintf("  Permutation %d / %d\n", i, n_perm))
  shuffled_labels <- sample(true_labels)
  null_aucs[i] <- run_loocv_auc(shuffled_labels, me_pre)
}

perm_pvalue <- mean(null_aucs >= observed_auc)

cat(sprintf("  Observed AUC = %.3f, Permutation p = %.4f (%d permutations)\n",
            observed_auc, perm_pvalue, n_perm))

perm_df <- tibble(
  observed_auc = observed_auc,
  perm_pvalue  = perm_pvalue,
  null_auc_mean = mean(null_aucs),
  null_auc_sd   = sd(null_aucs),
  n_permutations = n_perm
)
write_csv(perm_df, file.path(DAT_DIR, "fig6_panel_A_permutation.csv"))

# Histogram of null AUCs
p_perm <- ggplot(tibble(auc = null_aucs), aes(x = auc)) +
  geom_histogram(bins = 40, fill = "grey70", color = "grey50", linewidth = 0.2) +
  geom_vline(xintercept = observed_auc, color = "#D6604D",
             linewidth = 1, linetype = "solid") +
  annotate("text", x = observed_auc, y = Inf, vjust = -0.5, hjust = -0.1,
           label = sprintf("Observed AUC = %.2f\nPermutation p = %.3f",
                           observed_auc, perm_pvalue),
           size = 3, fontface = "bold", color = "#D6604D") +
  labs(title = "Permutation Test for Age Classification AUC",
       subtitle = sprintf("Null distribution from %d label permutations", n_perm),
       x = "AUC", y = "Count") +
  THEME_PUB +
  theme(plot.title = element_text(size = 10, face = "bold"))

ggsave(file.path(RPT_DIR, "panel_A_permutation.pdf"), p_perm,
       width = 180, height = 120, units = "mm", device = pdf)

# ============================================================================
# PHASE 5c: PRE vs POST AUC COMPARISON
# Tests whether age classification is weakened after training (reversal
# evidence: if training narrows the age gap, Post AUC should be lower).
# Reference: DeLong et al. 1988, Biometrics; Robin et al. 2011, BMC Bioinformatics
# ============================================================================

cat("\n--- Pre vs Post AUC Comparison ---\n")

# Run LOOCV on Pre-only samples (already done above -> roc_obj)
auc_pre <- observed_auc
roc_pre <- roc_obj

# Run LOOCV on Post-only samples
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

# DeLong test comparing Pre vs Post ROC curves
# Note: these are unpaired (different observations), so use unpaired test
roc_compare <- roc.test(roc_pre, roc_post, method = "bootstrap",
                         paired = FALSE, boot.n = 10000)

cat(sprintf("  AUC Pre = %.3f, AUC Post = %.3f\n", auc_pre, auc_post))
cat(sprintf("  Bootstrap test p = %.4g\n", roc_compare$p.value))

pre_post_df <- tibble(
  auc_pre  = round(auc_pre, 3),
  auc_post = round(auc_post, 3),
  auc_diff = round(auc_pre - auc_post, 3),
  bootstrap_p = roc_compare$p.value,
  interpretation = ifelse(auc_post < auc_pre,
    "Post AUC lower: training narrows age gap (reversal evidence)",
    "Post AUC not lower: no evidence of age gap narrowing")
)
write_csv(pre_post_df, file.path(DAT_DIR, "fig6_panel_A_pre_vs_post.csv"))

# Overlay both ROC curves
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
           size = 2.5, fontface = "bold", color = "grey25") +
  labs(title = "Pre vs Post Age Classification",
       subtitle = "ROC curves from LOOCV logistic regression",
       x = "1 \u2013 Specificity", y = "Sensitivity") +
  THEME_PUB +
  coord_equal() +
  theme(legend.position = c(0.75, 0.25),
        legend.background = element_rect(fill = scales::alpha("white", 0.85), color = NA))

ggsave(file.path(RPT_DIR, "panel_A_roc_comparison.pdf"), p_roc_compare,
       width = 180, height = 180, units = "mm", device = pdf)

cat("  Saved: panel_A_permutation.pdf, panel_A_roc_comparison.pdf\n")
cat("Panel A done\n")
