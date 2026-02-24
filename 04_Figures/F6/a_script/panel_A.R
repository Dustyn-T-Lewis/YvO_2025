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

cat("Panel A done\n")
