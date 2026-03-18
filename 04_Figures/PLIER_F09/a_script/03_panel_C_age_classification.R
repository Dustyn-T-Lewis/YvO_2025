# Panel C — Age classification from baseline pathway activity.
# Q2: Can baseline LV scores classify Young vs Old?
# Method: PCA + LOOCV logistic regression + ROC/AUC.

setwd(rprojroot::find_rstudio_root_file())
source("04_Figures/shared/style.R")

library(readr)
library(dplyr)
library(ggplot2)
library(pROC)

DAT <- "04_Figures/PLIER_F09/c_data"
RPT <- "04_Figures/PLIER_F09/b_results"

# --- Load data ---
lv_scores <- read_csv(file.path(DAT, "02_lv_scores.csv"), show_col_types = FALSE)
dal <- readRDS("02_Imputation/c_data/01_DAList_imputed.rds")
meta <- as.data.frame(dal$metadata)

lv_meta <- lv_scores |>
  inner_join(meta, by = c("sample_id" = "Col_ID"))

pre <- lv_meta |> filter(Timepoint == "Pre")
lv_cols <- grep("^LV\\d+$", names(pre), value = TRUE)
lv_mat <- as.matrix(pre[, lv_cols])

n_young <- sum(pre$Group == "Young")
n_old <- sum(pre$Group == "Old")
cat(sprintf("Classification: Young n=%d, Old n=%d\n", n_young, n_old))

# --- PCA on baseline LV scores ---
pca_res <- prcomp(lv_mat, center = TRUE, scale. = TRUE)
pca_df <- data.frame(
  PC1 = pca_res$x[, 1],
  PC2 = pca_res$x[, 2],
  Group = pre$Group,
  Group_Time = pre$Group_Time
)
var_exp <- summary(pca_res)$importance[2, 1:2] * 100

pC_pca <- ggplot(pca_df, aes(x = PC1, y = PC2, color = Group, shape = Group)) +
  stat_ellipse(aes(fill = Group), geom = "polygon", alpha = 0.10,
               level = 0.80, linewidth = 0.4, linetype = "dashed") +
  geom_point(size = 2.5, alpha = 0.85) +
  scale_color_manual(values = AGE_COLORS) +
  scale_fill_manual(values = AGE_COLORS) +
  labs(title = "C  PCA of Baseline Pathway Activity",
       subtitle = "LV scores, Pre-training only",
       x = sprintf("PC1 (%.1f%%)", var_exp[1]),
       y = sprintf("PC2 (%.1f%%)", var_exp[2])) +
  FIG_THEME +
  theme(legend.position = c(0.85, 0.15),
        legend.background = element_rect(fill = alpha("white", 0.85), color = NA))

ggsave(file.path(RPT, "13_panel_C_pca.pdf"), pC_pca,
       width = 140, height = 130, units = "mm")
ggsave(file.path(RPT, "13_panel_C_pca.png"), pC_pca,
       width = 140, height = 130, units = "mm", dpi = 300)

# --- LOOCV logistic regression ---
n <- nrow(pre)
pred_probs <- numeric(n)
true_labels <- as.integer(pre$Group == "Old")

for (i in seq_len(n)) {
  train_df <- data.frame(y = true_labels[-i], lv_mat[-i, , drop = FALSE])
  test_df <- data.frame(lv_mat[i, , drop = FALSE])

  # Use top PCs to avoid overfitting with many LVs
  pca_train <- prcomp(as.matrix(train_df[, -1]), center = TRUE, scale. = TRUE)
  n_pc <- min(5, ncol(pca_train$x))  # cap at 5 PCs
  train_pcs <- data.frame(y = train_df$y, pca_train$x[, 1:n_pc])
  test_scaled <- scale(as.matrix(test_df), center = pca_train$center,
                       scale = pca_train$scale)
  test_pcs <- data.frame(test_scaled %*% pca_train$rotation[, 1:n_pc])

  fit <- glm(y ~ ., data = train_pcs, family = binomial)
  pred_probs[i] <- predict(fit, newdata = test_pcs, type = "response")
}

# --- ROC/AUC ---
roc_obj <- roc(true_labels, pred_probs, quiet = TRUE)
auc_val <- as.numeric(auc(roc_obj))
ci_auc <- ci.auc(roc_obj)

roc_df <- data.frame(
  specificity = 1 - roc_obj$specificities,
  sensitivity = roc_obj$sensitivities
)

pC_roc <- ggplot(roc_df, aes(x = specificity, y = sensitivity)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey60") +
  geom_line(color = AGE_COLORS["Old"], linewidth = 1) +
  annotate("text", x = 0.55, y = 0.15,
           label = sprintf("AUC = %.2f [%.2f, %.2f]\nWGCNA AUC = 0.97",
                           auc_val, ci_auc[1], ci_auc[3]),
           size = 3.5, fontface = "bold", color = "grey25", hjust = 0) +
  labs(title = "C  LOOCV Age Classification",
       subtitle = sprintf("Logistic regression on LV-derived PCs (n=%d)", n),
       x = "1 - Specificity", y = "Sensitivity") +
  FIG_THEME

ggsave(file.path(RPT, "13_panel_C_roc.pdf"), pC_roc,
       width = 130, height = 120, units = "mm")
ggsave(file.path(RPT, "13_panel_C_roc.png"), pC_roc,
       width = 130, height = 120, units = "mm", dpi = 300)

# --- Export ---
class_summary <- tibble(
  method = "LOOCV logistic regression (top 5 PCs of LV scores)",
  n_young = n_young,
  n_old = n_old,
  AUC = auc_val,
  AUC_95CI_lo = as.numeric(ci_auc[1]),
  AUC_95CI_hi = as.numeric(ci_auc[3]),
  WGCNA_AUC_reference = 0.97
)

pred_df <- tibble(
  sample_id = pre$sample_id,
  true_group = pre$Group,
  true_label = true_labels,
  predicted_prob_old = pred_probs,
  predicted_class = ifelse(pred_probs > 0.5, "Old", "Young")
)

write_csv(bind_rows(class_summary, tibble()), file.path(DAT, "13_panel_C_age_classification.csv"))
write_csv(pred_df, file.path(DAT, "13_panel_C_predictions.csv"))

cat(sprintf("Panel C: AUC = %.3f [%.3f, %.3f] (WGCNA ref: 0.97)\n",
            auc_val, ci_auc[1], ci_auc[3]))
