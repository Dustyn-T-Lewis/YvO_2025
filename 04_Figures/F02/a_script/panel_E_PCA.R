# Figure 2 — Panel E: PCA Biplot + PERMANOVA
# Outputs: pE (ggplot object), panel_E_pca.pdf/.png

setwd(rprojroot::find_rstudio_root_file())
source("04_Figures/F02/a_script/style.R")

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(readr)
  library(ggplot2)
  library(vegan)
})

PE_W <- 145; PE_H <- 100

RPT_DIR <- "04_Figures/F02/b_reports"
DAT_DIR <- "04_Figures/F02/c_data"
dir.create(RPT_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(DAT_DIR, recursive = TRUE, showWarnings = FALSE)

imp_df <- read_csv("02_Imputation/c_data/01_imputed.csv",
                   show_col_types = FALSE)

ann_cols   <- c("uniprot_id", "protein", "gene", "description")
samp_names <- setdiff(names(imp_df), ann_cols)

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

imp_mat <- as.matrix(imp_df[, samp_names])
rownames(imp_mat) <- imp_df$gene

pdf_device <- get_pdf_device()

pca <- prcomp(t(imp_mat), center = TRUE, scale. = TRUE)
var_pct <- round(100 * summary(pca)$importance[2, 1:2], 1)

set.seed(42)
n_boot <- 1000
n_prot <- nrow(imp_mat)
boot_var <- matrix(NA_real_, nrow = n_boot, ncol = 2)
for (b in seq_len(n_boot)) {
  idx <- sample(n_prot, replace = TRUE)
  bp  <- prcomp(t(imp_mat[idx, ]), center = TRUE, scale. = TRUE)
  boot_var[b, ] <- 100 * summary(bp)$importance[2, 1:2]
}
var_ci <- data.frame(
  PC      = c("PC1", "PC2"),
  var_pct = var_pct,
  ci_lo   = apply(boot_var, 2, quantile, 0.025),
  ci_hi   = apply(boot_var, 2, quantile, 0.975)
)

pca_df <- as.data.frame(pca$x[, 1:2]) |>
  mutate(sample_id = rownames(pca$x)) |>
  left_join(meta, by = "sample_id")

# PERMANOVA — paired permutations within subjects
dist_mat <- dist(scale(t(imp_mat)))
set.seed(42)
perm_res <- adonis2(dist_mat ~ age * time, data = meta,
                    permutations = how(nperm = 999, blocks = meta$subject),
                    by = "terms")

perm_terms <- c("age", "time", "age:time")
perm_r2 <- perm_res[perm_terms, "R2"]
perm_pv <- perm_res[perm_terms, "Pr(>F)"]
perm_label <- sprintf(
  "PERMANOVA\nAge  R\u00b2 = %.3f,  %s\nTime  R\u00b2 = %.3f,  %s\nAge\u00d7Time  R\u00b2 = %.3f,  %s",
  perm_r2[1], fmt_p(perm_pv[1]),
  perm_r2[2], fmt_p(perm_pv[2]),
  perm_r2[3], fmt_p(perm_pv[3]))

bd_age  <- betadisper(dist_mat, meta$age)
bd_time <- betadisper(dist_mat, meta$time)
bd_grp  <- betadisper(dist_mat, meta$group)
bd_age_p  <- permutest(bd_age,  pairwise = FALSE, permutations = 999)$tab$`Pr(>F)`[1]
bd_time_p <- permutest(bd_time, pairwise = FALSE, permutations = 999)$tab$`Pr(>F)`[1]
bd_grp_p  <- permutest(bd_grp,  pairwise = FALSE, permutations = 999)$tab$`Pr(>F)`[1]
if (bd_age_p < 0.05 || bd_time_p < 0.05)
  warning("Heterogeneous dispersions detected — interpret PERMANOVA with caution")

pE <- ggplot(pca_df, aes(x = PC1, y = PC2, color = group, shape = group)) +
  stat_ellipse(aes(fill = group), geom = "polygon",
               alpha = 0.10, level = 0.80, show.legend = FALSE) +
  stat_ellipse(aes(group = group), level = 0.80, linewidth = 0.4,
               linetype = "dashed", show.legend = FALSE) +
  geom_point(size = 2.5, alpha = 0.85) +
  annotate("text", x = Inf, y = Inf, label = perm_label,
           hjust = 1.05, vjust = 1.15,
           size = scale_text(BASE_STAT - 1.2, PE_W), color = "grey30",
           fontface = "bold") +
  scale_color_manual(values = PCA_COLORS,
                     labels = c("Young Pre", "Young Post", "Old Pre", "Old Post"),
                     guide = guide_legend(override.aes = list(size = 2.5))) +
  scale_fill_manual(values = PCA_COLORS, guide = "none") +
  scale_shape_manual(values = PCA_SHAPES,
                     labels = c("Young Pre", "Young Post", "Old Pre", "Old Post")) +
  labs(title = "Principal Component Analysis (PCA)",
       subtitle = sprintf("PERMANOVA R\u00b2 = %.2f, %s | %s proteins, %d samples",
                          sum(perm_r2), fmt_p(min(perm_pv)),
                          format(nrow(imp_df), big.mark = ","), nrow(meta)),
       x = sprintf("PC1 (%.1f%% [%.1f, %.1f])", var_pct[1], var_ci$ci_lo[1], var_ci$ci_hi[1]),
       y = sprintf("PC2 (%.1f%% [%.1f, %.1f])", var_pct[2], var_ci$ci_lo[2], var_ci$ci_hi[2]),
       tag = "E") +
  FIG_THEME + theme(legend.position = c(0.88, 0.12),
                    legend.background = element_rect(fill = alpha("white", 0.8),
                                                     color = NA),
                    legend.title = element_blank(),
                    legend.text = element_text(size = FIG_LEGEND_TEXT),
                    legend.key.size = unit(3, "mm"),
                    legend.spacing.y = unit(0.5, "mm"))

write.csv(var_ci, file.path(DAT_DIR, "audit_panel_E_pca_variance_ci.csv"),
          row.names = FALSE)
betadisper_results <- data.frame(
  factor      = c("age", "time", "group"),
  p_value     = c(bd_age_p, bd_time_p, bd_grp_p),
  significant = c(bd_age_p < 0.05, bd_time_p < 0.05, bd_grp_p < 0.05)
)
write.csv(betadisper_results, file.path(DAT_DIR, "audit_panel_E_betadisper.csv"),
          row.names = FALSE)

ggsave(file.path(RPT_DIR, "panel_E_pca.pdf"), pE,
       width = PE_W, height = PE_H, units = "mm", device = pdf_device)
ggsave(file.path(RPT_DIR, "panel_E_pca.png"), pE,
       width = PE_W, height = PE_H, units = "mm", dpi = 300)
