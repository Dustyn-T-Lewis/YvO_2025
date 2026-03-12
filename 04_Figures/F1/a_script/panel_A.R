# Figure 1 — Panel A: CV% Violins (Inter-Individual Variability)
# Faceted by Age (Young | Old), Pre/Post x-axis. Median labels with bootstrap CIs.
# Outputs: pA (ggplot object), panel_A_cv.pdf/.png

setwd(rprojroot::find_rstudio_root_file())
source("04_Figures/shared/style.R")

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(readr)
  library(ggplot2)
  library(ggbeeswarm)
})

PA_W <- 110; PA_H <- 120

RPT_DIR <- "04_Figures/F1/b_reports"
DAT_DIR <- "04_Figures/F1/c_data"
dir.create(RPT_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(DAT_DIR, recursive = TRUE, showWarnings = FALSE)

norm_df <- read_csv("01_normalization/c_data/02_normalized.csv",
                    show_col_types = FALSE)

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

pdf_device <- get_pdf_device()

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

cv_df$age  <- factor(ifelse(grepl("Young", cv_df$group), "Young", "Old"),
                     levels = c("Young", "Old"))
cv_df$time <- factor(ifelse(grepl("Pre", cv_df$group), "Pre", "Post"),
                     levels = c("Pre", "Post"))

cv_med <- cv_df |> group_by(age, time, group) |>
  summarise(med = median(cv), .groups = "drop")

# Bootstrap 95% CI on median CV per group
set.seed(42)
boot_median_ci <- function(x, R = 2000, conf = 0.95) {
  meds <- replicate(R, median(sample(x, replace = TRUE)))
  qs   <- quantile(meds, c((1 - conf) / 2, (1 + conf) / 2))
  c(lower = unname(qs[1]), upper = unname(qs[2]))
}

# Pairwise Wilcoxon tests, BH corrected (audit only — not on figure)
bracket_comps <- list(c("Young_Pre", "Young_Post"), c("Old_Pre", "Old_Post"),
                      c("Young_Pre", "Old_Pre"),    c("Young_Post", "Old_Post"))
bracket_pvals_raw <- sapply(bracket_comps, function(pair)
  wilcox.test(cv_df$cv[cv_df$group == pair[1]],
              cv_df$cv[cv_df$group == pair[2]])$p.value)
bracket_pvals <- p.adjust(bracket_pvals_raw, method = "BH")

cliffs_delta <- function(x, y) {
  nx <- length(x); ny <- length(y)
  d <- outer(x, y, function(a, b) sign(a - b))
  sum(d) / (nx * ny)
}
cliff_results <- data.frame(
  comparison = sapply(bracket_comps, paste, collapse = " vs "),
  p_raw      = bracket_pvals_raw,
  p_bh       = bracket_pvals,
  cliffs_d   = sapply(bracket_comps, function(pair)
    cliffs_delta(cv_df$cv[cv_df$group == pair[1]],
                 cv_df$cv[cv_df$group == pair[2]]))
)

y_cap <- ceiling(max(cv_df$cv, na.rm = TRUE) / 10) * 10

cv_ci <- cv_df |>
  group_by(age, time, group) |>
  summarise(
    med   = median(cv),
    ci_lo = boot_median_ci(cv)[["lower"]],
    ci_hi = boot_median_ci(cv)[["upper"]],
    .groups = "drop"
  )

pA <- ggplot(cv_df, aes(x = time, y = cv, fill = group)) +
  geom_violin(alpha = 0.5, linewidth = 0.3, color = "black", scale = "width") +
  geom_quasirandom(aes(color = group), alpha = 0.15, size = 0.3,
                   width = 0.25, groupOnX = TRUE, show.legend = FALSE) +
  geom_boxplot(width = 0.15, outlier.shape = NA, linewidth = 0.3,
               color = "black", fill = "white", coef = 0) +
  geom_hline(yintercept = 25, linetype = "dashed", color = "grey50",
             linewidth = 0.4) +
  # Median labels with bootstrap CI below violins
  geom_label(data = cv_ci,
             aes(x = time, y = -2,
                 label = sprintf("%.0f%% [%.0f\u2013%.0f]", med, ci_lo, ci_hi)),
             size = scale_text(BASE_COUNT - 0.5, PA_W / 2),
             fontface = "bold", fill = alpha("white", 0.8),
             linewidth = 0.2, label.padding = unit(1.2, "pt")) +
  facet_wrap(~ age, nrow = 1) +
  scale_fill_manual(values = GROUP_FILL) +
  scale_color_manual(values = GROUP_FILL) +
  coord_cartesian(ylim = c(-5, y_cap + 5)) +
  labs(title = "Inter-Individual Variability (CV%)",
       subtitle = "Protein-level CV% (cycloess-normalized)",
       x = NULL, y = "CV (%)",
       tag = "A") +
  FIG_THEME + theme(legend.position = "none",
                    panel.spacing = unit(8, "mm"))

write.csv(as.data.frame(cv_ci),
          file.path(DAT_DIR, "audit_panel_A_median_cv_ci.csv"), row.names = FALSE)
write.csv(cliff_results,
          file.path(DAT_DIR, "audit_panel_A_wilcoxon_effects.csv"), row.names = FALSE)

ggsave(file.path(RPT_DIR, "panel_A_cv.pdf"), pA,
       width = PA_W, height = PA_H, units = "mm", device = pdf_device)
ggsave(file.path(RPT_DIR, "panel_A_cv.png"), pA,
       width = PA_W, height = PA_H, units = "mm", dpi = 300)
