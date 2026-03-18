# Figure 2 — Panel A (ICC): Intraclass Correlation Coefficient
# ICC(3,1) per protein using Pre/Post as two measurements per subject.
# High ICC = stable subject trait; low ICC = responsive to training.
# Outputs: pA_ICC (ggplot object), panel_A_ICC.pdf/.png

setwd(rprojroot::find_rstudio_root_file())
source("04_Figures/F02/a_script/style.R")

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(readr)
  library(ggplot2)
  library(psych)
})

PA_ICC_W <- 110; PA_ICC_H <- 120

RPT_DIR <- "04_Figures/F02/b_reports"
DAT_DIR <- "04_Figures/F02/c_data"
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
    subject  = paste0(prefix, "_", subj_num)
  )
meta$age  <- factor(meta$age, levels = c("Young", "Old"))
meta$time <- factor(meta$time, levels = c("Pre", "Post"))

log_mat <- as.matrix(norm_df[, samp_names])
rownames(log_mat) <- norm_df$gene

pdf_device <- get_pdf_device()

# Compute ICC(3,1) per protein within each age group
# ICC type 3 = two-way mixed, single measures, consistency
compute_icc_per_age <- function(age_group) {
  age_meta <- meta |> filter(age == age_group)
  subjects <- unique(age_meta$subject)

  # Build n_subjects x 2 matrix (Pre, Post) per protein
  pre_ids  <- age_meta$sample_id[age_meta$time == "Pre"]
  post_ids <- age_meta$sample_id[age_meta$time == "Post"]

  # Match subjects for paired structure
  pre_subj  <- age_meta$subject[age_meta$time == "Pre"]
  post_subj <- age_meta$subject[age_meta$time == "Post"]
  common    <- intersect(pre_subj, post_subj)

  pre_idx  <- match(common, pre_subj)
  post_idx <- match(common, post_subj)

  pre_mat  <- log_mat[, pre_ids[pre_idx], drop = FALSE]
  post_mat <- log_mat[, post_ids[post_idx], drop = FALSE]

  n_prot <- nrow(log_mat)
  icc_vals <- numeric(n_prot)

  for (i in seq_len(n_prot)) {
    pre_vals  <- pre_mat[i, ]
    post_vals <- post_mat[i, ]

    # Skip proteins with too many NAs
    ok <- !is.na(pre_vals) & !is.na(post_vals)
    if (sum(ok) < 3) { icc_vals[i] <- NA_real_; next }

    rating_mat <- cbind(pre_vals[ok], post_vals[ok])
    icc_res <- tryCatch(
      psych::ICC(rating_mat, missing = FALSE, lmer = FALSE),
      error = function(e) NULL
    )
    if (is.null(icc_res)) { icc_vals[i] <- NA_real_; next }
    # ICC3 = two-way mixed, consistency, single measures
    icc_vals[i] <- icc_res$results$ICC[3]
  }

  tibble(gene = norm_df$gene, age = age_group, icc = icc_vals)
}

icc_young <- compute_icc_per_age("Young")
icc_old   <- compute_icc_per_age("Old")
icc_df    <- bind_rows(icc_young, icc_old) |> filter(!is.na(icc))
icc_df$age <- factor(icc_df$age, levels = c("Young", "Old"))

# Bootstrap 95% CI on median ICC per group
set.seed(42)
boot_median_ci <- function(x, R = 2000, conf = 0.95) {
  meds <- replicate(R, median(sample(x, replace = TRUE)))
  qs   <- quantile(meds, c((1 - conf) / 2, (1 + conf) / 2))
  c(lower = unname(qs[1]), upper = unname(qs[2]))
}

icc_summary <- icc_df |>
  group_by(age) |>
  summarise(
    n        = n(),
    med_icc  = median(icc),
    ci_lo    = boot_median_ci(icc)[["lower"]],
    ci_hi    = boot_median_ci(icc)[["upper"]],
    pct_high = 100 * mean(icc > 0.75),
    pct_low  = 100 * mean(icc < 0.40),
    .groups  = "drop"
  )

wt <- wilcox.test(icc ~ age, data = icc_df)

sub_txt <- sprintf(
  "ICC(3,1) | Young: %.2f [%.2f, %.2f] | Old: %.2f [%.2f, %.2f] | Wilcoxon %s",
  icc_summary$med_icc[1], icc_summary$ci_lo[1], icc_summary$ci_hi[1],
  icc_summary$med_icc[2], icc_summary$ci_lo[2], icc_summary$ci_hi[2],
  fmt_p(wt$p.value)
)

pA_ICC <- ggplot(icc_df, aes(x = age, y = icc, fill = age)) +
  geom_violin(alpha = 0.5, linewidth = 0.3, color = "black", scale = "width") +
  geom_boxplot(width = 0.15, outlier.shape = NA, linewidth = 0.3,
               color = "black", fill = "white", coef = 0) +
  geom_hline(yintercept = 0.75, linetype = "dashed", color = "grey50",
             linewidth = 0.4) +
  annotate("text", x = 2.4, y = 0.77, label = "Good (>0.75)",
           hjust = 0, size = scale_text(BASE_STAT - 1, PA_ICC_W),
           color = "grey50", fontface = "italic") +
  geom_label(data = icc_summary,
             aes(x = age, y = 1.05,
                 label = sprintf("%.2f [%.2f, %.2f]", med_icc, ci_lo, ci_hi)),
             size = scale_text(BASE_COUNT + 0.5, PA_ICC_W),
             fontface = "bold", fill = scales::alpha("white", 0.8),
             linewidth = 0.2, label.padding = unit(1.5, "pt")) +
  scale_fill_manual(values = AGE_COLORS) +
  coord_cartesian(ylim = c(-0.2, 1.1)) +
  labs(title = "Test-Retest Reliability (ICC)",
       subtitle = sub_txt,
       x = NULL, y = "ICC(3,1)",
       tag = "A'") +
  FIG_THEME + theme(legend.position = "none")

write.csv(icc_summary, file.path(DAT_DIR, "audit_panel_A_ICC.csv"),
          row.names = FALSE)

ggsave(file.path(RPT_DIR, "panel_A_ICC.pdf"), pA_ICC,
       width = PA_ICC_W, height = PA_ICC_H, units = "mm", device = pdf_device)
ggsave(file.path(RPT_DIR, "panel_A_ICC.png"), pA_ICC,
       width = PA_ICC_W, height = PA_ICC_H, units = "mm", dpi = 300)
