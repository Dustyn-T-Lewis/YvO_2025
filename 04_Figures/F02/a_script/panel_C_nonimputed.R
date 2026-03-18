# Figure 2 — Panel C (Non-Imputed): Intra-Individual Proteomic Variability
# Complete-case logFC only (both Pre and Post observed). No imputation.
# Annotated with n_complete and %MNAR excluded per subject.
# Outputs: pC_nonimp (ggplot object), panel_C_nonimputed.pdf/.png

setwd(rprojroot::find_rstudio_root_file())
source("04_Figures/F02/a_script/style.R")

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(readr)
  library(ggplot2)
})

PC_NI_W <- 160; PC_NI_H <- 100

RPT_DIR <- "04_Figures/F02/b_reports"
DAT_DIR <- "04_Figures/F02/c_data"
dir.create(RPT_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(DAT_DIR, recursive = TRUE, showWarnings = FALSE)

# Non-imputed normalized data
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

norm_mat <- as.matrix(norm_df[, samp_names])
n_total  <- nrow(norm_mat)

# Load MAR/MNAR classification for annotation
mnar_df <- read_csv("02_Imputation/c_data/02_mar_mnar_classification.csv",
                    show_col_types = FALSE)
mnar_genes <- mnar_df$gene[mnar_df$classification == "MNAR"]

pdf_device <- get_pdf_device()
subjects <- unique(meta$subject)

lfc_list <- lapply(subjects, function(s) {
  pre_id  <- meta$sample_id[meta$subject == s & meta$time == "Pre"]
  post_id <- meta$sample_id[meta$subject == s & meta$time == "Post"]
  if (length(pre_id) != 1 || length(post_id) != 1) return(NULL)

  pre_vals  <- norm_mat[, pre_id]
  post_vals <- norm_mat[, post_id]

  # Complete case: both observed (not NA)
  complete <- !is.na(pre_vals) & !is.na(post_vals)
  lfc <- post_vals[complete] - pre_vals[complete]

  # Count MNAR proteins excluded (missing in at least one timepoint)
  missing_genes <- norm_df$gene[!complete]
  n_mnar_excluded <- sum(missing_genes %in% mnar_genes)

  tibble(
    subject         = s,
    age             = as.character(meta$age[meta$subject == s][1]),
    subj_num        = meta$subj_num[meta$subject == s][1],
    lfc             = as.numeric(lfc),
    n_complete      = sum(complete),
    n_missing       = sum(!complete),
    n_mnar_excluded = n_mnar_excluded,
    pct_mnar_excl   = 100 * n_mnar_excluded / n_total
  )
})

lfc_long <- bind_rows(lfc_list)
lfc_long$age <- factor(lfc_long$age, levels = c("Young", "Old"))

subj_summary <- lfc_long |>
  group_by(subject, age, subj_num, n_complete, n_missing,
           n_mnar_excluded, pct_mnar_excl) |>
  summarise(
    median_lfc = median(lfc, na.rm = TRUE),
    mad_lfc    = mad(lfc, na.rm = TRUE),
    sd_lfc     = sd(lfc, na.rm = TRUE),
    iqr_lfc    = IQR(lfc, na.rm = TRUE),
    .groups    = "drop"
  )

subj_summary <- subj_summary |>
  arrange(age, median_lfc) |>
  mutate(subj_order = factor(subject, levels = unique(subject)))

lfc_long <- lfc_long |>
  mutate(subj_order = factor(subject, levels = levels(subj_summary$subj_order)))

group_summary <- subj_summary |>
  group_by(age) |>
  summarise(
    mean_median    = mean(median_lfc),
    sd_median      = sd(median_lfc),
    mean_complete  = mean(n_complete),
    mean_pct_mnar  = mean(pct_mnar_excl),
    n              = n(),
    .groups        = "drop"
  )

wt <- wilcox.test(median_lfc ~ age, data = subj_summary)
n1 <- sum(subj_summary$age == "Young")
n2 <- sum(subj_summary$age == "Old")
r_rb <- 1 - 2 * wt$statistic / (n1 * n2)

subtitle_text <- sprintf(
  "Complete-case logFC (non-imputed) | n varies by subject (mean %.0f) | Wilcoxon %s",
  mean(subj_summary$n_complete), fmt_p(wt$p.value)
)

# Per-subject n_complete label at bottom of each box
n_label_df <- subj_summary |>
  mutate(label = as.character(n_complete))

pC_nonimp <- ggplot(lfc_long, aes(x = subj_order, y = lfc, fill = age)) +
  geom_boxplot(width = 0.5, linewidth = 0.3, color = "black",
               outlier.shape = NA, alpha = 0.5) +
  geom_text(data = n_label_df,
            aes(x = subj_order, y = -1.45, label = label),
            inherit.aes = FALSE, size = scale_text(BASE_COUNT - 1.5, PC_NI_W),
            color = "grey50", angle = 0) +
  facet_grid(~ age, scales = "free_x", space = "free_x") +
  coord_cartesian(ylim = c(-1.5, 1.5)) +
  scale_fill_manual(values = AGE_COLORS) +
  labs(x = "Subject",
       y = expression(bold(Delta~log[2]*"FC (Post/Pre)")),
       title = "Intra-Individual Variability (Non-Imputed)",
       subtitle = subtitle_text,
       tag = "C") +
  FIG_THEME +
  theme(legend.position = "none",
        panel.spacing = unit(3, "mm"),
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5,
                                   size = FIG_AXIS_TEXT - 1.5))

write.csv(subj_summary |>
            select(subject, age, subj_num, n_complete, n_missing,
                   n_mnar_excluded, pct_mnar_excl, median_lfc,
                   mad_lfc, sd_lfc, iqr_lfc),
          file.path(DAT_DIR, "audit_panel_C_nonimputed.csv"),
          row.names = FALSE)

ggsave(file.path(RPT_DIR, "panel_C_nonimputed.pdf"), pC_nonimp,
       width = PC_NI_W, height = PC_NI_H, units = "mm", device = pdf_device)
ggsave(file.path(RPT_DIR, "panel_C_nonimputed.png"), pC_nonimp,
       width = PC_NI_W, height = PC_NI_H, units = "mm", dpi = 300)
