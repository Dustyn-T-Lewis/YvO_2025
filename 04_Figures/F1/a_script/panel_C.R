# Figure 1 — Panel C: Intra-Individual Proteomic Variability
# One boxplot per subject, faceted by Young/Old, ordered by median log2FC.
# Outputs: pC (ggplot object), panel_C_intra_variability.pdf/.png

setwd(rprojroot::find_rstudio_root_file())
source("04_Figures/shared/style.R")

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(readr)
  library(ggplot2)
})

PC_W <- 160; PC_H <- 90

RPT_DIR <- "04_Figures/F1/b_reports"
DAT_DIR <- "04_Figures/F1/c_data"
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
    subject  = paste0(prefix, "_", subj_num)
  )
meta$age  <- factor(meta$age, levels = c("Young", "Old"))
meta$time <- factor(meta$time, levels = c("Pre", "Post"))

imp_mat <- as.matrix(imp_df[, samp_names])
n_proteins <- nrow(imp_mat)

pdf_device <- get_pdf_device()
subjects <- unique(meta$subject)

lfc_list <- lapply(subjects, function(s) {
  pre_id  <- meta$sample_id[meta$subject == s & meta$time == "Pre"]
  post_id <- meta$sample_id[meta$subject == s & meta$time == "Post"]
  if (length(pre_id) != 1 || length(post_id) != 1) return(NULL)

  lfc <- imp_mat[, post_id] - imp_mat[, pre_id]  # log2(Post) - log2(Pre)
  tibble(
    subject  = s,
    age      = as.character(meta$age[meta$subject == s][1]),
    subj_num = meta$subj_num[meta$subject == s][1],
    lfc      = as.numeric(lfc)
  )
})

lfc_long <- bind_rows(lfc_list)
lfc_long$age <- factor(lfc_long$age, levels = c("Young", "Old"))

subj_summary <- lfc_long |>
  group_by(subject, age, subj_num) |>
  summarise(
    median_lfc = median(lfc, na.rm = TRUE),
    mad_lfc    = mad(lfc, na.rm = TRUE),
    sd_lfc     = sd(lfc, na.rm = TRUE),
    iqr_lfc    = IQR(lfc, na.rm = TRUE),
    q25        = quantile(lfc, 0.25, na.rm = TRUE),
    q75        = quantile(lfc, 0.75, na.rm = TRUE),
    n_proteins = n(),
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
    mean_median = mean(median_lfc),
    sd_median   = sd(median_lfc),
    n           = n(),
    .groups     = "drop"
  )

wt <- wilcox.test(median_lfc ~ age, data = subj_summary)
wt_label <- fmt_p(wt$p.value)

subtitle_text <- sprintf(
  "%d proteins per subject | Wilcoxon %s",
  n_proteins, wt_label
)

pC <- ggplot(lfc_long, aes(x = subj_order, y = lfc, fill = age)) +
  geom_boxplot(width = 0.5, linewidth = 0.3, color = "black",
               outlier.shape = NA, alpha = 0.5) +
  facet_grid(~ age, scales = "free_x", space = "free_x") +
  coord_cartesian(ylim = c(-1.5, 1.5)) +
  scale_fill_manual(values = AGE_COLORS) +
  labs(x = "Subject",
       y = expression(Delta~log[2]*"FC (Post/Pre)"),
       title = "Intra-Individual Proteomic Variability",
       subtitle = subtitle_text,
       tag = "C") +
  FIG_THEME +
  theme(legend.position = "none",
        panel.spacing = unit(8, "mm"),
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5,
                                   size = FIG_AXIS_TEXT - 1.5))

write.csv(subj_summary |>
            select(subject, age, subj_num, median_lfc, mad_lfc, sd_lfc,
                   iqr_lfc, q25, q75, n_proteins),
          file.path(DAT_DIR, "audit_panel_C_intra_variability.csv"),
          row.names = FALSE)

write.csv(group_summary,
          file.path(DAT_DIR, "audit_panel_C_wilcoxon.csv"),
          row.names = FALSE)

ggsave(file.path(RPT_DIR, "panel_C_intra_variability.pdf"), pC,
       width = PC_W, height = PC_H, units = "mm", device = pdf_device)
ggsave(file.path(RPT_DIR, "panel_C_intra_variability.png"), pC,
       width = PC_W, height = PC_H, units = "mm", dpi = 300)
