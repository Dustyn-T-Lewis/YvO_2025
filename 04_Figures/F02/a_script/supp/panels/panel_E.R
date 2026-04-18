# F02 Supp Panel E: Intra-Individual Proteomic Variability (Imputed)
# One boxplot per subject, faceted by Young/Old, ordered by median log2FC.
# Annotated with per-subject imputation fractions.
# Outputs: pE (ggplot object), SUPP_panel_E_imputed.{pdf,png}

setwd(rprojroot::find_rstudio_root_file())
source("04_Figures/F02/a_script/style.R")

library(dplyr)
library(tidyr)
library(stringr)
library(readr)
library(ggplot2)

PE_W <- 160
PE_H <- 90

RPT_PNG <- "04_Figures/F02/b_reports/supp/png/panels"
RPT_PDF <- "04_Figures/F02/b_reports/supp/pdf/panels"
DAT_DIR <- "04_Figures/F02/c_data"
dir.create(RPT_PNG, recursive = TRUE, showWarnings = FALSE)
dir.create(RPT_PDF, recursive = TRUE, showWarnings = FALSE)
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
meta$age  <- factor(meta$age,  levels = c("Young", "Old"))
meta$time <- factor(meta$time, levels = c("Pre", "Post"))

imp_mat    <- as.matrix(imp_df[, samp_names])
n_proteins <- nrow(imp_mat)

# Imputation mask (TRUE = was missing, i.e. imputed)
mask_df  <- read_csv("02_Imputation/c_data/03_imputation_mask.csv",
                     show_col_types = FALSE)
mask_mat <- as.matrix(mask_df[, intersect(samp_names, names(mask_df))])
rownames(mask_mat) <- mask_df$gene

# MAR/MNAR classification
mnar_df    <- read_csv("02_Imputation/c_data/02_mar_mnar_classification.csv",
                       show_col_types = FALSE)
mnar_genes <- mnar_df$gene[mnar_df$classification == "MNAR"]

pdf_device <- get_pdf_device()
subjects   <- unique(meta$subject)

lfc_list <- lapply(subjects, function(s) {
  pre_id  <- meta$sample_id[meta$subject == s & meta$time == "Pre"]
  post_id <- meta$sample_id[meta$subject == s & meta$time == "Post"]
  if (length(pre_id) != 1 || length(post_id) != 1) return(NULL)

  lfc <- imp_mat[, post_id] - imp_mat[, pre_id]  # log2(Post) - log2(Pre)

  pre_imputed    <- if (pre_id  %in% colnames(mask_mat)) mask_mat[, pre_id]  else rep(FALSE, n_proteins)
  post_imputed   <- if (post_id %in% colnames(mask_mat)) mask_mat[, post_id] else rep(FALSE, n_proteins)
  either_imputed <- pre_imputed | post_imputed
  n_imputed_lfc  <- sum(either_imputed)
  pct_imputed    <- 100 * n_imputed_lfc / n_proteins

  genes            <- imp_df$gene
  mnar_and_imputed <- either_imputed & (genes %in% mnar_genes)
  n_mnar_lfc       <- sum(mnar_and_imputed)
  pct_mnar         <- 100 * n_mnar_lfc / n_proteins

  tibble(
    subject       = s,
    age           = as.character(meta$age[meta$subject == s][1]),
    subj_num      = meta$subj_num[meta$subject == s][1],
    lfc           = as.numeric(lfc),
    n_imputed_lfc = n_imputed_lfc,
    pct_imputed   = pct_imputed,
    n_mnar_lfc    = n_mnar_lfc,
    pct_mnar      = pct_mnar
  )
})

lfc_long     <- bind_rows(lfc_list)
lfc_long$age <- factor(lfc_long$age, levels = c("Young", "Old"))

subj_summary <- lfc_long |>
  group_by(subject, age, subj_num, n_imputed_lfc, pct_imputed,
           n_mnar_lfc, pct_mnar) |>
  summarise(
    median_lfc = median(lfc, na.rm = TRUE),
    mad_lfc    = mad(lfc, na.rm = TRUE),
    sd_lfc     = sd(lfc, na.rm = TRUE),
    iqr_lfc    = IQR(lfc, na.rm = TRUE),
    q25        = quantile(lfc, 0.25, na.rm = TRUE),
    q75        = quantile(lfc, 0.75, na.rm = TRUE),
    n_proteins = n(),
    .groups    = "drop"
  ) |>
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

wt   <- wilcox.test(median_lfc ~ age, data = subj_summary)
n1   <- sum(subj_summary$age == "Young")
n2   <- sum(subj_summary$age == "Old")
r_rb <- 1 - 2 * wt$statistic / (n1 * n2)

mean_pct_imp  <- mean(subj_summary$pct_imputed)
mean_pct_mnar <- mean(subj_summary$pct_mnar)
subtitle_text <- sprintf(
  paste0("Per-subject delta-log2FC (Post - Pre) | ",
         "%s missForest-imputed proteins\n",
         "%.0f%% imputed (%.0f%% MNAR-derived) | ",
         "age effect Wilcoxon %s"),
  format(n_proteins, big.mark = ","), mean_pct_imp, mean_pct_mnar, fmt_p(wt$p.value)
)

pE <- ggplot(lfc_long, aes(x = subj_order, y = lfc, fill = age)) +
  geom_boxplot(width = 0.5, linewidth = 0.3, color = "black",
               outlier.shape = NA, alpha = 0.5) +
  facet_grid(~ age, scales = "free_x", space = "free_x") +
  coord_cartesian(ylim = c(-1.5, 1.5)) +
  scale_fill_manual(values = AGE_COLORS) +
  labs(x = "Subject",
       y = expression(bold(Delta~log[2]*"FC (Post/Pre)")),
       title = "Intra-Individual Proteomic Variability",
       subtitle = subtitle_text,
       tag = "E") +
  FIG_THEME +
  theme(legend.position = "none",
        panel.spacing   = unit(3, "mm"),
        axis.text.x     = element_text(angle = 90, hjust = 1, vjust = 0.5,
                                       size = FIG_AXIS_TEXT - 1.5),
        plot.title      = element_text(margin = margin(b = 0)),
        plot.subtitle   = element_text(size = FIG_SUBTITLE_SIZE - 1.0,
                                       face = "bold.italic", color = "grey40",
                                       margin = margin(t = 0, b = 1)),
        strip.text      = element_text(face = "bold", size = FIG_STRIP_SIZE,
                                       margin = margin(t = 0, b = 1)),
        plot.margin     = margin(t = 0, r = 5.5, b = 5.5, l = 5.5))

write.csv(subj_summary |>
            select(subject, age, subj_num, median_lfc, mad_lfc, sd_lfc,
                   iqr_lfc, q25, q75, n_proteins,
                   n_imputed_lfc, pct_imputed, n_mnar_lfc, pct_mnar),
          file.path(DAT_DIR, "audit_panel_E_imputed.csv"),
          row.names = FALSE)

write.csv(group_summary,
          file.path(DAT_DIR, "audit_panel_E_wilcoxon.csv"),
          row.names = FALSE)

ggsave(file.path(RPT_PNG, "SUPP_panel_E_imputed.png"), pE,
       width = PE_W, height = PE_H, units = "mm", dpi = 300)
ggsave(file.path(RPT_PDF, "SUPP_panel_E_imputed.pdf"), pE,
       width = PE_W, height = PE_H, units = "mm", device = pdf_device)

cat("F02 Supp Panel E done\n")
