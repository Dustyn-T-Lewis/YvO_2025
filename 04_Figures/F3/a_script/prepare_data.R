# F3 prepare_data.R — run BEFORE panel scripts
# Melov permutation, reversal contingency + Fisher, signed reversal score, fGSEA cache

setwd(rprojroot::find_rstudio_root_file())
source("04_Figures/shared/style.R")
source("04_Figures/shared/pathway_utils.R")

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(tidyr)
  library(msigdbr)
  library(fgsea)
})

set.seed(42)

DAT <- "04_Figures/F3/c_data"
RPT <- "04_Figures/F3/b_reports"
dir.create(DAT, recursive = TRUE, showWarnings = FALSE)
dir.create(RPT, recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(DAT, "reversal_tests"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(DAT, "shared"), recursive = TRUE, showWarnings = FALSE)

dep_df <- read_csv("03_DEP/c_data/03_combined_results.csv", show_col_types = FALSE)
stopifnot(nrow(dep_df) > 2000)
stopifnot("logFC_Reversal" %in% names(dep_df))

imp_data <- read_csv("02_Imputation/c_data/01_imputed.csv", show_col_types = FALSE)
imp_ann_cols <- c("uniprot_id", "protein", "gene", "description")
imp_samp_cols <- setdiff(names(imp_data), imp_ann_cols)
imp_mat_rev <- as.matrix(imp_data[, imp_samp_cols])
rownames(imp_mat_rev) <- imp_data$uniprot_id

dal_meta_rev <- as.data.frame(
  readRDS("02_Imputation/c_data/01_DAList_imputed.rds")$metadata)
meta_rev <- tibble(
  sample_id = dal_meta_rev$Col_ID,
  subject   = sub("_(Pre|Post)$", "", dal_meta_rev$Col_ID),
  age       = dal_meta_rev$Group,
  time      = dal_meta_rev$Timepoint,
  group     = dal_meta_rev$Group_Time
)

aging_sig <- dep_df %>%
  filter(!is.na(P.Value_Aging) & P.Value_Aging < 0.05) %>%
  pull(uniprot_id)
aging_sig <- intersect(aging_sig, rownames(imp_mat_rev))
n_aging <- length(aging_sig)
message(sprintf("  Aging signature: %d proteins (nominal P < 0.05)", n_aging))

young_pre_ids  <- meta_rev$sample_id[meta_rev$age == "Young" & meta_rev$time == "Pre"]
old_pre_ids    <- meta_rev$sample_id[meta_rev$age == "Old"   & meta_rev$time == "Pre"]
old_post_ids   <- meta_rev$sample_id[meta_rev$age == "Old"   & meta_rev$time == "Post"]

young_pre_mean <- rowMeans(imp_mat_rev[aging_sig, intersect(young_pre_ids, colnames(imp_mat_rev)), drop = FALSE], na.rm = TRUE)
old_pre_mean   <- rowMeans(imp_mat_rev[aging_sig, intersect(old_pre_ids,   colnames(imp_mat_rev)), drop = FALSE], na.rm = TRUE)
old_post_mean  <- rowMeans(imp_mat_rev[aging_sig, intersect(old_post_ids,  colnames(imp_mat_rev)), drop = FALSE], na.rm = TRUE)

d_pre  <- sqrt(sum((old_pre_mean  - young_pre_mean)^2))
d_post <- sqrt(sum((old_post_mean - young_pre_mean)^2))
reversal_pct <- (d_pre - d_post) / d_pre * 100

set.seed(42)
n_perm <- 10000
perm_deltas <- numeric(n_perm)

old_subjects <- unique(meta_rev$subject[meta_rev$age == "Old"])
old_pre_meta  <- meta_rev %>% filter(age == "Old", time == "Pre")
old_post_meta <- meta_rev %>% filter(age == "Old", time == "Post")

for (i in seq_len(n_perm)) {
  swap <- sample(c(TRUE, FALSE), length(old_subjects), replace = TRUE)
  perm_pre_ids  <- character(0)
  perm_post_ids <- character(0)

  for (j in seq_along(old_subjects)) {
    subj <- old_subjects[j]
    pre_id  <- old_pre_meta$sample_id[old_pre_meta$subject == subj]
    post_id <- old_post_meta$sample_id[old_post_meta$subject == subj]
    if (swap[j]) {
      perm_pre_ids  <- c(perm_pre_ids, post_id)
      perm_post_ids <- c(perm_post_ids, pre_id)
    } else {
      perm_pre_ids  <- c(perm_pre_ids, pre_id)
      perm_post_ids <- c(perm_post_ids, post_id)
    }
  }

  perm_pre_mean  <- rowMeans(imp_mat_rev[aging_sig, intersect(perm_pre_ids, colnames(imp_mat_rev)), drop = FALSE], na.rm = TRUE)
  perm_post_mean <- rowMeans(imp_mat_rev[aging_sig, intersect(perm_post_ids, colnames(imp_mat_rev)), drop = FALSE], na.rm = TRUE)
  d_pre_perm  <- sqrt(sum((perm_pre_mean  - young_pre_mean)^2))
  d_post_perm <- sqrt(sum((perm_post_mean - young_pre_mean)^2))
  perm_deltas[i] <- d_pre_perm - d_post_perm
}

observed_delta <- d_pre - d_post
perm_pvalue <- mean(perm_deltas >= observed_delta)

n_exceed <- sum(perm_deltas >= observed_delta)
perm_pval_ci <- binom.test(n_exceed, n_perm)$conf.int

set.seed(42)
n_boot_rev <- 2000
boot_rev_pct <- replicate(n_boot_rev, {
  idx <- sample(seq_along(aging_sig), replace = TRUE)
  boot_d_pre  <- sqrt(sum((old_pre_mean[idx]  - young_pre_mean[idx])^2))
  boot_d_post <- sqrt(sum((old_post_mean[idx] - young_pre_mean[idx])^2))
  (boot_d_pre - boot_d_post) / boot_d_pre * 100
})
rev_pct_ci <- quantile(boot_rev_pct, c(0.025, 0.975))

message(sprintf("  Permutation p-value: %.4f [%.4f, %.4f]",
                perm_pvalue, perm_pval_ci[1], perm_pval_ci[2]))
message(sprintf("  Reversal %%: %.1f%% [%.1f%%, %.1f%%]",
                reversal_pct, rev_pct_ci[1], rev_pct_ci[2]))

melov_df <- tibble(
  d_pre = d_pre, d_post = d_post,
  reversal_pct = round(reversal_pct, 2),
  reversal_pct_ci_lower = round(rev_pct_ci[1], 2),
  reversal_pct_ci_upper = round(rev_pct_ci[2], 2),
  observed_delta = observed_delta,
  p_value = perm_pvalue,
  p_value_ci_lower = round(perm_pval_ci[1], 6),
  p_value_ci_upper = round(perm_pval_ci[2], 6),
  n_aging_proteins = n_aging,
  n_permutations = n_perm,
  n_boot_reversal_pct = n_boot_rev
)
write_csv(melov_df, file.path(DAT, "reversal_tests", "melov_permutation.csv"))

p_melov <- ggplot(tibble(delta = perm_deltas), aes(x = delta)) +
  geom_histogram(bins = 50, fill = "grey70", color = "grey50", linewidth = 0.2) +
  geom_vline(xintercept = observed_delta, color = "#D6604D",
             linewidth = 1, linetype = "solid") +
  annotate("text", x = observed_delta, y = Inf, vjust = 1.5, hjust = -0.1,
           label = sprintf("Observed = %.3f\np = %.4f [%.4f, %.4f]\nReversal = %.1f%% [%.1f, %.1f]",
                           observed_delta, perm_pvalue,
                           perm_pval_ci[1], perm_pval_ci[2],
                           reversal_pct, rev_pct_ci[1], rev_pct_ci[2]),
           size = 2.8, fontface = "bold", color = "#D6604D") +
  labs(title = "Melov Reversal Permutation Test",
       subtitle = sprintf("d(Old_Pre, Young_Pre) - d(Old_Post, Young_Pre) | %d aging-signature proteins",
                          n_aging),
       x = expression(d[pre] - d[post]), y = "Count (10,000 permutations)") +
  FIG_THEME +
  theme(plot.title = element_text(size = 10, face = "bold"))

ggsave(file.path(RPT, "melov_reversal_permutation.pdf"), p_melov,
       width = 180, height = 120, units = "mm", device = get_pdf_device())

aging_proteins <- dep_df %>%
  filter(!is.na(P.Value_Aging) & P.Value_Aging < 0.05 &
         !is.na(logFC_Aging) & !is.na(logFC_Training_Old))

contingency <- aging_proteins %>%
  mutate(
    aging_dir    = ifelse(logFC_Aging > 0, "Aging_Up", "Aging_Down"),
    training_dir = ifelse(logFC_Training_Old > 0, "Training_Up", "Training_Down"),
    pattern      = case_when(
      abs(logFC_Training_Old) <= 0.2 ~ "Negligible",
      sign(logFC_Aging) != sign(logFC_Training_Old) ~ "Reversed",
      TRUE ~ "Exacerbated"
    )
  )

ct <- table(contingency$aging_dir, contingency$training_dir)
fisher_res <- fisher.test(ct)
fisher_or_ci <- fisher_res$conf.int

n_rev_total <- sum(contingency$pattern == "Reversed")
rev_binom <- binom.test(n_rev_total, nrow(contingency))
rev_pct_binom_ci <- rev_binom$conf.int * 100

contingency_summary <- tibble(
  aging_up_training_down   = sum(contingency$aging_dir == "Aging_Up" &
                                  contingency$training_dir == "Training_Down"),
  aging_up_training_up     = sum(contingency$aging_dir == "Aging_Up" &
                                  contingency$training_dir == "Training_Up"),
  aging_down_training_up   = sum(contingency$aging_dir == "Aging_Down" &
                                  contingency$training_dir == "Training_Up"),
  aging_down_training_down = sum(contingency$aging_dir == "Aging_Down" &
                                  contingency$training_dir == "Training_Down"),
  n_reversed    = n_rev_total,
  n_exacerbated = sum(contingency$pattern == "Exacerbated"),
  n_negligible  = sum(contingency$pattern == "Negligible"),
  pct_reversed  = round(mean(contingency$pattern == "Reversed") * 100, 1),
  pct_reversed_ci_lower = round(rev_pct_binom_ci[1], 1),
  pct_reversed_ci_upper = round(rev_pct_binom_ci[2], 1),
  fisher_or     = round(fisher_res$estimate, 3),
  fisher_or_ci_lower = round(fisher_or_ci[1], 3),
  fisher_or_ci_upper = round(fisher_or_ci[2], 3),
  fisher_p      = fisher_res$p.value,
  n_aging_proteins = nrow(contingency)
)
write_csv(contingency_summary, file.path(DAT, "reversal_tests", "reversal_contingency.csv"))

reversal_cor <- dep_df %>%
  filter(!is.na(logFC_Aging) & !is.na(logFC_Training_Old))

cor_res <- cor.test(reversal_cor$logFC_Aging, reversal_cor$logFC_Training_Old,
                    method = "pearson")

signed_reversal <- tibble(
  r = round(cor_res$estimate, 4),
  ci_lower = round(cor_res$conf.int[1], 4),
  ci_upper = round(cor_res$conf.int[2], 4),
  p_value = cor_res$p.value,
  n_proteins = nrow(reversal_cor),
  interpretation = ifelse(cor_res$estimate < 0,
    "Negative: training opposes aging globally",
    "Positive: training reinforces aging direction")
)
write_csv(signed_reversal, file.path(DAT, "reversal_tests", "signed_reversal_score.csv"))

fgsea_cache <- file.path(DAT, "shared", "fgsea_tstat_all_v2.csv")
if (!file.exists(fgsea_cache)) {
  f2_cache <- "04_Figures/F2/c_data/shared/fgsea_tstat_all_v2.csv"
  f1_cache <- "04_Figures/F1/c_data/06_panel_F_fgsea_results.csv"
  if (file.exists(f2_cache)) {
    file.copy(f2_cache, fgsea_cache)
  } else if (file.exists(f1_cache)) {
    file.copy(f1_cache, fgsea_cache)
  } else {
    stop("fGSEA cache not found")
  }
}
fgsea_all <- read_csv(fgsea_cache, show_col_types = FALSE)

if (!"Reversal" %in% unique(fgsea_all$contrast)) {
  message("Computing fGSEA for Reversal contrast (unified H + C2:CP + GO:BP)...")

  reversal_stats <- dep_df %>%
    filter(!is.na(t_Reversal)) %>%
    distinct(gene, .keep_all = TRUE) %>%
    { setNames(.$t_Reversal, .$gene) } %>%
    sort(decreasing = TRUE)

  pw_collection_rev <- build_pathway_collection(min_size = 10, max_size = 500)

  set.seed(42)
  fgsea_reversal <- run_fgsea_deduplicated(
    ranks          = reversal_stats,
    pathways       = pw_collection_rev,
    jaccard_cutoff = 0.5,
    nperm          = 10000,
    min_size       = 10,
    max_size       = 500
  ) %>%
    mutate(contrast = "Reversal",
           leadingEdge = sapply(leadingEdge, paste, collapse = ";")) %>%
    dplyr::select(pathway, pval, padj, log2err, ES, NES, size,
                  leadingEdge, contrast, database)

  fgsea_all <- bind_rows(fgsea_all, fgsea_reversal)
  write_csv(fgsea_all, fgsea_cache)
  message(sprintf("  Added %d Reversal fGSEA results, updated cache",
                  nrow(fgsea_reversal)))
}

message("F3 prepare_data.R complete")
