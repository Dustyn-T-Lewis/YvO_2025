################################################################################
#   Figure 3 — Shared Setup
#   Common packages, style constants, helpers, and data loading used by all
#   per-panel scripts (panel_A.R through panel_E.R).
#
#   Central question: Does resistance training reverse aging?
#   Contrasts of interest: Aging vs Training_Old (and the Reversal contrast)
#
#   This script is idempotent: sourcing it multiple times has no side effects.
################################################################################

# === 1. PACKAGES ==============================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(patchwork)
  library(ggrepel)
  library(scales)
  library(grid)
  library(msigdbr)
  library(fgsea)
  library(rrvgo)
  library(GOSemSim)
  library(org.Hs.eg.db)
  library(clusterProfiler)
  library(ggforce)
  library(openxlsx)
  library(boot)
})

# === 2. SETUP =================================================================

set.seed(42)
setwd(rprojroot::find_rstudio_root_file())

source("05_Figures/shared/style.R")

RPT_DIR <- "05_Figures/F3/b_reports"
DAT_DIR <- "05_Figures/F3/c_data"
dir.create(RPT_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(DAT_DIR, recursive = TRUE, showWarnings = FALSE)
for (pnl in c("panel_A", "panel_B", "panel_C", "panel_D",
              "panel_E", "shared", "reversal_tests"))
  dir.create(file.path(DAT_DIR, pnl), recursive = TRUE, showWarnings = FALSE)

# === 3. F3-SPECIFIC PALETTES =================================================

SIG_COLORS <- c(
  "Reversal"           = "#7B5EA7",
  "Sig Both"           = "#2E7D32",
  "Sig Aging only"     = "#E05A4E",
  "Sig Training only"  = "#5DA5DA",
  "NS"                 = "grey70"
)

SIG_LABEL_FILL <- c(
  "Reversal"           = alpha("#7B5EA7", 0.75),
  "Sig Both"           = alpha("#2E7D32", 0.75),
  "Sig Aging only"     = alpha("#E05A4E", 0.75),
  "Sig Training only"  = alpha("#5DA5DA", 0.75),
  "NS"                 = alpha("grey70",  0.75)
)
SIG_LABEL_TEXT <- c(
  "Reversal"           = "white",
  "Sig Both"           = "white",
  "Sig Aging only"     = "white",
  "Sig Training only"  = "white",
  "NS"                 = "white"
)

# === 4. HELPERS ===============================================================

classify_proteins_f3 <- function(pi_aging, pi_training_old, pi_reversal,
                                 threshold = 0.05) {
  case_when(
    pi_reversal < threshold                          ~ "Reversal",
    pi_aging < threshold & pi_training_old < threshold ~ "Sig Both",
    pi_aging < threshold                             ~ "Sig Aging only",
    pi_training_old < threshold                      ~ "Sig Training only",
    TRUE                                             ~ "NS"
  ) %>%
    factor(levels = c("Reversal", "Sig Both",
                       "Sig Aging only", "Sig Training only", "NS"))
}

classify_reversal_pattern <- function(logFC_Aging, logFC_Training_Old,
                                      pi_aging, threshold = 0.05) {
  sig_A <- pi_aging < threshold
  case_when(
    sign(logFC_Aging) != sign(logFC_Training_Old) ~ "Reversed",
    sig_A & abs(logFC_Training_Old) < 0.5 * abs(logFC_Aging) ~ "Attenuated",
    sig_A & abs(logFC_Training_Old) > abs(logFC_Aging) ~ "Amplified",
    TRUE ~ "Concordant"
  ) %>%
    factor(levels = c("Reversed", "Attenuated", "Amplified", "Concordant"))
}

PATTERN_COLORS <- c(
  Reversed   = "#E05A4E",
  Attenuated = "#5DA5DA",
  Amplified  = "#F4A460",
  Concordant = "#78909C"
)

# === 5. DATA LOADING ==========================================================

message("Loading data...")
dep_df <- read_csv("03_DEP/c_data/03_combined_results.csv", show_col_types = FALSE)
stopifnot(nrow(dep_df) > 2000)
stopifnot("logFC_Reversal" %in% names(dep_df))

imputation_df <- read_csv("02_Imputation/c_data/02_mar_mnar_classification.csv",
                          show_col_types = FALSE) %>%
  transmute(gene, imputed = classification != "Complete")
message(sprintf("  %d proteins with imputation status (%d imputed)",
                nrow(imputation_df), sum(imputation_df$imputed)))

fgsea_cache <- "04_Figures/F2/c_data/shared/fgsea_tstat_all_v2.csv"
if (!file.exists(fgsea_cache)) stop("fGSEA cache not found at ", fgsea_cache)
fgsea_all <- read_csv(fgsea_cache, show_col_types = FALSE)

message(sprintf("Loaded %d proteins, %d fGSEA results", nrow(dep_df), nrow(fgsea_all)))

# === 6. REVERSAL fGSEA (compute if not in F2 cache) ==========================

if (!"Reversal" %in% unique(fgsea_all$contrast)) {
  message("Computing fGSEA for Reversal contrast...")

  reversal_stats <- dep_df %>%
    filter(!is.na(t_Reversal)) %>%
    distinct(gene, .keep_all = TRUE) %>%
    { setNames(.$t_Reversal, .$gene) } %>%
    sort(decreasing = TRUE)

  hallmark_gs <- msigdbr(species = "Homo sapiens", collection = "H") %>%
    split(x = .$gene_symbol, f = .$gs_name)
  gobp_gs <- msigdbr(species = "Homo sapiens", collection = "C5",
                      subcollection = "GO:BP") %>%
    split(x = .$gene_symbol, f = .$gs_name)

  fgsea_rev_h <- fgsea(pathways = hallmark_gs, stats = reversal_stats,
                        minSize = 10, maxSize = 500) %>%
    as_tibble() %>% mutate(contrast = "Reversal", database = "Hallmark")
  fgsea_rev_g <- fgsea(pathways = gobp_gs, stats = reversal_stats,
                        minSize = 10, maxSize = 500) %>%
    as_tibble() %>% mutate(contrast = "Reversal", database = "GO:BP")

  fgsea_reversal <- bind_rows(fgsea_rev_h, fgsea_rev_g) %>%
    mutate(leadingEdge = sapply(leadingEdge, paste, collapse = ";")) %>%
    dplyr::select(pathway, pval, padj, log2err, ES, NES, size,
                  leadingEdge, contrast, database)

  fgsea_all <- bind_rows(fgsea_all, fgsea_reversal)
  message(sprintf("  Added %d Reversal fGSEA results (%d sig at padj < 0.05)",
                  nrow(fgsea_reversal),
                  sum(fgsea_reversal$padj < 0.05, na.rm = TRUE)))
}

# === 7. REVERSAL TESTS (Melov permutation, contingency, signed score) ========
# References:
#   - Melov et al. 2007 (PMID 17520024): Euclidean distance reversal test
#   - Fisher's exact test on 2x2 reversal contingency table
#   - Signed reversal score: Pearson r(logFC_Aging, logFC_Training_Old)

message("Computing reversal statistical tests...")

# --- 7a. Melov Permutation Reversal Test ------------------------------------

message("  Melov permutation reversal test...")

imp_data <- read_csv("02_Imputation/c_data/01_imputed.csv", show_col_types = FALSE)
imp_ann_cols <- c("uniprot_id", "protein", "gene", "description")
imp_samp_cols <- setdiff(names(imp_data), imp_ann_cols)
imp_mat_rev <- as.matrix(imp_data[, imp_samp_cols])
rownames(imp_mat_rev) <- imp_data$uniprot_id

# Load metadata for group assignments
dal <- readRDS("02_Imputation/c_data/01_DAList_imputed.rds")
dal_meta_rev <- as.data.frame(dal$metadata)
meta_rev <- tibble(
  sample_id = dal_meta_rev$Col_ID,
  subject   = sub("_(Pre|Post)$", "", dal_meta_rev$Col_ID),
  age       = dal_meta_rev$Group,
  time      = dal_meta_rev$Timepoint,
  group     = dal_meta_rev$Group_Time
)

# Aging signature: proteins with nominal P < 0.05 from Aging contrast
aging_sig <- dep_df %>%
  filter(!is.na(P.Value_Aging) & P.Value_Aging < 0.05) %>%
  pull(uniprot_id)
aging_sig <- intersect(aging_sig, rownames(imp_mat_rev))
n_aging <- length(aging_sig)
message(sprintf("    Aging signature: %d proteins (nominal P < 0.05)", n_aging))

# Compute group means per protein
young_pre_ids  <- meta_rev$sample_id[meta_rev$age == "Young" & meta_rev$time == "Pre"]
old_pre_ids    <- meta_rev$sample_id[meta_rev$age == "Old"   & meta_rev$time == "Pre"]
old_post_ids   <- meta_rev$sample_id[meta_rev$age == "Old"   & meta_rev$time == "Post"]

young_pre_mean <- rowMeans(imp_mat_rev[aging_sig, intersect(young_pre_ids, colnames(imp_mat_rev)), drop = FALSE], na.rm = TRUE)
old_pre_mean   <- rowMeans(imp_mat_rev[aging_sig, intersect(old_pre_ids,   colnames(imp_mat_rev)), drop = FALSE], na.rm = TRUE)
old_post_mean  <- rowMeans(imp_mat_rev[aging_sig, intersect(old_post_ids,  colnames(imp_mat_rev)), drop = FALSE], na.rm = TRUE)

# Euclidean distances
d_pre  <- sqrt(sum((old_pre_mean  - young_pre_mean)^2))
d_post <- sqrt(sum((old_post_mean - young_pre_mean)^2))
reversal_pct <- (d_pre - d_post) / d_pre * 100

# Permutation test: shuffle Old_Pre vs Old_Post labels within old subjects
set.seed(42)
n_perm <- 10000
perm_deltas <- numeric(n_perm)

old_subjects  <- unique(meta_rev$subject[meta_rev$age == "Old"])
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

# 95% CI on permutation p-value (Clopper-Pearson exact binomial)
n_exceed <- sum(perm_deltas >= observed_delta)
perm_pval_ci <- binom.test(n_exceed, n_perm)$conf.int

# Bootstrap 95% CI on reversal % via resampling aging-signature proteins
set.seed(42)
n_boot_rev <- 2000
boot_rev_pct <- replicate(n_boot_rev, {
  idx <- sample(seq_along(aging_sig), replace = TRUE)
  boot_d_pre  <- sqrt(sum((old_pre_mean[idx]  - young_pre_mean[idx])^2))
  boot_d_post <- sqrt(sum((old_post_mean[idx] - young_pre_mean[idx])^2))
  (boot_d_pre - boot_d_post) / boot_d_pre * 100
})
rev_pct_ci <- quantile(boot_rev_pct, c(0.025, 0.975))

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
write_csv(melov_df, file.path(DAT_DIR, "reversal_tests", "melov_permutation.csv"))

# Histogram of permuted deltas
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

ggsave(file.path(RPT_DIR, "melov_reversal_permutation.pdf"), p_melov,
       width = 180, height = 120, units = "mm", device = pdf)

message(sprintf("    Melov: d_pre=%.4f, d_post=%.4f, reversal=%.1f%% [%.1f, %.1f], p=%.4f [%.4f, %.4f]",
                d_pre, d_post, reversal_pct, rev_pct_ci[1], rev_pct_ci[2],
                perm_pvalue, perm_pval_ci[1], perm_pval_ci[2]))

# --- 7b. Reversal Contingency Table + Fisher's Exact Test -------------------

message("  Reversal contingency table + Fisher's exact test...")

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
write_csv(contingency_summary, file.path(DAT_DIR, "reversal_tests", "reversal_contingency.csv"))

message(sprintf("    Fisher p=%.4g, OR=%.2f [%.2f, %.2f], Reversed=%d/%d (%.1f%%)",
                fisher_res$p.value, fisher_res$estimate,
                fisher_or_ci[1], fisher_or_ci[2],
                n_rev_total, nrow(contingency), contingency_summary$pct_reversed))

# --- 7c. Signed Reversal Score -----------------------------------------------

message("  Signed reversal score...")

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
write_csv(signed_reversal, file.path(DAT_DIR, "reversal_tests", "signed_reversal_score.csv"))

message(sprintf("    Signed reversal: r=%.3f [%.3f, %.3f], p=%.4g",
                cor_res$estimate, cor_res$conf.int[1], cor_res$conf.int[2],
                cor_res$p.value))

message("Reversal tests complete. F3 setup done.")
