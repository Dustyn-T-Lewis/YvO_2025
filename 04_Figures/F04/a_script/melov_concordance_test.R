# F04 Melov-style concordance tests — Option A (cosine) + Option B (Melov direct)
# Standalone exploratory script — run from project root

setwd(rprojroot::find_rstudio_root_file())

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(tidyr)
})

set.seed(42)

# ── Load data ──────────────────────────────────────────────────────────────
dep_df <- read_csv("03_DEP/c_data/03_combined_results.csv", show_col_types = FALSE)

imp_data <- read_csv("02_Imputation/c_data/01_imputed.csv", show_col_types = FALSE)
imp_ann_cols <- c("uniprot_id", "protein", "gene", "description")
imp_samp_cols <- setdiff(names(imp_data), imp_ann_cols)
imp_mat <- as.matrix(imp_data[, imp_samp_cols])
rownames(imp_mat) <- imp_data$uniprot_id

dal_meta <- as.data.frame(
  readRDS("02_Imputation/c_data/01_DAList_imputed.rds")$metadata)
meta <- tibble(
  sample_id = dal_meta$Col_ID,
  subject   = sub("_(Pre|Post)$", "", dal_meta$Col_ID),
  age       = dal_meta$Group,
  time      = dal_meta$Timepoint,
  group     = dal_meta$Group_Time
)

# ── Sample IDs ─────────────────────────────────────────────────────────────
young_pre_ids  <- meta$sample_id[meta$age == "Young" & meta$time == "Pre"]
young_post_ids <- meta$sample_id[meta$age == "Young" & meta$time == "Post"]
old_pre_ids    <- meta$sample_id[meta$age == "Old"   & meta$time == "Pre"]
old_post_ids   <- meta$sample_id[meta$age == "Old"   & meta$time == "Post"]

cat(sprintf("Subjects: %d Young, %d Old\n",
            length(unique(meta$subject[meta$age == "Young"])),
            length(unique(meta$subject[meta$age == "Old"]))))

# ── Define training-responsive signature (Young P < 0.05) ─────────────────
training_sig <- dep_df %>%
  filter(!is.na(P.Value_Training_Young) & P.Value_Training_Young < 0.05) %>%
  pull(uniprot_id)
training_sig <- intersect(training_sig, rownames(imp_mat))
n_sig <- length(training_sig)
cat(sprintf("Training signature: %d proteins (Training_Young nominal P < 0.05)\n\n", n_sig))

# ── Compute group centroids on training signature ──────────────────────────
young_pre_mean  <- rowMeans(imp_mat[training_sig, intersect(young_pre_ids,  colnames(imp_mat)), drop = FALSE], na.rm = TRUE)
young_post_mean <- rowMeans(imp_mat[training_sig, intersect(young_post_ids, colnames(imp_mat)), drop = FALSE], na.rm = TRUE)
old_pre_mean    <- rowMeans(imp_mat[training_sig, intersect(old_pre_ids,    colnames(imp_mat)), drop = FALSE], na.rm = TRUE)
old_post_mean   <- rowMeans(imp_mat[training_sig, intersect(old_post_ids,   colnames(imp_mat)), drop = FALSE], na.rm = TRUE)

# Training response vectors
young_delta <- young_post_mean - young_pre_mean
old_delta   <- old_post_mean   - old_pre_mean

# ══════════════════════════════════════════════════════════════════════════════
# OPTION A: Cosine Similarity Permutation Test
# Question: Is the Old training response directionally aligned with Young's?
# ══════════════════════════════════════════════════════════════════════════════
cosine_sim <- function(a, b) sum(a * b) / (sqrt(sum(a^2)) * sqrt(sum(b^2)))

observed_cosine <- cosine_sim(young_delta, old_delta)
cat("═══ OPTION A: Cosine Similarity Test ═══\n")
cat(sprintf("Observed cosine similarity: %.4f\n", observed_cosine))

# Permutation: shuffle Pre/Post within Old subjects
n_perm <- 10000
set.seed(42)

old_subjects  <- unique(meta$subject[meta$age == "Old"])
old_pre_meta  <- meta %>% filter(age == "Old", time == "Pre")
old_post_meta <- meta %>% filter(age == "Old", time == "Post")

perm_cosines <- numeric(n_perm)
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

  perm_pre_mean  <- rowMeans(imp_mat[training_sig, intersect(perm_pre_ids,  colnames(imp_mat)), drop = FALSE], na.rm = TRUE)
  perm_post_mean <- rowMeans(imp_mat[training_sig, intersect(perm_post_ids, colnames(imp_mat)), drop = FALSE], na.rm = TRUE)
  perm_old_delta <- perm_post_mean - perm_pre_mean
  perm_cosines[i] <- cosine_sim(young_delta, perm_old_delta)
}

cosine_p <- mean(perm_cosines >= observed_cosine)
cosine_p_ci <- binom.test(sum(perm_cosines >= observed_cosine), n_perm)$conf.int

cat(sprintf("Permutation p-value: %.4f [%.4f, %.4f]\n", cosine_p, cosine_p_ci[1], cosine_p_ci[2]))
cat(sprintf("Null cosine mean: %.4f, SD: %.4f\n", mean(perm_cosines), sd(perm_cosines)))
cat(sprintf("Observed is %.1f SDs above null mean\n\n",
            (observed_cosine - mean(perm_cosines)) / sd(perm_cosines)))

# ══════════════════════════════════════════════════════════════════════════════
# OPTION B: Direct Melov Adaptation — Distance to Young_Pre Reference
# Question: Does training in Old move the proteome toward the Young_Pre state?
#           (same as F05, but on the TRAINING signature instead of aging signature)
# ══════════════════════════════════════════════════════════════════════════════
cat("═══ OPTION B: Melov Distance-to-Reference Test ═══\n")

d_old_pre  <- sqrt(sum((old_pre_mean  - young_pre_mean)^2))
d_old_post <- sqrt(sum((old_post_mean - young_pre_mean)^2))
observed_delta_B <- d_old_pre - d_old_post
reversal_pct_B <- (d_old_pre - d_old_post) / d_old_pre * 100

cat(sprintf("d(Old_Pre,  Young_Pre) = %.4f\n", d_old_pre))
cat(sprintf("d(Old_Post, Young_Pre) = %.4f\n", d_old_post))
cat(sprintf("Observed delta (d_pre - d_post) = %.4f\n", observed_delta_B))
cat(sprintf("Magnitude reversal: %.1f%%\n", reversal_pct_B))

set.seed(42)
perm_deltas_B <- numeric(n_perm)

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

  perm_pre_mean  <- rowMeans(imp_mat[training_sig, intersect(perm_pre_ids,  colnames(imp_mat)), drop = FALSE], na.rm = TRUE)
  perm_post_mean <- rowMeans(imp_mat[training_sig, intersect(perm_post_ids, colnames(imp_mat)), drop = FALSE], na.rm = TRUE)
  d_pre_perm  <- sqrt(sum((perm_pre_mean  - young_pre_mean)^2))
  d_post_perm <- sqrt(sum((perm_post_mean - young_pre_mean)^2))
  perm_deltas_B[i] <- d_pre_perm - d_post_perm
}

melov_p_B <- mean(perm_deltas_B >= observed_delta_B)
melov_p_ci_B <- binom.test(sum(perm_deltas_B >= observed_delta_B), n_perm)$conf.int

# Bootstrap CI on reversal %
set.seed(42)
boot_rev_B <- replicate(2000, {
  idx <- sample(seq_along(training_sig), replace = TRUE)
  b_d_pre  <- sqrt(sum((old_pre_mean[idx]  - young_pre_mean[idx])^2))
  b_d_post <- sqrt(sum((old_post_mean[idx] - young_pre_mean[idx])^2))
  (b_d_pre - b_d_post) / b_d_pre * 100
})
rev_ci_B <- quantile(boot_rev_B, c(0.025, 0.975))

cat(sprintf("Permutation p-value: %.4f [%.4f, %.4f]\n", melov_p_B, melov_p_ci_B[1], melov_p_ci_B[2]))
cat(sprintf("Reversal %% = %.1f%% [%.1f, %.1f]\n\n", reversal_pct_B, rev_ci_B[1], rev_ci_B[2]))

# ══════════════════════════════════════════════════════════════════════════════
# BONUS: Also run Option A on the AGING signature (same space as F05)
# for direct comparison with the existing F05 Melov test
# ══════════════════════════════════════════════════════════════════════════════
cat("═══ BONUS: Cosine on AGING signature (for comparison with F05) ═══\n")

aging_sig <- dep_df %>%
  filter(!is.na(P.Value_Aging) & P.Value_Aging < 0.05) %>%
  pull(uniprot_id)
aging_sig <- intersect(aging_sig, rownames(imp_mat))
cat(sprintf("Aging signature: %d proteins\n", length(aging_sig)))

young_pre_aging  <- rowMeans(imp_mat[aging_sig, intersect(young_pre_ids,  colnames(imp_mat)), drop = FALSE], na.rm = TRUE)
young_post_aging <- rowMeans(imp_mat[aging_sig, intersect(young_post_ids, colnames(imp_mat)), drop = FALSE], na.rm = TRUE)
old_pre_aging    <- rowMeans(imp_mat[aging_sig, intersect(old_pre_ids,    colnames(imp_mat)), drop = FALSE], na.rm = TRUE)
old_post_aging   <- rowMeans(imp_mat[aging_sig, intersect(old_post_ids,   colnames(imp_mat)), drop = FALSE], na.rm = TRUE)

# For reversal: aging vector = Old_Pre - Young_Pre; training vector = Old_Post - Old_Pre
aging_vector   <- old_pre_aging - young_pre_aging
old_training_vector <- old_post_aging - old_pre_aging
cos_aging <- cosine_sim(-aging_vector, old_training_vector)  # negate aging = reversal direction
cat(sprintf("Cosine(−Aging, Training_Old) on aging sig: %.4f\n", cos_aging))
cat("  (positive = training opposes aging direction = reversal)\n")

# Young training direction on aging signature
young_training_aging <- young_post_aging - young_pre_aging
cos_young_old_aging <- cosine_sim(young_training_aging, old_training_vector)
cat(sprintf("Cosine(Training_Young, Training_Old) on aging sig: %.4f\n", cos_young_old_aging))
cat("  (positive = concordant training response on aging-affected proteins)\n\n")

# ══════════════════════════════════════════════════════════════════════════════
# Summary table
# ══════════════════════════════════════════════════════════════════════════════
cat("═══ SUMMARY ═══\n")
summary_df <- tibble(
  test = c("Option A: Cosine concordance (training sig)",
           "Option B: Melov distance-to-Young_Pre (training sig)",
           "F05 existing: Melov reversal (aging sig)"),
  signature = c(sprintf("%d proteins (TY P<0.05)", n_sig),
                sprintf("%d proteins (TY P<0.05)", n_sig),
                sprintf("%d proteins (Aging P<0.05)", length(aging_sig))),
  statistic = c(sprintf("cos = %.4f", observed_cosine),
                sprintf("delta = %.4f (%.1f%%)", observed_delta_B, reversal_pct_B),
                "delta = 1.404 (9.8%)"),
  p_value = c(sprintf("%.4f", cosine_p),
              sprintf("%.4f", melov_p_B),
              "0.1455"),
  interpretation = c(
    ifelse(cosine_p < 0.05, "SIGNIFICANT — Old training aligns with Young", "NS — directional alignment not distinguishable from noise"),
    ifelse(melov_p_B < 0.05, "SIGNIFICANT — Old moves toward Young_Pre", "NS — distance reduction not distinguishable from noise"),
    "NS — reversal magnitude not distinguishable from noise"
  )
)

for (r in seq_len(nrow(summary_df))) {
  cat(sprintf("\n%s\n  Signature: %s\n  Statistic: %s\n  p = %s\n  → %s\n",
              summary_df$test[r], summary_df$signature[r],
              summary_df$statistic[r], summary_df$p_value[r],
              summary_df$interpretation[r]))
}

cat("\nDone.\n")
