# F10 prepare_data.R — Baseline Proteome Predicts Training Response
# Loads ME, LV, phenotype data. Computes age gap and within-group correlations.
# Run BEFORE panel scripts.

setwd(rprojroot::find_rstudio_root_file())
source("04_Figures/F10/a_script/style.R")

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(tibble)
  library(tidyr)
})

set.seed(42)

DAT <- "04_Figures/F10/c_data"
RPT <- "04_Figures/F10/b_reports"
dir.create(DAT, recursive = TRUE, showWarnings = FALSE)
dir.create(RPT, recursive = TRUE, showWarnings = FALSE)

# --- Input file guards ---
wgcna_files <- c(
  "04_Figures/WGCNA_F08/c_data/me_pre.rds",
  "04_Figures/WGCNA_F08/c_data/pheno_wide.csv",
  "04_Figures/WGCNA_F08/c_data/subj_age.csv"
)
plier_files <- "04_Figures/PLIER_F09/c_data/02_lv_scores.csv"
meta_file   <- "02_Imputation/c_data/01_DAList_imputed.rds"
all_inputs  <- c(wgcna_files, plier_files, meta_file)
missing     <- all_inputs[!file.exists(all_inputs)]
if (length(missing) > 0) {
  stop("Missing input files for F10 prepare_data.R:\n  ",
       paste(missing, collapse = "\n  "),
       "\nRun WGCNA_F08 and PLIER_F09 first.", call. = FALSE)
}

# --- Load WGCNA data from F6 ---
me_pre    <- readRDS("04_Figures/WGCNA_F08/c_data/me_pre.rds")
pheno     <- read_csv("04_Figures/WGCNA_F08/c_data/pheno_wide.csv", show_col_types = FALSE)
subj_age  <- read_csv("04_Figures/WGCNA_F08/c_data/subj_age.csv", show_col_types = FALSE)

# --- Load PLIER LV scores ---
lv_scores <- read_csv("04_Figures/PLIER_F09/c_data/02_lv_scores.csv", show_col_types = FALSE)

# --- Load metadata for subject pairing ---
dal <- readRDS("02_Imputation/c_data/01_DAList_imputed.rds")
dal_meta <- as.data.frame(dal$metadata)
meta <- tibble(
  sample_id   = dal_meta$Col_ID,
  age         = dal_meta$Group,
  time        = dal_meta$Timepoint,
  subject_key = sub("_(Pre|Post)$", "", dal_meta$Col_ID)
)

# --- Build LV baseline matrix (Pre samples only) ---
lv_pre <- lv_scores %>%
  filter(sample_id %in% meta$sample_id[meta$time == "Pre"]) %>%
  mutate(subject_key = sub("_(Pre|Post)$", "", sample_id)) %>%
  filter(subject_key %in% rownames(me_pre))

lv_mat <- as.matrix(lv_pre[, grep("^LV", names(lv_pre))])
rownames(lv_mat) <- lv_pre$subject_key

common_subj <- intersect(rownames(me_pre), rownames(lv_mat))
common_subj <- intersect(common_subj, pheno$subject_key)
me_pre_use  <- me_pre[common_subj, , drop = FALSE]
lv_mat_use  <- lv_mat[common_subj, , drop = FALSE]
pheno_use   <- pheno %>% filter(subject_key %in% common_subj)
age_use     <- subj_age %>% filter(subject_key %in% common_subj)

# Align phenotype to common_subj order
pheno_use <- pheno_use[match(common_subj, pheno_use$subject_key), ]
age_use   <- age_use[match(common_subj, age_use$subject_key), ]

message(sprintf("Common subjects: %d (%d Young, %d Old)",
                length(common_subj),
                sum(age_use$age == "Young"),
                sum(age_use$age == "Old")))

# --- 1. Age Gap: Mahalanobis distance from Young centroid ---
# Use non-grey MEs
non_grey <- grep("^MEgrey$", colnames(me_pre_use), invert = TRUE)
me_agegap <- me_pre_use[, non_grey, drop = FALSE]

young_idx <- which(age_use$age == "Young")
young_me  <- me_agegap[young_idx, , drop = FALSE]
young_cov <- cov(young_me)
young_center <- colMeans(young_me)

maha_dist <- mahalanobis(me_agegap, center = young_center, cov = young_cov)
eucl_dist <- sqrt(rowSums(sweep(me_agegap, 2, young_center)^2))

age_gap_df <- tibble(
  subject_key = common_subj,
  age         = age_use$age,
  mahalanobis = maha_dist,
  euclidean   = eucl_dist,
  delta_VL    = pheno_use$delta_VL,
  delta_LBM   = pheno_use$delta_LBM
)

# Age gap correlations
agegap_cors <- bind_rows(
  {
    ct <- cor.test(age_gap_df$mahalanobis, age_gap_df$delta_VL)
    ci <- fisher_z_ci(ct$estimate, nrow(age_gap_df))
    tibble(outcome = "delta_VL", distance = "Mahalanobis",
           r = ct$estimate, p = ct$p.value, ci_lo = ci[1], ci_hi = ci[2],
           n = nrow(age_gap_df))
  },
  {
    ct <- cor.test(age_gap_df$mahalanobis, age_gap_df$delta_LBM)
    ci <- fisher_z_ci(ct$estimate, nrow(age_gap_df))
    tibble(outcome = "delta_LBM", distance = "Mahalanobis",
           r = ct$estimate, p = ct$p.value, ci_lo = ci[1], ci_hi = ci[2],
           n = nrow(age_gap_df))
  },
  {
    ct <- cor.test(age_gap_df$euclidean, age_gap_df$delta_VL)
    ci <- fisher_z_ci(ct$estimate, nrow(age_gap_df))
    tibble(outcome = "delta_VL", distance = "Euclidean",
           r = ct$estimate, p = ct$p.value, ci_lo = ci[1], ci_hi = ci[2],
           n = nrow(age_gap_df))
  },
  {
    ct <- cor.test(age_gap_df$euclidean, age_gap_df$delta_LBM)
    ci <- fisher_z_ci(ct$estimate, nrow(age_gap_df))
    tibble(outcome = "delta_LBM", distance = "Euclidean",
           r = ct$estimate, p = ct$p.value, ci_lo = ci[1], ci_hi = ci[2],
           n = nrow(age_gap_df))
  }
)

message("Age gap correlations:")
print(agegap_cors)

# --- 2. Within-group ME correlations ---
# Exclude grey (unassigned proteins)
me_cols <- grep("^MEgrey$", colnames(me_pre_use), invert = TRUE, value = TRUE)

within_me <- expand_grid(
  me = me_cols,
  outcome = c("delta_VL", "delta_LBM"),
  group = c("Young", "Old")
) %>%
  rowwise() %>%
  mutate({
    idx <- which(age_use$age == group)
    if (length(idx) < 5) {
      tibble(r = NA_real_, p = NA_real_, ci_lo = NA_real_, ci_hi = NA_real_, n = length(idx))
    } else {
      x_vals <- me_pre_use[idx, me]
      y_vals <- pheno_use[[outcome]][idx]
      ct <- cor.test(x_vals, y_vals)
      ci <- fisher_z_ci(ct$estimate, length(idx))
      tibble(r = ct$estimate, p = ct$p.value, ci_lo = ci[1], ci_hi = ci[2], n = length(idx))
    }
  }) %>%
  ungroup() %>%
  mutate(r = as.numeric(r), p = as.numeric(p),
         ci_lo = as.numeric(ci_lo), ci_hi = as.numeric(ci_hi),
         n = as.integer(n))

# Select top ME per group x outcome (by |r|)
top_me <- within_me %>%
  group_by(group, outcome) %>%
  slice_max(abs(r), n = 1, with_ties = FALSE) %>%
  ungroup()

message("Top within-group ME predictors:")
print(top_me)

# --- 3. Within-group LV correlations ---
within_lv <- expand_grid(
  lv = colnames(lv_mat_use),
  outcome = c("delta_VL", "delta_LBM"),
  group = c("Young", "Old")
) %>%
  rowwise() %>%
  mutate({
    idx <- which(age_use$age == group)
    if (length(idx) < 5) {
      tibble(r = NA_real_, p = NA_real_, ci_lo = NA_real_, ci_hi = NA_real_, n = length(idx))
    } else {
      x_vals <- lv_mat_use[idx, lv]
      y_vals <- pheno_use[[outcome]][idx]
      ct <- cor.test(x_vals, y_vals)
      ci <- fisher_z_ci(ct$estimate, length(idx))
      tibble(r = ct$estimate, p = ct$p.value, ci_lo = ci[1], ci_hi = ci[2], n = length(idx))
    }
  }) %>%
  ungroup() %>%
  mutate(r = as.numeric(r), p = as.numeric(p),
         ci_lo = as.numeric(ci_lo), ci_hi = as.numeric(ci_hi),
         n = as.integer(n))

# Select top LV per group x outcome (by |r|)
top_lv <- within_lv %>%
  group_by(group, outcome) %>%
  slice_max(abs(r), n = 1, with_ties = FALSE) %>%
  ungroup()

message("Top within-group LV predictors:")
print(top_lv)

# --- 4. Method comparison summary ---
# Pooled age gap
method_rows <- bind_rows(
  agegap_cors %>%
    filter(distance == "Mahalanobis") %>%
    mutate(method = paste0("Age Gap (", outcome, ")"),
           type = "nominal", predictor = "Mahalanobis ME dist.") %>%
    dplyr::select(method, predictor, r, p, ci_lo, ci_hi, n, type),

  top_me %>%
    mutate(method = paste0("Within-", group, " ME (", outcome, ")"),
           type = "nominal (within-group)",
           predictor = me) %>%
    dplyr::select(method, predictor, r, p, ci_lo, ci_hi, n, type),

  top_lv %>%
    mutate(method = paste0("Within-", group, " LV (", outcome, ")"),
           type = "nominal (within-group)",
           predictor = lv) %>%
    dplyr::select(method, predictor, r, p, ci_lo, ci_hi, n, type)
)

message("Method comparison:")
print(method_rows)

# --- Save intermediates ---
write_csv(age_gap_df,   file.path(DAT, "01_age_gap.csv"))
write_csv(agegap_cors,  file.path(DAT, "01_age_gap_cors.csv"))
write_csv(within_me,    file.path(DAT, "02_within_group_me.csv"))
write_csv(top_me,       file.path(DAT, "02_top_me.csv"))
write_csv(within_lv,    file.path(DAT, "03_within_group_lv.csv"))
write_csv(top_lv,       file.path(DAT, "03_top_lv.csv"))
write_csv(method_rows,  file.path(DAT, "04_method_comparison.csv"))

# Save objects for panel scripts
saveRDS(me_pre_use,  file.path(DAT, "me_pre.rds"))
saveRDS(lv_mat_use,  file.path(DAT, "lv_mat.rds"))
saveRDS(pheno_use,   file.path(DAT, "pheno.rds"))
saveRDS(age_use,     file.path(DAT, "age.rds"))

message("F10 prepare_data.R complete")
