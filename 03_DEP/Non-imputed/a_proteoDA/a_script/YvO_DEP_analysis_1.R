################################################################################
#
#   YvO Differential Expression of Proteins — proteoDA Pipeline (Non-imputed)
#
#   Pure proteoDA workflow:
#     0.  Setup & paths
#     1.  Load data & build metadata
#     2.  Create DAList
#     3.  Statistical design (repeated measures)
#     4.  Contrasts
#     5.  Fit model & extract results
#     6.  Write tables
#     7.  Write plots
#     8.  Design validation & summary
#
#   Design: 2x2 factorial (Age x Time) with repeated measures on subject
#   Samples: 62 (16 Young_Pre, 16 Young_Post, 15 Old_Pre, 15 Old_Post)
#   Proteins: 2124 (log2-transformed, median-normalized, NOT imputed — NAs present)
#
################################################################################

# ==============================================================================
# 0.  SETUP & PATHS
# ==============================================================================

cat("=== YvO proteoDA Pipeline (Non-imputed) ===\n\n")
cat(">> 0 -- Setup\n")

suppressPackageStartupMessages({
  library(tidyverse)
  library(proteoDA)
  library(openxlsx)
})

# Paths derived from a_script/ working directory (4 levels up to project root)
base_dir   <- normalizePath(file.path(dirname(getwd()), "..", "..", ".."), mustWork = TRUE)
DATA_FILE  <- file.path(base_dir, "01_normalization", "c_data", "01_normalized.csv")
REPORT_DIR <- file.path(base_dir, "03_DEP", "Non-imputed", "a_proteoDA", "b_reports")
DATA_DIR   <- file.path(base_dir, "03_DEP", "Non-imputed", "a_proteoDA", "c_data")

dir.create(REPORT_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(DATA_DIR,   recursive = TRUE, showWarnings = FALSE)

cat("   Base directory:", base_dir, "\n")

# ==============================================================================
# 1.  LOAD DATA & BUILD METADATA
# ==============================================================================

cat("\n>> 1 -- Load Data\n")

df <- read_csv(DATA_FILE, show_col_types = FALSE)
cat(sprintf("   Loaded: %d proteins x %d columns\n", nrow(df), ncol(df)))

# Separate annotation from intensity matrix
ann_cols   <- c("uniprot_id", "protein", "gene", "description")
ann        <- df[, ann_cols]
samp_names <- setdiff(names(df), ann_cols)
mat        <- as.matrix(df[, samp_names])
rownames(mat) <- ann$uniprot_id

n_missing  <- sum(is.na(mat))
pct_missing <- round(100 * n_missing / length(mat), 1)
cat(sprintf("   Matrix: %d proteins x %d samples\n", nrow(mat), ncol(mat)))
cat(sprintf("   Missing values: %d (%.1f%%)\n", n_missing, pct_missing))

# Parse sample IDs: {Prefix}_{SubjNum}_{Time}
# O/OP = Old, Y/YP = Young
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
                     levels = c("Young_Pre", "Young_Post",
                                "Old_Pre",   "Old_Post"))

cat("\n   Sample distribution:\n")
print(table(meta$age, meta$time))
stopifnot(identical(colnames(mat), meta$sample_id))

# ==============================================================================
# 2.  CREATE DAList
# ==============================================================================

cat("\n>> 2 -- Create DAList\n")

meta_df <- as.data.frame(meta)
rownames(meta_df) <- meta$sample_id

ann_df <- as.data.frame(ann)

dal <- DAList(
  data       = mat,
  annotation = ann_df,
  metadata   = meta_df,
  tags       = list(norm_method = "median")
)

cat(sprintf("   DAList: %d proteins x %d samples\n", nrow(dal$data), ncol(dal$data)))
cat("   Tag: norm_method = median (pre-normalized input, NAs retained)\n")

# ==============================================================================
# 3.  STATISTICAL DESIGN
# ==============================================================================

cat("\n>> 3 -- Statistical Design\n")

dal <- add_design(dal, "~ 0 + group + (1 | subject)")

cat("   Design: ~ 0 + group + (1 | subject)\n")
cat(sprintf("   Design matrix: %d samples x %d groups\n",
            nrow(dal$design$design_matrix), ncol(dal$design$design_matrix)))
cat(sprintf("   Random factor: subject (%d levels)\n",
            length(unique(meta$subject))))

# Strip "group" prefix from design matrix column names so contrast names are clean
colnames(dal$design$design_matrix) <- gsub("^group", "", colnames(dal$design$design_matrix))
cat("   Design columns:", paste(colnames(dal$design$design_matrix), collapse = ", "), "\n")

# ==============================================================================
# 4.  CONTRASTS
# ==============================================================================

cat("\n>> 4 -- Contrasts\n")

dal <- add_contrasts(dal, contrasts_vector = c(
  "Training_Young = Young_Post - Young_Pre",
  "Training_Old = Old_Post - Old_Pre",
  "Aging = Old_Pre - Young_Pre",
  "Interaction = (Old_Post - Old_Pre) - (Young_Post - Young_Pre)"
))

cat("   Training_Young = Young_Post - Young_Pre\n")
cat("   Training_Old   = Old_Post - Old_Pre\n")
cat("   Aging           = Old_Pre - Young_Pre\n")
cat("   Interaction     = (Old_Post - Old_Pre) - (Young_Post - Young_Pre)\n")

# ==============================================================================
# 5.  FIT MODEL & EXTRACT RESULTS
# ==============================================================================

cat("\n>> 5 -- Fit Model & Extract Results\n")

dal <- fit_limma_model(dal)

# Report within-subject correlation if available
if (!is.null(dal$eBayes_fit$correlation)) {
  cat(sprintf("   Within-subject correlation: %.3f\n", dal$eBayes_fit$correlation))
} else if (!is.null(dal$tags$duplicate_correlation)) {
  cat(sprintf("   Within-subject correlation: %.3f\n", dal$tags$duplicate_correlation))
}

dal <- extract_DA_results(dal,
                          pval_thresh  = 0.05,
                          lfc_thresh   = 0,
                          adj_method   = "BH")

cat("   Results extracted (FDR < 0.05, no LFC threshold)\n")

# ==============================================================================
# 6.  WRITE TABLES
# ==============================================================================

cat("\n>> 6 -- Write Tables\n")

write_limma_tables(dal,
                   output_dir = DATA_DIR,
                   overwrite  = TRUE)

cat("   Saved: DA_summary.csv, combined_results.csv, results.xlsx\n")
cat("   Saved: per_contrast_results/ (one CSV per contrast)\n")

# --- 6b: Add pi-scores, sig_pi direction & update outputs --------------------
cat("\n>> 6b -- Add Pi-scores, Direction Flags & Update Summary\n")

contrast_names <- names(dal$results)
per_contrast_dir <- file.path(DATA_DIR, "per_contrast_results")

# Per-contrast CSVs: add pi_score + sig_pi (1/-1/0), drop sig.PVal and sig.FDR
for (cname in contrast_names) {
  fpath <- file.path(per_contrast_dir, paste0(cname, ".csv"))
  if (file.exists(fpath)) {
    res <- read_csv(fpath, show_col_types = FALSE)
    res <- res %>%
      mutate(
        pi_score = P.Value ^ abs(logFC),
        sig_pi = case_when(
          pi_score < 0.05 & logFC > 0 ~  1L,
          pi_score < 0.05 & logFC < 0 ~ -1L,
          TRUE ~ 0L
        ),
        sig_nominal_01 = case_when(
          P.Value < 0.01 & logFC > 0 ~  1L,
          P.Value < 0.01 & logFC < 0 ~ -1L,
          TRUE ~ 0L
        )
      ) %>%
      dplyr::select(-any_of(c("sig.PVal", "sig.FDR")))
    write_csv(res, fpath)
  }
}
cat(sprintf("   Added pi_score + sig_pi to %d per-contrast CSVs\n",
            length(contrast_names)))

# Combined results: add pi_score + sig_pi per contrast, drop old sig columns
combined_path <- file.path(DATA_DIR, "combined_results.csv")
if (file.exists(combined_path)) {
  comb <- read_csv(combined_path, show_col_types = FALSE)
  for (cname in contrast_names) {
    lfc_col  <- paste0("logFC_", cname)
    pval_col <- paste0("P.Value_", cname)
    pi_col   <- paste0("pi_score_", cname)
    sig_col  <- paste0("sig_pi_", cname)
    if (all(c(lfc_col, pval_col) %in% names(comb))) {
      comb[[pi_col]]  <- comb[[pval_col]] ^ abs(comb[[lfc_col]])
      comb[[sig_col]] <- case_when(
        comb[[pi_col]] < 0.05 & comb[[lfc_col]] > 0 ~  1L,
        comb[[pi_col]] < 0.05 & comb[[lfc_col]] < 0 ~ -1L,
        TRUE ~ 0L
      )
      sig01_col <- paste0("sig_nominal_01_", cname)
      comb[[sig01_col]] <- case_when(
        comb[[pval_col]] < 0.01 & comb[[lfc_col]] > 0 ~  1L,
        comb[[pval_col]] < 0.01 & comb[[lfc_col]] < 0 ~ -1L,
        TRUE ~ 0L
      )
    }
  }
  drop_cols <- grep("^(sig\\.PVal_|sig\\.FDR_)", names(comb), value = TRUE)
  if (length(drop_cols) > 0) comb <- comb[, !names(comb) %in% drop_cols]
  write_csv(comb, combined_path)
  cat("   Added pi_score + sig_pi to combined_results.csv\n")
}

# DA summary: add sig.Pi counts alongside existing sig.PVal and sig.FDR
da_path <- file.path(DATA_DIR, "DA_summary.csv")
if (file.exists(da_path)) {
  da <- read_csv(da_path, show_col_types = FALSE)
  pi_rows <- list()
  for (cname in contrast_names) {
    fpath <- file.path(per_contrast_dir, paste0(cname, ".csv"))
    if (file.exists(fpath)) {
      res <- read_csv(fpath, show_col_types = FALSE)
      n_up   <- sum(res$sig_pi ==  1, na.rm = TRUE)
      n_down <- sum(res$sig_pi == -1, na.rm = TRUE)
      n_ns   <- sum(res$sig_pi ==  0, na.rm = TRUE)
      n_up_nom   <- sum(res$sig_nominal_01 ==  1, na.rm = TRUE)
      n_down_nom <- sum(res$sig_nominal_01 == -1, na.rm = TRUE)
      n_ns_nom   <- sum(res$sig_nominal_01 ==  0, na.rm = TRUE)
      pi_rows <- c(pi_rows, list(
        tibble(contrast = cname, type = "up",     sig.Pi = n_up,   sig.Nominal.01 = n_up_nom),
        tibble(contrast = cname, type = "down",   sig.Pi = n_down, sig.Nominal.01 = n_down_nom),
        tibble(contrast = cname, type = "nonsig", sig.Pi = n_ns,   sig.Nominal.01 = n_ns_nom)
      ))
    }
  }
  pi_df <- bind_rows(pi_rows)
  da <- da %>% left_join(pi_df, by = c("contrast", "type"))
  write_csv(da, da_path)
  cat("   Added sig.Pi counts to DA_summary.csv\n")
}

# Rebuild results.xlsx with pi_score + sig_pi (one sheet per contrast + summary)
xlsx_path <- file.path(DATA_DIR, "results.xlsx")
wb <- createWorkbook()
for (cname in contrast_names) {
  fpath <- file.path(per_contrast_dir, paste0(cname, ".csv"))
  if (file.exists(fpath)) {
    res <- read_csv(fpath, show_col_types = FALSE)
    addWorksheet(wb, cname)
    writeData(wb, cname, res)
  }
}
if (file.exists(da_path)) {
  da <- read_csv(da_path, show_col_types = FALSE)
  addWorksheet(wb, "DA_Summary")
  writeData(wb, "DA_Summary", da)
}
saveWorkbook(wb, xlsx_path, overwrite = TRUE)
cat("   Rebuilt results.xlsx with pi_score + sig_pi\n")

# ==============================================================================
# 7.  WRITE PLOTS
# ==============================================================================

cat("\n>> 7 -- Write Plots\n")

write_limma_plots(dal,
                  grouping_column = "group",
                  output_dir      = REPORT_DIR,
                  table_columns   = c("uniprot_id", "gene", "protein"),
                  title_column    = "gene",
                  overwrite       = TRUE)

# Reorganize outputs into per-contrast subdirectories
contrast_names <- names(dal$results)
static_dir <- file.path(REPORT_DIR, "static_plots")

for (cname in contrast_names) {
  contrast_dir <- file.path(REPORT_DIR, cname)
  dir.create(contrast_dir, showWarnings = FALSE)

  # Move HTML report into contrast subdirectory
  html_file <- file.path(REPORT_DIR, paste0(cname, "_DA_report.html"))
  if (file.exists(html_file)) {
    file.rename(html_file, file.path(contrast_dir, basename(html_file)))
  }

  # Move static PDFs into contrast subdirectory
  pdfs <- list.files(static_dir, pattern = paste0("^", cname, "-"), full.names = TRUE)
  for (pdf in pdfs) {
    # Strip contrast prefix from filename: "Training_Young-volcano-raw-pval.pdf" -> "volcano-raw-pval.pdf"
    clean_name <- sub(paste0("^", cname, "-"), "", basename(pdf))
    file.rename(pdf, file.path(contrast_dir, clean_name))
  }

  cat(sprintf("   %s/: HTML report + %d static PDFs\n", cname, length(pdfs)))
}

# Remove empty static_plots/ directory
if (length(list.files(static_dir)) == 0) {
  unlink(static_dir, recursive = TRUE)
}

# ==============================================================================
# 8.  DESIGN VALIDATION & SUMMARY
# ==============================================================================

cat("\n>> 8 -- Summary\n")

cat("\n   Design matrix dimensions:", paste(dim(dal$design$design_matrix), collapse = " x "), "\n")
cat("   Design matrix columns:", paste(colnames(dal$design$design_matrix), collapse = ", "), "\n")

cat("\n   Contrast matrix:\n")
print(dal$design$contrast_matrix)

cat("\n   DA summary (FDR < 0.05):\n")
da_summary <- read_csv(file.path(DATA_DIR, "DA_summary.csv"), show_col_types = FALSE)
print(da_summary)

cat("\n")
sessionInfo()

cat("\n=== YvO proteoDA Pipeline (Non-imputed) Complete ===\n")
