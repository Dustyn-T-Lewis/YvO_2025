################################################################################
#
#   YvO Differential Expression of Proteins — DEqMS Pipeline (Non-imputed)
#
#   Replaces proteoDA with limma + DEqMS:
#     - limma handles design, contrasts, duplicateCorrelation (repeated measures)
#     - DEqMS replaces eBayes() with spectraCounteBayes() to estimate
#       peptide-count-dependent prior variance per protein
#
#   Workflow:
#     0.  Setup & paths
#     1.  Load data & merge n_seq from raw xlsx
#     2.  Build design matrix & contrasts (manual limma, no proteoDA)
#     3.  Estimate within-subject correlation (duplicateCorrelation)
#     4.  Fit linear model (lmFit → contrasts.fit)
#     5.  DEqMS: spectraCounteBayes (replaces eBayes)
#     6.  Extract results & write tables
#     7.  Diagnostic plots (variance, volcanos, t-stat comparison, p-value histograms)
#     8.  Summary
#
#   Design: 2×2 factorial (Age × Time) with repeated measures on subject
#   Samples: 62 (16 Young_Pre, 16 Young_Post, 15 Old_Pre, 15 Old_Post)
#   Proteins: 2124 (log2-transformed, median-normalized, NOT imputed — NAs present)
#   DE method: DEqMS (peptide-count-weighted empirical Bayes)
#
################################################################################

# ==============================================================================
# 0.  SETUP & PATHS
# ==============================================================================

cat("=== YvO DEqMS Pipeline (Non-imputed) ===\n\n")
cat(">> 0 — Setup\n")

suppressPackageStartupMessages({
  library(tidyverse)
  library(limma)
  library(DEqMS)
  library(readxl)
  library(ggrepel)
  library(patchwork)
})

# Paths derived from a_script/ working directory (3 levels up to project root)
base_dir   <- normalizePath(file.path(dirname(getwd()), "..", ".."), mustWork = TRUE)
DATA_FILE  <- file.path(base_dir, "01_normalization", "c_data", "01_normalized.csv")
RAW_FILE   <- file.path(base_dir, "00_input", "YvO_raw.xlsx")
REPORT_DIR <- file.path(base_dir, "03_DEP", "DEqMS_Non-imputed", "b_reports")
DATA_DIR   <- file.path(base_dir, "03_DEP", "DEqMS_Non-imputed", "c_data")

dir.create(REPORT_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(DATA_DIR,   recursive = TRUE, showWarnings = FALSE)

cat("   Base directory:", base_dir, "\n")

# ==============================================================================
# 1.  LOAD DATA & MERGE n_seq
# ==============================================================================

cat("\n>> 1 — Load Data & Merge n_seq\n")

df <- read_csv(DATA_FILE, show_col_types = FALSE)
cat(sprintf("   Loaded: %d proteins x %d columns\n", nrow(df), ncol(df)))

# Separate annotation from intensity matrix (4 annotation cols, no n_seq)
ann_cols   <- c("uniprot_id", "protein", "gene", "description")
ann        <- df[, ann_cols]
samp_names <- setdiff(names(df), ann_cols)
mat        <- as.matrix(df[, samp_names])
rownames(mat) <- ann$uniprot_id

n_missing   <- sum(is.na(mat))
pct_missing <- round(100 * n_missing / length(mat), 1)
cat(sprintf("   Matrix: %d proteins x %d samples\n", nrow(mat), ncol(mat)))
cat(sprintf("   Missing values: %d (%.1f%%)\n", n_missing, pct_missing))

# Read n_seq from raw xlsx (2837 proteins; all 2124 normalized are a subset)
raw_ann <- read_xlsx(RAW_FILE) %>%
  dplyr::select(uniprot_id, n_seq) %>%
  distinct(uniprot_id, .keep_all = TRUE)
cat(sprintf("   Raw data: %d proteins with n_seq\n", nrow(raw_ann)))

# Left-join n_seq onto annotation
ann <- ann %>%
  left_join(raw_ann, by = "uniprot_id")

n_matched <- sum(!is.na(ann$n_seq))
n_total   <- nrow(ann)
cat(sprintf("   n_seq matched: %d / %d proteins\n", n_matched, n_total))
stopifnot(n_matched == n_total)  # all proteins must have n_seq

# Peptide count vector (aligned to protein rows) — this is what DEqMS needs
pep_count <- ann$n_seq
names(pep_count) <- ann$uniprot_id

cat(sprintf("   Peptide counts: range %d–%d, median %d\n",
            min(pep_count, na.rm = TRUE),
            max(pep_count, na.rm = TRUE),
            median(pep_count, na.rm = TRUE)))

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
# 2.  DESIGN MATRIX & CONTRASTS (manual limma — replaces proteoDA)
# ==============================================================================

cat("\n>> 2 — Design Matrix & Contrasts\n")

# Cell-means model (no intercept) — same structure as proteoDA version
design <- model.matrix(~ 0 + group, data = meta)
colnames(design) <- gsub("^group", "", colnames(design))

cat(sprintf("   Design matrix: %d samples x %d groups\n",
            nrow(design), ncol(design)))
cat("   Design columns:", paste(colnames(design), collapse = ", "), "\n")

# Same contrasts as proteoDA version
contrast_mat <- makeContrasts(
  Training_Young = Young_Post - Young_Pre,
  Training_Old   = Old_Post   - Old_Pre,
  Aging          = Old_Pre    - Young_Pre,
  Interaction    = (Old_Post - Old_Pre) - (Young_Post - Young_Pre),
  levels = design
)

cat("\n   Contrasts defined:\n")
cat("   Training_Young = Young_Post - Young_Pre\n")
cat("   Training_Old   = Old_Post - Old_Pre\n")
cat("   Aging          = Old_Pre - Young_Pre\n")
cat("   Interaction    = (Old_Post - Old_Pre) - (Young_Post - Young_Pre)\n")

# ==============================================================================
# 3.  WITHIN-SUBJECT CORRELATION (duplicateCorrelation)
# ==============================================================================

cat("\n>> 3 — Within-Subject Correlation\n")

block <- meta$subject

# duplicateCorrelation estimates consensus within-subject correlation
# This is the limma equivalent of proteoDA's (1 | subject) random effect
dupcor <- duplicateCorrelation(mat, design = design, block = block)

cat(sprintf("   Within-subject correlation: %.3f\n", dupcor$consensus.correlation))

# ==============================================================================
# 4.  FIT LINEAR MODEL
# ==============================================================================

cat("\n>> 4 — Fit Linear Model\n")

# lmFit with blocking (repeated measures) — equivalent to proteoDA's fit_limma_model
fit <- lmFit(mat,
             design      = design,
             block       = block,
             correlation = dupcor$consensus.correlation)

# Apply contrasts
fit2 <- contrasts.fit(fit, contrast_mat)

cat(sprintf("   Model fit: %d proteins x %d contrasts\n",
            nrow(fit2$coefficients), ncol(fit2$coefficients)))

# ==============================================================================
# 5.  DEqMS: spectraCounteBayes (replaces eBayes)
# ==============================================================================

cat("\n>> 5 — DEqMS: Peptide-Count-Weighted Empirical Bayes\n")

# Standard eBayes first (DEqMS builds on top of it)
fit2 <- eBayes(fit2)

# Attach peptide count to the fit object — DEqMS requires this
fit2$count <- pep_count[rownames(fit2$coefficients)]

# Replace standard eBayes variance with peptide-count-aware variance
# This is the core DEqMS step: fits a smooth variance ~ count curve
# so each protein gets a count-appropriate prior variance
fit2 <- spectraCounteBayes(fit2)

cat("   spectraCounteBayes complete\n")
cat(sprintf("   Limma prior variance (flat):  %.4f\n", fit2$s2.prior))
cat("   DEqMS prior variance: peptide-count-dependent (see diagnostic plot)\n")

# ==============================================================================
# 6.  EXTRACT RESULTS & WRITE TABLES
# ==============================================================================

cat("\n>> 6 — Extract Results & Write Tables\n")

contrast_names <- colnames(contrast_mat)

# Extract DEqMS results for each contrast
all_results <- list()

for (i in seq_along(contrast_names)) {
  cname <- contrast_names[i]

  # DEqMS::outputResult returns a data.frame with:
  #   logFC, count, sca.t (DEqMS t-stat), sca.P.Value, sca.adj.pval
  #   plus the standard limma columns (t, P.Value, adj.P.Val) for comparison
  res <- outputResult(fit2, coef_col = i)

  # Add gene annotation (drop DEqMS's 'gene' column which is just the rowname)
  res$gene <- NULL
  res$uniprot_id <- rownames(res)
  res <- res %>%
    left_join(as.data.frame(ann) %>% dplyr::select(uniprot_id, gene, protein, description),
              by = "uniprot_id") %>%
    arrange(sca.P.Value)

  all_results[[cname]] <- res

  # Write per-contrast CSV
  contrast_dir <- file.path(DATA_DIR, "per_contrast_results")
  dir.create(contrast_dir, showWarnings = FALSE, recursive = TRUE)
  write_csv(res, file.path(contrast_dir, paste0(cname, "_DEqMS.csv")))

  # Summary counts
  n_sig_deqms <- sum(res$sca.adj.pval < 0.05, na.rm = TRUE)
  n_sig_limma <- sum(res$adj.P.Val    < 0.05, na.rm = TRUE)
  cat(sprintf("   %s: %d DEqMS sig / %d limma sig (FDR < 0.05)\n",
              cname, n_sig_deqms, n_sig_limma))
}

# Combined results table (all contrasts)
combined <- bind_rows(
  lapply(names(all_results), function(cname) {
    all_results[[cname]] %>%
      mutate(contrast = cname, .before = 1)
  })
)
write_csv(combined, file.path(DATA_DIR, "combined_results_DEqMS.csv"))

# DA summary table (matches proteoDA output format, with both DEqMS and limma counts)
da_summary <- tibble(
  contrast     = contrast_names,
  n_tested     = sapply(all_results, nrow),
  n_sig_deqms  = sapply(all_results, function(r) sum(r$sca.adj.pval < 0.05, na.rm = TRUE)),
  n_up_deqms   = sapply(all_results, function(r) sum(r$sca.adj.pval < 0.05 & r$logFC > 0, na.rm = TRUE)),
  n_down_deqms = sapply(all_results, function(r) sum(r$sca.adj.pval < 0.05 & r$logFC < 0, na.rm = TRUE)),
  n_sig_limma  = sapply(all_results, function(r) sum(r$adj.P.Val < 0.05, na.rm = TRUE)),
  n_up_limma   = sapply(all_results, function(r) sum(r$adj.P.Val < 0.05 & r$logFC > 0, na.rm = TRUE)),
  n_down_limma = sapply(all_results, function(r) sum(r$adj.P.Val < 0.05 & r$logFC < 0, na.rm = TRUE))
)
write_csv(da_summary, file.path(DATA_DIR, "DA_summary_DEqMS.csv"))

cat("\n   DA Summary (FDR < 0.05):\n")
print(da_summary)

# Save fit object and results list for downstream use
saveRDS(fit2, file.path(DATA_DIR, "DEqMS_fit.rds"))
saveRDS(all_results, file.path(DATA_DIR, "DEqMS_results_list.rds"))

# ==============================================================================
# 7.  DIAGNOSTIC PLOTS
# ==============================================================================

cat("\n>> 7 — Diagnostic Plots\n")

# --- 7a: DEqMS variance diagnostic -------------------------------------------
# Core DEqMS plot: variance vs peptide count showing the fitted curve
# Red = DEqMS count-dependent prior; Green = limma flat prior

pdf(file.path(REPORT_DIR, "07a_DEqMS_variance_diagnostic.pdf"), width = 10, height = 8)

par(mfrow = c(1, 2))

# Scatter: variance vs log2(peptide count)
VarianceScatterplot(fit2, xlab = "log2(Peptide count)")
abline(h = log(fit2$s2.prior), col = "green", lwd = 2)
legend("topright",
       legend = c("DEqMS prior variance", "Limma prior variance"),
       col = c("red", "green"), lwd = 2, cex = 0.8)

# Box: variance by peptide count bin
VarianceBoxplot(fit2, n = 20, xlab = "Peptide count", main = "Variance by peptide count")

dev.off()
cat("   Saved: 07a_DEqMS_variance_diagnostic.pdf\n")

# --- 7b: Volcano plots (DEqMS p-values) per contrast -------------------------
# Saved into per-contrast subdirectories to match proteoDA structure

for (cname in contrast_names) {
  res <- all_results[[cname]]

  # Create per-contrast report subdirectory
  contrast_report_dir <- file.path(REPORT_DIR, cname)
  dir.create(contrast_report_dir, showWarnings = FALSE, recursive = TRUE)

  # Use DEqMS-specific p-values (sca.adj.pval)
  res <- res %>%
    mutate(
      neg_log10_p = -log10(sca.P.Value),
      sig = case_when(
        sca.adj.pval < 0.05 & logFC >  0 ~ "Up",
        sca.adj.pval < 0.05 & logFC <  0 ~ "Down",
        TRUE ~ "NS"
      ),
      label = if_else(sca.adj.pval < 0.05, gene, NA_character_)
    )

  n_up   <- sum(res$sig == "Up",   na.rm = TRUE)
  n_down <- sum(res$sig == "Down", na.rm = TRUE)

  p <- ggplot(res, aes(x = logFC, y = neg_log10_p, color = sig)) +
    geom_point(alpha = 0.5, size = 1.5) +
    geom_text_repel(aes(label = label), size = 2.2, max.overlaps = 20,
                    show.legend = FALSE, na.rm = TRUE) +
    scale_color_manual(
      values = c(Up = "#B2182B", Down = "#2166AC", NS = "gray70"),
      labels = c(
        Up   = sprintf("Up (%d)", n_up),
        Down = sprintf("Down (%d)", n_down),
        NS   = "NS"
      )
    ) +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed", alpha = 0.4) +
    labs(
      x = "log2 Fold Change",
      y = "-log10(DEqMS p-value)",
      title = sprintf("%s (DEqMS)", cname),
      color = NULL
    ) +
    theme_minimal(base_size = 12) +
    theme(legend.position = "bottom")

  ggsave(file.path(contrast_report_dir, "volcano_DEqMS.pdf"),
         p, width = 8, height = 6)
  cat(sprintf("   Saved: %s/volcano_DEqMS.pdf\n", cname))
}

# --- 7c: DEqMS vs limma t-stat comparison (4-panel) --------------------------
# Scatter: DEqMS t-stat vs limma t-stat to show impact of count weighting

p_list <- list()
for (cname in contrast_names) {
  res <- all_results[[cname]]

  p_list[[cname]] <- ggplot(res, aes(x = t, y = sca.t, color = log2(count + 1))) +
    geom_point(alpha = 0.5, size = 1) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "red") +
    scale_color_viridis_c(name = "log2(n_seq)") +
    labs(
      x = "limma t-statistic",
      y = "DEqMS t-statistic (sca.t)",
      title = cname
    ) +
    theme_minimal(base_size = 10)
}

p_compare <- wrap_plots(p_list, ncol = 2) +
  plot_annotation(
    title = "DEqMS vs limma t-statistics",
    subtitle = "Points colored by log2(peptide count); dashed = identity line"
  )

ggsave(file.path(REPORT_DIR, "07c_DEqMS_vs_limma_tstat.pdf"),
       p_compare, width = 12, height = 10)
cat("   Saved: 07c_DEqMS_vs_limma_tstat.pdf\n")

# --- 7d: P-value histograms per contrast (DEqMS) -----------------------------
# Saved into per-contrast subdirectories to match proteoDA structure

for (cname in contrast_names) {
  res <- all_results[[cname]]

  contrast_report_dir <- file.path(REPORT_DIR, cname)
  dir.create(contrast_report_dir, showWarnings = FALSE, recursive = TRUE)

  p <- ggplot(res, aes(x = sca.P.Value)) +
    geom_histogram(bins = 50, fill = "steelblue", color = "white", boundary = 0) +
    geom_vline(xintercept = 0.05, linetype = "dashed", color = "red", alpha = 0.6) +
    labs(
      x = "DEqMS p-value (sca.P.Value)",
      y = "Count",
      title = sprintf("%s — P-value Distribution (DEqMS)", cname)
    ) +
    theme_minimal(base_size = 12)

  ggsave(file.path(contrast_report_dir, "pval_histogram_DEqMS.pdf"),
         p, width = 7, height = 5)
  cat(sprintf("   Saved: %s/pval_histogram_DEqMS.pdf\n", cname))
}

# ==============================================================================
# 8.  SUMMARY
# ==============================================================================

cat("\n>> 8 — Summary\n")

cat("\n   Design matrix dimensions:", paste(dim(design), collapse = " x "), "\n")
cat("   Design matrix columns:", paste(colnames(design), collapse = ", "), "\n")
cat(sprintf("   Within-subject correlation: %.3f\n", dupcor$consensus.correlation))
cat("   DE method: DEqMS (spectraCounteBayes)\n")
cat(sprintf("   Peptide count covariate: n_seq (range %d–%d, median %d)\n",
            min(pep_count, na.rm = TRUE),
            max(pep_count, na.rm = TRUE),
            median(pep_count, na.rm = TRUE)))

cat("\n   Contrast matrix:\n")
print(contrast_mat)

cat("\n   DA summary (FDR < 0.05):\n")
print(da_summary)

cat("\n   Output files:\n")
cat("   - DA_summary_DEqMS.csv\n")
cat("   - combined_results_DEqMS.csv\n")
cat("   - per_contrast_results/ (one CSV per contrast)\n")
cat("   - DEqMS_fit.rds (full fit object for downstream)\n")
cat("   - DEqMS_results_list.rds (named list of result tables)\n")

cat("\n")
sessionInfo()

cat("\n=== YvO DEqMS Pipeline (Non-imputed) Complete ===\n")
