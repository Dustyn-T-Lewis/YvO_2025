# --- YvO Differential Expression — limma pipeline (proteoDA) -----------------
# Design: 2x2 factorial (Age x Time), repeated measures on subject
# Input:  cycloess-normalized, non-imputed (limma handles NAs per-protein)
#
# References:
#   Ritchie et al. 2015, Nucleic Acids Res 43(7):e47 — limma
#   Smyth 2005, Stat Appl Genet Mol Biol 3(1):3 — duplicateCorrelation
#   Karpievitch et al. 2012, BMC Bioinform 13(S16):S5 — non-imputed limma
#   Li, Cobbold & Smyth 2025, Bioinformatics — limpa (future direction)
#   Xiao et al. 2014, Bioinformatics 30(6):801-807 (PMID 24395753) — Pi-score
#     Pi = p^|logFC|; threshold Pi < 0.05 <-> original pi > 1.3
# ---------------------------------------------------------------------------

# --- SETUP ---

library(dplyr)
library(tidyr)
library(tibble)
library(readr)
library(purrr)
library(proteoDA)
library(openxlsx)

set.seed(42)

setwd(rprojroot::find_rstudio_root_file())

cfg <- list(
  norm_csv    = "01_normalization/c_data/02_normalized.csv",
  norm_rds    = "01_normalization/c_data/03_DAList_normalized.rds",
  data_dir    = "03_DEP/c_data",
  report_dir  = "03_DEP/b_reports",
  per_dir     = "03_DEP/c_data/04_per_contrast_results",
  proteoDA_dir = "03_DEP/b_reports/01_proteoDA",
  pval_thresh = 0.10,
  lfc_thresh  = 0,
  adj_method  = "BH",
  pi_thresh   = 0.05
)

dir.create(cfg$data_dir,     recursive = TRUE, showWarnings = FALSE)
dir.create(cfg$per_dir,      recursive = TRUE, showWarnings = FALSE)
dir.create(cfg$proteoDA_dir, recursive = TRUE, showWarnings = FALSE)

# --- LOAD DATA & BUILD METADATA ---

df <- read_csv(cfg$norm_csv, show_col_types = FALSE)

ann_cols   <- c("uniprot_id", "protein", "gene", "description")
ann        <- df[, ann_cols]
samp_names <- setdiff(names(df), ann_cols)
mat        <- as.matrix(df[, samp_names])
rownames(mat) <- ann$uniprot_id

cat(sprintf("Loaded: %d proteins x %d samples | missing: %d (%.1f%%)\n",
            nrow(mat), ncol(mat), sum(is.na(mat)),
            100 * sum(is.na(mat)) / length(mat)))

# Canonical metadata from normalisation DAList (not regex-derived)
dal_norm <- readRDS(cfg$norm_rds)
dal_meta <- as.data.frame(dal_norm$metadata)
meta <- tibble(
  sample_id = dal_meta$Col_ID,
  age       = dal_meta$Group,
  time      = dal_meta$Timepoint,
  group     = dal_meta$Group_Time,
  # Col_ID prefix as subject key (Subject_ID is NOT globally unique —
  # 8 IDs are shared between Young/Old cohorts, different individuals)
  subject   = sub("_(Pre|Post)$", "", dal_meta$Col_ID)
)
meta$age   <- factor(meta$age,  levels = c("Young", "Old"))
meta$time  <- factor(meta$time, levels = c("Pre", "Post"))
meta$group <- factor(meta$group,
                     levels = c("Young_Pre", "Young_Post", "Old_Pre", "Old_Post"))

print(table(meta$age, meta$time))
stopifnot(setequal(colnames(mat), meta$sample_id))

# --- CREATE DAList ---

meta_df <- as.data.frame(meta)
rownames(meta_df) <- meta$sample_id

dal <- DAList(
  data       = mat,
  annotation = as.data.frame(ann),
  metadata   = meta_df,
  tags       = list(norm_method = "cycloess")
)

# --- STATISTICAL DESIGN ---

dal <- add_design(dal, "~ 0 + group + (1 | subject)")
colnames(dal$design$design_matrix) <- gsub("^group", "",
                                            colnames(dal$design$design_matrix))

# --- CONTRASTS ---

# NOTE: Reversal = Training_Old - Aging = Old_Post - 2*Old_Pre + Young_Pre.
# Tests whether training in old subjects reverses the aging deficit.
# Non-standard but biologically motivated (Melov et al. 2007, PLOS ONE).
# Acknowledge in methods: not orthogonal to Training_Old and Aging.
#
# NOTE: Interaction yields <=1 FDR-significant protein (n ~ 15/group).
# Training_Old yields 0 FDR hits — consistent with attenuated training
# response in older adults (Robinson et al. 2022; Bechshoft et al. 2024).

dal <- add_contrasts(dal, contrasts_vector = c(
  "Training_Young = Young_Post - Young_Pre",
  "Training_Old = Old_Post - Old_Pre",
  "Aging = Old_Pre - Young_Pre",
  "Interaction = (Old_Post - Old_Pre) - (Young_Post - Young_Pre)",
  "Reversal = (Old_Post - Old_Pre) - (Old_Pre - Young_Pre)"
))

# --- FIT MODEL & EXTRACT RESULTS ---

dal <- fit_limma_model(dal)

within_cor <- dal$eBayes_fit$correlation %||%
  dal$tags$duplicate_correlation %||% NA_real_
if (!is.na(within_cor)) cat(sprintf("Within-subject correlation: %.3f\n", within_cor))

# FDR 0.10: appropriate for exploratory proteomics with small n
# (Benjamini & Hochberg 1995 used 0.10 in examples; Choi et al. 2008,
# BMC Bioinform 9:43, recommend FDR 0-15% for label-free proteomics).
# Pi-score (Pi < 0.05) provides a secondary effect-size-weighted filter.
dal <- extract_DA_results(dal,
                          pval_thresh = cfg$pval_thresh,
                          lfc_thresh  = cfg$lfc_thresh,
                          adj_method  = cfg$adj_method)

# --- SAVE FITTED DAList ---

saveRDS(dal, file.path(cfg$data_dir, "01_limma_DAList.rds"))

# --- GENERATE proteoDA REPORTS ---

tryCatch(
  write_limma_reports(dal, output_dir = cfg$proteoDA_dir, overwrite = TRUE),
  error = function(e) cat(sprintf("write_limma_reports: %s\n", conditionMessage(e)))
)

tryCatch(
  write_limma_plots(dal,
                    grouping_column = "group",
                    output_dir      = cfg$proteoDA_dir,
                    table_columns   = c("uniprot_id", "gene", "protein"),
                    title_column    = "gene",
                    overwrite       = TRUE),
  error = function(e) cat(sprintf("write_limma_plots: %s\n", conditionMessage(e)))
)

# --- BUILD RESULTS FROM dal$results (no disk round-trip) ---

contrast_names <- names(dal$results)
ann_df <- as.data.frame(dal$annotation)

# dal$results stores stats with uniprot_id as rownames; convert to column
# and join annotation (gene, protein, description) for per-contrast CSVs
results_list <- lapply(contrast_names, function(cname) {
  dal$results[[cname]] |>
    rownames_to_column("uniprot_id") |>
    left_join(ann_df, by = "uniprot_id") |>
    mutate(
      pi_score = P.Value ^ abs(logFC),
      sig_pi = case_when(
        pi_score < cfg$pi_thresh & logFC > 0 ~  1L,
        pi_score < cfg$pi_thresh & logFC < 0 ~ -1L,
        TRUE ~ 0L)
    ) |>
    select(-any_of(c("sig.PVal", "sig.FDR"))) |>
    mutate(contrast = cname)
})
names(results_list) <- contrast_names

# --- WRITE PER-CONTRAST CSVs ---

# Match column order of existing per-contrast files: annotation, sample data, stats
# The old pipeline wrote via write_limma_tables which includes sample columns
# We replicate the same structure: ann + sample_data + stats + pi
data_df <- as.data.frame(dal$data)

for (cname in contrast_names) {
  res <- results_list[[cname]]
  # Build per-contrast CSV with same structure as old pipeline:
  # annotation cols, sample intensity cols, stat cols
  out <- ann_df |>
    bind_cols(data_df) |>
    left_join(
      res |> select(uniprot_id, logFC, CI.L, CI.R, average_intensity,
                     t, B, P.Value, adj.P.Val, pi_score, sig_pi),
      by = "uniprot_id"
    )
  write_csv(out, file.path(cfg$per_dir, paste0(cname, ".csv")))
}

# --- BUILD COMBINED RESULTS (wide format) ---

base_df <- bind_cols(
  ann_df |> select(any_of(c("uniprot_id", "protein", "gene", "description"))),
  data_df
)

# Add per-contrast stats as suffixed columns
for (cname in contrast_names) {
  res <- results_list[[cname]]
  stat_cols <- res |>
    select(uniprot_id, logFC, CI.L, CI.R, average_intensity, t, B,
           P.Value, adj.P.Val, pi_score, sig_pi)
  names(stat_cols)[-1] <- paste0(names(stat_cols)[-1], "_", cname)
  base_df <- left_join(base_df, stat_cols, by = "uniprot_id")
}

write_csv(base_df, file.path(cfg$data_dir, "03_combined_results.csv"))

# --- BUILD DA SUMMARY ---

da_summary <- map_dfr(contrast_names, function(cname) {
  res <- results_list[[cname]]
  bind_rows(
    tibble(contrast = cname, type = "up",
           sig.PVal   = sum(res$P.Value < cfg$pval_thresh & res$logFC > 0, na.rm = TRUE),
           sig.FDR    = sum(res$adj.P.Val < cfg$pval_thresh & res$logFC > 0, na.rm = TRUE),
           pval_thresh = cfg$pval_thresh, lfc_thresh = cfg$lfc_thresh,
           p_adj_method = cfg$adj_method,
           sig.Pi     = sum(res$sig_pi == 1, na.rm = TRUE),
           sig.FDR.05 = sum(res$adj.P.Val < 0.05 & res$logFC > 0, na.rm = TRUE),
           sig.FDR.10 = sum(res$adj.P.Val < 0.10 & res$logFC > 0, na.rm = TRUE)),
    tibble(contrast = cname, type = "down",
           sig.PVal   = sum(res$P.Value < cfg$pval_thresh & res$logFC < 0, na.rm = TRUE),
           sig.FDR    = sum(res$adj.P.Val < cfg$pval_thresh & res$logFC < 0, na.rm = TRUE),
           pval_thresh = cfg$pval_thresh, lfc_thresh = cfg$lfc_thresh,
           p_adj_method = cfg$adj_method,
           sig.Pi     = sum(res$sig_pi == -1, na.rm = TRUE),
           sig.FDR.05 = sum(res$adj.P.Val < 0.05 & res$logFC < 0, na.rm = TRUE),
           sig.FDR.10 = sum(res$adj.P.Val < 0.10 & res$logFC < 0, na.rm = TRUE)),
    tibble(contrast = cname, type = "nonsig",
           sig.PVal   = sum(res$P.Value >= cfg$pval_thresh, na.rm = TRUE),
           sig.FDR    = sum(res$adj.P.Val >= cfg$pval_thresh, na.rm = TRUE),
           pval_thresh = cfg$pval_thresh, lfc_thresh = cfg$lfc_thresh,
           p_adj_method = cfg$adj_method,
           sig.Pi     = sum(res$sig_pi == 0, na.rm = TRUE),
           sig.FDR.05 = sum(!res$adj.P.Val < 0.05, na.rm = TRUE),
           sig.FDR.10 = sum(!res$adj.P.Val < 0.10, na.rm = TRUE))
  )
})

write_csv(da_summary, file.path(cfg$data_dir, "02_DA_summary.csv"))

# --- BUILD RESULTS EXCEL ---

wb_results <- createWorkbook()
for (cname in contrast_names) {
  addWorksheet(wb_results, cname)
  writeData(wb_results, cname, results_list[[cname]])
}
addWorksheet(wb_results, "DA_Summary")
writeData(wb_results, "DA_Summary", da_summary)
saveWorkbook(wb_results, file.path(cfg$data_dir, "05_results.xlsx"), overwrite = TRUE)

# --- SUMMARY ---

print(dal$design$contrast_matrix)
print(da_summary)

cat(sprintf("Done: 01_run_dep.R — %d contrasts -> %s/\n",
            length(contrast_names), cfg$data_dir))
