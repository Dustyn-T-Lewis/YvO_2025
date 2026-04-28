#!/usr/bin/env Rscript
# Stage 03: Differential Expression
# limma + duplicateCorrelation, 4 contrasts, Pi-score
# Input: cycloess-normalized non-imputed DAList (limma handles NAs per-protein)
#
# Outputs:
#   c_data/01_limma_DAList.rds     fitted proteoDA object
#   c_data/03_combined_results.csv wide-format results for figures
#   c_data/03_DEP_results.xlsx     multi-sheet workbook (core sheets)
#   b_reports/01_proteoDA/         proteoDA HTML reports + static plots

library(dplyr)
library(tidyr)
library(tibble)
library(purrr)
library(proteoDA)
library(openxlsx)

set.seed(42)

DAT <- here::here("03_DEP", "c_data")
RPT <- here::here("03_DEP", "b_reports")
PDA <- file.path(RPT, "01_proteoDA")
dir.create(DAT, recursive = TRUE, showWarnings = FALSE)
dir.create(PDA, recursive = TRUE, showWarnings = FALSE)

PVAL_THRESH <- 0.10
PI_THRESH   <- 0.05

# ── 1. Load from Stage 01 ─────────────────────────────────────────────────────
# Read numeric matrix from CSV for cross-pipeline float reproducibility
# (RDS binary doubles differ at ~1e-15 from CSV round-trip).

df <- readr::read_csv(here::here("01_normalization", "c_data", "02_normalized.csv"),
                      show_col_types = FALSE)
mat <- as.matrix(df[, -(1:4)])
rownames(mat) <- df$gene

dal_norm <- readRDS(here::here("01_normalization", "c_data",
                               "03_DAList_normalized.rds"))
ann  <- as.data.frame(dal_norm$annotation)
meta <- tibble(
  sample_id = dal_norm$metadata$Col_ID,
  age       = factor(dal_norm$metadata$Group,    levels = c("Young", "Old")),
  time      = factor(dal_norm$metadata$Timepoint, levels = c("Pre", "Post")),
  group     = factor(dal_norm$metadata$Group_Time,
                     levels = c("Young_Pre", "Young_Post", "Old_Pre", "Old_Post")),
  subject   = sub("_(Pre|Post)$", "", dal_norm$metadata$Col_ID))  # group-prefixed: Young_S01 ≠ Old_S01

stopifnot(setequal(colnames(mat), meta$sample_id))
message(sprintf("Loaded: %d proteins x %d samples | %.1f%% missing",
                nrow(mat), ncol(mat), 100 * mean(is.na(mat))))

# ── 2. Build DAList ───────────────────────────────────────────────────────────

meta_df <- as.data.frame(meta)
rownames(meta_df) <- meta$sample_id

dal <- DAList(data = mat, annotation = ann, metadata = meta_df,
              tags = list(norm_method = "cycloess"))

# ── 3. Design + contrasts + fit ───────────────────────────────────────────────

dal <- add_design(dal, "~ 0 + group + (1 | subject)")
colnames(dal$design$design_matrix) <- gsub("^group", "",
                                            colnames(dal$design$design_matrix))

dal <- add_contrasts(dal, contrasts_vector = c(
  "Training_Young = Young_Post - Young_Pre",
  "Training_Old   = Old_Post - Old_Pre",
  "Aging          = Old_Pre - Young_Pre",
  "Interaction    = (Old_Post - Old_Pre) - (Young_Post - Young_Pre)"))

# proteoDA runs duplicateCorrelation once. Phipson 2016 recommends iterating
# when array weights are also estimated; for block-only designs (our case)
# duplicateCorrelation is deterministic given the same data + design, so a
# second pass returns the same value. Validated below.
dal <- fit_limma_model(dal)

cor_pass1 <- dal$eBayes_fit$correlation %||%
  dal$tags$duplicate_correlation %||% NA_real_
if (!is.na(cor_pass1)) {
  cor_pass2 <- limma::duplicateCorrelation(
    mat, dal$design$design_matrix, block = meta$subject
  )$consensus.correlation
  message(sprintf("dupCor validation: pass1=%.4f, pass2=%.4f (delta=%.1e)",
                  cor_pass1, cor_pass2, abs(cor_pass1 - cor_pass2)))
}

dal <- extract_DA_results(dal, pval_thresh = PVAL_THRESH,
                          lfc_thresh = 0, adj_method = "BH")

within_cor <- dal$eBayes_fit$correlation %||%
  dal$tags$duplicate_correlation %||% NA_real_
if (!is.na(within_cor)) message(sprintf("Within-subject correlation: %.3f", within_cor))

saveRDS(dal, file.path(DAT, "01_limma_DAList.rds"))

# ── 4. proteoDA reports ──────────────────────────────────────────────────────

tryCatch(
  write_limma_plots(dal, grouping_column = "group", output_dir = PDA,
                    table_columns = c("uniprot_id", "gene", "protein"),
                    title_column = "gene", overwrite = TRUE),
  error = \(e) message("write_limma_plots: ", conditionMessage(e)))

# ── 5. Extract results + Pi-score ─────────────────────────────────────────────

contrast_names <- names(dal$results)
ann_df <- as.data.frame(dal$annotation)

results_list <- lapply(contrast_names, \(cname) {
  dal$results[[cname]] |>
    rownames_to_column("uniprot_id") |>
    left_join(ann_df |> select(uniprot_id, gene, protein, description),
              by = join_by(uniprot_id)) |>
    mutate(pi_score = P.Value ^ abs(logFC),  # Pi-score: Xiao 2014 PMID 22321699
           sig_pi = case_when(
             pi_score < PI_THRESH & logFC > 0 ~  1L,
             pi_score < PI_THRESH & logFC < 0 ~ -1L,
             TRUE ~ 0L),
           contrast = cname) |>
    select(-any_of(c("sig.PVal", "sig.FDR")))
})
names(results_list) <- contrast_names

# ── 6. Combined results (wide format) ────────────────────────────────────────

data_df <- as.data.frame(dal$data)
base_df <- bind_cols(
  ann_df |> select(any_of(c("uniprot_id", "protein", "gene", "description"))),
  data_df)

for (cname in contrast_names) {
  stat_cols <- results_list[[cname]] |>
    select(uniprot_id, logFC, CI.L, CI.R, average_intensity,
           t, B, P.Value, adj.P.Val, pi_score, sig_pi)
  names(stat_cols)[-1] <- paste0(names(stat_cols)[-1], "_", cname)
  base_df <- left_join(base_df, stat_cols, by = join_by(uniprot_id))
}

readr::write_csv(base_df, file.path(DAT, "03_combined_results.csv"))

# ── 7. DA summary ────────────────────────────────────────────────────────────

da_summary <- list_rbind(lapply(contrast_names, \(cname) {
  res <- results_list[[cname]]
  bind_rows(
    tibble(contrast = cname, type = "up",
           sig.PVal = sum(res$P.Value < PVAL_THRESH & res$logFC > 0, na.rm = TRUE),
           sig.FDR  = sum(res$adj.P.Val < PVAL_THRESH & res$logFC > 0, na.rm = TRUE),
           sig.Pi   = sum(res$sig_pi == 1, na.rm = TRUE),
           sig.FDR.05 = sum(res$adj.P.Val < 0.05 & res$logFC > 0, na.rm = TRUE)),
    tibble(contrast = cname, type = "down",
           sig.PVal = sum(res$P.Value < PVAL_THRESH & res$logFC < 0, na.rm = TRUE),
           sig.FDR  = sum(res$adj.P.Val < PVAL_THRESH & res$logFC < 0, na.rm = TRUE),
           sig.Pi   = sum(res$sig_pi == -1, na.rm = TRUE),
           sig.FDR.05 = sum(res$adj.P.Val < 0.05 & res$logFC < 0, na.rm = TRUE)),
    tibble(contrast = cname, type = "nonsig",
           sig.PVal = sum(res$P.Value >= PVAL_THRESH, na.rm = TRUE),
           sig.FDR  = sum(res$adj.P.Val >= PVAL_THRESH, na.rm = TRUE),
           sig.Pi   = sum(res$sig_pi == 0, na.rm = TRUE),
           sig.FDR.05 = sum(!res$adj.P.Val < 0.05, na.rm = TRUE)))
}))

# ── 8. Build xlsx ─────────────────────────────────────────────────────────────

write_sheet <- function(wb, name, data) {
  addWorksheet(wb, name)
  writeData(wb, name, data,
    headerStyle = createStyle(textDecoration = "bold", fgFill = "#DCE6F1"))
  freezePane(wb, name, firstRow = TRUE)
  setColWidths(wb, name, cols = seq_len(ncol(data)), widths = "auto")
}

wb <- createWorkbook()
write_sheet(wb, "combined_results", base_df)
write_sheet(wb, "DA_summary", da_summary)
for (cname in contrast_names) {
  res <- results_list[[cname]] |> arrange(pi_score)
  write_sheet(wb, cname, res)
}
saveWorkbook(wb, file.path(DAT, "03_DEP_results.xlsx"), overwrite = TRUE)

message(sprintf("Done: %d contrasts | %d proteins -> %s/",
                length(contrast_names), nrow(dal$data), DAT))
