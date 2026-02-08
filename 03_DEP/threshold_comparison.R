################################################################################
#
#   DEqMS Threshold Comparison: BH FDR vs Pi-value (adjusted & raw)
#
#   6 criteria compared:
#     1. BH FDR < 0.05
#     2. BH FDR < 0.10
#     3. Pi_adj < 0.05   where Pi_adj = p_adj^|log2FC|   (BH-adjusted p)
#     4. Pi_adj < 0.10
#     5. Pi_raw < 0.05   where Pi_raw = p_raw^|log2FC|   (unadjusted p)
#     6. Pi_raw < 0.10
#
#   Applied to both DEqMS pipelines (Non-imputed & Imputed)
#
################################################################################

suppressPackageStartupMessages(library(tidyverse))

cat("=== DEqMS Threshold Comparison (6 methods) ===\n\n")

dep_dir <- getwd()  # should be 03_DEP/

# Load both pipelines
pipelines <- list(
  `Non-imputed` = readRDS(file.path(dep_dir, "DEqMS_Non-imputed", "c_data", "DEqMS_results_list.rds")),
  Imputed       = readRDS(file.path(dep_dir, "DEqMS_Imputed",     "c_data", "DEqMS_results_list.rds"))
)

methods <- c("BH_0.05", "BH_0.10", "PiAdj_0.05", "PiAdj_0.10", "PiRaw_0.05", "PiRaw_0.10")

# Build comparison table
rows <- list()

for (pipe_name in names(pipelines)) {
  res_list <- pipelines[[pipe_name]]

  for (cname in names(res_list)) {
    r <- res_list[[cname]]

    # Pi-value with BH-adjusted p:  Pi_adj = sca.adj.pval ^ |logFC|
    r$pi_adj <- r$sca.adj.pval ^ abs(r$logFC)

    # Pi-value with raw (unadjusted) p:  Pi_raw = sca.P.Value ^ |logFC|
    r$pi_raw <- r$sca.P.Value ^ abs(r$logFC)

    for (method in methods) {
      sig <- switch(method,
        BH_0.05     = r$sca.adj.pval < 0.05,
        BH_0.10     = r$sca.adj.pval < 0.10,
        PiAdj_0.05  = r$pi_adj < 0.05,
        PiAdj_0.10  = r$pi_adj < 0.10,
        PiRaw_0.05  = r$pi_raw < 0.05,
        PiRaw_0.10  = r$pi_raw < 0.10
      )
      sig[is.na(sig)] <- FALSE

      rows[[length(rows) + 1]] <- tibble(
        pipeline = pipe_name,
        contrast = cname,
        method   = method,
        n_sig    = sum(sig),
        n_up     = sum(sig & r$logFC > 0),
        n_down   = sum(sig & r$logFC < 0)
      )
    }
  }
}

comp <- bind_rows(rows)

# Wide format
wide <- comp %>%
  pivot_wider(
    id_cols     = c(pipeline, contrast),
    names_from  = method,
    values_from = c(n_sig, n_up, n_down),
    names_glue  = "{method}_{.value}"
  ) %>%
  dplyr::select(
    pipeline, contrast,
    BH_0.05_n_sig, BH_0.05_n_up, BH_0.05_n_down,
    BH_0.10_n_sig, BH_0.10_n_up, BH_0.10_n_down,
    PiAdj_0.05_n_sig, PiAdj_0.05_n_up, PiAdj_0.05_n_down,
    PiAdj_0.10_n_sig, PiAdj_0.10_n_up, PiAdj_0.10_n_down,
    PiRaw_0.05_n_sig, PiRaw_0.05_n_up, PiRaw_0.05_n_down,
    PiRaw_0.10_n_sig, PiRaw_0.10_n_up, PiRaw_0.10_n_down
  )

# Print per pipeline
for (pipe_name in c("Non-imputed", "Imputed")) {
  cat(sprintf("\n--- %s ---\n", pipe_name))
  sub <- wide %>% filter(pipeline == pipe_name) %>% dplyr::select(-pipeline)

  cat(sprintf("\n  %-16s  %s  %s  %s  %s  %s  %s\n",
              "", "BH<0.05", "BH<0.10", "PiAdj<.05", "PiAdj<.10", "PiRaw<.05", "PiRaw<.10"))
  cat(sprintf("  %-16s  %s  %s  %s  %s  %s  %s\n",
              "", "sig up dn", "sig up dn", "sig  up dn", "sig  up dn", "sig  up dn", "sig  up dn"))
  cat("  ", strrep("-", 106), "\n")

  for (i in seq_len(nrow(sub))) {
    cat(sprintf("  %-16s  %3d %3d %3d  %3d %3d %3d  %4d %3d %3d  %4d %3d %3d  %4d %3d %3d  %4d %3d %3d\n",
                sub$contrast[i],
                sub$BH_0.05_n_sig[i],    sub$BH_0.05_n_up[i],    sub$BH_0.05_n_down[i],
                sub$BH_0.10_n_sig[i],    sub$BH_0.10_n_up[i],    sub$BH_0.10_n_down[i],
                sub$PiAdj_0.05_n_sig[i], sub$PiAdj_0.05_n_up[i], sub$PiAdj_0.05_n_down[i],
                sub$PiAdj_0.10_n_sig[i], sub$PiAdj_0.10_n_up[i], sub$PiAdj_0.10_n_down[i],
                sub$PiRaw_0.05_n_sig[i], sub$PiRaw_0.05_n_up[i], sub$PiRaw_0.05_n_down[i],
                sub$PiRaw_0.10_n_sig[i], sub$PiRaw_0.10_n_up[i], sub$PiRaw_0.10_n_down[i]))
  }
}

# Save
write_csv(comp, file.path(dep_dir, "threshold_comparison_long.csv"))
write_csv(wide, file.path(dep_dir, "threshold_comparison_wide.csv"))
cat("\n\nSaved: threshold_comparison_long.csv, threshold_comparison_wide.csv\n")

cat("\n=== Done ===\n")
