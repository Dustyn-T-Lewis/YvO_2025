#!/usr/bin/env Rscript
# 02_imputation_reports.R — Visualization (no recomputation)
# Input:  c_data/00_report_intermediates.rds (from 01_run_imputation.R)
# Output: b_reports/01_missingness_report.pdf, 02_imputation_report.pdf

# --- Libraries ---------------------------------------------------------------
library(ggplot2)
library(patchwork)
library(gridExtra)
library(dplyr)
library(tidyr)
library(tibble)
library(scales)
library(MsCoreUtils)

# --- Configuration -----------------------------------------------------------
setwd(rprojroot::find_rstudio_root_file())

REPORT_DIR <- "02_Imputation/b_reports"
DATA_DIR   <- "02_Imputation/c_data"
dir.create(REPORT_DIR, showWarnings = FALSE, recursive = TRUE)

THM    <- theme_minimal(base_size = 11)
THM_SM <- theme_minimal(base_size = 8)

# --- Load intermediates ------------------------------------------------------
d <- readRDS(file.path(DATA_DIR, "00_report_intermediates.rds"))

# Unpack for readability
mat          <- d$mat
mat_imp      <- d$mat_imp
was_na       <- d$was_na
meta         <- d$meta
miss_class   <- d$miss_class
miss_by_group <- d$miss_by_group
prot_pct     <- d$prot_pct
pct_miss     <- d$pct_miss
mnar_genes   <- d$mnar_genes
mnar_audit   <- d$mnar_audit
bench_sum    <- d$bench_sum
bin_sum      <- d$bin_sum
ext_sum      <- d$ext_sum
top5         <- d$top5
best         <- d$best
METHOD_TYPE  <- d$METHOD_TYPE
N_ITER       <- d$N_ITER
MASK_FRAC    <- d$MASK_FRAC
mar_miss_vals  <- d$mar_miss_vals
mnar_miss_vals <- d$mnar_miss_vals
total_miss_vals <- d$total_miss_vals
PAL_GT       <- d$PAL_GT
PAL_MAR      <- d$PAL_MAR
PAL_CLASS    <- d$PAL_CLASS
PAL_MTYPE    <- d$PAL_MTYPE
PAL_BIN      <- d$PAL_BIN

# --- Report 1: Missingness & MAR/MNAR classification ------------------------

# A: Per-protein missingness histogram
p_miss_hist <- ggplot(tibble(x = prot_pct), aes(x)) +
  geom_histogram(binwidth = 2, fill = "#4393C3", color = "white", alpha = 0.8) +
  geom_vline(xintercept = c(20, 50), linetype = "dashed", color = "red", alpha = 0.6) +
  annotate("text", x = 22, y = Inf, vjust = 2, hjust = 0, size = 3,
           label = sprintf("Complete: %d\n1-20%%: %d\n20-50%%: %d\n>50%%: %d",
                           sum(prot_pct == 0), sum(prot_pct > 0 & prot_pct <= 20),
                           sum(prot_pct > 20 & prot_pct <= 50), sum(prot_pct > 50))) +
  labs(x = "% missing per protein", y = "Count",
       title = "A: Per-Protein Missingness") + THM

# B: Classification donut (protein counts)
class_counts <- miss_class |> count(classification)
n_total <- sum(class_counts$n)
donut_df <- class_counts |>
  mutate(frac = n / n_total,
         label = sprintf("%s\n%d (%.1f%%)", classification, n, frac * 100),
         ymax = cumsum(frac), ymin = lag(ymax, default = 0),
         ymid = (ymin + ymax) / 2)

p_donut <- ggplot(donut_df, aes(ymax = ymax, ymin = ymin, xmax = 4, xmin = 2.5,
                                 fill = classification)) +
  geom_rect(color = "white", linewidth = 0.6) +
  geom_text(aes(x = 3.25, y = ymid, label = label), size = 3, lineheight = 0.9) +
  annotate("text", x = 0, y = 0,
           label = sprintf("%d proteins\n%.1f%% missing", n_total, pct_miss),
           size = 4, fontface = "bold", lineheight = 1.1) +
  coord_polar(theta = "y") + xlim(c(0, 4.5)) +
  scale_fill_manual(values = PAL_CLASS, guide = "none") +
  theme_void(base_size = 11) +
  labs(title = "B: Protein Classification",
       subtitle = sprintf("MAR values: %d (%.0f%%) | MNAR values: %d (%.0f%%)",
                          mar_miss_vals, mar_miss_vals / total_miss_vals * 100,
                          mnar_miss_vals, mnar_miss_vals / total_miss_vals * 100))

# C: Missingness vs abundance
mc <- miss_class |> filter(classification != "Complete")
p_scatter <- ggplot(mc, aes(mean_intensity, pct_miss, color = classification)) +
  geom_point(alpha = 0.5, size = 1.2) +
  scale_color_manual(values = PAL_MAR) +
  labs(x = "Mean log2 intensity", y = "% missing",
       title = "C: Missingness vs Abundance") +
  THM + theme(legend.position = "bottom")

# D: Per-sample missingness
p_sample <- tibble(Col_ID = colnames(mat),
                   pct = colSums(is.na(mat)) / nrow(mat) * 100) |>
  left_join(meta, by = "Col_ID") |>
  ggplot(aes(reorder(Col_ID, pct), pct, fill = Group_Time)) +
  geom_col(alpha = 0.85) + scale_fill_manual(values = PAL_GT) + coord_flip() +
  labs(x = NULL, y = "% missing proteins", title = "D: Per-Sample Missingness") +
  THM_SM + theme(axis.text.y = element_text(size = 4))

# E: Per-group heatmap (top 50 most-missing proteins)
top_idx <- order(prot_pct, decreasing = TRUE)[1:min(50, sum(prot_pct > 0))]
p_heat <- as_tibble(miss_by_group[top_idx, ], rownames = "gene") |>
  pivot_longer(-gene, names_to = "Group_Time", values_to = "pct") |>
  mutate(gene = factor(gene, levels = rev(rownames(mat)[top_idx]))) |>
  ggplot(aes(Group_Time, gene, fill = pct)) +
  geom_tile(color = "white", linewidth = 0.3) +
  scale_fill_gradient2(low = "white", mid = "#FDDBC7", high = "#B2182B",
                       midpoint = 50, name = "% Miss") +
  labs(x = NULL, y = NULL, title = "E: Per-Group Missingness (top 50)") +
  theme_minimal(base_size = 9) +
  theme(axis.text.y = element_text(size = 5),
        axis.text.x = element_text(angle = 45, hjust = 1))

pdf(file.path(REPORT_DIR, "01_missingness_report.pdf"), width = 14, height = 12)
print(
  (p_miss_hist | p_donut) / (p_scatter | p_sample) +
    plot_annotation(
      title = "Missingness & MAR/MNAR Classification",
      subtitle = sprintf("%d proteins x %d samples | %.2f%% missing overall",
                         nrow(mat), ncol(mat), pct_miss)))
print(
  p_heat +
    plot_annotation(title = "Per-Group Missingness Detail"))
dev.off()

# --- Report 2: Imputation benchmark & QC ------------------------------------

# Page 1: Benchmark ranking (NRMSE + PSS)
visible  <- bench_sum |> filter(median_nrmse < 5)
excluded <- setdiff(bench_sum$method, visible$method)
excl_note <- if (length(excluded) > 0)
  paste0(" | Off-scale: ", paste(excluded, collapse = ", ")) else ""

p_bench_nrmse <- visible |>
  left_join(METHOD_TYPE, by = "method") |>
  ggplot(aes(reorder(method, median_nrmse), median_nrmse, color = type)) +
  geom_point(size = 3.5) +
  geom_errorbar(aes(ymin = median_nrmse - sd_nrmse,
                    ymax = median_nrmse + sd_nrmse),
                width = 0.3, linewidth = 0.5) +
  geom_hline(yintercept = visible$median_nrmse[1], linetype = "dashed",
             color = "#B2182B", alpha = 0.5) +
  scale_color_manual(values = PAL_MTYPE, name = "Type") + coord_flip() +
  labs(x = NULL, y = "NRMSE (lower = better)",
       title = "Imputation Benchmark",
       subtitle = sprintf("%d iter x %d%% masked | Best: %s (%.4f)%s",
                          N_ITER, round(MASK_FRAC * 100), best,
                          bench_sum$mean_nrmse[1], excl_note)) + THM

visible_pss <- bench_sum |> filter(!is.na(mean_pss))
p_bench_pss <- visible_pss |>
  left_join(METHOD_TYPE, by = "method") |>
  ggplot(aes(reorder(method, mean_pss), mean_pss, color = type)) +
  geom_point(size = 3.5) +
  geom_errorbar(aes(ymin = mean_pss - sd_pss, ymax = mean_pss + sd_pss),
                width = 0.3, linewidth = 0.5) +
  scale_color_manual(values = PAL_MTYPE, name = "Type") + coord_flip() +
  labs(x = NULL, y = "Procrustes SS (lower = better)",
       title = "PCA Structure Preservation",
       subtitle = "Lower PSS = less geometric distortion of sample relationships") +
  THM

# Page 2: Observed vs imputed distributions
hist_df <- bind_rows(
  tibble(value = as.vector(mat[!was_na]), source = "Observed"),
  tibble(value = mat_imp[was_na],         source = "Imputed"))

p_hist <- ggplot(hist_df, aes(value, fill = source)) +
  geom_histogram(binwidth = 0.2, alpha = 0.6, position = "identity",
                 color = "white", linewidth = 0.1) +
  scale_fill_manual(values = c(Observed = "#2166AC", Imputed = "#B2182B")) +
  labs(x = "log2 intensity", y = "Count",
       title = "Observed vs Imputed Distributions",
       subtitle = sprintf("%s observed | %s imputed",
                          comma(sum(!was_na)), comma(sum(was_na)))) +
  THM + theme(legend.position = "bottom", legend.title = element_blank())

grp_df <- do.call(rbind, lapply(unique(meta$Group_Time), function(g) {
  cols <- meta$Col_ID[meta$Group_Time == g]
  sm <- mat[, cols, drop = FALSE]
  si <- mat_imp[, cols, drop = FALSE]
  sn <- was_na[, cols, drop = FALSE]
  bind_rows(
    tibble(value = as.vector(sm[!is.na(sm)]), stage = "Observed", Group_Time = g),
    tibble(value = as.vector(si[sn]),         stage = "Imputed",  Group_Time = g))
}))

p_grp <- ggplot(grp_df, aes(value, color = stage, linetype = stage)) +
  geom_density(linewidth = 0.7) + facet_wrap(~Group_Time, scales = "free_y") +
  scale_color_manual(values = c(Observed = "#2166AC", Imputed = "#B2182B")) +
  scale_linetype_manual(values = c(Observed = "solid", Imputed = "dashed")) +
  labs(x = "log2 intensity", title = "Observed vs Imputed by Group") +
  theme_minimal(base_size = 10) +
  theme(legend.position = "bottom", legend.title = element_blank())

# Page 3: MNAR audit + sample-level shift
p_mnar_a <- ggplot(mnar_audit, aes(pre_mean, post_mean, color = pct_miss)) +
  geom_point(alpha = 0.7, size = 1.5) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey40") +
  scale_color_gradient(low = "#FDDBC7", high = "#B2182B", name = "% missing") +
  labs(x = "Pre-imputation mean", y = "Post-imputation mean",
       title = "MNAR Mean Shift") + THM

p_mnar_b <- ggplot(mnar_audit, aes(pct_miss, shift)) +
  geom_point(alpha = 0.6, size = 1.2, color = "#D6604D") +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey40") +
  geom_smooth(method = "loess", se = TRUE, color = "#2166AC", linewidth = 0.8) +
  labs(x = "% missing", y = "Mean shift (post - pre)",
       title = "Shift vs Missingness (MNAR)") + THM

samp_means <- tibble(
  Col_ID = colnames(mat),
  Pre  = colMeans(mat, na.rm = TRUE),
  Post = colMeans(mat_imp)) |>
  left_join(meta |> select(Col_ID, Group_Time), by = "Col_ID")

p_samp <- ggplot(samp_means, aes(Pre, Post, color = Group_Time)) +
  geom_point(size = 3, alpha = 0.7) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
  scale_color_manual(values = PAL_GT) +
  labs(x = "Pre-imputation sample mean", y = "Post-imputation sample mean",
       title = "Sample-Level Mean Shift") +
  THM + theme(legend.position = "bottom")

# Page 4: Extended diagnostics (per-intensity bin + summary table)
p_ext_bin <- bin_sum |>
  filter(method %in% top5) |>
  mutate(bin = factor(bin, levels = c("low", "mid", "high"))) |>
  ggplot(aes(reorder(method, mean_nrmse), mean_nrmse, fill = bin)) +
  geom_col(position = position_dodge(0.7), width = 0.6, alpha = 0.85) +
  geom_errorbar(aes(ymin = mean_nrmse - sd_nrmse, ymax = mean_nrmse + sd_nrmse),
                position = position_dodge(0.7), width = 0.2, linewidth = 0.4) +
  scale_fill_manual(values = PAL_BIN, name = "Intensity\ntertile") +
  coord_flip() +
  labs(x = NULL, y = "NRMSE", title = "A: Per-Intensity-Bin NRMSE (Top 5)") + THM

rank_tbl <- ext_sum |> mutate(across(where(is.numeric), ~ round(.x, 4)))
p_ext_tbl <- patchwork::wrap_elements(
  gridExtra::tableGrob(rank_tbl, rows = NULL,
    theme = gridExtra::ttheme_minimal(base_size = 9,
      core = list(fg_params = list(hjust = 0, x = 0.05)),
      colhead = list(fg_params = list(hjust = 0, x = 0.05, fontface = "bold"))))) +
  labs(title = "B: Extended Summary") +
  theme(plot.title = element_text(size = 11, face = "bold"))

pdf(file.path(REPORT_DIR, "02_imputation_report.pdf"), width = 12, height = 10)
print(
  p_bench_nrmse / p_bench_pss +
    plot_layout(heights = c(2, 1)) +
    plot_annotation(title = "Imputation Benchmark"))
print(
  p_hist / p_grp +
    plot_annotation(
      title = sprintf("Post-Imputation Quality (%s)", best),
      subtitle = sprintf("%d proteins | %d values imputed (%.1f%%)",
                         nrow(mat), sum(was_na), sum(was_na) / length(mat) * 100)))
print(
  (p_mnar_a | p_mnar_b) / p_samp +
    plot_annotation(
      title = sprintf("MNAR Audit (%d proteins) & Sample Shift",
                      length(mnar_genes))))
print(
  (p_ext_bin | p_ext_tbl) +
    plot_annotation(
      title = "Extended Diagnostics (Top 5 Methods)",
      subtitle = "Per-intensity accuracy and method summary"))
# Page 5: Per-Sample LOOCV -----------------------------------------------------
cat("Computing per-sample LOOCV (bpca, 20 iter x 10% mask)...\n")

samp_miss <- tibble(
  Col_ID   = colnames(mat),
  pct_miss = colSums(is.na(mat)) / nrow(mat) * 100) |>
  left_join(meta |> select(Col_ID, Group_Time), by = "Col_ID")

high5 <- samp_miss |> slice_max(pct_miss, n = 5, with_ties = FALSE)
low5  <- samp_miss |> slice_min(pct_miss, n = 5, with_ties = FALSE)
loocv_samples <- bind_rows(
  high5 |> mutate(group = "high_miss"),
  low5  |> mutate(group = "low_miss"))

set.seed(42)
N_LOOCV    <- 20L
LOOCV_FRAC <- 0.10

loocv_res <- do.call(rbind, lapply(seq_len(nrow(loocv_samples)), function(i) {
  sid     <- loocv_samples$Col_ID[i]
  col_idx <- which(colnames(mat) == sid)
  obs_rows <- which(!is.na(mat[, col_idx]))
  obs_sd   <- sd(mat[obs_rows, col_idx])

  nrmse_vals <- vapply(seq_len(N_LOOCV), function(iter) {
    n_mask   <- max(1L, round(length(obs_rows) * LOOCV_FRAC))
    mask_idx <- sample(obs_rows, n_mask)

    mat_masked <- mat
    mat_masked[mask_idx, col_idx] <- NA
    mat_re <- impute_matrix(mat_masked, method = "bpca")

    truth   <- mat[mask_idx, col_idx]
    imputed <- mat_re[mask_idx, col_idx]
    sqrt(mean((truth - imputed)^2)) / obs_sd
  }, numeric(1))

  tibble(Col_ID     = sid,
         Group_Time = loocv_samples$Group_Time[i],
         pct_miss   = loocv_samples$pct_miss[i],
         mean_nrmse = mean(nrmse_vals),
         sd_nrmse   = sd(nrmse_vals),
         group      = loocv_samples$group[i])
}))

write.csv(loocv_res, file.path(DATA_DIR, "11_per_sample_loocv.csv"), row.names = FALSE)
cat(sprintf("  Saved: %s/11_per_sample_loocv.csv\n", DATA_DIR))

ref_nrmse <- bench_sum$mean_nrmse[bench_sum$method == best]
ratio_val <- mean(loocv_res$mean_nrmse[loocv_res$group == "high_miss"]) /
             mean(loocv_res$mean_nrmse[loocv_res$group == "low_miss"])

p_loocv <- loocv_res |>
  mutate(Col_ID = factor(Col_ID, levels = Col_ID[order(pct_miss)])) |>
  ggplot(aes(mean_nrmse, Col_ID, color = group)) +
  geom_point(size = 3) +
  geom_errorbarh(aes(xmin = mean_nrmse - sd_nrmse,
                     xmax = mean_nrmse + sd_nrmse),
                 height = 0.3, linewidth = 0.5) +
  geom_vline(xintercept = ref_nrmse, linetype = "dashed", color = "grey40") +
  scale_color_manual(values = c(high_miss = "#B2182B", low_miss = "#2166AC"),
                     name = NULL) +
  labs(x = "NRMSE (normalized by sample SD)", y = NULL,
       title = "Per-Sample LOOCV (bpca, 20 iter \u00d7 10% mask)",
       subtitle = sprintf("High/Low mean NRMSE ratio: %.2f | Ref line: %s benchmark (%.4f)",
                          ratio_val, best, ref_nrmse)) +
  THM + theme(legend.position = "bottom")

print(
  p_loocv +
    plot_annotation(title = "Per-Sample Imputation Accuracy"))
dev.off()

cat(sprintf("Reports written to %s/\n", REPORT_DIR))
