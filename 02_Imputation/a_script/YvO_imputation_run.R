#!/usr/bin/env Rscript
###############################################################################
#   YvO Imputation — Missingness, MAR/MNAR classification, benchmarking
#   Input:  cycloess-normalized log2 protein matrix (02_normalized.csv)
#
#   Reports:
#     01_missingness_report.pdf  — missingness + MAR/MNAR classification
#     02_imputation_report.pdf   — benchmark, post-imputation QC, MNAR audit
#
#   Methodological references:
#     - Lazar et al. 2016, J Proteome Res 15(4):1116-1125
#       Hybrid MAR/MNAR imputation outperforms single-method approaches
#     - Webb-Robertson et al. 2015, J Proteome Res 14(3):920-930
#       Comprehensive review of imputation strategies for MS proteomics
#     - Hediyeh-zadeh et al. 2023, Mol Cell Proteomics 22(2):100477
#       msImpute EBM for entropy-based MAR/MNAR classification
#     - Webel et al. 2024, J Proteome Res 23(3):1108-1119
#       Imputation evaluation; DE performance degrades >50% missingness
#     - Jin et al. 2021, Sci Rep 11:1556
#       NRMSE benchmarking with masking for imputation evaluation
#
#   Exercise proteomics context:
#     - Hostrup et al. 2022, eLife 11:e69802
#       HIT remodels skeletal muscle proteome; used imputation + DE pipeline
#     - Deshmukh et al. 2021, Nat Commun 12:304
#       Deep muscle proteomics with fiber-type resolution; complete-case analysis
#     - Ubaida-Mohien et al. 2019, eLife 8:e49874
#       Discovery proteomics in aging muscle; >4000 proteins quantified
###############################################################################

if (!requireNamespace("pacman", quietly = TRUE)) install.packages("pacman")
pacman::p_load(
  MsCoreUtils, msImpute, pcaMethods, imputeLCMD, missForest, missMDA,
  vegan, gridExtra, openxlsx,
  ggplot2, patchwork, dplyr, tidyr, tibble, readr, stringr, scales
)

# ─── Paths ───────────────────────────────────────────────────────────────────
setwd(rprojroot::find_rstudio_root_file())
INPUT_CSV  <- "01_normalization/c_data/02_normalized.csv"
REPORT_DIR <- "02_Imputation/b_reports"
DATA_DIR   <- "02_Imputation/c_data"
dir.create(REPORT_DIR, showWarnings = FALSE, recursive = TRUE)
dir.create(DATA_DIR,   showWarnings = FALSE, recursive = TRUE)

# ─── Benchmark methods ───────────────────────────────────────────────────────
# 12 methods: 4 MNAR, 5 MAR, 3 Hybrid (one per MAR backbone).
# Only one MNAR variant per hybrid because the benchmark masks MAR proteins
# only — the MNAR component does not affect NRMSE under this design.
METHODS <- list(
  # MNAR (left-censored / Gaussian sampling from low tail)
  MinProb        = list(method = "MinProb"),
  MinDet         = list(method = "MinDet"),
  QRILC          = list(method = "QRILC"),
  zero           = list(method = "zero"),
  # MAR (global structure-based)
  knn            = list(method = "knn"),
  bpca           = list(method = "bpca"),
  RF             = list(method = "RF"),
  SVD            = list(method = "SVD"),
  imputePCA      = list(method = "imputePCA"),
  # Hybrid (MAR method for MAR proteins, QRILC for MNAR proteins)
  mix_bpca_QRILC = list(method = "mixed", mar = "bpca", mnar = "QRILC"),
  mix_knn_QRILC  = list(method = "mixed", mar = "knn",  mnar = "QRILC"),
  mix_RF_QRILC   = list(method = "mixed", mar = "RF",   mnar = "QRILC")
)

N_ITER <- 20

# Method type lookup (auto-classified from METHODS list)
mtype <- tibble(method = names(METHODS)) %>%
  mutate(type = case_when(
    method %in% c("MinProb", "MinDet", "QRILC", "zero") ~ "MNAR",
    str_starts(method, "mix_") ~ "Hybrid",
    TRUE ~ "MAR"))

# ─── Shared aesthetics ──────────────────────────────────────────────────────
source("04_Figures/shared/palettes.R")
pal_gt    <- GROUP_FILL
pal_mar   <- c(MAR = "#4393C3", MNAR = "#D6604D")
pal_class <- c(Complete = "#4DAF4A", MAR = "#4393C3", MNAR = "#D6604D")
pal_mtyp  <- c(MNAR = "#D6604D", MAR = "#4393C3", Hybrid = "#5AAE61")
thm       <- theme_minimal(base_size = 11)
thm_sm    <- theme_minimal(base_size = 8)

# ─── Helpers ─────────────────────────────────────────────────────────────────

# Dispatch imputation by method specification
run_impute <- function(m, mat, randna) {
  if (m$method == "mixed")
    MsCoreUtils::impute_matrix(mat, method = "mixed", randna = randna,
                               mar = m$mar, mnar = m$mnar)
  else if (m$method == "imputePCA")
    missMDA::imputePCA(mat, ncp = 2, method = "Regularized")$completeObs
  else if (m$method == "SVD")
    pcaMethods::pca(mat, method = "svdImpute", nPcs = 2, verbose = FALSE)@completeObs
  else
    MsCoreUtils::impute_matrix(mat, method = m$method)
}

# Normalized RMSE (SD-normalized, per Webb-Robertson et al. 2015)
nrmse <- function(true_vals, imp_vals) {
  sqrt(mean((true_vals - imp_vals)^2)) / sd(true_vals)
}

# Procrustes Sum of Squares: PCA structure preservation after imputation.
# Uses complete-case proteins from the *original* matrix as reference features,
# then compares PCA scores between original and imputed (which may differ at
# positions that were artificially masked during the benchmark).
pss_metric <- function(original, imputed, n_pc = 5) {
  complete_prots <- complete.cases(original)
  if (sum(complete_prots) < n_pc + 1) return(NA_real_)
  n_pc_use <- min(n_pc, ncol(original) - 1)
  pca_orig <- prcomp(t(original[complete_prots, ]), center = TRUE, scale. = TRUE)
  pca_imp  <- prcomp(t(imputed[complete_prots, ]),  center = TRUE, scale. = TRUE)
  vegan::procrustes(pca_orig$x[, seq_len(n_pc_use)],
                    pca_imp$x[, seq_len(n_pc_use)])$ss
}

###############################################################################
# 1: LOAD DATA
###############################################################################
df  <- read_csv(INPUT_CSV, show_col_types = FALSE)
ann <- df %>% select(uniprot_id, gene, protein, description)
mat <- as.matrix(df[, -(1:4)])
rownames(mat) <- df$gene
if (any(duplicated(df$gene))) {
  dup_genes <- unique(df$gene[duplicated(df$gene)])
  warning(sprintf("%d duplicate gene names found (%s...); using uniprot_id as rownames",
                  length(dup_genes), paste(head(dup_genes, 3), collapse = ", ")))
  rownames(mat) <- df$uniprot_id
}
cat(sprintf("%d proteins x %d samples\n", nrow(mat), ncol(mat)))

# Load canonical metadata from normalisation DAList (not regex-derived)
dal_norm <- readRDS("01_normalization/c_data/03_DAList_normalized.rds")
dal_meta <- as.data.frame(dal_norm$metadata)
meta <- tibble(
  Col_ID     = dal_meta$Col_ID,
  Group      = dal_meta$Group,
  Timepoint  = dal_meta$Timepoint,
  Group_Time = dal_meta$Group_Time
)
stopifnot(setequal(meta$Col_ID, colnames(mat)))
print(count(meta, Group, Timepoint))

###############################################################################
# 2: MISSINGNESS & MAR/MNAR CLASSIFICATION
###############################################################################
prot_miss <- rowSums(is.na(mat))
prot_pct  <- prot_miss / ncol(mat) * 100
obs_means <- rowMeans(mat, na.rm = TRUE)
pct_miss  <- round(sum(is.na(mat)) / length(mat) * 100, 2)

cat(sprintf("Missing: %d / %d (%.2f%%) | Complete: %d\n",
            sum(is.na(mat)), length(mat), pct_miss, sum(prot_miss == 0)))

miss_by_group <- sapply(unique(meta$Group_Time), function(g) {
  cols <- meta$Col_ID[meta$Group_Time == g]
  rowSums(is.na(mat[, cols, drop = FALSE])) / length(cols) * 100
})

has_na <- which(prot_miss > 0 & prot_miss < ncol(mat))
miss_class <- tibble(gene = rownames(mat), n_miss = prot_miss,
                     pct_miss = prot_pct, mean_intensity = obs_means)

# MAR/MNAR classification: msImpute EBM (primary) → k-means (fallback)
mar_result <- tryCatch({
  feat <- msImpute::selectFeatures(mat[has_na, ], method = "ebm",
                                   group = meta$Group_Time)
  list(mar_genes = feat, method = "msImpute_ebm")
}, error = function(e) {
  cat(sprintf("msImpute EBM failed: %s\n", conditionMessage(e)))
  NULL
})

if (is.null(mar_result)) {
  # K-means on standardized (intensity, missingness) — data-driven boundary
  mc_sub <- miss_class %>% filter(n_miss > 0, n_miss < ncol(mat))
  km <- kmeans(scale(cbind(mc_sub$mean_intensity, mc_sub$pct_miss)),
               centers = 2, nstart = 25)
  cl_means <- tapply(mc_sub$mean_intensity, km$cluster, mean)
  mnar_cl <- which.min(cl_means)
  mar_result <- list(
    mar_genes = mc_sub$gene[km$cluster != mnar_cl],
    method = sprintf("k-means (cluster means: %.1f vs %.1f)",
                     cl_means[mnar_cl], cl_means[-mnar_cl]))
}

cat(sprintf("Classification: %s\n", mar_result$method))
miss_class <- miss_class %>%
  mutate(
    classification = case_when(
      n_miss == 0 ~ "Complete",
      gene %in% mar_result$mar_genes ~ "MAR",
      TRUE ~ "MNAR"),
    # Flag proteins with >50% missingness as unreliable for downstream
    # interpretation (Webel et al. 2024: imputation can degrade DE at high MNAR)
    imputation_reliable = classification == "Complete" | pct_miss < 50
  )
print(count(miss_class, classification))

# Group-stratified missingness diagnostic (Fisher exact test per MNAR protein)
# Tests whether missingness is uniformly distributed across the 4 groups.
mnar_genes_for_test <- miss_class$gene[miss_class$classification == "MNAR"]
group_miss_pval <- rep(NA_real_, nrow(miss_class))
names(group_miss_pval) <- miss_class$gene

if (length(mnar_genes_for_test) > 0) {
  for (g in mnar_genes_for_test) {
    ct <- sapply(unique(meta$Group_Time), function(gt) {
      cols <- meta$Col_ID[meta$Group_Time == gt]
      c(missing = sum(is.na(mat[g, cols])), observed = sum(!is.na(mat[g, cols])))
    })
    group_miss_pval[g] <- tryCatch(
      fisher.test(ct, simulate.p.value = TRUE, B = 2000)$p.value,
      error = function(e) NA_real_
    )
  }
  n_sig <- sum(group_miss_pval[mnar_genes_for_test] < 0.05, na.rm = TRUE)
  cat(sprintf("Group-stratified missingness: %d/%d MNAR proteins with group bias (p < 0.05)\n",
              n_sig, length(mnar_genes_for_test)))
}
miss_class$group_miss_pval <- group_miss_pval[miss_class$gene]

###############################################################################
# 3: REPORT 1 — MISSINGNESS
###############################################################################

# A: Per-protein missingness histogram
p_miss_hist <- ggplot(tibble(x = prot_pct), aes(x)) +
  geom_histogram(binwidth = 2, fill = "#4393C3", color = "white", alpha = 0.8) +
  geom_vline(xintercept = c(20, 50), linetype = "dashed", color = "red", alpha = 0.6) +
  annotate("text", x = 22, y = Inf, vjust = 2, hjust = 0, size = 3,
           label = sprintf("Complete: %d\n1-20%%: %d\n20-50%%: %d\n>50%%: %d",
                           sum(prot_pct == 0), sum(prot_pct > 0 & prot_pct <= 20),
                           sum(prot_pct > 20 & prot_pct <= 50), sum(prot_pct > 50))) +
  labs(x = "% missing per protein", y = "Count",
       title = "A: Per-Protein Missingness") + thm

# B: Classification donut (Complete / MAR / MNAR)
class_counts <- miss_class %>% count(classification)
n_total <- sum(class_counts$n)
donut_df <- class_counts %>%
  mutate(frac = n / n_total,
         label = sprintf("%s\n%d (%.1f%%)", classification, n, frac * 100),
         ymax = cumsum(frac), ymin = lag(ymax, default = 0),
         ymid = (ymin + ymax) / 2)

p_class_donut <- ggplot(donut_df, aes(ymax = ymax, ymin = ymin, xmax = 4, xmin = 2.5,
                                       fill = classification)) +
  geom_rect(color = "white", linewidth = 0.6) +
  geom_text(aes(x = 3.25, y = ymid, label = label), size = 3, lineheight = 0.9) +
  annotate("text", x = 0, y = 0,
           label = sprintf("%d\nproteins\n%.1f%% missing", n_total, pct_miss),
           size = 4, fontface = "bold", lineheight = 1.1) +
  coord_polar(theta = "y") + xlim(c(0, 4.5)) +
  scale_fill_manual(values = pal_class, guide = "none") +
  theme_void(base_size = 11) +
  labs(title = "B: Missingness Classification",
       subtitle = "MAR: intensity-independent | MNAR: low-abundance driven")

# C: MAR vs MNAR scatter (intensity vs % missing)
mc <- miss_class %>% filter(classification != "Complete")
p_mar_scatter <- ggplot(mc, aes(mean_intensity, pct_miss, color = classification)) +
  geom_point(alpha = 0.5, size = 1.2) + scale_color_manual(values = pal_mar) +
  labs(x = "Mean log2 intensity", y = "% missing",
       title = "C: Missingness vs Abundance by Type") +
  thm + theme(legend.position = "bottom")

# D: Per-sample missingness
p_samp_miss <- tibble(Col_ID = colnames(mat),
                      pct = colSums(is.na(mat)) / nrow(mat) * 100) %>%
  left_join(meta, by = "Col_ID") %>%
  ggplot(aes(reorder(Col_ID, pct), pct, fill = Group_Time)) +
  geom_col(alpha = 0.85) + scale_fill_manual(values = pal_gt) + coord_flip() +
  labs(x = NULL, y = "% missing proteins", title = "D: Per-Sample Missingness") +
  thm_sm + theme(axis.text.y = element_text(size = 4))

# E: Intensity histogram by type with Complete as reference
p_int_hist <- ggplot(miss_class, aes(mean_intensity, fill = classification)) +
  geom_histogram(binwidth = 0.3, alpha = 0.6, position = "identity", color = "white",
                 linewidth = 0.1) +
  scale_fill_manual(values = pal_class, name = "Type") +
  labs(x = "Mean log2 intensity", y = "Count",
       title = "E: Intensity Distribution by Classification",
       subtitle = "Complete proteins shown as reference") +
  thm + theme(legend.position = "bottom")

# F: Per-group missingness heatmap (top 50)
top_idx <- order(prot_pct, decreasing = TRUE)[1:min(50, sum(prot_pct > 0))]
p_group_heat <- as_tibble(miss_by_group[top_idx, ], rownames = "gene") %>%
  pivot_longer(-gene, names_to = "Group_Time", values_to = "pct") %>%
  mutate(gene = factor(gene, levels = rev(rownames(mat)[top_idx]))) %>%
  ggplot(aes(Group_Time, gene, fill = pct)) +
  geom_tile(color = "white", linewidth = 0.3) +
  scale_fill_gradient2(low = "white", mid = "#FDDBC7", high = "#B2182B",
                       midpoint = 50, name = "% Miss") +
  labs(x = NULL, y = NULL, title = "F: Per-Group Missingness (top 50)") +
  theme_minimal(base_size = 9) +
  theme(axis.text.y = element_text(size = 5),
        axis.text.x = element_text(angle = 45, hjust = 1))

pdf(file.path(REPORT_DIR, "01_missingness_report.pdf"), width = 14, height = 12)
print((p_miss_hist | p_class_donut) / (p_mar_scatter | p_samp_miss) +
        plot_annotation(
          title = "Missingness Diagnostics & MAR/MNAR Classification",
          subtitle = sprintf("%d proteins x %d samples | %.2f%% missing overall",
                             nrow(mat), ncol(mat), pct_miss)))
print((p_int_hist | p_group_heat) +
        plot_annotation(title = "Missingness Detail"))
dev.off()

###############################################################################
# 4: BENCHMARKING
###############################################################################
set.seed(42)
# MsCoreUtils::impute_matrix(method="mixed") defines randna=TRUE as MAR
# (i.e., "random NA"), so TRUE = apply MAR method, FALSE = apply MNAR method.
randna <- setNames(miss_class$classification != "MNAR", miss_class$gene)

cat(sprintf("Benchmarking: %d methods x %d iterations\n", length(METHODS), N_ITER))
res <- vector("list", length(METHODS) * N_ITER); k <- 0L
fail_res <- list(); k_fail <- 0L

for (iter in seq_len(N_ITER)) {
  if (iter %% 5 == 0) cat(sprintf("  Iter %d/%d\n", iter, N_ITER))
  # Mask 10% of observed values from MAR/Complete proteins
  mar_obs_idx <- which(!is.na(mat) & randna[row(mat)])
  mask_idx    <- sample(mar_obs_idx, round(length(mar_obs_idx) * 0.10))
  true_v      <- mat[mask_idx]
  mm <- mat; mm[mask_idx] <- NA

  for (nm in names(METHODS)) {
    result <- tryCatch(
      list(val = run_impute(METHODS[[nm]], mm, randna)),
      error = function(e) list(val = NULL, err = conditionMessage(e)))
    if (is.null(result$val)) {
      cat(sprintf("  %s failed: %s\n", nm, result$err))
      k_fail <- k_fail + 1L
      fail_res[[k_fail]] <- tibble(method = nm, iter = iter, error = result$err)
      next
    }
    k <- k + 1L
    nrmse_val <- nrmse(true_v, result$val[mask_idx])
    pss_val   <- tryCatch(pss_metric(mat, result$val), error = function(e) NA_real_)
    res[[k]]  <- tibble(method = nm, iter = iter, nrmse = nrmse_val, pss = pss_val)
  }
}

bench_df <- bind_rows(res)
bench_sum <- bench_df %>%
  group_by(method) %>%
  summarise(mean_nrmse   = mean(nrmse),
            sd_nrmse     = sd(nrmse),
            median_nrmse = median(nrmse),
            mean_pss     = mean(pss, na.rm = TRUE),
            sd_pss       = sd(pss, na.rm = TRUE),
            .groups = "drop") %>%
  arrange(mean_nrmse)
print(bench_sum)

best <- bench_sum$method[1]
# NOTE: If the winner is a pure MAR method (e.g., bpca), the MAR/MNAR
# classification does not affect imputed values — all proteins get the same
# method. The classification is retained for reporting and the MNAR audit.
# At low overall missingness (~12%), even MNAR-classified proteins retain
# enough observed values for EM-based methods to fit a reliable PC model,
# which may explain why structure-preserving methods outperform left-censored
# approaches here.
cat(sprintf("Best: %s (mean NRMSE = %.4f)\n", best, bench_sum$mean_nrmse[1]))

###############################################################################
# 4b: EXTENDED DIAGNOSTICS — Per-intensity-bin NRMSE (top 5 methods)
###############################################################################
top5 <- bench_sum$method[1:min(5, nrow(bench_sum))]
cat(sprintf("Extended diagnostics for top %d: %s\n", length(top5),
            paste(top5, collapse = ", ")))

# Use all benchmarked proteins (Complete + MAR, i.e. randna == TRUE) for binning.
# The benchmark masks observed values from these proteins, so the per-intensity
# analysis should cover the same set — not just MAR-classified proteins.
bench_genes <- miss_class$gene[randna[miss_class$gene]]

if (length(bench_genes) > 0) {
  bench_means <- rowMeans(mat[bench_genes, , drop = FALSE], na.rm = TRUE)
  bin_breaks  <- unique(quantile(bench_means, probs = c(0, 1/3, 2/3, 1)))
  bin_labels  <- c("low", "mid", "high")[seq_len(length(bin_breaks) - 1)]
  bench_bins  <- cut(bench_means, breaks = bin_breaks, labels = bin_labels,
                     include.lowest = TRUE)
  names(bench_bins) <- bench_genes

  set.seed(42)
  bin_res <- vector("list", length(top5) * N_ITER)
  kb <- 0L
  for (iter in seq_len(N_ITER)) {
    mar_obs_idx <- which(!is.na(mat) & randna[row(mat)])
    mask_idx    <- sample(mar_obs_idx, round(length(mar_obs_idx) * 0.10))
    true_v      <- mat[mask_idx]
    mm <- mat; mm[mask_idx] <- NA
    mask_rows <- rownames(mat)[((mask_idx - 1) %% nrow(mat)) + 1]
    mask_bin  <- bench_bins[mask_rows]

    for (nm in top5) {
      imp_mat_tmp <- tryCatch(run_impute(METHODS[[nm]], mm, randna),
                              error = function(e) NULL)
      if (is.null(imp_mat_tmp)) next
      imp_v <- imp_mat_tmp[mask_idx]
      for (b in bin_labels) {
        sel <- which(mask_bin == b & !is.na(mask_bin))
        if (length(sel) < 5) next
        kb <- kb + 1L
        bin_res[[kb]] <- tibble(method = nm, iter = iter, bin = b,
                                nrmse = nrmse(true_v[sel], imp_v[sel]))
      }
    }
  }
  bin_df  <- bind_rows(bin_res)
  bin_sum <- bin_df %>%
    group_by(method, bin) %>%
    summarise(mean_nrmse = mean(nrmse), sd_nrmse = sd(nrmse), .groups = "drop")
} else {
  cat("  Skipping per-intensity-bin NRMSE: no benchmarkable proteins found.\n")
  bin_df  <- tibble(method = character(), iter = integer(),
                    bin = character(), nrmse = double())
  bin_sum <- tibble(method = character(), bin = character(),
                    mean_nrmse = double(), sd_nrmse = double())
}

# Extended summary table (NRMSE + PSS + per-bin breakdown)
ext_sum <- bench_sum %>%
  filter(method %in% top5) %>%
  select(method, nrmse = mean_nrmse, pss = mean_pss) %>%
  left_join(
    bin_sum %>% select(method, bin, mean_nrmse) %>%
      pivot_wider(names_from = bin, values_from = mean_nrmse,
                  names_prefix = "nrmse_"),
    by = "method"
  )
cat("\n=== Extended Diagnostics Summary (top 5) ===\n")
print(as.data.frame(ext_sum), digits = 4, row.names = FALSE)

###############################################################################
# 5: APPLY BEST METHOD
###############################################################################
set.seed(42)
mat_imp <- run_impute(METHODS[[best]], mat, randna)
stopifnot(sum(is.na(mat_imp)) == 0)

###############################################################################
# 6: REPORT 2 — IMPUTATION
###############################################################################

was_na <- is.na(mat)
top3 <- bench_sum$method[1:3]

# ── Page 1: Benchmark (NRMSE dot plot + PSS) ──────────────────────────────
visible <- bench_sum %>% filter(median_nrmse < 5)
excluded <- setdiff(bench_sum$method, visible$method)
excl_note <- if (length(excluded) > 0) {
  paste0(" | Off scale: ", paste(sprintf("%s (%.1f)", excluded,
    bench_sum$mean_nrmse[match(excluded, bench_sum$method)]), collapse = ", "))
} else ""

p_bench_main <- visible %>%
  left_join(mtype, by = "method") %>%
  ggplot(aes(reorder(method, median_nrmse), median_nrmse, color = type)) +
  geom_point(size = 3.5) +
  geom_errorbar(aes(ymin = median_nrmse - sd_nrmse,
                    ymax = median_nrmse + sd_nrmse),
                width = 0.3, linewidth = 0.5) +
  geom_hline(yintercept = visible$median_nrmse[1], linetype = "dashed",
             color = "#B2182B", alpha = 0.5) +
  scale_color_manual(values = pal_mtyp, name = "Type") + coord_flip() +
  labs(x = NULL, y = "NRMSE (lower = better)", title = "Imputation Benchmark",
       subtitle = sprintf(
         "%d iter x 10%% masked | Best: %s (%.4f)%s\nLower = more accurate point-wise imputation of masked MAR values",
         N_ITER, best, bench_sum$mean_nrmse[1], excl_note)) + thm

visible_pss <- bench_sum %>% filter(!is.na(mean_pss))
p_bench_pss <- visible_pss %>%
  left_join(mtype, by = "method") %>%
  ggplot(aes(reorder(method, mean_pss), mean_pss, color = type)) +
  geom_point(size = 3.5) +
  geom_errorbar(aes(ymin = mean_pss - sd_pss, ymax = mean_pss + sd_pss),
                width = 0.3, linewidth = 0.5) +
  scale_color_manual(values = pal_mtyp, name = "Type") + coord_flip() +
  labs(x = NULL, y = "Procrustes SS (lower = better)",
       title = "PCA Structure Preservation",
       subtitle = "Lower PSS = better preservation of global PCA geometry") +
  thm

# ── Page 2: Post-imputation quality ───────────────────────────────────────
hist_df <- bind_rows(
  tibble(value = as.vector(mat[!was_na]), source = "Observed"),
  tibble(value = mat_imp[was_na],         source = "Imputed"))

p_hist <- ggplot(hist_df, aes(value, fill = source)) +
  geom_histogram(binwidth = 0.2, alpha = 0.6, position = "identity",
                 color = "white", linewidth = 0.1) +
  scale_fill_manual(values = c(Observed = "#2166AC", Imputed = "#B2182B")) +
  labs(x = "log2 intensity", y = "Count",
       title = "Observed vs Imputed Value Distributions",
       subtitle = sprintf("%s observed | %s imputed",
                          scales::comma(sum(!was_na)), scales::comma(sum(was_na)))) +
  thm + theme(legend.position = "bottom", legend.title = element_blank())

grp_df <- do.call(rbind, lapply(unique(meta$Group_Time), function(g) {
  cols <- meta$Col_ID[meta$Group_Time == g]
  sm <- mat[, cols, drop = FALSE]
  si <- mat_imp[, cols, drop = FALSE]
  sn <- was_na[, cols, drop = FALSE]
  bind_rows(
    tibble(value = as.vector(sm[!is.na(sm)]), stage = "Observed", Group_Time = g),
    tibble(value = as.vector(si[sn]),         stage = "Imputed",  Group_Time = g))
}))

p_grp_dens <- ggplot(grp_df, aes(value, color = stage, linetype = stage)) +
  geom_density(linewidth = 0.7) + facet_wrap(~Group_Time, scales = "free_y") +
  scale_color_manual(values = c(Observed = "#2166AC", Imputed = "#B2182B")) +
  scale_linetype_manual(values = c(Observed = "solid", Imputed = "dashed")) +
  labs(x = "log2 intensity", title = "Observed vs Imputed by Group") +
  theme_minimal(base_size = 10) +
  theme(legend.position = "bottom", legend.title = element_blank())

# ── Page 3: MNAR audit + sample-level shift ───────────────────────────────
mnar_genes <- miss_class$gene[miss_class$classification == "MNAR"]
mnar_audit <- tibble(
  gene      = mnar_genes,
  pre_mean  = rowMeans(mat[mnar_genes, ], na.rm = TRUE),
  post_mean = rowMeans(mat_imp[mnar_genes, ]),
  pre_sd    = apply(mat[mnar_genes, ], 1, sd, na.rm = TRUE),
  pct_miss  = prot_pct[mnar_genes],
  shift     = rowMeans(mat_imp[mnar_genes, ]) - rowMeans(mat[mnar_genes, ], na.rm = TRUE),
  effect_d  = shift / pre_sd,
  imputation_reliable = pct_miss < 50
)

p_mnar_a <- ggplot(mnar_audit, aes(pre_mean, post_mean, color = pct_miss)) +
  geom_point(alpha = 0.7, size = 1.5) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey40") +
  scale_color_gradient(low = "#FDDBC7", high = "#B2182B", name = "% missing") +
  labs(x = "Pre-imputation mean log2", y = "Post-imputation mean log2",
       title = "MNAR Mean Shift After Imputation") + thm

p_mnar_b <- ggplot(mnar_audit, aes(pct_miss, shift)) +
  geom_point(alpha = 0.6, size = 1.2, color = "#D6604D") +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey40") +
  geom_smooth(method = "loess", se = TRUE, color = "#2166AC", linewidth = 0.8) +
  labs(x = "% missing", y = "Mean shift (post - pre)",
       title = "Imputation Shift vs Missingness (MNAR)") + thm

samp_means <- tibble(
  Col_ID = colnames(mat),
  Pre    = colMeans(mat, na.rm = TRUE),
  Post   = colMeans(mat_imp)) %>%
  left_join(meta %>% select(Col_ID, Group_Time), by = "Col_ID")

p_samp_mean <- ggplot(samp_means, aes(Pre, Post, color = Group_Time)) +
  geom_point(size = 3, alpha = 0.7) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
  scale_color_manual(values = pal_gt) +
  labs(x = "Pre-imputation sample mean (log2)",
       y = "Post-imputation sample mean (log2)",
       title = "Sample-level mean shift") +
  thm + theme(legend.position = "bottom")

# ── Page 4: Extended diagnostics (per-intensity + summary table) ──────────
p_ext_a <- bin_sum %>%
  mutate(bin = factor(bin, levels = c("low", "mid", "high"))) %>%
  ggplot(aes(reorder(method, mean_nrmse), mean_nrmse, fill = bin)) +
  geom_col(position = position_dodge(0.7), width = 0.6, alpha = 0.85) +
  geom_errorbar(aes(ymin = mean_nrmse - sd_nrmse, ymax = mean_nrmse + sd_nrmse),
                position = position_dodge(0.7), width = 0.2, linewidth = 0.4) +
  scale_fill_manual(values = c(low = "#D6604D", mid = "#F4A582", high = "#92C5DE"),
                    name = "Intensity\ntertile") +
  coord_flip() +
  labs(x = NULL, y = "NRMSE", title = "A: Per-Intensity-Bin NRMSE") + thm

rank_tbl <- ext_sum %>%
  mutate(across(where(is.numeric), ~ round(.x, 4))) %>%
  arrange(nrmse)
p_ext_d <- gridExtra::tableGrob(rank_tbl,
              rows = NULL,
              theme = gridExtra::ttheme_minimal(base_size = 9,
                        core = list(fg_params = list(hjust = 0, x = 0.05)),
                        colhead = list(fg_params = list(hjust = 0, x = 0.05,
                                                         fontface = "bold"))))
p_ext_d <- patchwork::wrap_elements(p_ext_d) +
  labs(title = "B: Extended Summary Ranks") +
  theme(plot.title = element_text(size = 11, face = "bold"))

# ── Write report 2 ────────────────────────────────────────────────────────
pdf(file.path(REPORT_DIR, "02_imputation_report.pdf"), width = 12, height = 10)
print(p_bench_main / p_bench_pss +
      plot_layout(heights = c(2, 1)) +
      plot_annotation(title = "Imputation Benchmark"))
print(p_hist / p_grp_dens +
        plot_annotation(title = sprintf("Post-Imputation Quality (%s)", best),
                        subtitle = sprintf("%d proteins x %d samples | %d values imputed (%.1f%%)",
                                           nrow(mat), ncol(mat), sum(was_na),
                                           sum(was_na) / length(mat) * 100)))
print((p_mnar_a | p_mnar_b) / p_samp_mean +
        plot_annotation(title = sprintf("MNAR Post-Imputation Audit (%d proteins)",
                                        length(mnar_genes))))
print((p_ext_a | p_ext_d) +
        plot_annotation(title = "Extended Diagnostics (Top 5 Methods)",
                        subtitle = "Per-intensity accuracy and method summary"))
dev.off()

cat(sprintf("MNAR shift: %.3f (median), range [%.3f, %.3f]\n",
            median(mnar_audit$shift), min(mnar_audit$shift), max(mnar_audit$shift)))

###############################################################################
# 7: EXPORT — all data outputs consolidated here, numbered by pipeline order
###############################################################################
stopifnot(identical(ann$gene, rownames(mat_imp)))

# 01: Primary deliverables (imputed matrix + DAList)
write_csv(bind_cols(ann, as_tibble(mat_imp)), file.path(DATA_DIR, "01_imputed.csv"))

dal <- readRDS("01_normalization/c_data/03_DAList_normalized.rds")
rownames(mat_imp) <- ann$uniprot_id
dal$data <- mat_imp
imp_meta <- miss_class %>%
  select(gene, n_miss, pct_miss, miss_classification = classification,
         imputation_reliable)
dal$annotation <- merge(dal$annotation, imp_meta, by = "gene",
                        all.x = TRUE, sort = FALSE)
saveRDS(dal, file.path(DATA_DIR, "01_DAList_imputed.rds"))

# 02: MAR/MNAR classification (Section 2 product)
write_csv(miss_class, file.path(DATA_DIR, "02_mar_mnar_classification.csv"))

# 03–06: Benchmark outputs (Section 4 products)
write_csv(bench_sum, file.path(DATA_DIR, "03_benchmark_summary.csv"))
write_csv(bench_df,  file.path(DATA_DIR, "04_benchmark_raw_iterations.csv"))
write_csv(ext_sum,   file.path(DATA_DIR, "05_benchmark_extended.csv"))
write_csv(bin_sum,   file.path(DATA_DIR, "06_benchmark_per_intensity.csv"))

# 07: Imputation mask (Section 6 product)
write.csv(was_na, file.path(DATA_DIR, "07_imputation_mask.csv"), row.names = TRUE)

# 08: MNAR audit (Section 6 product)
write_csv(mnar_audit, file.path(DATA_DIR, "08_mnar_imputation_audit.csv"))

# 09: Pipeline summary
info <- list(
  n_proteins   = nrow(mat), n_samples = ncol(mat), pct_missing = pct_miss,
  n_mar        = sum(miss_class$classification == "MAR"),
  n_mnar       = sum(miss_class$classification == "MNAR"),
  n_complete   = sum(miss_class$classification == "Complete"),
  n_unreliable = sum(!miss_class$imputation_reliable),
  best_method  = best, best_nrmse = bench_sum$mean_nrmse[1])
writeLines(paste(names(info), info, sep = " = "),
           file.path(DATA_DIR, "09_imputation_summary.txt"))

# Conditional: benchmark failures (only if any methods crashed)
if (k_fail > 0) {
  write_csv(bind_rows(fail_res), file.path(DATA_DIR, "benchmark_failures.csv"))
}

# ── Supplementary Excel workbook ──────────────────────────────────────────
add_sheet <- function(wb, sheet_name, title_text, notes_text, data_df) {
  addWorksheet(wb, sheet_name)
  title_style <- createStyle(fontSize = 13, textDecoration = "bold")
  notes_style <- createStyle(fontSize = 10, fontColour = "#555555", wrapText = TRUE)
  header_style <- createStyle(textDecoration = "bold", border = "Bottom",
                              fgFill = "#DCE6F1")
  writeData(wb, sheet_name, x = title_text, startRow = 1, startCol = 1)
  mergeCells(wb, sheet_name, cols = 1:ncol(data_df), rows = 1)
  addStyle(wb, sheet_name, title_style, rows = 1, cols = 1)
  writeData(wb, sheet_name, x = notes_text, startRow = 2, startCol = 1)
  mergeCells(wb, sheet_name, cols = 1:ncol(data_df), rows = 2)
  addStyle(wb, sheet_name, notes_style, rows = 2, cols = 1)
  writeData(wb, sheet_name, x = data_df, startRow = 4, headerStyle = header_style)
  freezePane(wb, sheet_name, firstActiveRow = 5, firstActiveCol = 1)
  setColWidths(wb, sheet_name, cols = 1:ncol(data_df), widths = "auto")
}

wb <- createWorkbook()

readme_df <- tibble(
  Sheet = c("Benchmark_Summary", "Extended_Diagnostics", "MAR_MNAR_Classification",
            "MNAR_Audit", "Benchmark_Iterations", "Pipeline_Summary"),
  Description = c(
    "Imputation benchmark: 12 methods x 20 iterations ranked by mean NRMSE",
    "Extended diagnostics for top 5 methods: per-intensity-tertile NRMSE",
    "Per-protein missingness classification (Complete/MAR/MNAR) with reliability flag",
    "MNAR imputation audit: pre/post means, shift direction, Cohen's d",
    "Raw benchmark iterations (20 x 12 methods): per-iteration NRMSE and PSS",
    "Key pipeline statistics as parameter/value pairs"))
add_sheet(wb, "README",
  "YvO Imputation Supplementary Workbook",
  "Sheet-by-sheet index. See individual sheets for column definitions.",
  readme_df)

add_sheet(wb, "Benchmark_Summary",
  "Imputation Benchmark: 12 Methods x 20 Iterations",
  "mean_nrmse/sd_nrmse/median_nrmse: NRMSE across iterations (lower=better) | mean_pss/sd_pss: Procrustes SS (lower=better PCA preservation) | Ranked by mean NRMSE",
  bench_sum)

add_sheet(wb, "Extended_Diagnostics",
  "Extended Diagnostics: Top 5 Methods",
  "nrmse: overall mean NRMSE | pss: mean Procrustes SS | nrmse_low/nrmse_mid/nrmse_high: NRMSE by intensity tertile (low/mid/high abundance proteins)",
  ext_sum)

add_sheet(wb, "MAR_MNAR_Classification",
  "Per-Protein Missingness Classification",
  "classification: Complete (0% missing) / MAR (intensity-independent) / MNAR (low-abundance driven) | imputation_reliable: TRUE if Complete or <50% missing (Webel et al. 2024) | group_miss_pval: Fisher exact test for group-stratified missingness (MNAR only)",
  miss_class)

add_sheet(wb, "MNAR_Audit",
  "MNAR Imputation Audit",
  "pre_mean/post_mean: mean log2 intensity before/after imputation | pre_sd: SD of observed values | shift: post_mean - pre_mean | effect_d: Cohen's d (shift / pre_sd) | imputation_reliable: TRUE if <50% missing",
  mnar_audit)

add_sheet(wb, "Benchmark_Iterations",
  "Raw Benchmark Iterations (20 x 12)",
  "method: imputation method name | iter: iteration number (1-20) | nrmse: normalized RMSE for this iteration | pss: Procrustes SS for this iteration",
  bench_df)

info_df <- tibble(
  Parameter = names(info),
  Value     = as.character(info))
add_sheet(wb, "Pipeline_Summary",
  "Imputation Pipeline Summary",
  "Key statistics as parameter/value pairs",
  info_df)

saveWorkbook(wb, file.path(DATA_DIR, "imputation_supp.xlsx"), overwrite = TRUE)
cat("Supplementary workbook: imputation_supp.xlsx\n")

cat(sprintf("\n=== YvO Imputation complete === %s | NRMSE %.4f | %d unreliable ===\n",
            best, info$best_nrmse, info$n_unreliable))
