#!/usr/bin/env Rscript
################################################################################
#   YvO Sensitivity Analysis — Top-3 norms x (top-3 imps + none) = 12 combos
#
#   Tests whether downstream DEP results are robust to the choice of
#   normalization and imputation method by running all 12 pipeline
#   combinations through the same limma model and comparing results.
#
#   INPUT:  01_normalization/c_data/01_DAList_prenorm.rds
#           02_Imputation/c_data/03_benchmark_summary.csv
#           02_Imputation/c_data/02_mar_mnar_classification.csv
#
#   OUTPUT: 03_DEP/sensitivity_analysis/c_data/01_sensitivity_comparison.csv
#           03_DEP/sensitivity_analysis/c_data/02_pairwise_spearman.csv
#           03_DEP/sensitivity_analysis/c_data/03_dep_counts.csv
#           03_DEP/sensitivity_analysis/c_data/04_cat_curves.csv
#           03_DEP/sensitivity_analysis/c_data/05_sensitivity_results.xlsx
#           03_DEP/sensitivity_analysis/c_data/06_sensitivity_supplementary.xlsx
#           03_DEP/sensitivity_analysis/b_reports/01_sensitivity_comparison.pdf
#
#   Methodological references:
#     - Irizarry et al. 2005, Nat Methods 2(5):345-350  — CAT curves
#     - Xiao et al. 2014, Bioinformatics 30(6):801-807  — Pi-score
#     - Valikangas et al. 2018, Brief Bioinform 19(1):1-11
#       Systematic evaluation of normalization for label-free proteomics
#     - Lazar et al. 2016, J Proteome Res 15(4):1116-1125
#       Hybrid MAR/MNAR imputation benchmarking
#     - Webb-Robertson et al. 2015, J Proteome Res 14(3):920-930
#       Review of imputation strategies for MS proteomics
#
#   NOTE: The Reversal contrast is excluded because it is a linear
#   combination of Aging and Training_Old, making it non-independent.
#   Including it would inflate concordance metrics without adding
#   information about pipeline sensitivity.
################################################################################

if (!requireNamespace("pacman", quietly = TRUE)) install.packages("pacman")
pacman::p_load(proteoDA, MsCoreUtils, pcaMethods, imputeLCMD, missForest, missMDA,
               limma, openxlsx, readr, dplyr, tidyr, tibble, stringr,
               ggplot2, patchwork, scales)

setwd(rprojroot::find_rstudio_root_file())

# Local theme/palette definitions (avoid backwards dependency on 04_Figures)
DIR_COLORS <- c(Up = "#D6604D", Down = "#4393C3", NS = "grey70")
THEME_PUB <- theme_bw(base_size = 10) +
  theme(plot.title       = element_text(face = "bold", size = 12),
        plot.subtitle    = element_text(size = 9, color = "grey30", face = "bold.italic"),
        strip.background = element_blank(),
        strip.text       = element_text(face = "bold", size = 10),
        legend.key.size  = unit(3, "mm"))

# --- Paths ---
report_dir <- "03_DEP/sensitivity_analysis/b_reports"
data_dir   <- "03_DEP/sensitivity_analysis/c_data"
dir.create(report_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(data_dir,   recursive = TRUE, showWarnings = FALSE)

bench_file    <- "02_Imputation/c_data/03_benchmark_summary.csv"
pre_norm_file <- "01_normalization/c_data/01_DAList_prenorm.rds"
mar_mnar_file <- "02_Imputation/c_data/02_mar_mnar_classification.csv"

stopifnot("Benchmark summary not found — run 02_Imputation first"     = file.exists(bench_file))
stopifnot("Pre-norm DAList not found"    = file.exists(pre_norm_file))
stopifnot("MAR/MNAR classification not found — run 02_Imputation first" = file.exists(mar_mnar_file))

# === Helpers ==================================================================
# Imputation dispatcher — mirrors 02_Imputation/YvO_imputation_run.R::run_impute
run_impute_raw <- function(spec, mat, randna) {
  if (spec$method == "mixed")
    MsCoreUtils::impute_matrix(mat, method = "mixed", randna = randna,
                               mar = spec$mar, mnar = spec$mnar)
  else if (spec$method == "imputePCA")
    missMDA::imputePCA(mat, ncp = 2, method = "Regularized")$completeObs
  else if (spec$method == "SVD")
    pcaMethods::pca(mat, method = "svdImpute", nPcs = 2, verbose = FALSE)@completeObs
  else if (spec$method == "GaussianSample") {
    for (i in seq_len(nrow(mat))) {
      na_idx <- is.na(mat[i, ])
      if (any(na_idx) && sum(!na_idx) >= 2)
        mat[i, na_idx] <- rnorm(sum(na_idx), mean(mat[i,], na.rm=TRUE),
                                sd(mat[i,], na.rm=TRUE))
    }
    mat
  } else
    MsCoreUtils::impute_matrix(mat, method = spec$method)
}

# === 1. READ UPSTREAM RANKINGS ================================================

# Normalization methods: canonical (cycloess) + two diverse alternatives
# cycloess = local regression, quantile = distributional, vsn = variance stabilization
top_norms <- c("cycloess", "quantile", "vsn")
cat(sprintf("Norm methods: %s\n", paste(top_norms, collapse = ", ")))

bench <- read_csv(bench_file, show_col_types = FALSE) %>% arrange(mean_nrmse)
top_imps <- bench$method[1:3]
cat(sprintf("Top 3 imps:  %s\n", paste(top_imps, collapse = ", ")))

# Build imputation specs from method names
IMP_SPECS <- list(none = NULL)
for (imp_name in top_imps) {
  if (grepl("^mix_", imp_name)) {
    parts <- strsplit(sub("^mix_", "", imp_name), "_")[[1]]
    IMP_SPECS[[imp_name]] <- list(method = "mixed", mar = parts[1], mnar = parts[2])
  } else {
    IMP_SPECS[[imp_name]] <- list(method = imp_name)
  }
}

# === 2. LOAD DATA =============================================================

dal_raw <- readRDS(pre_norm_file)
dal_raw$metadata$group <- factor(dal_raw$metadata$Group_Time,
                                  levels = c("Young_Pre","Young_Post","Old_Pre","Old_Post"))
meta <- dal_raw$metadata
cat(sprintf("%d proteins x %d samples\n", nrow(dal_raw$data), ncol(dal_raw$data)))

mar_mnar <- read_csv(mar_mnar_file, show_col_types = FALSE)
# Key randna by uniprot_id (matrix rownames) not gene name
# mar_mnar has 'gene' column; map to uniprot_id via DAList annotation
ann_map <- dal_raw$annotation[, c("uniprot_id", "gene")]
mar_mnar_uid <- merge(mar_mnar, ann_map, by = "gene", all.x = TRUE)
randna <- setNames(mar_mnar_uid$classification != "MNAR", mar_mnar_uid$uniprot_id)
randna <- randna[!is.na(names(randna))]  # drop any unmatched

CONTRASTS <- c(
  Training_Young = "Young_Post - Young_Pre",
  Training_Old   = "Old_Post - Old_Pre",
  Aging          = "Old_Pre - Young_Pre",
  Interaction    = "(Old_Post - Old_Pre) - (Young_Post - Young_Pre)"
)

# === 3. RUN 12-COMBO GRID ====================================================
# Each combination: normalize → optionally impute → limma (same model as 03_DEP)

run_pipeline <- function(dal, norm_method, imp_spec, randna_vec) {
  dal_n <- tryCatch(normalize_data(dal, norm_method = norm_method),
                    error = function(e) NULL)
  if (is.null(dal_n)) return(NULL)

  mat <- as.matrix(dal_n$data)

  if (!is.null(imp_spec)) {
    rn <- randna_vec[rownames(mat)]
    rn[is.na(rn)] <- FALSE
    mat <- tryCatch(run_impute_raw(imp_spec, mat, rn), error = function(e) NULL)
    if (is.null(mat)) return(NULL)
  }

  meta   <- dal_n$metadata
  design <- model.matrix(~ 0 + group, data = meta)
  colnames(design) <- gsub("^group", "", colnames(design))

  contrast_mat <- makeContrasts(
    Training_Young = Young_Post - Young_Pre,
    Training_Old   = Old_Post - Old_Pre,
    Aging          = Old_Pre - Young_Pre,
    Interaction    = (Old_Post - Old_Pre) - (Young_Post - Young_Pre),
    levels = design)

  aw <- arrayWeights(mat, design)

  dupcor <- tryCatch(
    duplicateCorrelation(mat, design,
                        block = sub("_(Pre|Post)$", "", meta$Col_ID),
                        weights = aw),
    error = function(e) {
      warning(sprintf("[%s] duplicateCorrelation failed: %s -> correlation = 0",
                      norm_method, e$message), call. = FALSE)
      list(consensus.correlation = 0)
    })

  fit <- lmFit(mat, design, block = sub("_(Pre|Post)$", "", meta$Col_ID),
               correlation = dupcor$consensus.correlation, weights = aw) |>
    contrasts.fit(contrast_mat) |> eBayes(robust = TRUE, trend = TRUE)

  lapply(setNames(seq_along(CONTRASTS), names(CONTRASTS)), function(i) {
    tt <- topTable(fit, coef = i, number = Inf, sort.by = "none")
    tt$gene <- rownames(mat)
    tt
  })
}

# Build combo grid with norm/imp stored as proper columns (not parsed from label)
combos <- expand.grid(norm = top_norms, imp = names(IMP_SPECS), stringsAsFactors = FALSE)
combos$label <- paste0(combos$norm, "_", combos$imp)
imp_clustered_order <- combos %>%
  arrange(factor(imp, levels = names(IMP_SPECS)), norm) %>%
  pull(label)

set.seed(42)
all_results <- list()
for (i in seq_len(nrow(combos))) {
  lab <- combos$label[i]
  cat(sprintf("  [%d/%d] %s\n", i, nrow(combos), lab))
  res <- tryCatch(
    run_pipeline(dal_raw, combos$norm[i], IMP_SPECS[[combos$imp[i]]], randna),
    error = function(e) { cat(sprintf("  FAILED: %s\n", e$message)); NULL })
  if (!is.null(res)) all_results[[lab]] <- res
}

cat(sprintf("Completed: %d / %d combinations\n", length(all_results), nrow(combos)))
method_names <- names(all_results)
n_methods    <- length(method_names)

# === 4. PAIRWISE SPEARMAN RHO ================================================

all_genes <- unique(unlist(lapply(all_results, function(res) res[[1]]$gene)))

spearman_by_contrast <- list()
for (cname in names(CONTRASTS)) {
  t_mat <- sapply(method_names, function(lab) {
    tt <- all_results[[lab]][[cname]]
    t_vals <- setNames(rep(NA_real_, length(all_genes)), all_genes)
    t_vals[tt$gene] <- tt$t
    t_vals
  })
  spearman_mat <- cor(t_mat, method = "spearman", use = "pairwise.complete.obs")
  spearman_by_contrast[[cname]] <- spearman_mat

  mean_rho <- mean(spearman_mat[lower.tri(spearman_mat)], na.rm = TRUE)
  min_rho  <- min(spearman_mat[lower.tri(spearman_mat)], na.rm = TRUE)
  cat(sprintf("  %s: mean rho = %.4f, min rho = %.4f\n", cname, mean_rho, min_rho))
}

spearman_df <- bind_rows(lapply(names(spearman_by_contrast), function(cname) {
  sm <- spearman_by_contrast[[cname]]
  pairs <- which(lower.tri(sm), arr.ind = TRUE)
  tibble(contrast = cname,
         method_a = rownames(sm)[pairs[, 1]],
         method_b = colnames(sm)[pairs[, 2]],
         spearman_rho = sm[pairs])
}))
write_csv(spearman_df, file.path(data_dir, "02_pairwise_spearman.csv"))

# === 5. CAT CURVES ============================================================
# Concordance At the Top (Irizarry et al. 2005): fraction of overlap between
# each method's top-k list and a consensus ranking (mean rank across pipelines)

compute_cat <- function(ranks_a, ranks_b, max_k = 300) {
  shared <- intersect(names(ranks_a), names(ranks_b))
  if (length(shared) < 2) return(rep(NA_real_, max_k))
  ord_a <- names(sort(ranks_a[shared]))
  ord_b <- names(sort(ranks_b[shared]))
  ks <- seq_len(min(max_k, length(shared)))
  sapply(ks, function(k) length(intersect(ord_a[1:k], ord_b[1:k])) / k)
}

# Vectorized consensus ranking: build matrix of ranks, take row means
cat_data_all <- list()
for (cname in names(CONTRASTS)) {
  rank_mat <- sapply(all_results, function(res) {
    tt <- res[[cname]]
    r <- setNames(rep(NA_real_, length(all_genes)), all_genes)
    r[tt$gene] <- rank(-abs(tt$t))
    r
  })
  consensus_ranks <- rowMeans(rank_mat, na.rm = TRUE)
  consensus_ranks <- consensus_ranks[!is.na(consensus_ranks)]

  cat_data <- lapply(method_names, function(lab) {
    tt <- all_results[[lab]][[cname]]
    method_ranks <- setNames(rank(-abs(tt$t)), tt$gene)
    cat_vals <- compute_cat(method_ranks, consensus_ranks, max_k = 300)
    tibble(k = seq_along(cat_vals), concordance = cat_vals, method = lab, contrast = cname)
  }) %>% bind_rows()

  cat_data_all[[cname]] <- cat_data
}
write_csv(bind_rows(cat_data_all), file.path(data_dir, "04_cat_curves.csv"))

# === 6. DEP COUNTS ============================================================
# Pi-score uses Xiao et al. 2014 original formulation: Pi = p^|logFC|
# Threshold Pi < 0.05 (equivalent to -log10(Pi) > 1.3)

dep_summary <- list()
for (lab in method_names) {
  for (cname in names(CONTRASTS)) {
    tt <- all_results[[lab]][[cname]]
    tt$pi_score <- tt$P.Value ^ abs(tt$logFC)

    criteria <- list(
      nominal_p05  = tt %>% filter(!is.na(P.Value)   & P.Value   < 0.05),
      FDR_05       = tt %>% filter(!is.na(adj.P.Val) & adj.P.Val < 0.05),
      FDR_10       = tt %>% filter(!is.na(adj.P.Val) & adj.P.Val < 0.10),
      pi_score_05  = tt %>% filter(!is.na(pi_score)  & pi_score  < 0.05)
    )

    for (crit_name in names(criteria)) {
      sig <- criteria[[crit_name]]
      dep_summary <- c(dep_summary, list(tibble(
        method = lab, contrast = cname, criterion = crit_name,
        n_up = sum(sig$logFC > 0), n_down = sum(sig$logFC < 0),
        n_total = nrow(sig)
      )))
    }
  }
}
dep_df <- bind_rows(dep_summary)
write_csv(dep_df, file.path(data_dir, "03_dep_counts.csv"))

# === 7. SUMMARY TABLE =========================================================

# Use combos data frame for norm/imp (not fragile string parsing)
combo_lookup <- combos %>% select(label, norm, imp)

summary_rows <- list()
for (lab in method_names) {
  mean_rhos <- sapply(names(spearman_by_contrast), function(cname) {
    sm <- spearman_by_contrast[[cname]]
    mean(sm[lab, setdiff(method_names, lab)], na.rm = TRUE)
  })

  dep_fdr05 <- dep_df %>% filter(method == lab, criterion == "FDR_05") %>%
    summarise(
      n_DEP_Training_Y  = n_total[contrast == "Training_Young"],
      n_DEP_Training_O  = n_total[contrast == "Training_Old"],
      n_DEP_Aging       = n_total[contrast == "Aging"],
      n_DEP_Interaction = n_total[contrast == "Interaction"]
    )

  dep_extra <- dep_df %>%
    filter(method == lab, contrast == "Aging") %>%
    summarise(
      n_Aging_nominal_p05 = n_total[criterion == "nominal_p05"],
      n_Aging_FDR_10      = n_total[criterion == "FDR_10"],
      n_Aging_pi_score    = n_total[criterion == "pi_score_05"]
    )

  info <- combo_lookup %>% filter(label == lab)

  summary_rows <- c(summary_rows, list(tibble(
    method             = lab,
    norm               = info$norm,
    imputation         = info$imp,
    mean_spearman      = round(mean(mean_rhos), 4),
    rho_Training_Y     = round(mean_rhos["Training_Young"], 4),
    rho_Training_O     = round(mean_rhos["Training_Old"], 4),
    rho_Aging          = round(mean_rhos["Aging"], 4),
    rho_Interaction    = round(mean_rhos["Interaction"], 4),
    n_DEP_Training_Y   = dep_fdr05$n_DEP_Training_Y,
    n_DEP_Training_O   = dep_fdr05$n_DEP_Training_O,
    n_DEP_Aging        = dep_fdr05$n_DEP_Aging,
    n_DEP_Interaction  = dep_fdr05$n_DEP_Interaction,
    n_Aging_nominal_p05 = dep_extra$n_Aging_nominal_p05,
    n_Aging_FDR_10      = dep_extra$n_Aging_FDR_10,
    n_Aging_pi_score    = dep_extra$n_Aging_pi_score
  )))
}
summary_df <- bind_rows(summary_rows) %>% arrange(desc(mean_spearman))
print(summary_df)
write_csv(summary_df, file.path(data_dir, "01_sensitivity_comparison.csv"))

# === 8. EXCEL WORKBOOK ========================================================

wb <- createWorkbook()

addWorksheet(wb, "Overview")
writeData(wb, "Overview", summary_df)

addWorksheet(wb, "Pairwise Concordance")
writeData(wb, "Pairwise Concordance", spearman_df)

addWorksheet(wb, "DEP Counts")
writeData(wb, "DEP Counts", dep_df)

for (cname in names(spearman_by_contrast)) {
  sheet_name <- paste0("Rho — ", cname)
  addWorksheet(wb, sheet_name)
  sm <- spearman_by_contrast[[cname]]
  sm_df <- as.data.frame(sm) %>% tibble::rownames_to_column("method")
  writeData(wb, sheet_name, sm_df)
}

saveWorkbook(wb, file.path(data_dir, "05_sensitivity_results.xlsx"), overwrite = TRUE)

# === 9. SUPPLEMENTARY EXCEL ===================================================
# Manuscript-ready supplementary with presentation-friendly column names

wb_supp <- createWorkbook()

# Sheet 1: Pipeline grid — which methods were tested
grid_info <- combos %>%
  select(Pipeline = label, Normalization = norm, Imputation = imp) %>%
  mutate(Imputation = ifelse(Imputation == "none", "None (limma handles NAs)", Imputation))
addWorksheet(wb_supp, "S1 — Pipeline Grid")
writeData(wb_supp, "S1 — Pipeline Grid", grid_info)

# Sheet 2: Concordance summary — one row per pipeline, ranked
supp_summary <- summary_df %>%
  select(Pipeline = method, Normalization = norm, Imputation = imputation,
         `Mean Spearman rho` = mean_spearman,
         `rho Training (Young)` = rho_Training_Y,
         `rho Training (Old)` = rho_Training_O,
         `rho Aging` = rho_Aging,
         `rho Interaction` = rho_Interaction)
addWorksheet(wb_supp, "S2 — Concordance Summary")
writeData(wb_supp, "S2 — Concordance Summary", supp_summary)

# Sheet 3: DEP counts (FDR < 0.05) — wide format for readability
dep_fdr05_wide <- dep_df %>%
  filter(criterion == "FDR_05") %>%
  select(Pipeline = method, Contrast = contrast,
         Up = n_up, Down = n_down, Total = n_total) %>%
  arrange(Contrast, Pipeline)
addWorksheet(wb_supp, "S3 — DEP Counts (FDR 0.05)")
writeData(wb_supp, "S3 — DEP Counts (FDR 0.05)", dep_fdr05_wide)

# Sheet 4: DEP counts (all criteria)
dep_all_clean <- dep_df %>%
  mutate(criterion = factor(criterion,
    levels = c("nominal_p05", "FDR_10", "FDR_05", "pi_score_05"),
    labels = c("Nominal p<0.05", "FDR<0.10", "FDR<0.05", "Pi<0.05"))) %>%
  select(Pipeline = method, Contrast = contrast, Criterion = criterion,
         Up = n_up, Down = n_down, Total = n_total) %>%
  arrange(Contrast, Criterion, Pipeline)
addWorksheet(wb_supp, "S4 — DEP Counts (All Criteria)")
writeData(wb_supp, "S4 — DEP Counts (All Criteria)", dep_all_clean)

# Sheet 5: Pairwise Spearman — all pairs
supp_spearman <- spearman_df %>%
  select(Contrast = contrast, `Pipeline A` = method_a,
         `Pipeline B` = method_b, `Spearman rho` = spearman_rho) %>%
  arrange(Contrast, desc(`Spearman rho`))
addWorksheet(wb_supp, "S5 — Pairwise Spearman rho")
writeData(wb_supp, "S5 — Pairwise Spearman rho", supp_spearman)

# Sheet 6: CAT curve summary — concordance at k = 50, 100, 200
cat_summary <- bind_rows(cat_data_all) %>%
  filter(k %in% c(50, 100, 200)) %>%
  select(Pipeline = method, Contrast = contrast, k, Concordance = concordance) %>%
  pivot_wider(names_from = k, values_from = Concordance,
              names_prefix = "CAT_k") %>%
  arrange(Contrast, desc(CAT_k50))
addWorksheet(wb_supp, "S6 — CAT Concordance Summary")
writeData(wb_supp, "S6 — CAT Concordance Summary", cat_summary)

saveWorkbook(wb_supp, file.path(data_dir, "06_sensitivity_supplementary.xlsx"),
             overwrite = TRUE)

# === 10. VISUALS ==============================================================

pdf(file.path(report_dir, "01_sensitivity_comparison.pdf"), width = 14, height = 10)

# Page 1: Spearman heatmaps
heat_df <- bind_rows(lapply(names(spearman_by_contrast), function(cname) {
  sm <- spearman_by_contrast[[cname]]
  expand.grid(method_a = factor(rownames(sm), levels = imp_clustered_order),
              method_b = factor(colnames(sm), levels = imp_clustered_order),
              stringsAsFactors = FALSE) %>%
    mutate(rho = as.vector(sm), contrast = cname)
}))

p_heat <- ggplot(heat_df, aes(x = method_a, y = method_b, fill = rho)) +
  geom_tile(color = "white", linewidth = 0.3) +
  geom_text(aes(label = sprintf("%.3f", rho)), size = 2) +
  facet_wrap(~ contrast, ncol = 2) +
  scale_fill_gradient2(low = "#B2182B", mid = "#F7F7F7", high = "#2166AC",
                       midpoint = 0.95,
                       limits = range(heat_df$rho, na.rm = TRUE),
                       name = "Spearman\nrho") +
  labs(title = "Pairwise Spearman rho of t-statistics",
       subtitle = "Higher = methods rank proteins more similarly",
       x = NULL, y = NULL) +
  THEME_PUB +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        panel.grid = element_blank())
print(p_heat)

# Page 2: CAT curves
cat_df <- bind_rows(cat_data_all)
cat_plots <- lapply(names(CONTRASTS), function(cname) {
  ggplot(cat_df %>% filter(contrast == cname),
         aes(x = k, y = concordance,
             color = factor(method, levels = imp_clustered_order))) +
    geom_line(linewidth = 0.8) +
    geom_hline(yintercept = 0.7, linetype = "dashed", alpha = 0.4) +
    labs(x = "Top k proteins", y = "Concordance", title = cname, color = "Pipeline") +
    THEME_PUB +
    scale_color_brewer(palette = "Paired") +
    ylim(0, 1)
})
p_cat <- wrap_plots(cat_plots, ncol = 1) +
  plot_annotation(
    title = "Concordance At the Top (CAT) curves",
    subtitle = "Each method vs consensus ranking (Irizarry et al. 2005) | Dashed = 0.70",
    theme = theme(plot.title = element_text(face = "bold")))
print(p_cat)

# Page 3: DEP counts by direction (FDR < 0.05)
dep_fdr05_long <- dep_df %>%
  filter(criterion == "FDR_05") %>%
  mutate(method = factor(method, levels = imp_clustered_order)) %>%
  pivot_longer(c(n_up, n_down), names_to = "direction", values_to = "count") %>%
  mutate(direction = ifelse(direction == "n_up", "Up", "Down"),
         count = ifelse(direction == "Down", -count, count))

p_dep <- ggplot(dep_fdr05_long, aes(x = method, y = count, fill = direction)) +
  geom_col(width = 0.7) +
  geom_hline(yintercept = 0, color = "grey30") +
  facet_wrap(~ contrast, scales = "free_y", ncol = 1) +
  scale_fill_manual(values = DIR_COLORS[c("Up", "Down")]) +
  labs(x = NULL, y = "DEP count (up / down)",
       title = "DEP counts per pipeline (FDR < 0.05)", fill = "Direction") +
  THEME_PUB +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
print(p_dep)

# Page 4: DEP counts by significance criterion
dep_criteria <- dep_df %>%
  mutate(method = factor(method, levels = imp_clustered_order),
         criterion = factor(criterion,
    levels = c("nominal_p05", "FDR_10", "FDR_05", "pi_score_05"),
    labels = c("Nominal p<0.05", "FDR<0.10", "FDR<0.05", "Pi<0.05")))

p_criteria <- ggplot(dep_criteria, aes(x = method, y = n_total, fill = criterion)) +
  geom_col(position = "dodge", width = 0.7) +
  facet_wrap(~ contrast, scales = "free_y", ncol = 1) +
  scale_fill_manual(values = c("Nominal p<0.05" = "#92C5DE",
                                "FDR<0.10"       = "#F4A582",
                                "FDR<0.05"       = "#D6604D",
                                "Pi<0.05"        = "#5AAE61")) +
  labs(x = NULL, y = "Total DEPs",
       title = "DEP counts across significance criteria",
       subtitle = "Pi = p^|logFC| (Xiao et al. 2014); threshold Pi < 0.05",
       fill = "Criterion") +
  THEME_PUB +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
print(p_criteria)

dev.off()

# === INTERPRETATION ===========================================================

overall_rho <- mean(spearman_df$spearman_rho, na.rm = TRUE)
cat(sprintf("\nOverall mean Spearman rho: %.4f\n", overall_rho))
if (overall_rho > 0.95) {
  cat("Methods agree closely (rho > 0.95). Choice matters little.\n")
} else if (overall_rho > 0.85) {
  cat("Methods mostly agree (rho > 0.85). Some variation in extremes.\n")
} else {
  cat("Methods diverge meaningfully (rho < 0.85). Inspect contrasts individually.\n")
}

cat(sprintf("Top pipeline: %s (rho = %.4f)\n",
            summary_df$method[1], summary_df$mean_spearman[1]))
cat(sprintf("Done: %d pipelines x %d contrasts\n", n_methods, length(CONTRASTS)))
