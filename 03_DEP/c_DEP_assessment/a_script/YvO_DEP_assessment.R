################################################################################
#
#   YvO DEP Assessment — Imputed vs Non-imputed Comparison
#
#   Focus: limma results across Imputed vs Non-imputed data (2 cells per
#          contrast = 8 total result sets). DEqMS concordance is a brief
#          supplementary section.
#
#   Sections:
#     0.  Setup & load data
#     1.  Threshold sensitivity
#     2.  Data concordance (Imputed vs Non-imputed)
#     3.  Imputation-dependency flags
#     4.  Robustness classification
#     5.  Normalization sensitivity analysis
#     6.  Sensitivity power analysis
#     7.  Confident effect sizes (topconfects)
#     8.  DEqMS concordance (supplementary)
#     9.  Summary workbook
#    10.  Session info
#
#   Design: 2x2 factorial (Age x Time) with repeated measures on subject
#   Contrasts: Training_Young, Training_Old, Aging, Interaction
#
################################################################################

# ==============================================================================
# 0.  SETUP & LOAD DATA
# ==============================================================================

cat("=== YvO DEP Assessment ===\n\n")
cat(">> 0 -- Setup & Load Data\n")

suppressPackageStartupMessages({
  library(tidyverse)
  library(patchwork)
  library(ggrepel)
  library(eulerr)
  library(pwr)
  library(topconfects)
  library(openxlsx)
  library(limma)
  library(proteoDA)
  library(fgsea)
  library(msigdbr)
})

# Paths
base_dir   <- normalizePath(file.path(dirname(getwd()), "..", ".."), mustWork = TRUE)
REPORT_DIR <- file.path(base_dir, "03_DEP", "c_DEP_assessment", "b_reports")
DATA_DIR   <- file.path(base_dir, "03_DEP", "c_DEP_assessment", "c_data")

dir.create(REPORT_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(DATA_DIR,   recursive = TRUE, showWarnings = FALSE)

cat("   Base directory:", base_dir, "\n")

# Contrast names
contrast_names <- c("Training_Young", "Training_Old", "Aging", "Interaction")

# Load 8 limma result sets: 4 contrasts x 2 data types
limma_results <- list()
for (dtype in c("Imputed", "Non-imputed")) {
  per_dir <- file.path(base_dir, "03_DEP", dtype, "a_limma", "c_data", "per_contrast_results")
  for (cname in contrast_names) {
    fpath <- file.path(per_dir, paste0(cname, ".csv"))
    key   <- paste(dtype, cname, sep = ".")
    limma_results[[key]] <- read_csv(fpath, show_col_types = FALSE)
    cat(sprintf("   Loaded: %s (%d proteins)\n", key, nrow(limma_results[[key]])))
  }
}

cat(sprintf("   Total result sets loaded: %d\n", length(limma_results)))

# ==============================================================================
# 1.  THRESHOLD SENSITIVITY
# ==============================================================================

cat("\n>> 1 -- Threshold Sensitivity\n")

fdr_thresholds <- c(0.01, 0.05, 0.10)
lfc_thresholds <- c(0, 0.26, 0.58, 1.0)

thresh_rows <- list()
for (dtype in c("Imputed", "Non-imputed")) {
  for (cname in contrast_names) {
    key <- paste(dtype, cname, sep = ".")
    res <- limma_results[[key]]
    for (fdr_t in fdr_thresholds) {
      for (lfc_t in lfc_thresholds) {
        n_sig <- sum(res$adj.P.Val < fdr_t & abs(res$logFC) > lfc_t, na.rm = TRUE)
        thresh_rows <- c(thresh_rows, list(tibble(
          data_type = dtype, contrast = cname,
          FDR_threshold = fdr_t, logFC_threshold = lfc_t,
          n_DEP = n_sig
        )))
      }
    }
  }
}

thresh_df <- bind_rows(thresh_rows)
write_csv(thresh_df, file.path(DATA_DIR, "threshold_comparison.csv"))

# Faceted grouped bar chart
thresh_df <- thresh_df %>%
  mutate(
    FDR_label = paste0("FDR < ", FDR_threshold),
    LFC_label = paste0("|log2FC| > ", logFC_threshold)
  )

p_thresh <- ggplot(thresh_df, aes(x = LFC_label, y = n_DEP, fill = data_type)) +
  geom_col(position = position_dodge(width = 0.8), width = 0.7) +
  facet_grid(FDR_label ~ contrast, scales = "free_y") +
  scale_fill_manual(values = c("Imputed" = "#4A90D9", "Non-imputed" = "#D94A4A")) +
  labs(x = "log2FC threshold", y = "DEP count", fill = "Data type",
       title = "Threshold sensitivity: Imputed vs Non-imputed") +
  theme_bw(base_size = 10) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "bottom")

ggsave(file.path(REPORT_DIR, "01_threshold_sensitivity.pdf"), p_thresh,
       width = 12, height = 8)
cat("   Saved: threshold_comparison.csv, 01_threshold_sensitivity.pdf\n")

# ==============================================================================
# 2.  DATA CONCORDANCE (Imputed vs Non-imputed)
# ==============================================================================

cat("\n>> 2 -- Data Concordance (Imputed vs Non-imputed)\n")

concord_rows <- list()
euler_plots  <- list()
scatter_plots <- list()

for (cname in contrast_names) {
  imp_key  <- paste("Imputed", cname, sep = ".")
  nimp_key <- paste("Non-imputed", cname, sep = ".")

  imp  <- limma_results[[imp_key]]
  nimp <- limma_results[[nimp_key]]

  # Inner join on uniprot_id
  merged <- inner_join(imp, nimp, by = "uniprot_id", suffix = c(".imp", ".nimp"))

  # Significant sets (FDR < 0.05)
  sig_imp  <- imp$uniprot_id[imp$adj.P.Val < 0.05]
  sig_nimp <- nimp$uniprot_id[nimp$adj.P.Val < 0.05]

  # Jaccard
  union_n   <- length(union(sig_imp, sig_nimp))
  intersect_n <- length(intersect(sig_imp, sig_nimp))
  jaccard <- if (union_n > 0) intersect_n / union_n else NA_real_

  # Spearman correlation of logFC
  spearman_r <- cor(merged$logFC.imp, merged$logFC.nimp,
                    method = "spearman", use = "complete.obs")

  concord_rows <- c(concord_rows, list(tibble(
    contrast = cname,
    n_shared_proteins = nrow(merged),
    n_sig_imputed = length(sig_imp),
    n_sig_nonimputed = length(sig_nimp),
    n_sig_both = intersect_n,
    jaccard_FDR05 = round(jaccard, 3),
    spearman_logFC = round(spearman_r, 3)
  )))

  # Euler diagram
  euler_input <- c(
    "Imputed"         = length(sig_imp) - intersect_n,
    "Non-imputed"     = length(sig_nimp) - intersect_n,
    "Imputed&Non-imputed" = intersect_n
  )
  euler_fit <- euler(euler_input)
  euler_plots[[cname]] <- plot(euler_fit,
    fills = list(fill = c("#4A90D9", "#D94A4A"), alpha = 0.5),
    labels = TRUE,
    quantities = TRUE,
    main = cname
  )

  # logFC scatter
  scatter_plots[[cname]] <- ggplot(merged, aes(x = logFC.imp, y = logFC.nimp)) +
    geom_point(alpha = 0.3, size = 0.8) +
    geom_abline(slope = 1, intercept = 0, color = "red", linetype = "dashed") +
    labs(x = "logFC (Imputed)", y = "logFC (Non-imputed)",
         title = sprintf("%s (rho = %.3f)", cname, spearman_r)) +
    theme_bw(base_size = 10)
}

concord_df <- bind_rows(concord_rows)
write_csv(concord_df, file.path(DATA_DIR, "data_concordance.csv"))

# Save euler diagrams
pdf(file.path(REPORT_DIR, "02a_data_euler.pdf"), width = 10, height = 8)
par(mfrow = c(2, 2))
for (cname in contrast_names) {
  print(euler_plots[[cname]])
}
dev.off()

# Save logFC scatter plots
p_scatter <- wrap_plots(scatter_plots, ncol = 2) +
  plot_annotation(title = "logFC concordance: Imputed vs Non-imputed")
ggsave(file.path(REPORT_DIR, "02b_data_logFC_scatter.pdf"), p_scatter,
       width = 10, height = 8)

cat("   Saved: data_concordance.csv, 02a_data_euler.pdf, 02b_data_logFC_scatter.pdf\n")

# ==============================================================================
# 3.  IMPUTATION-DEPENDENCY FLAGS
# ==============================================================================

cat("\n>> 3 -- Imputation-Dependency Flags\n")

flag_rows <- list()
for (cname in contrast_names) {
  imp_key  <- paste("Imputed", cname, sep = ".")
  nimp_key <- paste("Non-imputed", cname, sep = ".")

  imp  <- limma_results[[imp_key]]
  nimp <- limma_results[[nimp_key]]

  sig_imp  <- imp$uniprot_id[imp$adj.P.Val < 0.05]
  sig_nimp <- nimp$uniprot_id[nimp$adj.P.Val < 0.05]

  all_ids <- union(imp$uniprot_id, nimp$uniprot_id)

  for (uid in all_ids) {
    in_imp  <- uid %in% sig_imp
    in_nimp <- uid %in% sig_nimp
    flag <- case_when(
      in_imp & in_nimp   ~ "robust",
      in_imp & !in_nimp  ~ "imputation_dependent",
      !in_imp & in_nimp  ~ "imputation_masked",
      TRUE               ~ "never_sig"
    )
    gene_name <- NA_character_
    if (uid %in% imp$uniprot_id) {
      gene_name <- imp$gene[imp$uniprot_id == uid][1]
    } else if (uid %in% nimp$uniprot_id) {
      gene_name <- nimp$gene[nimp$uniprot_id == uid][1]
    }
    flag_rows <- c(flag_rows, list(tibble(
      uniprot_id = uid, gene = gene_name, contrast = cname, flag = flag
    )))
  }
}

flag_df <- bind_rows(flag_rows)
write_csv(flag_df, file.path(DATA_DIR, "imputation_dependency_flags.csv"))
cat(sprintf("   Flagged %d protein-contrast pairs\n", nrow(flag_df)))

flag_summary <- flag_df %>% count(contrast, flag) %>% pivot_wider(names_from = flag, values_from = n, values_fill = 0)
print(flag_summary)
cat("   Saved: imputation_dependency_flags.csv\n")

# ==============================================================================
# 4.  ROBUSTNESS CLASSIFICATION
# ==============================================================================

cat("\n>> 4 -- Robustness Classification\n")

robust_df <- flag_df %>%
  mutate(classification = case_when(
    flag == "robust"               ~ "Robust",
    flag %in% c("imputation_dependent", "imputation_masked") ~ "Conditional",
    TRUE                           ~ "NS"
  ))

write_csv(robust_df, file.path(DATA_DIR, "robustness_classification.csv"))

p_robust <- ggplot(robust_df, aes(x = contrast, fill = classification)) +
  geom_bar(position = "stack") +
  scale_fill_manual(values = c("Robust" = "#2E8B57", "Conditional" = "#DAA520", "NS" = "grey70")) +
  labs(x = "Contrast", y = "Protein count", fill = "Classification",
       title = "Robustness: Imputed vs Non-imputed agreement") +
  theme_bw(base_size = 11) +
  theme(axis.text.x = element_text(angle = 30, hjust = 1))

ggsave(file.path(REPORT_DIR, "03_robustness_barplot.pdf"), p_robust,
       width = 8, height = 5)

cat("   Saved: robustness_classification.csv, 03_robustness_barplot.pdf\n")

# ==============================================================================
# 5.  NORMALIZATION SENSITIVITY ANALYSIS
# ==============================================================================

cat("\n>> 5 -- Normalization Sensitivity Analysis\n")

# Load pre-normalization DAList
prenorm_path <- file.path(base_dir, "01_normalization", "c_data", "00_DAList_prenorm.rds")
if (!file.exists(prenorm_path)) {
  cat("   WARNING: 00_DAList_prenorm.rds not found. Skipping normalization sensitivity.\n")
  cat("   Run 01_normalization/a_script/YvO_normalization_run.R first to generate this file.\n")
} else {
  dal_prenorm <- readRDS(prenorm_path)
  cat(sprintf("   Pre-norm DAList: %d proteins x %d samples\n",
              nrow(dal_prenorm$data), ncol(dal_prenorm$data)))

  # Normalization methods to test
  norm_methods <- c("cycloess", "quantile", "median", "loess", "vsn", "rlr", "global", "none")

  # Hallmark gene sets
  hallmark_df <- msigdbr(species = "Homo sapiens", category = "H")
  hallmark_sets <- split(hallmark_df$gene_symbol, hallmark_df$gs_name)

  # Store rankings per method per contrast
  all_rankings <- list()
  all_pathways <- list()

  for (method in norm_methods) {
    cat(sprintf("   Testing normalization: %s\n", method))

    tryCatch({
      # Normalize
      dal_norm <- normalize_data(dal_prenorm, norm_method = method)

      # Build limma model (same design as DEP scripts)
      dal_norm <- add_design(dal_norm, "~ 0 + group + (1 | subject)")
      colnames(dal_norm$design$design_matrix) <- gsub("^group", "",
        colnames(dal_norm$design$design_matrix))

      dal_norm <- add_contrasts(dal_norm, contrasts_vector = c(
        "Training_Young = Young_Post - Young_Pre",
        "Training_Old = Old_Post - Old_Pre",
        "Aging = Old_Pre - Young_Pre",
        "Interaction = (Old_Post - Old_Pre) - (Young_Post - Young_Pre)"
      ))

      dal_norm <- fit_limma_model(dal_norm)
      dal_norm <- extract_DA_results(dal_norm, pval_thresh = 0.05,
                                      lfc_thresh = 0, adj_method = "BH")

      # Extract per-contrast rankings
      for (cname in contrast_names) {
        res <- dal_norm$results[[cname]]
        if (is.null(res)) next

        # Rank by |moderated t|
        res <- res %>%
          mutate(abs_t = abs(t)) %>%
          arrange(desc(abs_t)) %>%
          mutate(rank = row_number())

        key <- paste(method, cname, sep = ".")
        all_rankings[[key]] <- res %>%
          dplyr::select(uniprot_id, gene, logFC, t, abs_t, rank)

        # fgsea with moderated t-statistic
        gene_ranks <- setNames(res$t, res$gene)
        gene_ranks <- gene_ranks[!is.na(names(gene_ranks)) & !duplicated(names(gene_ranks))]

        fgsea_res <- fgsea(pathways = hallmark_sets, stats = gene_ranks,
                           minSize = 15, maxSize = 500)
        top10 <- fgsea_res %>% arrange(pval) %>% head(10)
        all_pathways[[key]] <- top10 %>% dplyr::select(pathway, pval, NES)
      }
    }, error = function(e) {
      cat(sprintf("   WARNING: %s normalization failed: %s\n", method, e$message))
    })
  }

  # Protein ranking concordance: top-50 overlap matrix
  ranking_rows <- list()
  for (cname in contrast_names) {
    ref_key <- paste("cycloess", cname, sep = ".")
    if (is.null(all_rankings[[ref_key]])) next
    ref_top50 <- all_rankings[[ref_key]]$uniprot_id[1:50]

    for (method in norm_methods) {
      key <- paste(method, cname, sep = ".")
      if (is.null(all_rankings[[key]])) next
      test_top50 <- all_rankings[[key]]$uniprot_id[1:50]
      overlap <- length(intersect(ref_top50, test_top50))
      jaccard <- overlap / length(union(ref_top50, test_top50))
      ranking_rows <- c(ranking_rows, list(tibble(
        contrast = cname, method = method,
        top50_overlap_with_cycloess = overlap,
        jaccard = round(jaccard, 3)
      )))
    }
  }

  rank_concord <- bind_rows(ranking_rows)
  write_csv(rank_concord, file.path(DATA_DIR, "norm_sensitivity_rankings.csv"))

  # Normalization-sensitive flags: rank shift > 50 positions between any pair
  flag_norm_rows <- list()
  for (cname in contrast_names) {
    method_ranks <- list()
    for (method in norm_methods) {
      key <- paste(method, cname, sep = ".")
      if (!is.null(all_rankings[[key]])) {
        r <- all_rankings[[key]] %>% dplyr::select(uniprot_id, rank)
        names(r)[2] <- method
        method_ranks[[method]] <- r
      }
    }
    if (length(method_ranks) < 2) next

    # Merge all method ranks
    rank_mat <- method_ranks[[1]]
    for (i in 2:length(method_ranks)) {
      rank_mat <- full_join(rank_mat, method_ranks[[i]], by = "uniprot_id")
    }

    rank_cols <- setdiff(names(rank_mat), "uniprot_id")
    rank_values <- rank_mat[, rank_cols]

    max_rank_shift <- apply(rank_values, 1, function(x) {
      x <- x[!is.na(x)]
      if (length(x) < 2) return(NA_real_)
      max(x) - min(x)
    })
    min_rank <- apply(rank_values, 1, min, na.rm = TRUE)
    max_rank <- apply(rank_values, 1, max, na.rm = TRUE)

    # Get gene names
    ref_key <- paste("cycloess", cname, sep = ".")
    gene_lookup <- all_rankings[[ref_key]] %>% dplyr::select(uniprot_id, gene)

    flag_norm_rows <- c(flag_norm_rows, list(
      tibble(
        uniprot_id = rank_mat$uniprot_id,
        contrast = cname,
        max_rank_shift = max_rank_shift,
        min_rank = min_rank,
        max_rank = max_rank,
        sensitive = max_rank_shift > 50
      ) %>% left_join(gene_lookup, by = "uniprot_id")
    ))
  }

  norm_flags <- bind_rows(flag_norm_rows)
  write_csv(norm_flags, file.path(DATA_DIR, "norm_sensitivity_flags.csv"))

  n_sensitive <- sum(norm_flags$sensitive, na.rm = TRUE)
  cat(sprintf("   Normalization-sensitive protein-contrast pairs: %d\n", n_sensitive))

  # --- Figures for 04_norm_sensitivity.pdf ---
  pdf(file.path(REPORT_DIR, "04_norm_sensitivity.pdf"), width = 14, height = 12)

  # Panel A: Heatmap of protein ranks (top-50 from cycloess, rows = proteins, cols = methods)
  for (cname in contrast_names) {
    ref_key <- paste("cycloess", cname, sep = ".")
    if (is.null(all_rankings[[ref_key]])) next
    top50_ids <- all_rankings[[ref_key]]$uniprot_id[1:50]
    top50_genes <- all_rankings[[ref_key]]$gene[1:50]

    heatmap_mat <- matrix(NA, nrow = 50, ncol = length(norm_methods),
                          dimnames = list(top50_genes, norm_methods))
    for (method in norm_methods) {
      key <- paste(method, cname, sep = ".")
      if (!is.null(all_rankings[[key]])) {
        idx <- match(top50_ids, all_rankings[[key]]$uniprot_id)
        heatmap_mat[, method] <- all_rankings[[key]]$rank[idx]
      }
    }

    heatmap(heatmap_mat, Colv = NA, scale = "none",
            col = colorRampPalette(c("darkblue", "white", "darkred"))(50),
            main = paste("Protein ranks across norms:", cname),
            margins = c(8, 10))
  }

  # Panel B: Rank-rank scatter (cycloess vs each other method, top 200)
  for (cname in contrast_names) {
    ref_key <- paste("cycloess", cname, sep = ".")
    if (is.null(all_rankings[[ref_key]])) next
    ref <- all_rankings[[ref_key]] %>% filter(rank <= 200)

    par(mfrow = c(2, 4))
    for (method in setdiff(norm_methods, "cycloess")) {
      key <- paste(method, cname, sep = ".")
      if (is.null(all_rankings[[key]])) next
      test <- all_rankings[[key]]
      merged_rr <- inner_join(
        ref %>% dplyr::select(uniprot_id, rank),
        test %>% dplyr::select(uniprot_id, rank),
        by = "uniprot_id", suffix = c(".cyc", ".test")
      )
      plot(merged_rr$rank.cyc, merged_rr$rank.test,
           xlab = "Rank (cycloess)", ylab = paste0("Rank (", method, ")"),
           main = paste(cname, "-", method),
           pch = 16, cex = 0.5, col = rgb(0, 0, 0, 0.3))
      abline(0, 1, col = "red", lty = 2)
    }
    par(mfrow = c(1, 1))
  }

  # Panel C: Bar chart of normalization-sensitive protein counts per contrast
  sens_counts <- norm_flags %>%
    filter(sensitive) %>%
    count(contrast)
  if (nrow(sens_counts) > 0) {
    barplot(sens_counts$n, names.arg = sens_counts$contrast,
            main = "Normalization-sensitive proteins per contrast",
            ylab = "Count", col = "#DAA520", las = 2)
  }

  # Panel D: Top-10 pathway concordance tile plot
  pathway_long <- list()
  for (method in norm_methods) {
    for (cname in contrast_names) {
      key <- paste(method, cname, sep = ".")
      if (!is.null(all_pathways[[key]])) {
        pathway_long <- c(pathway_long, list(
          all_pathways[[key]] %>% mutate(method = method, contrast = cname, present = 1)
        ))
      }
    }
  }
  if (length(pathway_long) > 0) {
    pw_df <- bind_rows(pathway_long)
    pw_tile <- ggplot(pw_df, aes(x = method, y = pathway, fill = present)) +
      geom_tile(color = "white") +
      facet_wrap(~ contrast, scales = "free_y") +
      scale_fill_gradient(low = "white", high = "#4A90D9", guide = "none") +
      labs(title = "Top-10 pathway concordance across normalization methods") +
      theme_bw(base_size = 8) +
      theme(axis.text.x = element_text(angle = 45, hjust = 1),
            axis.text.y = element_text(size = 6))
    print(pw_tile)
  }

  dev.off()
  cat("   Saved: norm_sensitivity_rankings.csv, norm_sensitivity_flags.csv, 04_norm_sensitivity.pdf\n")
}

# ==============================================================================
# 6.  SENSITIVITY POWER ANALYSIS
# ==============================================================================

cat("\n>> 6 -- Sensitivity Power Analysis\n")

# Load imputed data to estimate residual SD from limma fit
imp_data_path <- file.path(base_dir, "02_Imputation", "c_data", "01_imputed.csv")
df_imp <- read_csv(imp_data_path, show_col_types = FALSE)

ann_cols   <- c("uniprot_id", "protein", "gene", "description")
samp_names <- setdiff(names(df_imp), ann_cols)
mat_imp    <- as.matrix(df_imp[, samp_names])
rownames(mat_imp) <- df_imp$uniprot_id

# Build metadata
meta <- tibble(sample_id = samp_names) %>%
  mutate(
    prefix   = str_extract(sample_id, "^[A-Z]+"),
    subj_num = str_extract(sample_id, "S\\d+"),
    time     = str_extract(sample_id, "(Pre|Post)$"),
    age      = if_else(str_detect(prefix, "^O"), "Old", "Young"),
    subject  = paste0(prefix, "_", subj_num),
    group    = paste(age, time, sep = "_")
  )
meta$group <- factor(meta$group, levels = c("Young_Pre", "Young_Post", "Old_Pre", "Old_Post"))

# Fit limma to estimate sigma
design_mat <- model.matrix(~ 0 + group, data = meta)
colnames(design_mat) <- gsub("^group", "", colnames(design_mat))

corr_fit <- duplicateCorrelation(mat_imp, design_mat, block = meta$subject)
fit <- lmFit(mat_imp, design_mat, block = meta$subject,
             correlation = corr_fit$consensus)
fit <- eBayes(fit)

# Observed residual SD (median sigma)
median_sigma <- median(sqrt(fit$s2.post))
cat(sprintf("   Median residual SD (sigma): %.3f\n", median_sigma))

# Sensitivity power curves for n = 15, 16 (group sizes in this study)
effect_sizes <- seq(0.1, 2.0, by = 0.05)
power_rows <- list()
for (n in c(15, 16)) {
  for (d in effect_sizes) {
    pw <- pwr.t.test(n = n, d = d / median_sigma, sig.level = 0.05,
                     type = "two.sample", alternative = "two.sided")
    power_rows <- c(power_rows, list(tibble(
      n = n, effect_size_log2FC = d, cohen_d = d / median_sigma,
      power = pw$power
    )))
  }
}

power_df <- bind_rows(power_rows)
write_csv(power_df, file.path(DATA_DIR, "power_sensitivity.csv"))

p_power <- ggplot(power_df, aes(x = effect_size_log2FC, y = power,
                                 color = factor(n))) +
  geom_line(linewidth = 1) +
  geom_hline(yintercept = 0.8, linetype = "dashed", color = "grey50") +
  scale_color_manual(values = c("15" = "#D94A4A", "16" = "#4A90D9")) +
  labs(x = "Effect size (log2FC)", y = "Power",
       color = "n per group",
       title = "Sensitivity power curves",
       subtitle = sprintf("Median sigma = %.3f", median_sigma)) +
  theme_bw(base_size = 11)

ggsave(file.path(REPORT_DIR, "05_power_curves.pdf"), p_power,
       width = 7, height = 5)
cat("   Saved: power_sensitivity.csv, 05_power_curves.pdf\n")

# ==============================================================================
# 7.  CONFIDENT EFFECT SIZES (topconfects)
# ==============================================================================

cat("\n>> 7 -- Confident Effect Sizes (topconfects)\n")

# Build limma fit in-script (self-contained, no dependency on DEP workspace)
contrast_matrix <- makeContrasts(
  Training_Young = Young_Post - Young_Pre,
  Training_Old   = Old_Post - Old_Pre,
  Aging          = Old_Pre - Young_Pre,
  Interaction    = (Old_Post - Old_Pre) - (Young_Post - Young_Pre),
  levels = design_mat
)

fit_contrasts <- contrasts.fit(fit, contrast_matrix)
fit_contrasts <- eBayes(fit_contrasts)

confect_rows <- list()
for (cname in contrast_names) {
  confects <- limma_confects(fit_contrasts, coef = cname, fdr = 0.05)
  top_df <- as.data.frame(confects$table) %>%
    head(50) %>%
    mutate(contrast = cname)

  # Add gene names
  if ("name" %in% names(top_df)) {
    gene_lookup <- df_imp %>% dplyr::select(uniprot_id, gene)
    top_df <- top_df %>%
      left_join(gene_lookup, by = c("name" = "uniprot_id"))
  }

  confect_rows <- c(confect_rows, list(top_df))
}

confect_df <- bind_rows(confect_rows)
write_csv(confect_df, file.path(DATA_DIR, "confident_effect_sizes.csv"))

# Confect ranking plot
p_confect_list <- list()
for (cname in contrast_names) {
  cdata <- confect_df %>%
    filter(contrast == cname) %>%
    head(20)
  if (nrow(cdata) == 0) next

  # Use gene name if available, otherwise use 'name'
  label_col <- if ("gene" %in% names(cdata) && any(!is.na(cdata$gene))) "gene" else "name"

  cdata <- cdata %>%
    mutate(label = if (label_col == "gene") coalesce(gene, name) else name) %>%
    mutate(label = fct_reorder(label, confect, .na_rm = TRUE))

  p_confect_list[[cname]] <- ggplot(cdata, aes(x = confect, y = label)) +
    geom_point(size = 2) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "grey50") +
    labs(x = "Confident effect size (log2FC)", y = NULL, title = cname) +
    theme_bw(base_size = 9)
}

if (length(p_confect_list) > 0) {
  p_confect <- wrap_plots(p_confect_list, ncol = 2) +
    plot_annotation(title = "Top confident effect sizes (topconfects, FDR < 0.05)")
  ggsave(file.path(REPORT_DIR, "06_confect_ranking.pdf"), p_confect,
         width = 12, height = 10)
}
cat("   Saved: confident_effect_sizes.csv, 06_confect_ranking.pdf\n")

# ==============================================================================
# 8.  DEqMS CONCORDANCE (SUPPLEMENTARY)
# ==============================================================================

cat("\n>> 8 -- DEqMS Concordance (Supplementary)\n")

deqms_concord_rows <- list()
deqms_scatter_list <- list()

for (cname in contrast_names) {
  # Load limma (imputed) result
  limma_res <- limma_results[[paste("Imputed", cname, sep = ".")]]

  # Load DEqMS (imputed) result
  deqms_path <- file.path(base_dir, "03_DEP", "Imputed", "b_DEqMS", "c_data",
                           "per_contrast_results", paste0(cname, "_DEqMS.csv"))
  if (!file.exists(deqms_path)) {
    cat(sprintf("   WARNING: DEqMS result not found: %s\n", deqms_path))
    next
  }
  deqms_res <- read_csv(deqms_path, show_col_types = FALSE)

  # Merge
  merged <- inner_join(
    limma_res %>% dplyr::select(uniprot_id, logFC, P.Value, adj.P.Val),
    deqms_res %>% dplyr::select(uniprot_id, logFC, sca.P.Value, sca.adj.pval),
    by = "uniprot_id", suffix = c(".limma", ".deqms")
  )

  # Jaccard (FDR < 0.05)
  sig_limma <- merged$uniprot_id[merged$adj.P.Val < 0.05]
  sig_deqms <- merged$uniprot_id[merged$sca.adj.pval < 0.05]
  union_n   <- length(union(sig_limma, sig_deqms))
  intersect_n <- length(intersect(sig_limma, sig_deqms))
  jaccard <- if (union_n > 0) intersect_n / union_n else NA_real_

  deqms_concord_rows <- c(deqms_concord_rows, list(tibble(
    contrast = cname,
    n_sig_limma = length(sig_limma),
    n_sig_deqms = length(sig_deqms),
    n_sig_both = intersect_n,
    jaccard_FDR05 = round(jaccard, 3),
    logFC_cor = round(cor(merged$logFC.limma, merged$logFC.deqms, use = "complete.obs"), 4)
  )))

  # -log10(p) scatter
  deqms_scatter_list[[cname]] <- ggplot(merged, aes(
    x = -log10(P.Value), y = -log10(sca.P.Value))) +
    geom_point(alpha = 0.3, size = 0.8) +
    geom_abline(slope = 1, intercept = 0, color = "red", linetype = "dashed") +
    labs(x = "-log10(p) limma", y = "-log10(p) DEqMS", title = cname) +
    theme_bw(base_size = 10)
}

if (length(deqms_concord_rows) > 0) {
  deqms_concord <- bind_rows(deqms_concord_rows)
  write_csv(deqms_concord, file.path(DATA_DIR, "deqms_concordance.csv"))

  p_deqms <- wrap_plots(deqms_scatter_list, ncol = 2) +
    plot_annotation(title = "limma vs DEqMS p-value concordance (Imputed data)")
  ggsave(file.path(REPORT_DIR, "07_deqms_concordance.pdf"), p_deqms,
         width = 10, height = 8)
  cat("   Saved: deqms_concordance.csv, 07_deqms_concordance.pdf\n")
}

# ==============================================================================
# 9.  SUMMARY WORKBOOK
# ==============================================================================

cat("\n>> 9 -- Summary Workbook\n")

# Multi-panel summary PDF
pdf(file.path(REPORT_DIR, "08_DEP_assessment_summary.pdf"), width = 14, height = 10)

# Page 1: Threshold sensitivity + robustness
print(p_thresh)
print(p_robust)

# Page 2: logFC scatter + power
print(p_scatter)
print(p_power)

# Page 3: Confects
if (exists("p_confect")) print(p_confect)

# Page 4: DEqMS
if (exists("p_deqms")) print(p_deqms)

dev.off()
cat("   Saved: 08_DEP_assessment_summary.pdf\n")

# Multi-sheet Excel workbook
wb <- createWorkbook()
addWorksheet(wb, "Threshold_Sensitivity")
writeData(wb, "Threshold_Sensitivity", thresh_df %>% dplyr::select(-FDR_label, -LFC_label))

addWorksheet(wb, "Data_Concordance")
writeData(wb, "Data_Concordance", concord_df)

addWorksheet(wb, "Imputation_Flags")
writeData(wb, "Imputation_Flags", flag_df)

addWorksheet(wb, "Robustness")
writeData(wb, "Robustness", robust_df)

if (exists("rank_concord")) {
  addWorksheet(wb, "Norm_Sensitivity_Rankings")
  writeData(wb, "Norm_Sensitivity_Rankings", rank_concord)
}

if (exists("norm_flags")) {
  addWorksheet(wb, "Norm_Sensitivity_Flags")
  writeData(wb, "Norm_Sensitivity_Flags", norm_flags)
}

addWorksheet(wb, "Power_Sensitivity")
writeData(wb, "Power_Sensitivity", power_df)

addWorksheet(wb, "Confident_Effects")
writeData(wb, "Confident_Effects", confect_df)

if (exists("deqms_concord")) {
  addWorksheet(wb, "DEqMS_Concordance")
  writeData(wb, "DEqMS_Concordance", deqms_concord)
}

saveWorkbook(wb, file.path(DATA_DIR, "DEP_assessment_workbook.xlsx"), overwrite = TRUE)
cat("   Saved: DEP_assessment_workbook.xlsx\n")

# ==============================================================================
# 10.  SESSION INFO
# ==============================================================================

cat("\n>> 10 -- Session Info\n\n")
sessionInfo()

cat("\n=== YvO DEP Assessment Complete ===\n")
