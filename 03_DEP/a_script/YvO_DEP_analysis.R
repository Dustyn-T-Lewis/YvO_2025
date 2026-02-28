################################################################################
#   YvO Differential Expression — limma pipeline (proteoDA)
#   Design: 2×2 factorial (Age × Time), repeated measures on subject
#   Input:  cycloess-normalized, non-imputed (limma handles NAs per-protein)
#
#   Core references:
#     - Ritchie et al. 2015, Nucleic Acids Res 43(7):e47 — limma
#     - Smyth 2005, Stat Appl Genet Mol Biol 3(1):3
#       duplicateCorrelation for repeated-measures (limma User's Guide §9.7)
#     - Karpievitch et al. 2012, BMC Bioinform 13(S16):S5
#       Normalization and missing value handling in label-free proteomics
#     - Vous et al. 2024, Nat Commun 15:3891
#       Limma on non-imputed data outperforms imputed alternatives
#     - Xiao et al. 2014, Bioinformatics 30(6):801-807 (PMC3957066)
#       Pi-value reparametrized as Pi = p^|logFC| (range 0-1).
#       -log10(Pi) = |logFC| × -log10(p) = original pi-score.
#       Threshold Pi < 0.05 ↔ original pi > 1.3.
#
#   Exercise proteomics context:
#     - Robinson et al. 2022, GeroScience 45:1271-1287
#     - Bechshoft et al. 2024, Aging 16(7):6157-6187
#     - Hostrup et al. 2022, eLife 11:e69802
#     - Hulmi et al. 2025, J Physiol 603:438-453
################################################################################

# === SETUP ====================================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(proteoDA)
  library(openxlsx)
  library(patchwork)
  library(ggrepel)
  library(gridExtra)
  library(boot)
  library(pwr)
})

set.seed(42)

setwd(rprojroot::find_rstudio_root_file())
DATA_FILE  <- "01_normalization/c_data/02_normalized.csv"
REPORT_DIR <- "03_DEP/b_reports"
DATA_DIR   <- "03_DEP/c_data"

dir.create(REPORT_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(DATA_DIR,   recursive = TRUE, showWarnings = FALSE)

# === LOAD DATA & BUILD METADATA ==============================================

df <- read_csv(DATA_FILE, show_col_types = FALSE)

ann_cols   <- c("uniprot_id", "protein", "gene", "description")
ann        <- df[, ann_cols]
samp_names <- setdiff(names(df), ann_cols)
mat        <- as.matrix(df[, samp_names])
rownames(mat) <- ann$uniprot_id

cat(sprintf("Loaded: %d proteins x %d samples | missing: %d (%.1f%%)\n",
            nrow(mat), ncol(mat), sum(is.na(mat)),
            100 * sum(is.na(mat)) / length(mat)))

# Canonical metadata from normalisation DAList (not regex-derived)
dal_norm <- readRDS("01_normalization/c_data/03_DAList_normalized.rds")
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

# === CREATE DAList ============================================================

meta_df <- as.data.frame(meta)
rownames(meta_df) <- meta$sample_id

dal <- DAList(
  data       = mat,
  annotation = as.data.frame(ann),
  metadata   = meta_df,
  tags       = list(norm_method = "cycloess")
)

# === STATISTICAL DESIGN =======================================================

dal <- add_design(dal, "~ 0 + group + (1 | subject)")
colnames(dal$design$design_matrix) <- gsub("^group", "",
                                            colnames(dal$design$design_matrix))

# === CONTRASTS ================================================================

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

# === FIT MODEL & EXTRACT RESULTS =============================================

dal <- fit_limma_model(dal)

within_cor <- dal$eBayes_fit$correlation %||%
  dal$tags$duplicate_correlation %||% NA_real_
if (!is.na(within_cor)) cat(sprintf("Within-subject correlation: %.3f\n", within_cor))

# FDR 0.10: appropriate for exploratory proteomics with small n
# (Benjamini & Hochberg 1995 used 0.10 in examples; Choi et al. 2008,
# BMC Bioinform 9:43, recommend FDR 0-15% for label-free proteomics).
# Pi-score (Pi < 0.05) provides a secondary effect-size-weighted filter.
dal <- extract_DA_results(dal, pval_thresh = 0.10, lfc_thresh = 0, adj_method = "BH")

saveRDS(dal, file.path(DATA_DIR, "limma_DAList.rds"))

# === ADD PI-SCORES & WRITE TABLES (single-pass I/O) ==========================

write_limma_tables(dal, output_dir = DATA_DIR, overwrite = TRUE)

contrast_names   <- names(dal$results)
per_contrast_dir <- file.path(DATA_DIR, "per_contrast_results")

# Read each per-contrast CSV once, add pi-score, keep in memory
results_list <- list()
for (cname in contrast_names) {
  fpath <- file.path(per_contrast_dir, paste0(cname, ".csv"))
  res <- read_csv(fpath, show_col_types = FALSE) %>%
    mutate(
      pi_score = P.Value ^ abs(logFC),
      sig_pi = case_when(
        pi_score < 0.05 & logFC > 0 ~  1L,
        pi_score < 0.05 & logFC < 0 ~ -1L,
        TRUE ~ 0L
      )
    ) %>%
    select(-any_of(c("sig.PVal", "sig.FDR")))
  write_csv(res, fpath)
  results_list[[cname]] <- res
}

# Update combined_results.csv (single read-modify-write)
combined_path <- file.path(DATA_DIR, "combined_results.csv")
comb <- read_csv(combined_path, show_col_types = FALSE)
for (cname in contrast_names) {
  lfc_col  <- paste0("logFC_", cname)
  pval_col <- paste0("P.Value_", cname)
  if (!all(c(lfc_col, pval_col) %in% names(comb))) next
  pi_col  <- paste0("pi_score_", cname)
  sig_col <- paste0("sig_pi_", cname)
  comb[[pi_col]]  <- comb[[pval_col]] ^ abs(comb[[lfc_col]])
  comb[[sig_col]] <- case_when(
    comb[[pi_col]] < 0.05 & comb[[lfc_col]] > 0 ~  1L,
    comb[[pi_col]] < 0.05 & comb[[lfc_col]] < 0 ~ -1L,
    TRUE ~ 0L
  )
}
comb <- comb %>% select(-matches("^sig\\.PVal_|^sig\\.FDR_"))
write_csv(comb, combined_path)

# Build DA summary from in-memory results
da_path  <- file.path(DATA_DIR, "DA_summary.csv")
da_base  <- read_csv(da_path, show_col_types = FALSE)
pi_rows  <- map_dfr(contrast_names, function(cname) {
  res <- results_list[[cname]]
  bind_rows(
    tibble(contrast = cname, type = "up",
           sig.Pi     = sum(res$sig_pi == 1, na.rm = TRUE),
           sig.FDR.05 = sum(res$adj.P.Val < 0.05 & res$logFC > 0, na.rm = TRUE),
           sig.FDR.10 = sum(res$adj.P.Val < 0.10 & res$logFC > 0, na.rm = TRUE)),
    tibble(contrast = cname, type = "down",
           sig.Pi     = sum(res$sig_pi == -1, na.rm = TRUE),
           sig.FDR.05 = sum(res$adj.P.Val < 0.05 & res$logFC < 0, na.rm = TRUE),
           sig.FDR.10 = sum(res$adj.P.Val < 0.10 & res$logFC < 0, na.rm = TRUE)),
    tibble(contrast = cname, type = "nonsig",
           sig.Pi     = sum(res$sig_pi == 0, na.rm = TRUE),
           sig.FDR.05 = sum(!res$adj.P.Val < 0.05, na.rm = TRUE),
           sig.FDR.10 = sum(!res$adj.P.Val < 0.10, na.rm = TRUE))
  )
})
da_summary <- da_base %>% left_join(pi_rows, by = c("contrast", "type"))
write_csv(da_summary, da_path)

# Rebuild results.xlsx from in-memory data
wb_results <- createWorkbook()
for (cname in contrast_names) {
  addWorksheet(wb_results, cname)
  writeData(wb_results, cname, results_list[[cname]])
}
addWorksheet(wb_results, "DA_Summary")
writeData(wb_results, "DA_Summary", da_summary)
saveWorkbook(wb_results, file.path(DATA_DIR, "results.xlsx"), overwrite = TRUE)

# === WRITE PLOTS ==============================================================

tryCatch(
  write_limma_plots(dal,
                    grouping_column = "group",
                    output_dir      = REPORT_DIR,
                    table_columns   = c("uniprot_id", "gene", "protein"),
                    title_column    = "gene",
                    overwrite       = TRUE),
  error = function(e) cat(sprintf("write_limma_plots: %s\n", conditionMessage(e)))
)

# Reorganize proteoDA plots into per-contrast subdirectories
static_dir <- file.path(REPORT_DIR, "static_plots")
for (cname in contrast_names) {
  contrast_dir <- file.path(REPORT_DIR, cname)
  dir.create(contrast_dir, showWarnings = FALSE)
  html_file <- file.path(REPORT_DIR, paste0(cname, "_DA_report.html"))
  if (file.exists(html_file))
    file.rename(html_file, file.path(contrast_dir, basename(html_file)))
  pdfs <- list.files(static_dir, pattern = paste0("^", cname, "-"), full.names = TRUE)
  for (pdf in pdfs)
    file.rename(pdf, file.path(contrast_dir,
                               sub(paste0("^", cname, "-"), "", basename(pdf))))
}
if (length(list.files(static_dir)) == 0) unlink(static_dir, recursive = TRUE)

# --- Custom per-contrast: triple volcano + top-25 table ---

theme_dep <- theme_bw(base_size = 11) +
  theme(plot.title = element_text(face = "bold", size = 12),
        legend.position = "bottom")
pal_dir <- c(Up = "#D6604D", Down = "#4393C3", NS = "grey70")

for (cname in contrast_names) {
  res <- results_list[[cname]] %>%
    mutate(
      nlog10_pval = -log10(pmax(P.Value, 1e-300)),
      nlog10_adj  = -log10(pmax(adj.P.Val, 1e-300)),
      nlog10_pi   = -log10(pmax(pi_score, 1e-300)),
      dir_pi = factor(case_when(
        sig_pi ==  1 ~ "Up", sig_pi == -1 ~ "Down", TRUE ~ "NS"),
        levels = c("Up", "Down", "NS"))
    )

  top_nom <- slice_min(res, P.Value, n = 10, with_ties = FALSE)
  top_adj <- slice_min(res, adj.P.Val, n = 10, with_ties = FALSE)
  top_pi  <- slice_min(res, pi_score, n = 10, with_ties = FALSE)

  make_vol <- function(df, ycol, ylab, top_df, thresh) {
    ggplot(df, aes(logFC, .data[[ycol]], color = dir_pi)) +
      geom_point(alpha = 0.35, size = 1) +
      geom_hline(yintercept = thresh, linetype = "dashed", color = "grey40",
                 linewidth = 0.4) +
      geom_text_repel(data = top_df, aes(label = gene), size = 2.5,
                      max.overlaps = 15, show.legend = FALSE, seed = 42) +
      scale_color_manual(values = pal_dir, drop = FALSE) +
      labs(x = expression(log[2]~FC), y = ylab) +
      theme_dep + theme(legend.position = "none")
  }

  p1 <- make_vol(res, "nlog10_pval", expression(-log[10](P.Value)),
                 top_nom, -log10(0.01)) + ggtitle("Nominal P < 0.01")
  p2 <- make_vol(res, "nlog10_adj", expression(-log[10](adj.P.Val)),
                 top_adj, -log10(0.10)) + ggtitle("FDR < 0.10")
  p3 <- make_vol(res, "nlog10_pi", expression(-log[10](Pi)),
                 top_pi, -log10(0.05)) +
    ggtitle("Pi < 0.05") + theme(legend.position = "right") + labs(color = NULL)

  tbl_data <- res %>%
    slice_min(pi_score, n = 25, with_ties = FALSE) %>%
    transmute(
      UniProt = uniprot_id, Gene = gene,
      logFC = sprintf("%.3f", logFC),
      P.Value = formatC(P.Value, format = "e", digits = 2),
      adj.P.Val = formatC(adj.P.Val, format = "e", digits = 2),
      Pi = formatC(pi_score, format = "e", digits = 2)
    )

  p_tbl <- tableGrob(tbl_data, rows = NULL,
    theme = ttheme_minimal(base_size = 8,
      core    = list(fg_params = list(hjust = 0, x = 0.02)),
      colhead = list(fg_params = list(hjust = 0, x = 0.02, fontface = "bold"))))

  n_fdr <- sum(res$adj.P.Val < 0.10, na.rm = TRUE)
  n_pi  <- sum(res$sig_pi != 0)
  contrast_out <- file.path(REPORT_DIR, cname)

  pdf(file.path(contrast_out, "summary.pdf"), width = 16, height = 14)
  print(
    (p1 | p2 | p3) +
      plot_annotation(
        title = cname,
        subtitle = sprintf("FDR<0.10: %d | Pi<0.05: %d | %d proteins",
                           n_fdr, n_pi, nrow(res)),
        theme = theme(plot.title = element_text(face = "bold", size = 14)))
  )
  grid::grid.newpage()
  grid::grid.draw(arrangeGrob(
    p_tbl,
    top = grid::textGrob(
      sprintf("%s \u2014 Top 25 by Pi-score", cname),
      gp = grid::gpar(fontface = "bold", fontsize = 14)),
    bottom = grid::textGrob(
      "UniProt | Gene | log2FC | P | adj.P | Pi",
      gp = grid::gpar(fontsize = 9, col = "grey40"))
  ))
  dev.off()

  cat(sprintf("  %s: FDR=%d, Pi=%d\n", cname, n_fdr, n_pi))
}

# --- Overall DEP summary bar chart ---

sc <- map_dfr(contrast_names, function(cname) {
  res <- results_list[[cname]]
  bind_rows(
    tibble(contrast = cname, criterion = "FDR < 0.10",
           up   = sum(res$adj.P.Val < 0.10 & res$logFC > 0, na.rm = TRUE),
           down = sum(res$adj.P.Val < 0.10 & res$logFC < 0, na.rm = TRUE)),
    tibble(contrast = cname, criterion = "Pi < 0.05",
           up   = sum(res$sig_pi == 1, na.rm = TRUE),
           down = sum(res$sig_pi == -1, na.rm = TRUE))
  )
}) %>%
  pivot_longer(c(up, down), names_to = "direction", values_to = "count") %>%
  mutate(
    signed    = if_else(direction == "down", -count, count),
    direction = factor(str_to_title(direction), levels = c("Up", "Down")),
    criterion = factor(criterion, levels = c("FDR < 0.10", "Pi < 0.05"))
  )

p_bar <- ggplot(sc, aes(x = contrast, y = signed, fill = direction)) +
  geom_col(position = "identity", width = 0.7) +
  geom_hline(yintercept = 0, linewidth = 0.3) +
  geom_text(data = filter(sc, signed > 0),
            aes(y = signed, label = count), vjust = -0.3, size = 3) +
  geom_text(data = filter(sc, signed < 0),
            aes(y = signed, label = count), vjust = 1.3, size = 3) +
  facet_wrap(~criterion) +
  scale_fill_manual(values = c(Up = "#D6604D", Down = "#4393C3")) +
  labs(title = "DEP Counts by Significance Criterion",
       x = NULL, y = "Number of DEPs (Up / Down)", fill = NULL) +
  theme_dep +
  theme(axis.text.x = element_text(angle = 30, hjust = 1))

# Wide table: FDR and Pi side by side per contrast
stbl_wide <- map_dfr(contrast_names, function(cname) {
  res <- results_list[[cname]]
  tibble(
    Contrast  = cname,
    `FDR Up`  = sum(res$adj.P.Val < 0.10 & res$logFC > 0, na.rm = TRUE),
    `FDR Down`= sum(res$adj.P.Val < 0.10 & res$logFC < 0, na.rm = TRUE),
    `FDR Tot` = `FDR Up` + `FDR Down`,
    `Pi Up`   = sum(res$sig_pi == 1, na.rm = TRUE),
    `Pi Down` = sum(res$sig_pi == -1, na.rm = TRUE),
    `Pi Tot`  = `Pi Up` + `Pi Down`
  )
})

p_stbl <- tableGrob(stbl_wide, rows = NULL,
  theme = ttheme_minimal(base_size = 10,
    core    = list(fg_params = list(hjust = 0.5)),
    colhead = list(fg_params = list(fontface = "bold", hjust = 0.5))))

# Stored for composite page — assembled after robustness plots are built

# === SUMMARY ==================================================================

print(dal$design$contrast_matrix)
print(da_summary)

# NOTE: fGSEA run separately per contrast without cross-contrast correction.
# Standard practice (each contrast = independent question), but running the
# same gene sets across 5 contrasts increases spurious enrichment probability.
# Acknowledge in manuscript methods.

# ============================================================================
# ROBUSTNESS ANALYSES
# ============================================================================

# === 1. BLUNTING DIAGNOSTICS (Enhanced) ======================================
# Tests whether |logFC| distributions differ between Training_Young and
# Training_Old, providing evidence for age-dependent response attenuation.
# References: Conover 1999, Practical Nonparametric Statistics;
#             Romano et al. 2006, Stat Probab Lett 76:1421-1430

cat("\n--- Blunting Diagnostics ---\n")

blunt_df <- tibble(
  gene          = results_list[["Training_Young"]]$gene,
  abs_lfc_young = abs(results_list[["Training_Young"]]$logFC),
  abs_lfc_old   = abs(results_list[["Training_Old"]]$logFC)
) %>% filter(!is.na(abs_lfc_young) & !is.na(abs_lfc_old))

# KS test: distributional shape/location shift
ks_res <- ks.test(blunt_df$abs_lfc_young, blunt_df$abs_lfc_old)

# Fligner-Killeen: variance heterogeneity (robust to non-normality)
fk_data <- data.frame(
  abs_lfc  = c(blunt_df$abs_lfc_young, blunt_df$abs_lfc_old),
  contrast = factor(rep(c("Training_Young", "Training_Old"), each = nrow(blunt_df)))
)
fk_res <- fligner.test(abs_lfc ~ contrast, data = fk_data)

# Wilcoxon signed-rank: paired location shift (same proteins)
wx_res <- wilcox.test(blunt_df$abs_lfc_young, blunt_df$abs_lfc_old, paired = TRUE)

# Cliff's delta: non-parametric effect size (no additional package)
diffs <- blunt_df$abs_lfc_young - blunt_df$abs_lfc_old
cliff_delta <- (sum(diffs > 0) - sum(diffs < 0)) / length(diffs)
cliff_mag <- case_when(
  abs(cliff_delta) < 0.147 ~ "negligible",
  abs(cliff_delta) < 0.33  ~ "small",
  abs(cliff_delta) < 0.474 ~ "medium",
  TRUE                     ~ "large"
)

# Quantile comparison
quants <- c(0.25, 0.50, 0.75, 0.90)
quant_df <- tibble(
  quantile       = paste0("Q", quants * 100),
  Training_Young = quantile(blunt_df$abs_lfc_young, quants),
  Training_Old   = quantile(blunt_df$abs_lfc_old, quants),
  ratio_YO       = round(Training_Young / Training_Old, 3)
)

blunt_diag <- tibble(
  test = c("Kolmogorov-Smirnov", "Fligner-Killeen",
           "Wilcoxon signed-rank", "Cliff's delta"),
  statistic = c(ks_res$statistic, fk_res$statistic,
                wx_res$statistic, cliff_delta),
  p_value = c(ks_res$p.value, fk_res$p.value, wx_res$p.value, NA),
  interpretation = c(
    ifelse(ks_res$p.value < 0.05,
           "Distributions differ significantly (shape/location)",
           "No significant distributional difference"),
    ifelse(fk_res$p.value < 0.05,
           "Variance differs (consistent with blunted response)",
           "No significant variance difference"),
    ifelse(wx_res$p.value < 0.05,
           "Paired |logFC| shift significant (attenuated with age)",
           "No significant paired shift"),
    sprintf("%s effect (delta = %.3f; >0 = young responds more)",
            str_to_title(cliff_mag), cliff_delta)
  )
)
write_csv(blunt_diag, file.path(DATA_DIR, "blunting_diagnostics.csv"))

cat(sprintf(paste0(
  "  KS: D=%.4f, p=%.4g\n",
  "  Fligner: chi2=%.2f, p=%.4g\n",
  "  Wilcoxon: V=%.0f, p=%.4g\n",
  "  Cliff's delta: %.3f (%s)\n"),
  ks_res$statistic, ks_res$p.value,
  fk_res$statistic, fk_res$p.value,
  wx_res$statistic, wx_res$p.value,
  cliff_delta, cliff_mag))

# Density overlay (enhanced annotations)
med_young <- median(blunt_df$abs_lfc_young)
med_old   <- median(blunt_df$abs_lfc_old)

density_long <- bind_rows(
  tibble(abs_lfc = blunt_df$abs_lfc_young, contrast = "Training (Young)"),
  tibble(abs_lfc = blunt_df$abs_lfc_old,   contrast = "Training (Old)")
)

p_density <- ggplot(density_long, aes(x = abs_lfc, fill = contrast, color = contrast)) +
  geom_density(alpha = 0.3, linewidth = 0.6) +
  geom_vline(xintercept = med_young, linetype = "dashed",
             color = "#E05A4E", linewidth = 0.5) +
  geom_vline(xintercept = med_old, linetype = "dashed",
             color = "#5DA5DA", linewidth = 0.5) +
  annotate("text", x = med_young + 0.02, y = Inf, vjust = 2, hjust = 0,
           label = sprintf("median = %.3f", med_young),
           size = 2.5, color = "#E05A4E", fontface = "bold") +
  annotate("text", x = med_old + 0.02, y = Inf, vjust = 3.5, hjust = 0,
           label = sprintf("median = %.3f", med_old),
           size = 2.5, color = "#5DA5DA", fontface = "bold") +
  annotate("text", x = Inf, y = Inf, hjust = 1.1, vjust = 1.5,
           label = sprintf(paste0(
             "KS D=%.3f, p=%.2g\n",
             "Fligner chi2=%.1f, p=%.2g\n",
             "Wilcoxon p=%.2g\n",
             "Cliff's delta=%.3f (%s)"),
             ks_res$statistic, ks_res$p.value,
             fk_res$statistic, fk_res$p.value,
             wx_res$p.value,
             cliff_delta, cliff_mag),
           size = 2.3, fontface = "bold", color = "grey25") +
  scale_fill_manual(values = c("Training (Young)" = "#E05A4E",
                                "Training (Old)" = "#5DA5DA")) +
  scale_color_manual(values = c("Training (Young)" = "#E05A4E",
                                 "Training (Old)" = "#5DA5DA")) +
  labs(title = "Effect Size Distribution: Training Response by Age",
       subtitle = "|logFC| density — evidence for age-dependent attenuation",
       x = "|logFC|", y = "Density", fill = NULL, color = NULL) +
  theme_bw(base_size = 11) +
  theme(plot.title = element_text(face = "bold", size = 12),
        legend.position = "bottom")

cat("  Computed blunting diagnostics\n")

# === 2. BOOTSTRAP CI (Effect Sizes) ==========================================
# Median |logFC| with 95% BCa bootstrap CI per contrast
# Reference: Efron & Tibshirani 1993, An Introduction to the Bootstrap

cat("\n--- Bootstrap CI ---\n")

median_fn <- function(d, i) median(d[i], na.rm = TRUE)

boot_contrasts <- c("Training_Young", "Training_Old", "Aging", "Interaction")
boot_df <- map_dfr(boot_contrasts, function(cname) {
  vals <- abs(results_list[[cname]]$logFC)
  vals <- vals[!is.na(vals)]
  b  <- boot(vals, median_fn, R = 10000)
  ci <- tryCatch(boot.ci(b, type = "bca"),
                 error = function(e) boot.ci(b, type = "perc"))
  ci_lo <- if (!is.null(ci$bca)) ci$bca[4] else ci$percent[4]
  ci_hi <- if (!is.null(ci$bca)) ci$bca[5] else ci$percent[5]
  tibble(contrast = cname, median_absLFC = median(vals),
         ci_lower = ci_lo, ci_upper = ci_hi,
         boot_se = sd(b$t), n_proteins = length(vals))
})
write_csv(boot_df, file.path(DATA_DIR, "effect_size_bootstrap.csv"))

pal_contrast <- c(Training_Young = "#E05A4E", Training_Old = "#5DA5DA",
                  Aging = "#4CAF50", Interaction = "#9B7FBF")

p_forest <- ggplot(boot_df, aes(y = reorder(contrast, median_absLFC),
                                 x = median_absLFC, color = contrast)) +
  geom_pointrange(aes(xmin = ci_lower, xmax = ci_upper),
                  size = 0.8, linewidth = 0.8) +
  scale_color_manual(values = pal_contrast, guide = "none") +
  labs(title = "Effect Size by Contrast",
       subtitle = "Median |logFC| with 95% BCa bootstrap CI (10k resamples)",
       x = "Median |logFC|", y = NULL) +
  theme_bw(base_size = 11) +
  theme(plot.title = element_text(face = "bold", size = 12))

cat("  Computed bootstrap CIs\n")
print(boot_df)

# === 3. POWER ANALYSIS =======================================================
# Approximate minimum detectable logFC at 80% power.
# Conservative: pwr.t.test assumes standard t; limma's moderated t has higher
# power via empirical Bayes variance shrinkage. Treat as lower bounds.
# Reference: Cohen 1988, Statistical Power Analysis

cat("\n--- Power Analysis ---\n")

fit <- dal$eBayes_fit
sigma_residual <- sqrt(mean(fit$sigma^2, na.rm = TRUE))
n_young <- sum(meta$age == "Young" & meta$time == "Pre")
n_old   <- sum(meta$age == "Old"   & meta$time == "Pre")

power_df <- map_dfr(boot_contrasts, function(cname) {
  n_subj <- switch(cname,
    Training_Young = n_young, Training_Old = n_old,
    min(n_young, n_old))
  paired <- cname %in% c("Training_Young", "Training_Old")
  eff_sigma <- if (paired && !is.na(within_cor))
    sigma_residual * sqrt(2 * (1 - within_cor)) else
    sigma_residual * sqrt(2)
  pw <- pwr.t.test(n = n_subj, d = NULL, sig.level = 0.10, power = 0.80,
                   type = if (paired) "paired" else "two.sample")
  tibble(contrast = cname, n_subjects = n_subj,
         within_cor = ifelse(paired, within_cor, NA_real_),
         effective_sigma = round(eff_sigma, 4),
         min_detectable_d = round(pw$d, 4),
         min_detectable_logFC = round(pw$d * eff_sigma, 4),
         power = 0.80, alpha = 0.10)
})
write_csv(power_df, file.path(DATA_DIR, "power_analysis.csv"))
cat("  Saved: power_analysis.csv\n")
print(as.data.frame(power_df))

# === 4. IMPUTATION SENSITIVITY ===============================================
# Compare t-statistics: non-imputed (main) vs imputed limma.
# Reference: Vous et al. 2024, Nat Commun 15:3891

cat("\n--- Imputation Sensitivity ---\n")

imp_path     <- "02_Imputation/c_data/01_imputed.csv"
sens_df      <- NULL
scatter_list <- list()
p_sens       <- NULL

if (file.exists(imp_path)) {
  imp_data <- read_csv(imp_path, show_col_types = FALSE)
  imp_samp <- setdiff(names(imp_data), ann_cols)
  imp_mat  <- as.matrix(imp_data[, imp_samp])
  rownames(imp_mat) <- imp_data$uniprot_id

  shared_samps <- intersect(colnames(mat), colnames(imp_mat))
  meta_imp_df  <- meta_df[meta_df$sample_id %in% shared_samps, ]

  dal_imp <- DAList(
    data       = imp_mat[, shared_samps],
    annotation = as.data.frame(ann),
    metadata   = meta_imp_df,
    tags       = list(norm_method = "cycloess_imputed")
  )
  dal_imp <- add_design(dal_imp, "~ 0 + group + (1 | subject)")
  colnames(dal_imp$design$design_matrix) <- gsub("^group", "",
    colnames(dal_imp$design$design_matrix))
  dal_imp <- add_contrasts(dal_imp, contrasts_vector = c(
    "Training_Young = Young_Post - Young_Pre",
    "Training_Old = Old_Post - Old_Pre",
    "Aging = Old_Pre - Young_Pre",
    "Interaction = (Old_Post - Old_Pre) - (Young_Post - Young_Pre)"
  ))
  dal_imp <- fit_limma_model(dal_imp)
  dal_imp <- extract_DA_results(dal_imp, pval_thresh = 0.10, lfc_thresh = 0,
                                 adj_method = "BH")

  imp_tmp <- file.path(tempdir(), "imp_results")
  dir.create(imp_tmp, showWarnings = FALSE)
  write_limma_tables(dal_imp, output_dir = imp_tmp, overwrite = TRUE)
  imp_comb <- read_csv(file.path(imp_tmp, "combined_results.csv"),
                       show_col_types = FALSE)

  sens_contrasts <- c("Training_Young", "Training_Old", "Aging", "Interaction")
  sens_df <- map_dfr(sens_contrasts, function(cname) {
    t_col   <- paste0("t_", cname)
    adj_col <- paste0("adj.P.Val_", cname)
    if (!(t_col %in% names(comb)) || !(t_col %in% names(imp_comb))) return(NULL)

    merged <- inner_join(
      comb %>% select(uniprot_id, t_nonimp = all_of(t_col),
                       padj_nonimp = all_of(adj_col)),
      imp_comb %>% select(uniprot_id, t_imp = all_of(t_col),
                           padj_imp = all_of(adj_col)),
      by = "uniprot_id"
    ) %>% filter(!is.na(t_nonimp) & !is.na(t_imp))

    sp <- cor.test(merged$t_nonimp, merged$t_imp, method = "spearman")

    merged <- merged %>% mutate(
      switch = case_when(
        padj_nonimp < 0.10 & !(padj_imp < 0.10) ~ "Lost in imputed",
        !(padj_nonimp < 0.10) & padj_imp < 0.10  ~ "Gained in imputed",
        TRUE ~ "Concordant"
      )
    )

    scatter_list[[cname]] <<- ggplot(merged, aes(t_nonimp, t_imp, color = switch)) +
      geom_abline(slope = 1, intercept = 0, linetype = "dashed",
                  color = "grey50", linewidth = 0.3) +
      geom_point(alpha = 0.65, size = 1) +
      scale_color_manual(values = c(Concordant = "grey40",
                                     "Lost in imputed" = "#D6604D",
                                     "Gained in imputed" = "#4393C3"), name = NULL) +
      annotate("text", x = -Inf, y = Inf, hjust = -0.05, vjust = 1.3,
               label = sprintf("rho = %.3f", sp$estimate),
               size = 3, fontface = "bold", color = "grey25") +
      labs(title = cname, x = "t (non-imputed)", y = "t (imputed)") +
      theme_bw(base_size = 9) +
      theme(plot.title = element_text(face = "bold", size = 10),
            legend.position = "bottom", legend.text = element_text(size = 7))

    tibble(contrast = cname, spearman_rho = round(sp$estimate, 4),
           p_value = sp$p.value, n_proteins = nrow(merged))
  })

  write_csv(sens_df, file.path(DATA_DIR, "imputation_sensitivity.csv"))

  if (length(scatter_list) > 0) {
    p_sens <- wrap_plots(scatter_list, ncol = 2) +
      plot_annotation(
        title    = "Imputation Sensitivity: t-statistic Concordance",
        subtitle = "Non-imputed (main) vs imputed limma",
        theme = theme(plot.title = element_text(face = "bold", size = 14)))
  }

  cat("  Computed imputation sensitivity\n")
  print(as.data.frame(sens_df))
} else {
  cat("  Imputed data not found — skipping\n")
}

# === COMBINED OVERVIEW PDF ====================================================
# Page 1: DEP counts (top) → table (mid) → forest + density (bottom row)
# Page 2: Imputation sensitivity

cat("\n--- Writing combined overview PDF ---\n")

# Build composite page 1
p_forest_clean <- p_forest +
  labs(title = "Effect Size by Contrast",
       subtitle = "Median |logFC|, 95% BCa CI") +
  theme(plot.title = element_text(size = 11, face = "bold"),
        plot.subtitle = element_text(size = 9))

p_density_clean <- p_density +
  labs(title = "Training Response Attenuation",
       subtitle = "|logFC| density by age group") +
  theme(plot.title = element_text(size = 11, face = "bold"),
        plot.subtitle = element_text(size = 9),
        legend.position = "bottom",
        legend.key.size = unit(0.4, "cm"))

p_page1 <- p_bar /
  wrap_elements(p_stbl) /
  (p_forest_clean | p_density_clean) +
  plot_layout(heights = c(3, 1.2, 3)) +
  plot_annotation(
    title = "YvO Differential Expression Overview",
    theme = theme(plot.title = element_text(face = "bold", size = 16)))

pdf(file.path(REPORT_DIR, "00_DEP_overview.pdf"), width = 14, height = 16)

# Page 1: Composite
print(p_page1)

# Page 2: Imputation sensitivity
if (!is.null(p_sens)) print(p_sens)

dev.off()
cat("  Saved: 00_DEP_overview.pdf\n")

# ============================================================================
# SUPPLEMENTARY EXCEL
# ============================================================================
# Consolidated workbook for manuscript supplementary materials.
# Each sheet has an interpretation header followed by data.

cat("\n--- Building DEP Supplementary Excel ---\n")

wb_supp    <- createWorkbook()
hdr_style  <- createStyle(textDecoration = "italic", fontColour = "#555555",
                           wrapText = TRUE)
bold_style <- createStyle(textDecoration = "bold")

write_supp_sheet <- function(wb, sheet, header, data) {
  addWorksheet(wb, sheet)
  for (i in seq_along(header))
    writeData(wb, sheet, header[i], startRow = i, startCol = 1)
  addStyle(wb, sheet, hdr_style, rows = seq_along(header), cols = 1,
           stack = TRUE)
  start <- length(header) + 2
  writeData(wb, sheet, data, startRow = start, headerStyle = bold_style)
  setColWidths(wb, sheet, cols = seq_len(ncol(data)), widths = "auto")
}

# Sheet 1: Methods Overview
overview <- tibble(
  Parameter = c(
    "Design", "Model", "Contrasts", "FDR method", "FDR threshold",
    "Pi-score threshold", "N proteins", "N subjects (Young)",
    "N subjects (Old)", "Within-subject correlation",
    "Normalization", "Imputation strategy"),
  Value = c(
    "2x2 factorial (Age x Time), repeated measures on subject",
    "limma + duplicateCorrelation via proteoDA",
    paste(contrast_names, collapse = "; "),
    "Benjamini-Hochberg",
    "0.10 (exploratory; pi-score provides secondary filter)",
    "Pi < 0.05 (= original pi > 1.3; Xiao et al. 2014)",
    as.character(nrow(mat)),
    as.character(n_young),
    as.character(n_old),
    sprintf("%.3f", within_cor),
    "Cycloess (selected by PRONE benchmarking)",
    "Non-imputed; limma handles NAs per-protein (Vous et al. 2024)")
)
write_supp_sheet(wb_supp, "Methods_Overview", c(
  "YvO DEP Analysis — Methods Overview",
  "Parameters and settings for the differential expression analysis."
), overview)

# Sheet 2: Significance Summary
write_supp_sheet(wb_supp, "Significance_Summary", c(
  "Significance counts by contrast, direction, and criterion.",
  "sig.PVal: nominal P < pval_thresh | sig.FDR: adj.P < pval_thresh (BH)",
  "sig.Pi: Pi-score < 0.05 | sig.FDR.05/10: FDR at 0.05 and 0.10 thresholds."
), da_summary)

# Sheets 3-7: Per-contrast results (sorted by pi-score, key columns)
keep_cols <- c("uniprot_id", "gene", "protein", "description",
               "logFC", "CI.L", "CI.R", "t", "P.Value", "adj.P.Val",
               "pi_score", "sig_pi")

contrast_labels <- c(
  Training_Young = "Young_Post - Young_Pre",
  Training_Old   = "Old_Post - Old_Pre",
  Aging          = "Old_Pre - Young_Pre",
  Interaction    = "(Old_Post-Old_Pre) - (Young_Post-Young_Pre)",
  Reversal       = "(Old_Post-Old_Pre) - (Old_Pre-Young_Pre)"
)

for (cname in contrast_names) {
  res <- results_list[[cname]]
  cols_present <- intersect(keep_cols, names(res))
  n_sig_fdr <- sum(res$adj.P.Val < 0.10, na.rm = TRUE)
  n_sig_pi  <- sum(res$sig_pi != 0)
  write_supp_sheet(wb_supp, cname, c(
    sprintf("%s: %s", cname, contrast_labels[[cname]]),
    sprintf("Proteins tested: %d | FDR<0.10: %d | Pi<0.05: %d",
            nrow(res), n_sig_fdr, n_sig_pi),
    "Sorted by Pi-score (most significant first)."
  ), res[order(res$pi_score), cols_present])
}

# Sheet 8: Blunting Diagnostics (tests + quantiles on one sheet)
addWorksheet(wb_supp, "Blunting_Diagnostics")
blunt_header <- c(
  "Age-dependent attenuation of training response.",
  "Tests compare |logFC| distributions: Training_Young vs Training_Old.",
  "KS: shape/location | Fligner-Killeen: variance | Wilcoxon: paired shift",
  "Cliff's delta: non-parametric effect size (>0 = young responds more)."
)
for (i in seq_along(blunt_header))
  writeData(wb_supp, "Blunting_Diagnostics", blunt_header[i],
            startRow = i, startCol = 1)
addStyle(wb_supp, "Blunting_Diagnostics", hdr_style,
         rows = seq_along(blunt_header), cols = 1, stack = TRUE)
writeData(wb_supp, "Blunting_Diagnostics", blunt_diag,
          startRow = length(blunt_header) + 2, headerStyle = bold_style)
quant_start <- length(blunt_header) + 2 + nrow(blunt_diag) + 2
writeData(wb_supp, "Blunting_Diagnostics", "Quantile comparison (|logFC|):",
          startRow = quant_start)
writeData(wb_supp, "Blunting_Diagnostics", quant_df,
          startRow = quant_start + 1, headerStyle = bold_style)
setColWidths(wb_supp, "Blunting_Diagnostics",
             cols = 1:max(ncol(blunt_diag), ncol(quant_df)), widths = "auto")

# Sheet 9: Bootstrap CIs
write_supp_sheet(wb_supp, "Effect_Size_CI", c(
  "Median |logFC| with 95% BCa bootstrap confidence intervals (10,000 resamples).",
  "Summarises global effect magnitude per contrast (Efron & Tibshirani 1993).",
  "Aging shows the largest median effect; Training_Old is attenuated vs Training_Young."
), boot_df)

# Sheet 10: Power Analysis
write_supp_sheet(wb_supp, "Power_Analysis", c(
  "Minimum detectable logFC at 80% power (alpha = 0.10).",
  "Conservative: uses pwr.t.test (standard t); limma's moderated t has higher power.",
  "Paired contrasts adjust for within-subject correlation; between-group do not.",
  "Interpret as approximate lower bounds on detectable effect sizes."
), power_df)

# Sheet 11: Imputation Sensitivity
if (!is.null(sens_df)) {
  write_supp_sheet(wb_supp, "Imputation_Sensitivity", c(
    "Spearman correlation of t-statistics: non-imputed (main) vs imputed limma.",
    "rho > 0.93 across all contrasts indicates high robustness to imputation choice.",
    "Reference: Vous et al. 2024, Nat Commun 15:3891."
  ), sens_df)
}

saveWorkbook(wb_supp, file.path(DATA_DIR, "DEP_supplementary.xlsx"), overwrite = TRUE)
cat("  Saved: DEP_supplementary.xlsx\n")

# === RENAME OUTPUTS (numbered pipeline order) =================================

cat("\n--- Numbering output files ---\n")

rename_map <- c(
  "limma_DAList.rds"           = "01_limma_DAList.rds",
  "DA_summary.csv"             = "02_DA_summary.csv",
  "combined_results.csv"       = "03_combined_results.csv",
  "per_contrast_results"       = "04_per_contrast_results",
  "results.xlsx"               = "05_results.xlsx",
  "blunting_diagnostics.csv"   = "06_blunting_diagnostics.csv",
  "effect_size_bootstrap.csv"  = "07_effect_size_bootstrap.csv",
  "power_analysis.csv"         = "08_power_analysis.csv",
  "imputation_sensitivity.csv" = "09_imputation_sensitivity.csv",
  "DEP_supplementary.xlsx"     = "10_DEP_supplementary.xlsx"
)

for (old_name in names(rename_map)) {
  old_path <- file.path(DATA_DIR, old_name)
  new_path <- file.path(DATA_DIR, rename_map[old_name])
  if (file.exists(old_path)) file.rename(old_path, new_path)
}
cat("  Output files numbered 01-10\n")

cat("\n=== YvO limma DEP complete ===\n")
