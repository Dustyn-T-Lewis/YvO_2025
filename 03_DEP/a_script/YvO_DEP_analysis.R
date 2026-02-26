################################################################################
#   YvO Differential Expression — limma pipeline (proteoDA)
#   Design: 2x2 factorial (Age x Time), repeated measures on subject
#   Input:  cycloess-normalized, non-imputed (limma handles NAs per-protein)
#
#   Methodological references:
#     - Ritchie et al. 2015, Nucleic Acids Res 43(7):e47 — limma
#     - Smyth 2005, Stat Appl Genet Mol Biol 3(1):3
#       duplicateCorrelation for multi-level/repeated-measures designs
#       (limma User's Guide Section 9.7)
#     - Karpievitch et al. 2019, BMC Bioinform 20:391
#       Limma on non-imputed data avoids imputation artifacts
#     - Xiao et al. 2014, Bioinformatics 30(6):801-807 (PMC3957066)
#       Pi-value = p^|logFC| (equivalent to -log10(Pi) = |logFC|*-log10(p))
#       Threshold Pi < 0.05 corresponds to original pi > 1.3
#
#   Exercise proteomics context:
#     - Robinson et al. 2022, GeroScience 45:1271-1287
#       Proteomic features of skeletal muscle adaptation to RET as a function
#       of age; 2x2 factorial (young/old x pre/post 20wk RET)
#     - Bechshoft et al. 2024, Aging 16(7):6157-6187
#       Deep proteomics of aging + resistance training in skeletal muscle
#     - Hostrup et al. 2022, eLife 11:e69802
#       HIT remodels skeletal muscle proteome; limma-based DE
#     - Hulmi et al. 2025, J Physiol 603:438-453
#       dia-PASEF proteomics with repeated measures training design;
#       150 altered proteins during first training period
################################################################################

# === SETUP ====================================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(proteoDA)
  library(openxlsx)
  library(patchwork)
  library(ggrepel)
  library(gridExtra)
})

set.seed(42)

setwd(rprojroot::find_rstudio_root_file())
DATA_FILE  <- "01_normalization/c_data/01_normalized.csv"
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

# Load canonical metadata from normalisation DAList (not regex-derived)
dal_norm <- readRDS("01_normalization/c_data/01_DAList_normalized.rds")
dal_meta <- as.data.frame(dal_norm$metadata)
meta <- tibble(
  sample_id = dal_meta$Col_ID,
  age       = dal_meta$Group,
  time      = dal_meta$Timepoint,
  group     = dal_meta$Group_Time,
  # Use Col_ID prefix (not Subject_ID) as subject key:
  # Subject_ID is NOT globally unique — 8 IDs are shared between Young/Old
  # cohorts (different individuals). The Col_ID prefix uniquely identifies

  # each person: e.g., "OP_S06" ≠ "O_S06" ≠ "Y_S07" ≠ "O_S07".
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

# Clean column names for readable contrast definitions
colnames(dal$design$design_matrix) <- gsub("^group", "", colnames(dal$design$design_matrix))

# === CONTRASTS ================================================================

# MANUSCRIPT NOTE: The Reversal contrast = (Old_Post - Old_Pre) - (Old_Pre - Young_Pre)
# is mathematically valid but uncommon. Justify clearly in the methods section.
#
# MANUSCRIPT NOTE: Interaction yields only 1 FDR-significant protein (n ~ 15/group).
# Acknowledge the inherent power limitation for interaction effects.
# Training_Old yields 0 FDR-significant proteins — this is NOT a bug but reflects
# attenuated training response with age (Robinson et al. 2022; Bechshoft et al. 2024).

dal <- add_contrasts(dal, contrasts_vector = c(
  "Training_Young = Young_Post - Young_Pre",
  "Training_Old = Old_Post - Old_Pre",
  "Aging = Old_Pre - Young_Pre",
  "Interaction = (Old_Post - Old_Pre) - (Young_Post - Young_Pre)",
  "Reversal = (Old_Post - Old_Pre) - (Old_Pre - Young_Pre)"
))

# === FIT MODEL & EXTRACT RESULTS =============================================

dal <- fit_limma_model(dal)

if (!is.null(dal$eBayes_fit$correlation)) {
  cat(sprintf("Within-subject correlation: %.3f\n", dal$eBayes_fit$correlation))
} else if (!is.null(dal$tags$duplicate_correlation)) {
  cat(sprintf("Within-subject correlation: %.3f\n", dal$tags$duplicate_correlation))
}

# NOTE: FDR threshold set to 0.10 (rather than 0.05) because:
#   (a) small sample size (~15/group) limits statistical power
#   (b) pi-score (Xiao et al. 2014) provides secondary effect-size filter
#   (c) consistent with exploratory proteomics discovery (Tyanova et al. 2016)
dal <- extract_DA_results(dal, pval_thresh = 0.10, lfc_thresh = 0, adj_method = "BH")

# Save DAList for downstream figure scripts (F2B/F3B complex analysis)
saveRDS(dal, file.path(DATA_DIR, "limma_DAList.rds"))

# === WRITE TABLES =============================================================

write_limma_tables(dal, output_dir = DATA_DIR, overwrite = TRUE)

# --- Add pi-scores, direction flags & update outputs -------------------------

contrast_names   <- names(dal$results)
per_contrast_dir <- file.path(DATA_DIR, "per_contrast_results")

# Per-contrast CSVs: add Pi-value + sig_pi + sig_nominal_01, drop old sig cols
# Pi-value (Xiao et al. 2014, Bioinformatics 30(6):801-807):
#   Pi = p^|logFC|  (range 0-1; smaller = more significant)
#   Equivalent to: -log10(Pi) = |logFC| * -log10(p) = original pi-score
#   Threshold Pi < 0.05 <==> pi-score > 1.3 (recommended in Xiao et al.)
for (cname in contrast_names) {
  fpath <- file.path(per_contrast_dir, paste0(cname, ".csv"))
  if (file.exists(fpath)) {
    res <- read_csv(fpath, show_col_types = FALSE) %>%
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

# Combined results: add pi_score + sig_pi per contrast
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
}

# DA summary: add sig.Pi and sig.Nominal.01 counts
da_path <- file.path(DATA_DIR, "DA_summary.csv")
if (file.exists(da_path)) {
  da <- read_csv(da_path, show_col_types = FALSE)
  pi_rows <- list()
  for (cname in contrast_names) {
    fpath <- file.path(per_contrast_dir, paste0(cname, ".csv"))
    if (file.exists(fpath)) {
      res <- read_csv(fpath, show_col_types = FALSE)
      pi_rows <- c(pi_rows, list(
        tibble(contrast = cname, type = "up",
               sig.Pi = sum(res$sig_pi == 1, na.rm = TRUE),
               sig.Nominal.01 = sum(res$sig_nominal_01 == 1, na.rm = TRUE)),
        tibble(contrast = cname, type = "down",
               sig.Pi = sum(res$sig_pi == -1, na.rm = TRUE),
               sig.Nominal.01 = sum(res$sig_nominal_01 == -1, na.rm = TRUE)),
        tibble(contrast = cname, type = "nonsig",
               sig.Pi = sum(res$sig_pi == 0, na.rm = TRUE),
               sig.Nominal.01 = sum(res$sig_nominal_01 == 0, na.rm = TRUE))
      ))
    }
  }
  da <- da %>% left_join(bind_rows(pi_rows), by = c("contrast", "type"))
  write_csv(da, da_path)
}

# Rebuild results.xlsx (one sheet per contrast + summary)
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
  addWorksheet(wb, "DA_Summary")
  writeData(wb, "DA_Summary", read_csv(da_path, show_col_types = FALSE))
}
saveWorkbook(wb, xlsx_path, overwrite = TRUE)

# === WRITE PLOTS ==============================================================

tryCatch(
  write_limma_plots(dal,
                    grouping_column = "group",
                    output_dir      = REPORT_DIR,
                    table_columns   = c("uniprot_id", "gene", "protein"),
                    title_column    = "gene",
                    overwrite       = TRUE),
  error = function(e) cat(sprintf("write_limma_plots warning: %s\n", conditionMessage(e)))
)

# Reorganize into per-contrast subdirectories
static_dir <- file.path(REPORT_DIR, "static_plots")

for (cname in contrast_names) {
  contrast_dir <- file.path(REPORT_DIR, cname)
  dir.create(contrast_dir, showWarnings = FALSE)

  html_file <- file.path(REPORT_DIR, paste0(cname, "_DA_report.html"))
  if (file.exists(html_file)) {
    file.rename(html_file, file.path(contrast_dir, basename(html_file)))
  }

  pdfs <- list.files(static_dir, pattern = paste0("^", cname, "-"), full.names = TRUE)
  for (pdf in pdfs) {
    clean_name <- sub(paste0("^", cname, "-"), "", basename(pdf))
    file.rename(pdf, file.path(contrast_dir, clean_name))
  }
}

if (length(list.files(static_dir)) == 0) unlink(static_dir, recursive = TRUE)

# === CUSTOM SUMMARY FIGURES ===================================================

theme_dep <- theme_bw(base_size = 11) +
  theme(plot.title = element_text(face = "bold", size = 12),
        legend.position = "bottom")
pal_dir <- c(Up = "#D6604D", Down = "#4393C3", NS = "grey70")

# --- Per-contrast: triple volcano + top hits table ---

for (cname in contrast_names) {
  fpath <- file.path(per_contrast_dir, paste0(cname, ".csv"))
  if (!file.exists(fpath)) next

  res <- read_csv(fpath, show_col_types = FALSE) %>%
    mutate(
      nlog10_pval = -log10(pmax(P.Value, 1e-300)),
      nlog10_adj  = -log10(pmax(adj.P.Val, 1e-300)),
      nlog10_pi   = -log10(pmax(pi_score, 1e-300)),
      dir_pi = factor(case_when(
        sig_pi ==  1 ~ "Up", sig_pi == -1 ~ "Down", TRUE ~ "NS"),
        levels = c("Up", "Down", "NS"))
    )

  top_nom <- res %>% arrange(P.Value) %>% slice_head(n = 10)
  top_adj <- res %>% arrange(adj.P.Val) %>% slice_head(n = 10)
  top_pi  <- res %>% arrange(pi_score) %>% slice_head(n = 10)

  make_vol <- function(df, ycol, ylab, top_df, thresh) {
    ggplot(df, aes(logFC, .data[[ycol]], color = dir_pi)) +
      geom_point(alpha = 0.35, size = 1) +
      geom_hline(yintercept = thresh, linetype = "dashed", color = "grey40",
                 linewidth = 0.4) +
      ggrepel::geom_text_repel(
        data = top_df, aes(label = gene), size = 2.5,
        max.overlaps = 15, show.legend = FALSE, seed = 42) +
      scale_color_manual(values = pal_dir, drop = FALSE) +
      labs(x = expression(log[2]~FC), y = ylab) +
      theme_dep + theme(legend.position = "none")
  }

  p1 <- make_vol(res, "nlog10_pval", expression(-log[10](P.Value)),
                 top_nom, -log10(0.01)) +
    ggtitle("Nominal P-value (P < 0.01)")

  p2 <- make_vol(res, "nlog10_adj", expression(-log[10](adj.P.Val)),
                 top_adj, -log10(0.10)) +
    ggtitle("Adjusted P-value (FDR < 0.10)")

  p3 <- make_vol(res, "nlog10_pi", expression(-log[10](pi~score)),
                 top_pi, -log10(0.05)) +
    ggtitle("Pi-score (Pi < 0.05)") +
    theme(legend.position = "right") + labs(color = NULL)

  # Top 25 hits table
  tbl_data <- res %>%
    arrange(pi_score) %>%
    slice_head(n = 25) %>%
    transmute(
      UniProt = uniprot_id, Gene = gene,
      logFC = sprintf("%.3f", logFC),
      `P.Value` = formatC(P.Value, format = "e", digits = 2),
      `adj.P.Val` = formatC(adj.P.Val, format = "e", digits = 2),
      `Pi-score` = formatC(pi_score, format = "e", digits = 2)
    )

  p_tbl <- tableGrob(tbl_data, rows = NULL,
    theme = ttheme_minimal(base_size = 8,
      core    = list(fg_params = list(hjust = 0, x = 0.02)),
      colhead = list(fg_params = list(hjust = 0, x = 0.02, fontface = "bold"))))

  n_fdr <- sum(res$adj.P.Val < 0.10, na.rm = TRUE)
  n_nom <- sum(res$P.Value < 0.01, na.rm = TRUE)
  n_pi  <- sum(res$sig_pi != 0)

  contrast_out <- file.path(REPORT_DIR, cname)

  pdf(file.path(contrast_out, "summary.pdf"), width = 16, height = 14)
  # Page 1: Triple volcano
  print(
    (p1 | p2 | p3) +
      plot_annotation(
        title = cname,
        subtitle = sprintf("FDR<0.10: %d | P<0.01: %d | Pi<0.05: %d | %d proteins tested",
                           n_fdr, n_nom, n_pi, nrow(res)),
        theme = theme(plot.title = element_text(face = "bold", size = 14)))
  )
  # Page 2: Top hits table
  grid::grid.newpage()
  grid::grid.draw(arrangeGrob(
    p_tbl,
    top = grid::textGrob(
      sprintf("%s \u2014 Top 25 Proteins by Pi-score", cname),
      gp = grid::gpar(fontface = "bold", fontsize = 14)),
    bottom = grid::textGrob(
      "UniProt ID | Gene | log2FC | Nominal P | BH-adjusted P | Pi-score",
      gp = grid::gpar(fontsize = 9, col = "grey40"))
  ))
  dev.off()

  cat(sprintf("  %s: FDR=%d, Nominal=%d, Pi=%d\n", cname, n_fdr, n_nom, n_pi))
}

# --- Overall DEP summary ---

all_counts <- list()
for (cname in contrast_names) {
  fpath <- file.path(per_contrast_dir, paste0(cname, ".csv"))
  if (!file.exists(fpath)) next
  res <- read_csv(fpath, show_col_types = FALSE)
  all_counts <- c(all_counts, list(
    tibble(contrast = cname, criterion = "FDR < 0.10",
           up = sum(res$adj.P.Val < 0.10 & res$logFC > 0, na.rm = TRUE),
           down = sum(res$adj.P.Val < 0.10 & res$logFC < 0, na.rm = TRUE)),
    tibble(contrast = cname, criterion = "P < 0.01",
           up = sum(res$sig_nominal_01 == 1, na.rm = TRUE),
           down = sum(res$sig_nominal_01 == -1, na.rm = TRUE)),
    tibble(contrast = cname, criterion = "Pi < 0.05",
           up = sum(res$sig_pi == 1, na.rm = TRUE),
           down = sum(res$sig_pi == -1, na.rm = TRUE))
  ))
}

sc <- bind_rows(all_counts) %>%
  pivot_longer(c(up, down), names_to = "direction", values_to = "count") %>%
  mutate(
    signed    = if_else(direction == "down", -count, count),
    direction = factor(str_to_title(direction), levels = c("Up", "Down")),
    criterion = factor(criterion, levels = c("FDR < 0.10", "P < 0.01", "Pi < 0.05"))
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

stbl <- bind_rows(all_counts) %>%
  mutate(total = up + down) %>%
  arrange(contrast, criterion)

p_stbl <- tableGrob(stbl, rows = NULL,
  theme = ttheme_minimal(base_size = 10,
    colhead = list(fg_params = list(fontface = "bold"))))

pdf(file.path(REPORT_DIR, "DEP_summary.pdf"), width = 14, height = 12)
print(
  p_bar / wrap_elements(p_stbl) +
    plot_layout(heights = c(2, 1)) +
    plot_annotation(
      title = "YvO Differential Expression Summary",
      theme = theme(plot.title = element_text(face = "bold", size = 16)))
)
dev.off()

# === SUMMARY ==================================================================

print(dal$design$contrast_matrix)

da_summary <- read_csv(file.path(DATA_DIR, "DA_summary.csv"), show_col_types = FALSE)
print(da_summary)

# ============================================================================
# NOTE ON PATHWAY ENRICHMENT (fGSEA)
# ============================================================================
# fGSEA is run independently per contrast without cross-contrast correction.
# This is standard practice in the field (each contrast is treated as an
# independent biological question), but running the same gene sets across
# 5 contrasts increases the chance of finding spurious enrichment in at least
# one contrast. This limitation should be acknowledged in the manuscript
# methods section.
#
# References:
#   - Subramanian et al. (2005) PNAS: Original GSEA methodology
#   - Korotkevich et al. (2021) bioRxiv: fGSEA implementation
# ============================================================================

# ============================================================================
# PHASE 1: ROBUSTNESS ENHANCEMENTS
# ============================================================================

# === 1a. KS + LEVENE'S DISTRIBUTIONAL TEST (H1 — Blunting Diagnostics) ======
# Tests whether Training_Young and Training_Old |logFC| distributions differ,
# providing evidence for age-dependent response attenuation.
# Reference: Conover 1999, Practical Nonparametric Statistics
# Reference: Levene 1960; Brown-Forsythe variant

cat("\n--- Blunting Diagnostics ---\n")

comb <- read_csv(combined_path, show_col_types = FALSE)

blunt_df <- comb %>%
  filter(!is.na(logFC_Training_Young) & !is.na(logFC_Training_Old)) %>%
  transmute(
    abs_lfc_young = abs(logFC_Training_Young),
    abs_lfc_old   = abs(logFC_Training_Old)
  )

# KS test: do the two |logFC| distributions differ in shape?
ks_res <- ks.test(blunt_df$abs_lfc_young, blunt_df$abs_lfc_old)

# Fligner-Killeen test for homogeneity of variances (robust to non-normality)
fk_data <- data.frame(
  abs_lfc = c(blunt_df$abs_lfc_young, blunt_df$abs_lfc_old),
  contrast = factor(rep(c("Training_Young", "Training_Old"),
                        each = nrow(blunt_df)))
)
fk_res <- fligner.test(abs_lfc ~ contrast, data = fk_data)

blunt_diag <- tibble(
  test = c("Kolmogorov-Smirnov", "Fligner-Killeen"),
  statistic = c(ks_res$statistic, fk_res$statistic),
  p_value = c(ks_res$p.value, fk_res$p.value),
  interpretation = c(
    ifelse(ks_res$p.value < 0.05,
           "Distributions differ significantly",
           "No significant distributional difference"),
    ifelse(fk_res$p.value < 0.05,
           "Variances differ significantly (blunting evidence)",
           "No significant variance difference")
  )
)
write_csv(blunt_diag, file.path(DATA_DIR, "blunting_diagnostics.csv"))

cat(sprintf("  KS test: D = %.4f, p = %.4g\n", ks_res$statistic, ks_res$p.value))
cat(sprintf("  Fligner-Killeen: chi-sq = %.4f, p = %.4g\n",
            fk_res$statistic, fk_res$p.value))

# Density overlay plot
density_long <- bind_rows(
  tibble(abs_lfc = blunt_df$abs_lfc_young, contrast = "Training (Young)"),
  tibble(abs_lfc = blunt_df$abs_lfc_old,   contrast = "Training (Old)")
)
med_young <- median(blunt_df$abs_lfc_young)
med_old   <- median(blunt_df$abs_lfc_old)

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
           label = sprintf("KS D = %.3f, p = %.3g\nFligner chi2 = %.2f, p = %.3g",
                           ks_res$statistic, ks_res$p.value,
                           fk_res$statistic, fk_res$p.value),
           size = 2.5, fontface = "bold", color = "grey25") +
  scale_fill_manual(values = c("Training (Young)" = "#E05A4E",
                                "Training (Old)" = "#5DA5DA")) +
  scale_color_manual(values = c("Training (Young)" = "#E05A4E",
                                 "Training (Old)" = "#5DA5DA")) +
  labs(title = "Effect Size Distribution: Training Response by Age",
       subtitle = "Density of |logFC| — attenuation evidence",
       x = "|logFC|", y = "Density", fill = NULL, color = NULL) +
  theme_bw(base_size = 11) +
  theme(plot.title = element_text(face = "bold", size = 12),
        legend.position = "bottom")

pdf(file.path(REPORT_DIR, "blunting_density.pdf"), width = 7, height = 5.5)
print(p_density)
dev.off()
cat("  Saved: blunting_density.pdf, blunting_diagnostics.csv\n")

# === 1b. BOOTSTRAP CI FOREST PLOT (H1, supplementary) =======================
# Median |logFC| with 95% BCa bootstrap CI per contrast
# Reference: Efron & Tibshirani 1993, An Introduction to the Bootstrap

cat("\n--- Bootstrap CI Forest Plot ---\n")

suppressPackageStartupMessages(library(boot))

median_fn <- function(d, i) median(d[i], na.rm = TRUE)

boot_contrasts <- c("Training_Young", "Training_Old", "Aging", "Interaction")
boot_results <- list()

for (cname in boot_contrasts) {
  lfc_col <- paste0("logFC_", cname)
  if (!(lfc_col %in% names(comb))) next
  vals <- abs(comb[[lfc_col]])
  vals <- vals[!is.na(vals)]

  b <- boot(vals, median_fn, R = 10000)
  ci <- tryCatch(boot.ci(b, type = "bca"),
                 error = function(e) boot.ci(b, type = "perc"))
  ci_lo <- if (!is.null(ci$bca)) ci$bca[4] else ci$percent[4]
  ci_hi <- if (!is.null(ci$bca)) ci$bca[5] else ci$percent[5]

  boot_results <- c(boot_results, list(tibble(
    contrast      = cname,
    median_absLFC = median(vals),
    ci_lower      = ci_lo,
    ci_upper      = ci_hi,
    boot_se       = sd(b$t),
    n_proteins    = length(vals)
  )))
}

boot_df <- bind_rows(boot_results)
write_csv(boot_df, file.path(DATA_DIR, "effect_size_bootstrap.csv"))

pal_contrast <- c(Training_Young = "#E05A4E", Training_Old = "#5DA5DA",
                  Aging = "#4CAF50", Interaction = "#9B7FBF")

p_forest <- ggplot(boot_df, aes(y = reorder(contrast, median_absLFC),
                                 x = median_absLFC, color = contrast)) +
  geom_pointrange(aes(xmin = ci_lower, xmax = ci_upper),
                  size = 0.8, linewidth = 0.8) +
  scale_color_manual(values = pal_contrast, guide = "none") +
  labs(title = "Effect Size by Contrast",
       subtitle = "Median |logFC| with 95% BCa bootstrap CI (10,000 resamples)",
       x = "Median |logFC|", y = NULL) +
  theme_bw(base_size = 11) +
  theme(plot.title = element_text(face = "bold", size = 12))

pdf(file.path(REPORT_DIR, "effect_size_forest.pdf"), width = 7, height = 4)
print(p_forest)
dev.off()
cat("  Saved: effect_size_forest.pdf, effect_size_bootstrap.csv\n")
print(boot_df)

# === 1c. POWER ANALYSIS (H1) ================================================
# Minimum detectable logFC at 80% power for each contrast
# Reference: Cohen 1988, Statistical Power Analysis

cat("\n--- Power Analysis ---\n")

suppressPackageStartupMessages(library(pwr))

# Get residual SD from limma fit
fit <- dal$eBayes_fit
sigma_residual <- sqrt(mean(fit$sigma^2, na.rm = TRUE))

# Within-subject correlation (from duplicateCorrelation)
within_cor <- if (!is.null(fit$correlation)) fit$correlation else
  if (!is.null(dal$tags$duplicate_correlation)) dal$tags$duplicate_correlation else NA

# Sample sizes per group
n_young <- sum(meta$age == "Young" & meta$time == "Pre")
n_old   <- sum(meta$age == "Old"   & meta$time == "Pre")

power_results <- list()
power_contrasts <- c("Training_Young", "Training_Old", "Aging", "Interaction")
for (cname in power_contrasts) {
  n_subj <- switch(cname,
    Training_Young = n_young,
    Training_Old   = n_old,
    Aging          = min(n_young, n_old),
    Interaction    = min(n_young, n_old)
  )

  # For paired contrasts, effective SD is reduced by within-subject correlation
  paired <- cname %in% c("Training_Young", "Training_Old")
  if (paired && !is.na(within_cor)) {
    effective_sigma <- sigma_residual * sqrt(2 * (1 - within_cor))
  } else {
    effective_sigma <- sigma_residual * sqrt(2)
  }

  # Cohen's d = logFC / effective_sigma; solve for d at 80% power
  pw <- pwr.t.test(n = n_subj, d = NULL, sig.level = 0.10,
                   power = 0.80,
                   type = if (paired) "paired" else "two.sample")
  min_lfc <- pw$d * effective_sigma

  power_results <- c(power_results, list(tibble(
    contrast          = cname,
    n_subjects        = n_subj,
    within_cor        = ifelse(paired, within_cor, NA_real_),
    effective_sigma   = round(effective_sigma, 4),
    min_detectable_d  = round(pw$d, 4),
    min_detectable_logFC = round(min_lfc, 4),
    power             = 0.80,
    alpha             = 0.10
  )))
}

power_df <- bind_rows(power_results)
write_csv(power_df, file.path(DATA_DIR, "power_analysis.csv"))
cat("  Power analysis results:\n")
print(as.data.frame(power_df))

# === 1d. IMPUTATION SENSITIVITY ANALYSIS ====================================
# Compare t-statistics from non-imputed limma (main) vs imputed limma
# Reference: Karpievitch et al. 2019, BMC Bioinform 20:391

cat("\n--- Imputation Sensitivity Analysis ---\n")

imp_path <- "02_Imputation/c_data/01_imputed.csv"
if (file.exists(imp_path)) {
  imp_data <- read_csv(imp_path, show_col_types = FALSE)
  imp_ann  <- imp_data[, ann_cols]
  imp_samp <- setdiff(names(imp_data), ann_cols)
  imp_mat  <- as.matrix(imp_data[, imp_samp])
  rownames(imp_mat) <- imp_ann$uniprot_id

  # Only use samples present in both datasets
  shared_samps <- intersect(colnames(mat), colnames(imp_mat))
  imp_mat_sub  <- imp_mat[, shared_samps]

  # Build same design/model on imputed data
  meta_imp <- meta %>% filter(sample_id %in% shared_samps)
  meta_imp_df <- as.data.frame(meta_imp)
  rownames(meta_imp_df) <- meta_imp$sample_id

  dal_imp <- DAList(
    data       = imp_mat_sub,
    annotation = as.data.frame(imp_ann),
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

  # Compare t-statistics
  sens_results <- list()
  imp_combined_dir <- file.path(tempdir(), "imp_results")
  dir.create(imp_combined_dir, showWarnings = FALSE)
  write_limma_tables(dal_imp, output_dir = imp_combined_dir, overwrite = TRUE)
  imp_combined <- read_csv(file.path(imp_combined_dir, "combined_results.csv"),
                           show_col_types = FALSE)

  sens_contrasts <- c("Training_Young", "Training_Old", "Aging", "Interaction")
  scatter_list <- list()

  for (cname in sens_contrasts) {
    t_col <- paste0("t_", cname)
    adj_col <- paste0("adj.P.Val_", cname)

    if (!(t_col %in% names(comb)) || !(t_col %in% names(imp_combined))) next

    merged <- inner_join(
      comb %>% dplyr::select(uniprot_id, t_nonimp = all_of(t_col),
                              padj_nonimp = all_of(adj_col)),
      imp_combined %>% dplyr::select(uniprot_id, t_imp = all_of(t_col),
                                      padj_imp = all_of(adj_col)),
      by = "uniprot_id"
    ) %>% filter(!is.na(t_nonimp) & !is.na(t_imp))

    sp <- cor.test(merged$t_nonimp, merged$t_imp, method = "spearman")

    sens_results <- c(sens_results, list(tibble(
      contrast    = cname,
      spearman_rho = round(sp$estimate, 4),
      p_value     = sp$p.value,
      n_proteins  = nrow(merged)
    )))

    # Flag proteins that switch significance
    merged <- merged %>%
      mutate(
        sig_nonimp = padj_nonimp < 0.10,
        sig_imp    = padj_imp < 0.10,
        switch     = case_when(
          sig_nonimp & !sig_imp ~ "Lost in imputed",
          !sig_nonimp & sig_imp ~ "Gained in imputed",
          TRUE ~ "Concordant"
        )
      )

    scatter_list[[cname]] <- ggplot(merged, aes(x = t_nonimp, y = t_imp,
                                                 color = switch)) +
      geom_abline(slope = 1, intercept = 0, linetype = "dashed",
                  color = "grey50", linewidth = 0.3) +
      geom_point(alpha = 0.4, size = 0.8) +
      scale_color_manual(values = c("Concordant" = "grey60",
                                     "Lost in imputed" = "#D6604D",
                                     "Gained in imputed" = "#4393C3"),
                         name = NULL) +
      annotate("text", x = -Inf, y = Inf, hjust = -0.05, vjust = 1.3,
               label = sprintf("rho = %.3f", sp$estimate),
               size = 3, fontface = "bold", color = "grey25") +
      labs(title = cname, x = "t (non-imputed)", y = "t (imputed)") +
      theme_bw(base_size = 9) +
      theme(plot.title = element_text(face = "bold", size = 10),
            legend.position = "bottom",
            legend.text = element_text(size = 7))
  }

  sens_df <- bind_rows(sens_results)
  write_csv(sens_df, file.path(DATA_DIR, "imputation_sensitivity.csv"))

  if (length(scatter_list) > 0) {
    p_sens <- wrap_plots(scatter_list, ncol = 2) +
      plot_annotation(
        title    = "Imputation Sensitivity: t-statistic Concordance",
        subtitle = "Non-imputed (limma, main) vs imputed limma",
        theme = theme(plot.title = element_text(face = "bold", size = 14),
                      plot.subtitle = element_text(size = 10, color = "grey30")))

    pdf(file.path(REPORT_DIR, "imputation_sensitivity.pdf"), width = 11, height = 9.5)
    print(p_sens)
    dev.off()
  }

  cat("  Sensitivity results:\n")
  print(as.data.frame(sens_df))
  cat("  Saved: imputation_sensitivity.pdf, imputation_sensitivity.csv\n")
} else {
  cat("  Imputed data not found — skipping sensitivity analysis\n")
}

cat("\n=== YvO limma DEP complete ===\n")
