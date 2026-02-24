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

setwd("/Users/dtl0018/Desktop/A_Proteomics_Analysis/A_YvO_2025")
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
  subject   = dal_meta$Subject_ID
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

write_limma_plots(dal,
                  grouping_column = "group",
                  output_dir      = REPORT_DIR,
                  table_columns   = c("uniprot_id", "gene", "protein"),
                  title_column    = "gene",
                  overwrite       = TRUE)

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

cat("\n=== YvO limma DEP complete ===\n")
