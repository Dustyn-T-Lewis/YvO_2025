#!/usr/bin/env Rscript
# Stage 03: DEP reports
# Per-contrast volcano + top-25 table PDFs, overview bar chart, outlier sensitivity

library(dplyr)
library(tidyr)
library(tibble)
library(purrr)
library(readxl)
library(ggplot2)
library(ggrepel)
library(patchwork)
library(gridExtra)
library(limma)

DAT <- here::here("03_DEP", "c_data")
RPT <- here::here("03_DEP", "b_reports")
SUM <- file.path(RPT, "03_contrast_summaries")
dir.create(SUM, recursive = TRUE, showWarnings = FALSE)

XLSX <- file.path(DAT, "03_DEP_results.xlsx")

theme_dep <- theme_bw(base_size = 11) +
  theme(plot.title = element_text(face = "bold", size = 12),
        legend.position = "bottom")
pal_dir <- c(Up = "#D6604D", Down = "#4393C3", NS = "grey70")

# ── Load ──────────────────────────────────────────────────────────────────────

dal <- readRDS(file.path(DAT, "01_limma_DAList.rds"))
contrast_names <- names(dal$results)

results_list <- lapply(contrast_names, \(cn) {
  as.data.frame(read_excel(XLSX, sheet = cn))
})
names(results_list) <- contrast_names

da_summary <- as.data.frame(read_excel(XLSX, sheet = "DA_summary"))

# ── Per-contrast volcano + top-25 table ──────────────────────────────────────

for (cname in contrast_names) {
  res <- results_list[[cname]] |>
    mutate(
      nlog10_pval = -log10(pmax(P.Value, 1e-300)),
      nlog10_adj  = -log10(pmax(adj.P.Val, 1e-300)),
      nlog10_pi   = -log10(pmax(pi_score, 1e-300)),
      dir_pi = factor(case_when(
        sig_pi ==  1 ~ "Up", sig_pi == -1 ~ "Down", TRUE ~ "NS"),
        levels = c("Up", "Down", "NS")))

  make_vol <- function(df, ycol, ylab, top_n, thresh) {
    top <- slice_min(df, .data[[ycol]] * -1, n = top_n, with_ties = FALSE)
    ggplot(df, aes(logFC, .data[[ycol]], color = dir_pi)) +
      geom_point(alpha = 0.35, size = 1) +
      geom_hline(yintercept = thresh, linetype = "dashed", color = "grey40",
                 linewidth = 0.4) +
      geom_text_repel(data = top, aes(label = gene), size = 2.5,
                      max.overlaps = 15, show.legend = FALSE, seed = 42) +
      scale_color_manual(values = pal_dir, drop = FALSE) +
      labs(x = expression(log[2]~FC), y = ylab) +
      theme_dep + theme(legend.position = "none")
  }

  p1 <- make_vol(res, "nlog10_pval", expression(-log[10](P)), 10,
                 -log10(0.01)) + ggtitle("Nominal P < 0.01")
  p2 <- make_vol(res, "nlog10_adj", expression(-log[10](FDR)), 10,
                 -log10(0.10)) + ggtitle("FDR < 0.10")
  p3 <- make_vol(res, "nlog10_pi", expression(-log[10](Pi)), 10,
                 -log10(0.05)) +
    ggtitle("Pi < 0.05") + theme(legend.position = "right") + labs(color = NULL)

  tbl_data <- res |>
    slice_min(pi_score, n = 25, with_ties = FALSE) |>
    transmute(Gene = gene, logFC = sprintf("%.3f", logFC),
              P = formatC(P.Value, format = "e", digits = 2),
              FDR = formatC(adj.P.Val, format = "e", digits = 2),
              Pi = formatC(pi_score, format = "e", digits = 2))

  p_tbl <- tableGrob(tbl_data, rows = NULL,
    theme = ttheme_minimal(base_size = 8,
      core    = list(fg_params = list(hjust = 0, x = 0.02)),
      colhead = list(fg_params = list(hjust = 0, x = 0.02, fontface = "bold"))))

  n_fdr <- sum(res$adj.P.Val < 0.10, na.rm = TRUE)
  n_pi  <- sum(res$sig_pi != 0)
  cdir  <- file.path(SUM, cname)
  dir.create(cdir, showWarnings = FALSE)

  pdf(file.path(cdir, "summary.pdf"), width = 16, height = 14)
  print((p1 | p2 | p3) + plot_annotation(
    title = cname,
    subtitle = sprintf("FDR<0.10: %d | Pi<0.05: %d | %d proteins",
                       n_fdr, n_pi, nrow(res)),
    theme = theme(plot.title = element_text(face = "bold", size = 14))))
  grid::grid.newpage()
  grid::grid.draw(arrangeGrob(p_tbl,
    top = grid::textGrob(sprintf("%s \u2014 Top 25 by Pi", cname),
                         gp = grid::gpar(fontface = "bold", fontsize = 14))))
  dev.off()
  message(sprintf("  %s: FDR<0.10=%d, Pi<0.05=%d", cname, n_fdr, n_pi))
}

# ── Overview bar chart ────────────────────────────────────────────────────────

sc <- list_rbind(lapply(contrast_names, \(cname) {
  res <- results_list[[cname]]
  bind_rows(
    tibble(contrast = cname, criterion = "FDR < 0.10",
           up = sum(res$adj.P.Val < 0.10 & res$logFC > 0, na.rm = TRUE),
           down = sum(res$adj.P.Val < 0.10 & res$logFC < 0, na.rm = TRUE)),
    tibble(contrast = cname, criterion = "Pi < 0.05",
           up = sum(res$sig_pi == 1), down = sum(res$sig_pi == -1)))
})) |>
  pivot_longer(c(up, down), names_to = "direction", values_to = "count") |>
  mutate(signed = if_else(direction == "down", -count, count),
         direction = factor(stringr::str_to_title(direction),
                           levels = c("Up", "Down")),
         criterion = factor(criterion))

p_bar <- ggplot(sc, aes(contrast, signed, fill = direction)) +
  geom_col(position = "identity", width = 0.7) +
  geom_hline(yintercept = 0, linewidth = 0.3) +
  geom_text(data = filter(sc, signed > 0),
            aes(label = count), vjust = -0.3, size = 3) +
  geom_text(data = filter(sc, signed < 0),
            aes(label = count), vjust = 1.3, size = 3) +
  facet_wrap(~criterion) +
  scale_fill_manual(values = c(Up = "#D6604D", Down = "#4393C3")) +
  labs(title = "DEP Counts", x = NULL, y = "DEPs (Up / Down)", fill = NULL) +
  theme_dep + theme(axis.text.x = element_text(angle = 30, hjust = 1))

# ── Outlier sensitivity (64 vs 62 samples) ───────────────────────────────────

run_limma_sens <- function(mat, meta) {
  meta$Group_Time <- factor(meta$Group_Time,
    levels = c("Young_Pre", "Young_Post", "Old_Pre", "Old_Post"))
  design <- model.matrix(~ 0 + Group_Time, data = meta)
  colnames(design) <- gsub("^Group_Time", "", colnames(design))
  block_id <- sub("_(Pre|Post)$", "", meta$Col_ID)  # group-prefixed: Young_S01 ≠ Old_S01
  corfit <- duplicateCorrelation(mat, design, block = block_id)
  fit <- lmFit(mat, design, block = block_id,
               correlation = corfit$consensus.correlation)
  cm <- makeContrasts(
    Training_Young = Young_Post - Young_Pre,
    Training_Old   = Old_Post - Old_Pre,
    Aging          = Old_Pre - Young_Pre,
    Interaction    = (Old_Post - Old_Pre) - (Young_Post - Young_Pre),
    Reversal       = (Old_Post - Old_Pre) - (Old_Pre - Young_Pre),
    levels = design)
  eBayes(contrasts.fit(fit, cm), robust = TRUE)
}

int <- readRDS(here::here("01_normalization", "c_data",
                          "00_report_intermediates.rds"))
full_mat <- normalizeBetweenArrays(log2(int$data_pre_outlier + 1),
                                    method = "cyclicloess")
full_fit <- run_limma_sens(full_mat, int$meta_pre_outlier)

dal_norm <- readRDS(here::here("01_normalization", "c_data",
                               "03_DAList_normalized.rds"))
red_mat  <- as.matrix(dal_norm$data)
red_fit  <- run_limma_sens(red_mat, dal_norm$metadata)

sens_compare <- list_rbind(lapply(contrast_names, \(cname) {
  full_res <- topTable(full_fit, coef = cname, number = Inf, sort.by = "none")
  red_res  <- topTable(red_fit,  coef = cname, number = Inf, sort.by = "none")
  shared   <- intersect(rownames(full_res), rownames(red_res))
  tibble(Contrast = cname,
         FDR_full    = sum(full_res$adj.P.Val < 0.10, na.rm = TRUE),
         FDR_reduced = sum(red_res$adj.P.Val < 0.10, na.rm = TRUE),
         Pi_full     = sum(full_res$P.Value^abs(full_res$logFC) < 0.05, na.rm = TRUE),
         Pi_reduced  = sum(red_res$P.Value^abs(red_res$logFC) < 0.05, na.rm = TRUE),
         Pearson_r   = cor(full_res[shared, "logFC"], red_res[shared, "logFC"],
                          use = "complete.obs"),
         Spearman_rho = cor(full_res[shared, "logFC"], red_res[shared, "logFC"],
                           use = "complete.obs", method = "spearman"))
}))

# Add sensitivity sheet to xlsx
wb <- loadWorkbook(XLSX)
write_sheet_simple <- function(wb, name, data) {
  addWorksheet(wb, name)
  writeData(wb, name, data,
    headerStyle = createStyle(textDecoration = "bold", fgFill = "#DCE6F1"))
  freezePane(wb, name, firstRow = TRUE)
  setColWidths(wb, name, cols = seq_len(ncol(data)), widths = "auto")
}
write_sheet_simple(wb, "outlier_sensitivity", sens_compare)
saveWorkbook(wb, XLSX, overwrite = TRUE)

# ── Assemble overview PDF ─────────────────────────────────────────────────────

sens_display <- sens_compare |>
  mutate(across(c(Pearson_r, Spearman_rho), \(x) sprintf("%.3f", x)))
p_sens <- tableGrob(sens_display, rows = NULL,
  theme = ttheme_minimal(base_size = 10,
    core    = list(fg_params = list(hjust = 0.5)),
    colhead = list(fg_params = list(fontface = "bold", hjust = 0.5))))

pdf(file.path(RPT, "02_dep_overview.pdf"), width = 14, height = 10)
print(p_bar + plot_annotation(
  title = "YvO Differential Expression Overview",
  theme = theme(plot.title = element_text(face = "bold", size = 16))))
print(wrap_elements(p_sens) + plot_annotation(
  title = "Outlier-Removal Sensitivity (64 vs 62 samples)",
  subtitle = sprintf("Removed: %s", paste(int$outlier_ids, collapse = ", ")),
  theme = theme(plot.title = element_text(face = "bold", size = 14),
                plot.subtitle = element_text(size = 11))))
dev.off()
message("Saved: ", file.path(RPT, "02_dep_overview.pdf"))
