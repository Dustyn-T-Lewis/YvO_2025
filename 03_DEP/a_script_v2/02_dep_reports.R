# --- YvO DEP — Visualization & Reporting ------------------------------------
# Reads fitted DAList and per-contrast results from c_data_v2/
# Produces per-contrast volcano + top-25 table PDFs and overview PDF
# Depends on: 01_run_dep.R outputs
# ---------------------------------------------------------------------------

# --- SETUP ---

library(dplyr)
library(tidyr)
library(tibble)
library(stringr)
library(readr)
library(purrr)
library(ggplot2)
library(ggrepel)
library(patchwork)
library(gridExtra)
library(limma)

setwd(rprojroot::find_rstudio_root_file())

cfg <- list(
  data_dir     = "03_DEP/c_data_v2",
  per_dir      = "03_DEP/c_data_v2/04_per_contrast_results",
  report_dir   = "03_DEP/b_reports_v2",
  summary_dir  = "03_DEP/b_reports_v2/03_contrast_summaries"
)

dir.create(cfg$summary_dir, recursive = TRUE, showWarnings = FALSE)

# --- THEME & PALETTE ---

theme_dep <- theme_bw(base_size = 11) +
  theme(plot.title = element_text(face = "bold", size = 12),
        legend.position = "bottom")
pal_dir <- c(Up = "#D6604D", Down = "#4393C3", NS = "grey70")

# --- LOAD DATA ---

dal <- readRDS(file.path(cfg$data_dir, "01_limma_DAList.rds"))
da_summary <- read_csv(file.path(cfg$data_dir, "02_DA_summary.csv"),
                        show_col_types = FALSE)

contrast_names <- names(dal$results)

# Read per-contrast CSVs into a list
results_list <- lapply(contrast_names, function(cname) {
  read_csv(file.path(cfg$per_dir, paste0(cname, ".csv")), show_col_types = FALSE)
})
names(results_list) <- contrast_names

# --- PER-CONTRAST VOLCANO + TOP-25 TABLE ---

for (cname in contrast_names) {
  res <- results_list[[cname]] |>
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

  tbl_data <- res |>
    slice_min(pi_score, n = 25, with_ties = FALSE) |>
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
  contrast_out <- file.path(cfg$summary_dir, cname)
  dir.create(contrast_out, showWarnings = FALSE)

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

# --- OVERALL DEP SUMMARY BAR CHART ---

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
}) |>
  pivot_longer(c(up, down), names_to = "direction", values_to = "count") |>
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

# --- WRITE OVERVIEW PDF ---

pdf(file.path(cfg$report_dir, "02_dep_overview.pdf"), width = 14, height = 10)
print(
  p_bar / wrap_elements(p_stbl) +
    plot_layout(heights = c(3, 1.2)) +
    plot_annotation(
      title = "YvO Differential Expression Overview",
      theme = theme(plot.title = element_text(face = "bold", size = 16)))
)

# Page 2: Outlier-Removal Sensitivity

run_limma_sens <- function(mat, meta) {
  meta$Group_Time <- factor(meta$Group_Time,
    levels = c("Young_Pre", "Young_Post", "Old_Pre", "Old_Post"))
  design <- model.matrix(~ 0 + Group_Time, data = meta)
  colnames(design) <- gsub("^Group_Time", "", colnames(design))
  corfit <- duplicateCorrelation(mat, design, block = meta$Subject_ID)
  aw     <- arrayWeights(mat, design)
  fit    <- lmFit(mat, design, block = meta$Subject_ID,
                  correlation = corfit$consensus.correlation, weights = aw)
  cm <- makeContrasts(
    Training_Young = Young_Post - Young_Pre,
    Training_Old   = Old_Post - Old_Pre,
    Aging          = Old_Pre - Young_Pre,
    Interaction    = (Old_Post - Old_Pre) - (Young_Post - Young_Pre),
    Reversal       = (Old_Post - Old_Pre) - (Old_Pre - Young_Pre),
    levels = design)
  eBayes(contrasts.fit(fit, cm), robust = TRUE, trend = TRUE)
}

# Full dataset (all 64 samples, no outlier removal)
int <- readRDS("01_normalization/c_data_v2/00_report_intermediates.rds")
full_mat  <- normalizeBetweenArrays(log2(int$data_pre_outlier + 1),
                                     method = "cyclicloess")
full_fit  <- run_limma_sens(full_mat, int$meta_pre_outlier)

# Reduced dataset (62 samples, outliers removed — v2 normalized)
norm_v2   <- read_csv("01_normalization/c_data_v2/02_normalized.csv",
                       show_col_types = FALSE)
ann_cols  <- c("uniprot_id", "protein", "gene", "description")
red_mat   <- as.matrix(norm_v2[, setdiff(names(norm_v2), ann_cols)])
rownames(red_mat) <- norm_v2$uniprot_id
red_meta  <- readRDS("01_normalization/c_data_v2/03_DAList_normalized.rds")$metadata
red_fit   <- run_limma_sens(red_mat, red_meta)

pi_thresh <- 0.05
sens_compare <- do.call(rbind, lapply(contrast_names, function(cname) {
  full_res <- topTable(full_fit, coef = cname, number = Inf, sort.by = "none")
  red_res  <- topTable(red_fit, coef = cname, number = Inf, sort.by = "none")

  full_pi <- full_res$P.Value ^ abs(full_res$logFC)
  red_pi  <- red_res$P.Value ^ abs(red_res$logFC)

  shared  <- intersect(rownames(full_res), rownames(red_res))
  tibble(
    Contrast      = cname,
    FDR_full      = sum(full_res$adj.P.Val < 0.10, na.rm = TRUE),
    FDR_reduced   = sum(red_res$adj.P.Val < 0.10, na.rm = TRUE),
    Pi_full       = sum(full_pi < pi_thresh, na.rm = TRUE),
    Pi_reduced    = sum(red_pi < pi_thresh, na.rm = TRUE),
    Pearson_r     = cor(full_res[shared, "logFC"], red_res[shared, "logFC"],
                        use = "complete.obs"),
    Spearman_rho  = cor(full_res[shared, "logFC"], red_res[shared, "logFC"],
                        use = "complete.obs", method = "spearman"))
}))

write.csv(sens_compare, file.path(cfg$data_dir, "11_outlier_sensitivity.csv"),
          row.names = FALSE)

sens_display <- sens_compare |>
  mutate(across(c(Pearson_r, Spearman_rho), ~ sprintf("%.3f", .x)))

p_sens_tbl <- tableGrob(sens_display, rows = NULL,
  theme = ttheme_minimal(base_size = 10,
    core    = list(fg_params = list(hjust = 0.5)),
    colhead = list(fg_params = list(fontface = "bold", hjust = 0.5))))

outlier_str <- paste(int$outlier_ids, collapse = ", ")
print(
  wrap_elements(p_sens_tbl) +
    plot_annotation(
      title = "Outlier-Removal Sensitivity (64 vs 62 samples)",
      subtitle = sprintf("Removed: %s (4-method consensus, >=3/4)", outlier_str),
      theme = theme(
        plot.title = element_text(face = "bold", size = 14),
        plot.subtitle = element_text(size = 11))))
dev.off()

cat("Done: 02_dep_reports.R\n")
