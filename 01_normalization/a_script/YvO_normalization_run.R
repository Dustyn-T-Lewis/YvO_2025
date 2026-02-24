#!/usr/bin/env Rscript
###############################################################################
#   YvO Normalization — DIA-MS proteomics (Young vs Old, skeletal muscle)
#   HPA tissue filter, 66% group missingness, consensus outlier detection,
#   cycloess normalization via proteoDA
#
#   Reports:
#     01_normalization_report.pdf  — proteoDA norm method comparison
#     02_qc_report_pre.pdf         — proteoDA QC before normalization
#     03_qc_report_post.pdf        — proteoDA QC after normalization
#     04_filtering_report.pdf      — HPA + missingness filtering diagnostics
#     05_diagnostics_report.pdf    — variability, outliers, post-norm PCA
#
#   Methodological references:
#     - Thurman et al. 2023, JOSS 8(85):5184 — proteoDA package
#     - Bolstad et al. 2003, Bioinformatics 19(2):185-193 — cyclic loess
#     - Valikangas et al. 2018, Brief Bioinform 19(1):1-11
#       Systematic evaluation of normalization methods for label-free proteomics
#     - PRONE: Arend et al. 2025, Brief Bioinform 26(3):bbaf201
#       Normalization evaluation showing dataset-specific method selection
#     - Brenes 2024, J Proteome Res, PMC11629372
#       CV calculation on linear (not log) scale for DIA proteomics
#     - SEAOP: Wen et al. 2024, Brief Bioinform 25(3):bbae129
#       Statistical ensemble outlier detection in proteomics
#
#   Exercise proteomics context:
#     - Hostrup et al. 2022, eLife 11:e69802
#       HIT remodels skeletal muscle proteome & acetylome; 3168 proteins
#     - Deshmukh et al. 2021, Nat Commun 12:304
#       Deep muscle-proteomic analysis of freeze-dried biopsies
#     - Hulmi et al. 2025, J Physiol 603:438-453
#       Skeletal muscle proteomic memory after repeated resistance training;
#       dia-PASEF proteomics with repeated measures
#     - O'Leary et al. 2025, Exp Physiol 110:438-453
#       Skeletal muscle proteomics in young/older women after 8wk RET
###############################################################################

if (!requireNamespace("pacman", quietly = TRUE)) install.packages("pacman")
pacman::p_load(
  proteoDA,
  readxl, readr, dplyr, tidyr, stringr,
  org.Hs.eg.db, GO.db, AnnotationDbi,
  ggplot2, ggrepel, ggnewscale, patchwork
)

set.seed(42)

setwd("/Users/dtl0018/Desktop/A_Proteomics_Analysis/A_YvO_2025")
input_dir  <- "00_input"
report_dir <- "01_normalization/b_reports"
data_dir   <- "01_normalization/c_data"

dir.create(report_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(data_dir,   recursive = TRUE, showWarnings = FALSE)

source("04_Figures/shared/palettes.R")
pal_group_time <- GROUP_FILL
shape_tp       <- SHAPE_TP

# --- Helpers -----------------------------------------------------------------

impute_median <- function(mat) {
  for (j in seq_len(ncol(mat)))
    mat[is.na(mat[, j]), j] <- median(mat[, j], na.rm = TRUE)
  mat
}

run_pca <- function(mat, metadata, log_transform = TRUE) {
  if (log_transform) mat <- log2(impute_median(mat) + 1)
  else mat <- impute_median(mat)
  pca <- prcomp(t(mat), center = TRUE, scale. = TRUE)
  var_exp <- round(summary(pca)$importance[2, 1:3] * 100, 1)
  pc_df <- as.data.frame(pca$x[, 1:3]) %>%
    mutate(Col_ID = rownames(.)) %>%
    left_join(metadata, by = "Col_ID")
  list(pca = pca, pc_df = pc_df, var_exp = var_exp)
}

make_gocc_donut <- function(gene_vec, title_text, n_top = 8) {
  if (length(gene_vec) == 0) return(NULL)

  go_map <- AnnotationDbi::select(
    org.Hs.eg.db, keys = gene_vec,
    keytype = "SYMBOL", columns = c("SYMBOL", "GO", "ONTOLOGY")
  ) %>%
    filter(ONTOLOGY == "CC") %>%
    dplyr::select(gene = SYMBOL, GO) %>%
    left_join(
      AnnotationDbi::select(GO.db, keys = unique(.$GO), columns = "TERM") %>%
        dplyr::select(GO = GOID, term = TERM),
      by = "GO")

  term_freq <- go_map %>% count(term, sort = TRUE)
  top_terms <- term_freq %>% slice_head(n = n_top) %>% pull(term)

  go_primary <- go_map %>%
    left_join(term_freq, by = "term") %>%
    group_by(gene) %>% slice_max(n, n = 1, with_ties = FALSE) %>% ungroup()

  donut_data <- go_primary %>%
    mutate(category = if_else(term %in% top_terms, term, "Other")) %>%
    bind_rows(tibble(gene = setdiff(gene_vec, go_primary$gene),
                     category = "No GO:CC")) %>%
    count(category) %>% arrange(desc(n)) %>%
    mutate(frac = n / sum(n),
           label = sprintf("%s\n(%d)", str_wrap(category, 20), n),
           ymax = cumsum(frac), ymin = lag(ymax, default = 0),
           ymid = (ymin + ymax) / 2)

  ggplot(donut_data, aes(ymax = ymax, ymin = ymin, xmax = 4, xmin = 2.5,
                         fill = category)) +
    geom_rect(color = "white", linewidth = 0.5) +
    geom_text(aes(x = 3.25, y = ymid, label = label), size = 2.5, lineheight = 0.9) +
    annotate("text", x = 0, y = 0,
             label = paste0(length(gene_vec), "\nproteins"),
             size = 5, fontface = "bold") +
    coord_polar(theta = "y") + xlim(c(0, 4.5)) +
    scale_fill_brewer(palette = "Set2", guide = "none") +
    theme_void() + labs(title = title_text)
}

plot_pca_biplot <- function(pca_out, color_pal, annotation_df = NULL,
                            n_loadings = 10, n_top_cc = 6,
                            title_text = "PCA Biplot") {
  pc_df   <- pca_out$pc_df
  pca     <- pca_out$pca
  var_exp <- pca_out$var_exp

  loadings <- as.data.frame(pca$rotation[, 1:2])
  loadings$magnitude <- sqrt(loadings$PC1^2 + loadings$PC2^2)
  loadings$id <- rownames(loadings)
  top_load <- loadings %>% slice_max(magnitude, n = n_loadings)

  if (!is.null(annotation_df)) {
    id_to_gene <- setNames(annotation_df$gene, annotation_df$uniprot_id)
    top_load$label <- ifelse(top_load$id %in% names(id_to_gene),
                             id_to_gene[top_load$id], top_load$id)
  } else {
    top_load$label <- top_load$id
  }

  go_assigned <- FALSE
  if (!is.null(annotation_df)) {
    go_map <- tryCatch(suppressWarnings({
      raw_map <- AnnotationDbi::select(
        org.Hs.eg.db, keys = unique(top_load$label),
        keytype = "SYMBOL", columns = c("SYMBOL", "GO", "ONTOLOGY")
      ) %>% filter(ONTOLOGY == "CC", !is.na(GO)) %>%
        dplyr::select(gene = SYMBOL, GO)
      if (nrow(raw_map) > 0) {
        raw_map %>% left_join(
          AnnotationDbi::select(GO.db, keys = unique(raw_map$GO),
                                columns = "TERM") %>%
            dplyr::select(GO = GOID, term = TERM), by = "GO")
      } else NULL
    }), error = function(e) NULL)

    if (!is.null(go_map) && nrow(go_map) > 0) {
      term_freq <- go_map %>% count(term, sort = TRUE)
      top_terms <- term_freq %>% slice_head(n = n_top_cc) %>% pull(term)
      go_primary <- go_map %>%
        left_join(term_freq, by = "term") %>%
        group_by(gene) %>% slice_max(n, n = 1, with_ties = FALSE) %>% ungroup() %>%
        mutate(compartment = if_else(term %in% top_terms, term, "Other")) %>%
        dplyr::select(label = gene, compartment) %>%
        distinct(label, .keep_all = TRUE)
      top_load <- top_load %>%
        left_join(go_primary, by = "label") %>%
        mutate(compartment = if_else(is.na(compartment), "No GO:CC", compartment))
      go_assigned <- TRUE
    }
  }
  if (!go_assigned) top_load$compartment <- "No GO:CC"

  score_range <- max(abs(c(pc_df$PC1, pc_df$PC2)))
  load_range  <- max(top_load$magnitude)
  sf <- score_range / load_range * 0.8
  top_load <- top_load %>% mutate(PC1s = PC1 * sf, PC2s = PC2 * sf)

  ggplot() +
    geom_point(data = pc_df,
               aes(x = PC1, y = PC2, color = Group_Time, shape = Timepoint),
               size = 3.5, alpha = 0.85) +
    scale_color_manual(values = color_pal, name = "Group") +
    scale_shape_manual(values = shape_tp) +
    ggnewscale::new_scale_color() +
    geom_segment(data = top_load,
                 aes(x = 0, y = 0, xend = PC1s, yend = PC2s, color = compartment),
                 arrow = arrow(length = unit(0.2, "cm")),
                 alpha = 0.7, linewidth = 0.5) +
    geom_text_repel(data = top_load,
                    aes(x = PC1s, y = PC2s, label = label, color = compartment),
                    size = 2.5, max.overlaps = 15) +
    scale_color_brewer(palette = "Set2", name = "GO:CC") +
    labs(x = sprintf("PC1 (%.1f%%)", var_exp[1]),
         y = sprintf("PC2 (%.1f%%)", var_exp[2]),
         title = title_text) +
    theme_minimal(base_size = 12) +
    theme(legend.position = "bottom", legend.box = "vertical")
}

# =============================================================================
# 1. LOAD DATA
# =============================================================================

raw <- read_excel(file.path(input_dir, "YvO_raw.xlsx"))
annot_cols <- c("uniprot_id", "protein", "gene", "description", "n_seq")
annotation <- raw[, annot_cols]
intensity  <- raw[, setdiff(names(raw), annot_cols)]
cat(sprintf("Raw: %d proteins x %d samples\n", nrow(raw), ncol(intensity)))

metadata <- as.data.frame(read_excel(file.path(input_dir, "YvO_meta.xlsx")))
rownames(metadata) <- metadata$Col_ID
stopifnot("Sample mismatch" = setequal(colnames(intensity), metadata$Col_ID))
intensity <- intensity[, metadata$Col_ID]

n_raw <- nrow(annotation)
filter_log <- tibble(step = "Raw input", n_before = NA_integer_,
                     n_after = n_raw, n_removed = NA_integer_)

# =============================================================================
# 2. HPA TISSUE FILTER
# =============================================================================

hpa <- read_tsv(file.path(input_dir, "HPA_skeletal_muscle_annotations.tsv"),
                show_col_types = FALSE) %>%
  dplyr::select(Gene, Ensembl, Evidence,
                Protein_class    = `Protein class`,
                Subcellular_main = `Subcellular main location`,
                Interactions) %>%
  distinct(Gene, .keep_all = TRUE)

n_before <- nrow(annotation)
keep_mask <- annotation$gene %in% hpa$Gene
intensity  <- intensity[keep_mask, ]
annotation <- annotation[keep_mask, ] %>%
  left_join(hpa, by = c("gene" = "Gene"))
n_after <- nrow(annotation)

filter_log <- bind_rows(filter_log, tibble(
  step = "HPA non-muscle removal", n_before = n_before,
  n_after = n_after, n_removed = n_before - n_after))
cat(sprintf("HPA filter: %d -> %d (%d removed)\n", n_before, n_after, n_before - n_after))

removed_genes <- setdiff(raw$gene, annotation$gene)
p_hpa_donut <- make_gocc_donut(removed_genes, "GO:CC — Proteins removed by HPA filter")

# HPA sensitivity audit
aging_patterns <- c("^COL\\d", "^MYH", "^HLA-", "^IGF", "^MSTN", "^PAX7",
                    "^FOXO", "^TRIM63", "^FBXO32", "^MT-", "^SOD", "^CAT$",
                    "^GPX", "^TNF", "^IL\\d")
hpa_audit <- tibble(gene = removed_genes) %>%
  mutate(aging_relevant = sapply(gene, function(g)
    any(sapply(aging_patterns, function(p) grepl(p, g)))))
write_csv(hpa_audit, file.path(data_dir, "02c_hpa_sensitivity_audit.csv"))

# =============================================================================
# 3. ASSEMBLE DAList
# =============================================================================

if (any(duplicated(annotation$uniprot_id))) {
  n_before_dup <- nrow(annotation)
  annotation$row_mean <- rowMeans(data.matrix(intensity), na.rm = TRUE)
  keep_idx <- annotation %>%
    mutate(row_idx = row_number()) %>%
    group_by(uniprot_id) %>%
    slice_max(row_mean, n = 1, with_ties = FALSE) %>%
    pull(row_idx)
  annotation <- annotation[keep_idx, ]; intensity <- intensity[keep_idx, ]
  annotation$row_mean <- NULL
  filter_log <- bind_rows(filter_log, tibble(
    step = "Deduplication", n_before = n_before_dup,
    n_after = nrow(annotation), n_removed = n_before_dup - nrow(annotation)))
  cat(sprintf("Deduplicated: %d unique proteins\n", nrow(annotation)))
}

intensity_mat <- as.data.frame(data.matrix(intensity))
rownames(intensity_mat) <- annotation$uniprot_id
annot_df <- as.data.frame(annotation); rownames(annot_df) <- annotation$uniprot_id
meta_df  <- as.data.frame(metadata);   rownames(meta_df)  <- metadata$Col_ID

dal <- DAList(data = intensity_mat, annotation = annot_df, metadata = meta_df)
cat(sprintf("DAList: %d proteins x %d samples\n", nrow(dal$data), ncol(dal$data)))

# =============================================================================
# 4. QUALITY FILTERING
# =============================================================================

dal <- zero_to_missing(dal)
n_before <- nrow(dal$data)

group_prop <- dal$metadata %>%
  split(.$Group_Time) %>%
  lapply(function(g) rowMeans(!is.na(dal$data[, g$Col_ID, drop = FALSE]))) %>%
  bind_cols()
keep <- apply(group_prop >= 0.66, 1, any)
dal <- filter_proteins_by_annotation(dal, keep)

n_after <- nrow(dal$data)
filter_log <- bind_rows(filter_log, tibble(
  step = "Missingness (>=66% in >=1 group)", n_before = n_before,
  n_after = n_after, n_removed = n_before - n_after))
cat(sprintf("Missingness filter: %d -> %d (%d removed)\n",
            n_before, n_after, n_before - n_after))

miss_removed_genes <- annot_df %>%
  filter(!uniprot_id %in% rownames(dal$data), !gene %in% removed_genes) %>%
  pull(gene) %>% unique()
p_miss_donut <- make_gocc_donut(miss_removed_genes,
                                "GO:CC — Proteins removed by missingness filter")

filter_log <- filter_log %>% mutate(pct_of_raw = round(n_after / n_raw * 100, 1))
write_csv(filter_log, file.path(data_dir, "03_filtering_effects.csv"))
write_csv(
  annot_df %>% filter(!uniprot_id %in% rownames(dal$data)) %>%
    dplyr::select(uniprot_id, gene, description),
  file.path(data_dir, "03_filtered_proteins.csv"))
print(filter_log)

# Filtering stacked bar
filter_plot_data <- filter_log %>%
  filter(!is.na(n_removed)) %>%
  mutate(step = factor(step, levels = step)) %>%
  pivot_longer(c(n_after, n_removed), names_to = "status", values_to = "n") %>%
  mutate(status = recode(status, n_after = "Retained", n_removed = "Removed"))

p_filter_bar <- ggplot(filter_plot_data, aes(x = step, y = n, fill = status)) +
  geom_col(width = 0.6) +
  geom_text(aes(label = n), position = position_stack(vjust = 0.5), size = 3.5) +
  scale_fill_manual(values = c(Retained = "#2166AC", Removed = "#B2182B")) +
  labs(x = NULL, y = "Proteins", fill = NULL,
       title = "Protein retention through filtering pipeline") +
  theme_minimal(base_size = 13) +
  theme(axis.text.x = element_text(angle = 25, hjust = 1), legend.position = "top")

# Per-sample missingness bars
miss_profile <- dal$metadata %>%
  dplyr::select(Col_ID, Group_Time) %>%
  mutate(n_detected = colSums(!is.na(dal$data[, Col_ID])),
         n_missing  = nrow(dal$data) - n_detected) %>%
  pivot_longer(c(n_detected, n_missing), names_to = "status", values_to = "n") %>%
  mutate(status = recode(status, n_detected = "Detected", n_missing = "Missing"))

p_miss_bars <- ggplot(miss_profile, aes(x = reorder(Col_ID, -n * (status == "Detected")),
                                         y = n, fill = status)) +
  geom_col(width = 0.8) +
  scale_fill_manual(values = c(Detected = "#2166AC", Missing = "#D6604D")) +
  facet_grid(~ Group_Time, scales = "free_x", space = "free_x") +
  labs(x = NULL, y = "Proteins", fill = NULL,
       title = "Per-sample detection vs missingness") +
  theme_minimal(base_size = 11) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 4),
        legend.position = "top", strip.text = element_text(face = "bold"))

# =============================================================================
# 5. CV VARIABILITY
# =============================================================================

# Back-transform to linear scale for CV calculation (Brenes 2024,
# J Proteome Res, PMC11629372). The base CV formula (SD/mean * 100) should
# only be applied to non-log-transformed intensities.
lin_dat <- 2^dal$data

cv_by_group <- dal$metadata %>%
  split(.$Group_Time) %>%
  lapply(function(g) {
    sub <- lin_dat[, g$Col_ID, drop = FALSE]
    tibble(uniprot_id = rownames(sub),
           cv = apply(sub, 1, function(x) {
             x <- x[!is.na(x)]
             if (length(x) < 2) return(NA_real_)
             sd(x) / mean(x)
           }),
           Group_Time = g$Group_Time[1])
  }) %>% bind_rows()

p_cv_violin <- ggplot(cv_by_group, aes(x = Group_Time, y = cv, fill = Group_Time)) +
  geom_violin(alpha = 0.7, draw_quantiles = c(0.25, 0.5, 0.75)) +
  scale_fill_manual(values = pal_group_time, guide = "none") +
  coord_cartesian(ylim = c(0, quantile(cv_by_group$cv, 0.99, na.rm = TRUE))) +
  labs(x = NULL, y = "CV (linear scale)",
       title = "Inter-individual protein variability by group") +
  theme_minimal(base_size = 12)

cv_wide <- cv_by_group %>% pivot_wider(names_from = Group_Time, values_from = cv)

plot_cv_scatter <- function(df, x_col, y_col) {
  ggplot(df, aes(x = .data[[x_col]], y = .data[[y_col]])) +
    geom_point(alpha = 0.3, size = 1, color = "gray30") +
    geom_abline(slope = 1, linetype = "dashed", color = "red", alpha = 0.5) +
    geom_smooth(method = "lm", se = FALSE, linewidth = 0.6, color = "#2166AC") +
    coord_cartesian(xlim = c(0, quantile(df[[x_col]], 0.99, na.rm = TRUE)),
                    ylim = c(0, quantile(df[[y_col]], 0.99, na.rm = TRUE))) +
    labs(x = paste("CV:", x_col), y = paste("CV:", y_col)) +
    theme_minimal(base_size = 11)
}

p_cv1 <- plot_cv_scatter(cv_wide, "Young_Pre", "Young_Post")
p_cv2 <- plot_cv_scatter(cv_wide, "Old_Pre",   "Old_Post")
p_cv3 <- plot_cv_scatter(cv_wide, "Young_Pre",  "Old_Pre")

# =============================================================================
# 6. OUTLIER DETECTION & REMOVAL
# =============================================================================

pct_missing <- colMeans(is.na(dal$data)) * 100

delta_missing <- dal$metadata %>%
  dplyr::select(Col_ID, Subject_ID, Group, Timepoint) %>%
  mutate(pct_missing = pct_missing[Col_ID],
         prefix = str_remove(Col_ID, "_(Pre|Post)$")) %>%
  arrange(prefix, match(Timepoint, c("Pre", "Post"))) %>%
  group_by(prefix) %>%
  mutate(delta_missing = abs(pct_missing - dplyr::first(pct_missing))) %>%
  ungroup()

miss_threshold  <- quantile(pct_missing, 0.75) + 1.5 * IQR(pct_missing)
delta_vals      <- delta_missing$delta_missing[delta_missing$Timepoint != "Pre"]
delta_threshold <- quantile(delta_vals, 0.75) + 1.5 * IQR(delta_vals)

delta_missing <- delta_missing %>%
  mutate(miss_flag = pct_missing > miss_threshold |
           (Timepoint != "Pre" & delta_missing > delta_threshold))

zero_miss_mat <- dal$data[rowSums(is.na(dal$data)) == 0, ]
cat(sprintf("Pre-norm PCA: %d complete proteins (of %d)\n",
            nrow(zero_miss_mat), nrow(dal$data)))
pca_out   <- run_pca(zero_miss_mat, dal$metadata, log_transform = TRUE)
pc_scores <- pca_out$pca$x[, 1:3]
mahal_dist <- mahalanobis(pc_scores, colMeans(pc_scores), cov(pc_scores))
pca_flags <- tibble(Col_ID = colnames(dal$data), mahal_dist = mahal_dist,
                    pca_flag = mahal_dist > qchisq(0.99, df = 3))

sample_medians <- apply(log2(dal$data + 1), 2, median, na.rm = TRUE)
global_median  <- median(sample_medians)
mad_val        <- mad(sample_medians)
mad_flags <- tibble(Col_ID = names(sample_medians), sample_median = sample_medians,
                    mad_flag = abs(sample_medians - global_median) > 3 * mad_val)

outlier_diag <- delta_missing %>%
  left_join(pca_flags, by = "Col_ID") %>%
  left_join(mad_flags, by = "Col_ID") %>%
  mutate(n_flags = miss_flag + pca_flag + mad_flag,
         consensus_outlier = n_flags >= 3)

write_csv(outlier_diag, file.path(data_dir, "04_outlier_diagnostics.csv"))
n_outliers <- sum(outlier_diag$consensus_outlier)
cat(sprintf("Outlier consensus: %d sample(s) flagged (3/3 methods)\n", n_outliers))

# Outlier plots
p_out_a <- ggplot(outlier_diag, aes(pct_missing, delta_missing)) +
  geom_point(aes(color = consensus_outlier, shape = Timepoint), size = 3) +
  geom_text_repel(data = . %>% filter(consensus_outlier),
                  aes(label = prefix), size = 2.5, color = "red", max.overlaps = 20) +
  geom_hline(yintercept = delta_threshold, linetype = "dashed", color = "red", alpha = 0.5) +
  geom_vline(xintercept = miss_threshold, linetype = "dashed", color = "red", alpha = 0.5) +
  scale_color_manual(values = c("FALSE" = "gray40", "TRUE" = "red")) +
  scale_shape_manual(values = shape_tp) +
  labs(x = "% Missing", y = "|Delta Missingness|", title = "A: Paired Missingness") +
  theme_minimal(base_size = 11) + theme(legend.position = "bottom")

p_out_b <- plot_pca_biplot(pca_out, color_pal = pal_group_time,
                           annotation_df = dal$annotation,
                           title_text = "B: PCA Biplot (pre-normalization)")

p_out_c <- ggplot(outlier_diag, aes(reorder(prefix, sample_median), sample_median)) +
  geom_point(aes(color = consensus_outlier), size = 2.5) +
  geom_text_repel(data = . %>% filter(consensus_outlier),
                  aes(label = prefix), size = 2.5, color = "red") +
  geom_hline(yintercept = global_median) +
  geom_hline(yintercept = global_median + c(-3, 3) * mad_val,
             linetype = "dashed", color = "red", alpha = 0.5) +
  scale_color_manual(values = c("FALSE" = "gray40", "TRUE" = "red")) +
  labs(x = "Sample", y = "Median log2 intensity", title = "C: MAD Outliers") +
  theme_minimal(base_size = 11) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 4),
        legend.position = "none")

# Remove outliers (paired)
if (n_outliers > 0) {
  flagged_ids <- outlier_diag %>% filter(consensus_outlier) %>% pull(Col_ID)
  flagged_prefixes <- unique(str_remove(flagged_ids, "_(Pre|Post)$"))
  remove_ids <- dal$metadata$Col_ID[
    str_remove(dal$metadata$Col_ID, "_(Pre|Post)$") %in% flagged_prefixes]
  cat(sprintf("Removing (with pairs): %s\n", paste(flagged_prefixes, collapse = ", ")))
  dal <- filter_samples(dal, !(Col_ID %in% remove_ids))
  cat(sprintf("%d samples remain\n", ncol(dal$data)))
}

# =============================================================================
# 7. NORMALIZATION & proteoDA REPORTS
# =============================================================================

# Report 1: normalization method comparison
write_norm_report(dal, grouping_column = "Group_Time", output_dir = report_dir,
                  filename = "01_normalization_report.pdf", overwrite = TRUE)
saveRDS(dal, file.path(data_dir, "00_DAList_prenorm.rds"))

# Report 2: QC pre-normalization
write_qc_report(dal, color_column = "Group_Time", label_column = "Col_ID",
                output_dir = report_dir, filename = "02_qc_report_pre.pdf",
                overwrite = TRUE)

# Normalize
dal <- normalize_data(dal, norm_method = "cycloess")
cat(sprintf("Cycloess: %d proteins x %d samples\n", nrow(dal$data), ncol(dal$data)))

# Report 3: QC post-normalization
write_qc_report(dal, color_column = "Group_Time", label_column = "Col_ID",
                output_dir = report_dir, filename = "03_qc_report_post.pdf",
                overwrite = TRUE)

# =============================================================================
# 8. REPORT 4: FILTERING REPORT (combined)
# =============================================================================

pdf(file.path(report_dir, "04_filtering_report.pdf"), width = 10, height = 10)

if (!is.null(p_hpa_donut)) print(p_hpa_donut)
print(p_filter_bar)
if (!is.null(p_miss_donut)) print(p_miss_donut)
print(p_miss_bars)

dev.off()

# =============================================================================
# 9. REPORT 5: DIAGNOSTICS REPORT (combined)
# =============================================================================

# Post-norm PCA (Group x Timepoint only)
pca_post <- run_pca(dal$data, dal$metadata, log_transform = FALSE)

p_pca_gt <- ggplot(pca_post$pc_df, aes(PC1, PC2, color = Group_Time, shape = Timepoint)) +
  geom_point(size = 3.5, alpha = 0.85) +
  stat_ellipse(aes(group = Group_Time), type = "norm", level = 0.68, linewidth = 0.7) +
  scale_color_manual(values = pal_group_time) + scale_shape_manual(values = shape_tp) +
  labs(x = sprintf("PC1 (%.1f%%)", pca_post$var_exp[1]),
       y = sprintf("PC2 (%.1f%%)", pca_post$var_exp[2]),
       title = "Post-normalization PCA (Group x Timepoint)") +
  theme_minimal(base_size = 12) + theme(legend.position = "bottom")

pdf(file.path(report_dir, "05_diagnostics_report.pdf"), width = 14, height = 10)

# Page 1: CV variability
print((p_cv_violin | plot_spacer()) / (p_cv1 | p_cv2 | p_cv3) +
        plot_annotation(title = "Inter-Individual Variability"))

# Page 2: Outlier diagnostics
print(p_out_a / p_out_b / p_out_c +
        plot_annotation(title = "Outlier Diagnostics (3-method consensus)"))

# Page 3: Post-norm PCA
print(p_pca_gt)

dev.off()

# =============================================================================
# 10. EXPORT
# =============================================================================

export_df <- bind_cols(
  as_tibble(dal$annotation) %>% dplyr::select(uniprot_id, protein, gene, description),
  as_tibble(dal$data))

write_csv(export_df, file.path(data_dir, "01_normalized.csv"))
saveRDS(dal, file.path(data_dir, "01_DAList_normalized.rds"))

cat(sprintf("Exported: %d proteins x %d samples\n", nrow(dal$data), ncol(dal$data)))
cat("\n=== YvO Normalization complete ===\n")
