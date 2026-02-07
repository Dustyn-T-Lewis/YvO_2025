#!/usr/bin/env Rscript
# ============================================================================
# YvO Proteomics — Normalization Pipeline
# Young vs. Old Skeletal Muscle DIA-MS Proteomics —
# Resistance Training Intervention
#
# Run from: 01_normalization/a_script/
# Usage:    source("YvO_normalization_run.R")  OR  Rscript YvO_normalization_run.R
# ============================================================================

cat("=== YvO Normalization Pipeline ===\n\n")

# --- 0: Setup ----------------------------------------------------------------
cat(">> 0 — Setup\n")

if (!requireNamespace("pacman", quietly = TRUE)) install.packages("pacman")
pacman::p_load(
  proteoDA,
  readxl, readr, dplyr, tidyr, stringr,
  org.Hs.eg.db, GO.db, AnnotationDbi,
  ggplot2, ggrepel, patchwork
)

base_dir   <- normalizePath(file.path(dirname(getwd()), ".."), mustWork = TRUE)
input_dir  <- file.path(base_dir, "00_input")
report_dir <- file.path(base_dir, "01_normalization", "b_reports")
data_dir   <- file.path(base_dir, "01_normalization", "c_data")

dir.create(report_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(data_dir,   recursive = TRUE, showWarnings = FALSE)

pal_group      <- c(Young = "#2166AC", Old = "#B2182B")
pal_group_time <- c(Young_Pre = "#92C5DE", Young_Post = "#2166AC",
                    Old_Pre   = "#F4A582", Old_Post   = "#B2182B")
shape_tp       <- c(Pre = 16, Post = 17)

# --- Helper functions ---------------------------------------------------------

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

plot_gocc_donut <- function(gene_vec, title_text, output_path, n_top = 8) {
  if (length(gene_vec) == 0) { cat("   No genes to plot.\n"); return(invisible(NULL)) }

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

  p <- ggplot(donut_data, aes(ymax = ymax, ymin = ymin, xmax = 4, xmin = 2.5,
                              fill = category)) +
    geom_rect(color = "white", linewidth = 0.5) +
    geom_text(aes(x = 3.25, y = ymid, label = label), size = 2.5, lineheight = 0.9) +
    annotate("text", x = 0, y = 0,
             label = paste0(length(gene_vec), "\nproteins"),
             size = 5, fontface = "bold") +
    coord_polar(theta = "y") + xlim(c(0, 4.5)) +
    scale_fill_brewer(palette = "Set2", guide = "none") +
    theme_void() + labs(title = title_text)

  ggsave(output_path, p, width = 8, height = 8)
  cat(sprintf("   Saved: %s\n", basename(output_path)))
}

plot_pca_biplot <- function(pca_out, color_col, color_pal, shape_col = "Timepoint",
                            n_loadings = 10, title_text = "PCA Biplot") {
  pc_df   <- pca_out$pc_df
  pca     <- pca_out$pca
  var_exp <- pca_out$var_exp

  # Top loadings by L2 norm on PC1:2
  loadings <- as.data.frame(pca$rotation[, 1:2])
  loadings$magnitude <- sqrt(loadings$PC1^2 + loadings$PC2^2)
  loadings$label <- rownames(loadings)
  top_load <- loadings %>% slice_max(magnitude, n = n_loadings)

  # Scale loadings to score range
  score_range <- max(abs(c(pc_df$PC1, pc_df$PC2)))
  load_range  <- max(top_load$magnitude)
  scale_factor <- score_range / load_range * 0.8

  top_load <- top_load %>% mutate(PC1s = PC1 * scale_factor,
                                   PC2s = PC2 * scale_factor)

  ggplot() +
    geom_point(data = pc_df,
               aes(x = PC1, y = PC2, color = .data[[color_col]],
                   shape = .data[[shape_col]]),
               size = 3.5, alpha = 0.85) +
    geom_segment(data = top_load,
                 aes(x = 0, y = 0, xend = PC1s, yend = PC2s),
                 arrow = arrow(length = unit(0.2, "cm")),
                 color = "firebrick", alpha = 0.6, linewidth = 0.5) +
    geom_text_repel(data = top_load,
                    aes(x = PC1s, y = PC2s, label = label),
                    size = 2.5, color = "firebrick", max.overlaps = 15) +
    scale_color_manual(values = color_pal) +
    scale_shape_manual(values = shape_tp) +
    labs(x = sprintf("PC1 (%.1f%%)", var_exp[1]),
         y = sprintf("PC2 (%.1f%%)", var_exp[2]),
         title = title_text) +
    theme_minimal(base_size = 12) + theme(legend.position = "bottom")
}

cat("   Base directory:", base_dir, "\n")

# --- 1: Load & Validate Data -------------------------------------------------
cat("\n>> 1 — Load & Validate Data\n")

raw <- read_excel(file.path(input_dir, "YvO_raw.xlsx"))

annot_cols <- c("uniprot_id", "protein", "gene", "description", "n_seq")
annotation <- raw[, annot_cols]
intensity  <- raw[, setdiff(names(raw), annot_cols)]

cat(sprintf("   Raw data: %d proteins x %d samples\n", nrow(raw), ncol(intensity)))

metadata <- as.data.frame(read_excel(file.path(input_dir, "YvO_meta.xlsx")))
rownames(metadata) <- metadata$Col_ID

stopifnot("Sample mismatch" = setequal(colnames(intensity), metadata$Col_ID))
intensity <- intensity[, metadata$Col_ID]

cat(sprintf("   Metadata: %d samples aligned.\n", nrow(metadata)))
cat("\n   Sample distribution:\n")
print(metadata %>% count(Group, Timepoint, Group_Time))

# Initialize filter log
n_raw <- nrow(annotation)
filter_log <- tibble(step = "Raw input", n_before = NA_integer_,
                     n_after = n_raw, n_removed = NA_integer_)

# --- 2: HPA Annotation & Non-Muscle Filtering --------------------------------
cat("\n>> 2 — HPA Annotation & Non-Muscle Filtering\n")

hpa <- read_tsv(file.path(input_dir, "HPA_skeletal_muscle_annotations.tsv"),
                show_col_types = FALSE) %>%
  dplyr::select(Gene, Ensembl, Evidence,
                Protein_class    = `Protein class`,
                Subcellular_main = `Subcellular main location`,
                Interactions) %>%
  distinct(Gene, .keep_all = TRUE)

n_before <- nrow(annotation)
annotation <- annotation %>% inner_join(hpa, by = c("gene" = "Gene"))
n_after <- nrow(annotation)
intensity <- intensity[seq_len(n_after), ]

filter_log <- bind_rows(filter_log, tibble(
  step = "HPA non-muscle removal", n_before = n_before,
  n_after = n_after, n_removed = n_before - n_after))

cat(sprintf("   HPA join: %d -> %d proteins (%d removed)\n",
            n_before, n_after, n_before - n_after))

# --- 2b: GO:CC donut for removed proteins ------------------------------------
removed_genes <- setdiff(raw$gene, annotation$gene)
plot_gocc_donut(removed_genes, "GO:CC — Proteins removed by HPA filter",
                file.path(report_dir, "02b_removed_proteins_gocc.pdf"))

# --- 2c: HPA filter sensitivity audit ----------------------------------------
cat("\n>> 2c — HPA filter sensitivity audit\n")

aging_patterns <- c(
  "^COL\\d",   # collagens (ECM remodeling)
  "^MYH",      # myosin heavy chains (fiber type)
  "^HLA-",     # MHC class I/II (immune)
  "^IGF",      # insulin-like growth factors
  "^MSTN",     # myostatin
  "^PAX7",     # satellite cell marker
  "^FOXO",     # atrophy signaling
  "^TRIM63",   # MuRF1 (ubiquitin ligase)
  "^FBXO32",   # atrogin-1 (ubiquitin ligase)
  "^MT-",      # mitochondrial-encoded
  "^SOD",      # superoxide dismutases
  "^CAT$",     # catalase
  "^GPX",      # glutathione peroxidases
  "^TNF",      # tumor necrosis factor
  "^IL\\d"     # interleukins
)

hpa_audit <- tibble(gene = removed_genes) %>%
  mutate(aging_relevant = sapply(gene, function(g) {
    any(sapply(aging_patterns, function(p) grepl(p, g)))
  }))

n_flagged <- sum(hpa_audit$aging_relevant)
cat(sprintf("   %d of %d removed proteins match aging-relevant patterns\n",
            n_flagged, length(removed_genes)))
if (n_flagged > 0) {
  cat("   Flagged genes:", paste(hpa_audit$gene[hpa_audit$aging_relevant], collapse = ", "), "\n")
}

write_csv(hpa_audit, file.path(data_dir, "02c_hpa_sensitivity_audit.csv"))

# --- 3: Assemble DAList -------------------------------------------------------
cat("\n>> 3 — Assemble DAList\n")

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
  cat(sprintf("   Deduplicated: %d unique proteins\n", nrow(annotation)))
} else {
  cat("   No duplicate uniprot_ids.\n")
}

intensity_mat <- as.data.frame(data.matrix(intensity))
rownames(intensity_mat) <- annotation$uniprot_id
annot_df <- as.data.frame(annotation); rownames(annot_df) <- annotation$uniprot_id
meta_df  <- as.data.frame(metadata);   rownames(meta_df)  <- metadata$Col_ID

dal <- DAList(data = intensity_mat, annotation = annot_df, metadata = meta_df)
cat(sprintf("   DAList: %d proteins x %d samples\n", nrow(dal$data), ncol(dal$data)))

# --- 4: Quality Filtering ----------------------------------------------------
cat("\n>> 4 — Quality Filtering\n")

dal <- zero_to_missing(dal)

# Missingness filter: >=66% present in at least 1 Group_Time
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

cat(sprintf("   Missingness filter: %d -> %d (%d removed)\n",
            n_before, n_after, n_before - n_after))

# GO:CC donut for missingness-filtered proteins
miss_removed_genes <- annot_df %>%
  filter(!uniprot_id %in% rownames(dal$data),
         !annot_df$gene %in% removed_genes) %>%  # exclude HPA-removed
  pull(gene) %>% unique()
plot_gocc_donut(miss_removed_genes, "GO:CC — Proteins removed by missingness filter",
                file.path(report_dir, "04a_missingness_filtered_gocc.pdf"))

# Filter log finalize
filter_log <- filter_log %>%
  mutate(pct_of_raw = round(n_after / n_raw * 100, 1))
write_csv(filter_log, file.path(data_dir, "03_filtering_effects.csv"))

filtered_proteins <- annot_df %>%
  filter(!uniprot_id %in% rownames(dal$data)) %>%
  dplyr::select(uniprot_id, gene, description)
write_csv(filtered_proteins, file.path(data_dir, "03_filtered_proteins.csv"))

cat("\n   Filtering summary:\n"); print(filter_log)

# Stacked bar: filtering pipeline
filter_plot_data <- filter_log %>%
  filter(!is.na(n_removed)) %>%
  mutate(step = factor(step, levels = step)) %>%
  pivot_longer(c(n_after, n_removed), names_to = "status", values_to = "n") %>%
  mutate(status = recode(status, n_after = "Retained", n_removed = "Removed"))

p_filter <- ggplot(filter_plot_data, aes(x = step, y = n, fill = status)) +
  geom_col(width = 0.6) +
  geom_text(aes(label = n), position = position_stack(vjust = 0.5), size = 3.5) +
  scale_fill_manual(values = c(Retained = "#2166AC", Removed = "#B2182B")) +
  labs(x = NULL, y = "Proteins", fill = NULL,
       title = "Protein retention through filtering pipeline") +
  theme_minimal(base_size = 13) +
  theme(axis.text.x = element_text(angle = 25, hjust = 1), legend.position = "top")
ggsave(file.path(report_dir, "03_filtering_effects.pdf"), p_filter, width = 8, height = 5)

# --- 4b: Subcellular composition & Missingness profile ------------------------
cat("\n>> 4b — Composition & Missingness Profiles\n")

# Subcellular composition per Group_Time (detected proteins only)
sc_annot <- dal$annotation %>%
  mutate(compartment = str_trim(str_extract(Subcellular_main, "^[^,]+")),
         compartment = if_else(is.na(compartment) | compartment == "",
                               "Unannotated", compartment))

# Bin rare compartments into "Other"
top_compartments <- sc_annot %>% count(compartment, sort = TRUE) %>%
  slice_head(n = 10) %>% pull(compartment)
sc_annot <- sc_annot %>%
  mutate(compartment = if_else(compartment %in% top_compartments,
                               compartment, "Other"))

# Per group_time: fraction of detected proteins per compartment
sc_group <- dal$metadata %>%
  split(.$Group_Time) %>%
  lapply(function(g) {
    detected <- rownames(dal$data)[rowSums(!is.na(dal$data[, g$Col_ID, drop = FALSE])) > 0]
    sc_annot %>% filter(uniprot_id %in% detected) %>%
      count(compartment) %>%
      mutate(frac = n / sum(n), Group_Time = g$Group_Time[1])
  }) %>% bind_rows()

p_cc <- ggplot(sc_group, aes(x = Group_Time, y = frac, fill = compartment)) +
  geom_col(position = "stack", width = 0.7) +
  scale_fill_brewer(palette = "Set3") +
  scale_y_continuous(labels = scales::percent_format()) +
  labs(x = NULL, y = "Fraction of detected proteins", fill = "Compartment",
       title = "Subcellular composition of detected proteome") +
  theme_minimal(base_size = 12) +
  theme(legend.position = "right", axis.text.x = element_text(angle = 25, hjust = 1))

# Per-sample missingness stacked bars
miss_profile <- dal$metadata %>%
  dplyr::select(Col_ID, Group_Time) %>%
  mutate(n_detected = colSums(!is.na(dal$data[, Col_ID])),
         n_missing  = nrow(dal$data) - n_detected) %>%
  pivot_longer(c(n_detected, n_missing), names_to = "status", values_to = "n") %>%
  mutate(status = recode(status, n_detected = "Detected", n_missing = "Missing"))

p_miss <- ggplot(miss_profile, aes(x = reorder(Col_ID, -n * (status == "Detected")),
                                    y = n, fill = status)) +
  geom_col(width = 0.8) +
  scale_fill_manual(values = c(Detected = "#2166AC", Missing = "#D6604D")) +
  facet_grid(~ Group_Time, scales = "free_x", space = "free_x") +
  labs(x = NULL, y = "Proteins", fill = NULL,
       title = "Per-sample detection vs. missingness") +
  theme_minimal(base_size = 11) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 4),
        legend.position = "top", strip.text = element_text(face = "bold"))

ggsave(file.path(report_dir, "04b_composition_missingness.pdf"),
       p_cc / p_miss, width = 12, height = 10)
cat("   Saved: 04b_composition_missingness.pdf\n")

# --- 4c: Inter-individual variability & CV scatters ---------------------------
cat("\n>> 4c — Variability & CV Profiles\n")

log2_dat <- log2(dal$data + 1)

# Compute per-protein CV within each Group_Time
cv_by_group <- dal$metadata %>%
  split(.$Group_Time) %>%
  lapply(function(g) {
    sub <- log2_dat[, g$Col_ID, drop = FALSE]
    mu <- rowMeans(sub, na.rm = TRUE)
    sd <- apply(sub, 1, sd, na.rm = TRUE)
    tibble(uniprot_id = rownames(sub), cv = sd / mu, Group_Time = g$Group_Time[1])
  }) %>% bind_rows()

# Violin: CV distributions per group
p_cv_violin <- ggplot(cv_by_group, aes(x = Group_Time, y = cv, fill = Group_Time)) +
  geom_violin(alpha = 0.7, draw_quantiles = c(0.25, 0.5, 0.75)) +
  scale_fill_manual(values = pal_group_time, guide = "none") +
  coord_cartesian(ylim = c(0, quantile(cv_by_group$cv, 0.99, na.rm = TRUE))) +
  labs(x = NULL, y = "Coefficient of Variation",
       title = "Inter-individual protein variability by group") +
  theme_minimal(base_size = 12)

# CV scatter helper
cv_wide <- cv_by_group %>%
  pivot_wider(names_from = Group_Time, values_from = cv)

plot_cv_scatter <- function(df, x_col, y_col) {
  ggplot(df, aes(x = .data[[x_col]], y = .data[[y_col]])) +
    geom_point(alpha = 0.3, size = 1, color = "gray30") +
    geom_abline(slope = 1, linetype = "dashed", color = "red", alpha = 0.5) +
    geom_smooth(method = "lm", se = FALSE, linewidth = 0.6, color = "#2166AC") +
    coord_cartesian(xlim = c(0, quantile(df[[x_col]], 0.99, na.rm = TRUE)),
                    ylim = c(0, quantile(df[[y_col]], 0.99, na.rm = TRUE))) +
    labs(x = paste("CV:", x_col), y = paste("CV:", y_col),
         title = sprintf("%s vs %s", x_col, y_col)) +
    theme_minimal(base_size = 11)
}

p_cv1 <- plot_cv_scatter(cv_wide, "Young_Pre", "Young_Post")
p_cv2 <- plot_cv_scatter(cv_wide, "Old_Pre",   "Old_Post")
p_cv3 <- plot_cv_scatter(cv_wide, "Young_Pre",  "Old_Pre")

ggsave(file.path(report_dir, "04c_variability_cv.pdf"),
       (p_cv_violin | plot_spacer()) / (p_cv1 | p_cv2 | p_cv3),
       width = 14, height = 10)
cat("   Saved: 04c_variability_cv.pdf\n")

# --- 5: Outlier Detection & Removal ------------------------------------------
cat("\n>> 5 — Outlier Detection & Removal\n")

pct_missing <- colMeans(is.na(dal$data)) * 100

# Paired missingness deltas
delta_missing <- dal$metadata %>%
  dplyr::select(Col_ID, Subject_ID, Group, Timepoint) %>%
  mutate(pct_missing = pct_missing[Col_ID],
         prefix = str_remove(Col_ID, "_(Pre|Post)$")) %>%
  arrange(prefix, match(Timepoint, c("Pre", "Post"))) %>%
  group_by(prefix) %>%
  mutate(delta_missing = pct_missing - dplyr::first(pct_missing)) %>%
  ungroup()

# IQR thresholds
miss_threshold  <- quantile(pct_missing, 0.75) + 1.5 * IQR(pct_missing)
delta_vals      <- delta_missing$delta_missing[delta_missing$Timepoint != "Pre"]
delta_threshold <- quantile(delta_vals, 0.75) + 1.5 * IQR(delta_vals)
delta_lower     <- quantile(delta_vals, 0.25) - 1.5 * IQR(delta_vals)

delta_missing <- delta_missing %>%
  mutate(miss_flag = pct_missing > miss_threshold |
           (Timepoint != "Pre" &
              (delta_missing > delta_threshold | delta_missing < delta_lower)))

# PCA outlier detection (Mahalanobis on PC1:3)
pca_out   <- run_pca(dal$data, dal$metadata, log_transform = TRUE)
pc_scores <- pca_out$pca$x[, 1:3]
mahal_dist <- mahalanobis(pc_scores, colMeans(pc_scores), cov(pc_scores))

pca_flags <- tibble(Col_ID = colnames(dal$data), mahal_dist = mahal_dist,
                    pca_flag = mahal_dist > qchisq(0.99, df = 3))

# MAD intensity outlier detection
sample_medians <- apply(log2(dal$data + 1), 2, median, na.rm = TRUE)
global_median  <- median(sample_medians)
mad_val        <- mad(sample_medians)

mad_flags <- tibble(Col_ID = names(sample_medians), sample_median = sample_medians,
                    mad_flag = abs(sample_medians - global_median) > 3 * mad_val)

# Consensus (all 3 methods must agree)
outlier_diag <- delta_missing %>%
  left_join(pca_flags, by = "Col_ID") %>%
  left_join(mad_flags, by = "Col_ID") %>%
  mutate(n_flags = miss_flag + pca_flag + mad_flag,
         consensus_outlier = n_flags >= 3)

write_csv(outlier_diag, file.path(data_dir, "04_outlier_diagnostics.csv"))
n_outliers <- sum(outlier_diag$consensus_outlier)
cat(sprintf("   Outlier consensus: %d sample(s) flagged (3/3 methods)\n", n_outliers))

if (n_outliers > 0) {
  print(outlier_diag %>% filter(consensus_outlier) %>%
          dplyr::select(prefix, Group, Timepoint, pct_missing,
                        delta_missing, mahal_dist, n_flags))
}

# Diagnostic plots
var_explained <- pca_out$var_exp

# A: Paired missingness
p1 <- ggplot(outlier_diag, aes(pct_missing, delta_missing)) +
  geom_point(aes(color = consensus_outlier, shape = Timepoint), size = 3) +
  geom_text_repel(data = . %>% filter(consensus_outlier),
                  aes(label = prefix), size = 2.5, color = "red", max.overlaps = 20) +
  geom_hline(yintercept = c(delta_lower, delta_threshold), linetype = "dashed",
             color = "red", alpha = 0.5) +
  geom_vline(xintercept = miss_threshold, linetype = "dashed",
             color = "red", alpha = 0.5) +
  scale_color_manual(values = c("FALSE" = "gray40", "TRUE" = "red")) +
  scale_shape_manual(values = shape_tp) +
  labs(x = "% Missing", y = "Delta Missing (Post-Pre)",
       title = "A: Paired Missingness") +
  theme_minimal(base_size = 11) + theme(legend.position = "bottom")

# B: PCA biplot with top loadings
p2 <- plot_pca_biplot(pca_out, color_col = "Group_Time",
                      color_pal = pal_group_time,
                      title_text = "B: PCA Biplot (pre-normalization)")

# C: MAD outliers
p3 <- ggplot(outlier_diag, aes(reorder(prefix, sample_median), sample_median)) +
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

ggsave(file.path(report_dir, "05_outlier_diagnostics.pdf"),
       p1 / p2 / p3, width = 10, height = 16)
cat("   Saved: 05_outlier_diagnostics.pdf\n")

# Remove outliers (paired: both timepoints for flagged subjects)
if (n_outliers > 0) {
  flagged_ids <- outlier_diag %>% filter(consensus_outlier) %>% pull(Col_ID)
  flagged_prefixes <- unique(str_remove(flagged_ids, "_(Pre|Post)$"))
  remove_ids <- dal$metadata$Col_ID[
    str_remove(dal$metadata$Col_ID, "_(Pre|Post)$") %in% flagged_prefixes]

  cat(sprintf("   Flagged: %s\n", paste(flagged_prefixes, collapse = ", ")))
  cat(sprintf("   Removing (with pairs): %s\n", paste(remove_ids, collapse = ", ")))
  dal <- filter_samples(dal, !(Col_ID %in% remove_ids))
  cat(sprintf("   %d samples remain\n", ncol(dal$data)))
} else {
  cat("   No outliers removed.\n")
}

# --- 6: Normalization ---------------------------------------------------------
cat("\n>> 6 — Normalization\n")

write_norm_report(dal, grouping_column = "Group_Time", output_dir = report_dir,
                  filename = "06_normalization_report.pdf", overwrite = TRUE)
dal <- normalize_data(dal, norm_method = "cycloess")
cat(sprintf("   Cycloess normalization: %d proteins x %d samples\n",
            nrow(dal$data), ncol(dal$data)))

# --- 7: Post-Normalization QC -------------------------------------------------
cat("\n>> 7 — Post-Normalization QC\n")

write_qc_report(dal, color_column = "Group_Time", label_column = "Col_ID",
                output_dir = report_dir, filename = "07_qc_report.pdf",
                overwrite = TRUE)

# Post-norm PCA
pca_post <- run_pca(dal$data, dal$metadata, log_transform = FALSE)

# PCA by Group
p_pca1 <- ggplot(pca_post$pc_df, aes(PC1, PC2, color = Group, shape = Timepoint)) +
  geom_point(size = 3.5, alpha = 0.85) +
  stat_ellipse(aes(group = Group), type = "norm", level = 0.68, linewidth = 0.7) +
  scale_color_manual(values = pal_group) + scale_shape_manual(values = shape_tp) +
  labs(x = sprintf("PC1 (%.1f%%)", pca_post$var_exp[1]),
       y = sprintf("PC2 (%.1f%%)", pca_post$var_exp[2]),
       title = "A: By Group") +
  theme_minimal(base_size = 12) + theme(legend.position = "bottom")

# PCA by Group_Time
p_pca2 <- ggplot(pca_post$pc_df, aes(PC1, PC2, color = Group_Time, shape = Timepoint)) +
  geom_point(size = 3.5, alpha = 0.85) +
  stat_ellipse(aes(group = Group_Time), type = "norm", level = 0.68, linewidth = 0.7) +
  scale_color_manual(values = pal_group_time) + scale_shape_manual(values = shape_tp) +
  labs(x = sprintf("PC1 (%.1f%%)", pca_post$var_exp[1]),
       y = sprintf("PC2 (%.1f%%)", pca_post$var_exp[2]),
       title = "B: By Group x Timepoint") +
  theme_minimal(base_size = 12) + theme(legend.position = "bottom")

# PCA biplot post-norm
p_pca3 <- plot_pca_biplot(pca_post, color_col = "Group_Time",
                          color_pal = pal_group_time,
                          title_text = "C: Post-normalization biplot")

ggsave(file.path(report_dir, "07_post_norm_pca.pdf"),
       p_pca1 / p_pca2 / p_pca3, width = 10, height = 16)
cat("   Saved: 07_post_norm_pca.pdf\n")

# --- 8: Export ----------------------------------------------------------------
cat("\n>> 8 — Export\n")

export_df <- bind_cols(
  as_tibble(dal$annotation) %>% dplyr::select(uniprot_id, protein, gene, description),
  as_tibble(dal$data))

write_csv(export_df, file.path(data_dir, "01_normalized.csv"))
saveRDS(dal, file.path(data_dir, "01_DAList_normalized.rds"))

cat(sprintf("   Exported: %d proteins x %d samples\n", nrow(dal$data), ncol(dal$data)))

# --- 9: Session Info ----------------------------------------------------------
cat("\n>> 9 — Session Info\n"); sessionInfo()
cat("\n=== YvO Normalization Pipeline Complete ===\n")
