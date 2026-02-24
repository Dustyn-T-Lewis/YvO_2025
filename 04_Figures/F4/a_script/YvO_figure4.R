################################################################################
#   Figure 4: Proteomic Response Archetypes to Resistance Training (Mfuzz FCM)
#
#   Study: YvO — DIA-MS Proteomics of Skeletal Muscle
#   Design: 2x2 mixed (Age x Time) with repeated measures
#
#   Central question: Are there distinct proteomic response archetypes to
#   resistance training, and do young and old adults differ in which
#   archetypes they express?
#
#   Approach:
#     1. Compute per-subject delta matrix (Post - Pre) across all proteins
#     2. Z-score standardize per protein
#     3. Fuzzy c-means (Mfuzz) clustering on the delta matrix
#     4. Characterize each cluster's age-specific trajectory and pathway
#        enrichment
#
#   Panels (redesigned 2026-02-23):
#     A — FCM cluster profiles (spaghetti + ribbon)
#     B — PCA highlight-on-grey scatter
#     C — Per-cluster Sankey triptych (Hallmark | GO:BP | GO:CC)
#     D — Cluster synthesis Sankey
#
#   References:
#     Kumar & Futschik 2007, Bioinformatics 23:1418 (Mfuzz)
#     Schwammle & Jensen 2010, Bioinformatics 26:2841 (mestimate)
################################################################################

# === 1. PACKAGES ==============================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(patchwork)
  library(Biobase)
  library(e1071)
  library(clusterProfiler)
  library(msigdbr)
  library(org.Hs.eg.db)
  library(rrvgo)
  library(GOSemSim)
  library(ggrepel)
  library(scales)
  library(grid)
})

# Load Mfuzz core functions without tcltk/X11 dependency
# (Mfuzz imports tcltk for its GUI; we only need the clustering functions)
.mfuzz_env <- new.env()
suppressWarnings(
  lazyLoad(file.path(find.package("Mfuzz"), "R", "Mfuzz"), envir = .mfuzz_env)
)
mestimate   <- .mfuzz_env$mestimate
mfuzz       <- .mfuzz_env$mfuzz
standardise <- .mfuzz_env$standardise

# === 2. SEED ==================================================================

set.seed(42)

# === 3. PATH RESOLUTION =======================================================

.script_dir <- tryCatch(dirname(normalizePath(sys.frame(1)$ofile)),
                         error = function(e) {
                           args <- commandArgs(trailingOnly = FALSE)
                           f <- grep("^--file=", args, value = TRUE)
                           if (length(f)) dirname(normalizePath(sub("^--file=", "", f[1])))
                           else normalizePath(".")
                         })
BASE_DIR <- normalizePath(file.path(.script_dir, "..", "..", ".."))
FIG_DIR  <- normalizePath(file.path(.script_dir, ".."))
RPT_DIR  <- file.path(FIG_DIR, "b_reports")
DAT_DIR  <- file.path(FIG_DIR, "c_data")

# === 4. DIRECTORY CREATION ====================================================

dir.create(RPT_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(DAT_DIR, recursive = TRUE, showWarnings = FALSE)

# === 5. CANONICAL CONSTANTS ===================================================

CONTRAST_COLORS <- c(Aging = "#4CAF50", Training_Young = "#E05A4E",
                     Training_Old = "#5DA5DA", Interaction = "#9B7FBF")
AGE_COLORS <- c(Young = "#4393C3", Old = "#D6604D")
DIR_COLORS <- c(Up = "#D6604D", Down = "#4393C3")
DB_COLORS  <- c(Hallmark = "#AA336A", "GO:BP" = "#00796B",
                "GO:CC" = "#26A69A", "GO:MF" = "#CD5C5C")
CLUSTER_COLORS <- c(C1 = "#E74C3C", C2 = "#3498DB", C3 = "#2ECC71",
                     C4 = "#F39C12", C5 = "#9B59B6", C6 = "#1ABC9C",
                     C7 = "#E67E22", C8 = "#34495E", C9 = "#D35400",
                     C10 = "#7F8C8D")
YOUNG_COL <- "#4393C3"
OLD_COL   <- "#D6604D"
AGING_GAP_LINE <- "#7B68EE"
KEY_TEXT  <- 2.2
KEY_TITLE <- 2.3

GROUP_COLS <- c("Young_Pre", "Young_Post", "Old_Pre", "Old_Post")
GROUP_LABS <- c("Young\nPre", "Young\nPost", "Old\nPre", "Old\nPost")

GROUP_FILL <- c(
  Young_Pre  = scales::alpha("#4393C3", 0.5),
  Young_Post = "#4393C3",
  Old_Pre    = scales::alpha("#D6604D", 0.5),
  Old_Post   = "#D6604D"
)

THEME_PUB <- theme_bw(base_size = 8) +
  theme(plot.title       = element_text(face = "bold", size = 9),
        plot.subtitle    = element_text(size = 6.5, color = "grey30", face = "italic"),
        strip.background = element_blank(),
        strip.text       = element_text(face = "bold", size = 6.5),
        legend.key.size  = unit(3, "mm"))

# --- Output dimensions (mm) ---
FIG_W <- 650
FIG_H <- 320
COL_WIDTHS <- c(0.12, 0.10, 0.45, 0.33)   # A, B, C, D

# === 6. HELPER FUNCTIONS ======================================================

clean_pathway_name <- function(name) {
  name |>
    str_remove("^HALLMARK_") |>
    str_remove("^GOBP_") |>
    str_remove("^GOCC_") |>
    str_remove("^GOMF_") |>
    str_replace_all("_", " ") |>
    str_to_title() |>
    str_replace("Mtorc1", "mTORC1") |>
    str_replace("Myc ", "MYC ") |>
    str_replace("E2f ", "E2F ") |>
    str_replace("Dna ", "DNA ") |>
    str_replace("Rna ", "RNA ") |>
    str_replace("Tnfa ", "TNFa ") |>
    str_replace("Uv ", "UV ") |>
    str_replace("G2m ", "G2M ") |>
    str_replace("Il6 ", "IL6 ") |>
    str_replace("Il2 ", "IL2 ") |>
    str_replace("Kras ", "KRAS ") |>
    str_replace("P53 ", "p53 ") |>
    str_replace("Tgf ", "TGF ") |>
    str_replace("Nf Kb", "NF-kB") |>
    str_replace("Atp ", "ATP ") |>
    str_replace("Nadh ", "NADH ") |>
    str_trunc(45, ellipsis = "...")
}

sig_stars <- function(padj) {
  dplyr::case_when(
    padj < 0.001 ~ "***",
    padj < 0.01  ~ "**",
    padj < 0.05  ~ "*",
    TRUE         ~ ""
  )
}

make_sigmoid_ribbon <- function(x0, x1, y0_top, y0_bot, y1_top, y1_bot,
                                n_pts = 50, ribbon_id) {
  t <- seq(0, 1, length.out = n_pts)
  blend <- (1 - cos(pi * t)) / 2
  tibble(
    x = c(x0 + (x1 - x0) * t, rev(x0 + (x1 - x0) * t)),
    y = c(y0_top + (y1_top - y0_top) * blend,
          rev(y0_bot + (y1_bot - y0_bot) * blend)),
    ribbon_id = ribbon_id
  )
}

reorder_within <- function(x, by, within, fun = mean, sep = "___", ...) {
  new_x <- paste(x, within, sep = sep)
  stats::reorder(new_x, by, FUN = fun, ...)
}

scale_y_reordered <- function(..., sep = "___") {
  reg <- paste0(sep, ".+$")
  ggplot2::scale_y_discrete(labels = function(x) gsub(reg, "", x), ...)
}

# === 7. LOAD IMPUTED MATRIX ===================================================

cat("Loading imputed matrix...\n")
imputed <- read_csv(file.path(BASE_DIR, "02_Imputation", "c_data", "01_imputed.csv"),
                    show_col_types = FALSE)
cat(sprintf("  Loaded: %d proteins x %d columns\n", nrow(imputed), ncol(imputed)))

# === 8. PARSE SAMPLE METADATA =================================================

# Identify annotation vs sample columns
annot_cols  <- c("uniprot_id", "protein", "gene", "description")
sample_cols <- setdiff(colnames(imputed), annot_cols)

# Build sample metadata
sample_meta <- tibble(sample = sample_cols) |>
  mutate(
    # Age group: O_ and OP_ prefixes = Old; Y_ and YP_ = Young
    age  = ifelse(str_detect(sample, "^(O_|OP_)"), "Old", "Young"),
    # Time point: _Pre or _Post suffix
    time = ifelse(str_detect(sample, "_Post$"), "Post", "Pre"),
    # Subject ID: strip _Pre/_Post suffix
    subject = str_remove(sample, "_(Pre|Post)$")
  )

cat(sprintf("  Samples: %d total (%d Young, %d Old)\n",
            nrow(sample_meta),
            sum(sample_meta$age == "Young"),
            sum(sample_meta$age == "Old")))

# === 9. COMPUTE DELTA MATRIX (Post - Pre) ====================================

cat("Computing delta matrix (Post - Pre)...\n")

# Identify Pre/Post pairs by subject
subjects <- unique(sample_meta$subject)
# Verify all subjects have both Pre and Post
pre_subjects  <- sample_meta$subject[sample_meta$time == "Pre"]
post_subjects <- sample_meta$subject[sample_meta$time == "Post"]
paired_subjects <- intersect(pre_subjects, post_subjects)
cat(sprintf("  Paired subjects: %d (of %d unique)\n",
            length(paired_subjects), length(subjects)))

# Sort subjects: Young first (sorted), then Old (sorted)
young_subjects <- sort(paired_subjects[str_detect(paired_subjects, "^(Y_|YP_)")])
old_subjects   <- sort(paired_subjects[str_detect(paired_subjects, "^(O_|OP_)")])
ordered_subjects <- c(young_subjects, old_subjects)

cat(sprintf("  Young subjects: %d | Old subjects: %d\n",
            length(young_subjects), length(old_subjects)))

# Prepare gene identifiers: use gene column, handle NA/duplicates
gene_ids <- imputed$gene
uniprot_ids <- imputed$uniprot_id

# Handle NA/empty gene names: use uniprot_id as fallback
na_mask <- is.na(gene_ids) | gene_ids == ""
if (any(na_mask)) {
  gene_ids[na_mask] <- uniprot_ids[na_mask]
  cat(sprintf("  Replaced %d missing gene names with uniprot IDs\n", sum(na_mask)))
}

# Handle duplicate gene names: append _2, _3, etc.
dup_counts <- table(gene_ids)
dup_genes  <- names(dup_counts[dup_counts > 1])
if (length(dup_genes) > 0) {
  cat(sprintf("  Resolving %d duplicate gene names\n", length(dup_genes)))
  for (dg in dup_genes) {
    idx <- which(gene_ids == dg)
    for (j in seq_along(idx)[-1]) {
      gene_ids[idx[j]] <- paste0(dg, "_", j)
    }
  }
}

# Build abundance matrix with gene symbols as row names
abund_mat <- imputed |>
  dplyr::select(all_of(sample_cols)) |>
  as.matrix()
rownames(abund_mat) <- gene_ids

# Compute delta for each paired subject
delta_mat <- matrix(NA_real_,
                    nrow = nrow(abund_mat),
                    ncol = length(ordered_subjects),
                    dimnames = list(gene_ids, ordered_subjects))

for (subj in ordered_subjects) {
  pre_col  <- paste0(subj, "_Pre")
  post_col <- paste0(subj, "_Post")
  delta_mat[, subj] <- abund_mat[, post_col] - abund_mat[, pre_col]
}

cat(sprintf("  Delta matrix: %d proteins x %d subjects\n",
            nrow(delta_mat), ncol(delta_mat)))

# === 10. Z-SCORE STANDARDIZE PER PROTEIN =====================================

cat("Z-score standardizing per protein (row-wise)...\n")
delta_z <- t(scale(t(delta_mat)))

# Check for any proteins with zero variance (constant delta across subjects)
zero_var <- apply(delta_mat, 1, sd, na.rm = TRUE) == 0
if (any(zero_var)) {
  cat(sprintf("  Removing %d zero-variance proteins\n", sum(zero_var)))
  delta_z <- delta_z[!zero_var, , drop = FALSE]
}

# Remove any rows with NaN (from zero-variance scaling)
nan_rows <- apply(delta_z, 1, function(x) any(is.nan(x)))
if (any(nan_rows)) {
  cat(sprintf("  Removing %d rows with NaN after scaling\n", sum(nan_rows)))
  delta_z <- delta_z[!nan_rows, , drop = FALSE]
}

cat(sprintf("  Z-scored matrix: %d proteins x %d subjects\n",
            nrow(delta_z), ncol(delta_z)))

# === 11. CREATE EXPRESSIONSET FOR MFUZZ ======================================

cat("Creating ExpressionSet for Mfuzz...\n")

# Mfuzz requires an ExpressionSet (rows = genes, cols = conditions/samples)
eset <- new("ExpressionSet",
            exprs = delta_z)

# Standardize for Mfuzz (already z-scored, but Mfuzz standardize() centers/scales)
eset <- standardise(eset)

cat(sprintf("  ExpressionSet: %d features x %d samples\n",
            nrow(exprs(eset)), ncol(exprs(eset))))

# === 12. MFUZZ CLUSTERING =====================================================

# Estimate fuzzifier
m_est <- mestimate(eset)
cat(sprintf("Estimated fuzzifier m = %.3f\n", m_est))

# Optimal k selection via Dmin
# Run for k=2 to 6, manually inspect
# Expected: 3-4 clusters based on prior analyses
set.seed(42)
# Use Dmin approach: compute Dmin for each k, pick elbow
dmin_vals <- numeric(5)
for (ki in 2:6) {
  cl_tmp <- mfuzz(eset, c = ki, m = m_est)
  dmin_vals[ki - 1] <- min(dist(cl_tmp$centers))
}
# Print Dmin values for inspection
cat("Dmin values (k=2-6):", paste(round(dmin_vals, 4), collapse = ", "), "\n")

# Select optimal k (use 4 as default, can be adjusted)
optimal_k <- 4
cat(sprintf("Using k = %d clusters\n", optimal_k))

# Final clustering
set.seed(42)
cl <- mfuzz(eset, c = optimal_k, m = m_est)
membership <- cl$membership
centroids  <- cl$centers

# Assign proteins to clusters (hard assignment = max membership)
cluster_assign <- tibble(
  gene = rownames(membership),
  cluster = paste0("C", max.col(membership)),
  membership = apply(membership, 1, max)
)

# Print cluster sizes
cat("Cluster sizes:\n")
print(table(cluster_assign$cluster))
cat(sprintf("Core proteins (membership >= 0.7): %d\n", sum(cluster_assign$membership >= 0.7)))

# === 13. SAVE CLUSTER ASSIGNMENTS AND EXPORT ==================================

write_csv(cluster_assign, file.path(DAT_DIR, "mfuzz_assignments.csv"))

centroid_export <- as.data.frame(centroids) |>
  rownames_to_column("cluster") |>
  mutate(cluster = paste0("C", row_number()))
write_csv(centroid_export, file.path(DAT_DIR, "mfuzz_centroids.csv"))

cat(sprintf("  Saved: %s\n", file.path(DAT_DIR, "mfuzz_assignments.csv")))
cat(sprintf("  Saved: %s\n", file.path(DAT_DIR, "mfuzz_centroids.csv")))

cat("=== Data pipeline complete: Setup + Mfuzz clustering ===\n")

# --- Core protein filter (membership >= 0.5) ---
CORE_THRESH <- 0.5
core_proteins <- cluster_assign %>%
  filter(membership >= CORE_THRESH)
cat(sprintf("Core proteins (membership >= %.1f): %d / %d (%.1f%%)\n",
            CORE_THRESH, nrow(core_proteins), nrow(cluster_assign),
            100 * nrow(core_proteins) / nrow(cluster_assign)))

# ╔══════════════════════════════════════════════════════════════════════════════╗
# ║                    PANEL BUILDING (REDESIGNED 2026-02-23)                   ║
# ║  Design: docs/plans/2026-02-23-figure4-redesign-design.md                  ║
# ║  Layout: A (Profiles) | B (PCA) | C (Sankey Triptych) | D (Cluster Synth)  ║
# ╚══════════════════════════════════════════════════════════════════════════════╝

# ============================================================================
# PANEL A — Cluster profiles with subject spaghetti
# ============================================================================

cat("=== Building Panel A: cluster profiles with subject spaghetti ===\n")

# --- A1. Group-level z-scores from raw abundance matrix ---
# For each of 4 groups (Young_Pre, Young_Post, Old_Pre, Old_Post), compute
# row means across matching samples, then z-score per protein (row-wise).
# Use only core_proteins$gene rows.

# Subset to core proteins
core_genes <- core_proteins$gene
core_abund <- abund_mat[core_genes, , drop = FALSE]

# Identify samples per group
group_samples <- list(
  Young_Pre  = sample_meta$sample[sample_meta$age == "Young" & sample_meta$time == "Pre"],
  Young_Post = sample_meta$sample[sample_meta$age == "Young" & sample_meta$time == "Post"],
  Old_Pre    = sample_meta$sample[sample_meta$age == "Old"   & sample_meta$time == "Pre"],
  Old_Post   = sample_meta$sample[sample_meta$age == "Old"   & sample_meta$time == "Post"]
)

cat(sprintf("  Group sample sizes: Young_Pre=%d, Young_Post=%d, Old_Pre=%d, Old_Post=%d\n",
            length(group_samples$Young_Pre), length(group_samples$Young_Post),
            length(group_samples$Old_Pre),   length(group_samples$Old_Post)))

# Compute group means (row means per group)
group_means <- sapply(group_samples, function(samps) {
  rowMeans(core_abund[, samps, drop = FALSE], na.rm = TRUE)
})
# group_means: proteins x 4 groups

# Z-score per protein (row-wise)
group_z <- t(scale(t(group_means)))
colnames(group_z) <- names(group_samples)

cat(sprintf("  Group z-score matrix: %d proteins x %d groups\n",
            nrow(group_z), ncol(group_z)))

# --- A2. Reshape to long format ---
panel_a_long <- as.data.frame(group_z) %>%
  rownames_to_column("gene") %>%
  pivot_longer(cols = all_of(GROUP_COLS),
               names_to = "group",
               values_to = "z_score") %>%
  left_join(core_proteins %>% dplyr::select(gene, cluster, membership), by = "gene") %>%
  mutate(
    age      = ifelse(str_detect(group, "^Young"), "Young", "Old"),
    time     = ifelse(str_detect(group, "_Post$"), "Post", "Pre"),
    time_num = ifelse(time == "Pre", 1, 2)
  )

cat(sprintf("  Long format: %d rows\n", nrow(panel_a_long)))

# --- A3. Per-cluster summary: mean_z and se_z ---
panel_a_summary <- panel_a_long %>%
  group_by(cluster, age, time, time_num) %>%
  summarise(
    mean_z = mean(z_score, na.rm = TRUE),
    se_z   = sd(z_score, na.rm = TRUE) / sqrt(n()),
    .groups = "drop"
  )

# --- A4. Shared y-axis range across all clusters ---
y_range <- panel_a_summary %>%
  summarise(
    y_lo = min(mean_z - 1.96 * se_z, na.rm = TRUE),
    y_hi = max(mean_z + 1.96 * se_z, na.rm = TRUE)
  )
y_pad   <- (y_range$y_hi - y_range$y_lo) * 0.1
y_limits <- c(y_range$y_lo - y_pad, y_range$y_hi + y_pad)

cat(sprintf("  Shared y-axis range: [%.2f, %.2f]\n", y_limits[1], y_limits[2]))

# --- A5. Build one plot per cluster ---
cluster_ids <- paste0("C", seq_len(optimal_k))
n_clusters  <- length(cluster_ids)

panels_A <- lapply(seq_along(cluster_ids), function(i) {
  cid <- cluster_ids[i]

  # Subset data for this cluster
  cl_data    <- panel_a_long %>% filter(cluster == cid)
  cl_summary <- panel_a_summary %>% filter(cluster == cid)

  n_total <- n_distinct(cl_data$gene)
  n_core  <- n_total  # already filtered to core_proteins

  # Determine if this is the first (top) or last (bottom) cluster
  is_first <- (i == 1)
  is_last  <- (i == n_clusters)

  # Build plot
  p <- ggplot() +
    # Subject spaghetti: ultra-thin lines
    geom_line(
      data = cl_data,
      aes(x = time_num, y = z_score, group = interaction(gene, age),
          colour = age, alpha = membership),
      linewidth = 0.15
    ) +
    # SE ribbon
    geom_ribbon(
      data = cl_summary,
      aes(x = time_num, ymin = mean_z - se_z, ymax = mean_z + se_z,
          fill = age, group = age),
      alpha = 0.15
    ) +
    # Centroid lines
    geom_line(
      data = cl_summary,
      aes(x = time_num, y = mean_z, colour = age, group = age),
      linewidth = 1.2
    ) +
    # Centroid points
    geom_point(
      data = cl_summary,
      aes(x = time_num, y = mean_z, colour = age),
      size = 2.5
    ) +
    # Scales
    scale_colour_manual(values = AGE_COLORS, guide = "none") +
    scale_fill_manual(values = AGE_COLORS, guide = "none") +
    scale_alpha_continuous(range = c(0.02, 0.15), guide = "none") +
    scale_x_continuous(
      breaks = c(1, 2),
      labels = if (is_last) c("Pre", "Post") else NULL,
      limits = c(0.7, 2.3),
      expand = expansion(0)
    ) +
    scale_y_continuous(
      limits = y_limits,
      name   = if (is_first) "Z-score (group means)" else NULL
    ) +
    # Labels
    labs(
      title    = paste0("Cluster ", i),
      subtitle = sprintf("(n = %d)", n_core),
      x        = NULL
    ) +
    # Theme
    THEME_PUB +
    theme(
      plot.title       = element_text(colour = CLUSTER_COLORS[cid],
                                      face = "bold", size = 8, hjust = 0.5),
      plot.subtitle    = element_text(colour = "grey30", face = "italic",
                                      size = 6.5, hjust = 0.5),
      panel.border     = element_rect(colour = CLUSTER_COLORS[cid],
                                      linewidth = 0.6, fill = NA),
      axis.title.y     = if (is_first) element_text(size = 7) else element_blank(),
      axis.text.y      = element_text(size = 6),
      axis.text.x      = if (is_last) element_text(size = 7) else element_blank(),
      axis.ticks.x     = if (is_last) element_line() else element_blank(),
      plot.margin      = margin(t = 2, r = 2, b = if (is_last) 4 else 1, l = 2)
    )

  p
})

cat("  Panel A: built", length(panels_A), "cluster profile plots\n")

# --- A6. Verification: save temporary PDF ---
pA_stack <- wrap_plots(panels_A, ncol = 1)
ggsave(file.path(RPT_DIR, "panel_A_profiles.pdf"), pA_stack,
       width = 80, height = 280, units = "mm")
cat(sprintf("  Saved verification PDF: %s\n", file.path(RPT_DIR, "panel_A_profiles.pdf")))

# ============================================================================
# PANEL B — PCA highlight-on-grey scatter
# ============================================================================

cat("=== Building Panel B: PCA highlight-on-grey scatter ===\n")

# --- B1. Run PCA on the delta matrix ---
# delta_mat: proteins (rows) x subjects (cols)
# t(delta_mat) => subjects x proteins; prcomp returns:
#   $rotation  = protein loadings (proteins x PCs) — one point per protein
pca_result <- prcomp(t(delta_mat), center = TRUE, scale. = FALSE)

# Protein scores for plotting (each dot = one protein)
pca_coords <- as.data.frame(pca_result$rotation[, 1:2])
pca_coords$gene <- rownames(pca_coords)

# Variance explained
var_explained <- summary(pca_result)$importance[2, 1:2] * 100
pc1_lab <- sprintf("PC1 (%.1f%%)", var_explained[1])
pc2_lab <- sprintf("PC2 (%.1f%%)", var_explained[2])

cat(sprintf("  PCA: %d proteins projected, PC1=%.1f%%, PC2=%.1f%%\n",
            nrow(pca_coords), var_explained[1], var_explained[2]))

# --- B2. Create pca_scores data frame ---
pca_scores <- pca_coords %>%
  left_join(core_proteins %>% dplyr::select(gene, cluster, membership),
            by = "gene") %>%
  mutate(
    cluster    = ifelse(is.na(cluster), NA_character_, cluster),
    membership = ifelse(is.na(membership), 0, membership)
  )

cat(sprintf("  pca_scores: %d rows (%d with cluster assignment)\n",
            nrow(pca_scores), sum(!is.na(pca_scores$cluster))))

# --- B3. Export PCA scores ---
write_csv(pca_scores, file.path(DAT_DIR, "fig4_panel_B_pca.csv"))
cat(sprintf("  Saved: %s\n", file.path(DAT_DIR, "fig4_panel_B_pca.csv")))

# --- B4. Build one highlight-on-grey plot per cluster ---
panels_B <- lapply(seq_along(cluster_ids), function(i) {
  cid <- cluster_ids[i]

  is_first <- (i == 1)
  is_last  <- (i == n_clusters)

  # Core proteins for this cluster
  highlight_data <- pca_scores %>%
    filter(cluster == cid)

  p <- ggplot() +
    # Background: ALL proteins as grey points
    geom_point(
      data = pca_scores,
      aes(x = PC1, y = PC2),
      color = "grey80", size = 0.3, alpha = 0.3
    ) +
    # Highlighted cluster: core proteins in cluster color
    geom_point(
      data = highlight_data,
      aes(x = PC1, y = PC2, alpha = membership),
      color = CLUSTER_COLORS[cid], size = 0.6
    ) +
    scale_alpha_continuous(range = c(0.3, 1.0), guide = "none") +
    # Cluster label: top-left corner
    annotate("text",
             x = min(pca_scores$PC1), y = max(pca_scores$PC2),
             label = cid,
             color = CLUSTER_COLORS[cid], fontface = "bold", size = 2.8,
             hjust = 0, vjust = 1) +
    # Axis labels: only on bottom row
    labs(
      x = if (is_last) pc1_lab else NULL,
      y = if (is_first) pc2_lab else NULL
    ) +
    THEME_PUB +
    theme(
      panel.border = element_rect(colour = CLUSTER_COLORS[cid],
                                  linewidth = 0.6, fill = NA),
      axis.title.x = if (is_last) element_text(size = 6) else element_blank(),
      axis.title.y = if (is_first) element_text(size = 6) else element_blank(),
      axis.text.x  = if (is_last) element_text(size = 5) else element_blank(),
      axis.text.y  = element_text(size = 5),
      axis.ticks.x = if (is_last) element_line() else element_blank(),
      plot.margin  = margin(t = 2, r = 2, b = if (is_last) 4 else 1, l = 2)
    )

  p
})

cat("  Panel B: built", length(panels_B), "PCA highlight-on-grey plots\n")

# --- B5. Verification: save temporary PDF ---
pB_stack <- wrap_plots(panels_B, ncol = 1)
ggsave(file.path(RPT_DIR, "panel_B_pca.pdf"), pB_stack,
       width = 65, height = 280, units = "mm")
cat(sprintf("  Saved verification PDF: %s\n", file.path(RPT_DIR, "panel_B_pca.pdf")))

# === ENRICHMENT ANALYSIS (shared by Panels C and D) ==========================

cat("=== Enrichment analysis: ORA + rrvgo + 1:1 greedy assignment ===\n")

# --- Step 1: Gene set loading ------------------------------------------------

hallmark_t2g <- msigdbr(species = "Homo sapiens", category = "H") %>%
  dplyr::select(gs_name, gene_symbol) %>%
  dplyr::rename(term = gs_name, gene = gene_symbol)

gobp_full <- msigdbr(species = "Homo sapiens", category = "C5",
                      subcategory = "GO:BP")
gobp_t2g <- gobp_full %>%
  dplyr::select(gs_name, gene_symbol) %>%
  dplyr::rename(term = gs_name, gene = gene_symbol)
# Build name -> GO ID map for rrvgo
gobp_id_map <- gobp_full %>%
  dplyr::select(gs_name, gs_exact_source) %>%
  dplyr::distinct()

gocc_full <- msigdbr(species = "Homo sapiens", category = "C5",
                      subcategory = "GO:CC")
gocc_t2g <- gocc_full %>%
  dplyr::select(gs_name, gene_symbol) %>%
  dplyr::rename(term = gs_name, gene = gene_symbol)
gocc_id_map <- gocc_full %>%
  dplyr::select(gs_name, gs_exact_source) %>%
  dplyr::distinct()

universe_genes <- unique(cluster_assign$gene)

cat(sprintf("  Gene sets loaded: Hallmark=%d terms, GO:BP=%d terms, GO:CC=%d terms\n",
            n_distinct(hallmark_t2g$term),
            n_distinct(gobp_t2g$term),
            n_distinct(gocc_t2g$term)))
cat(sprintf("  Universe: %d genes\n", length(universe_genes)))

# --- Step 2: Per-cluster ORA -------------------------------------------------

cat("Running per-cluster ORA...\n")

enrich_list <- list()

for (cl_id in seq_len(optimal_k)) {
  cl_label <- paste0("C", cl_id)
  cl_genes <- core_proteins$gene[core_proteins$cluster == cl_label]
  cat(sprintf("  %s: %d core genes\n", cl_label, length(cl_genes)))

  # Run enricher for each database
  res_hall <- enricher(cl_genes, TERM2GENE = hallmark_t2g,
                       universe = universe_genes,
                       pAdjustMethod = "BH",
                       pvalueCutoff = 0.05, qvalueCutoff = 1)

  res_gobp <- enricher(cl_genes, TERM2GENE = gobp_t2g,
                        universe = universe_genes,
                        pAdjustMethod = "BH",
                        pvalueCutoff = 0.05, qvalueCutoff = 1)

  res_gocc <- enricher(cl_genes, TERM2GENE = gocc_t2g,
                        universe = universe_genes,
                        pAdjustMethod = "BH",
                        pvalueCutoff = 0.05, qvalueCutoff = 1)

  # Combine results with database column
  combined <- bind_rows(
    if (!is.null(res_hall) && nrow(as.data.frame(res_hall)) > 0)
      as.data.frame(res_hall) %>% mutate(database = "Hallmark") else NULL,
    if (!is.null(res_gobp) && nrow(as.data.frame(res_gobp)) > 0)
      as.data.frame(res_gobp) %>% mutate(database = "GO:BP") else NULL,
    if (!is.null(res_gocc) && nrow(as.data.frame(res_gocc)) > 0)
      as.data.frame(res_gocc) %>% mutate(database = "GO:CC") else NULL
  )

  # Filter to p.adjust < 0.05
  if (nrow(combined) > 0) {
    combined <- combined %>% filter(p.adjust < 0.05)
  }

  enrich_list[[cl_label]] <- combined %>% mutate(cluster = cl_label)
  cat(sprintf("    %s: %d significant terms (H=%d, BP=%d, CC=%d)\n",
              cl_label, nrow(combined),
              sum(combined$database == "Hallmark"),
              sum(combined$database == "GO:BP"),
              sum(combined$database == "GO:CC")))
}

# --- Step 3: rrvgo reduction -------------------------------------------------

cat("Applying rrvgo redundancy reduction for GO terms...\n")

# Pre-compute semantic data objects
bp_semdata <- tryCatch(
  godata("org.Hs.eg.db", ont = "BP", computeIC = TRUE),
  error = function(e) { cat("  Warning: could not compute BP semdata\n"); NULL }
)
cc_semdata <- tryCatch(
  godata("org.Hs.eg.db", ont = "CC", computeIC = TRUE),
  error = function(e) { cat("  Warning: could not compute CC semdata\n"); NULL }
)

reduce_go_terms <- function(enrich_df, ont, id_map, semdata) {
  # If no terms or no semdata, return as-is
  if (is.null(semdata) || nrow(enrich_df) == 0) return(enrich_df)

  db_label <- paste0("GO:", ont)
  go_terms <- enrich_df %>% filter(database == db_label)
  other_terms <- enrich_df %>% filter(database != db_label)

  if (nrow(go_terms) < 2) return(enrich_df)

  # Map MSigDB names to GO IDs
  go_terms <- go_terms %>%
    left_join(id_map, by = c("ID" = "gs_name")) %>%
    dplyr::rename(go_id = gs_exact_source)

  # Remove terms without GO ID mapping
  go_terms <- go_terms %>% filter(!is.na(go_id))
  if (nrow(go_terms) < 2) {
    go_terms <- go_terms %>% dplyr::select(-go_id)
    return(bind_rows(other_terms, go_terms))
  }

  # Build named p-value vector (GO IDs as names)
  scores <- setNames(go_terms$p.adjust, go_terms$go_id)

  reduced <- tryCatch({
    sim_mat <- calculateSimMatrix(go_terms$go_id,
                                   orgdb = "org.Hs.eg.db",
                                   ont = ont,
                                   method = "Rel",
                                   semdata = semdata)
    red <- reduceSimMatrix(sim_mat, scores = scores, threshold = 0.85,
                               orgdb = "org.Hs.eg.db")
    # Keep only parent terms (use 'parent' column which contains GO IDs)
    parent_go_ids <- unique(red$parent)
    go_terms %>% filter(go_id %in% parent_go_ids) %>% dplyr::select(-go_id)
  }, error = function(e) {
    cat(sprintf("    rrvgo warning (%s): %s — keeping all terms\n", ont, e$message))
    go_terms %>% dplyr::select(-go_id)
  })

  bind_rows(other_terms, reduced)
}

# Apply rrvgo to each cluster's results
for (cl_label in names(enrich_list)) {
  before_n <- nrow(enrich_list[[cl_label]])
  enrich_list[[cl_label]] <- reduce_go_terms(enrich_list[[cl_label]], "BP",
                                              gobp_id_map, bp_semdata)
  enrich_list[[cl_label]] <- reduce_go_terms(enrich_list[[cl_label]], "CC",
                                              gocc_id_map, cc_semdata)
  after_n <- nrow(enrich_list[[cl_label]])
  cat(sprintf("  %s: %d -> %d terms after rrvgo\n", cl_label, before_n, after_n))
}

# --- Step 4: Top term selection + 1:1 greedy assignment -----------------------

cat("Selecting top terms and performing 1:1 greedy assignment...\n")

# Select top 3 per database per cluster
enrich_top <- bind_rows(enrich_list) %>%
  group_by(cluster, database) %>%
  slice_min(p.adjust, n = 3, with_ties = FALSE) %>%
  ungroup()

cat(sprintf("  Top terms: %d total across %d clusters\n",
            nrow(enrich_top), n_distinct(enrich_top$cluster)))

# 1:1 greedy assignment: for each cluster, assign each core protein to its
# best pathway (lowest p.adjust that contains the gene)
protein_pathway_links <- bind_rows(lapply(seq_len(optimal_k), function(cl_id) {
  cl_label <- paste0("C", cl_id)
  cl_genes <- core_proteins$gene[core_proteins$cluster == cl_label]
  cl_top   <- enrich_top %>% filter(cluster == cl_label) %>% arrange(p.adjust)

  if (nrow(cl_top) == 0) {
    return(tibble(gene = cl_genes, pathway = "Unmapped",
                  database = NA_character_, cluster = cl_label))
  }

  # Parse geneID column (slash-separated) into a list
  cl_top$gene_list <- strsplit(cl_top$geneID, "/")

  # Greedy assignment: iterate through genes, assign to best pathway containing it
  assigned <- tibble(gene = character(), pathway = character(),
                     database = character(), cluster = character())

  for (g in cl_genes) {
    best_idx <- NA
    for (j in seq_len(nrow(cl_top))) {
      if (g %in% cl_top$gene_list[[j]]) {
        best_idx <- j
        break  # Already sorted by p.adjust, so first match is best
      }
    }
    if (!is.na(best_idx)) {
      assigned <- bind_rows(assigned, tibble(
        gene     = g,
        pathway  = cl_top$Description[best_idx],
        database = cl_top$database[best_idx],
        cluster  = cl_label
      ))
    } else {
      assigned <- bind_rows(assigned, tibble(
        gene     = g,
        pathway  = "Unmapped",
        database = NA_character_,
        cluster  = cl_label
      ))
    }
  }

  assigned
}))

# Print summary
for (cl_label in cluster_ids) {
  cl_links <- protein_pathway_links %>% filter(cluster == cl_label)
  n_mapped   <- sum(cl_links$pathway != "Unmapped")
  n_unmapped <- sum(cl_links$pathway == "Unmapped")
  n_pathways <- n_distinct(cl_links$pathway[cl_links$pathway != "Unmapped"])
  cat(sprintf("  Cluster %s: %d pathways, %d proteins mapped, %d unmapped\n",
              cl_label, n_pathways, n_mapped, n_unmapped))
}

# --- Step 5: Export -----------------------------------------------------------

write_csv(enrich_top, file.path(DAT_DIR, "fig4_panel_C_enrichment.csv"))
write_csv(protein_pathway_links, file.path(DAT_DIR, "fig4_panel_C_sankey_links.csv"))

cat(sprintf("  Saved: %s\n", file.path(DAT_DIR, "fig4_panel_C_enrichment.csv")))
cat(sprintf("  Saved: %s\n", file.path(DAT_DIR, "fig4_panel_C_sankey_links.csv")))

# --- Step 6: Top Hallmark term per cluster (for Panel A subtitles later) ------

top_hallmark <- enrich_top %>%
  filter(database == "Hallmark") %>%
  group_by(cluster) %>%
  slice_min(p.adjust, n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  mutate(label = clean_pathway_name(Description))

cat("Top Hallmark terms per cluster:\n")
for (i in seq_len(nrow(top_hallmark))) {
  cat(sprintf("  %s: %s (p.adj = %.2e)\n",
              top_hallmark$cluster[i],
              top_hallmark$label[i],
              top_hallmark$p.adjust[i]))
}

cat("=== Enrichment analysis complete ===\n")
