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

# ============================================================================
# PANEL C — Per-cluster Sankey triptych (heatmap | Sankey | enrichment bars)
# ============================================================================

cat("=== Building Panel C: per-cluster Sankey triptych ===\n")

Z_CAP <- 2

# --- Shared x-axis max for enrichment bars across all clusters ---
shared_x_max_C <- enrich_top %>%
  mutate(neg_log10_p = -log10(p.adjust)) %>%
  pull(neg_log10_p) %>%
  max(na.rm = TRUE) * 1.15

cat(sprintf("  Shared x-axis max (enrichment bars): %.2f\n", shared_x_max_C))

# --- build_panel_C function ---
build_panel_C <- function(cl_id, show_xlab, shared_x_max) {

  cl_label <- paste0("C", cl_id)
  cl_col   <- CLUSTER_COLORS[cl_label]

  # -- Genes for this cluster, sorted by descending membership --
  cl_core <- core_proteins %>%
    filter(cluster == cl_label) %>%
    arrange(desc(membership))
  cl_genes <- cl_core$gene

  n_genes <- length(cl_genes)
  cat(sprintf("  Panel C — %s: %d core genes\n", cl_label, n_genes))

  # -- Links and enrichment for this cluster --
  cl_links <- protein_pathway_links %>%
    filter(cluster == cl_label)
  cl_enrich <- enrich_top %>%
    filter(cluster == cl_label) %>%
    arrange(p.adjust)

  # Mapped pathways (exclude Unmapped)
  mapped_links <- cl_links %>% filter(pathway != "Unmapped")
  mapped_pathways <- unique(mapped_links$pathway)
  n_pw <- length(mapped_pathways)

  cat(sprintf("    Mapped pathways: %d, Mapped proteins: %d, Unmapped: %d\n",
              n_pw, nrow(mapped_links), sum(cl_links$pathway == "Unmapped")))

  # ===========================================================================
  # Sub-panel C1: Heatmap Strip
  # ===========================================================================

  # Get z-scores for this cluster's genes from group_z
  ht_genes <- intersect(cl_genes, rownames(group_z))
  ht_mat   <- group_z[ht_genes, , drop = FALSE]

  # Cap at +/- Z_CAP
  ht_mat[ht_mat >  Z_CAP] <-  Z_CAP
  ht_mat[ht_mat < -Z_CAP] <- -Z_CAP

  # Reshape to long
  ht_long <- as.data.frame(ht_mat) %>%
    rownames_to_column("gene") %>%
    pivot_longer(cols = all_of(GROUP_COLS),
                 names_to = "group", values_to = "z") %>%
    mutate(
      gene  = factor(gene, levels = rev(ht_genes)),
      group = factor(group, levels = GROUP_COLS)
    )

  # Diverging color palette
  ht_colors <- c("#2166AC", "#92C5DE", "white", "#F4A582", "#B2182B")

  p_ht <- ggplot(ht_long, aes(x = group, y = gene, fill = z)) +
    geom_raster() +
    scale_fill_gradientn(
      colours = ht_colors,
      limits  = c(-Z_CAP, Z_CAP),
      oob     = scales::squish,
      name    = "Z",
      guide   = guide_colorbar(barwidth = unit(2, "mm"),
                                barheight = unit(12, "mm"),
                                title.position = "top", title.hjust = 0.5)
    ) +
    scale_x_discrete(
      labels = if (show_xlab) GROUP_LABS else NULL,
      expand = expansion(0)
    ) +
    scale_y_discrete(expand = expansion(0)) +
    labs(x = NULL, y = NULL) +
    THEME_PUB +
    theme(
      axis.text.y     = element_blank(),
      axis.ticks.y    = element_blank(),
      axis.text.x     = if (show_xlab) element_text(size = 5, lineheight = 0.85)
                         else element_blank(),
      axis.ticks.x    = if (show_xlab) element_line() else element_blank(),
      panel.border    = element_rect(colour = cl_col, linewidth = 0.6, fill = NA),
      legend.position = "none",
      plot.margin     = margin(t = 2, r = 1, b = if (show_xlab) 4 else 1, l = 2)
    )

  # ===========================================================================
  # Sub-panel C2: Sigmoid Sankey
  # ===========================================================================

  if (n_pw == 0) {
    # No mapped pathways — return placeholder
    p_sankey <- ggplot() +
      annotate("text", x = 0.5, y = 0.5, label = "No mapped\npathways",
               size = 2.5, color = "grey50") +
      theme_void() +
      theme(plot.margin = margin(t = 2, r = 1, b = if (show_xlab) 4 else 1, l = 1))
  } else {

    Y_SPAN  <- n_genes
    X_GENE  <- 0.5
    X_PW    <- 2.5
    BAR_W   <- 0.06

    # -- Gene bar positions (top to bottom) --
    gene_h   <- Y_SPAN / (n_genes * 1.15)
    gene_gap <- if (n_genes > 1) (Y_SPAN - n_genes * gene_h) / (n_genes - 1) else 0

    gene_bars <- tibble(
      gene  = cl_genes,
      idx   = seq_along(cl_genes),
      y_top = Y_SPAN - (idx - 1) * (gene_h + gene_gap),
      y_bot = y_top - gene_h,
      y_ctr = (y_top + y_bot) / 2,
      x_left  = X_GENE - BAR_W / 2,
      x_right = X_GENE + BAR_W / 2,
      fill  = cl_col
    )

    # -- Pathway bar positions (top to bottom) --
    pw_h   <- Y_SPAN / (n_pw * 1.4)
    pw_gap <- if (n_pw > 1) (Y_SPAN - n_pw * pw_h) / (n_pw - 1) else 0

    # Order pathways to match enrichment order (lowest p.adjust first at top)
    pw_order <- cl_enrich %>%
      filter(Description %in% mapped_pathways) %>%
      pull(Description)
    # Add any mapped_pathways not in cl_enrich at the end
    pw_order <- c(pw_order, setdiff(mapped_pathways, pw_order))

    pw_bars <- tibble(
      pathway = pw_order,
      idx     = seq_along(pw_order),
      y_top   = Y_SPAN - (idx - 1) * (pw_h + pw_gap),
      y_bot   = y_top - pw_h,
      y_ctr   = (y_top + y_bot) / 2,
      x_left  = X_PW - BAR_W / 2,
      x_right = X_PW + BAR_W / 2
    ) %>%
      left_join(
        mapped_links %>% dplyr::select(pathway, database) %>% distinct(),
        by = "pathway"
      ) %>%
      mutate(fill = DB_COLORS[database])

    # -- Slot allocation within each pathway bar --
    # Each pathway bar is subdivided into slots (one per gene mapping to it)
    slot_data <- mapped_links %>%
      group_by(pathway) %>%
      mutate(slot_idx = row_number(), n_slots = n()) %>%
      ungroup() %>%
      left_join(pw_bars %>% dplyr::select(pathway, y_top, y_bot), by = "pathway") %>%
      mutate(
        slot_h   = (y_top - y_bot) / n_slots,
        slot_top = y_top - (slot_idx - 1) * slot_h,
        slot_bot = slot_top - slot_h,
        slot_ctr = (slot_top + slot_bot) / 2
      )

    # -- Build ribbon data --
    ribbon_input <- slot_data %>%
      left_join(gene_bars %>% dplyr::select(gene, y_top_gene = y_top,
                                             y_bot_gene = y_bot, y_ctr_gene = y_ctr),
                by = "gene") %>%
      mutate(
        ribbon_id = paste0(cl_label, "_", gene, "_", pathway),
        cross_dist = abs(y_ctr_gene - slot_ctr)
      ) %>%
      arrange(desc(cross_dist))  # far-crossing ribbons drawn first (behind)

    ribbon_polys <- pmap_dfr(ribbon_input, function(gene, pathway, database, cluster,
                                                     slot_idx, n_slots,
                                                     y_top, y_bot,
                                                     slot_h, slot_top, slot_bot, slot_ctr,
                                                     y_top_gene, y_bot_gene, y_ctr_gene,
                                                     ribbon_id, cross_dist) {
      make_sigmoid_ribbon(
        x0 = X_GENE + BAR_W / 2,
        x1 = X_PW   - BAR_W / 2,
        y0_top = y_top_gene,
        y0_bot = y_bot_gene,
        y1_top = slot_top,
        y1_bot = slot_bot,
        ribbon_id = ribbon_id
      ) %>%
        mutate(fill = DB_COLORS[database])
    })

    # -- Build the plot --
    # Gene bar rectangles as data frame for geom_rect
    gene_rect <- gene_bars %>%
      dplyr::select(x_left, x_right, y_bot, y_top, fill)
    pw_rect <- pw_bars %>%
      dplyr::select(x_left, x_right, y_bot, y_top, fill)

    p_sankey <- ggplot() +
      # Ribbons (drawn first = behind bars)
      geom_polygon(
        data = ribbon_polys,
        aes(x = x, y = y, group = ribbon_id, fill = fill),
        alpha = 0.25, color = NA
      ) +
      # Gene bars
      geom_rect(
        data = gene_rect,
        aes(xmin = x_left, xmax = x_right, ymin = y_bot, ymax = y_top, fill = fill),
        color = NA
      ) +
      # Pathway bars
      geom_rect(
        data = pw_rect,
        aes(xmin = x_left, xmax = x_right, ymin = y_bot, ymax = y_top, fill = fill),
        color = NA
      ) +
      scale_fill_identity() +
      coord_cartesian(xlim = c(0, 3), ylim = c(0, Y_SPAN), expand = FALSE) +
      theme_void() +
      theme(
        plot.margin = margin(t = 2, r = 1, b = if (show_xlab) 4 else 1, l = 1)
      )
  }

  # ===========================================================================
  # Sub-panel C3: Enrichment Bars
  # ===========================================================================

  if (nrow(cl_enrich) == 0) {
    p_bars <- ggplot() +
      annotate("text", x = 0.5, y = 0.5, label = "No significant\nenrichment",
               size = 2.5, color = "grey50") +
      theme_void() +
      theme(plot.margin = margin(t = 2, r = 2, b = if (show_xlab) 4 else 1, l = 1))
  } else {

    bar_data <- cl_enrich %>%
      mutate(
        neg_log10_p = -log10(p.adjust),
        clean_name  = clean_pathway_name(Description),
        gene_count  = as.numeric(sub("/.*", "", GeneRatio)),
        db_fill     = DB_COLORS[database]
      )

    # Order pathways by -log10(p.adjust) within cluster
    bar_data <- bar_data %>%
      mutate(clean_name = reorder_within(clean_name, neg_log10_p, cluster))

    p_bars <- ggplot(bar_data, aes(x = neg_log10_p, y = clean_name)) +
      geom_col(aes(fill = db_fill), color = "black", linewidth = 0.3) +
      geom_text(aes(label = gene_count),
                hjust = -0.2, size = KEY_TEXT, fontface = "bold") +
      geom_vline(xintercept = -log10(0.05), linetype = "dashed",
                 color = "grey40", linewidth = 0.3) +
      scale_fill_identity() +
      scale_x_continuous(
        limits = c(0, shared_x_max),
        expand = expansion(mult = c(0, 0.05)),
        name   = if (show_xlab) expression(-log[10](p[adj])) else NULL
      ) +
      scale_y_reordered() +
      labs(y = NULL) +
      THEME_PUB +
      theme(
        axis.text.y  = element_text(size = 5),
        axis.text.x  = if (show_xlab) element_text(size = 5) else element_blank(),
        axis.ticks.x = if (show_xlab) element_line() else element_blank(),
        axis.title.x = if (show_xlab) element_text(size = 6) else element_blank(),
        panel.border = element_rect(colour = "grey70", linewidth = 0.3, fill = NA),
        plot.margin  = margin(t = 2, r = 2, b = if (show_xlab) 4 else 1, l = 1)
      )
  }

  # ===========================================================================
  # Assemble the triptych
  # ===========================================================================

  (p_ht | p_sankey | p_bars) + plot_layout(widths = c(0.15, 0.45, 0.40))
}

# --- Build all clusters ---
panels_C <- lapply(1:optimal_k, function(i) {
  build_panel_C(i, show_xlab = (i == optimal_k), shared_x_max = shared_x_max_C)
})

cat("  Panel C: built", length(panels_C), "triptych rows\n")

# --- Verification: save temporary PDF ---
pC_stack <- wrap_plots(panels_C, ncol = 1)
ggsave(file.path(RPT_DIR, "panel_C_triptych.pdf"), pC_stack,
       width = 280, height = 300, units = "mm")
cat(sprintf("  Saved verification PDF: %s\n", file.path(RPT_DIR, "panel_C_triptych.pdf")))

# ============================================================================
# PANEL D — Cluster synthesis Sankey
# ============================================================================

cat("=== Building Panel D: Cluster synthesis Sankey ===\n")

# === PANEL D: CLUSTER SYNTHESIS SANKEY ========================================

# --- Step 1: Theme curation ---------------------------------------------------

assign_theme <- function(pathway_name) {
  pw <- tolower(pathway_name)
  dplyr::case_when(
    stringr::str_detect(pw, "mitochon|oxidative|respiratory|electron|tca|citrate|nadh|atp") ~
      "Mitochondrial Metabolism",
    stringr::str_detect(pw, "myogen|muscle|contract|actin|myofib|sarco") ~
      "Muscle Structure",
    stringr::str_detect(pw, "riboso|translat|proteasom|ubiquitin|protein fold") ~
      "Protein Homeostasis",
    stringr::str_detect(pw, "immun|inflam|complement|cytokine|interferon") ~
      "Immune & Inflammatory",
    stringr::str_detect(pw, "vesicle|transport|endosom|golgi|lysosom") ~
      "Intracellular Transport",
    stringr::str_detect(pw, "mtorc|signal|kinase|cell cycle|proliferat") ~
      "Cell Signaling",
    stringr::str_detect(pw, "extracellular|matrix|collagen|adhesion|integrin") ~
      "ECM & Adhesion",
    stringr::str_detect(pw, "glycol|lipid|fatty|metabol|xenobiot|bile") ~
      "Metabolic Regulation",
    TRUE ~ "Other"
  )
}

theme_links <- protein_pathway_links %>%
  filter(pathway != "Unmapped") %>%
  mutate(theme = assign_theme(pathway))

cat(sprintf("  Theme assignment: %d mapped proteins across themes\n", nrow(theme_links)))
cat("  Theme distribution:\n")
print(table(theme_links$theme))

# Compute summaries
cluster_theme_counts <- theme_links %>%
  count(cluster, theme, name = "n_proteins") %>%
  arrange(cluster, desc(n_proteins))

theme_totals <- cluster_theme_counts %>%
  group_by(theme) %>%
  summarise(total = sum(n_proteins), .groups = "drop") %>%
  arrange(desc(total))

theme_order <- theme_totals$theme

cat("  Theme totals (ordered):\n")
for (i in seq_len(nrow(theme_totals))) {
  cat(sprintf("    %s: %d proteins\n", theme_totals$theme[i], theme_totals$total[i]))
}

# Export
write_csv(cluster_theme_counts, file.path(DAT_DIR, "fig4_panel_D_cluster_themes.csv"))
cat(sprintf("  Saved: %s\n", file.path(DAT_DIR, "fig4_panel_D_cluster_themes.csv")))

# --- Step 2: Cluster-to-theme Sankey coordinate system ------------------------

# Exclude "Other" from the Sankey diagram — it adds no biological insight
# and compresses the meaningful themes.  Keep it in the exported CSV.
sankey_theme_counts <- cluster_theme_counts %>% filter(theme != "Other")
sankey_theme_totals <- sankey_theme_counts %>%
  group_by(theme) %>%
  summarise(total = sum(n_proteins), .groups = "drop") %>%
  arrange(desc(total))
sankey_theme_order <- sankey_theme_totals$theme

cat(sprintf("  Sankey themes (excl. Other): %d themes, %d proteins\n",
            length(sankey_theme_order), sum(sankey_theme_totals$total)))

D_Y_SPAN <- 100
D_X_CL   <- 1.0
D_X_TH   <- 3.0
D_BAR_W  <- 0.20   # wider bars so cluster labels are legible

# Cluster mapped protein counts (only non-Other themed proteins)
cluster_mapped <- sankey_theme_counts %>%
  group_by(cluster) %>%
  summarise(n_mapped = sum(n_proteins), .groups = "drop") %>%
  arrange(cluster)

# Ensure all active clusters are present
active_clusters <- sort(unique(sankey_theme_counts$cluster))
n_active <- length(active_clusters)

cat(sprintf("  Active clusters with mapped proteins: %d\n", n_active))

# --- Cluster bars (left side) ---
# Stack from top to bottom, sized proportional to mapped protein count
cl_total_mapped <- sum(cluster_mapped$n_mapped)
cl_gap_frac <- 0.05  # 5% of Y_SPAN for gaps between bars
cl_usable   <- D_Y_SPAN * (1 - cl_gap_frac * (n_active - 1) / n_active)
cl_gap_size <- if (n_active > 1) (D_Y_SPAN - cl_usable) / (n_active - 1) else 0

cl_bars <- cluster_mapped %>%
  mutate(
    bar_h   = (n_mapped / cl_total_mapped) * cl_usable,
    # Stack from top
    y_top   = D_Y_SPAN - cumsum(c(0, head(bar_h + cl_gap_size, -1))),
    y_bot   = y_top - bar_h,
    y_ctr   = (y_top + y_bot) / 2,
    x_left  = D_X_CL - D_BAR_W / 2,
    x_right = D_X_CL + D_BAR_W / 2,
    fill    = CLUSTER_COLORS[cluster]
  )

# --- Theme bars (right side) ---
# One bar per theme (excl. Other), ordered by total protein count (descending)
n_themes <- length(sankey_theme_order)
th_gap_frac <- 0.03
th_usable   <- D_Y_SPAN * (1 - th_gap_frac * (n_themes - 1) / n_themes)
th_gap_size <- if (n_themes > 1) (D_Y_SPAN - th_usable) / (n_themes - 1) else 0
th_total    <- sum(sankey_theme_totals$total)

th_bars <- sankey_theme_totals %>%
  mutate(
    theme = factor(theme, levels = sankey_theme_order),
    bar_h   = (total / th_total) * th_usable,
    y_top   = D_Y_SPAN - cumsum(c(0, head(bar_h + th_gap_size, -1))),
    y_bot   = y_top - bar_h,
    y_ctr   = (y_top + y_bot) / 2,
    x_left  = D_X_TH - D_BAR_W / 2,
    x_right = D_X_TH + D_BAR_W / 2
  ) %>%
  arrange(theme)

# --- Ribbon construction: cluster side tracking ---
# For each cluster, track cumulative allocation of its bar space across themes
# For each theme, track cumulative allocation of its bar space across clusters

# Initialize cumulative trackers
cl_cum <- setNames(rep(0, n_active), active_clusters)   # fraction used per cluster
th_cum <- setNames(rep(0, n_themes), sankey_theme_order) # fraction used per theme

# Build ribbons by iterating themes (top to bottom) then clusters within each theme
ribbon_list <- list()
ribbon_idx  <- 0

for (th in sankey_theme_order) {
  # Get cluster contributions to this theme
  th_contribs <- sankey_theme_counts %>%
    filter(theme == th) %>%
    arrange(cluster)

  for (r in seq_len(nrow(th_contribs))) {
    cl <- th_contribs$cluster[r]
    n  <- th_contribs$n_proteins[r]
    if (n == 0) next

    ribbon_idx <- ribbon_idx + 1

    # --- Cluster bar side: slice proportional to this cluster's mapped count ---
    cl_row   <- cl_bars %>% filter(cluster == cl)
    cl_n     <- cl_row$n_mapped
    cl_h     <- cl_row$bar_h
    cl_y_top <- cl_row$y_top

    # This ribbon occupies fraction n/cl_n of the cluster bar
    frac_cl  <- n / cl_n
    # Start from cumulative position (top of bar going down)
    y0_top <- cl_y_top - cl_cum[cl] * cl_h
    y0_bot <- y0_top - frac_cl * cl_h
    cl_cum[cl] <- cl_cum[cl] + frac_cl

    # --- Theme bar side: slice proportional to theme total ---
    th_row   <- th_bars %>% filter(theme == th)
    th_n     <- th_row$total
    th_h     <- th_row$bar_h
    th_y_top <- th_row$y_top

    frac_th  <- n / th_n
    y1_top <- th_y_top - th_cum[th] * th_h
    y1_bot <- y1_top - frac_th * th_h
    th_cum[th] <- th_cum[th] + frac_th

    # Build sigmoid ribbon polygon
    ribbon_poly <- make_sigmoid_ribbon(
      x0 = D_X_CL + D_BAR_W / 2,
      x1 = D_X_TH - D_BAR_W / 2,
      y0_top = y0_top,
      y0_bot = y0_bot,
      y1_top = y1_top,
      y1_bot = y1_bot,
      ribbon_id = paste0("D_", cl, "_", th)
    ) %>%
      mutate(
        cluster = cl,
        theme   = th,
        fill    = CLUSTER_COLORS[cl]
      )

    ribbon_list[[ribbon_idx]] <- ribbon_poly
  }
}

D_ribbons <- bind_rows(ribbon_list)
cat(sprintf("  Built %d ribbons for Panel D Sankey\n", ribbon_idx))

# --- Step 3: Build stacked bar geometry in the same coordinate system ---------
# Layout zones (x-coordinates):
#   D_X_CL (1.0)              — cluster bars
#   D_X_TH (3.0)              — theme bars
#   D_X_LABEL (3.2 - 5.2)     — theme name labels (right-aligned text zone)
#   D_X_BAR_START (5.4)       — stacked bars begin
#   D_X_BAR_END (~7.0)        — stacked bars end + total labels

D_X_BAR_START <- 5.4
D_MAX_BAR_LEN <- 1.5
max_theme_count <- max(sankey_theme_totals$total)

# Stacked bar height: use a fixed bar height centered on each theme bar center
D_SBAR_H <- 2.2   # fixed height for all stacked bars

# Build stacked bar rectangles: for each theme, one rect per cluster contribution
stacked_rects <- list()
stacked_idx <- 0

for (th in sankey_theme_order) {
  th_row <- th_bars %>% filter(theme == th)
  th_contribs <- sankey_theme_counts %>%
    filter(theme == th) %>%
    arrange(cluster)

  # Center stacked bar on theme bar center
  bar_y_ctr <- th_row$y_ctr
  bar_y_top <- bar_y_ctr + D_SBAR_H / 2
  bar_y_bot <- bar_y_ctr - D_SBAR_H / 2

  # Stack cluster segments left to right
  x_cursor <- D_X_BAR_START

  for (r in seq_len(nrow(th_contribs))) {
    cl <- th_contribs$cluster[r]
    n  <- th_contribs$n_proteins[r]
    if (n == 0) next

    stacked_idx <- stacked_idx + 1
    seg_w <- (n / max_theme_count) * D_MAX_BAR_LEN

    stacked_rects[[stacked_idx]] <- tibble(
      xmin = x_cursor,
      xmax = x_cursor + seg_w,
      ymin = bar_y_bot,
      ymax = bar_y_top,
      cluster = cl,
      theme   = th,
      fill    = CLUSTER_COLORS[cl]
    )

    x_cursor <- x_cursor + seg_w
  }
}

D_stacked <- bind_rows(stacked_rects)

# Total labels at bar tips
D_bar_totals <- sankey_theme_totals %>%
  mutate(
    x_tip = D_X_BAR_START + (total / max_theme_count) * D_MAX_BAR_LEN
  ) %>%
  left_join(th_bars %>% dplyr::select(theme, y_ctr), by = "theme")

cat(sprintf("  Built %d stacked bar segments\n", stacked_idx))

# --- Build the combined Sankey + stacked bars plot ---
p_D_sankey <- ggplot() +
  # Ribbons (drawn first = behind bars)
  geom_polygon(
    data = D_ribbons,
    aes(x = x, y = y, group = ribbon_id, fill = fill),
    alpha = 0.30, color = NA
  ) +
  # Cluster bars (left)
  geom_rect(
    data = cl_bars,
    aes(xmin = x_left, xmax = x_right, ymin = y_bot, ymax = y_top),
    fill = cl_bars$fill, color = "black", linewidth = 0.3
  ) +
  # Cluster labels (inside bars, white text)
  geom_text(
    data = cl_bars,
    aes(x = (x_left + x_right) / 2, y = y_ctr, label = cluster),
    color = "white", fontface = "bold", size = 2.8
  ) +
  # Theme bars (right side of Sankey)
  geom_rect(
    data = th_bars,
    aes(xmin = x_left, xmax = x_right, ymin = y_bot, ymax = y_top),
    fill = "grey85", color = "black", linewidth = 0.3
  ) +
  # Theme labels (in dedicated zone between theme bars and stacked bars)
  geom_text(
    data = th_bars,
    aes(x = x_right + 0.12, y = y_ctr, label = theme),
    hjust = 0, size = 2.3, fontface = "bold",
    color = "grey20"
  ) +
  # Stacked bars (aligned to theme bar y-centers, fixed height)
  geom_rect(
    data = D_stacked,
    aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
    fill = D_stacked$fill, color = "black", linewidth = 0.2
  ) +
  # Total count at bar tips
  geom_text(
    data = D_bar_totals,
    aes(x = x_tip + 0.06, y = y_ctr, label = total),
    hjust = 0, size = KEY_TEXT, fontface = "bold", color = "grey30"
  ) +
  # X-axis label for stacked bars
  annotate("text",
           x = D_X_BAR_START + D_MAX_BAR_LEN / 2,
           y = -3,
           label = "Protein count", size = 2.2, color = "grey40") +
  scale_fill_identity() +
  coord_cartesian(
    xlim = c(0.5, D_X_BAR_START + D_MAX_BAR_LEN + 0.7),
    ylim = c(-5, D_Y_SPAN + 2),
    expand = FALSE
  ) +
  theme_void() +
  theme(plot.margin = margin(t = 4, r = 4, b = 4, l = 2))

cat("  Panel D combined plot built\n")

# --- Step 4: Panel D is now a single unified plot ---
panel_D <- p_D_sankey

cat("  Panel D composed\n")

# --- Verification: save PDF ---
ggsave(file.path(RPT_DIR, "panel_D_synthesis.pdf"), panel_D,
       width = 200, height = 280, units = "mm")
cat(sprintf("  Saved verification PDF: %s\n", file.path(RPT_DIR, "panel_D_synthesis.pdf")))
