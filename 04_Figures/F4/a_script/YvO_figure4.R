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
#   Panels (to be added in later tasks):
#     A — FCM cluster centroid profiles (spaghetti + ribbon)
#     B — Z-score heatmap stack per cluster
#     C — Pre/Post trajectory plots per cluster (Young vs Old)
#     D — Pathway enrichment per cluster (Hallmark + GO:BP)
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

cat("=== Task 1 complete: Setup + Mfuzz clustering ===\n")

# === 14. PANEL A — FCM CLUSTER PROFILES (Interaction Plot) ================

cat("Building Panel A: FCM cluster profiles (interaction plot)...\n")

n_young <- length(young_subjects)
n_old   <- length(old_subjects)
n_subj  <- n_young + n_old

# Compute group-mean expression per protein across 4 conditions
group_means_A <- matrix(NA_real_, nrow = nrow(abund_mat), ncol = 4,
                        dimnames = list(rownames(abund_mat), GROUP_COLS))
for (grp in GROUP_COLS) {
  grp_samples <- sample_meta$sample[paste0(sample_meta$age, "_", sample_meta$time) == grp]
  group_means_A[, grp] <- rowMeans(abund_mat[, grp_samples, drop = FALSE], na.rm = TRUE)
}

# Z-score per protein across the 4 group means
group_z_A <- t(scale(t(group_means_A)))

# Only keep proteins that survived clustering
keep_genes_A <- intersect(rownames(group_z_A), cluster_assign$gene)
group_z_A <- group_z_A[keep_genes_A, ]

# Reshape for interaction plot: age x time
profile_interaction <- as.data.frame(group_z_A) |>
  rownames_to_column("gene") |>
  pivot_longer(-gene, names_to = "group", values_to = "z_score") |>
  left_join(cluster_assign |> dplyr::select(gene, cluster, membership),
            by = "gene") |>
  mutate(
    age  = ifelse(str_detect(group, "^Young"), "Young", "Old"),
    time = ifelse(str_detect(group, "_Post$"), "Post", "Pre"),
    time_num = ifelse(time == "Pre", 1, 2)
  )

# Compute mean +/- SE per cluster x age x time
profile_summary <- profile_interaction |>
  group_by(cluster, age, time, time_num) |>
  summarise(
    mean_z = mean(z_score, na.rm = TRUE),
    se_z   = sd(z_score, na.rm = TRUE) / sqrt(n()),
    .groups = "drop"
  )

# Shared y-axis range
y_range_A <- range(
  profile_summary$mean_z - profile_summary$se_z,
  profile_summary$mean_z + profile_summary$se_z
) * c(1.05, 1.05)

# Build per-cluster interaction plots
panel_A_list <- lapply(seq_len(optimal_k), function(ci) {
  cl_id  <- paste0("C", ci)
  cl_sum <- profile_summary |> filter(cluster == cl_id)
  n_total <- sum(cluster_assign$cluster == cl_id)
  n_core  <- sum(cluster_assign$membership[cluster_assign$cluster == cl_id] >= 0.7)

  # Aging gap at Pre
  pre_vals  <- cl_sum |> filter(time == "Pre")
  young_pre <- pre_vals$mean_z[pre_vals$age == "Young"]
  old_pre   <- pre_vals$mean_z[pre_vals$age == "Old"]

  sub_text <- paste0("(n = ", n_total, ", core = ", n_core, ")")

  p <- ggplot(cl_sum, aes(x = time_num, y = mean_z, color = age,
                           fill = age, group = age)) +
    geom_ribbon(aes(ymin = mean_z - se_z, ymax = mean_z + se_z),
                alpha = 0.15, color = NA) +
    geom_line(linewidth = 1.2) +
    geom_point(size = 2.5) +
    annotate("segment", x = 1, xend = 1,
             y = min(young_pre, old_pre), yend = max(young_pre, old_pre),
             linetype = "dashed", color = AGING_GAP_LINE, linewidth = 0.5) +
    annotate("text", x = 0.85, y = (young_pre + old_pre) / 2,
             label = expression(Delta * "age"), size = 2, color = AGING_GAP_LINE,
             fontface = "italic") +
    scale_color_manual(values = AGE_COLORS, guide = "none") +
    scale_fill_manual(values = AGE_COLORS, guide = "none") +
    scale_x_continuous(breaks = c(1, 2), labels = c("Pre", "Post"),
                       limits = c(0.7, 2.3)) +
    coord_cartesian(ylim = y_range_A) +
    labs(title    = paste0("Cluster ", ci),
         subtitle = sub_text,
         x = NULL,
         y = if (ci == 1) "Z-score (group means)" else NULL) +
    THEME_PUB +
    theme(plot.title = element_text(color = CLUSTER_COLORS[cl_id]))

  # Suppress x-axis on all except bottom cluster
  if (ci < optimal_k) {
    p <- p + theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())
  } else {
    p <- p + theme(axis.text.x = element_text(size = 5))
  }

  # Only show y-axis on first panel
  if (ci > 1) {
    p <- p + theme(axis.text.y = element_blank(), axis.ticks.y = element_blank())
  }

  p
})

cat("  Panel A complete\n")
cat("=== Task 2 complete: Panel A — FCM Cluster Profiles ===\n")

# === 15. PANEL B — CONTRAST VIOLIN DISTRIBUTIONS ============================

cat("Building Panel B: Contrast violin distributions...\n")

# Compute group means from raw imputed matrix for each protein
sample_groups <- sample_meta |>
  mutate(group = paste0(age, "_", time))

group_means <- matrix(NA_real_, nrow = nrow(abund_mat), ncol = 4,
                      dimnames = list(rownames(abund_mat), GROUP_COLS))

for (grp in GROUP_COLS) {
  grp_samples <- sample_meta$sample[paste0(sample_meta$age, "_", sample_meta$time) == grp]
  group_means[, grp] <- rowMeans(abund_mat[, grp_samples, drop = FALSE], na.rm = TRUE)
}

# Compute 4 contrasts from group means (raw logFC, NOT z-scored)
CONTRAST_COLS <- c("Aging", "Training_Young", "Training_Old", "Interaction")
CONTRAST_LABS <- c("Aging", "Training\nYoung", "Training\nOld", "Interaction")

contrast_mat <- matrix(NA_real_, nrow = nrow(group_means), ncol = 4,
                       dimnames = list(rownames(group_means), CONTRAST_COLS))
contrast_mat[, "Aging"]          <- group_means[, "Old_Pre"]    - group_means[, "Young_Pre"]
contrast_mat[, "Training_Young"] <- group_means[, "Young_Post"] - group_means[, "Young_Pre"]
contrast_mat[, "Training_Old"]   <- group_means[, "Old_Post"]   - group_means[, "Old_Pre"]
contrast_mat[, "Interaction"]    <- (group_means[, "Old_Post"] - group_means[, "Old_Pre"]) -
                                    (group_means[, "Young_Post"] - group_means[, "Young_Pre"])

# Only keep proteins that are in cluster_assign
keep_genes <- intersect(rownames(contrast_mat), cluster_assign$gene)

# Order: by cluster, then by descending membership within cluster
protein_order <- cluster_assign |>
  filter(gene %in% keep_genes) |>
  arrange(cluster, desc(membership))

# Build long format for violins
violin_long <- as.data.frame(contrast_mat[protein_order$gene, ]) |>
  rownames_to_column("gene") |>
  pivot_longer(all_of(CONTRAST_COLS), names_to = "contrast", values_to = "logFC") |>
  left_join(protein_order |> dplyr::select(gene, cluster), by = "gene") |>
  mutate(contrast = factor(contrast, levels = rev(CONTRAST_COLS)))

# Shared x-axis range
x_range_B <- range(violin_long$logFC, na.rm = TRUE) * c(1.05, 1.05)

# Build per-cluster violin sub-plots
panel_B_list <- lapply(seq_len(optimal_k), function(ci) {
  cl_id <- paste0("C", ci)
  cl_viol <- violin_long |> filter(cluster == cl_id)

  p <- ggplot(cl_viol, aes(y = contrast, x = logFC, fill = contrast)) +
    geom_violin(scale = "width", alpha = 0.7, color = NA, linewidth = 0.2) +
    geom_boxplot(width = 0.15, outlier.shape = NA, fill = "white",
                 alpha = 0.6, linewidth = 0.3) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "grey40",
               linewidth = 0.3) +
    scale_fill_manual(values = CONTRAST_COLORS, guide = "none") +
    scale_y_discrete(labels = rev(CONTRAST_LABS)) +
    coord_cartesian(xlim = x_range_B) +
    labs(x = NULL, y = NULL) +
    THEME_PUB +
    theme(
      panel.border = element_rect(color = CLUSTER_COLORS[cl_id],
                                  linewidth = 1.2, fill = NA)
    )

  # Suppress x-axis on all except bottom cluster
  if (ci < optimal_k) {
    p <- p + theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())
  } else {
    p <- p + theme(axis.text.x = element_text(size = 5))
  }

  # Only show contrast labels on first cluster
  if (ci > 1) {
    p <- p + theme(axis.text.y = element_blank(), axis.ticks.y = element_blank())
  }

  p
})

cat("  Panel B complete (per-cluster violin list)\n")

# === 16. PANEL C — GROUP TRAJECTORY STACK ====================================

cat("Building Panel C: Group trajectories...\n")

# For each cluster, compute per-sample cluster-mean expression,
# then group means +/- SE

# Per-sample mean for each cluster
traj_data_list <- lapply(paste0("C", seq_len(optimal_k)), function(cl_id) {
  cl_genes <- cluster_assign$gene[cluster_assign$cluster == cl_id]
  cl_genes_in_mat <- intersect(cl_genes, rownames(abund_mat))

  # Per-sample mean across cluster genes
  sample_means <- colMeans(abund_mat[cl_genes_in_mat, , drop = FALSE], na.rm = TRUE)

  tibble(
    sample = names(sample_means),
    value = sample_means,
    cluster = cl_id
  ) |>
    left_join(sample_meta, by = "sample") |>
    mutate(group = paste0(age, "_", time))
})

traj_data <- bind_rows(traj_data_list)

# Compute group means +/- SE per cluster
traj_summary <- traj_data |>
  group_by(cluster, age, time) |>
  summarise(
    mean_val = mean(value, na.rm = TRUE),
    se_val = sd(value, na.rm = TRUE) / sqrt(n()),
    .groups = "drop"
  ) |>
  mutate(
    time_num = ifelse(time == "Pre", 1, 2),
    cluster = factor(cluster, levels = paste0("C", seq_len(optimal_k)))
  )

# Export trajectory data
write_csv(traj_summary, file.path(DAT_DIR, "fig4_group_trajectories.csv"))
cat(sprintf("  Saved: %s\n", file.path(DAT_DIR, "fig4_group_trajectories.csv")))

# Shared y range across all trajectory panels
y_range <- range(traj_summary$mean_val - traj_summary$se_val,
                 traj_summary$mean_val + traj_summary$se_val) * c(0.98, 1.02)

# Build per-cluster individual trajectory panels
panel_C_list <- lapply(paste0("C", seq_len(optimal_k)), function(cl_id) {
  cl_traj <- traj_data |> filter(cluster == cl_id) |>
    mutate(time_num = ifelse(time == "Pre", 1, 2))
  cl_sum  <- traj_summary |> filter(cluster == cl_id)

  # Per-cluster auto-range y-axis (tight fit for visual drama)
  y_range_cl <- range(cl_traj$value, na.rm = TRUE)
  y_pad      <- diff(y_range_cl) * 0.02
  y_range_cl <- y_range_cl + c(-y_pad, y_pad)

  # Aging gap at Pre (from group means)
  pre_vals  <- cl_sum |> filter(time == "Pre")
  young_pre <- pre_vals$mean_val[pre_vals$age == "Young"]
  old_pre   <- pre_vals$mean_val[pre_vals$age == "Old"]

  ci <- as.integer(sub("C", "", cl_id))

  p <- ggplot() +
    # Individual subject lines
    geom_line(data = cl_traj,
              aes(x = time_num, y = value, group = subject, color = age),
              linewidth = 0.4, alpha = 0.5) +
    # Group mean overlay
    geom_line(data = cl_sum,
              aes(x = time_num, y = mean_val, color = age, group = age),
              linewidth = 1.4) +
    geom_point(data = cl_sum,
               aes(x = time_num, y = mean_val, color = age),
               size = 2.5) +
    # Aging gap bracket at Pre
    annotate("segment", x = 1, xend = 1,
             y = min(young_pre, old_pre), yend = max(young_pre, old_pre),
             linetype = "dashed", color = AGING_GAP_LINE, linewidth = 0.5) +
    annotate("text", x = 0.85, y = (young_pre + old_pre) / 2,
             label = expression(Delta * "age"), size = 2, color = AGING_GAP_LINE,
             fontface = "italic") +
    scale_color_manual(values = AGE_COLORS, guide = "none") +
    scale_x_continuous(breaks = c(1, 2), labels = c("Pre", "Post"),
                       limits = c(0.7, 2.3)) +
    coord_cartesian(ylim = y_range_cl) +
    labs(x = NULL,
         y = if (cl_id == "C1") "Mean Intensity\n(cluster proteins)" else NULL) +
    THEME_PUB

  # Suppress x-axis on all except bottom cluster
  if (ci < optimal_k) {
    p <- p + theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())
  } else {
    p <- p + theme(axis.text.x = element_text(size = 5))
  }

  # Only show y-axis labels on first panel
  if (cl_id != "C1") {
    p <- p + theme(axis.text.y = element_blank(), axis.ticks.y = element_blank())
  }

  p
})

# Keep as list for per-cluster row assembly (do NOT stack vertically)

cat("  Panel C complete\n")
cat("=== Task 3 complete: Panel B + Panel C — Z-Score Heatmap + Group Trajectories ===\n")

# === 17. PANEL D — PER-CLUSTER PATHWAY ENRICHMENT DOTPLOT ====================

cat("Building Panel D: Per-cluster pathway enrichment...\n")

# ------ Step 1: Prepare gene lists -------------------------------------------

universe_genes <- unique(cluster_assign$gene)

# ------ Step 2: Hallmark enrichment via enricher() ---------------------------

cat("  Loading Hallmark gene sets...\n")
hallmark_t2g <- msigdbr(species = "Homo sapiens", collection = "H") |>
  dplyr::select(gs_name, gene_symbol)
cat(sprintf("    Hallmark: %d gene-set-gene pairs\n", nrow(hallmark_t2g)))

# ------ Step 3: GO:BP enrichment via enrichGO() + rrvgo reduction ------------
# ------ Step 4: Per-cluster enrichment loop ----------------------------------

cat("  Running per-cluster enrichment (Hallmark + GO:BP + GO:CC)...\n")

# Pre-compute GO semantic similarity data for rrvgo (reuse across clusters)
cat("    Pre-computing GO:BP semantic similarity data for rrvgo...\n")
semdata_bp <- tryCatch(
  GOSemSim::godata("org.Hs.eg.db", ont = "BP", keytype = "ENTREZID"),
  error = function(e) {
    cat("    WARNING: Could not pre-compute BP semdata:", conditionMessage(e), "\n")
    NULL
  }
)

cat("    Pre-computing GO:CC semantic similarity data for rrvgo...\n")
semdata_cc <- tryCatch(
  GOSemSim::godata("org.Hs.eg.db", ont = "CC", keytype = "ENTREZID"),
  error = function(e) {
    cat("    WARNING: Could not pre-compute CC semdata:", conditionMessage(e), "\n")
    NULL
  }
)

# Map entire universe to ENTREZID once (reuse across clusters)
entrez_universe <- tryCatch(
  bitr(universe_genes, fromType = "SYMBOL", toType = "ENTREZID",
       OrgDb = org.Hs.eg.db),
  error = function(e) {
    cat("    WARNING: bitr universe mapping failed:", conditionMessage(e), "\n")
    tibble(SYMBOL = character(), ENTREZID = character())
  }
)
cat(sprintf("    Universe: %d genes mapped to ENTREZID (of %d)\n",
            nrow(entrez_universe), length(universe_genes)))

enrich_results_list <- lapply(paste0("C", seq_len(optimal_k)), function(cl_id) {
  cat(sprintf("    Cluster %s ...\n", cl_id))

  # Genes with membership >= 0.5 in this cluster
  cl_genes <- cluster_assign$gene[cluster_assign$cluster == cl_id &
                                    cluster_assign$membership >= 0.5]
  cl_genes <- unique(cl_genes[!is.na(cl_genes) & cl_genes != ""])
  cat(sprintf("      %d genes with membership >= 0.5\n", length(cl_genes)))

  if (length(cl_genes) < 5) {
    cat("      Skipping: too few genes\n")
    return(tibble())
  }

  # --- Hallmark enrichment ---
  h_res <- tryCatch(
    enricher(cl_genes, TERM2GENE = hallmark_t2g, universe = universe_genes,
             pvalueCutoff = 0.05, pAdjustMethod = "BH"),
    error = function(e) {
      cat(sprintf("      Hallmark enricher error: %s\n", conditionMessage(e)))
      NULL
    }
  )

  h_df <- if (!is.null(h_res) && nrow(as.data.frame(h_res)) > 0) {
    as.data.frame(h_res) |>
      mutate(database = "Hallmark", cluster = cl_id) |>
      head(5)
  } else {
    cat("      No significant Hallmark terms\n")
    tibble()
  }

  # --- GO:BP enrichment ---
  entrez_map <- tryCatch(
    bitr(cl_genes, fromType = "SYMBOL", toType = "ENTREZID",
         OrgDb = org.Hs.eg.db),
    error = function(e) {
      cat(sprintf("      bitr mapping error: %s\n", conditionMessage(e)))
      tibble(SYMBOL = character(), ENTREZID = character())
    }
  )

  bp_res <- tryCatch(
    enrichGO(gene = entrez_map$ENTREZID,
             universe = entrez_universe$ENTREZID,
             OrgDb = org.Hs.eg.db, ont = "BP",
             pAdjustMethod = "BH", pvalueCutoff = 0.05,
             readable = TRUE),
    error = function(e) {
      cat(sprintf("      enrichGO error: %s\n", conditionMessage(e)))
      NULL
    }
  )

  bp_df <- if (!is.null(bp_res) && nrow(as.data.frame(bp_res)) > 0) {
    bp_full <- as.data.frame(bp_res)
    cat(sprintf("      GO:BP raw hits: %d\n", nrow(bp_full)))

    # Reduce with rrvgo (threshold 0.85)
    if (nrow(bp_full) > 1 && !is.null(semdata_bp)) {
      sim_matrix <- tryCatch({
        calculateSimMatrix(bp_full$ID, orgdb = "org.Hs.eg.db",
                           semdata = semdata_bp, ont = "BP", method = "Rel")
      }, error = function(e) {
        cat(sprintf("      rrvgo simMatrix error: %s\n", conditionMessage(e)))
        NULL
      })

      if (!is.null(sim_matrix) && nrow(sim_matrix) > 0) {
        scores_vec <- setNames(-log10(bp_full$p.adjust), bp_full$ID)
        # Only keep IDs that appear in sim_matrix rows
        scores_vec <- scores_vec[intersect(names(scores_vec),
                                           rownames(sim_matrix))]
        reduced <- tryCatch(
          reduceSimMatrix(sim_matrix, scores = scores_vec,
                          threshold = 0.85, orgdb = "org.Hs.eg.db"),
          error = function(e) {
            cat(sprintf("      rrvgo reduce error: %s\n", conditionMessage(e)))
            NULL
          }
        )

        if (!is.null(reduced) && nrow(reduced) > 0) {
          parent_terms <- unique(reduced$parentTerm)
          bp_full <- bp_full |> filter(Description %in% parent_terms |
                                         ID %in% parent_terms)
          cat(sprintf("      GO:BP after rrvgo: %d parent terms\n",
                      nrow(bp_full)))
        }
      }
    }

    bp_full |>
      mutate(database = "GO:BP", cluster = cl_id) |>
      head(5)
  } else {
    cat("      No significant GO:BP terms\n")
    tibble()
  }

  # --- GO:CC enrichment ---
  cc_res <- tryCatch(
    enrichGO(gene = entrez_map$ENTREZID,
             universe = entrez_universe$ENTREZID,
             OrgDb = org.Hs.eg.db, ont = "CC",
             pAdjustMethod = "BH", pvalueCutoff = 0.05,
             readable = TRUE),
    error = function(e) {
      cat(sprintf("      enrichGO CC error: %s\n", conditionMessage(e)))
      NULL
    }
  )

  cc_df <- if (!is.null(cc_res) && nrow(as.data.frame(cc_res)) > 0) {
    cc_full <- as.data.frame(cc_res)
    cat(sprintf("      GO:CC raw hits: %d\n", nrow(cc_full)))

    # Reduce with rrvgo (threshold 0.85)
    if (nrow(cc_full) > 1 && !is.null(semdata_cc)) {
      sim_matrix_cc <- tryCatch({
        calculateSimMatrix(cc_full$ID, orgdb = "org.Hs.eg.db",
                           semdata = semdata_cc, ont = "CC", method = "Rel")
      }, error = function(e) {
        cat(sprintf("      rrvgo CC simMatrix error: %s\n", conditionMessage(e)))
        NULL
      })

      if (!is.null(sim_matrix_cc) && nrow(sim_matrix_cc) > 0) {
        scores_vec_cc <- setNames(-log10(cc_full$p.adjust), cc_full$ID)
        scores_vec_cc <- scores_vec_cc[intersect(names(scores_vec_cc),
                                                 rownames(sim_matrix_cc))]
        reduced_cc <- tryCatch(
          reduceSimMatrix(sim_matrix_cc, scores = scores_vec_cc,
                          threshold = 0.85, orgdb = "org.Hs.eg.db"),
          error = function(e) {
            cat(sprintf("      rrvgo CC reduce error: %s\n", conditionMessage(e)))
            NULL
          }
        )

        if (!is.null(reduced_cc) && nrow(reduced_cc) > 0) {
          parent_terms_cc <- unique(reduced_cc$parentTerm)
          cc_full <- cc_full |> filter(Description %in% parent_terms_cc |
                                         ID %in% parent_terms_cc)
          cat(sprintf("      GO:CC after rrvgo: %d parent terms\n",
                      nrow(cc_full)))
        }
      }
    }

    cc_full |>
      mutate(database = "GO:CC", cluster = cl_id) |>
      head(5)
  } else {
    cat("      No significant GO:CC terms\n")
    tibble()
  }

  bind_rows(h_df, bp_df, cc_df)
})

enrich_combined <- bind_rows(enrich_results_list)
cat(sprintf("  Combined enrichment: %d rows\n", nrow(enrich_combined)))

# Export enrichment data
write_csv(enrich_combined, file.path(DAT_DIR, "mfuzz_enrichment.csv"))
cat(sprintf("  Saved: %s\n", file.path(DAT_DIR, "mfuzz_enrichment.csv")))

# ------ Step 5: Update Panel A subtitles with top Hallmark term per cluster --

cat("  Updating Panel A subtitles with top Hallmark enrichment labels...\n")

# Extract top Hallmark term per cluster
top_hallmark_per_cluster <- enrich_combined |>
  filter(database == "Hallmark") |>
  group_by(cluster) |>
  slice_min(p.adjust, n = 1, with_ties = FALSE) |>
  ungroup() |>
  mutate(label = clean_pathway_name(Description))

# Build a lookup: cluster -> biology label
bio_labels <- setNames(
  rep("", optimal_k),
  paste0("C", seq_len(optimal_k))
)
for (i in seq_len(nrow(top_hallmark_per_cluster))) {
  bio_labels[top_hallmark_per_cluster$cluster[i]] <- top_hallmark_per_cluster$label[i]
}

cat("  Biology labels per cluster:\n")
for (cl_id in paste0("C", seq_len(optimal_k))) {
  cat(sprintf("    %s: %s\n", cl_id, bio_labels[cl_id]))
}

# Rebuild Panel A with biology labels (interaction plot)
panel_A_list <- lapply(seq_len(optimal_k), function(ci) {
  cl_id  <- paste0("C", ci)
  cl_sum <- profile_summary |> filter(cluster == cl_id)
  n_total <- sum(cluster_assign$cluster == cl_id)
  n_core  <- sum(cluster_assign$membership[cluster_assign$cluster == cl_id] >= 0.7)

  # Aging gap at Pre
  pre_vals  <- cl_sum |> filter(time == "Pre")
  young_pre <- pre_vals$mean_z[pre_vals$age == "Young"]
  old_pre   <- pre_vals$mean_z[pre_vals$age == "Old"]

  bio_lbl <- bio_labels[cl_id]
  if (nchar(bio_lbl) > 0) {
    sub_text <- paste0("(n = ", n_total, ", core = ", n_core, ")  ", bio_lbl)
  } else {
    sub_text <- paste0("(n = ", n_total, ", core = ", n_core, ")")
  }

  p <- ggplot(cl_sum, aes(x = time_num, y = mean_z, color = age,
                           fill = age, group = age)) +
    geom_ribbon(aes(ymin = mean_z - se_z, ymax = mean_z + se_z),
                alpha = 0.15, color = NA) +
    geom_line(linewidth = 1.2) +
    geom_point(size = 2.5) +
    annotate("segment", x = 1, xend = 1,
             y = min(young_pre, old_pre), yend = max(young_pre, old_pre),
             linetype = "dashed", color = AGING_GAP_LINE, linewidth = 0.5) +
    annotate("text", x = 0.85, y = (young_pre + old_pre) / 2,
             label = expression(Delta * "age"), size = 2, color = AGING_GAP_LINE,
             fontface = "italic") +
    scale_color_manual(values = AGE_COLORS, guide = "none") +
    scale_fill_manual(values = AGE_COLORS, guide = "none") +
    scale_x_continuous(breaks = c(1, 2), labels = c("Pre", "Post"),
                       limits = c(0.7, 2.3)) +
    coord_cartesian(ylim = y_range_A) +
    labs(title    = paste0("Cluster ", ci),
         subtitle = sub_text,
         x = NULL,
         y = if (ci == 1) "Z-score (group means)" else NULL) +
    THEME_PUB +
    theme(plot.title = element_text(color = CLUSTER_COLORS[cl_id]))

  if (ci < optimal_k) {
    p <- p + theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())
  } else {
    p <- p + theme(axis.text.x = element_text(size = 5))
  }

  if (ci > 1) {
    p <- p + theme(axis.text.y = element_blank(), axis.ticks.y = element_blank())
  }

  p
})

cat("  Panel A updated with biology labels\n")

# ------ Build Panel D — Separate per-cluster enrichment plots ---------------

if (nrow(enrich_combined) > 0) {
  # Clean pathway names and compute -log10(padj)
  enrich_plot_df <- enrich_combined |>
    mutate(
      term_clean = clean_pathway_name(Description),
      neg_log10_padj = -log10(p.adjust),
      # Parse GeneRatio
      gene_ratio = sapply(strsplit(GeneRatio, "/"), function(x) {
        as.numeric(x[1]) / as.numeric(x[2])
      }),
      cluster = factor(cluster, levels = paste0("C", seq_len(optimal_k)))
    )

  # Build separate per-cluster panels
  panel_D_list <- list()

  for (ci in seq_len(optimal_k)) {
    cl_id <- paste0("C", ci)
    cl_color <- CLUSTER_COLORS[cl_id]

    cl_enrich <- enrich_plot_df |> filter(cluster == cl_id)
    if (nrow(cl_enrich) == 0) next

    # Build facet label: "C1: mTORC1 Signaling"
    bio_lbl <- bio_labels[cl_id]
    facet_label <- if (nchar(bio_lbl) > 0) {
      paste0(cl_id, ": ", bio_lbl)
    } else {
      cl_id
    }
    cl_enrich <- cl_enrich |>
      mutate(
        label = facet_label,
        term_clean = fct_reorder(term_clean, neg_log10_padj)
      )

    p <- ggplot(cl_enrich, aes(x = neg_log10_padj, y = term_clean)) +
      geom_point(aes(size = gene_ratio, fill = database),
                 shape = 21, stroke = 0.8, color = cl_color) +
      geom_vline(xintercept = -log10(0.05), linetype = "dashed",
                 color = "grey40", linewidth = 0.3) +
      scale_fill_manual(values = c("Hallmark" = "#AA336A",
                                   "GO:BP" = "#00796B",
                                   "GO:CC" = "#26A69A"),
                        name = "Database") +
      scale_size_continuous(range = c(1.5, 5), name = "Gene Ratio") +
      facet_wrap(~ label) +
      labs(x = expression(-log[10](p.adjust)), y = NULL) +
      THEME_PUB +
      theme(
        strip.background = element_rect(fill = cl_color, color = cl_color),
        strip.text = element_text(color = "white", face = "bold", size = 7),
        panel.border = element_rect(color = alpha(cl_color, 0.4),
                                    linewidth = 0.8, fill = NA),
        panel.background = element_rect(fill = alpha(cl_color, 0.04)),
        legend.position = "none"
      )

    panel_D_list[[cl_id]] <- p
  }

  if (length(panel_D_list) > 0) {
    panel_D <- Reduce(`|`, panel_D_list) +
      plot_layout(guides = "collect") +
      plot_annotation(
        title = "D  Per-Cluster Pathway Enrichment",
        subtitle = "ORA: Hallmark + GO:BP + GO:CC (rrvgo-reduced, threshold = 0.85)",
        theme = THEME_PUB
      ) &
      theme(legend.position = "bottom")
  } else {
    cat("  WARNING: All clusters empty — creating placeholder Panel D\n")
    panel_D <- ggplot() +
      annotate("text", x = 0.5, y = 0.5,
               label = "No significant enrichment terms", size = 4) +
      labs(title = "D  Per-Cluster Pathway Enrichment") +
      THEME_PUB + theme_void()
  }
} else {
  # Fallback: empty Panel D if no enrichment results
  cat("  WARNING: No enrichment results — creating placeholder Panel D\n")
  panel_D <- ggplot() +
    annotate("text", x = 0.5, y = 0.5,
             label = "No significant enrichment terms", size = 4) +
    labs(title = "D  Per-Cluster Pathway Enrichment") +
    THEME_PUB +
    theme_void()
}

cat("  Panel D complete\n")
cat("=== Task 4 complete: Panel D + Enrichment ===\n")

# === 18. CLUSTER COLOR BAR STRIPS =============================================

cat("Building cluster color bar strips...\n")

# Get cluster protein counts for labels
cluster_ns <- protein_order |>
  group_by(cluster) |>
  summarise(n = n(), .groups = "drop") |>
  arrange(cluster) |>
  pull(n)

color_bar_list <- lapply(seq_len(optimal_k), function(ci) {
  cl_id <- paste0("C", ci)
  ggplot() +
    annotate("rect", xmin = 0, xmax = 1, ymin = 0, ymax = 1,
             fill = CLUSTER_COLORS[cl_id]) +
    annotate("text", x = 0.5, y = 0.5,
             label = paste0(cl_id, "\n(", cluster_ns[ci], ")"),
             color = "white", fontface = "bold", size = 2.5, lineheight = 0.9) +
    theme_void() +
    theme(plot.margin = margin(0, 0, 0, 0))
})

cat("  Color bars complete\n")

# === 19. FINAL FIGURE ASSEMBLY ================================================

cat("Assembling Figure 4 (cluster-aligned layout)...\n")

# Column widths: Panel A (22%) | color bar (3%) | Panel B (45%) | Panel C (30%)
col_widths <- c(0.22, 0.03, 0.45, 0.30)

# Header row with panel titles
header_A   <- ggplot() + annotate("text", x = 0.5, y = 0.5,
               label = "A  Cluster Profiles", fontface = "bold", size = 3.5) +
               theme_void()
header_bar <- ggplot() + theme_void()
header_B   <- ggplot() + annotate("text", x = 0.5, y = 0.5,
               label = "B  Contrast Heatmap", fontface = "bold", size = 3.5) +
               theme_void()
header_C   <- ggplot() + annotate("text", x = 0.5, y = 0.5,
               label = "C  Trajectories", fontface = "bold", size = 3.5) +
               theme_void()

header_row <- (header_A | header_bar | header_B | header_C) +
  plot_layout(widths = col_widths)

# Per-cluster rows: Panel A | color bar | Panel B | Panel C
cluster_rows <- lapply(seq_len(optimal_k), function(ci) {
  (panel_A_list[[ci]] | color_bar_list[[ci]] | panel_B_list[[ci]] | panel_C_list[[ci]]) +
    plot_layout(widths = col_widths)
})

# Row heights proportional to cluster protein counts
cluster_heights <- cluster_ns / sum(cluster_ns) * 0.60

# Stack: header + cluster rows + Panel D
fig4 <- Reduce(`/`, c(list(header_row), cluster_rows, list(panel_D))) +
  plot_layout(heights = c(0.03, cluster_heights, 0.37))

ggsave(file.path(RPT_DIR, "Figure_4.pdf"), fig4,
       width = 380, height = 500, units = "mm", device = pdf, bg = "white")
cat(sprintf("  Saved: %s\n", file.path(RPT_DIR, "Figure_4.pdf")))

ggsave(file.path(RPT_DIR, "Figure_4.png"), fig4,
       width = 380, height = 500, units = "mm", dpi = 300, bg = "white")
cat(sprintf("  Saved: %s\n", file.path(RPT_DIR, "Figure_4.png")))

cat("=== Figure 4 complete ===\n")
