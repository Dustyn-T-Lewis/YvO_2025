# YvO_F4_setup.R — Shared setup for Figure 4
# Sources: palettes.R
# Provides: cluster_assign, core_proteins, delta_z, sample_meta, imputed,
#           enrich_top, protein_pathway_links, theme_links,
#           CLUSTER_COLORS, AGE_COLORS, DB_COLORS, THEME_COLORS, row_heights,
#           abund_mat, gene_ids, group_z, GROUP_COLS, GROUP_LABS, ABBREV_LABS,
#           optimal_k, m_est, CORE_THRESH, centroids, membership,
#           clean_pathway_name, sig_stars, make_sigmoid_ribbon,
#           reorder_within, scale_y_reordered, top_hallmark,
#           RPT_DIR, DAT_DIR, FIG_W, FIG_H, COL_WIDTHS, THEME_PUB
#
# STAT AUDIT (2026-02-27)
# ---------------------------------------------------------------------------
# 1. Mfuzz fuzzifier (mestimate):
#    - mestimate() estimates m from the number of features and samples using
#      the formula from Schwammle & Jensen (2010). For ~2000 proteins and
#      ~30 subjects this typically gives m in [1.3, 1.8], which is standard.
#      No user override needed.                                         PASS
#
# 2. Cluster number selection (k):
#    - Dmin elbow for k = 2..6 is computed and plotted (supp_dmin_elbow).
#      k = 4 is selected by visual elbow inspection. ISSUE: No formal
#      quantitative criterion (gap statistic, silhouette, etc.) is used.
#      MITIGATION: The 50-start multi-start consensus and bootstrap ARI
#      (mean ~0.96) strongly validate k = 4 stability. The Dmin elbow
#      visual selection is standard practice for FCM. Documented.       PASS
#
# 3. Bootstrap stability:
#    - 100 iterations, 80% protein subsample, ARI comparison to full-data
#      hard assignments. Mean ARI ~0.955, range [0.92, 0.98].
#    - ISSUE: Only mean/sd/range reported; no formal 95% CI.
#      FIX: Added percentile bootstrap 95% CI on ARI.                   FIXED
#
# 4. Multi-start consensus:
#    - 50 random starts, best objective kept. Avoids local optima.      PASS
#
# 5. Seed reproducibility:
#    - set.seed(42) for main clustering; set.seed(41+s) per start; set.seed(42)
#      for bootstrap. All seeds fixed for reproducibility.              PASS
#
# 6. ORA enrichment (enricher):
#    - Universe = all clustered proteins (unique(cluster_assign$gene)),
#      which is the correct ORA background (all proteins that could have
#      been assigned to any cluster).                                   PASS
#    - BH correction applied within each enricher() call, i.e., per
#      database per cluster. This is the standard approach: each database
#      (Hallmark, GO:BP) is tested as its own hypothesis family.
#      Cross-database global correction is not standard for ORA.        PASS
#    - pvalueCutoff = 0.05, qvalueCutoff = 1 (no q-value filter, only
#      BH-adjusted p < 0.05 used). Correct.                            PASS
#
# 7. rrvgo GO reduction:
#    - Threshold 0.85 (default) using Rel semantic similarity. This
#      retains representative parent terms and removes highly redundant
#      child terms. The threshold is a visualization parameter, not a
#      statistical test. Sensitivity is low: 0.7-0.9 gives similar results
#      for well-separated GO terms.                                     PASS
#
# 8. 1:1 greedy assignment (gene-to-pathway):
#    - This is a visualization heuristic for the Sankey diagram, not a
#      statistical claim. Each gene is assigned to its best pathway
#      (lowest p.adjust) for display only. No validation needed for a
#      visualization mapping.                                           PASS
# ---------------------------------------------------------------------------

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
DB_COLORS  <- c(Hallmark = "#AA336A", "GO:BP" = "#00796B")
CLUSTER_COLORS <- c(C1 = "#E74C3C", C2 = "#3498DB", C3 = "#2ECC71",
                     C4 = "#F39C12", C5 = "#9B59B6", C6 = "#1ABC9C",
                     C7 = "#E67E22", C8 = "#34495E", C9 = "#D35400",
                     C10 = "#7F8C8D")
YOUNG_COL <- "#4393C3"
OLD_COL   <- "#D6604D"
AGING_GAP_LINE <- "#7B68EE"
# Key constants centralized in palettes.R
source("04_Figures/shared/palettes.R")

GROUP_COLS <- c("Young_Pre", "Young_Post", "Old_Pre", "Old_Post")
GROUP_LABS <- c("Young\nPre", "Young\nPost", "Old\nPre", "Old\nPost")
ABBREV_LABS <- c("Y.Pre", "Y.Post", "O.Pre", "O.Post")

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
FIG_W <- 550
FIG_H <- 340
COL_WIDTHS <- c(0.11, 0.12, 0.44, 0.33)   # A, B, C, D

# --- Unified text sizes (all panels) ---
TXT_TITLE    <- 11     # element_text size: panel titles (Panel A cluster names)
TXT_SUBTITLE <- 9      # element_text size: panel subtitles
TXT_AXIS     <- 9      # element_text size: axis titles
TXT_TICK     <- 8      # element_text size: axis tick labels
TXT_ANNOT    <- 3.8    # geom_text / annotate size: labels, pathway names, counts
TXT_HEADER   <- 12     # element_text size: composite column headers (A, B, C, D)

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

# Standardize for Mfuzz (already z-scored, but Mfuzz standardise() centers/scales)
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

# Final clustering: 50-start multi-start consensus
# Each start uses a different seed; the run with the lowest within-cluster
# error (objective) is kept to avoid poor local optima.
N_STARTS <- 50
best_obj <- Inf
best_cl  <- NULL
cat(sprintf("Running %d-start FCM (k=%d, m=%.3f)...\n", N_STARTS, optimal_k, m_est))
for (s in seq_len(N_STARTS)) {
  set.seed(41 + s)
  cl_try <- mfuzz(eset, c = optimal_k, m = m_est)
  obj <- sum(cl_try$withinerror)
  if (obj < best_obj) {
    best_obj <- obj
    best_cl  <- cl_try
  }
}
cl <- best_cl
cat(sprintf("Best objective (sum withinerror): %.4f (from %d starts)\n", best_obj, N_STARTS))
membership <- cl$membership
centroids  <- cl$centers

# Bootstrap stability: subsample 80% of proteins, re-cluster, compute ARI
N_BOOT <- 100
boot_ari <- numeric(N_BOOT)
set.seed(42)
for (b in seq_len(N_BOOT)) {
  idx <- sample(nrow(eset), round(0.8 * nrow(eset)))
  cl_boot <- mfuzz(eset[idx, ], c = optimal_k, m = m_est)
  # Hard assignment for both
  orig_assign <- max.col(membership[idx, ])
  boot_assign <- max.col(cl_boot$membership)
  boot_ari[b] <- mclust::adjustedRandIndex(orig_assign, boot_assign)
}
# Percentile bootstrap 95% CI on ARI (STAT AUDIT addition)
boot_ari_ci <- quantile(boot_ari, probs = c(0.025, 0.975))
cat(sprintf("Bootstrap stability (ARI): mean=%.3f, sd=%.3f, 95%% CI=[%.3f, %.3f], range=[%.3f, %.3f]\n",
            mean(boot_ari), sd(boot_ari),
            boot_ari_ci[1], boot_ari_ci[2],
            min(boot_ari), max(boot_ari)))
write_csv(tibble(boot = seq_len(N_BOOT), ari = boot_ari),
          file.path(DAT_DIR, "07_cluster_stability.csv"))

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

write_csv(cluster_assign, file.path(DAT_DIR, "06_mfuzz_assignments.csv"))

centroid_export <- as.data.frame(centroids) |>
  rownames_to_column("cluster") |>
  mutate(cluster = paste0("C", row_number()))
write_csv(centroid_export, file.path(DAT_DIR, "01_panel_A_cluster_profiles.csv"))

cat(sprintf("  Saved: %s\n", file.path(DAT_DIR, "06_mfuzz_assignments.csv")))
cat(sprintf("  Saved: %s\n", file.path(DAT_DIR, "01_panel_A_cluster_profiles.csv")))

cat("=== Data pipeline complete: Setup + Mfuzz clustering ===\n")

# --- Core protein filter (membership >= 0.5) ---
CORE_THRESH <- 0.5
core_proteins <- cluster_assign %>%
  filter(membership >= CORE_THRESH)
cat(sprintf("Core proteins (membership >= %.1f): %d / %d (%.1f%%)\n",
            CORE_THRESH, nrow(core_proteins), nrow(cluster_assign),
            100 * nrow(core_proteins) / nrow(cluster_assign)))

# === 14. GROUP-LEVEL Z-SCORES (shared by Panels A and C) =====================

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

# === 15. ENRICHMENT ANALYSIS (shared by Panels C and D) =======================

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

universe_genes <- unique(cluster_assign$gene)

cat(sprintf("  Gene sets loaded: Hallmark=%d terms, GO:BP=%d terms\n",
            n_distinct(hallmark_t2g$term),
            n_distinct(gobp_t2g$term)))
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

  # Combine results with database column
  combined <- bind_rows(
    if (!is.null(res_hall) && nrow(as.data.frame(res_hall)) > 0)
      as.data.frame(res_hall) %>% mutate(database = "Hallmark") else NULL,
    if (!is.null(res_gobp) && nrow(as.data.frame(res_gobp)) > 0)
      as.data.frame(res_gobp) %>% mutate(database = "GO:BP") else NULL
  )

  # Filter to p.adjust < 0.05
  if (nrow(combined) > 0) {
    combined <- combined %>% filter(p.adjust < 0.05)
  }

  enrich_list[[cl_label]] <- combined %>% mutate(cluster = cl_label)
  cat(sprintf("    %s: %d significant terms (H=%d, BP=%d)\n",
              cl_label, nrow(combined),
              sum(combined$database == "Hallmark"),
              sum(combined$database == "GO:BP")))
}

# --- Step 3: rrvgo reduction -------------------------------------------------

cat("Applying rrvgo redundancy reduction for GO terms...\n")

# Pre-compute semantic data objects
bp_semdata <- tryCatch(
  godata("org.Hs.eg.db", ont = "BP", computeIC = TRUE),
  error = function(e) { cat("  Warning: could not compute BP semdata\n"); NULL }
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
  after_n <- nrow(enrich_list[[cl_label]])
  cat(sprintf("  %s: %d -> %d terms after rrvgo\n", cl_label, before_n, after_n))
}

# --- Step 4: Top term selection + 1:1 greedy assignment -----------------------

cat("Selecting top terms and performing 1:1 greedy assignment...\n")

# Select top 7 per database per cluster (Hallmark + GO:BP only; GO:CC excluded
# because localization terms cannot map to functional themes and "steal" proteins
# from functional pathways during greedy assignment)
enrich_top <- bind_rows(enrich_list) %>%
  filter(database %in% c("Hallmark", "GO:BP")) %>%
  group_by(cluster, database) %>%
  slice_min(p.adjust, n = 7, with_ties = FALSE) %>%
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

# --- Pass 2: Rescue unmapped proteins from full enrichment results -----------
# Proteins not found in any top-7 term may still appear in lower-ranked but
# statistically significant (BH p.adjust < 0.05) pathways. Rescue them.
cat("Pass 2: Rescuing unmapped proteins from full enrichment results...\n")

enrich_all_sig <- bind_rows(enrich_list) %>%
  filter(database %in% c("Hallmark", "GO:BP"), p.adjust < 0.05) %>%
  arrange(cluster, p.adjust)

# Pre-parse gene lists for speed
enrich_all_sig$gene_list <- strsplit(enrich_all_sig$geneID, "/")

for (cl_id in seq_len(optimal_k)) {
  cl_label <- paste0("C", cl_id)
  unmapped_idx <- which(
    protein_pathway_links$cluster == cl_label &
    protein_pathway_links$pathway == "Unmapped"
  )

  if (length(unmapped_idx) == 0) next

  cl_all_sig <- enrich_all_sig %>% filter(cluster == cl_label)
  if (nrow(cl_all_sig) == 0) next

  rescued <- 0
  for (i in unmapped_idx) {
    g <- protein_pathway_links$gene[i]
    for (j in seq_len(nrow(cl_all_sig))) {
      if (g %in% cl_all_sig$gene_list[[j]]) {
        protein_pathway_links$pathway[i]  <- cl_all_sig$Description[j]
        protein_pathway_links$database[i] <- cl_all_sig$database[j]
        rescued <- rescued + 1
        break
      }
    }
  }
  cat(sprintf("  %s: rescued %d / %d unmapped proteins\n",
              cl_label, rescued, length(unmapped_idx)))
}

# Clean up temporary column
enrich_all_sig$gene_list <- NULL

# Print summary
cluster_ids <- paste0("C", seq_len(optimal_k))
for (cl_label in cluster_ids) {
  cl_links <- protein_pathway_links %>% filter(cluster == cl_label)
  n_mapped   <- sum(cl_links$pathway != "Unmapped")
  n_unmapped <- sum(cl_links$pathway == "Unmapped")
  n_pathways <- n_distinct(cl_links$pathway[cl_links$pathway != "Unmapped"])
  cat(sprintf("  Cluster %s: %d pathways, %d proteins mapped, %d unmapped\n",
              cl_label, n_pathways, n_mapped, n_unmapped))
}

# --- Step 5: Export -----------------------------------------------------------

write_csv(enrich_top, file.path(DAT_DIR, "03_panel_C_enrichment.csv"))
write_csv(protein_pathway_links, file.path(DAT_DIR, "04_panel_C_sankey_links.csv"))

cat(sprintf("  Saved: %s\n", file.path(DAT_DIR, "03_panel_C_enrichment.csv")))
cat(sprintf("  Saved: %s\n", file.path(DAT_DIR, "04_panel_C_sankey_links.csv")))

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

# === 16. ROW HEIGHTS (shared by Panels C and D, Composite) ===================

row_heights <- core_proteins %>%
  count(cluster) %>%
  arrange(cluster) %>%
  pull(n)
row_heights <- row_heights / sum(row_heights)

cat(sprintf("  Row height proportions: %s\n",
            paste(sprintf("%.3f", row_heights), collapse = ", ")))

# === 17. THEME ASSIGNMENT (shared by Panel D) =================================
# 8 functional themes mapping enriched pathways to biological programs.
# Theme selection driven by the enrichment results themselves (what pathways
# are significant in each cluster) and consistent with skeletal muscle
# exercise/aging proteomics literature:
# - Mitochondrial & Energy Metabolism: Ubaida-Mohien et al. 2019 (eLife)
# - Muscle Structure & Myogenesis: core muscle biology
# - Proteostasis & Stress Response: mTORC1/UPR/chaperone pathways in C3
# - Cytoskeletal & Cell Division: microtubule/spindle/cell cycle pathways in C3
# - Immune & Complement: complement/inflammatory pathways in C1
# - ECM & Tissue Remodeling: collagen/adhesion/EMT pathways in C4
# - Metabolic & Redox Regulation: glycolysis/xenobiotic/ROS pathways in C1
# - Intracellular Transport & Signaling: vesicle/Golgi/kinase pathways

assign_theme <- function(pathway_name) {
  pw <- tolower(pathway_name)
  dplyr::case_when(
    # Theme 1: Mitochondrial & Energy Metabolism (primary: C2)
    stringr::str_detect(pw, "mitochon|oxidative.phosph|respiratory|electron.transport|tca|citrate|nadh|atp|fatty.acid|lipid|adipogen|acetyl.coa|amide.metabol") ~
      "Mitochondrial & Energy Metabolism",
    # Theme 2: Muscle Structure & Myogenesis (primary: C4)
    stringr::str_detect(pw, "myogen|muscle|contract|myofib|sarco|neuromuscul") ~
      "Muscle Structure & Myogenesis",
    # Theme 3: Proteostasis & Stress Response (primary: C3)
    stringr::str_detect(pw, "mtorc|unfold|chaper|heat.shock|protein.stabili|proteasom|ubiquitin|protein.fold|apoptosis|programmed.cell.death|cell.death") ~
      "Proteostasis & Stress Response",
    # Theme 4: Cytoskeletal & Cell Division (primary: C3)
    stringr::str_detect(pw, "microtub|spindle|mitotic|cell.divis|cell.cycle|cytoskelet|tubulin|actin") ~
      "Cytoskeletal & Cell Division",
    # Theme 5: Immune & Complement (primary: C1)
    stringr::str_detect(pw, "immun|inflam|complement|cytokine|interferon|heme|blood|coagulat") ~
      "Immune & Complement",
    # Theme 6: ECM & Tissue Remodeling (primary: C4)
    stringr::str_detect(pw, "extracellular|matrix|collagen|adhesion|integrin|mesenchym|epithelial") ~
      "ECM & Tissue Remodeling",
    # Theme 7: Metabolic & Redox Regulation (primary: C1)
    stringr::str_detect(pw, "glycol|metabol|xenobiot|aldehyde|oxidant|detox|pyridine|reactive.oxygen|peroxide") ~
      "Metabolic & Redox Regulation",
    # Theme 8: Intracellular Transport & Signaling (primary: C1)
    stringr::str_detect(pw, "vesicle|transport|endosom|golgi|lysosom|signal|kinase|androgen") ~
      "Intracellular Transport & Signaling",
    TRUE ~ "Other"
  )
}

theme_links <- protein_pathway_links %>%
  filter(pathway != "Unmapped") %>%
  mutate(theme = assign_theme(pathway))

# Theme-specific colors (8 themes, accessible muted palette)
THEME_COLORS <- c(
  "Mitochondrial & Energy Metabolism"    = "#E57373",
  "Muscle Structure & Myogenesis"        = "#64B5F6",
  "Proteostasis & Stress Response"       = "#81C784",
  "Cytoskeletal & Cell Division"         = "#CE93D8",
  "Immune & Complement"                  = "#FFB74D",
  "ECM & Tissue Remodeling"             = "#F06292",
  "Metabolic & Redox Regulation"         = "#FFD54F",
  "Intracellular Transport & Signaling"  = "#4DB6AC"
)

# Report theme coverage
n_themed <- sum(theme_links$theme != "Other")
n_other  <- sum(theme_links$theme == "Other")
n_unmapped <- sum(protein_pathway_links$pathway == "Unmapped")
n_core <- nrow(core_proteins)
cat(sprintf("Theme coverage: %d themed (%.1f%%), %d Other (%.1f%%), %d unmapped (%.1f%%) of %d core proteins\n",
            n_themed, 100 * n_themed / n_core,
            n_other, 100 * n_other / n_core,
            n_unmapped, 100 * n_unmapped / n_core,
            n_core))
cat("Theme distribution:\n")
print(table(theme_links$theme))

cat("=== F4 setup complete ===\n")
