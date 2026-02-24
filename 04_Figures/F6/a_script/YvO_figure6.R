################################################################################
#   Figure 6: Supervised Integration — sPLS & DIABLO
#
#   Study: YvO — DIA-MS Proteomics of Skeletal Muscle
#   Design: 2x2 mixed (Age x Time) with repeated measures
#
#   Central question: Can baseline proteome predict training-induced changes,
#   and does multi-block integration separate age groups?
#
#   Panels:
#     A — sPLS Sample Score Plot (Component 1 x 2)
#     B — sPLS Loading Plot (Top 20 Baseline Predictors)
#     C — sPLS Feature Stability (CV Selection Frequency)
#     D — DIABLO Arrow/Sample Plot
#     E — DIABLO Circos Plot
#     F — Multi-Method Convergence (UpSet Plot)
#
#   Inputs:
#     - 02_Imputation/c_data/01_imputed.csv
#     - 02_Imputation/c_data/01_DAList_imputed.rds
#     - 03_DEP/c_data/combined_results.csv
#     - 05_WGCNA/c_data/wgcna_hub_proteins.csv
#     - 05_WGCNA/c_data/wgcna_module_assignments.csv
#     - 04_Figures/F4/c_data/mfuzz_assignments.csv
################################################################################

# === 1. PACKAGES ==============================================================

suppressPackageStartupMessages({
  library(mixOmics)
  library(tidyverse)
  library(patchwork)
  library(UpSetR)
  library(ggrepel)
  library(scales)
  library(grid)
  library(png)
})

# === 2. SEED ==================================================================

set.seed(42)

# === 3. PATH RESOLUTION =======================================================

setwd("/Users/dtl0018/Desktop/A_Proteomics_Analysis/A_YvO_2025")
RPT_DIR <- "04_Figures/F6/b_reports"
DAT_DIR <- "04_Figures/F6/c_data"
SUP_DIR <- file.path(RPT_DIR, "supplementary")
dir.create(RPT_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(DAT_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(SUP_DIR, recursive = TRUE, showWarnings = FALSE)

# === 5. FILE PATH CONSTANTS ===================================================

IMPUTED_FILE <- "02_Imputation/c_data/01_imputed.csv"
DALIST_RDS   <- "02_Imputation/c_data/01_DAList_imputed.rds"
DEP_FILE     <- "03_DEP/c_data/combined_results.csv"
WGCNA_DIR    <- "05_WGCNA/c_data"
MOD_ASSIGN   <- file.path(WGCNA_DIR, "wgcna_module_assignments.csv")
HUB_FILE     <- file.path(WGCNA_DIR, "wgcna_hub_proteins.csv")
FCM_FILE     <- "04_Figures/F4/c_data/mfuzz_assignments.csv"

# === 6. CANONICAL CONSTANTS ===================================================

# Palettes — canonical (must match all figures F1-F6)
CONTRAST_COLORS <- c(Aging          = "#4CAF50",
                     Training_Young = "#E05A4E",
                     Training_Old   = "#5DA5DA",
                     Interaction    = "#9B7FBF")
AGE_COLORS  <- c(Young = "#4393C3", Old = "#D6604D")
DIR_COLORS  <- c(Up = "#D6604D", Down = "#4393C3")
GROUP_FILL  <- c(
  Young_Pre  = scales::alpha("#4393C3", 0.5),
  Young_Post = "#4393C3",
  Old_Pre    = scales::alpha("#D6604D", 0.5),
  Old_Post   = "#D6604D"
)
SIG_COLORS  <- c(Interaction      = "#9B7FBF",
                 `Sig Both`       = "#2E7D32",
                 `Sig Young only` = "#E05A4E",
                 `Sig Old only`   = "#5DA5DA",
                 NS               = "grey55")
KEY_TEXT  <- 2.2
KEY_TITLE <- 2.3

THEME_PUB <- theme_bw(base_size = 8) +
  theme(plot.title       = element_text(face = "bold", size = 9),
        plot.subtitle    = element_text(size = 6.5, color = "grey30", face = "italic"),
        strip.background = element_blank(),
        strip.text       = element_text(face = "bold", size = 6.5),
        legend.key.size  = unit(3, "mm"))

# === 7. DATA LOADING ==========================================================

cat("Loading data...\n")

# Imputed expression matrix
imp_df <- read_csv(IMPUTED_FILE, show_col_types = FALSE)
ann_cols   <- c("uniprot_id", "protein", "gene", "description")
samp_names <- setdiff(names(imp_df), ann_cols)
mat <- as.matrix(imp_df[, samp_names])
rownames(mat) <- imp_df$uniprot_id

# Gene lookup: uniprot -> gene
gene_lookup <- setNames(imp_df$gene, imp_df$uniprot_id)
# Handle NA/empty gene names: use uniprot_id as fallback
na_mask <- is.na(gene_lookup) | gene_lookup == ""
gene_lookup[na_mask] <- names(gene_lookup)[na_mask]
# Handle duplicate gene names
dup_counts <- table(gene_lookup)
dup_genes  <- names(dup_counts[dup_counts > 1])
if (length(dup_genes) > 0) {
  for (dg in dup_genes) {
    idx <- which(gene_lookup == dg)
    for (j in seq_along(idx)[-1]) {
      gene_lookup[idx[j]] <- paste0(dg, "_", j)
    }
  }
}

cat(sprintf("  Imputed matrix: %d proteins x %d samples\n", nrow(mat), ncol(mat)))

# DAList (for phenotype metadata)
dal <- readRDS(DALIST_RDS)
dal_meta <- as.data.frame(dal$metadata)

# DEP results
dep_df <- read_csv(DEP_FILE, show_col_types = FALSE)
cat(sprintf("  DEP results: %d proteins\n", nrow(dep_df)))

# WGCNA module assignments
mod_df <- read_csv(MOD_ASSIGN, show_col_types = FALSE)
cat(sprintf("  WGCNA modules: %d proteins\n", nrow(mod_df)))

# WGCNA hub proteins
hub_df <- read_csv(HUB_FILE, show_col_types = FALSE)
cat(sprintf("  WGCNA hubs: %d proteins\n", nrow(hub_df)))

# FCM assignments
fcm_df <- read_csv(FCM_FILE, show_col_types = FALSE)
cat(sprintf("  FCM clusters: %d proteins, %d clusters\n",
            nrow(fcm_df), length(unique(fcm_df$cluster))))

# === 8. PARSE SAMPLE METADATA =================================================

meta <- tibble(sample_id = samp_names) |>
  mutate(
    age     = if_else(str_detect(sample_id, "^(O_|OP_)"), "Old", "Young"),
    time    = if_else(str_detect(sample_id, "_Post$"), "Post", "Pre"),
    subject = str_remove(sample_id, "_(Pre|Post)$"),
    group   = paste(age, time, sep = "_")
  )

cat(sprintf("  Samples: %d total (%d Young, %d Old)\n",
            nrow(meta), sum(meta$age == "Young"), sum(meta$age == "Old")))

# === 9. COMPUTE DELTA MATRIX (Post - Pre) PER SUBJECT ========================

cat("Computing delta matrix (Post - Pre)...\n")

pre_subjects  <- meta$subject[meta$time == "Pre"]
post_subjects <- meta$subject[meta$time == "Post"]
paired_subjects <- intersect(pre_subjects, post_subjects)

young_subjects <- sort(paired_subjects[str_detect(paired_subjects, "^(Y_|YP_)")])
old_subjects   <- sort(paired_subjects[str_detect(paired_subjects, "^(O_|OP_)")])
ordered_subjects <- c(young_subjects, old_subjects)

cat(sprintf("  Paired subjects: %d (%d Young, %d Old)\n",
            length(ordered_subjects), length(young_subjects), length(old_subjects)))

# Delta matrix: Post - Pre for each subject (proteins x subjects)
delta_mat <- matrix(NA_real_, nrow = nrow(mat), ncol = length(ordered_subjects),
                    dimnames = list(rownames(mat), ordered_subjects))
for (subj in ordered_subjects) {
  pre_col  <- paste0(subj, "_Pre")
  post_col <- paste0(subj, "_Post")
  delta_mat[, subj] <- mat[, post_col] - mat[, pre_col]
}

cat(sprintf("  Delta matrix: %d proteins x %d subjects\n",
            nrow(delta_mat), ncol(delta_mat)))

# === 10. PREPARE PHENOTYPE DELTAS PER SUBJECT =================================

cat("Computing phenotype deltas...\n")

# Phenotype columns of interest
pheno_cols <- c("VL_thick_cm", "DXA_LBM_kg", "Type_I_fCSA", "Type_II_fCSA")

# Ensure numeric
for (pc in pheno_cols) {
  if (pc %in% names(dal_meta)) {
    vals <- dal_meta[[pc]]
    if (is.character(vals)) {
      vals[vals %in% c("NA", "na", "")] <- NA
      vals <- as.numeric(vals)
    }
    dal_meta[[pc]] <- as.numeric(vals)
  }
}

# Build per-subject phenotype delta
pheno_delta <- tibble(subject = ordered_subjects)

for (pc in pheno_cols) {
  if (pc %in% names(dal_meta)) {
    pre_vals  <- dal_meta[[pc]][match(paste0(ordered_subjects, "_Pre"),  dal_meta$Col_ID)]
    post_vals <- dal_meta[[pc]][match(paste0(ordered_subjects, "_Post"), dal_meta$Col_ID)]
    pheno_delta[[paste0("delta_", pc)]] <- post_vals - pre_vals
  }
}

# Age group per subject
subj_age <- meta |>
  filter(time == "Pre") |>
  dplyr::select(subject, age) |>
  distinct()
pheno_delta <- pheno_delta |> left_join(subj_age, by = "subject")

cat(sprintf("  Phenotype deltas: %d subjects\n", nrow(pheno_delta)))
for (pc in pheno_cols) {
  dc <- paste0("delta_", pc)
  if (dc %in% names(pheno_delta)) {
    n_ok <- sum(complete.cases(pheno_delta[[dc]]))
    cat(sprintf("    %s: %d non-NA\n", dc, n_ok))
  }
}

# === 11. PREPARE sPLS BLOCKS ==================================================

cat("\n=== Preparing sPLS blocks ===\n")

# X block: baseline (Pre) expression, subjects x proteins
pre_samples <- paste0(ordered_subjects, "_Pre")
X_baseline  <- t(mat[, pre_samples])
rownames(X_baseline) <- ordered_subjects
# Use gene symbols as column names
colnames(X_baseline) <- gene_lookup[colnames(X_baseline)]

# Center and scale
X_baseline <- scale(X_baseline)

cat(sprintf("  X (baseline): %d subjects x %d proteins\n",
            nrow(X_baseline), ncol(X_baseline)))

# Y block: training delta for interaction proteins
# Interaction proteins: p < 0.05 from Interaction contrast
interaction_prots <- dep_df |>
  filter(P.Value_Interaction < 0.05) |>
  pull(uniprot_id)

cat(sprintf("  Interaction proteins (p < 0.05): %d\n", length(interaction_prots)))

# If fewer than 20, relax to p < 0.1
if (length(interaction_prots) < 20) {
  interaction_prots <- dep_df |>
    filter(P.Value_Interaction < 0.1) |>
    pull(uniprot_id)
  cat(sprintf("  Relaxed to p < 0.1: %d proteins\n", length(interaction_prots)))
}

# Filter delta_mat to interaction proteins
interaction_prots_avail <- intersect(interaction_prots, rownames(delta_mat))
Y_delta <- t(delta_mat[interaction_prots_avail, ])
colnames(Y_delta) <- gene_lookup[interaction_prots_avail]
Y_delta <- scale(Y_delta)

cat(sprintf("  Y (delta, interaction proteins): %d subjects x %d proteins\n",
            nrow(Y_delta), ncol(Y_delta)))

# Age class factor
Y_class <- factor(
  pheno_delta$age[match(ordered_subjects, pheno_delta$subject)],
  levels = c("Young", "Old")
)

cat(sprintf("  Y_class: %s\n", paste(table(Y_class), collapse = " / ")))

# === 12. sPLS ANALYSIS ========================================================

cat("\n=== sPLS Analysis ===\n")

# Step 1: Determine ncomp via non-sparse PLS
cat("  Fitting initial PLS (ncomp = 5)...\n")
max_ncomp <- min(5, nrow(X_baseline) - 1, ncol(Y_delta))
pls_init <- pls(X_baseline, Y_delta, ncomp = max_ncomp)

cat("  Running perf() for ncomp selection (5-fold CV x 20 repeats)...\n")
set.seed(42)
pls_perf <- perf(pls_init, validation = "Mfold", folds = 5, nrepeat = 20,
                 progressBar = FALSE)

# Extract Q2 per component from perf() output
# mixOmics stores Q2 in $measures$Q2.total (list with $values + $summary)
q2_avg <- NULL
q2_label <- "Q2"

# Try Q2.total first
q2_raw <- pls_perf$measures$Q2.total %||% pls_perf$Q2.total %||% pls_perf$Q2
if (!is.null(q2_raw)) {
  if (is.list(q2_raw) && "summary" %in% names(q2_raw)) {
    # mixOmics format: list with $summary dataframe (feature, comp, mean, sd)
    q2_summary <- q2_raw$summary
    q2_avg <- q2_summary$mean
    names(q2_avg) <- q2_summary$comp
    cat("  Q2 per component (mean):", paste(round(q2_avg, 4), collapse = ", "), "\n")
  } else if (is.matrix(q2_raw)) {
    q2_avg <- rowMeans(q2_raw, na.rm = TRUE)
  } else if (is.numeric(q2_raw)) {
    q2_avg <- q2_raw
  }
}

# Fallback: use cor.tpred if Q2 unavailable
if (is.null(q2_avg)) {
  cor_raw <- pls_perf$measures$cor.tpred %||% pls_perf$cor.tpred
  if (!is.null(cor_raw) && is.list(cor_raw) && "summary" %in% names(cor_raw)) {
    q2_avg <- cor_raw$summary$mean
    q2_label <- "cor"
    cat("  cor.tpred per component:", paste(round(q2_avg, 4), collapse = ", "), "\n")
  } else {
    cat("  No Q2 or cor metric found; defaulting to 1 component.\n")
    q2_avg <- rep(-1, max_ncomp)
  }
}

# Select ncomp: use components with Q2 > 0.0975 (standard threshold)
# If all Q2 < 0 (poor prediction), default to 1 component
ncomp_selected <- max(1, sum(q2_avg > 0.0975))
ncomp_selected <- min(ncomp_selected, 3)  # cap at 3

cat(sprintf("  Selected ncomp = %d\n", ncomp_selected))

# Step 2: Tune sparsity
cat("  Tuning sPLS sparsity...\n")
keepX_grid <- c(10, 25, 50, 100)
keepY_grid <- c(5, 10, 25)

# Filter grids to not exceed available features
keepX_grid <- keepX_grid[keepX_grid <= ncol(X_baseline)]
keepY_grid <- keepY_grid[keepY_grid <= ncol(Y_delta)]

set.seed(42)
tune_spls <- tune.spls(X_baseline, Y_delta,
                        ncomp = ncomp_selected,
                        test.keepX = keepX_grid,
                        test.keepY = keepY_grid,
                        validation = "Mfold",
                        folds = 5, nrepeat = 10,
                        measure = "cor",
                        progressBar = FALSE)

# Extract optimal keepX and keepY
optimal_keepX <- tune_spls$choice.keepX
optimal_keepY <- tune_spls$choice.keepY
cat(sprintf("  Optimal keepX: %s\n", paste(optimal_keepX, collapse = ", ")))
cat(sprintf("  Optimal keepY: %s\n", paste(optimal_keepY, collapse = ", ")))

# Step 3: Fit final sPLS
cat("  Fitting final sPLS model...\n")
spls_final <- spls(X_baseline, Y_delta,
                    ncomp = ncomp_selected,
                    keepX = optimal_keepX,
                    keepY = optimal_keepY)

cat(sprintf("  Final sPLS: ncomp = %d, keepX = %s, keepY = %s\n",
            ncomp_selected,
            paste(optimal_keepX, collapse = ","),
            paste(optimal_keepY, collapse = ",")))

# Step 4: Feature stability via custom CV
cat("  Running feature stability CV (50 iterations x 5 folds)...\n")

n_iter   <- 50
n_folds  <- 5
n_subj   <- nrow(X_baseline)
feature_counts <- rep(0, ncol(X_baseline))
names(feature_counts) <- colnames(X_baseline)

set.seed(42)
for (iter in seq_len(n_iter)) {
  fold_ids <- sample(rep(seq_len(n_folds), length.out = n_subj))
  for (f in seq_len(n_folds)) {
    train_idx <- which(fold_ids != f)
    # Need at least 2 samples per fold
    if (length(train_idx) < 5) next
    spls_cv <- tryCatch(
      spls(X_baseline[train_idx, , drop = FALSE],
           Y_delta[train_idx, , drop = FALSE],
           ncomp = ncomp_selected,
           keepX = optimal_keepX,
           keepY = optimal_keepY),
      error = function(e) NULL
    )
    if (is.null(spls_cv)) next
    # Get selected features from Component 1
    sel_x <- selectVar(spls_cv, comp = 1)$X$name
    feature_counts[sel_x] <- feature_counts[sel_x] + 1
  }
}

# Convert to frequency
total_fits <- n_iter * n_folds
feature_freq <- feature_counts / total_fits

stable_features <- names(feature_freq[feature_freq > 0.70])
cat(sprintf("  Feature stability: %d features selected in >70%% of folds\n",
            length(stable_features)))
cat(sprintf("  Total features ever selected: %d\n", sum(feature_freq > 0)))

# === 13. DIABLO ANALYSIS =====================================================

cat("\n=== DIABLO Analysis ===\n")

# Prepare blocks for DIABLO
# Block 1: baseline expression (top 500 by variance)
var_baseline <- apply(X_baseline, 2, var, na.rm = TRUE)
top500_baseline <- names(sort(var_baseline, decreasing = TRUE))[1:min(500, ncol(X_baseline))]
X_diablo_base <- X_baseline[, top500_baseline]

# Block 2: delta expression (top 500 by variance)
delta_all_scaled <- scale(t(delta_mat))
rownames(delta_all_scaled) <- ordered_subjects
colnames(delta_all_scaled) <- gene_lookup[colnames(delta_all_scaled)]
var_delta <- apply(delta_all_scaled, 2, var, na.rm = TRUE)
top500_delta <- names(sort(var_delta, decreasing = TRUE))[1:min(500, ncol(delta_all_scaled))]
X_diablo_delta <- delta_all_scaled[, top500_delta]

# Block 3: phenotype deltas
pheno_block_cols <- c("delta_VL_thick_cm", "delta_DXA_LBM_kg",
                      "delta_Type_I_fCSA", "delta_Type_II_fCSA")
pheno_block_avail <- intersect(pheno_block_cols, names(pheno_delta))

# Build phenotype matrix for DIABLO subjects
pheno_mat_diablo <- as.matrix(pheno_delta[match(ordered_subjects, pheno_delta$subject),
                                           pheno_block_avail])
rownames(pheno_mat_diablo) <- ordered_subjects

# Filter to complete subjects across all blocks
complete_mask <- complete.cases(X_diablo_base) &
                 complete.cases(X_diablo_delta) &
                 complete.cases(pheno_mat_diablo)

n_complete <- sum(complete_mask)
cat(sprintf("  Complete subjects across all blocks: %d\n", n_complete))

# Check if we have enough subjects for DIABLO with phenotype block
run_diablo_pheno <- n_complete >= 15

if (run_diablo_pheno) {
  diablo_subjects <- ordered_subjects[complete_mask]
  diablo_blocks <- list(
    baseline  = X_diablo_base[diablo_subjects, , drop = FALSE],
    delta     = X_diablo_delta[diablo_subjects, , drop = FALSE],
    phenotype = scale(pheno_mat_diablo[diablo_subjects, , drop = FALSE])
  )
  Y_diablo <- Y_class[complete_mask]
  cat(sprintf("  DIABLO with phenotype block: %d subjects\n", length(diablo_subjects)))
} else {
  # Fall back to 2-block DIABLO (baseline + delta only)
  cat("  Not enough complete subjects for phenotype block; using 2-block DIABLO\n")
  diablo_subjects <- ordered_subjects
  diablo_blocks <- list(
    baseline = X_diablo_base,
    delta    = X_diablo_delta
  )
  Y_diablo <- Y_class
}

cat(sprintf("  DIABLO blocks: %s\n",
            paste(names(diablo_blocks), sapply(diablo_blocks, ncol), sep = "=", collapse = ", ")))

# Design matrix: full design (off-diagonal = 1)
n_blocks <- length(diablo_blocks)
design_mat <- matrix(1, n_blocks, n_blocks,
                     dimnames = list(names(diablo_blocks), names(diablo_blocks)))
diag(design_mat) <- 0

# Determine ncomp (cap at min block ncol to avoid dimension mismatch)
min_block_ncol <- min(sapply(diablo_blocks, ncol))
max_diablo_ncomp <- min(5, length(diablo_subjects) - 2, min_block_ncol)
cat(sprintf("  Fitting initial block.plsda (ncomp = %d)...\n", max_diablo_ncomp))

diablo_init <- tryCatch(
  block.plsda(diablo_blocks, Y_diablo,
              ncomp = max_diablo_ncomp,
              design = design_mat),
  error = function(e) {
    cat(sprintf("  block.plsda error: %s\n", conditionMessage(e)))
    NULL
  }
)

diablo_ncomp <- 2  # default

if (!is.null(diablo_init)) {
  cat("  Running perf() for DIABLO ncomp selection...\n")
  set.seed(42)
  diablo_perf <- tryCatch(
    perf(diablo_init, validation = "Mfold", folds = 5, nrepeat = 10,
         progressBar = FALSE),
    error = function(e) {
      cat(sprintf("  perf() error: %s\n", conditionMessage(e)))
      NULL
    }
  )

  if (!is.null(diablo_perf)) {
    # Select ncomp by BER (balanced error rate)
    ber_vals <- diablo_perf$WeightedVote.error.rate
    # Use centroids.dist if available
    if ("centroids.dist" %in% names(ber_vals)) {
      ber_mat <- ber_vals$centroids.dist
    } else {
      ber_mat <- ber_vals[[1]]
    }
    cat("  BER values:\n")
    print(ber_mat)
    # Pick ncomp with min overall BER
    # BER matrix: rows = classes + Overall.ER/BER, cols = components
    if (is.matrix(ber_mat)) {
      if ("Overall.BER" %in% rownames(ber_mat)) {
        overall_ber <- ber_mat["Overall.BER", ]
      } else if ("Overall.ER" %in% rownames(ber_mat)) {
        overall_ber <- ber_mat["Overall.ER", ]
      } else {
        # Last row is typically overall
        overall_ber <- ber_mat[nrow(ber_mat), ]
      }
      diablo_ncomp <- which.min(overall_ber)
      cat(sprintf("  Overall BER: %s\n", paste(round(overall_ber, 4), collapse = ", ")))
    }
    diablo_ncomp <- max(1, min(diablo_ncomp, 3))
    cat(sprintf("  Selected DIABLO ncomp = %d\n", diablo_ncomp))
  }
}

# Tune sparsity
cat("  Tuning DIABLO sparsity...\n")
keepX_diablo <- list()
for (bn in names(diablo_blocks)) {
  nc <- ncol(diablo_blocks[[bn]])
  if (bn == "phenotype") {
    keepX_diablo[[bn]] <- rep(min(nc, 4), diablo_ncomp)
  } else {
    keepX_diablo[[bn]] <- rep(min(25, nc), diablo_ncomp)
  }
}

# Try tune.block.splsda; fall back to fixed keepX if it fails
diablo_tuned <- tryCatch({
  test_keepX_diablo <- list()
  for (bn in names(diablo_blocks)) {
    nc <- ncol(diablo_blocks[[bn]])
    if (bn == "phenotype") {
      test_keepX_diablo[[bn]] <- c(2, 3, min(4, nc))
    } else {
      test_keepX_diablo[[bn]] <- c(10, 25, min(50, nc))
    }
  }

  set.seed(42)
  tune_res <- tune.block.splsda(
    diablo_blocks, Y_diablo,
    ncomp = diablo_ncomp,
    test.keepX = test_keepX_diablo,
    design = design_mat,
    validation = "Mfold",
    folds = 5, nrepeat = 10,
    progressBar = FALSE
  )
  cat("  DIABLO tuning complete.\n")
  cat(sprintf("  Tuned keepX: %s\n",
              paste(sapply(tune_res$choice.keepX, paste, collapse = ","),
                    collapse = " | ")))
  tune_res$choice.keepX
}, error = function(e) {
  cat(sprintf("  tune.block.splsda error: %s — using defaults\n", conditionMessage(e)))
  keepX_diablo
})

# Fit final DIABLO
cat("  Fitting final DIABLO model...\n")
diablo_final <- tryCatch(
  block.splsda(diablo_blocks, Y_diablo,
               ncomp = diablo_ncomp,
               keepX = diablo_tuned,
               design = design_mat),
  error = function(e) {
    cat(sprintf("  block.splsda error: %s\n", conditionMessage(e)))
    NULL
  }
)

if (!is.null(diablo_final)) {
  cat(sprintf("  Final DIABLO: ncomp = %d, %d blocks, %d subjects\n",
              diablo_ncomp, length(diablo_blocks), length(diablo_subjects)))
} else {
  cat("  WARNING: DIABLO fitting failed — panels D/E will use fallbacks.\n")
}

# === DIABLO Permutation Test ==================================================
message("Running DIABLO permutation test (50 permutations)...")

# Determine the observed BER from the diablo_perf object
observed_ber <- NA_real_
if (exists("diablo_perf") && !is.null(diablo_perf)) {
  ber_vals_perm <- diablo_perf$WeightedVote.error.rate
  if ("centroids.dist" %in% names(ber_vals_perm)) {
    ber_mat_perm <- ber_vals_perm$centroids.dist
  } else {
    ber_mat_perm <- ber_vals_perm[[1]]
  }
  if (is.matrix(ber_mat_perm)) {
    if ("Overall.BER" %in% rownames(ber_mat_perm)) {
      observed_ber <- min(ber_mat_perm["Overall.BER", ])
    } else if ("Overall.ER" %in% rownames(ber_mat_perm)) {
      observed_ber <- min(ber_mat_perm["Overall.ER", ])
    } else {
      observed_ber <- min(ber_mat_perm[nrow(ber_mat_perm), ])
    }
  }
}

perm_pval <- NA_real_
perm_ber_valid <- numeric(0)

if (!is.na(observed_ber) && !is.null(diablo_final)) {
  set.seed(42)
  n_perm <- 50
  perm_ber <- numeric(n_perm)

  for (i in seq_len(n_perm)) {
    if (i %% 10 == 0) cat(sprintf("  Permutation %d/%d\n", i, n_perm))
    Y_perm <- sample(Y_diablo)  # shuffle class labels
    perm_fit <- tryCatch(
      block.splsda(X = diablo_blocks, Y = Y_perm, ncomp = diablo_ncomp,
                   keepX = diablo_tuned, design = design_mat),
      error = function(e) NULL
    )
    if (!is.null(perm_fit)) {
      perm_perf <- tryCatch(
        perf(perm_fit, validation = "Mfold", folds = 5, nrepeat = 5,
             progressBar = FALSE),
        error = function(e) NULL
      )
      if (!is.null(perm_perf)) {
        perm_ber_vals <- perm_perf$WeightedVote.error.rate
        if ("centroids.dist" %in% names(perm_ber_vals)) {
          perm_ber_m <- perm_ber_vals$centroids.dist
        } else {
          perm_ber_m <- perm_ber_vals[[1]]
        }
        if (is.matrix(perm_ber_m)) {
          if ("Overall.BER" %in% rownames(perm_ber_m)) {
            perm_ber[i] <- min(perm_ber_m["Overall.BER", ])
          } else if ("Overall.ER" %in% rownames(perm_ber_m)) {
            perm_ber[i] <- min(perm_ber_m["Overall.ER", ])
          } else {
            perm_ber[i] <- min(perm_ber_m[nrow(perm_ber_m), ])
          }
        } else {
          perm_ber[i] <- NA
        }
      } else {
        perm_ber[i] <- NA
      }
    } else {
      perm_ber[i] <- NA
    }
  }

  perm_ber_valid <- perm_ber[!is.na(perm_ber)]
  perm_pval <- (sum(perm_ber_valid <= observed_ber) + 1) / (length(perm_ber_valid) + 1)
  cat(sprintf("DIABLO permutation p-value: %.4f (observed BER: %.3f, null median: %.3f)\n",
              perm_pval, observed_ber, median(perm_ber_valid)))

  # Export permutation results
  perm_results <- tibble(
    observed_ber = observed_ber,
    null_median_ber = median(perm_ber_valid),
    null_mean_ber = mean(perm_ber_valid),
    perm_pvalue = perm_pval,
    n_permutations = length(perm_ber_valid)
  )
  write_csv(perm_results, file.path(DAT_DIR, "diablo_permutation_test.csv"))
  write_csv(tibble(permutation = seq_along(perm_ber), ber = perm_ber),
            file.path(DAT_DIR, "diablo_permutation_null.csv"))
} else {
  cat("  Skipping permutation test: observed BER not available or DIABLO failed.\n")
}

# === 14. PANEL A — sPLS Sample Score Plot =====================================

cat("\n=== Panel A: sPLS Sample Score Plot ===\n")

# Extract scores
spls_scores_X <- as.data.frame(spls_final$variates$X)
spls_scores_Y <- as.data.frame(spls_final$variates$Y)

score_df <- tibble(
  subject = rownames(spls_scores_X),
  comp1_X = spls_scores_X[, 1],
  comp1_Y = spls_scores_Y[, 1],
  age     = pheno_delta$age[match(rownames(spls_scores_X), pheno_delta$subject)]
)

if (ncomp_selected >= 2) {
  score_df$comp2_X <- spls_scores_X[, 2]
  score_df$comp2_Y <- spls_scores_Y[, 2]
}

# Explained variance per component
if (!is.null(spls_final$prop_expl_var) && "X" %in% names(spls_final$prop_expl_var)) {
  expl_var <- spls_final$prop_expl_var$X[seq_len(ncomp_selected)] * 100
} else {
  # Fallback: reconstruction-based variance explained
  expl_var <- sapply(seq_len(ncomp_selected), function(i) {
    X_hat_i <- spls_scores_X[, i, drop = FALSE] %*%
               t(spls_final$loadings$X[, i, drop = FALSE])
    sum(X_hat_i^2) / sum(X_baseline^2, na.rm = TRUE) * 100
  })
}

# Q2/cor annotation for panels
q2_annot <- paste(
  sprintf("%s comp%d = %.3f", q2_label,
          seq_len(min(ncomp_selected, length(q2_avg))),
          q2_avg[seq_len(min(ncomp_selected, length(q2_avg)))]),
  collapse = "\n"
)

n_sel_X <- sum(spls_final$loadings$X[, 1] != 0)
n_sel_Y <- sum(spls_final$loadings$Y[, 1] != 0)
sel_annot <- sprintf("Selected: %d baseline, %d delta features", n_sel_X, n_sel_Y)

if (ncomp_selected >= 2) {
  # 2D scatter: Comp1 vs Comp2 (X variates)
  panel_A <- ggplot(score_df, aes(x = comp1_X, y = comp2_X, color = age)) +
    stat_ellipse(aes(group = age), level = 0.80, linetype = "dashed",
                 linewidth = 0.5, show.legend = FALSE) +
    geom_point(size = 2.5, alpha = 0.8, shape = 16) +
    scale_color_manual(values = AGE_COLORS, name = "Age") +
    annotate("text", x = Inf, y = -Inf, label = q2_annot,
             hjust = 1.05, vjust = -0.5, size = 2.2, fontface = "italic",
             color = "grey30") +
    annotate("text", x = -Inf, y = Inf, label = sel_annot,
             hjust = -0.05, vjust = 1.5, size = 2.2, fontface = "italic",
             color = "grey30") +
    labs(
      title    = "A  sPLS Sample Score Plot",
      subtitle = "Baseline proteome predicting training response",
      x        = sprintf("Component 1 (%.1f%% var.)", expl_var[1]),
      y        = sprintf("Component 2 (%.1f%% var.)", expl_var[2])
    ) +
    THEME_PUB +
    theme(legend.position = c(0.02, 0.02),
          legend.justification = c(0, 0),
          legend.background = element_rect(fill = alpha("white", 0.85), color = NA))
} else {
  # Fallback: boxplot of Comp1 scores by age
  wt <- wilcox.test(comp1_X ~ age, data = score_df)
  wt_label <- if (wt$p.value < 0.001) {
    sprintf("Wilcoxon p = %.1e", wt$p.value)
  } else {
    sprintf("Wilcoxon p = %.3f", wt$p.value)
  }

  panel_A <- ggplot(score_df, aes(x = age, y = comp1_X, fill = age)) +
    geom_boxplot(width = 0.5, outlier.shape = NA, alpha = 0.7) +
    geom_jitter(aes(color = age), width = 0.15, size = 2, alpha = 0.6, shape = 16) +
    scale_fill_manual(values = AGE_COLORS) +
    scale_color_manual(values = AGE_COLORS) +
    annotate("text", x = 1.5, y = max(score_df$comp1_X) * 1.05,
             label = wt_label, size = 2.5, fontface = "italic", color = "grey30") +
    annotate("text", x = 1.5, y = min(score_df$comp1_X) - 0.3,
             label = q2_annot, size = 2.2, fontface = "italic", color = "grey30") +
    labs(
      title    = "A  sPLS Sample Scores",
      subtitle = sprintf("%s | %s", sel_annot, q2_annot),
      x        = NULL,
      y        = "Component 1 Score"
    ) +
    THEME_PUB +
    theme(legend.position = "none")
}

cat("  Panel A built.\n")

# Save Panel A individual
ggsave(file.path(RPT_DIR, "panel_A_spls_scores.pdf"),
       plot = panel_A, width = 127, height = 127, units = "mm", device = pdf)

# Save Panel A source data
write_csv(score_df, file.path(DAT_DIR, "fig6_panel_A_spls_scores.csv"))
cat("  Saved: panel_A_spls_scores.pdf + fig6_panel_A_spls_scores.csv\n")

# === 15. PANEL B — sPLS Loading Plot (Top 20 Baseline Predictors) =============

cat("\n=== Panel B: sPLS Loading Plot ===\n")

# Extract Component 1 loadings for X block
loadings_X <- spls_final$loadings$X[, 1]
loadings_df <- tibble(
  gene    = names(loadings_X),
  loading = as.numeric(loadings_X)
) |>
  filter(loading != 0) |>
  arrange(desc(abs(loading)))

# Get top 20
top20 <- head(loadings_df, 20) |>
  mutate(gene = factor(gene, levels = rev(gene)))

# Annotate with WGCNA module
wgcna_lookup <- setNames(mod_df$module_color, mod_df$gene)
fcm_lookup   <- setNames(fcm_df$cluster, fcm_df$gene)

top20 <- top20 |>
  mutate(
    wgcna_module = wgcna_lookup[as.character(gene)],
    fcm_cluster  = fcm_lookup[as.character(gene)],
    # Mark stable features
    stable       = as.character(gene) %in% stable_features,
    stable_label = if_else(stable, "*", "")
  )

# Compute median expression in Young vs Old for color
pre_young <- meta |> filter(time == "Pre", age == "Young") |> pull(sample_id)
pre_old   <- meta |> filter(time == "Pre", age == "Old")   |> pull(sample_id)

# Map gene symbols back to uniprot IDs for expression lookup
gene_to_uniprot <- setNames(names(gene_lookup), gene_lookup)

top20 <- top20 |>
  rowwise() |>
  mutate(
    uid = gene_to_uniprot[as.character(gene)],
    median_young = if (!is.na(uid) && uid %in% rownames(mat))
      median(mat[uid, pre_young], na.rm = TRUE) else NA_real_,
    median_old = if (!is.na(uid) && uid %in% rownames(mat))
      median(mat[uid, pre_old], na.rm = TRUE) else NA_real_,
    fc_young_vs_old = median_young - median_old  # log2 scale
  ) |>
  ungroup()

panel_B <- ggplot(top20, aes(x = loading, y = gene)) +
  geom_col(aes(fill = fc_young_vs_old), width = 0.7) +
  geom_text(aes(label = stable_label, x = loading + sign(loading) * 0.002),
            size = 3, color = "#DAA520", hjust = if_else(top20$loading > 0, 0, 1)) +
  scale_fill_gradient2(
    low      = AGE_COLORS["Old"],
    mid      = "grey90",
    high     = AGE_COLORS["Young"],
    midpoint = 0,
    name     = "Young-Old\n(median log2)",
    guide    = guide_colorbar(barwidth = unit(2.5, "cm"), barheight = unit(0.3, "cm"),
                              title.position = "top")
  ) +
  geom_vline(xintercept = 0, color = "grey40", linewidth = 0.3) +
  labs(
    title    = "B  Top Baseline Predictors (sPLS Loadings)",
    subtitle = sprintf("Component 1 | * = stable (>70%% CV selection) | %d/%d stable",
                       sum(top20$stable), nrow(loadings_df)),
    x        = "sPLS Loading (Component 1)",
    y        = NULL
  ) +
  THEME_PUB +
  theme(legend.position  = "bottom",
        axis.text.y      = element_text(size = 6, face = "bold"))

cat("  Panel B built.\n")

# Save Panel B
ggsave(file.path(RPT_DIR, "panel_B_spls_loadings.pdf"),
       plot = panel_B, width = 127, height = 127, units = "mm", device = pdf)

write_csv(top20 |> dplyr::select(gene, loading, wgcna_module, fcm_cluster,
                                   stable, fc_young_vs_old),
          file.path(DAT_DIR, "fig6_panel_B_spls_loadings.csv"))
cat("  Saved: panel_B_spls_loadings.pdf + fig6_panel_B_spls_loadings.csv\n")

# === 16. PANEL C — Feature Stability (CV Selection Frequency) =================

cat("\n=== Panel C: Feature Stability ===\n")

# Top 40 features by CV selection frequency
stab_df <- tibble(
  gene = names(feature_freq),
  freq = as.numeric(feature_freq)
) |>
  filter(freq > 0) |>
  arrange(desc(freq)) |>
  head(40) |>
  mutate(
    gene   = factor(gene, levels = rev(gene)),
    stable = freq >= 0.70
  )

n_total_selected <- sum(feature_freq > 0)
n_stable <- sum(feature_freq >= 0.70)

panel_C <- ggplot(stab_df, aes(x = freq, y = gene, fill = stable)) +
  geom_col(width = 0.7) +
  geom_vline(xintercept = 0.70, linetype = "dashed", color = "grey30", linewidth = 0.4) +
  scale_fill_manual(values = c("TRUE" = AGE_COLORS["Young"], "FALSE" = "grey70"),
                    guide = "none") +
  scale_x_continuous(labels = percent_format(), limits = c(0, 1.02), expand = c(0, 0)) +
  annotate("text", x = 0.72, y = nrow(stab_df) * 0.95,
           label = "70% threshold", hjust = 0, size = 2.2,
           color = "grey30", fontface = "italic") +
  labs(
    title    = "C  Feature Stability (Cross-Validation)",
    subtitle = sprintf("%d features selected, %d exceed 70%% threshold | 50 iter x 5-fold",
                       n_total_selected, n_stable),
    x        = "Selection Frequency",
    y        = NULL
  ) +
  THEME_PUB +
  theme(axis.text.y = element_text(size = 5.5))

cat("  Panel C built.\n")

# Save Panel C
ggsave(file.path(RPT_DIR, "panel_C_feature_stability.pdf"),
       plot = panel_C, width = 127, height = 152, units = "mm", device = pdf)

stab_export <- tibble(gene = names(feature_freq), freq = as.numeric(feature_freq)) |>
  filter(freq > 0) |>
  arrange(desc(freq))
write_csv(stab_export, file.path(DAT_DIR, "fig6_panel_C_feature_stability.csv"))
cat("  Saved: panel_C_feature_stability.pdf + fig6_panel_C_feature_stability.csv\n")

# === 17. PANEL D — DIABLO Arrow/Sample Plot ===================================

cat("\n=== Panel D: DIABLO Arrow/Sample Plot ===\n")

if (!is.null(diablo_final)) {

  # Extract per-block scores
  block_names <- names(diablo_final$variates)
  block_names <- setdiff(block_names, "Y")  # exclude response

  if (diablo_ncomp >= 2) {
    # Build arrow data: per subject, per block, Comp1 + Comp2
    arrow_data <- lapply(block_names, function(bn) {
      scores <- diablo_final$variates[[bn]]
      tibble(
        subject = rownames(scores),
        comp1   = scores[, 1],
        comp2   = scores[, 2],
        block   = bn
      )
    }) |> bind_rows()

    # Consensus centroid per subject (mean across blocks)
    consensus <- arrow_data |>
      group_by(subject) |>
      summarise(c1_mean = mean(comp1), c2_mean = mean(comp2), .groups = "drop")

    arrow_data <- arrow_data |>
      left_join(consensus, by = "subject") |>
      mutate(age = pheno_delta$age[match(subject, pheno_delta$subject)])

    # BER annotation
    ber_label <- ""
    if (exists("diablo_perf") && !is.null(diablo_perf)) {
      ber_vals <- diablo_perf$WeightedVote.error.rate
      if ("centroids.dist" %in% names(ber_vals)) {
        ber_mat_d <- ber_vals$centroids.dist
        if ("Overall.BER" %in% rownames(ber_mat_d)) {
          ber_overall <- ber_mat_d["Overall.BER", diablo_ncomp]
        } else if ("Overall.ER" %in% rownames(ber_mat_d)) {
          ber_overall <- ber_mat_d["Overall.ER", diablo_ncomp]
        } else {
          ber_overall <- ber_mat_d[nrow(ber_mat_d), diablo_ncomp]
        }
        ber_label <- sprintf("BER = %.3f", ber_overall)
      }
    }

    # Permutation p-value annotation
    perm_label <- ""
    if (!is.na(perm_pval) && length(perm_ber_valid) > 0) {
      perm_label <- sprintf(" | Permutation p = %.3f (n = %d)", perm_pval, length(perm_ber_valid))
    }

    panel_D <- ggplot(arrow_data, aes(color = age)) +
      # Arrows from block position to consensus centroid
      geom_segment(aes(x = comp1, y = comp2, xend = c1_mean, yend = c2_mean),
                   alpha = 0.3, linewidth = 0.3,
                   arrow = arrow(length = unit(1.5, "mm"), type = "closed")) +
      # Block points — all circles (shape = 16)
      geom_point(aes(x = comp1, y = comp2), size = 2, alpha = 0.7, shape = 16) +
      # Consensus points (larger)
      geom_point(data = consensus |>
                   mutate(age = pheno_delta$age[match(subject, pheno_delta$subject)]),
                 aes(x = c1_mean, y = c2_mean), size = 3, shape = 4, stroke = 1) +
      stat_ellipse(data = consensus |>
                     mutate(age = pheno_delta$age[match(subject, pheno_delta$subject)]),
                   aes(x = c1_mean, y = c2_mean, group = age),
                   level = 0.80, linetype = "dashed", linewidth = 0.4,
                   show.legend = FALSE) +
      scale_color_manual(values = AGE_COLORS, name = "Age") +
      annotate("text", x = Inf, y = -Inf, label = ber_label,
               hjust = 1.1, vjust = -0.5, size = 2.2, fontface = "italic",
               color = "grey30") +
      labs(
        title    = "D  DIABLO Multi-Block Integration",
        subtitle = sprintf("%d blocks, %d subjects | Arrows: block -> consensus%s",
                           length(block_names), length(unique(arrow_data$subject)),
                           perm_label),
        x        = "Component 1",
        y        = "Component 2"
      ) +
      THEME_PUB +
      theme(legend.position  = "bottom",
            legend.box       = "horizontal")
  } else {
    # Fallback: grouped boxplot of Comp1 per block
    box_data <- lapply(block_names, function(bn) {
      scores <- diablo_final$variates[[bn]]
      tibble(
        subject = rownames(scores),
        comp1   = scores[, 1],
        block   = bn,
        age     = pheno_delta$age[match(rownames(scores), pheno_delta$subject)]
      )
    }) |> bind_rows()

    # Permutation p-value annotation
    perm_label <- ""
    if (!is.na(perm_pval) && length(perm_ber_valid) > 0) {
      perm_label <- sprintf(" | Permutation p = %.3f (n = %d)", perm_pval, length(perm_ber_valid))
    }

    panel_D <- ggplot(box_data, aes(x = block, y = comp1, fill = age)) +
      geom_boxplot(width = 0.6, outlier.shape = NA, alpha = 0.7,
                   position = position_dodge(0.7)) +
      geom_point(aes(color = age), position = position_jitterdodge(
        jitter.width = 0.1, dodge.width = 0.7), size = 1.5, alpha = 0.5, shape = 16) +
      scale_fill_manual(values = AGE_COLORS, name = "Age") +
      scale_color_manual(values = AGE_COLORS, guide = "none") +
      labs(
        title    = "D  DIABLO Component 1 by Block",
        subtitle = sprintf("%d blocks, %d subjects%s", length(block_names),
                           length(unique(box_data$subject)), perm_label),
        x        = "Data Block",
        y        = "Component 1 Score"
      ) +
      THEME_PUB +
      theme(legend.position = "bottom")
  }

  cat("  Panel D built.\n")
} else {
  # DIABLO failed: placeholder
  panel_D <- ggplot() +
    annotate("text", x = 0.5, y = 0.5,
             label = "DIABLO analysis did not converge", size = 3) +
    labs(title = "D  DIABLO Multi-Block Integration") +
    THEME_PUB +
    theme(axis.text = element_blank(), axis.ticks = element_blank())
  cat("  Panel D placeholder built (DIABLO failed).\n")
}

# Save Panel D
ggsave(file.path(RPT_DIR, "panel_D_diablo_arrows.pdf"),
       plot = panel_D, width = 152, height = 127, units = "mm", device = pdf)

# Save Panel D source data
if (!is.null(diablo_final)) {
  if (diablo_ncomp >= 2 && exists("arrow_data")) {
    write_csv(arrow_data, file.path(DAT_DIR, "fig6_panel_D_diablo_arrows.csv"))
  } else if (exists("box_data")) {
    write_csv(box_data, file.path(DAT_DIR, "fig6_panel_D_diablo_scores.csv"))
  }
}
cat("  Saved: panel_D_diablo_arrows.pdf\n")

# === 18. PANEL E — DIABLO Circos Plot =========================================

cat("\n=== Panel E: DIABLO Circos Plot ===\n")

if (!is.null(diablo_final)) {
  # Render circos plot to temporary PNG
  circos_tmp <- tempfile(fileext = ".png")
  png(circos_tmp, width = 380 * 4, height = 380 * 4, units = "px", res = 300,
      bg = "white")

  # Colorblind-safe block colors (Okabe-Ito derived)
  block_colors_circos <- c(baseline  = "#0072B2",  # blue
                           delta     = "#009E73",  # teal
                           phenotype = "#E69F00")  # orange
  block_colors_circos <- block_colors_circos[names(diablo_blocks)]

  # Try different cutoffs until we get a non-empty plot
  circos_cutoff <- 0.7
  circos_ok <- FALSE

  for (co in c(0.7, 0.6, 0.5, 0.4)) {
    tryCatch({
      circosPlot(diablo_final, cutoff = co, line = TRUE,
                 color.blocks = unname(block_colors_circos),
                 size.labels = 1.2, size.variables = 0.5)
      circos_cutoff <- co
      circos_ok <- TRUE
      break
    }, error = function(e) {
      cat(sprintf("    circosPlot cutoff %.1f failed: %s\n", co, conditionMessage(e)))
    })
  }

  if (!circos_ok) {
    # Draw a placeholder message
    plot.new()
    text(0.5, 0.5, "Circos plot could not be generated\n(insufficient cross-block correlations)",
         cex = 1.2)
  }
  dev.off()

  # Read back as raster and embed in ggplot
  circos_img  <- readPNG(circos_tmp)
  circos_grob <- rasterGrob(circos_img, interpolate = TRUE)

  panel_E <- ggplot() +
    annotation_custom(circos_grob, xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf) +
    labs(
      title    = "E  Cross-Block Feature Correlations",
      subtitle = sprintf("DIABLO circos | cutoff = %.1f", circos_cutoff)
    ) +
    THEME_PUB +
    theme(axis.text = element_blank(), axis.ticks = element_blank(),
          axis.title = element_blank(), panel.grid = element_blank(),
          panel.border = element_blank(), plot.background = element_blank(),
          panel.background = element_blank())

  unlink(circos_tmp)
  cat(sprintf("  Panel E built (cutoff = %.1f).\n", circos_cutoff))
} else {
  panel_E <- ggplot() +
    annotate("text", x = 0.5, y = 0.5,
             label = "DIABLO analysis did not converge", size = 3) +
    labs(title = "E  Cross-Block Feature Correlations") +
    THEME_PUB +
    theme(axis.text = element_blank(), axis.ticks = element_blank())
  cat("  Panel E placeholder built.\n")
}

# Save Panel E
ggsave(file.path(RPT_DIR, "panel_E_circos.pdf"),
       plot = panel_E, width = 152, height = 152, units = "mm", device = pdf)
cat("  Saved: panel_E_circos.pdf\n")

# === 19. PANEL F — Multi-Method Convergence (UpSet Plot) ======================

cat("\n=== Panel F: Multi-Method Convergence ===\n")

# Define 5 feature sets (all as gene symbols)

# Set 1: Limma interaction DEPs (p < 0.05)
limma_interaction <- dep_df |>
  filter(P.Value_Interaction < 0.05) |>
  pull(gene) |>
  unique()
limma_interaction <- limma_interaction[!is.na(limma_interaction) & limma_interaction != ""]
cat(sprintf("  Limma interaction DEPs: %d\n", length(limma_interaction)))

# Set 2: FCM core members (membership > 0.5)
fcm_core <- fcm_df |>
  filter(membership > 0.5) |>
  pull(gene) |>
  unique()
fcm_core <- fcm_core[!is.na(fcm_core) & fcm_core != ""]
cat(sprintf("  FCM core members: %d\n", length(fcm_core)))

# Set 3: WGCNA hub proteins
wgcna_hubs <- hub_df$gene |> unique()
wgcna_hubs <- wgcna_hubs[!is.na(wgcna_hubs) & wgcna_hubs != ""]
cat(sprintf("  WGCNA hub proteins: %d\n", length(wgcna_hubs)))

# Set 4: sPLS stable features (>70% CV selection)
spls_stable <- stable_features
cat(sprintf("  sPLS stable features: %d\n", length(spls_stable)))

# Set 5: DIABLO selected features
diablo_selected <- character(0)
if (!is.null(diablo_final)) {
  for (bn in setdiff(names(diablo_final$loadings), "Y")) {
    for (comp_i in seq_len(diablo_ncomp)) {
      sel <- selectVar(diablo_final, block = bn, comp = comp_i)
      if (!is.null(sel) && !is.null(sel[[1]]$name)) {
        diablo_selected <- c(diablo_selected, sel[[1]]$name)
      }
    }
  }
  diablo_selected <- unique(diablo_selected)
}
cat(sprintf("  DIABLO selected features: %d\n", length(diablo_selected)))

# Build UpSet input: list of sets
upset_sets <- list(
  "Limma Interaction" = limma_interaction,
  "FCM Core"          = fcm_core,
  "WGCNA Hubs"        = wgcna_hubs,
  "sPLS Stable"       = spls_stable,
  "DIABLO"            = diablo_selected
)

# Remove empty sets
upset_sets <- upset_sets[sapply(upset_sets, length) > 0]
cat(sprintf("  Non-empty sets: %d\n", length(upset_sets)))

# Build binary membership matrix for UpSetR
all_genes <- unique(unlist(upset_sets))
upset_mat <- data.frame(row.names = all_genes)
for (sn in names(upset_sets)) {
  upset_mat[[sn]] <- as.integer(all_genes %in% upset_sets[[sn]])
}

# Render UpSet to temporary PNG
upset_tmp <- tempfile(fileext = ".png")
png(upset_tmp, width = 380 * 4, height = 300 * 4, units = "px", res = 300,
    bg = "white")

tryCatch({
  # UpSetR requires at least 2 sets
  if (length(upset_sets) >= 2) {
    upset(upset_mat, sets = names(upset_sets), order.by = "freq",
          nsets = length(upset_sets), nintersects = 25,
          main.bar.color = "grey30", sets.bar.color = "grey30",
          text.scale = c(1.3, 1.2, 1.1, 1.1, 1.3, 1.1),
          point.size = 2.5, line.size = 0.8,
          mb.ratio = c(0.6, 0.4))
  } else {
    plot.new()
    text(0.5, 0.5, "Insufficient sets for UpSet plot", cex = 1.2)
  }
}, error = function(e) {
  cat(sprintf("  UpSet rendering error: %s\n", conditionMessage(e)))
  plot.new()
  text(0.5, 0.5, paste("UpSet error:", conditionMessage(e)), cex = 0.8)
})
dev.off()

# Read back as raster
upset_img  <- readPNG(upset_tmp)
upset_grob <- rasterGrob(upset_img, interpolate = TRUE)

# Count consensus genes (in >= 3 methods)
method_counts <- rowSums(upset_mat)
consensus_genes <- names(method_counts[method_counts >= 3])
consensus_label <- sprintf("%d genes in >= 3 methods", length(consensus_genes))

panel_F <- ggplot() +
  annotation_custom(upset_grob, xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf) +
  labs(
    title    = "F  Multi-Method Convergence",
    subtitle = sprintf("%d total genes | %s", length(all_genes), consensus_label)
  ) +
  THEME_PUB +
  theme(axis.text = element_blank(), axis.ticks = element_blank(),
        axis.title = element_blank(), panel.grid = element_blank(),
        panel.border = element_blank(), plot.background = element_blank(),
        panel.background = element_blank())

unlink(upset_tmp)
cat("  Panel F built.\n")

# Save Panel F
ggsave(file.path(RPT_DIR, "panel_F_upset_convergence.pdf"),
       plot = panel_F, width = 178, height = 127, units = "mm", device = pdf)

# Save Panel F source data
write_csv(upset_mat |> rownames_to_column("gene"),
          file.path(DAT_DIR, "fig6_panel_F_upset_data.csv"))

# Save consensus genes
if (length(consensus_genes) > 0) {
  write_csv(tibble(gene = consensus_genes,
                   n_methods = method_counts[consensus_genes]),
            file.path(DAT_DIR, "fig6_consensus_genes.csv"))
  cat(sprintf("  Saved: fig6_consensus_genes.csv (%d genes)\n", length(consensus_genes)))
}
cat("  Saved: panel_F_upset_convergence.pdf + fig6_panel_F_upset_data.csv\n")

# === 20. FINAL ASSEMBLY =======================================================

cat("\n=== Assembling Figure 6 ===\n")

row1 <- (panel_A | panel_B) + plot_layout(widths = c(0.45, 0.55))
row2 <- (panel_C | panel_D) + plot_layout(widths = c(0.45, 0.55))
row3 <- (panel_E | panel_F) + plot_layout(widths = c(0.45, 0.55))

fig6 <- (row1 / row2 / row3) + plot_layout(heights = c(0.33, 0.33, 0.34))

# Save PDF and PNG
ggsave(file.path(RPT_DIR, "Figure_6.pdf"), fig6,
       width = 380, height = 520, units = "mm", device = pdf, bg = "white")
cat("  Saved: Figure_6.pdf\n")

ggsave(file.path(RPT_DIR, "Figure_6.png"), fig6,
       width = 380, height = 520, units = "mm", dpi = 300, bg = "white")
cat("  Saved: Figure_6.png\n")

cat("\n=== Figure 6 complete ===\n")
