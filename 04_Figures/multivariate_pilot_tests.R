# ──────────────────────────────────────────────────────────────────────────────
# multivariate_pilot_tests.R
# Pilot alternative multivariate tests for reversal (Q1) and concordance (Q2)
# Standalone script — numeric results only, no figures
#
# Q1 (Reversal):   Does training in Old move the aging-affected proteome
#                   back toward the Young state?
# Q2 (Concordance): Do Young and Old respond to training in the same
#                   multivariate direction?
#
# Groups:
#   1. Limma rotation tests (fry, camera) — exploit protein correlation structure
#   2. Distance-based tests (PERMANOVA, paired Euclidean) — geometry in protein space
#   3. Projection-based tests — scalar summary per subject, paired t-test
#   4. fGSEA enrichment — rank-based, nonparametric
# ──────────────────────────────────────────────────────────────────────────────

setwd(rprojroot::find_rstudio_root_file())
pacman::p_load(limma, vegan, fgsea, permute, dplyr)
set.seed(42)

# ── Load data ────────────────────────────────────────────────────────────────

dal     <- readRDS("03_DEP/c_data/01_limma_DAList.rds")
dep     <- read.csv("03_DEP/c_data/03_combined_results.csv", check.names = FALSE)
imp_csv <- read.csv("02_Imputation/c_data/01_imputed.csv", check.names = FALSE)
dal_imp <- readRDS("02_Imputation/c_data/01_DAList_imputed.rds")

# ── Build imputed matrix ─────────────────────────────────────────────────────

imp_mat <- as.matrix(imp_csv[, -(1:4)])
rownames(imp_mat) <- imp_csv$uniprot_id

# ── Metadata from imputed DAList ─────────────────────────────────────────────

meta <- dal_imp$metadata
# Standardize column names
col_map <- c(Col_ID = "sample_id", Group = "age", Timepoint = "time",
             Group_Time = "group")
for (orig in names(col_map)) {
  target <- col_map[orig]
  if (orig %in% names(meta) && !target %in% names(meta))
    meta[[target]] <- meta[[orig]]
}
meta$subject_key <- sub("_(Pre|Post)$", "", meta$sample_id)

stopifnot(all(meta$sample_id %in% colnames(imp_mat)))
imp_mat <- imp_mat[, meta$sample_id]

cat("Data loaded:\n")
cat(sprintf("  Limma matrix:   %d proteins x %d samples\n",
            nrow(dal$data), ncol(dal$data)))
cat(sprintf("  Imputed matrix: %d proteins x %d samples\n",
            nrow(imp_mat), ncol(imp_mat)))

# ── Define signatures (nominal P < 0.05) ─────────────────────────────────────

aging_up   <- dep$uniprot_id[!is.na(dep$P.Value_Aging) &
                               dep$P.Value_Aging < 0.05 & dep$logFC_Aging > 0]
aging_down <- dep$uniprot_id[!is.na(dep$P.Value_Aging) &
                               dep$P.Value_Aging < 0.05 & dep$logFC_Aging < 0]
aging_all  <- union(aging_up, aging_down)

ty_up   <- dep$uniprot_id[!is.na(dep$P.Value_Training_Young) &
                             dep$P.Value_Training_Young < 0.05 &
                             dep$logFC_Training_Young > 0]
ty_down <- dep$uniprot_id[!is.na(dep$P.Value_Training_Young) &
                             dep$P.Value_Training_Young < 0.05 &
                             dep$logFC_Training_Young < 0]
ty_all  <- union(ty_up, ty_down)

cat(sprintf("\nSignatures:\n"))
cat(sprintf("  Aging:          %d up + %d down = %d\n",
            length(aging_up), length(aging_down), length(aging_all)))
cat(sprintf("  Training_Young: %d up + %d down = %d\n",
            length(ty_up), length(ty_down), length(ty_all)))

# ── Results collector ────────────────────────────────────────────────────────

results <- list()

add_result <- function(test_name, grp, question, stat_name, stat_val,
                       p_val, n_feat, interp) {
  results[[length(results) + 1]] <<- data.frame(
    test_name       = test_name,
    group           = grp,
    question        = question,
    statistic_name  = stat_name,
    statistic_value = round(stat_val, 4),
    p_value         = signif(p_val, 4),
    n_features      = n_feat,
    interpretation  = interp,
    stringsAsFactors = FALSE
  )
}

# ══════════════════════════════════════════════════════════════════════════════
# GROUP 1: Limma Rotation Tests (non-imputed matrix + limma design)
# ══════════════════════════════════════════════════════════════════════════════

cat("\n-- Group 1: Limma rotation tests --\n")

design_mat   <- dal$design$design_matrix
contrast_mat <- dal$design$contrast_matrix

# fry/camera require a complete matrix (no NAs). The limma DAList's non-imputed
# matrix has NAs by design. Use the imputed matrix aligned to the design matrix.
limma_samples <- rownames(design_mat)
stopifnot(all(limma_samples %in% colnames(imp_mat)))
expr_rot <- imp_mat[, limma_samples]
cat(sprintf("  Rotation matrix: %d proteins x %d samples (imputed, aligned to design)\n",
            nrow(expr_rot), ncol(expr_rot)))

# Within-subject correlation from the fitted model
cor_val <- dal$eBayes_fit$correlation
if (is.null(cor_val)) cor_val <- 0

# Block variable from limma metadata
meta_limma <- dal$metadata
if ("subject" %in% names(meta_limma)) {
  block_var <- meta_limma$subject
} else {
  block_var <- sub("_(Pre|Post)$", "", limma_samples)
}

cat(sprintf("  Within-subject correlation: %.4f\n", cor_val))
cat(sprintf("  Unique subjects: %d\n", length(unique(block_var))))

# Row indices in expression matrix for each gene set
all_prot       <- rownames(expr_rot)
idx_aging_up   <- which(all_prot %in% aging_up)
idx_aging_down <- which(all_prot %in% aging_down)
idx_aging_all  <- which(all_prot %in% aging_all)
idx_ty_up      <- which(all_prot %in% ty_up)
idx_ty_down    <- which(all_prot %in% ty_down)
idx_ty_all     <- which(all_prot %in% ty_all)

# Training_Old contrast vector
to_contrast <- contrast_mat[, "Training_Old"]

# Helper: one-sided p-value for a specific expected direction from fry output
fry_dir_p <- function(fry_res, set_name, target_dir) {
  dir  <- as.character(fry_res[set_name, "Direction"])
  pval <- as.numeric(fry_res[set_name, "PValue"])
  if (dir == target_dir) pval else 1 - pval
}

# ── Test 1: fry Q1 (reversal) ───────────────────────────────────────────────
# Aging-up proteins should go DOWN in Training_Old = reversal
# Aging-down proteins should go UP in Training_Old = reversal
fry_q1 <- fry(
  expr_rot,
  index    = list(aging_up = idx_aging_up, aging_down = idx_aging_down),
  design   = design_mat,
  contrast = to_contrast,
  block    = block_var,
  correlation = cor_val
)

p_rev_up   <- fry_dir_p(fry_q1, "aging_up",   "Down")
p_rev_down <- fry_dir_p(fry_q1, "aging_down", "Up")

# Fisher's combined test (df = 4 for two independent one-sided p-values)
chi2_q1  <- -2 * (log(pmax(p_rev_up, 1e-300)) + log(pmax(p_rev_down, 1e-300)))
p_fry_q1 <- pchisq(chi2_q1, df = 4, lower.tail = FALSE)

cat(sprintf("  fry Q1: aging_up %s p=%.3g, aging_down %s p=%.3g -> Fisher p=%.3g\n",
            fry_q1["aging_up", "Direction"],   fry_q1["aging_up", "PValue"],
            fry_q1["aging_down", "Direction"], fry_q1["aging_down", "PValue"],
            p_fry_q1))

add_result("fry", "rotation", "Q1_reversal", "Fisher_chi2",
           chi2_q1, p_fry_q1, length(idx_aging_all),
           sprintf("aging_up->%s p=%.2g; aging_down->%s p=%.2g",
                   fry_q1["aging_up", "Direction"],   fry_q1["aging_up", "PValue"],
                   fry_q1["aging_down", "Direction"], fry_q1["aging_down", "PValue"]))

# ── Test 2: fry Q2 (concordance) ────────────────────────────────────────────
# TY-up should go UP in Training_Old = concordance
# TY-down should go DOWN in Training_Old = concordance
fry_q2 <- fry(
  expr_rot,
  index    = list(ty_up = idx_ty_up, ty_down = idx_ty_down),
  design   = design_mat,
  contrast = to_contrast,
  block    = block_var,
  correlation = cor_val
)

p_con_up   <- fry_dir_p(fry_q2, "ty_up",   "Up")
p_con_down <- fry_dir_p(fry_q2, "ty_down", "Down")

chi2_q2  <- -2 * (log(pmax(p_con_up, 1e-300)) + log(pmax(p_con_down, 1e-300)))
p_fry_q2 <- pchisq(chi2_q2, df = 4, lower.tail = FALSE)

cat(sprintf("  fry Q2: ty_up %s p=%.3g, ty_down %s p=%.3g -> Fisher p=%.3g\n",
            fry_q2["ty_up", "Direction"],   fry_q2["ty_up", "PValue"],
            fry_q2["ty_down", "Direction"], fry_q2["ty_down", "PValue"],
            p_fry_q2))

add_result("fry", "rotation", "Q2_concordance", "Fisher_chi2",
           chi2_q2, p_fry_q2, length(idx_ty_all),
           sprintf("ty_up->%s p=%.2g; ty_down->%s p=%.2g",
                   fry_q2["ty_up", "Direction"],   fry_q2["ty_up", "PValue"],
                   fry_q2["ty_down", "Direction"], fry_q2["ty_down", "PValue"]))

# ── Test 3: camera Q1 (reversal — competitive) ──────────────────────────────
# Is the aging signature MORE differentially expressed in Training_Old than
# random proteins? Two-sided test; direction indicates net shift.
# camera.default does not accept block/correlation — it handles inter-gene
# correlation via its own inter.gene.cor parameter (default 0.01)
cam_q1 <- camera(
  expr_rot,
  index    = list(aging_sig = idx_aging_all),
  design   = design_mat,
  contrast = to_contrast
)

cam_q1_dir <- as.character(cam_q1["aging_sig", "Direction"])
cam_q1_p2  <- as.numeric(cam_q1["aging_sig", "PValue"])  # two-sided

cat(sprintf("  camera Q1: %s p=%.3g (two-sided)\n", cam_q1_dir, cam_q1_p2))

add_result("camera", "rotation", "Q1_reversal", "Direction",
           ifelse(cam_q1_dir == "Down", -1, 1), cam_q1_p2,
           length(idx_aging_all),
           sprintf("aging sig net %s in Training_Old (2-sided)", cam_q1_dir))

# ── Test 4: camera Q2 (concordance — competitive) ───────────────────────────
cam_q2 <- camera(
  expr_rot,
  index    = list(ty_sig = idx_ty_all),
  design   = design_mat,
  contrast = to_contrast
)

cam_q2_dir <- as.character(cam_q2["ty_sig", "Direction"])
cam_q2_p2  <- as.numeric(cam_q2["ty_sig", "PValue"])

cat(sprintf("  camera Q2: %s p=%.3g (two-sided)\n", cam_q2_dir, cam_q2_p2))

add_result("camera", "rotation", "Q2_concordance", "Direction",
           ifelse(cam_q2_dir == "Up", 1, -1), cam_q2_p2,
           length(idx_ty_all),
           sprintf("TY sig net %s in Training_Old (2-sided)", cam_q2_dir))

# ══════════════════════════════════════════════════════════════════════════════
# GROUP 2: Distance-Based Tests (imputed matrix)
# ══════════════════════════════════════════════════════════════════════════════

cat("\n-- Group 2: Distance-based tests --\n")

# Subset imputed matrix to aging sig proteins present in both datasets
aging_imp <- intersect(aging_all, rownames(imp_mat))
imp_aging <- imp_mat[aging_imp, ]

# Sample group indices
yp_samples <- meta$sample_id[meta$age == "Young" & meta$time == "Pre"]
op_samples <- meta$sample_id[meta$age == "Old"   & meta$time == "Pre"]
oo_samples <- meta$sample_id[meta$age == "Old"   & meta$time == "Post"]

# Old subject pairing
old_meta    <- meta[meta$age == "Old", ]
old_subjects <- unique(old_meta$subject_key)
n_old       <- length(old_subjects)

# Young Pre centroid on aging signature
yp_centroid <- rowMeans(imp_aging[, yp_samples])

cat(sprintf("  Aging sig in imputed: %d proteins\n", length(aging_imp)))
cat(sprintf("  Old subjects: %d\n", n_old))

# ── Test 5: PERMANOVA (time effect in Old, within-subject permutations) ──────
old_all_samples <- meta$sample_id[meta$age == "Old"]
imp_old_t       <- t(imp_aging[, old_all_samples])  # samples x proteins
old_meta_perm   <- old_meta[match(rownames(imp_old_t), old_meta$sample_id), ]
old_meta_perm$time <- factor(old_meta_perm$time, levels = c("Pre", "Post"))

dist_old <- dist(imp_old_t, method = "euclidean")

# Generate explicit within-subject permutation matrix via shuffleSet
perm_ctrl <- how(blocks = factor(old_meta_perm$subject_key))
perm_mat  <- shuffleSet(n = nrow(old_meta_perm), nset = 9999, control = perm_ctrl)

perm_res <- adonis2(dist_old ~ time, data = old_meta_perm, permutations = perm_mat)

perm_f  <- perm_res["time", "F"]
perm_p  <- perm_res["time", "Pr(>F)"]
perm_r2 <- perm_res["time", "R2"]

cat(sprintf("  PERMANOVA: F=%.3f, R2=%.4f, p=%.4g\n", perm_f, perm_r2, perm_p))

add_result("PERMANOVA", "distance", "Q1_reversal", "F",
           perm_f, perm_p, length(aging_imp),
           sprintf("Time effect on aging sig in Old (R2=%.4f)", perm_r2))

# ── Test 6: Paired Euclidean t-test ──────────────────────────────────────────
# Per-subject distance to Young_Pre centroid: is d_Post < d_Pre?
d_pre  <- numeric(n_old)
d_post <- numeric(n_old)

for (i in seq_along(old_subjects)) {
  s <- old_subjects[i]
  pre_col  <- old_meta$sample_id[old_meta$subject_key == s & old_meta$time == "Pre"]
  post_col <- old_meta$sample_id[old_meta$subject_key == s & old_meta$time == "Post"]
  d_pre[i]  <- sqrt(sum((imp_aging[, pre_col]  - yp_centroid)^2))
  d_post[i] <- sqrt(sum((imp_aging[, post_col] - yp_centroid)^2))
}

tt_euc      <- t.test(d_post, d_pre, paired = TRUE, alternative = "less")
mean_delta  <- mean(d_post - d_pre)
pct_change  <- mean_delta / mean(d_pre) * 100

cat(sprintf("  Paired Euclidean: mean(d_post-d_pre)=%.2f (%.1f%%), t=%.3f, p=%.4g\n",
            mean_delta, pct_change, tt_euc$statistic, tt_euc$p.value))

add_result("paired_euclidean", "distance", "Q1_reversal", "t",
           as.numeric(tt_euc$statistic), tt_euc$p.value, length(aging_imp),
           sprintf("mean delta_d=%.1f (%.1f%% of d_pre), df=%d",
                   mean_delta, pct_change, n_old - 1))

# ══════════════════════════════════════════════════════════════════════════════
# GROUP 3: Projection-Based Tests (imputed matrix)
# ══════════════════════════════════════════════════════════════════════════════

cat("\n-- Group 3: Projection-based tests --\n")

# ── Test 7: Signed projection Q1 (reversal) ─────────────────────────────────
# Aging direction = Old_Pre centroid - Young_Pre centroid (young -> old)
# Project each Old subject's training vector onto this direction
# Negative projection = reversal (training moves AWAY from old, toward young)
aging_dir      <- rowMeans(imp_aging[, op_samples]) - yp_centroid
aging_dir_norm <- aging_dir / sqrt(sum(aging_dir^2))

proj_q1 <- numeric(n_old)
for (i in seq_along(old_subjects)) {
  s <- old_subjects[i]
  pre_col  <- old_meta$sample_id[old_meta$subject_key == s & old_meta$time == "Pre"]
  post_col <- old_meta$sample_id[old_meta$subject_key == s & old_meta$time == "Post"]
  train_vec   <- imp_aging[, post_col] - imp_aging[, pre_col]
  proj_q1[i]  <- sum(train_vec * aging_dir_norm)
}

tt_proj_q1 <- t.test(proj_q1, alternative = "less")

cat(sprintf("  Projection Q1: mean=%.3f, t=%.3f, p=%.4g\n",
            mean(proj_q1), tt_proj_q1$statistic, tt_proj_q1$p.value))

add_result("projection", "projection", "Q1_reversal", "t",
           as.numeric(tt_proj_q1$statistic), tt_proj_q1$p.value, length(aging_imp),
           sprintf("mean proj=%.3f (negative=reversal), df=%d",
                   mean(proj_q1), n_old - 1))

# ── Test 8: Signed projection Q2 (concordance) ──────────────────────────────
# Young training direction = Young_Post centroid - Young_Pre centroid
# Project each Old subject's training vector onto this
# Positive projection = concordance (Old training in same direction as Young)
ty_imp  <- intersect(ty_all, rownames(imp_mat))
imp_ty  <- imp_mat[ty_imp, ]

ypo_samples     <- meta$sample_id[meta$age == "Young" & meta$time == "Post"]
young_train_dir <- rowMeans(imp_ty[, ypo_samples]) - rowMeans(imp_ty[, yp_samples])
yt_dir_norm     <- young_train_dir / sqrt(sum(young_train_dir^2))

proj_q2 <- numeric(n_old)
for (i in seq_along(old_subjects)) {
  s <- old_subjects[i]
  pre_col  <- old_meta$sample_id[old_meta$subject_key == s & old_meta$time == "Pre"]
  post_col <- old_meta$sample_id[old_meta$subject_key == s & old_meta$time == "Post"]
  train_vec   <- imp_ty[, post_col] - imp_ty[, pre_col]
  proj_q2[i]  <- sum(train_vec * yt_dir_norm)
}

tt_proj_q2 <- t.test(proj_q2, alternative = "greater")

cat(sprintf("  Projection Q2: mean=%.3f, t=%.3f, p=%.4g\n",
            mean(proj_q2), tt_proj_q2$statistic, tt_proj_q2$p.value))

add_result("projection", "projection", "Q2_concordance", "t",
           as.numeric(tt_proj_q2$statistic), tt_proj_q2$p.value, length(ty_imp),
           sprintf("mean proj=%.3f (positive=concordance), df=%d",
                   mean(proj_q2), n_old - 1))

# ══════════════════════════════════════════════════════════════════════════════
# GROUP 4: fGSEA (rank-based enrichment)
# ══════════════════════════════════════════════════════════════════════════════

cat("\n-- Group 4: fGSEA --\n")

# Rank all proteins by t_Training_Old
ranks <- setNames(dep$t_Training_Old, dep$uniprot_id)
ranks <- ranks[!is.na(ranks)]

# Gene sets (intersect with ranked proteins)
gs <- list(
  aging_up   = intersect(aging_up,   names(ranks)),
  aging_down = intersect(aging_down, names(ranks)),
  ty_up      = intersect(ty_up,      names(ranks)),
  ty_down    = intersect(ty_down,    names(ranks))
)

fgsea_res <- fgseaMultilevel(gs, ranks, minSize = 5, maxSize = Inf,
                              nPermSimple = 10000)

# ── Test 9: fGSEA Q1 (reversal) ─────────────────────────────────────────────
# Reversal: aging_up should have negative NES, aging_down should have positive
nes_au <- fgsea_res$NES[fgsea_res$pathway == "aging_up"]
p_au   <- fgsea_res$pval[fgsea_res$pathway == "aging_up"]
nes_ad <- fgsea_res$NES[fgsea_res$pathway == "aging_down"]
p_ad   <- fgsea_res$pval[fgsea_res$pathway == "aging_down"]

# One-sided p-values for reversal direction
p_gsea_rev_up <- if (nes_au < 0) p_au / 2 else 1 - p_au / 2
p_gsea_rev_dn <- if (nes_ad > 0) p_ad / 2 else 1 - p_ad / 2

chi2_gsea_q1 <- -2 * (log(pmax(p_gsea_rev_up, 1e-300)) +
                       log(pmax(p_gsea_rev_dn, 1e-300)))
p_gsea_q1    <- pchisq(chi2_gsea_q1, df = 4, lower.tail = FALSE)

cat(sprintf("  fGSEA Q1: aging_up NES=%.3f p=%.3g, aging_down NES=%.3f p=%.3g -> Fisher p=%.3g\n",
            nes_au, p_au, nes_ad, p_ad, p_gsea_q1))

add_result("fGSEA", "enrichment", "Q1_reversal", "Fisher_chi2",
           chi2_gsea_q1, p_gsea_q1, length(aging_all),
           sprintf("aging_up NES=%.2f p=%.2g; aging_down NES=%.2f p=%.2g",
                   nes_au, p_au, nes_ad, p_ad))

# ── Test 10: fGSEA Q2 (concordance) ─────────────────────────────────────────
# Concordance: ty_up should have positive NES, ty_down should have negative
nes_tu <- fgsea_res$NES[fgsea_res$pathway == "ty_up"]
p_tu   <- fgsea_res$pval[fgsea_res$pathway == "ty_up"]
nes_td <- fgsea_res$NES[fgsea_res$pathway == "ty_down"]
p_td   <- fgsea_res$pval[fgsea_res$pathway == "ty_down"]

p_gsea_con_up <- if (nes_tu > 0) p_tu / 2 else 1 - p_tu / 2
p_gsea_con_dn <- if (nes_td < 0) p_td / 2 else 1 - p_td / 2

chi2_gsea_q2 <- -2 * (log(pmax(p_gsea_con_up, 1e-300)) +
                       log(pmax(p_gsea_con_dn, 1e-300)))
p_gsea_q2    <- pchisq(chi2_gsea_q2, df = 4, lower.tail = FALSE)

cat(sprintf("  fGSEA Q2: ty_up NES=%.3f p=%.3g, ty_down NES=%.3f p=%.3g -> Fisher p=%.3g\n",
            nes_tu, p_tu, nes_td, p_td, p_gsea_q2))

add_result("fGSEA", "enrichment", "Q2_concordance", "Fisher_chi2",
           chi2_gsea_q2, p_gsea_q2, length(ty_all),
           sprintf("ty_up NES=%.2f p=%.2g; ty_down NES=%.2f p=%.2g",
                   nes_tu, p_tu, nes_td, p_td))

# ══════════════════════════════════════════════════════════════════════════════
# Compile and export
# ══════════════════════════════════════════════════════════════════════════════

summary_df <- bind_rows(results)

cat("\n")
cat("================================================================\n")
cat("PILOT TEST SUMMARY\n")
cat("================================================================\n\n")

# Formatted console output
for (i in seq_len(nrow(summary_df))) {
  row <- summary_df[i, ]
  pv <- row$p_value
  sig <- if (is.na(pv)) "NA" else if (pv < 0.01) "***" else
    if (pv < 0.05) "**" else if (pv < 0.10) "*" else ""
  cat(sprintf("  [%2d] %-18s %-16s p = %-10s %s\n",
              i, row$test_name, row$question, format(row$p_value, digits = 3), sig))
  cat(sprintf("       %s\n", row$interpretation))
}

cat("\n  Significance: *** p<0.01, ** p<0.05, * p<0.10\n")

# Reference p-values from existing tests
cat("\n  Reference (existing tests):\n")
cat("    Melov permutation (F05 reversal):      p = 0.15\n")
cat("    Cosine similarity (F04 concordance):    p = 0.066\n")

out_path <- "04_Figures/multivariate_pilot_summary.csv"
write.csv(summary_df, out_path, row.names = FALSE)
cat(sprintf("\nSaved: %s (%d rows)\n", out_path, nrow(summary_df)))
