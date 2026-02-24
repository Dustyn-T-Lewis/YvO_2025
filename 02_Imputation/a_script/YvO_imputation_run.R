#!/usr/bin/env Rscript
###############################################################################
#   YvO Imputation — Missingness, MAR/MNAR classification, benchmarking
#   Input:  cycloess-normalized log2 protein matrix (01_normalized.csv)
#
#   Reports:
#     01_missingness_report.pdf  — missingness + MAR/MNAR classification
#     02_imputation_report.pdf   — benchmark, post-imputation QC, MNAR audit
#
#   Methodological references:
#     - Lazar et al. 2016, J Proteome Res 15(4):1116-1125
#       Hybrid MAR/MNAR imputation outperforms single-method approaches
#     - Webb-Robertson et al. 2015, J Proteome Res 14(3):920-930
#       Comprehensive review of imputation strategies for MS proteomics
#     - Hediyeh-zadeh et al. 2023, Mol Cell Proteomics 22(2):100477
#       msImpute EBM for entropy-based MAR/MNAR classification
#     - Borcherding et al. 2021, Brief Bioinform 22(3):bbaa112
#       SFI-hybrid (kNN for MAR + QRILC for MNAR) outperforms single methods
#     - Jin et al. 2021, Sci Rep 11:1556
#       NRMSE benchmarking with masking for imputation evaluation
#
#   Exercise proteomics context:
#     - Hostrup et al. 2022, eLife 11:e69802
#       HIT remodels skeletal muscle proteome; used imputation + DE pipeline
#     - Deshmukh et al. 2021, Nat Commun 12:304
#       Deep muscle proteomics with fiber-type resolution; complete-case analysis
#     - Ubaida-Mohien et al. 2019, eLife 8:e49874
#       Discovery proteomics in aging muscle; >4000 proteins quantified
###############################################################################

if (!requireNamespace("pacman", quietly = TRUE)) install.packages("pacman")
pacman::p_load(
  MsCoreUtils, msImpute, pcaMethods, imputeLCMD, missForest, missMDA,
  ggplot2, patchwork, dplyr, tidyr, tibble, readr, stringr, scales
)

# ─── Paths ───────────────────────────────────────────────────────────────────
setwd("/Users/dtl0018/Desktop/A_Proteomics_Analysis/A_YvO_2025")
INPUT_CSV  <- "01_normalization/c_data/01_normalized.csv"
REPORT_DIR <- "02_Imputation/b_reports"
DATA_DIR   <- "02_Imputation/c_data"
dir.create(REPORT_DIR, showWarnings = FALSE, recursive = TRUE)
dir.create(DATA_DIR,   showWarnings = FALSE, recursive = TRUE)

# ─── Benchmark methods ───────────────────────────────────────────────────────
# Edit this list to add/remove methods. Types: MNAR, MAR, Hybrid, Custom
METHODS <- list(
  # MNAR (left-censored / Gaussian sampling from low tail)
  MinProb         = list(method="MinProb"),
  MinDet          = list(method="MinDet"),
  QRILC           = list(method="QRILC"),
  zero            = list(method="zero"),
  # MAR (global structure-based)
  knn             = list(method="knn"),
  bpca            = list(method="bpca"),
  RF              = list(method="RF"),
  SVD             = list(method="SVD"),
  imputePCA       = list(method="imputePCA"),
  GaussianSample  = list(method="GaussianSample"),
  # Hybrid (MAR method for MAR proteins, MNAR method for MNAR proteins)
  mix_bpca_MinProb = list(method="mixed", mar="bpca",  mnar="MinProb"),
  mix_bpca_QRILC   = list(method="mixed", mar="bpca",  mnar="QRILC"),
  mix_bpca_MinDet  = list(method="mixed", mar="bpca",  mnar="MinDet"),
  mix_RF_MinProb   = list(method="mixed", mar="RF",    mnar="MinProb"),
  mix_RF_QRILC     = list(method="mixed", mar="RF",    mnar="QRILC")
)

N_ITER <- 20

# Method type lookup (auto-classified from METHODS list)
mtype <- tibble(method = names(METHODS)) %>%
  mutate(type = case_when(
    method %in% c("MinProb","MinDet","QRILC","zero") ~ "MNAR",
    str_starts(method, "mix_") ~ "Hybrid",
    TRUE ~ "MAR"))

# ─── Shared aesthetics ──────────────────────────────────────────────────────
source("04_Figures/shared/palettes.R")
pal_gt   <- GROUP_FILL
pal_mar   <- c(MAR="#4393C3", MNAR="#D6604D")
pal_class <- c(Complete="#4DAF4A", MAR="#4393C3", MNAR="#D6604D")
pal_mtyp  <- c(MNAR="#D6604D", MAR="#4393C3", Hybrid="#5AAE61")
thm      <- theme_minimal(base_size = 11)
thm_sm   <- theme_minimal(base_size = 8)

# ─── Helpers ─────────────────────────────────────────────────────────────────

# Dispatch imputation by method specification
run_impute <- function(m, mat, randna) {
  if (m$method == "mixed") {
    MsCoreUtils::impute_matrix(mat, method="mixed",
                               randna=randna, mar=m$mar, mnar=m$mnar)
  } else if (m$method == "imputePCA") {
    missMDA::imputePCA(mat, ncp=2, method="Regularized")$completeObs
  } else if (m$method == "SVD") {
    pcaMethods::pca(mat, method="svdImpute", nPcs=2, verbose=FALSE)@completeObs
  } else if (m$method == "GaussianSample") {
    impute_gaussian(mat)
  } else {
    MsCoreUtils::impute_matrix(mat, method=m$method)
  }
}

# Per-protein Gaussian sampling: draw from N(mu_obs, sigma_obs) for each row
impute_gaussian <- function(mat) {
  for (i in seq_len(nrow(mat))) {
    na_idx <- is.na(mat[i, ])
    if (any(na_idx) && sum(!na_idx) >= 2) {
      mat[i, na_idx] <- rnorm(sum(na_idx),
                               mean = mean(mat[i, ], na.rm = TRUE),
                               sd   = sd(mat[i, ], na.rm = TRUE))
    }
  }
  mat
}

###############################################################################
# 1: LOAD DATA
###############################################################################
df  <- read_csv(INPUT_CSV, show_col_types = FALSE)
ann <- df[, 1:4]
mat <- as.matrix(df[, -(1:4)])
rownames(mat) <- df$gene
cat(sprintf("%d proteins x %d samples\n", nrow(mat), ncol(mat)))

meta <- tibble(Col_ID = colnames(mat)) %>%
  mutate(Group     = if_else(str_detect(Col_ID, "^(Y_|YP_)"), "Young", "Old"),
         Timepoint = if_else(str_detect(Col_ID, "_Pre$"), "Pre", "Post"),
         Group_Time = paste(Group, Timepoint, sep = "_"))
print(count(meta, Group, Timepoint))

###############################################################################
# 2: MISSINGNESS & MAR/MNAR CLASSIFICATION
###############################################################################
prot_miss <- rowSums(is.na(mat))
prot_pct  <- prot_miss / ncol(mat) * 100
obs_means <- rowMeans(mat, na.rm = TRUE)
pct_miss  <- round(sum(is.na(mat)) / length(mat) * 100, 2)

cat(sprintf("Missing: %d / %d (%.2f%%) | Complete: %d\n",
            sum(is.na(mat)), length(mat), pct_miss, sum(prot_miss == 0)))

miss_by_group <- sapply(unique(meta$Group_Time), function(g) {
  cols <- meta$Col_ID[meta$Group_Time == g]
  rowSums(is.na(mat[, cols, drop=FALSE])) / length(cols) * 100
})

has_na <- which(prot_miss > 0 & prot_miss < ncol(mat))
miss_class <- tibble(gene=rownames(mat), n_miss=prot_miss,
                     pct_miss=prot_pct, mean_intensity=obs_means)

# MAR/MNAR classification: msImpute EBM (primary) → k-means (fallback)
mar_result <- tryCatch({
  feat <- msImpute::selectFeatures(mat[has_na,], method="ebm")
  list(mar_genes = feat, method = "msImpute_ebm")
}, error = function(e) {
  cat(sprintf("msImpute EBM failed: %s\n", conditionMessage(e)))
  NULL
})

if (is.null(mar_result)) {
  # K-means on standardized (intensity, missingness) — data-driven boundary
  mc_sub <- miss_class %>% filter(n_miss > 0, n_miss < ncol(mat))
  km <- kmeans(scale(cbind(mc_sub$mean_intensity, mc_sub$pct_miss)),
               centers = 2, nstart = 25)
  cl_means <- tapply(mc_sub$mean_intensity, km$cluster, mean)
  mnar_cl <- which.min(cl_means)
  mar_result <- list(
    mar_genes = mc_sub$gene[km$cluster != mnar_cl],
    method = sprintf("k-means (cluster means: %.1f vs %.1f)",
                     cl_means[mnar_cl], cl_means[-mnar_cl]))
}

cat(sprintf("Classification: %s\n", mar_result$method))
miss_class <- miss_class %>%
  mutate(classification = case_when(
    n_miss == 0 ~ "Complete",
    gene %in% mar_result$mar_genes ~ "MAR",
    TRUE ~ "MNAR"))
print(count(miss_class, classification))

###############################################################################
# 3: REPORT 1 -- MISSINGNESS
###############################################################################

# A: Per-protein missingness histogram
p_miss_hist <- ggplot(tibble(x = prot_pct), aes(x)) +
  geom_histogram(binwidth=2, fill="#4393C3", color="white", alpha=0.8) +
  geom_vline(xintercept=c(20,50), linetype="dashed", color="red", alpha=0.6) +
  annotate("text", x=22, y=Inf, vjust=2, hjust=0, size=3,
           label=sprintf("Complete: %d\n1-20%%: %d\n20-50%%: %d\n>50%%: %d",
                         sum(prot_pct==0), sum(prot_pct>0 & prot_pct<=20),
                         sum(prot_pct>20 & prot_pct<=50), sum(prot_pct>50))) +
  labs(x="% missing per protein", y="Count", title="A: Per-Protein Missingness") + thm

# B: Classification donut (Complete / MAR / MNAR)
class_counts <- miss_class %>% count(classification)
n_total <- sum(class_counts$n)
donut_df <- class_counts %>%
  mutate(frac = n / n_total,
         label = sprintf("%s\n%d (%.1f%%)", classification, n, frac * 100),
         ymax = cumsum(frac), ymin = lag(ymax, default = 0),
         ymid = (ymin + ymax) / 2)

p_class_donut <- ggplot(donut_df, aes(ymax=ymax, ymin=ymin, xmax=4, xmin=2.5,
                                       fill=classification)) +
  geom_rect(color="white", linewidth=0.6) +
  geom_text(aes(x=3.25, y=ymid, label=label), size=3, lineheight=0.9) +
  annotate("text", x=0, y=0,
           label=sprintf("%d\nproteins\n%.1f%% missing", n_total, pct_miss),
           size=4, fontface="bold", lineheight=1.1) +
  coord_polar(theta="y") + xlim(c(0, 4.5)) +
  scale_fill_manual(values=pal_class, guide="none") +
  theme_void(base_size=11) +
  labs(title="B: Missingness Classification",
       subtitle="MAR: intensity-independent | MNAR: low-abundance driven")

# C: MAR vs MNAR scatter (intensity vs % missing)
mc <- miss_class %>% filter(classification != "Complete")
p_mar_scatter <- ggplot(mc, aes(mean_intensity, pct_miss, color=classification)) +
  geom_point(alpha=0.5, size=1.2) + scale_color_manual(values=pal_mar) +
  labs(x="Mean log2 intensity", y="% missing",
       title="C: Missingness vs Abundance by Type") +
  thm + theme(legend.position="bottom")

# D: Per-sample missingness
p_samp_miss <- tibble(Col_ID=colnames(mat), pct=colSums(is.na(mat))/nrow(mat)*100) %>%
  left_join(meta, by="Col_ID") %>%
  ggplot(aes(reorder(Col_ID, pct), pct, fill=Group_Time)) +
  geom_col(alpha=0.85) + scale_fill_manual(values=pal_gt) + coord_flip() +
  labs(x=NULL, y="% missing proteins", title="D: Per-Sample Missingness") +
  thm_sm + theme(axis.text.y=element_text(size=4))

# E: Intensity histogram by type (counts) with Complete as reference
p_int_hist <- ggplot(miss_class, aes(mean_intensity, fill=classification)) +
  geom_histogram(binwidth=0.3, alpha=0.6, position="identity", color="white",
                 linewidth=0.1) +
  scale_fill_manual(values=pal_class, name="Type") +
  labs(x="Mean log2 intensity", y="Count",
       title="E: Intensity Distribution by Classification",
       subtitle="Complete proteins shown as reference") +
  thm + theme(legend.position="bottom")

# F: Per-group missingness heatmap (top 50)
top_idx <- order(prot_pct, decreasing=TRUE)[1:min(50, sum(prot_pct>0))]
p_group_heat <- as_tibble(miss_by_group[top_idx,], rownames="gene") %>%
  pivot_longer(-gene, names_to="Group_Time", values_to="pct") %>%
  mutate(gene = factor(gene, levels=rev(rownames(mat)[top_idx]))) %>%
  ggplot(aes(Group_Time, gene, fill=pct)) +
  geom_tile(color="white", linewidth=0.3) +
  scale_fill_gradient2(low="white", mid="#FDDBC7", high="#B2182B", midpoint=50, name="% Miss") +
  labs(x=NULL, y=NULL, title="F: Per-Group Missingness (top 50)") +
  theme_minimal(base_size=9) +
  theme(axis.text.y=element_text(size=5), axis.text.x=element_text(angle=45, hjust=1))

pdf(file.path(REPORT_DIR, "01_missingness_report.pdf"), width=14, height=12)
print((p_miss_hist | p_class_donut) / (p_mar_scatter | p_samp_miss) +
        plot_annotation(title="Missingness Diagnostics & MAR/MNAR Classification",
                        subtitle=sprintf("%d proteins x %d samples | %.2f%% missing overall",
                                         nrow(mat), ncol(mat), pct_miss)))
print((p_int_hist | p_group_heat) +
        plot_annotation(title="Missingness Detail"))
dev.off()
###############################################################################
# 4: BENCHMARKING
###############################################################################
set.seed(42)
# MsCoreUtils::impute_matrix(method="mixed") defines randna=TRUE as MAR
# (i.e., "random NA"), so TRUE = apply MAR method, FALSE = apply MNAR method.
# See: rformassspectrometry.github.io/MsCoreUtils/reference/imputation.html
randna <- setNames(miss_class$classification != "MNAR", miss_class$gene)
nrmse  <- function(t, i) sqrt(mean((t - i)^2)) / sd(t)

cat(sprintf("Benchmarking: %d methods x %d iterations\n", length(METHODS), N_ITER))
res <- vector("list", length(METHODS) * N_ITER); k <- 0L
fail_res <- list(); k_fail <- 0L

for (iter in seq_len(N_ITER)) {
  if (iter %% 5 == 0) cat(sprintf("  Iter %d/%d\n", iter, N_ITER))
  # Mask observed values from MAR proteins only (randna=TRUE means MAR)
  mar_obs_idx <- which(!is.na(mat) & randna[row(mat)])
  mask_idx    <- sample(mar_obs_idx, round(length(mar_obs_idx) * 0.10))
  true_v   <- mat[mask_idx]
  mm <- mat; mm[mask_idx] <- NA

  for (nm in names(METHODS)) {
    result <- tryCatch(
      list(val = run_impute(METHODS[[nm]], mm, randna)),
      error = function(e) list(val = NULL, err = conditionMessage(e)))
    if (is.null(result$val)) {
      cat(sprintf("  %s failed: %s\n", nm, result$err))
      k_fail <- k_fail + 1L
      fail_res[[k_fail]] <- tibble(method=nm, iter=iter, error=result$err)
      next
    }
    k <- k + 1L
    res[[k]] <- tibble(method=nm, iter=iter, nrmse=nrmse(true_v, result$val[mask_idx]))
  }
}

# Log benchmark failures
if (k_fail > 0) {
  write_csv(bind_rows(fail_res), file.path(DATA_DIR, "benchmark_failures.csv"))
}

bench_df <- bind_rows(res)
bench_sum <- bench_df %>%
  group_by(method) %>%
  summarise(mean_nrmse=mean(nrmse), sd_nrmse=sd(nrmse),
            median_nrmse=median(nrmse), .groups="drop") %>%
  arrange(mean_nrmse)
print(bench_sum)

best <- bench_sum$method[1]
cat(sprintf("Best: %s (mean NRMSE = %.4f)\n", best, bench_sum$mean_nrmse[1]))
write_csv(bench_sum, file.path(DATA_DIR, "benchmark_summary.csv"))

###############################################################################
# 5: APPLY BEST METHOD
###############################################################################
set.seed(42)
mat_imp <- run_impute(METHODS[[best]], mat, randna)
stopifnot(sum(is.na(mat_imp)) == 0)

###############################################################################
# 6: REPORT 2 -- IMPUTATION
###############################################################################

was_na <- is.na(mat)
top3 <- bench_sum$method[1:3]

# Page 1: Benchmark — points for all methods, boxplot inset for top 3
visible <- bench_sum %>% filter(median_nrmse < 5)
excluded <- setdiff(bench_sum$method, visible$method)
excl_note <- if (length(excluded) > 0) {
  paste0(" | Off scale: ", paste(sprintf("%s (%.1f)", excluded,
    bench_sum$mean_nrmse[match(excluded, bench_sum$method)]), collapse=", "))
} else ""

p_bench_main <- visible %>%
  left_join(mtype, by="method") %>%
  ggplot(aes(reorder(method, median_nrmse), median_nrmse, color=type)) +
  geom_point(size=3.5) +
  geom_errorbar(aes(ymin=median_nrmse - sd_nrmse, ymax=median_nrmse + sd_nrmse),
                width=0.3, linewidth=0.5) +
  geom_hline(yintercept=visible$median_nrmse[1], linetype="dashed",
             color="#B2182B", alpha=0.5) +
  scale_color_manual(values=pal_mtyp, name="Type") + coord_flip() +
  labs(x=NULL, y="NRMSE (lower = better)", title="Imputation Benchmark",
       subtitle=sprintf("%d iter x 10%% masked | Best: %s (%.4f)%s",
                        N_ITER, best, bench_sum$mean_nrmse[1], excl_note)) + thm

p_bench_inset <- bench_df %>%
  filter(method %in% top3) %>% left_join(mtype, by="method") %>%
  ggplot(aes(reorder(method, nrmse, FUN=median), nrmse, fill=type)) +
  geom_boxplot(alpha=0.7, outlier.size=0.8) +
  scale_fill_manual(values=pal_mtyp, guide="none") + coord_flip() +
  labs(x=NULL, y=NULL, title="Top 3 (zoomed)") +
  theme_minimal(base_size=7) +
  theme(plot.background=element_rect(fill="white", color="grey40", linewidth=0.4),
        plot.title=element_text(size=7, face="bold"),
        axis.text.x=element_text(size=6), axis.text.y=element_text(size=6),
        plot.margin=margin(4,4,4,4))

# Page 2: Post-imputation quality
hist_df <- bind_rows(
  tibble(value=as.vector(mat[!was_na]), source="Observed"),
  tibble(value=mat_imp[was_na],         source="Imputed"))

p_hist <- ggplot(hist_df, aes(value, fill=source)) +
  geom_histogram(binwidth=0.2, alpha=0.6, position="identity", color="white",
                 linewidth=0.1) +
  scale_fill_manual(values=c(Observed="#2166AC", Imputed="#B2182B")) +
  labs(x="log2 intensity", y="Count",
       title="Observed vs Imputed Value Distributions",
       subtitle=sprintf("%s observed | %s imputed",
                        scales::comma(sum(!was_na)), scales::comma(sum(was_na)))) +
  thm + theme(legend.position="bottom", legend.title=element_blank())

grp_df <- do.call(rbind, lapply(unique(meta$Group_Time), function(g) {
  cols <- meta$Col_ID[meta$Group_Time==g]
  sm <- mat[, cols, drop=FALSE]; si <- mat_imp[, cols, drop=FALSE]; sn <- was_na[, cols, drop=FALSE]
  bind_rows(
    tibble(value=as.vector(sm[!is.na(sm)]), stage="Observed", Group_Time=g),
    tibble(value=as.vector(si[sn]),         stage="Imputed",  Group_Time=g))
}))

p_grp_dens <- ggplot(grp_df, aes(value, color=stage, linetype=stage)) +
  geom_density(linewidth=0.7) + facet_wrap(~Group_Time, scales="free_y") +
  scale_color_manual(values=c(Observed="#2166AC", Imputed="#B2182B")) +
  scale_linetype_manual(values=c(Observed="solid", Imputed="dashed")) +
  labs(x="log2 intensity", title="Observed vs Imputed by Group") +
  theme_minimal(base_size=10) +
  theme(legend.position="bottom", legend.title=element_blank())

# Page 3: MNAR audit
mnar_genes <- miss_class$gene[miss_class$classification == "MNAR"]
mnar_audit <- tibble(
  gene = mnar_genes,
  pre_mean  = rowMeans(mat[mnar_genes, ], na.rm = TRUE),
  post_mean = rowMeans(mat_imp[mnar_genes, ]),
  pct_miss  = prot_pct[mnar_genes],
  shift     = rowMeans(mat_imp[mnar_genes, ]) - rowMeans(mat[mnar_genes, ], na.rm = TRUE),
  imputation_reliable = prot_pct[mnar_genes] < 50
)

p_mnar_a <- ggplot(mnar_audit, aes(pre_mean, post_mean, color = pct_miss)) +
  geom_point(alpha = 0.7, size = 1.5) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey40") +
  scale_color_gradient(low = "#FDDBC7", high = "#B2182B", name = "% missing") +
  labs(x = "Pre-imputation mean log2", y = "Post-imputation mean log2",
       title = "MNAR Mean Shift After Imputation") + thm

p_mnar_b <- ggplot(mnar_audit, aes(pct_miss, shift)) +
  geom_point(alpha = 0.6, size = 1.2, color = "#D6604D") +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey40") +
  geom_smooth(method = "loess", se = TRUE, color = "#2166AC", linewidth = 0.8) +
  labs(x = "% missing", y = "Mean shift (post - pre)",
       title = "Imputation Shift vs Missingness (MNAR)") + thm

pdf(file.path(REPORT_DIR, "02_imputation_report.pdf"), width=12, height=10)
print(p_bench_main +
        inset_element(p_bench_inset, left=0.55, bottom=0.02, right=0.98, top=0.40))
print(p_hist / p_grp_dens +
        plot_annotation(title=sprintf("Post-Imputation Quality (%s)", best),
                        subtitle=sprintf("%d proteins x %d samples | %d values imputed (%.1f%%)",
                                         nrow(mat), ncol(mat), sum(was_na),
                                         sum(was_na)/length(mat)*100)))
print((p_mnar_a | p_mnar_b) +
        plot_annotation(title=sprintf("MNAR Post-Imputation Audit (%d proteins)", length(mnar_genes))))
dev.off()

cat(sprintf("MNAR shift: %.3f (median), range [%.3f, %.3f]\n",
            median(mnar_audit$shift), min(mnar_audit$shift), max(mnar_audit$shift)))
write_csv(mnar_audit, file.path(DATA_DIR, "mnar_imputation_audit.csv"))

###############################################################################
# 7: EXPORT
###############################################################################
stopifnot(identical(ann$gene, rownames(mat_imp)))
write_csv(bind_cols(ann, as_tibble(mat_imp)), file.path(DATA_DIR, "01_imputed.csv"))
write_csv(miss_class, file.path(DATA_DIR, "mar_mnar_classification.csv"))

info <- list(n_proteins=nrow(mat), n_samples=ncol(mat), pct_missing=pct_miss,
             n_mar=sum(miss_class$classification=="MAR"),
             n_mnar=sum(miss_class$classification=="MNAR"),
             n_complete=sum(miss_class$classification=="Complete"),
             best_method=best, best_nrmse=bench_sum$mean_nrmse[1])
writeLines(paste(names(info), info, sep=" = "),
           file.path(DATA_DIR, "imputation_summary.txt"))

# Save DAList with imputed data
dal <- readRDS("01_normalization/c_data/01_DAList_normalized.rds")
rownames(mat_imp) <- ann$uniprot_id
dal$data <- mat_imp
saveRDS(dal, file.path(DATA_DIR, "01_DAList_imputed.rds"))

cat(sprintf("\n=== YvO Imputation complete === %s | NRMSE %.4f ===\n",
            best, info$best_nrmse))
