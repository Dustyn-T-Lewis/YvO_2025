#!/usr/bin/env Rscript
###############################################################################
# YvO_imputation.R — Missingness Diagnostics, Benchmarking & Imputation
#
# 1. Characterize missingness (overall, per-group, per-protein)
# 2. Classify MAR vs MNAR (msImpute)
# 3. Benchmark imputation methods (NRMSE on held-out observed values)
# 4. Apply best method, export imputed matrix + diagnostics
#
# Input:  01_normalized.csv  (log2 cycloess-normalized, protein-level)
# Output: 01_imputed.csv + 4 diagnostic PDFs + summary files
###############################################################################

cat("=== YvO Imputation Pipeline ===\n\n")

if (!requireNamespace("pacman", quietly = TRUE)) install.packages("pacman")
pacman::p_load(
  MsCoreUtils, msImpute, pcaMethods, imputeLCMD,
  ggplot2, patchwork, dplyr, tidyr, tibble, readr, stringr, scales
)

# ─── Paths ───────────────────────────────────────────────────────────────────
base_dir   <- normalizePath(file.path(dirname(getwd()), ".."), mustWork = TRUE)
INPUT_CSV  <- file.path(base_dir, "01_normalization", "c_data", "01_normalized.csv")
REPORT_DIR <- file.path(base_dir, "02_Imputation", "b_reports")
DATA_DIR   <- file.path(base_dir, "02_Imputation", "c_data")
dir.create(REPORT_DIR, showWarnings = FALSE, recursive = TRUE)
dir.create(DATA_DIR,   showWarnings = FALSE, recursive = TRUE)

# ─── Shared aesthetics ──────────────────────────────────────────────────────
pal_gt   <- c(Young_Pre="#92C5DE", Young_Post="#2166AC",
              Old_Pre="#F4A582",   Old_Post="#B2182B")
pal_mar  <- c(MAR="#4393C3", MNAR="#D6604D")
pal_mtyp <- c(MNAR="#D6604D", MAR="#4393C3", Hybrid="#5AAE61")
thm      <- theme_minimal(base_size = 11)
thm_sm   <- theme_minimal(base_size = 8)

# ─── Helper: dispatch imputation ────────────────────────────────────────────
run_impute <- function(m, mat, randna) {
  if (m$method == "mixed")
    MsCoreUtils::impute_matrix(mat, method="mixed",
                               randna=randna, mar=m$mar, mnar=m$mnar)
  else
    MsCoreUtils::impute_matrix(mat, method=m$method)
}

###############################################################################
# 1: LOAD DATA
###############################################################################
cat(">> 1 — Loading data\n")
df  <- read_csv(INPUT_CSV, show_col_types = FALSE)
ann <- df[, 1:4]
mat <- as.matrix(df[, -(1:4)])
rownames(mat) <- df$gene
cat(sprintf("   %d proteins × %d samples\n", nrow(mat), ncol(mat)))

meta <- tibble(Col_ID = colnames(mat)) %>%
  mutate(Group     = if_else(str_detect(Col_ID, "^(Y_|YP_)"), "Young", "Old"),
         Timepoint = if_else(str_detect(Col_ID, "_Pre$"), "Pre", "Post"),
         Group_Time = paste(Group, Timepoint, sep = "_"))
print(count(meta, Group, Timepoint))

###############################################################################
# 2: MISSINGNESS DIAGNOSTICS
###############################################################################
cat("\n>> 2 — Missingness diagnostics\n")

prot_miss <- rowSums(is.na(mat))
prot_pct  <- prot_miss / ncol(mat) * 100
obs_means <- rowMeans(mat, na.rm = TRUE)
pct_miss  <- round(sum(is.na(mat)) / length(mat) * 100, 2)

cat(sprintf("   Missing: %d / %d (%.2f%%) | Complete proteins: %d\n",
            sum(is.na(mat)), length(mat), pct_miss, sum(prot_miss == 0)))

miss_by_group <- sapply(unique(meta$Group_Time), function(g) {
  cols <- meta$Col_ID[meta$Group_Time == g]
  rowSums(is.na(mat[, cols, drop=FALSE])) / length(cols) * 100
})

pdf(file.path(REPORT_DIR, "01_missingness_diagnostics.pdf"), width=12, height=14)

p2a <- ggplot(tibble(x = prot_pct), aes(x)) +
  geom_histogram(binwidth=2, fill="#4393C3", color="white", alpha=0.8) +
  geom_vline(xintercept=c(20,50), linetype="dashed", color="red", alpha=0.6) +
  annotate("text", x=22, y=Inf, vjust=2, hjust=0, size=3,
           label=sprintf("Complete: %d\n1-20%%: %d\n20-50%%: %d\n>50%%: %d",
                         sum(prot_pct==0), sum(prot_pct>0 & prot_pct<=20),
                         sum(prot_pct>20 & prot_pct<=50), sum(prot_pct>50))) +
  labs(x="% missing per protein", y="Count", title="A: Per-Protein Missingness") + thm

top_idx <- order(prot_pct, decreasing=TRUE)[1:min(50, sum(prot_pct>0))]
p2b <- as_tibble(miss_by_group[top_idx,], rownames="gene") %>%
  pivot_longer(-gene, names_to="Group_Time", values_to="pct") %>%
  mutate(gene = factor(gene, levels=rev(rownames(mat)[top_idx]))) %>%
  ggplot(aes(Group_Time, gene, fill=pct)) +
  geom_tile(color="white", linewidth=0.3) +
  scale_fill_gradient2(low="white", mid="#FDDBC7", high="#B2182B", midpoint=50, name="% Miss") +
  labs(x=NULL, y=NULL, title="B: Per-Group Missingness (top 50)") +
  theme_minimal(base_size=9) +
  theme(axis.text.y=element_text(size=5), axis.text.x=element_text(angle=45, hjust=1))

p2c <- ggplot(tibble(int=obs_means, miss=prot_pct), aes(int, miss)) +
  geom_point(alpha=0.3, size=0.8, color="#4393C3") +
  geom_smooth(method="loess", se=TRUE, color="#B2182B", linewidth=0.8) +
  labs(x="Mean log2 intensity", y="% missing", title="C: Missingness vs Abundance") + thm

p2d <- tibble(Col_ID=colnames(mat), pct=colSums(is.na(mat))/nrow(mat)*100) %>%
  left_join(meta, by="Col_ID") %>%
  ggplot(aes(reorder(Col_ID, pct), pct, fill=Group_Time)) +
  geom_col(alpha=0.85) + scale_fill_manual(values=pal_gt) + coord_flip() +
  labs(x=NULL, y="% missing proteins", title="D: Per-Sample Missingness") +
  thm_sm + theme(axis.text.y=element_text(size=4))

print((p2a | p2c) / (p2b | p2d) + plot_annotation(
  title="Pre-Imputation Missingness Diagnostics",
  subtitle=sprintf("%d proteins × %d samples | %.2f%% missing", nrow(mat), ncol(mat), pct_miss)))
dev.off()
cat("   Saved: 01_missingness_diagnostics.pdf\n")

###############################################################################
# 3: MAR / MNAR CLASSIFICATION
###############################################################################
cat("\n>> 3 — MAR/MNAR classification\n")

has_na <- which(prot_miss > 0 & prot_miss < ncol(mat))
cat(sprintf("   Classifying %d proteins with missingness\n", length(has_na)))

miss_class <- tibble(gene=rownames(mat), n_miss=prot_miss,
                     pct_miss=prot_pct, mean_intensity=obs_means)

mar_result <- tryCatch({
  mar_feat <- msImpute::selectFeatures(mat[has_na,], method="ebm")
  list(success = TRUE, features = mar_feat)
}, error = function(e) {
  cat("   selectFeatures failed, using intensity heuristic\n")
  list(success = FALSE)
})

if (mar_result$success) {
  miss_class <- miss_class %>%
    mutate(classification = case_when(
      n_miss==0 ~ "Complete", gene %in% mar_result$features ~ "MAR", TRUE ~ "MNAR"))
} else {
  med <- median(obs_means, na.rm=TRUE)
  miss_class <- miss_class %>%
    mutate(classification = case_when(
      n_miss==0 ~ "Complete",
      pct_miss>30 & mean_intensity<med ~ "MNAR", pct_miss>50 ~ "MNAR",
      TRUE ~ "MAR"))
}
print(count(miss_class, classification))

pdf(file.path(REPORT_DIR, "02_mar_mnar_classification.pdf"), width=10, height=8)
mc <- miss_class %>% filter(classification != "Complete")

p3a <- ggplot(mc, aes(mean_intensity, pct_miss, color=classification)) +
  geom_point(alpha=0.5, size=1.2) + scale_color_manual(values=pal_mar) +
  labs(x="Mean log2 intensity", y="% missing", title="A: MAR vs MNAR") +
  thm + theme(legend.position="bottom")

p3b <- ggplot(mc, aes(mean_intensity, fill=classification)) +
  geom_density(alpha=0.5) + scale_fill_manual(values=pal_mar) +
  labs(x="Mean log2 intensity", title="B: Intensity by Missingness Type") +
  thm + theme(legend.position="bottom")

p3c <- ggplot(mc, aes(classification, pct_miss, fill=classification)) +
  geom_boxplot(alpha=0.6, outlier.size=0.5) + scale_fill_manual(values=pal_mar) +
  labs(x=NULL, y="% missing", title="C: Missingness by Type") +
  thm + theme(legend.position="none")

print((p3a | p3b) / (p3c | plot_spacer()) +
        plot_annotation(title="MAR / MNAR Classification"))
dev.off()
cat("   Saved: 02_mar_mnar_classification.pdf\n")

###############################################################################
# 4: BENCHMARKING
###############################################################################
cat("\n>> 4 — Imputation benchmarking\n")
set.seed(42); N_ITER <- 20

METHODS <- list(
  MinProb=list(method="MinProb"), MinDet=list(method="MinDet"),
  QRILC=list(method="QRILC"),    zero=list(method="zero"),
  knn=list(method="knn"),         bpca=list(method="bpca"),
  ppca=list(method="ppca"),       nipals=list(method="nipals"),
  MLE=list(method="MLE"),
  mix_knn_MinProb =list(method="mixed", mar="knn",  mnar="MinProb"),
  mix_knn_MinDet  =list(method="mixed", mar="knn",  mnar="MinDet"),
  mix_knn_QRILC   =list(method="mixed", mar="knn",  mnar="QRILC"),
  mix_bpca_MinProb=list(method="mixed", mar="bpca", mnar="MinProb"),
  mix_bpca_QRILC  =list(method="mixed", mar="bpca", mnar="QRILC")
)

randna <- setNames(miss_class$classification == "MNAR", miss_class$gene)
nrmse  <- function(t, i) sqrt(mean((t - i)^2)) / sd(t)

cat(sprintf("   %d methods × %d iterations\n", length(METHODS), N_ITER))
res <- vector("list", length(METHODS) * N_ITER); k <- 0L
fail_res <- list(); k_fail <- 0L

for (iter in seq_len(N_ITER)) {
  cat(sprintf("   Iter %d/%d\n", iter, N_ITER))
  obs_idx  <- which(!is.na(mat))
  mask_idx <- sample(obs_idx, round(length(obs_idx)*0.10))
  true_v   <- mat[mask_idx]
  mm <- mat; mm[mask_idx] <- NA

  # Correct randna for masked positions: treat masked MNAR proteins as MAR
  masked_rows <- unique((mask_idx - 1) %% nrow(mat) + 1)
  randna_iter <- randna
  randna_iter[masked_rows] <- FALSE

  for (nm in names(METHODS)) {
    result <- tryCatch(
      list(val = run_impute(METHODS[[nm]], mm, randna_iter)),
      error = function(e) list(val = NULL, err = conditionMessage(e)))
    if (is.null(result$val)) {
      cat(sprintf("      %s failed: %s\n", nm, result$err))
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
  cat(sprintf("   %d method failures logged to benchmark_failures.csv\n", k_fail))
}

bench_df <- bind_rows(res)
bench_sum <- bench_df %>%
  group_by(method) %>%
  summarise(mean_nrmse=mean(nrmse), sd_nrmse=sd(nrmse),
            median_nrmse=median(nrmse), .groups="drop") %>%
  arrange(mean_nrmse)
print(bench_sum)

best <- bench_sum$method[1]
cat(sprintf("\n   ★ Best: %s (mean NRMSE = %.4f)\n", best, bench_sum$mean_nrmse[1]))

mtype <- tibble(method=names(METHODS)) %>%
  mutate(type = case_when(
    method %in% c("MinProb","MinDet","QRILC","zero") ~ "MNAR",
    method %in% c("knn","bpca","ppca","nipals","MLE") ~ "MAR",
    TRUE ~ "Hybrid"))

pdf(file.path(REPORT_DIR, "03_imputation_benchmark.pdf"), width=10, height=7)
print(
  bench_df %>% filter(nrmse < 10) %>% left_join(mtype, by="method") %>%
    ggplot(aes(reorder(method, nrmse, FUN=median), nrmse, fill=type)) +
    geom_boxplot(alpha=0.6, outlier.size=1) +
    geom_jitter(width=0.15, size=1.2, alpha=0.5) +
    geom_hline(yintercept=bench_sum$mean_nrmse[1], linetype="dashed",
               color="#B2182B", alpha=0.6) +
    scale_fill_manual(values=pal_mtyp, name="Type") + coord_flip() +
    labs(x=NULL, y="NRMSE (lower = better)", title="Imputation Benchmark",
         subtitle=sprintf("%d iter × 10%% masked | Best: %s (%.4f)",
                          N_ITER, best, bench_sum$mean_nrmse[1])) + thm
)
dev.off()
cat("   Saved: 03_imputation_benchmark.pdf\n")
write_csv(bench_sum, file.path(DATA_DIR, "benchmark_summary.csv"))

###############################################################################
# 5: APPLY BEST METHOD
###############################################################################
cat("\n>> 5 — Applying best method\n")
set.seed(42)
mat_imp <- run_impute(METHODS[[best]], mat, randna)
cat(sprintf("   Remaining NAs: %d\n", sum(is.na(mat_imp))))

###############################################################################
# 6: POST-IMPUTATION DIAGNOSTICS
###############################################################################
cat("\n>> 6 — Post-imputation diagnostics\n")

was_na <- is.na(mat)
mat_med <- mat  # median-fill baseline for PCA comparison
for (j in seq_len(ncol(mat_med)))
  mat_med[is.na(mat_med[,j]), j] <- median(mat_med[,j], na.rm=TRUE)

pdf(file.path(REPORT_DIR, "04_pre_vs_post_imputation.pdf"), width=14, height=20)

# 6A: Global density
dens_df <- bind_rows(
  tibble(value=as.vector(mat),       stage="Observed"),
  tibble(value=as.vector(mat_imp),   stage="All (post)"),
  tibble(value=mat_imp[was_na],      stage="Imputed only"))

p6a <- ggplot(dens_df, aes(value, color=stage, linetype=stage)) +
  geom_density(linewidth=0.8) +
  scale_color_manual(values=c(Observed="#2166AC", `All (post)`="#4DAF4A",
                               `Imputed only`="#B2182B")) +
  scale_linetype_manual(values=c(Observed="solid", `All (post)`="solid",
                                  `Imputed only`="dashed")) +
  labs(x="log2 intensity", y="Density", title="A: Global Intensity Distributions") +
  thm + theme(legend.position="bottom", legend.title=element_blank())

# 6B: Observed vs imputed histogram (counts × intensity)
hist_df <- bind_rows(
  tibble(value=as.vector(mat[!was_na]), source="Observed"),
  tibble(value=mat_imp[was_na],         source="Imputed"))

p6b_hist <- ggplot(hist_df, aes(value, fill=source)) +
  geom_histogram(binwidth=0.2, alpha=0.6, position="identity", color="white",
                 linewidth=0.1) +
  scale_fill_manual(values=c(Observed="#2166AC", Imputed="#B2182B")) +
  labs(x="log2 intensity", y="Count",
       title="B: Observed vs Imputed Counts",
       subtitle=sprintf("%s observed | %s imputed",
                        scales::comma(sum(!was_na)), scales::comma(sum(was_na)))) +
  thm + theme(legend.position="bottom", legend.title=element_blank())

# 6C: Per-group density
grp_df <- do.call(rbind, lapply(unique(meta$Group_Time), function(g) {
  cols <- meta$Col_ID[meta$Group_Time==g]
  sm <- mat[, cols, drop=FALSE]; si <- mat_imp[, cols, drop=FALSE]; sn <- was_na[, cols, drop=FALSE]
  bind_rows(
    tibble(value=as.vector(sm[!is.na(sm)]), stage="Observed", Group_Time=g),
    tibble(value=as.vector(si[sn]),         stage="Imputed",  Group_Time=g))
}))

p6c <- ggplot(grp_df, aes(value, color=stage, linetype=stage)) +
  geom_density(linewidth=0.7) + facet_wrap(~Group_Time, scales="free_y") +
  scale_color_manual(values=c(Observed="#2166AC", Imputed="#B2182B")) +
  scale_linetype_manual(values=c(Observed="solid", Imputed="dashed")) +
  labs(x="log2 intensity", title="C: Observed vs Imputed by Group") +
  theme_minimal(base_size=10) +
  theme(legend.position="bottom", legend.title=element_blank())

# 6D: Per-sample boxplots
samp_df <- bind_rows(
  tibble(Col_ID=rep(colnames(mat), each=nrow(mat)),
         value=as.vector(mat), stage="Pre") %>% filter(!is.na(value)),
  tibble(Col_ID=rep(colnames(mat_imp), each=nrow(mat_imp)),
         value=as.vector(mat_imp), stage="Post")
) %>% left_join(meta, by="Col_ID")

p6d <- ggplot(samp_df, aes(Col_ID, value, fill=stage)) +
  geom_boxplot(outlier.size=0.2, alpha=0.6, linewidth=0.3) +
  scale_fill_manual(values=c(Pre="#92C5DE", Post="#F4A582")) + coord_flip() +
  labs(x=NULL, y="log2 intensity", title="D: Per-Sample Distributions (Pre vs Post)") +
  thm_sm + theme(axis.text.y=element_text(size=4), legend.position="bottom")

# 6E: Per-protein mean scatter
prot_df <- tibble(pre=rowMeans(mat, na.rm=TRUE), post=rowMeans(mat_imp),
                  had_na=rowSums(was_na) > 0)

p6e <- ggplot(prot_df, aes(pre, post, color=had_na)) +
  geom_point(alpha=0.4, size=0.8) +
  geom_abline(slope=1, intercept=0, linetype="dashed", color="grey40") +
  scale_color_manual(values=c("FALSE"="grey70","TRUE"="#B2182B"),
                     labels=c("Complete","Had missing"), name=NULL) +
  labs(x="Pre mean log2", y="Post mean log2", title="E: Per-Protein Mean Shift") +
  thm + theme(legend.position="bottom")

print(p6a / p6b_hist / p6c / p6d / p6e + plot_annotation(
  title=sprintf("Pre vs Post Imputation — %s", best),
  subtitle=sprintf("%d proteins × %d samples | %d imputed (%.1f%%)",
                   nrow(mat), ncol(mat), sum(was_na), sum(was_na)/length(mat)*100)))

# ─── Page 2: PCA + correlation ───────────────────────────────────────────────
pca_pre  <- prcomp(t(mat_med), center=TRUE, scale.=TRUE)
pca_post <- prcomp(t(mat_imp), center=TRUE, scale.=TRUE)
vp <- round(summary(pca_pre)$importance[2,1:2]*100, 1)
vi <- round(summary(pca_post)$importance[2,1:2]*100, 1)

pca_df <- bind_rows(
  as_tibble(pca_pre$x[,1:2])  %>% mutate(Col_ID=colnames(mat), stage="Pre"),
  as_tibble(pca_post$x[,1:2]) %>% mutate(Col_ID=colnames(mat), stage="Post")
) %>% left_join(meta, by="Col_ID")

make_pca_plot <- function(stg, v, ttl) {
  ggplot(filter(pca_df, stage==stg), aes(PC1, PC2, color=Group_Time)) +
    geom_point(size=2.5, alpha=0.8) + stat_ellipse(level=0.68, linewidth=0.6) +
    scale_color_manual(values=pal_gt) +
    labs(x=sprintf("PC1 (%.1f%%)", v[1]), y=sprintf("PC2 (%.1f%%)", v[2]), title=ttl) +
    theme_minimal(base_size=10) + theme(legend.position="bottom")
}

cor_diff <- cor(mat_imp) - cor(mat_med)
p_cor <- as_tibble(cor_diff, rownames="S1") %>%
  pivot_longer(-S1, names_to="S2", values_to="d") %>%
  ggplot(aes(S1, S2, fill=d)) + geom_tile() +
  scale_fill_gradient2(low="#2166AC", mid="white", high="#B2182B", name="Δ cor") +
  labs(title="H: Δ Sample Correlation (Post − Pre)") +
  thm_sm + theme(axis.text=element_text(size=3),
                 axis.text.x=element_text(angle=90, hjust=1), axis.title=element_blank())

n_pcs <- min(10, ncol(mat)-1)
p_var <- tibble(PC=paste0("PC",seq_len(n_pcs)),
                Pre=summary(pca_pre)$importance[2,seq_len(n_pcs)]*100,
                Post=summary(pca_post)$importance[2,seq_len(n_pcs)]*100) %>%
  pivot_longer(-PC, names_to="stage", values_to="v") %>%
  mutate(PC=factor(PC, levels=paste0("PC",seq_len(n_pcs)))) %>%
  ggplot(aes(PC, v, fill=stage)) +
  geom_col(position="dodge", alpha=0.7) +
  scale_fill_manual(values=c(Pre="#92C5DE", Post="#F4A582")) +
  labs(x=NULL, y="% variance", title="I: Variance Explained (Pre vs Post)") +
  theme_minimal(base_size=10) + theme(legend.position="bottom")

print(
  (make_pca_plot("Pre", vp, "F: PCA Pre (median fill)") |
   make_pca_plot("Post", vi, sprintf("G: PCA Post (%s)", best))) /
  (p_cor | p_var) + plot_annotation(title="PCA & Correlation Diagnostics"))

dev.off()
cat("   Saved: 04_pre_vs_post_imputation.pdf\n")

###############################################################################
# 6b: MNAR POST-IMPUTATION AUDIT
###############################################################################
cat("\n>> 6b — MNAR post-imputation audit\n")

mnar_genes <- miss_class$gene[miss_class$classification == "MNAR"]
mnar_audit <- tibble(
  gene = mnar_genes,
  pre_mean  = rowMeans(mat[mnar_genes, ], na.rm = TRUE),
  post_mean = rowMeans(mat_imp[mnar_genes, ]),
  pct_miss  = prot_pct[mnar_genes],
  shift     = rowMeans(mat_imp[mnar_genes, ]) - rowMeans(mat[mnar_genes, ], na.rm = TRUE),
  imputation_reliable = prot_pct[mnar_genes] < 50
)

pdf(file.path(REPORT_DIR, "05_mnar_imputation_audit.pdf"), width = 10, height = 8)
p_mnar_a <- ggplot(mnar_audit, aes(pre_mean, post_mean, color = pct_miss)) +
  geom_point(alpha = 0.7, size = 1.5) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey40") +
  scale_color_gradient(low = "#FDDBC7", high = "#B2182B", name = "% missing") +
  labs(x = "Pre-imputation mean log2", y = "Post-imputation mean log2",
       title = "A: MNAR Protein Mean Shift After Imputation") + thm

p_mnar_b <- ggplot(mnar_audit, aes(pct_miss, shift)) +
  geom_point(alpha = 0.6, size = 1.2, color = "#D6604D") +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey40") +
  geom_smooth(method = "loess", se = TRUE, color = "#2166AC", linewidth = 0.8) +
  labs(x = "% missing", y = "Mean shift (post - pre)",
       title = "B: Imputation Shift vs Missingness for MNAR Proteins") + thm

print((p_mnar_a | p_mnar_b) + plot_annotation(
  title = sprintf("MNAR Post-Imputation Audit (%d proteins, %s method)", length(mnar_genes), best)))
dev.off()
cat(sprintf("   MNAR mean shift: %.3f (median), range [%.3f, %.3f]\n",
            median(mnar_audit$shift), min(mnar_audit$shift), max(mnar_audit$shift)))
cat("   Saved: 05_mnar_imputation_audit.pdf\n")

write_csv(mnar_audit, file.path(DATA_DIR, "mnar_imputation_audit.csv"))

###############################################################################
# 7: EXPORT
###############################################################################
cat("\n>> 7 — Export\n")
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

# ─── Save DAList with imputed data ──────────────────────────────────────────
dal <- readRDS(file.path(base_dir, "01_normalization", "c_data", "01_DAList_normalized.rds"))
rownames(mat_imp) <- ann$uniprot_id
dal$data <- mat_imp
saveRDS(dal, file.path(DATA_DIR, "01_DAList_imputed.rds"))

cat(sprintf("\n=== Done === %s | NRMSE %.4f | %.2f%% missing ===\n",
            best, info$best_nrmse, pct_miss))
