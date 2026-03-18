# 00_run_plier.R — Fit PLIER on imputed proteomics data.
# Mao et al. 2019, Nat Methods 16:607-610 (PMID 31249421).
# Prior: Reactome + KEGG + Hallmark from msigdbr.
# Duplicate genes collapsed by highest mean abundance.

setwd(rprojroot::find_rstudio_root_file())

library(PLIER)
library(msigdbr)
library(readr)
library(dplyr)
library(tibble)

set.seed(42)

# --- Paths ---
DAT <- "04_Figures/PLIER_F09/c_data"
dir.create(DAT, recursive = TRUE, showWarnings = FALSE)

# --- Load expression data ---
imp_df <- read_csv("02_Imputation/c_data/01_imputed.csv", show_col_types = FALSE)
ann_cols <- c("uniprot_id", "protein", "gene", "description")
samp_names <- setdiff(names(imp_df), ann_cols)

expr_raw <- as.matrix(imp_df[, samp_names])
rownames(expr_raw) <- imp_df$gene
gene_means <- rowMeans(expr_raw, na.rm = TRUE)

# --- Resolve duplicate gene symbols ---
# Rule: keep row with highest mean abundance
dup_genes <- imp_df$gene[duplicated(imp_df$gene)]
dup_summary <- imp_df |>
  filter(gene %in% dup_genes) |>
  mutate(mean_abundance = gene_means[match(gene, imp_df$gene)]) |>
  select(uniprot_id, gene, protein, mean_abundance)

if (nrow(dup_summary) > 0) {
  # For each gene, keep the row with highest mean
  keep_idx <- imp_df |>
    mutate(row_idx = row_number(), mean_abund = gene_means) |>
    group_by(gene) |>
    slice_max(mean_abund, n = 1, with_ties = FALSE) |>
    ungroup() |>
    pull(row_idx)
  expr_mat <- expr_raw[keep_idx, ]
  rownames(expr_mat) <- imp_df$gene[keep_idx]
  cat(sprintf("Duplicate genes: %d symbols (%d rows removed)\n",
              length(unique(dup_genes)), nrow(expr_raw) - nrow(expr_mat)))
} else {
  expr_mat <- expr_raw
  cat("No duplicate gene symbols found\n")
}

write_csv(dup_summary, file.path(DAT, "06_gene_overlap_summary.csv"))

# --- Build prior matrix from msigdbr (v25+: collection/subcollection) ---
pathways <- bind_rows(
  msigdbr(species = "Homo sapiens", collection = "C2", subcollection = "CP:REACTOME"),
  msigdbr(species = "Homo sapiens", collection = "C2", subcollection = "CP:KEGG_MEDICUS"),
  msigdbr(species = "Homo sapiens", collection = "H")
)

genes_in_data <- rownames(expr_mat)
pathways_filt <- pathways |> filter(gene_symbol %in% genes_in_data)
all_sets <- unique(pathways_filt$gs_name)

C <- matrix(0L, nrow = length(genes_in_data), ncol = length(all_sets),
            dimnames = list(genes_in_data, all_sets))
for (i in seq_len(nrow(pathways_filt))) {
  g <- pathways_filt$gene_symbol[i]
  s <- pathways_filt$gs_name[i]
  C[g, s] <- 1L
}

# Filter pathways: require minGenes genes in common
# Use minGenes=5 for proteomics (fewer features than transcriptomics)
min_genes <- 5
pathway_sizes_pre <- colSums(C)
C <- C[, pathway_sizes_pre >= min_genes, drop = FALSE]
cat(sprintf("Prior matrix: %d genes x %d pathways (minGenes=%d)\n",
            nrow(C), ncol(C), min_genes))
cat(sprintf("Pathways dropped (< %d genes): %d\n",
            min_genes, sum(pathway_sizes_pre < min_genes)))

# --- Intersect genes ---
cm_genes <- PLIER::commonRows(C, expr_mat)
cat(sprintf("Common genes between data and prior: %d / %d\n",
            length(cm_genes), nrow(expr_mat)))

expr_sub <- expr_mat[cm_genes, ]
C_sub <- C[cm_genes, ]

# Drop any pathways that now have 0 genes
C_sub <- C_sub[, colSums(C_sub) >= min_genes, drop = FALSE]
cat(sprintf("Final prior: %d genes x %d pathways\n",
            nrow(C_sub), ncol(C_sub)))

# --- Estimate k ---
# num.pc() has a class() bug in R >= 4.0; use SVD elbow directly
expr_norm <- PLIER::rowNorm(expr_sub)
svd_res <- rsvd::rsvd(expr_norm, k = min(ncol(expr_norm) - 1, 50))
# Elbow: find where explained variance flattens (second derivative)
var_exp <- svd_res$d^2 / sum(svd_res$d^2)
diffs <- diff(var_exp)
diffs2 <- diff(diffs)
k_est <- which.min(diffs2[1:min(30, length(diffs2))]) + 1
k_use <- min(k_est * 2, ncol(expr_sub) - 1)
cat(sprintf("SVD elbow estimate: %d -> using k = %d (2x rule per Mao et al.)\n",
            k_est, k_use))

k_summary <- tibble(
  method = "num.pc (elbow) x2",
  num_pc_estimate = k_est,
  k_used = k_use,
  n_genes = nrow(expr_sub),
  n_samples = ncol(expr_sub),
  n_pathways = ncol(C_sub),
  rationale = "Mao et al. 2019: 'increase initial k by factor of 2'; robust to overestimation"
)
write_csv(k_summary, file.path(DAT, "07_k_selection_summary.csv"))

# --- Fit PLIER ---
cat("Fitting PLIER...\n")
plier_res <- PLIER(
  data     = expr_sub,
  priorMat = C_sub,
  k        = k_use,
  scale    = TRUE,
  trace    = TRUE,
  doCrossval = TRUE,
  minGenes = min_genes,
  seed     = 42
)

# --- Summary statistics ---
n_lv <- nrow(plier_res$B)
n_with_prior <- length(plier_res$withPrior)
cat(sprintf("\nPLIER fit: %d LVs, %d pathway-aligned (%.0f%%)\n",
            n_lv, n_with_prior, 100 * n_with_prior / n_lv))

# --- Export core outputs ---
# 01: full PLIER object
saveRDS(plier_res, file.path(DAT, "01_plier_fit.rds"))

# 02: LV scores (k x samples -> transpose to samples x LVs for convenience)
lv_scores <- as.data.frame(t(plier_res$B))
lv_scores$sample_id <- colnames(plier_res$B)
lv_scores <- lv_scores |> select(sample_id, everything())
colnames(lv_scores)[2:(n_lv + 1)] <- paste0("LV", seq_len(n_lv))
write_csv(lv_scores, file.path(DAT, "02_lv_scores.csv"))

# 03: gene loadings (genes x LVs)
gene_loadings <- as.data.frame(plier_res$Z)
colnames(gene_loadings) <- paste0("LV", seq_len(n_lv))
gene_loadings$gene <- rownames(plier_res$Z)
gene_loadings <- gene_loadings |> select(gene, everything())
write_csv(gene_loadings, file.path(DAT, "03_gene_loadings.csv"))

# 04: U matrix (pathway x LV)
u_mat <- as.data.frame(plier_res$U)
colnames(u_mat) <- paste0("LV", seq_len(n_lv))
u_mat$pathway <- rownames(plier_res$U)
u_mat <- u_mat |> select(pathway, everything())
write_csv(u_mat, file.path(DAT, "04_pathway_lv_u_matrix.csv"))

# 05: pathway-LV alignment summary (from cross-validation)
if (!is.null(plier_res$summary) && nrow(plier_res$summary) > 0) {
  align_summary <- plier_res$summary |>
    as_tibble()
  # Normalize column names (PLIER uses "LV index", "p-value")
  names(align_summary) <- gsub("LV index", "LV", names(align_summary))
  names(align_summary) <- gsub("p-value", "pval", names(align_summary))
  names(align_summary) <- gsub(" ", "_", names(align_summary))
  align_summary <- align_summary |> arrange(FDR, desc(AUC))
  write_csv(align_summary, file.path(DAT, "05_lv_pathway_alignment_summary.csv"))
  cat(sprintf("Pathway-LV associations (FDR<0.05): %d\n",
              sum(align_summary$FDR < 0.05, na.rm = TRUE)))
} else {
  # Fallback: extract from U and Uauc manually
  sig_pairs <- which(plier_res$U > 0, arr.ind = TRUE)
  if (nrow(sig_pairs) > 0) {
    align_summary <- tibble(
      pathway = rownames(plier_res$U)[sig_pairs[, 1]],
      LV = sig_pairs[, 2],
      U_coef = plier_res$U[sig_pairs],
      AUC = if (!is.null(plier_res$Uauc)) plier_res$Uauc[sig_pairs] else NA_real_,
      pval = if (!is.null(plier_res$Up)) plier_res$Up[sig_pairs] else NA_real_
    ) |>
      mutate(FDR = p.adjust(pval, method = "BH")) |>
      arrange(FDR, desc(AUC))
  } else {
    align_summary <- tibble(pathway = character(), LV = integer(),
                            U_coef = numeric(), AUC = numeric(),
                            pval = numeric(), FDR = numeric())
  }
  write_csv(align_summary, file.path(DAT, "05_lv_pathway_alignment_summary.csv"))
}

cat(sprintf("\n=== PLIER complete ===\n"))
cat(sprintf("  LVs: %d\n  Pathway-aligned: %d\n", n_lv, n_with_prior))
cat(sprintf("  Genes used: %d\n  Pathways in prior: %d\n",
            nrow(expr_sub), ncol(C_sub)))
cat(sprintf("  Outputs saved to: %s/\n", DAT))
