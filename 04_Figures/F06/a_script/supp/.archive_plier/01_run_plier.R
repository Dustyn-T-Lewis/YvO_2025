# Supplementary: PLIER — Fit PLIER on imputed proteomics data.
# Mao et al. 2019, Nat Methods 16:607-610 (PMID 31249421).
# Prior: Reactome + KEGG + Hallmark from msigdbr.
# Duplicate genes collapsed by highest mean abundance.

setwd(rprojroot::find_rstudio_root_file())

library(PLIER)
library(msigdbr)
# Guard: msigdbr v25+ uses collection/subcollection (not category/subcategory)
if (packageVersion("msigdbr") < "10.0" &&
    !"collection" %in% names(formals(msigdbr::msigdbr))) {
  stop("msigdbr version too old — requires v7.5.1+ with collection/subcollection API.",
       call. = FALSE)
}
library(readr)
library(dplyr)
library(tibble)

set.seed(42)

# --- Paths
DAT <- "04_Figures/F06/c_data/supp/plier"
dir.create(DAT, recursive = TRUE, showWarnings = FALSE)

# --- Load expression data
imp_df <- read_csv("02_Imputation/c_data/01_imputed.csv", show_col_types = FALSE)
ann_cols <- c("uniprot_id", "protein", "gene", "description")
samp_names <- setdiff(names(imp_df), ann_cols)

expr_raw <- as.matrix(imp_df[, samp_names])
rownames(expr_raw) <- imp_df$gene
gene_means <- rowMeans(expr_raw, na.rm = TRUE)

# --- Resolve duplicate gene symbols
# Rule: keep row with highest mean abundance
dup_genes <- imp_df$gene[duplicated(imp_df$gene)]
dup_summary <- imp_df |>
  mutate(mean_abundance = gene_means) |>
  filter(gene %in% dup_genes) |>
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
  message(sprintf("Duplicate genes: %d symbols (%d rows removed)",
              length(unique(dup_genes)), nrow(expr_raw) - nrow(expr_mat)))
} else {
  expr_mat <- expr_raw
  message("No duplicate gene symbols found")
}

write_csv(dup_summary, file.path(DAT, "06_gene_overlap_summary.csv"))

# --- Build prior matrix from msigdbr (v25+: collection/subcollection)
pathways <- bind_rows(
  msigdbr(species = "Homo sapiens", collection = "C2", subcollection = "CP:REACTOME"),
  msigdbr(species = "Homo sapiens", collection = "C2", subcollection = "CP:KEGG_MEDICUS"),
  msigdbr(species = "Homo sapiens", collection = "H")
)

genes_in_data <- rownames(expr_mat)
pathways_filt <- pathways |> filter(gene_symbol %in% genes_in_data)
all_sets <- unique(pathways_filt$gs_name)

# Build binary prior matrix via matrix indexing (vectorized)
C <- matrix(0L, nrow = length(genes_in_data), ncol = length(all_sets),
            dimnames = list(genes_in_data, all_sets))
idx <- cbind(match(pathways_filt$gene_symbol, genes_in_data),
             match(pathways_filt$gs_name, all_sets))
C[idx] <- 1L

# Filter pathways: require minGenes genes in common
# minGenes=5 for proteomics: ~2K proteins vs ~20K genes means fewer pathway members
# survive intersection. Default 10 drops too many pathways. Lavin et al. 2020
# (PMID 32625117) used 5884 genes; our 2K proteins warrant a lower floor.
min_genes <- 5
pathway_sizes_pre <- colSums(C)
C <- C[, pathway_sizes_pre >= min_genes, drop = FALSE]
message(sprintf("Prior matrix: %d genes x %d pathways (minGenes=%d, %d dropped)",
                nrow(C), ncol(C), min_genes, sum(pathway_sizes_pre < min_genes)))

# --- Intersect genes
cm_genes <- PLIER::commonRows(C, expr_mat)
message(sprintf("Common genes between data and prior: %d / %d",
                length(cm_genes), nrow(expr_mat)))

expr_sub <- expr_mat[cm_genes, ]
C_sub <- C[cm_genes, ]

# Drop any pathways that now have 0 genes
C_sub <- C_sub[, colSums(C_sub) >= min_genes, drop = FALSE]
message(sprintf("Final prior: %d genes x %d pathways", nrow(C_sub), ncol(C_sub)))

# --- Estimate k
# num.pc() has a class() bug in R >= 4.0; use SVD elbow directly
expr_norm <- PLIER::rowNorm(expr_sub)
svd_res <- rsvd::rsvd(expr_norm, k = min(ncol(expr_norm) - 1, 50))
# Elbow: find where explained variance flattens (second derivative)
var_exp <- svd_res$d^2 / sum(svd_res$d^2)
diffs <- diff(var_exp)
diffs2 <- diff(diffs)
k_est <- which.min(diffs2[1:min(30, length(diffs2))]) + 1
k_use <- min(k_est * 2, ncol(expr_sub) - 1)
message(sprintf("SVD elbow estimate: %d -> using k = %d (2x rule per Mao et al.)",
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

message("Fitting PLIER...")
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

n_lv <- nrow(plier_res$B)
n_with_prior <- length(plier_res$withPrior)
message(sprintf("PLIER fit: %d LVs, %d pathway-aligned (%.0f%%)",
                n_lv, n_with_prior, 100 * n_with_prior / n_lv))

# Export core outputs ----
saveRDS(plier_res, file.path(DAT, "01_plier_fit.rds"))

# Helper: convert named-row matrix to tidy CSV with ID column first
export_matrix <- function(mat, id_col, lv_names, path) {
  df <- as.data.frame(mat)
  colnames(df) <- lv_names
  df[[id_col]] <- rownames(mat)
  df <- df |> select(all_of(id_col), everything())
  write_csv(df, path)
}

lv_names <- paste0("LV", seq_len(n_lv))
export_matrix(t(plier_res$B), "sample_id", lv_names, file.path(DAT, "02_lv_scores.csv"))
export_matrix(plier_res$Z,    "gene",      lv_names, file.path(DAT, "03_gene_loadings.csv"))
export_matrix(plier_res$U,    "pathway",   lv_names, file.path(DAT, "04_pathway_lv_u_matrix.csv"))

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
  message(sprintf("Pathway-LV associations (FDR<0.05): %d",
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

message(sprintf("PLIER complete: %d LVs, %d pathway-aligned, %d genes, %d pathways",
                n_lv, n_with_prior, nrow(expr_sub), ncol(C_sub)))
