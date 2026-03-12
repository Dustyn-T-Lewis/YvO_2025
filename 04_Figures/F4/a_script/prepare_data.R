# F4 prepare_data.R: Mfuzz FCM clustering + ORA enrichment cache
# Run BEFORE panel scripts. 50-start FCM (k=4), bootstrap ARI,
# per-cluster ORA (Hallmark + GO:BP + rrvgo), theme assignment.

setwd(rprojroot::find_rstudio_root_file())
source("04_Figures/shared/style.R")

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(tibble)
  library(Biobase)
  library(e1071)
  library(clusterProfiler)
  library(msigdbr)
  library(org.Hs.eg.db)
  library(rrvgo)
  library(GOSemSim)
})

# Load Mfuzz core functions without tcltk/X11 dependency
.mfuzz_env <- new.env()
suppressWarnings(
  lazyLoad(file.path(find.package("Mfuzz"), "R", "Mfuzz"), envir = .mfuzz_env)
)
mestimate   <- .mfuzz_env$mestimate
mfuzz       <- .mfuzz_env$mfuzz
standardise <- .mfuzz_env$standardise

set.seed(42)

DAT <- "04_Figures/F4/c_data"
RPT <- "04_Figures/F4/b_reports"
dir.create(DAT, recursive = TRUE, showWarnings = FALSE)
dir.create(RPT, recursive = TRUE, showWarnings = FALSE)

optimal_k   <- 4
CORE_THRESH <- 0.5
N_STARTS    <- 50
N_BOOT      <- 100

imputed <- read_csv("02_Imputation/c_data/01_imputed.csv", show_col_types = FALSE)
annot_cols  <- c("uniprot_id", "protein", "gene", "description")
sample_cols <- setdiff(colnames(imputed), annot_cols)

sample_meta <- tibble(sample = sample_cols) |>
  mutate(
    age     = ifelse(str_detect(sample, "^(O_|OP_)"), "Old", "Young"),
    time    = ifelse(str_detect(sample, "_Post$"), "Post", "Pre"),
    subject = str_remove(sample, "_(Pre|Post)$")
  )
write_csv(sample_meta, file.path(DAT, "sample_meta.csv"))

pre_subjects  <- sample_meta$subject[sample_meta$time == "Pre"]
post_subjects <- sample_meta$subject[sample_meta$time == "Post"]
paired_subjects <- intersect(pre_subjects, post_subjects)

young_subjects <- sort(paired_subjects[str_detect(paired_subjects, "^(Y_|YP_)")])
old_subjects   <- sort(paired_subjects[str_detect(paired_subjects, "^(O_|OP_)")])
ordered_subjects <- c(young_subjects, old_subjects)

gene_ids <- imputed$gene
uniprot_ids <- imputed$uniprot_id

na_mask <- is.na(gene_ids) | gene_ids == ""
if (any(na_mask)) gene_ids[na_mask] <- uniprot_ids[na_mask]

dup_counts <- table(gene_ids)
dup_genes  <- names(dup_counts[dup_counts > 1])
if (length(dup_genes) > 0) {
  for (dg in dup_genes) {
    idx <- which(gene_ids == dg)
    for (j in seq_along(idx)[-1]) {
      gene_ids[idx[j]] <- paste0(dg, "_", j)
    }
  }
}

abund_mat <- imputed |>
  dplyr::select(all_of(sample_cols)) |>
  as.matrix()
rownames(abund_mat) <- gene_ids

delta_mat <- matrix(NA_real_,
                    nrow = nrow(abund_mat),
                    ncol = length(ordered_subjects),
                    dimnames = list(gene_ids, ordered_subjects))

for (subj in ordered_subjects) {
  pre_col  <- paste0(subj, "_Pre")
  post_col <- paste0(subj, "_Post")
  delta_mat[, subj] <- abund_mat[, post_col] - abund_mat[, pre_col]
}

# Z-score standardize per protein
delta_z <- t(scale(t(delta_mat)))

zero_var <- apply(delta_mat, 1, sd, na.rm = TRUE) == 0
if (any(zero_var)) delta_z <- delta_z[!zero_var, , drop = FALSE]

nan_rows <- apply(delta_z, 1, function(x) any(is.nan(x)))
if (any(nan_rows)) delta_z <- delta_z[!nan_rows, , drop = FALSE]

# Mfuzz FCM clustering
eset <- new("ExpressionSet", exprs = delta_z)
eset <- standardise(eset)

m_est <- mestimate(eset)
message(sprintf("Estimated fuzzifier m = %.3f", m_est))

set.seed(42)
dmin_vals <- numeric(5)
for (ki in 2:6) {
  cl_tmp <- mfuzz(eset, c = ki, m = m_est)
  dmin_vals[ki - 1] <- min(dist(cl_tmp$centers))
}
write_csv(tibble(k = 2:6, dmin = dmin_vals), file.path(DAT, "dmin_values.csv"))

best_obj <- Inf
best_cl  <- NULL
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
membership <- cl$membership
centroids  <- cl$centers

boot_ari <- numeric(N_BOOT)
set.seed(42)
for (b in seq_len(N_BOOT)) {
  idx <- sample(nrow(eset), round(0.8 * nrow(eset)))
  cl_boot <- mfuzz(eset[idx, ], c = optimal_k, m = m_est)
  orig_assign <- max.col(membership[idx, ])
  boot_assign <- max.col(cl_boot$membership)
  boot_ari[b] <- mclust::adjustedRandIndex(orig_assign, boot_assign)
}
write_csv(tibble(boot = seq_len(N_BOOT), ari = boot_ari),
          file.path(DAT, "07_cluster_stability.csv"))

cluster_assign <- tibble(
  gene = rownames(membership),
  cluster = paste0("C", max.col(membership)),
  membership = apply(membership, 1, max)
)

write_csv(cluster_assign, file.path(DAT, "06_mfuzz_assignments.csv"))

centroid_export <- as.data.frame(centroids) |>
  rownames_to_column("cluster") |>
  mutate(cluster = paste0("C", row_number()))
write_csv(centroid_export, file.path(DAT, "01_panel_A_cluster_profiles.csv"))

core_proteins <- cluster_assign %>% filter(membership >= CORE_THRESH)
core_genes <- core_proteins$gene
core_abund <- abund_mat[core_genes, , drop = FALSE]

group_samples <- list(
  Young_Pre  = sample_meta$sample[sample_meta$age == "Young" & sample_meta$time == "Pre"],
  Young_Post = sample_meta$sample[sample_meta$age == "Young" & sample_meta$time == "Post"],
  Old_Pre    = sample_meta$sample[sample_meta$age == "Old"   & sample_meta$time == "Pre"],
  Old_Post   = sample_meta$sample[sample_meta$age == "Old"   & sample_meta$time == "Post"]
)

group_means <- sapply(group_samples, function(samps) {
  rowMeans(core_abund[, samps, drop = FALSE], na.rm = TRUE)
})

group_z <- t(scale(t(group_means)))
colnames(group_z) <- names(group_samples)

saveRDS(delta_mat, file.path(DAT, "delta_mat.rds"))
saveRDS(group_z, file.path(DAT, "group_z.rds"))

# Enrichment analysis

hallmark_t2g <- msigdbr(species = "Homo sapiens", collection = "H") %>%
  dplyr::select(gs_name, gene_symbol) %>%
  dplyr::rename(term = gs_name, gene = gene_symbol)

universe_genes <- unique(cluster_assign$gene)

enrich_list <- list()

for (cl_id in seq_len(optimal_k)) {
  cl_label <- paste0("C", cl_id)
  cl_genes <- core_proteins$gene[core_proteins$cluster == cl_label]

  res_hall <- enricher(cl_genes, TERM2GENE = hallmark_t2g,
                       universe = universe_genes,
                       pAdjustMethod = "BH",
                       pvalueCutoff = 0.05, qvalueCutoff = 1,
                       minGSSize = 10, maxGSSize = 500)

  res_gobp <- enrichGO(gene = cl_genes,
                        universe = universe_genes,
                        OrgDb = org.Hs.eg.db,
                        keyType = "SYMBOL",
                        ont = "BP",
                        pAdjustMethod = "BH",
                        pvalueCutoff = 0.05, qvalueCutoff = 1,
                        minGSSize = 10, maxGSSize = 500,
                        readable = FALSE)

  combined <- bind_rows(
    if (!is.null(res_hall) && nrow(as.data.frame(res_hall)) > 0)
      as.data.frame(res_hall) %>% mutate(database = "Hallmark") else NULL,
    if (!is.null(res_gobp) && nrow(as.data.frame(res_gobp)) > 0)
      as.data.frame(res_gobp) %>% mutate(database = "GO:BP") else NULL
  )

  if (nrow(combined) > 0) combined <- combined %>% filter(p.adjust < 0.05)
  enrich_list[[cl_label]] <- combined %>% mutate(cluster = cl_label)
}

# rrvgo redundancy reduction (threshold=0.5)

bp_semdata <- tryCatch(
  godata("org.Hs.eg.db", ont = "BP", computeIC = TRUE),
  error = function(e) { message("  Warning: could not compute BP semdata"); NULL }
)

reduce_go_terms <- function(enrich_df, ont, semdata) {
  if (is.null(semdata) || nrow(enrich_df) == 0) return(enrich_df)
  db_label <- paste0("GO:", ont)
  go_terms <- enrich_df %>% filter(database == db_label)
  other_terms <- enrich_df %>% filter(database != db_label)
  if (nrow(go_terms) < 2) return(enrich_df)
  go_ids <- go_terms$ID
  scores <- setNames(go_terms$p.adjust, go_ids)
  reduced <- tryCatch({
    sim_mat <- calculateSimMatrix(go_ids, orgdb = "org.Hs.eg.db",
                                   ont = ont, method = "Rel", semdata = semdata)
    red <- reduceSimMatrix(sim_mat, scores = scores, threshold = 0.5,
                               orgdb = "org.Hs.eg.db")
    go_terms %>% filter(ID %in% unique(red$parent))
  }, error = function(e) {
    message(sprintf("    rrvgo warning (%s): %s -- keeping all", ont, e$message))
    go_terms
  })
  bind_rows(other_terms, reduced)
}

for (cl_label in names(enrich_list)) {
  enrich_list[[cl_label]] <- reduce_go_terms(enrich_list[[cl_label]], "BP", bp_semdata)
}

enrich_top <- bind_rows(enrich_list) %>%
  filter(database %in% c("Hallmark", "GO:BP")) %>%
  group_by(cluster, database) %>%
  slice_min(p.adjust, n = 7, with_ties = FALSE) %>%
  ungroup()

# 1:1 greedy gene-pathway assignment
protein_pathway_links <- bind_rows(lapply(seq_len(optimal_k), function(cl_id) {
  cl_label <- paste0("C", cl_id)
  cl_genes <- core_proteins$gene[core_proteins$cluster == cl_label]
  cl_top   <- enrich_top %>% filter(cluster == cl_label) %>% arrange(p.adjust)
  if (nrow(cl_top) == 0) {
    return(tibble(gene = cl_genes, pathway = "Unmapped",
                  database = NA_character_, cluster = cl_label))
  }
  cl_top$gene_list <- strsplit(cl_top$geneID, "/")
  pw_out <- rep("Unmapped", length(cl_genes))
  db_out <- rep(NA_character_, length(cl_genes))
  for (gi in seq_along(cl_genes)) {
    g <- cl_genes[gi]
    for (j in seq_len(nrow(cl_top))) {
      if (g %in% cl_top$gene_list[[j]]) {
        pw_out[gi] <- cl_top$Description[j]
        db_out[gi] <- cl_top$database[j]
        break
      }
    }
  }
  tibble(gene = cl_genes, pathway = pw_out, database = db_out, cluster = cl_label)
}))

enrich_all_sig <- bind_rows(enrich_list) %>%
  filter(database %in% c("Hallmark", "GO:BP"), p.adjust < 0.05) %>%
  arrange(cluster, p.adjust)
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
  for (i in unmapped_idx) {
    g <- protein_pathway_links$gene[i]
    for (j in seq_len(nrow(cl_all_sig))) {
      if (g %in% cl_all_sig$gene_list[[j]]) {
        protein_pathway_links$pathway[i]  <- cl_all_sig$Description[j]
        protein_pathway_links$database[i] <- cl_all_sig$database[j]
        break
      }
    }
  }
}
enrich_all_sig$gene_list <- NULL

write_csv(enrich_top, file.path(DAT, "03_panel_C_enrichment.csv"))
write_csv(protein_pathway_links, file.path(DAT, "04_panel_C_sankey_links.csv"))

top_hallmark <- enrich_top %>%
  filter(database == "Hallmark") %>%
  group_by(cluster) %>%
  slice_min(p.adjust, n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  mutate(label = clean_pathway_name(Description))

write_csv(top_hallmark, file.path(DAT, "top_hallmark.csv"))
theme_links <- protein_pathway_links %>%
  filter(pathway != "Unmapped") %>%
  mutate(theme = assign_theme(pathway))

write_csv(theme_links, file.path(DAT, "theme_links.csv"))
cluster_theme_counts <- theme_links %>%
  count(cluster, theme, name = "n") %>%
  arrange(cluster, desc(n))
write_csv(cluster_theme_counts, file.path(DAT, "05_panel_D_cluster_themes.csv"))

message("F4 prepare_data.R complete")
