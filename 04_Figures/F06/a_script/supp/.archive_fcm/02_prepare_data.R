# Supplementary: FCM — prepare_data.R: Mfuzz FCM clustering + ORA enrichment cache
# Run BEFORE panel scripts. 50-start FCM (k=4), bootstrap ARI,
# per-cluster ORA (Hallmark + GO:BP + rrvgo), theme assignment.

setwd(rprojroot::find_rstudio_root_file())
source("04_Figures/shared/style.R")

library(readr)
library(dplyr)
library(tidyr)
library(stringr)
library(tibble)
library(e1071)
source("04_Figures/shared/pathway_utils.R")

# FCM via e1071::cmeans directly (no Mfuzz/Biobase dependency needed)
# mestimate: Schwämmle & Jensen 2010 (DOI:10.1093/bioinformatics/btq534)
# Full formula from Mfuzz source — NOT the simplified approximation
mestimate <- function(mat) {
  N <- nrow(mat); D <- ncol(mat)
  1 + (1418/N + 22.05) * D^(-2) +
      (12.33/N + 0.243) * D^(-0.0406 * log(N) - 0.1134)
}

set.seed(42)

DAT <- "04_Figures/F06/c_data/supp/fcm"
RPT <- "04_Figures/F06/b_reports/supp/fcm"
dir.create(DAT, recursive = TRUE, showWarnings = FALSE)
dir.create(RPT, recursive = TRUE, showWarnings = FALSE)

# k=4 selected via Dmin elbow test (see dmin_elbow.R); Dmin values for k=2:6
# computed below. k=4 maximizes inter-cluster separation before plateau.
optimal_k   <- 4
CORE_THRESH <- 0.5
N_STARTS    <- 50
N_BOOT      <- 100

imputed <- read_csv("02_Imputation/c_data/01_imputed.csv", show_col_types = FALSE)
annot_cols  <- c("uniprot_id", "protein", "gene", "description")
sample_cols <- setdiff(colnames(imputed), annot_cols)

# Metadata from DAList (not regex) — canonical source of Group/Timepoint
dal_meta <- as.data.frame(
  readRDS("02_Imputation/c_data/01_DAList_imputed.rds")$metadata)
sample_meta <- tibble(
  sample  = dal_meta$Col_ID,
  age     = dal_meta$Group,
  time    = dal_meta$Timepoint,
  subject = sub("_(Pre|Post)$", "", dal_meta$Col_ID)
) |>
  filter(sample %in% sample_cols)
write_csv(sample_meta, file.path(DAT, "sample_meta.csv"))

pre_subjects  <- sample_meta$subject[sample_meta$time == "Pre"]
post_subjects <- sample_meta$subject[sample_meta$time == "Post"]
paired_subjects <- intersect(pre_subjects, post_subjects)

young_subjects <- sort(paired_subjects[paired_subjects %in%
  sample_meta$subject[sample_meta$age == "Young"]])
old_subjects   <- sort(paired_subjects[paired_subjects %in%
  sample_meta$subject[sample_meta$age == "Old"]])
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

# FCM clustering (e1071::cmeans)
# delta_z is already z-scored per protein; ensure clean numeric matrix for cmeans
fcm_mat <- matrix(as.numeric(delta_z), nrow = nrow(delta_z), ncol = ncol(delta_z),
                  dimnames = dimnames(delta_z))

m_est <- mestimate(fcm_mat)
message(sprintf("Estimated fuzzifier m = %.3f", m_est))

set.seed(42)
dmin_vals <- numeric(5)
for (ki in 2:6) {
  cl_tmp <- cmeans(fcm_mat, centers = ki, m = m_est)
  dmin_vals[ki - 1] <- min(dist(cl_tmp$centers))
}
write_csv(tibble(k = 2:6, dmin = dmin_vals), file.path(DAT, "dmin_values.csv"))

best_obj <- Inf
best_cl  <- NULL
for (s in seq_len(N_STARTS)) {
  set.seed(41 + s)
  cl_try <- cmeans(fcm_mat, centers = optimal_k, m = m_est)
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
  idx <- sample(nrow(fcm_mat), round(0.8 * nrow(fcm_mat)))
  cl_boot <- cmeans(fcm_mat[idx, , drop = FALSE], centers = optimal_k, m = m_est)
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

group_z <- t(scale(t(group_means)))
colnames(group_z) <- names(group_samples)

saveRDS(delta_mat, file.path(DAT, "delta_mat.rds"))
saveRDS(group_z, file.path(DAT, "group_z.rds"))

# Multi-database ORA enrichment (replaces Hallmark + GO:BP + rrvgo)
# Uses run_ora_deduplicated() from pathway_utils.R for consistency with F03-F04

universe_genes <- unique(cluster_assign$gene)
pw_collection <- build_pathway_collection(include_goslim = TRUE)

enrich_list <- list()

for (cl_id in seq_len(optimal_k)) {
  cl_label <- paste0("C", cl_id)
  cl_genes <- core_proteins$gene[core_proteins$cluster == cl_label]

  ora_res <- tryCatch(
    run_ora_deduplicated(
      genes          = cl_genes,
      universe       = universe_genes,
      pathways       = pw_collection,
      jaccard_cutoff = 0.5,
      min_size       = 10,
      max_size       = 500,
      padj_cutoff    = 0.05
    ),
    error = function(e) { warning(sprintf("ORA failed for '%s': %s", cl_label, e$message)); NULL }
  )

  if (!is.null(ora_res) && nrow(ora_res) > 0) {
    ora_res$cluster <- cl_label
    # Backward-compat columns for downstream panels
    ora_res$Description <- clean_pathway_name(ora_res$pathway)
    ora_res$geneID <- vapply(ora_res$overlapGenes, function(g) {
      paste(g, collapse = "/")
    }, character(1))
    ora_res$Count    <- ora_res$overlap
    ora_res$p.adjust <- ora_res$padj
    ora_res$ID       <- ora_res$pathway
    ora_res$overlapGenes <- NULL
    enrich_list[[cl_label]] <- ora_res
  } else {
    enrich_list[[cl_label]] <- tibble()
  }
  message(sprintf("  Cluster %s: %d enriched terms (padj < 0.05)", cl_label,
                  if (!is.null(ora_res)) nrow(ora_res) else 0L))
}

enrich_top <- bind_rows(enrich_list) %>%
  group_by(cluster, database) %>%
  slice_min(padj, n = 7, with_ties = FALSE) %>%
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
  filter(padj < 0.05) %>%
  arrange(cluster, padj)
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

message("FCM prepare_data.R complete")

})
