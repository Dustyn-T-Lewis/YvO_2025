################################################################################
#   Figure 6 — Panel F: Convergence UpSet
#   Multi-method convergence: WGCNA hubs, Aging DEPs, Interaction DEPs,
#   FCM core members, top GS proteins
################################################################################

source("04_Figures/F6/a_script/YvO_F6_setup.R")

# ---- Verify cross-figure dependencies ----
fcm_path <- "04_Figures/F4/c_data/fcm_assignments_aging_k4.csv"
if (!file.exists(fcm_path)) {
  stop("F6 Panel F requires F4 output. Please run F4 first.\n",
       "  Missing: ", fcm_path, call. = FALSE)
}

# ---- Load additional data sources ----
dep_results <- read_csv("03_DEP/c_data/combined_results.csv",
                         show_col_types = FALSE)

fcm_aging <- read_csv(fcm_path, show_col_types = FALSE)

# ---- 1. WGCNA hubs (kME > 0.7 from hub file; fallback to all hubs) ----
all_hub_genes <- hub_df %>% filter(kME > 0.7) %>% pull(gene) %>% unique()
if (length(all_hub_genes) < 20) {
  all_hub_genes <- hub_df %>% pull(gene) %>% unique()
}
all_hub_genes <- all_hub_genes[!is.na(all_hub_genes) & all_hub_genes != ""]

# ---- 2. Aging DEPs (FDR < 0.05) ----
aging_dep_genes <- dep_results %>%
  filter(adj.P.Val_Aging < 0.05) %>%
  pull(gene) %>% unique()
aging_dep_genes <- aging_dep_genes[!is.na(aging_dep_genes) &
                                     aging_dep_genes != ""]

# ---- 3. Interaction DEPs (relaxed: nominal p < 0.01) ----
interaction_dep_genes <- dep_results %>%
  filter(P.Value_Interaction < 0.01) %>%
  pull(gene) %>% unique()
interaction_dep_genes <- interaction_dep_genes[
  !is.na(interaction_dep_genes) & interaction_dep_genes != ""]

# ---- 4. FCM core members (any membership > 0.5) ----
mem_cols <- grep("^mem_", names(fcm_aging), value = TRUE)
if (length(mem_cols) > 0) {
  fcm_core_genes <- fcm_aging %>%
    filter(if_any(all_of(mem_cols), ~ . > 0.5)) %>%
    pull(gene) %>% unique()
} else {
  fcm_core_genes <- fcm_aging %>% pull(gene) %>% unique()
}
fcm_core_genes <- fcm_core_genes[!is.na(fcm_core_genes) &
                                   fcm_core_genes != ""]

# ---- 5. Top GS proteins for delta VL (top 100 by |GS|) ----
gs_sorted <- sort(abs(gs_vl[, 1]), decreasing = TRUE)
gs_top_ids <- names(gs_sorted)[seq_len(min(100, length(gs_sorted)))]
gs_top_genes <- module_df$gene[match(gs_top_ids, module_df$uniprot_id)]
gs_top_genes <- gs_top_genes[!is.na(gs_top_genes) & gs_top_genes != ""]

# ---- Assemble and filter empty sets ----
upset_input <- list(
  "WGCNA Hubs"        = all_hub_genes,
  "Aging DEPs"        = aging_dep_genes,
  "Interaction DEPs"  = interaction_dep_genes,
  "FCM Core"          = fcm_core_genes,
  "Top GS (delta VL)" = gs_top_genes
)
upset_input <- upset_input[lengths(upset_input) > 0]

# ---- Count consensus genes (>= 3 methods) ----
all_genes_union   <- unique(unlist(upset_input))
gene_method_count <- sapply(all_genes_union, function(g) {
  sum(vapply(upset_input, function(s) g %in% s, logical(1)))
})
consensus_genes_3 <- names(gene_method_count[gene_method_count >= 3])
consensus_genes_2 <- names(gene_method_count[gene_method_count >= 2])
cat(sprintf("Consensus genes: %d (>= 3 methods), %d (>= 2 methods)\n",
            length(consensus_genes_3), length(consensus_genes_2)))

# ---- Save standalone UpSet PDF ----
if (length(upset_input) >= 2) {
  pdf(file.path(RPT_DIR, "panel_F_upset.pdf"), width = 14, height = 10)
  upset(fromList(upset_input),
        sets        = rev(names(upset_input)),
        keep.order  = TRUE,
        order.by    = "freq",
        nsets       = length(upset_input),
        nintersects = 15,
        text.scale  = c(1.8, 1.3, 1.3, 1.3, 1.8, 1.3),
        main.bar.color = "grey30",
        sets.bar.color = "grey50",
        mb.ratio    = c(0.65, 0.35))
  dev.off()
}

# ---- Render to PNG ----
png(file.path(RPT_DIR, "panel_F_upset.png"),
    width = 3600, height = 2000, res = 300, bg = "white")
tryCatch({
  if (length(upset_input) >= 2) {
    upset(fromList(upset_input),
          sets        = rev(names(upset_input)),
          keep.order  = TRUE,
          order.by    = "freq",
          nsets       = length(upset_input),
          nintersects = 15,
          text.scale  = c(1.8, 1.3, 1.3, 1.3, 1.8, 1.3),
          main.bar.color = "grey30",
          sets.bar.color = "grey50",
          mb.ratio    = c(0.65, 0.35))
  } else {
    plot.new()
    text(0.5, 0.5, "Insufficient sets for UpSet plot", cex = 1.2)
  }
}, error = function(e) {
  plot.new()
  text(0.5, 0.5, paste("UpSet error:", conditionMessage(e)), cex = 0.8)
})
dev.off()

# ---- Export consensus genes CSV ----
if (length(consensus_genes_2) > 0) {
  write_csv(
    tibble(gene = consensus_genes_2,
           n_methods = gene_method_count[consensus_genes_2],
           tier = ifelse(gene_method_count[consensus_genes_2] >= 3,
                         "high", "moderate")),
    file.path(DAT_DIR, "fig6_consensus_genes.csv"))
}

cat("Panel F done\n")
