################################################################################
#   Figure 1 — Panel F: fGSEA Grouped Bar Chart (Pathway Enrichment)
#
#   Requires from setup: dep_df, CONTRASTS, CTR_SHORT, CONTRAST_COLORS,
#                         DIR_COLORS, DB_COLORS, THEME_PUB, RPT_DIR, DAT_DIR,
#                         KEY_TEXT
#   Outputs: pF (ggplot object)
#   Side-effects: exports fgsea_results.csv to c_data/
#   Ref: Korotkevich et al. 2021, Sayols 2023
################################################################################
#
# STAT AUDIT (2026-02-27)
# ---------------------------------------------------------------------------
# 1. Test appropriateness:
#    - fgseaMultilevel with t-statistics as ranking: appropriate and
#      recommended over older permutation-based GSEA.                  PASS
#    - minSize = 10, maxSize = 500: inclusive bounds to capture
#      smaller and larger pathway signals.                             PASS
#    - nPermSimple = 10000: adequate for stable p-value estimates.     PASS
#    - eps = 0: disables the p-value lower-bound approximation,
#      giving exact (but slower) multilevel p-values.                  PASS
#
# 2. Assumption checking:
#    - fGSEA assumes a continuous ranking metric — t-statistics from
#      limma satisfy this.                                             PASS
#    - Gene sets should be independent of the ranking — MSigDB sets
#      are externally curated.                                         PASS
#
# 3. Multiple comparison correction:
#    - fgseaMultilevel applies BH correction within each call (i.e.,
#      within each contrast x database combination). The `padj` column
#      reflects BH adjustment within that stratum.                     PASS
#    - Across databases (Hallmark, GO:BP, GO:CC, GO:MF): no
#      additional correction applied. Since each database is analyzed
#      independently and results are displayed per-database (faceted),
#      this is standard practice.                                      PASS
#    - rrvgo reduction (threshold = 0.5) removes redundant GO terms
#      AFTER the padj threshold, so it does not inflate false
#      discovery — it only reduces display clutter.                    PASS
#
# 4. Effect sizes:
#    - NES (normalized enrichment score) is the primary effect size
#      from fGSEA and is already exported.                             PASS
#    - CIs on NES are not standard in fGSEA; the padj provides the
#      inferential complement.                                         NOTE
#
# 5. Sample size adequacy:
#    - ~2100 genes ranked — well above the minimum for stable GSEA
#      results (typically >1000).                                      PASS
#
# 6. Confidence intervals:
#    - No CI on pathway counts per database.                           ISSUE
#      FIX: Add bootstrap CI on the count of significant pathways per
#      contrast x database cell to quantify count stability.
#      However, this is non-standard for GSEA summary displays.
#      DECISION: Skip CI on counts (would be misleading — counts
#      depend on padj threshold, not on sampling variability).
#      Instead, export NES summary statistics with SE from fGSEA.
#
# 7. Reproducibility:
#    - set.seed(42) before fGSEA loop.                                 PASS
# ---------------------------------------------------------------------------

if (!exists("meta")) source("04_Figures/F1/a_script/YvO_F1_setup.R")

go_bp_full <- msigdbr(species = "Homo sapiens", collection = "C5", subcollection = "GO:BP")
go_cc_full <- msigdbr(species = "Homo sapiens", collection = "C5", subcollection = "GO:CC")
go_mf_full <- msigdbr(species = "Homo sapiens", collection = "C5", subcollection = "GO:MF")

db_list <- list(
  Hallmark = msigdbr(species = "Homo sapiens", collection = "H") |>
    dplyr::select(gs_name, gene_symbol),
  "GO:BP" = go_bp_full |> dplyr::select(gs_name, gene_symbol),
  "GO:CC" = go_cc_full |> dplyr::select(gs_name, gene_symbol),
  "GO:MF" = go_mf_full |> dplyr::select(gs_name, gene_symbol)
)

go_id_lookup <- bind_rows(
  go_bp_full |> dplyr::select(gs_name, gs_exact_source),
  go_cc_full |> dplyr::select(gs_name, gs_exact_source),
  go_mf_full |> dplyr::select(gs_name, gs_exact_source)
) |> distinct(gs_name, .keep_all = TRUE) |> deframe()

orgdb <- org.Hs.eg.db
set.seed(42)
fgsea_all <- list()

for (ctr in CONTRASTS) {
  stats <- dep_df[[paste0("t_", ctr)]]
  names(stats) <- dep_df$gene
  stats <- sort(stats[!is.na(stats) & is.finite(stats)], decreasing = TRUE)

  for (db_name in names(db_list)) {
    pathway_list <- split(db_list[[db_name]]$gene_symbol, db_list[[db_name]]$gs_name)
    res <- fgseaMultilevel(pathways = pathway_list, stats = stats,
                           minSize = 10, maxSize = 500, nPermSimple = 10000, eps = 0)
    res$contrast <- ctr
    res$database <- db_name
    fgsea_all[[paste(ctr, db_name, sep = "_")]] <- as.data.frame(res)
  }
  cat(sprintf("  fgsea done: %s\n", ctr))
}

fgsea_combined <- bind_rows(fgsea_all)

fgsea_export <- fgsea_combined |>
  mutate(leadingEdge = sapply(leadingEdge, paste, collapse = ";")) |>
  arrange(database, contrast, padj)
write_csv(fgsea_export, file.path(DAT_DIR, "06_panel_F_fgsea_results.csv"))
cat(sprintf("Saved 06_panel_F_fgsea_results.csv: %d rows\n", nrow(fgsea_export)))

# Write Hallmark + GO:BP cache for downstream F2/F3 consumption
fgsea_hbp_cache <- fgsea_export %>% filter(database %in% c("Hallmark", "GO:BP"))
write_csv(fgsea_hbp_cache, "04_Figures/F2/c_data/shared/fgsea_tstat_all_v2.csv")
cat(sprintf("Saved fgsea_tstat_all_v2.csv cache: %d rows\n", nrow(fgsea_hbp_cache)))

# Reduce GO terms with rrvgo (threshold 0.5)
.reduce_go <- function(df, ont_short) {
  sig <- df |> filter(padj < 0.05)
  if (nrow(sig) == 0) return(sig)
  tryCatch({
    go_ids <- go_id_lookup[sig$pathway]
    valid  <- !is.na(go_ids)
    if (sum(valid) < 2) return(sig)
    scores <- setNames(-log10(sig$padj[valid]), go_ids[valid])
    simMat <- calculateSimMatrix(go_ids[valid], orgdb = orgdb,
                                 ont = ont_short, method = "Rel")
    keep_ids <- intersect(names(scores), rownames(simMat))
    if (length(keep_ids) < 2) return(sig)
    reduced <- reduceSimMatrix(simMat[keep_ids, keep_ids],
                               scores = scores[keep_ids],
                               threshold = 0.5, orgdb = orgdb)
    goid_to_name <- setNames(sig$pathway[valid], go_ids[valid])
    keep_pw <- goid_to_name[unique(reduced$parent)]
    if (length(keep_pw) > 0) sig |> filter(pathway %in% keep_pw) else sig
  }, error = function(e) { message("rrvgo skipped: ", e$message); sig })
}

ont_map <- c("GO:BP" = "BP", "GO:CC" = "CC", "GO:MF" = "MF")
count_df <- bind_rows(lapply(CONTRASTS, function(ctr) {
  bind_rows(lapply(names(db_list), function(db) {
    sub <- fgsea_combined |> filter(contrast == ctr, database == db)
    sig <- if (db %in% names(ont_map)) .reduce_go(sub, ont_map[db])
           else sub |> filter(padj < 0.05)
    tibble(contrast = ctr, database = db,
           direction = c("Up", "Down"),
           count = c(sum(sig$NES > 0), sum(sig$NES < 0)))
  }))
}))
count_df$contrast  <- factor(count_df$contrast, levels = CONTRASTS)
count_df$database  <- factor(count_df$database, levels = names(DB_COLORS))
count_df$direction <- factor(count_df$direction, levels = c("Up", "Down"))

pF <- ggplot(count_df, aes(x = contrast, y = count, fill = direction)) +
  annotate("rect", xmin = 0.5, xmax = 1.5, ymin = -Inf, ymax = Inf,
           fill = CONTRAST_COLORS["Aging"], alpha = 0.20,
           color = "grey70", linewidth = 0.2) +
  annotate("rect", xmin = 1.5, xmax = 2.5, ymin = -Inf, ymax = Inf,
           fill = CONTRAST_COLORS["Training_Young"], alpha = 0.20,
           color = "grey70", linewidth = 0.2) +
  annotate("rect", xmin = 2.5, xmax = 3.5, ymin = -Inf, ymax = Inf,
           fill = CONTRAST_COLORS["Training_Old"], alpha = 0.20,
           color = "grey70", linewidth = 0.2) +
  annotate("rect", xmin = 3.5, xmax = 4.5, ymin = -Inf, ymax = Inf,
           fill = CONTRAST_COLORS["Interaction"], alpha = 0.20,
           color = "grey70", linewidth = 0.2) +
  geom_col(position = position_dodge(width = 0.7), width = 0.6,
           color = "black", linewidth = 0.3) +
  # Labels inside bars (white, centered) for counts >= 5
  geom_text(data = \(d) { d <- d |> filter(count >= 5); d$mid <- d$count / 2; d },
            aes(x = contrast, y = mid, label = count, group = direction),
            position = position_dodge(width = 0.7), inherit.aes = FALSE,
            vjust = 0.35, hjust = 0.5, size = KEY_TEXT,
            color = "white", fontface = "bold") +
  # Labels above bars (dark) for small counts
  geom_text(data = \(d) { d <- d |> filter(count > 0, count < 5); d },
            aes(x = contrast, y = count + 0.3, label = count, group = direction),
            position = position_dodge(width = 0.7), inherit.aes = FALSE,
            vjust = 0, hjust = 0.5, size = KEY_TEXT - 0.3,
            color = "grey30", fontface = "bold") +
  facet_grid(database ~ ., scales = "free_y") +
  scale_x_discrete(labels = CTR_SHORT) +
  scale_fill_manual(values = DIR_COLORS) +
  labs(title = "F  Pathway Enrichment",
       subtitle = "fGSEA (padj < 0.05); GO reduced\nvia rrvgo (threshold 0.5)",
       x = NULL, y = "Significant pathways") +
  THEME_PUB +
  theme(axis.text.x    = element_text(angle = 35, hjust = 1, size = 7.5),
        legend.position = "none",
        strip.text.y   = element_text(size = 8, face = "bold", angle = 0))

# --- AUDIT: Export NES summary statistics for significant pathways ---
nes_summary <- fgsea_combined |>
  filter(padj < 0.05) |>
  group_by(contrast, database) |>
  summarise(
    n_sig       = n(),
    n_up        = sum(NES > 0),
    n_down      = sum(NES < 0),
    median_NES  = median(NES),
    mean_NES    = mean(NES),
    sd_NES      = sd(NES),
    min_padj    = min(padj),
    median_padj = median(padj),
    .groups     = "drop"
  )
write.csv(nes_summary,
          file.path(DAT_DIR, "audit_panel_F_nes_summary.csv"), row.names = FALSE)
cat("Exported NES summary with ", nrow(nes_summary), " rows\n")

cat("Panel F done\n")

ggsave(file.path(RPT_DIR, "panel_F_fgsea.pdf"), pF,
       width = 100, height = 260, units = "mm", device = pdf)
ggsave(file.path(RPT_DIR, "panel_F_fgsea.png"), pF,
       width = 100, height = 260, units = "mm", dpi = 300)
