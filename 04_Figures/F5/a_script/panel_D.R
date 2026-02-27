################################################################################
#   YvO Figure 5 — Panel D: GO Enrichment Dotplot
################################################################################
#
# ── STAT AUDIT (2026-02-27) ──────────────────────────────────────────────────
#
# BH CORRECTION SCOPE
#   GO enrichment is performed PER MODULE by clusterProfiler::enrichGO()
#   in YvO_WGCNA_run.R (line 341: pAdjustMethod = "BH").  Each module's
#   enrichment is an independent ORA with its own BH correction — this is
#   standard practice in WGCNA workflows (Langfelder & Horvath 2008).
#   No cross-module BH correction is applied, which is appropriate because
#   each module's gene list is a biologically distinct hypothesis.
#
# TOP-5 SELECTION
#   slice_head(n = 5) selects the top 5 terms per module by p.adjust.
#   This is post-hoc filtering for display purposes only — all enriched
#   terms are available in the full GO CSV.  The dashed vertical line at
#   -log10(0.05) provides visual reference for significance.
#
# BP PREFERENCE
#   BP ontology is preferred; if no BP terms exist for a module, all
#   ontologies are used.  This is documented in the subtitle.
#
# AUDIT VERDICT: BH per-module is standard for ORA.  Top-5 filtering is
# descriptive, not inferential.  No changes required.
# ─────────────────────────────────────────────────────────────────────────────

source("04_Figures/F5/a_script/YvO_F5_setup.R")

# Four biologically-curated key modules
key_mods <- c("MEblue", "MEgreen", "MEturquoise", "MEyellow")
key_mods <- intersect(key_mods, colnames(MEs))

go_plot <- go_df %>%
  filter(module %in% gsub("^ME", "", key_mods)) %>%
  group_by(module) %>%
  # Prefer BP; if no BP terms, fall back to best available ontology
  mutate(has_bp = any(ONTOLOGY == "BP")) %>%
  filter(if_else(has_bp, ONTOLOGY == "BP", TRUE)) %>%
  arrange(p.adjust) %>%
  slice_head(n = 5) %>%
  ungroup() %>%
  mutate(Description     = clean_pathway_name(Description),
         neg_log10_padj = -log10(p.adjust))

pD <- ggplot(go_plot,
             aes(x = neg_log10_padj,
                 y = reorder_within(Description, neg_log10_padj, module),
                 size = Count, color = module)) +
  geom_point(alpha = 0.8) +
  geom_vline(xintercept = -log10(0.05), linetype = "dashed",
             color = "grey40", linewidth = 0.3) +
  facet_wrap(~ module, scales = "free_y", ncol = 1) +
  scale_y_reordered() +
  scale_color_identity() +
  scale_size_continuous(range = c(1.5, 5), name = "Gene count") +
  labs(title    = "D  GO Enrichment",
       subtitle = "Top 5 GO terms per key module (BP preferred; fallback to CC/MF)",
       x = expression(-log[10](p[adj])), y = NULL) +
  THEME_PUB +
  LEGEND_THEME +
  theme(axis.text.y = element_text(size = 5.5))

write_csv(go_plot, file.path(DAT_DIR, "04_panel_D_enrichment_data.csv"))

ggsave(file.path(RPT_DIR, "panel_D_go_enrichment.pdf"), pD,
       width = 250, height = 300, units = "mm",
       device = pdf, limitsize = FALSE)
ggsave(file.path(RPT_DIR, "panel_D_go_enrichment.png"), pD,
       width = 250, height = 300, units = "mm",
       dpi = 300, limitsize = FALSE)

cat("Panel D saved\n")
