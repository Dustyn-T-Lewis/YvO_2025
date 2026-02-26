################################################################################
#   YvO Figure 5 — Panel D: GO Enrichment Dotplot
################################################################################

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

write_csv(go_plot, file.path(DAT_DIR, "fig5_panel_D_enrichment_data.csv"))

ggsave(file.path(RPT_DIR, "panel_D_go_enrichment.pdf"), pD,
       width = 250, height = 300, units = "mm",
       device = pdf, limitsize = FALSE)
ggsave(file.path(RPT_DIR, "panel_D_go_enrichment.png"), pD,
       width = 250, height = 300, units = "mm",
       dpi = 300, limitsize = FALSE)

cat("Panel D saved\n")
