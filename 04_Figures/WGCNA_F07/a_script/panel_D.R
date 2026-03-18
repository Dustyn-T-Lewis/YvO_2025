# Figure 5 — Panel D: GO Enrichment Dotplot
# Generates: panel_D_go_enrichment.pdf/png, c_data/04_panel_D_enrichment_data.csv

setwd(rprojroot::find_rstudio_root_file())
source("04_Figures/shared/style.R")

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
})

RPT <- "04_Figures/WGCNA_F07/b_reports"
DAT <- "04_Figures/WGCNA_F07/c_data"
dir.create(RPT, recursive = TRUE, showWarnings = FALSE)
dir.create(DAT, recursive = TRUE, showWarnings = FALSE)

pdf_device <- get_pdf_device()

go_df <- read_csv(file.path(DAT, "wgcna/wgcna_module_GO_enrichment.csv"), show_col_types = FALSE)
MEs   <- readRDS(file.path(DAT, "MEs.rds"))

message("Panel D: GO enrichment dotplot...")

PD_W <- 250  # panel width mm
PD_H <- 300  # panel height mm

txt_term <- scale_text(BASE_GENE, PD_W)
txt_axis <- scale_text(BASE_STAT, PD_W)

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
  labs(title    = "GO Enrichment",
       subtitle = "Top 5 GO terms per key module (BP preferred; fallback to CC/MF)",
       x = expression(-log[10](p[adj])), y = NULL) +
  FIG_THEME +
  theme(axis.text.y      = element_text(size = txt_term),
        legend.position  = "none")

write_csv(go_plot, file.path(DAT, "04_panel_D_enrichment_data.csv"))

ggsave(file.path(RPT, "panel_D_go_enrichment.pdf"), pD,
       width = PD_W, height = PD_H, units = "mm",
       device = pdf_device, limitsize = FALSE)
ggsave(file.path(RPT, "panel_D_go_enrichment.png"), pD,
       width = PD_W, height = PD_H, units = "mm",
       dpi = 300, limitsize = FALSE)

message("  Panel D saved")
