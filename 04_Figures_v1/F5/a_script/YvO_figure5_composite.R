################################################################################
#   Figure 5 — Composite Assembly (Main + Supplementary)
#
#   Main Figure 5 (~180 x 240 mm):
#     A — Module dendrogram with color bar
#     B — Module-trait correlation heatmap
#     C — Simplified eigengene slopes + trait strip (6 key modules)
#     D — Module preservation (compacted bar chart)
#
#   Supplementary Figure S5:
#     Full triptych (Panel C)
#     GO enrichment dotplot (Panel D)
#     Hub protein networks (Panel E)
#     Age gap closure
#     FCM-WGCNA overlap
################################################################################

# === 0. Ensure WGCNA outputs exist ==========================================

wgcna_cache <- "04_Figures/F5/c_data/wgcna/wgcna_network.rds"
if (!file.exists(wgcna_cache)) {
  message("WGCNA cache not found — running WGCNA pipeline...")
  source("04_Figures/F5/a_script/YvO_WGCNA_run.R")
}

# === 1. Source all panel scripts (F5 setup loaded automatically by each) =====

source("04_Figures/F5/a_script/panel_A.R")     # → pA
source("04_Figures/F5/a_script/panel_B.R")     # → pB
source("04_Figures/F5/a_script/panel_C.R")     # → panel_C, eigen_all, all_enrich
source("04_Figures/F5/a_script/panel_D.R")     # → pD, go_plot
source("04_Figures/F5/a_script/panel_E.R")     # → panel_E, all_node_data, all_periph_data
source("04_Figures/F5/a_script/panel_F.R")     # → pF, pres_df

message("All F5 panel scripts sourced — assembling composite...")

# === 2. Build simplified Panel C for main figure ============================
# Eigengene slope plots + trait correlation strip only, for 6 key modules

MAIN_MODULES <- c("blue", "brown", "turquoise", "green", "black", "pink")
GROUP_COLS <- c("Young_Pre", "Young_Post", "Old_Pre", "Old_Post")

# Trait correlation data from Panel B
trait_heatmap_df <- read_csv(file.path(DAT_DIR, "02_panel_B_heatmap_data.csv"),
                             show_col_types = FALSE)
KEY_TRAITS      <- c("age_num", "time_num", "VL_thick_cm", "Type_II_fCSA")
KEY_TRAIT_LABELS <- c("Age", "Training", "VL Thick.", "Type II fCSA")

build_simple_row <- function(mod, show_xlab) {
  me_col <- paste0("ME", mod)
  border_col <- darken_color(mod, 0.3)
  n_mod <- sum(module_df$module_color == mod)
  bio_lbl <- mod_display_label(mod)

  # Accent bar
  p_accent <- ggplot() +
    geom_rect(aes(xmin = 0, xmax = 1, ymin = 0, ymax = 1), fill = mod) +
    annotate("text", x = 0.5, y = 0.5,
             label = paste0(bio_lbl, "\nn=", n_mod),
             size = 2.8, fontface = "bold",
             color = if (mod %in% c("yellow", "pink")) "grey20" else "white",
             lineheight = 0.85, angle = 90) +
    coord_cartesian(expand = FALSE) +
    theme_void() +
    theme(plot.margin = margin(0, 0, 0, 0))

  # Trait tile strip
  trait_mod <- trait_heatmap_df %>%
    filter(module == paste0("ME", mod), trait %in% KEY_TRAITS) %>%
    mutate(trait_label = factor(trait_label, levels = rev(KEY_TRAIT_LABELS)))

  p_trait <- ggplot(trait_mod, aes(x = 1, y = trait_label, fill = cor)) +
    geom_tile(color = "grey60", linewidth = 0.3) +
    geom_text(aes(label = stars), size = 2.5, fontface = "bold", vjust = 0.75) +
    scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B",
                         midpoint = 0, limits = c(-1, 1), guide = "none") +
    scale_x_continuous(expand = expansion(0)) +
    labs(x = NULL, y = NULL) +
    THEME_PUB +
    theme(
      axis.text.y  = element_text(size = 6, face = "bold"),
      axis.text.x  = element_blank(),
      axis.ticks   = element_blank(),
      panel.border = element_rect(colour = "grey70", linewidth = 0.3, fill = NA),
      plot.margin  = margin(t = 1, r = 1, b = if (show_xlab) 4 else 0, l = 1)
    )

  # Eigengene slope plot
  eig_mod <- eigen_all %>% filter(module == mod) %>%
    mutate(age_fct  = factor(age, levels = c("Young", "Old")),
           time_fct = factor(time, levels = c("Pre", "Post")))

  eig_means <- eig_mod %>%
    group_by(age_fct, time_fct) %>%
    summarise(eigengene = mean(eigengene), .groups = "drop")

  p_eigen <- ggplot(eig_mod, aes(x = time_fct, y = eigengene)) +
    geom_line(aes(group = subject, color = age_fct),
              alpha = 0.25, linewidth = 0.3) +
    geom_point(aes(fill = age_fct), shape = 21, size = 1.0,
               alpha = 0.45, stroke = 0.2, color = "white") +
    geom_line(data = eig_means, aes(group = 1, color = age_fct),
              linewidth = 1.3) +
    geom_point(data = eig_means, aes(fill = age_fct), shape = 21,
               size = 2.5, stroke = 0.6, color = "white") +
    facet_wrap(~age_fct, nrow = 1) +
    scale_color_manual(values = AGE_COLORS, guide = "none") +
    scale_fill_manual(values = AGE_COLORS, guide = "none") +
    labs(x = NULL, y = if (show_xlab) "Eigengene" else NULL) +
    THEME_PUB +
    theme(
      strip.text       = element_text(size = 7, face = "bold", color = "grey15"),
      strip.background = element_rect(fill = "grey93", colour = NA),
      axis.text.x  = if (show_xlab) element_text(size = 7, face = "bold")
                      else element_blank(),
      axis.ticks.x = if (show_xlab) element_line() else element_blank(),
      axis.text.y  = element_text(size = 6),
      axis.title.y = if (show_xlab) element_text(size = 7, face = "bold") else element_blank(),
      panel.border = element_rect(colour = border_col, linewidth = 0.5, fill = NA),
      panel.spacing = unit(2, "pt"),
      plot.margin  = margin(t = 1, r = 1, b = if (show_xlab) 4 else 0, l = 1)
    )

  (p_accent | p_trait | p_eigen) +
    plot_layout(widths = c(0.08, 0.15, 0.77), nrow = 1)
}

simple_rows <- lapply(seq_along(MAIN_MODULES), function(i) {
  build_simple_row(MAIN_MODULES[i], show_xlab = (i == length(MAIN_MODULES)))
})

panel_C_simple <- wrap_plots(simple_rows, ncol = 1) +
  plot_annotation(
    title    = "C   Module Eigengene Dynamics",
    subtitle = "Trait correlation  |  Eigengene Pre\u2192Post (Young vs Old)",
    theme = theme(
      plot.title    = element_text(face = "bold", size = 9, color = "black"),
      plot.subtitle = element_text(size = 6.5, color = "grey25", face = "italic")
    )
  )

# === 3. Assemble MAIN Figure 5 ==============================================

# Row 1: A (dendrogram, full width)
# Row 2: B (heatmap) | simplified C (eigengene)
# Row 3: F (preservation bar chart, compact)

row1 <- pA

row2 <- (pB | panel_C_simple) +
  plot_layout(widths = c(0.55, 0.45))

row3 <- pF

fig5_main <- row1 / row2 / row3 +
  plot_layout(heights = c(0.22, 0.60, 0.18)) +
  plot_annotation(
    title = "WGCNA Co-Expression Network Architecture",
    subtitle = sprintf("Signed hybrid network | %d modules from %s proteins | %d samples",
                       n_modules, format(ncol(datExpr), big.mark = ","), nrow(datExpr)),
    theme = theme(
      plot.title    = element_text(face = "bold", size = 11, hjust = 0.5),
      plot.subtitle = element_text(size = 8, color = "grey30", hjust = 0.5)
    )
  )

ggsave(file.path(RPT_DIR, "Figure_5.pdf"), fig5_main,
       width = 180, height = 240, units = "mm",
       device = pdf, limitsize = FALSE)
ggsave(file.path(RPT_DIR, "Figure_5.png"), fig5_main,
       width = 180, height = 240, units = "mm",
       dpi = 300, limitsize = FALSE)

cat("Saved Figure_5.pdf (main, 180x240 mm)\n")

# === 4. Assemble SUPPLEMENTARY Figure S5 ====================================

# Row 1: Full triptych (Panel C)
# Row 2: GO dotplot (Panel D) | Hub networks (Panel E)
# Row 3: Age gap closure (from panel_C.R's p_closure)

supp_row1 <- wrap_elements(full = panel_C)

supp_row2 <- (pD | wrap_elements(full = panel_E)) +
  plot_layout(widths = c(0.30, 0.70))

fig5_supp <- supp_row1 / supp_row2 +
  plot_layout(heights = c(0.50, 0.50)) +
  plot_annotation(
    title = "Supplementary Figure S5: WGCNA Module Characterization",
    subtitle = sprintf("Full functional triptych, GO enrichment, and hub protein networks | %d modules",
                       n_modules),
    theme = theme(
      plot.title    = element_text(face = "bold", size = 11, hjust = 0.5),
      plot.subtitle = element_text(size = 8, color = "grey30", hjust = 0.5)
    )
  )

ggsave(file.path(RPT_DIR, "Figure_S5.pdf"), fig5_supp,
       width = 580, height = 650, units = "mm",
       device = pdf, limitsize = FALSE)
ggsave(file.path(RPT_DIR, "Figure_S5.png"), fig5_supp,
       width = 580, height = 650, units = "mm",
       dpi = 300, limitsize = FALSE)

cat("Saved Figure_S5.pdf (supplementary)\n")

# === 5. Supplementary Excel Workbook ========================================

library(openxlsx)

wb <- createWorkbook()

addWorksheet(wb, "A_dendrogram")
addWorksheet(wb, "B_heatmap")
addWorksheet(wb, "C_eigengenes")
addWorksheet(wb, "C_enrichment")
addWorksheet(wb, "D_go_enrichment")
addWorksheet(wb, "E_hub_network")
addWorksheet(wb, "F_preservation")
addWorksheet(wb, "_metadata")

# Panel A: dendrogram module assignments
writeData(wb, "A_dendrogram", dendro_data)

# Panel B: module-trait heatmap
writeData(wb, "B_heatmap", heat_df)

# Panel C: module eigengene data
writeData(wb, "C_eigengenes", eigen_all)

# Panel C: triptych GO enrichment
writeData(wb, "C_enrichment", all_enrich)

# Panel D: GO enrichment dotplot
writeData(wb, "D_go_enrichment", go_plot)

# Panel E: hub network nodes (hub + context ring combined)
hub_combined <- bind_rows(
  all_node_data  %>% mutate(node_type = "hub"),
  all_periph_data %>% mutate(node_type = "periphery")
)
writeData(wb, "E_hub_network", hub_combined)

# Panel F: module preservation
writeData(wb, "F_preservation", pres_df)

# Metadata sheet
meta_sheet <- tibble(
  sheet       = c("A_dendrogram", "B_heatmap", "C_eigengenes", "C_enrichment",
                   "D_go_enrichment", "E_hub_network", "F_preservation"),
  description = c(
    "Module assignments from WGCNA dendrogram (unmerged + merged colors)",
    "Module-trait correlation heatmap (Pearson r, BH-corrected p-values)",
    "Module eigengene values per sample for 6 key modules",
    "GO enrichment results for triptych panel (top terms per module)",
    "GO enrichment dotplot data (top 5 per key module, BP preferred)",
    "Hub protein network nodes (top 30 by kME) + periphery ring proteins",
    "Module preservation Zsummary (Pre -> Post training, 200 permutations)"
  ),
  source_panel = c("A", "B", "C", "C", "D", "E", "F"),
  n_rows       = c(nrow(dendro_data), nrow(heat_df), nrow(eigen_all),
                    nrow(all_enrich), nrow(go_plot), nrow(hub_combined),
                    nrow(pres_df))
)
writeData(wb, "_metadata", meta_sheet)

saveWorkbook(wb, file.path(DAT_DIR, "Figure_5_data.xlsx"), overwrite = TRUE)
cat("Saved Figure_5_data.xlsx\n")
