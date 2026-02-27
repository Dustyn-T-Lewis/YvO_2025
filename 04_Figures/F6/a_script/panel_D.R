################################################################################
#   Figure 6 — Panel D: kME x GS Scatter (Hub Proteins)
#   Module membership vs |gene significance| for the best predictive module
################################################################################

source("04_Figures/F6/a_script/YvO_F6_setup.R")

# ---- Best predictive module (highest |r| to either phenotype delta) ----
best_mod_name  <- pred_cor %>% arrange(desc(max_r)) %>%
  slice(1) %>% pull(module)
best_mod_color <- gsub("^ME", "", best_mod_name)
kme_col <- paste0("kME", best_mod_color)

mod_proteins <- module_df$uniprot_id[module_df$module_color == best_mod_color]

# Only keep proteins that exist in kME_all rows and gs_vl rows
mod_proteins <- intersect(mod_proteins, rownames(kME_all))
mod_proteins <- intersect(mod_proteins, rownames(gs_vl))

kme_gs_df <- tibble(
  uniprot_id = mod_proteins,
  kME = kME_all[mod_proteins, kme_col],
  GS  = abs(gs_vl[mod_proteins, 1])
) %>%
  left_join(module_df %>% dplyr::select(uniprot_id, gene),
            by = "uniprot_id") %>%
  filter(!is.na(kME), !is.na(GS))

kme_thresh <- 0.7
gs_thresh  <- 0.3
n_candidates <- sum(kme_gs_df$kME > kme_thresh & kme_gs_df$GS > gs_thresh)

kme_gs_df <- kme_gs_df %>%
  mutate(is_hub = kME > kme_thresh & GS > gs_thresh,
         label  = ifelse(is_hub | rank(-kME * GS) <= 15, gene, NA_character_))

# ---- Plot ----
pD <- ggplot(kme_gs_df, aes(x = kME, y = GS)) +
  # Module background tint
  annotate("rect", xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf,
           fill = best_mod_color, alpha = 0.05) +
  geom_vline(xintercept = kme_thresh, linetype = "dashed",
             color = "grey40", linewidth = 0.3) +
  geom_hline(yintercept = gs_thresh, linetype = "dashed",
             color = "grey40", linewidth = 0.3) +
  geom_point(aes(color = is_hub, alpha = is_hub), size = 1.5) +
  scale_color_manual(values = c("FALSE" = scales::alpha(best_mod_color, 0.4),
                                 "TRUE" = best_mod_color),
                     guide = "none") +
  scale_alpha_manual(values = c("FALSE" = 0.4, "TRUE" = 0.8), guide = "none") +
  geom_text_repel(aes(label = label), size = 2.0, max.overlaps = 20,
                  segment.size = 0.2, na.rm = TRUE) +
  annotate("text",
           x = max(kme_gs_df$kME, na.rm = TRUE),
           y = max(kme_gs_df$GS, na.rm = TRUE) * 0.95,
           label = sprintf("%d candidates\n(kME > %.1f & |GS| > %.1f)",
                           n_candidates, kme_thresh, gs_thresh),
           hjust = 1, size = 2.5, fontface = "bold", color = "grey25") +
  labs(title = sprintf("D  Hub Proteins (%s Module)",
                       str_to_title(best_mod_color)),
       subtitle = "Module Membership vs |Gene Significance| for delta VL",
       x = sprintf("Module Membership (kME, %s)", best_mod_color),
       y = "|Gene Significance| (delta VL)") +
  THEME_PUB

# ---- Save ----
ggsave(file.path(RPT_DIR, "panel_D_hub_scatter.pdf"), pD,
       width = 250, height = 250, units = "mm")
ggsave(file.path(RPT_DIR, "panel_D_hub_scatter.png"), pD,
       width = 250, height = 250, units = "mm", dpi = 300)

write_csv(kme_gs_df, file.path(DAT_DIR, "04_panel_D_kme_gs.csv"))

cat("Panel D done\n")
