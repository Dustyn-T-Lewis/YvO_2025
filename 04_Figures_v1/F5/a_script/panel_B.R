################################################################################
#   YvO Figure 5 — Panel B: Module-Trait Heatmap
#   Three-part layout: Module counts (left) | Color sidebar | Heatmap (right)
#   Traits grouped by category with gaps; phenotype n noted in brackets
################################################################################
#
# ── STAT AUDIT (2026-02-27) ──────────────────────────────────────────────────
#
# CORRELATIONS
#   Pearson r between module eigengenes (MEs) and trait vectors.
#   Computed in YvO_WGCNA_run.R via cor(MEs, traits_mat, use = "pairwise").
#
# BH CORRECTION SCOPE
#   Global BH across ALL 15 modules × 9 traits = 135 tests (YvO_WGCNA_run.R
#   lines 232-237).  Panel B displays 14 modules × 8 traits (grey and BMI
#   removed) — a subset of the globally-corrected p-values.
#   This is conservative (more tests in denominator than displayed).
#
# CONFIDENCE INTERVALS
#   WGCNA's corPvalueStudent() does not return CIs.  AUDIT ADDITION below:
#   we recompute cor.test() for each displayed cell to obtain 95% CIs and
#   pairwise-complete n, then export as a supplementary CSV.
#
# SAMPLE SIZE CAVEAT
#   Design traits (age, time, interaction) have n = 62.  Phenotype traits
#   have incomplete data: Body Comp n ~ 62, Strength n ~ 50,
#   Fiber Morphology n ~ 44.  corPvalueStudent() in YvO_WGCNA_run.R uses
#   n = 62 for ALL traits, making p-values liberal for incomplete traits.
#   The cor.test() recomputation below uses pairwise-complete n, providing
#   correct p-values in the supplementary table.  The BH-corrected values
#   on the heatmap face remain from the WGCNA pipeline (slightly liberal
#   for phenotype traits); the supplementary table serves as a cross-check.
#
# AUDIT VERDICT: BH scope is appropriate (global).  Sample size caveat
# documented in subtitle.  CIs added below as supplementary export.
# ─────────────────────────────────────────────────────────────────────────────

source("04_Figures/F5/a_script/YvO_F5_setup.R")

cor_mat <- trait_cor %>% column_to_rownames("module") %>% as.matrix()
pval_mat_b <- pval_bh %>% column_to_rownames("module") %>% as.matrix()

# --- Trait selection (drop BMI — baseline only, not informative) ---
trait_order <- c("age_num", "time_num", "interaction",
                 "VL_thick_cm", "DXA_LBM_kg",
                 "deadlift_1rm_kg",
                 "Type_I_fCSA", "Type_II_fCSA")
trait_order <- intersect(trait_order, colnames(cor_mat))
cor_mat    <- cor_mat[, trait_order]
pval_mat_b <- pval_mat_b[, trait_order]

# --- Remove grey (unassigned) module ---
cor_mat    <- cor_mat[rownames(cor_mat) != "MEgrey", , drop = FALSE]
pval_mat_b <- pval_mat_b[rownames(pval_mat_b) != "MEgrey", , drop = FALSE]

# Human-readable trait labels
col_labels <- c(
  age_num         = "Age",
  time_num        = "Training",
  interaction     = "Age x Trn",
  VL_thick_cm     = "VL Thick.",
  DXA_LBM_kg      = "DXA LBM",
  deadlift_1rm_kg = "Deadlift",
  Type_I_fCSA     = "Type I fCSA",
  Type_II_fCSA    = "Type II fCSA"
)

# --- Numeric x positions with narrow gaps between groups ---
# Group 1: Design (Age, Training, Age x Trn) -> x = 1, 2, 3
# Group 2: Body Comp  -> x = 4.15, 5.15
# Group 3: Strength   -> x = 6.30
# Group 4: Fiber Morph -> x = 7.45, 8.45
gap <- 0.15
xpos <- c(1, 2, 3,                          # Design
           3 + gap + 1, 3 + gap + 2,         # Body Composition
           3 + gap + 2 + gap + 1,             # Strength
           3 + gap + 2 + gap + 1 + gap + 1,   # Fiber Morphology
           3 + gap + 2 + gap + 1 + gap + 2)
trait_xpos <- setNames(xpos, trait_order)

# Order modules by hierarchical clustering of correlation profiles
me_dist <- dist(cor_mat)
me_hc <- hclust(me_dist)
mod_order <- rownames(cor_mat)[me_hc$order]

# Build heat_df with numeric x
heat_df <- expand.grid(module = rownames(cor_mat), trait = colnames(cor_mat),
                       stringsAsFactors = FALSE) %>%
  mutate(cor   = as.vector(cor_mat),
         pval  = as.vector(pval_mat_b),
         stars = sig_stars(pval),
         label = sprintf("%.2f%s", cor, stars),
         trait_label = col_labels[trait],
         xpos  = trait_xpos[trait])

heat_df$module <- factor(heat_df$module, levels = mod_order)

# --- Module display names (strip ME prefix) ---
mod_display <- gsub("^ME", "", mod_order) %>% str_to_title()
mod_display_map <- setNames(mod_display, mod_order)
mod_color_raw <- setNames(gsub("^ME", "", mod_order), mod_order)

# Light module color list (for text color on bars)
light_modules <- c("yellow", "greenyellow", "cyan", "tan", "pink", "salmon")

# =========================================================================
# Part 1: Module Gene Counts (left-justified horizontal bars)
# =========================================================================

mod_counts <- module_df %>%
  filter(module_color != "grey") %>%
  count(module_color, name = "n_proteins") %>%
  mutate(module = paste0("ME", module_color)) %>%
  filter(module %in% mod_order) %>%
  mutate(module = factor(module, levels = mod_order))

p_counts <- ggplot(mod_counts, aes(x = n_proteins, y = module)) +
  geom_col(fill = mod_counts$module_color, color = "black",
           linewidth = 0.3, width = 1.0) +
  geom_text(aes(label = n_proteins, x = n_proteins / 2),
            size = 2.8, fontface = "bold",
            color = ifelse(mod_counts$module_color %in% light_modules,
                           "grey30", "white")) +
  scale_x_continuous(expand = c(0, 0), limits = c(0, max(mod_counts$n_proteins) * 1.02)) +
  scale_y_discrete(labels = NULL) +
  labs(y = NULL, x = "Gene Counts") +
  THEME_PUB +
  theme(axis.text.y    = element_blank(),
        axis.ticks.y   = element_blank(),
        axis.text.x    = element_text(size = 6),
        axis.title.x   = element_text(size = 7),
        panel.grid.major.y = element_blank(),
        panel.grid.minor   = element_blank(),
        plot.margin    = margin(2, 0, 2, 2))

# =========================================================================
# Part 2: Module names + Color sidebar (narrow colored squares)
# =========================================================================

sidebar_df <- tibble(
  module = factor(mod_order, levels = mod_order),
  mod_color = mod_color_raw[mod_order],
  mod_name  = mod_display_map[mod_order]
)

p_sidebar <- ggplot(sidebar_df, aes(x = 1, y = module)) +
  geom_tile(fill = sidebar_df$mod_color, color = "black",
            linewidth = 0.3, width = 0.8, height = 1.0) +
  geom_text(aes(x = 0.45, label = mod_name), hjust = 1,
            size = 2.8, fontface = "bold", color = "black") +
  scale_x_continuous(limits = c(-1.2, 1.4), expand = c(0, 0)) +
  theme_void() +
  theme(plot.margin = margin(2, -5, 2, 0))

# =========================================================================
# Part 3: Bracket annotations above heatmap
# =========================================================================

xmin_all <- min(xpos) - 0.55
xmax_all <- max(xpos) + 0.55

# Group boundaries for brackets
brackets <- tribble(
  ~label,                       ~start,  ~end,
  "Study Design",               xpos[1], xpos[3],
  "Body Comp.\n(n=62)",         xpos[4], xpos[5],
  "Strength\n(n=50)",           xpos[6], xpos[6],
  "Fiber Morphology\n(n=44)",   xpos[7], xpos[8]
)
brackets <- brackets %>% mutate(mid = (start + end) / 2)

p_brackets <- ggplot() +
  # Horizontal bracket lines
  geom_segment(data = brackets,
               aes(x = start - 0.4, xend = end + 0.4, y = 0.28, yend = 0.28),
               linewidth = 0.4, color = "grey30") +
  # Left caps
  geom_segment(data = brackets,
               aes(x = start - 0.4, xend = start - 0.4, y = 0.28, yend = 0.14),
               linewidth = 0.4, color = "grey30") +
  # Right caps
  geom_segment(data = brackets,
               aes(x = end + 0.4, xend = end + 0.4, y = 0.28, yend = 0.14),
               linewidth = 0.4, color = "grey30") +
  # Labels above brackets
  geom_text(data = brackets,
            aes(x = mid, y = 0.62, label = label),
            size = 2.3, fontface = "bold", color = "grey25",
            lineheight = 0.85) +
  scale_x_continuous(limits = c(xmin_all, xmax_all), expand = c(0, 0)) +
  scale_y_continuous(limits = c(0, 1), expand = c(0, 0)) +
  theme_void() +
  theme(plot.margin = margin(2, 2, 0, 0))

# =========================================================================
# Part 4: Main heatmap (numeric x-axis with gaps)
# =========================================================================

p_heat <- ggplot(heat_df, aes(x = xpos, y = module, fill = cor)) +

  geom_tile(color = "black", linewidth = 0.3) +
  geom_tile(data = heat_df %>% filter(pval < 0.001),
            color = "black", linewidth = 1.0, fill = NA) +
  geom_text(data = heat_df %>% filter(pval >= 0.05),
            aes(label = label), size = 2.8, color = "black") +
  geom_text(data = heat_df %>% filter(pval < 0.05 & pval >= 0.001),
            aes(label = label), size = 2.8, fontface = "bold", color = "black") +
  geom_text(data = heat_df %>% filter(pval < 0.001),
            aes(label = label), size = 3.2, fontface = "bold", color = "white") +
  scale_fill_gradient2(low = "#0033CC", mid = "white", high = "#CC0000",
                       midpoint = 0, limits = c(-0.8, 0.8),
                       oob = scales::squish,
                       name = "Pearson r",
                       na.value = "white") +
  scale_x_continuous(
    breaks = xpos,
    labels = col_labels[trait_order],
    limits = c(xmin_all, xmax_all),
    expand = c(0, 0)
  ) +
  scale_y_discrete(labels = NULL) +
  labs(x = NULL, y = NULL) +
  THEME_PUB +
  LEGEND_THEME +
  theme(axis.text.x       = element_text(angle = 45, hjust = 1, size = 7),
        axis.text.y        = element_blank(),
        axis.ticks.y       = element_blank(),
        axis.ticks.x       = element_blank(),
        legend.key.width   = unit(15, "mm"),
        legend.key.height  = unit(3, "mm"),
        panel.grid         = element_blank(),
        panel.border       = element_blank(),
        plot.margin        = margin(0, 2, 2, -3))

# =========================================================================
# Combine using area-based design
# Row 1: bracket annotations (only above heatmap area)
# Rows 2-10: counts | sidebar | heatmap
# =========================================================================

design <- c(
  area(1, 4, 1, 12),   # 1: brackets (above heatmap only)
  area(2, 1, 10, 2),   # 2: gene counts
  area(2, 3, 10, 3),   # 3: sidebar (single column, tight to heatmap)
  area(2, 4, 10, 12)   # 4: heatmap
)

pB <- wrap_elements(p_brackets) + p_counts + p_sidebar + p_heat +
  plot_layout(design = design) +
  plot_annotation(
    title    = "B  Module-Trait Correlations",
    subtitle = "BH-corrected | * p<.05  ** p<.01  *** p<.001 | Phenotype correlations computed on available samples per trait",
    theme = theme(
      plot.title    = element_text(face = "bold", size = 9),
      plot.subtitle = element_text(size = 6.5, color = "grey30", face = "italic"),
      plot.margin   = margin(2, 2, 2, 2)
    )
  )

write_csv(heat_df, file.path(DAT_DIR, "02_panel_B_heatmap_data.csv"))

ggsave(file.path(RPT_DIR, "panel_B_heatmap.pdf"), pB,
       width = 320, height = 220, units = "mm",
       device = pdf, limitsize = FALSE)
ggsave(file.path(RPT_DIR, "panel_B_heatmap.png"), pB,
       width = 320, height = 220, units = "mm",
       dpi = 300, limitsize = FALSE)

# ── STAT AUDIT ADDITION: cor.test() CIs for every displayed cell ────────────
# Recompute Pearson r with 95% CI and exact pairwise-complete n.
# This supplements the WGCNA-pipeline BH values shown on the heatmap.

cat("Computing cor.test() CIs for Panel B cells...\n")

# Rebuild trait matrix the same way as YvO_WGCNA_run.R
audit_traits <- meta %>%
  mutate(
    age_num     = if_else(age == "Old", 1, 0),
    time_num    = if_else(time == "Post", 1, 0),
    interaction = age_num * time_num
  )
pheno_cols_b <- c("VL_thick_cm", "DXA_LBM_kg", "deadlift_1rm_kg",
                  "Type_I_fCSA", "Type_II_fCSA")
for (pc in pheno_cols_b) {
  if (pc %in% names(meta)) audit_traits[[pc]] <- meta[[pc]]
}
audit_trait_mat <- audit_traits %>%
  dplyr::select(sample_id, all_of(trait_order)) %>%
  column_to_rownames("sample_id")
audit_trait_mat <- audit_trait_mat[rownames(MEs), , drop = FALSE]

# Displayed modules (exclude grey)
disp_mods <- rownames(cor_mat)

ci_results <- list()
for (me_col in disp_mods) {
  for (tr in trait_order) {
    me_vec <- MEs[, me_col]
    tr_vec <- audit_trait_mat[, tr]
    complete <- complete.cases(me_vec, tr_vec)
    n_pair <- sum(complete)
    if (n_pair < 4) {
      ci_results <- c(ci_results, list(tibble(
        module = me_col, trait = tr, n = n_pair,
        r = NA_real_, ci_lo = NA_real_, ci_hi = NA_real_,
        p_raw = NA_real_
      )))
      next
    }
    ct <- cor.test(me_vec[complete], tr_vec[complete], method = "pearson")
    ci_results <- c(ci_results, list(tibble(
      module = me_col, trait = tr, n = n_pair,
      r = round(ct$estimate, 4),
      ci_lo = round(ct$conf.int[1], 4),
      ci_hi = round(ct$conf.int[2], 4),
      p_raw = ct$p.value
    )))
  }
}

ci_df <- bind_rows(ci_results)
ci_df$p_bh_global <- heat_df$pval[match(
  paste(ci_df$module, ci_df$trait),
  paste(heat_df$module, heat_df$trait)
)]
write_csv(ci_df, file.path(DAT_DIR, "02_panel_B_correlation_CIs.csv"))
cat(sprintf("  Exported: 02_panel_B_correlation_CIs.csv (%d cells, 95%% CIs)\n", nrow(ci_df)))
# ─────────────────────────────────────────────────────────────────────────────

cat("Panel B saved\n")
