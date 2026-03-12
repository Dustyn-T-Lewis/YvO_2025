################################################################################
#   Figure 5 — Panel B: Module-Trait Heatmap
#   Self-contained: loads own data from c_data/
#   Three-part layout: Module counts (left) | Color sidebar | Heatmap (right)
#   Traits grouped by category with gaps; phenotype n noted in brackets
#   Generates: panel_B_heatmap.pdf/png, c_data/02_panel_B_*.csv
################################################################################

setwd(rprojroot::find_rstudio_root_file())
source("04_Figures_v2/shared/style.R")

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(stringr)
  library(patchwork)
})

RPT <- "04_Figures_v2/F5/b_reports"
DAT <- "04_Figures_v2/F5/c_data"
dir.create(RPT, recursive = TRUE, showWarnings = FALSE)
dir.create(DAT, recursive = TRUE, showWarnings = FALSE)

pdf_device <- get_pdf_device()

# --- Load data ---
trait_cor      <- read_csv(file.path(DAT, "wgcna/wgcna_module_trait_correlations.csv"), show_col_types = FALSE)
pval_bh        <- read_csv(file.path(DAT, "wgcna/wgcna_module_trait_pvalues_bh.csv"), show_col_types = FALSE)
module_df      <- read_csv(file.path(DAT, "wgcna/wgcna_module_assignments.csv"), show_col_types = FALSE)
MEs            <- readRDS(file.path(DAT, "MEs.rds"))
meta           <- read_csv(file.path(DAT, "meta.csv"), show_col_types = FALSE)
mod_bio_labels <- read_csv(file.path(DAT, "mod_bio_labels.csv"), show_col_types = FALSE)

# --- mod_display_label helper ---
mod_bio_labels_vec <- setNames(mod_bio_labels$bio_label, mod_bio_labels$module_color)
mod_display_label <- function(color) {
  lbl <- mod_bio_labels_vec[color]
  lbl[is.na(lbl)] <- stringr::str_to_title(color[is.na(lbl)])
  paste0(lbl, " (", stringr::str_to_title(color), ")")
}

message("Panel B: module-trait heatmap...")

PB_W <- 320  # panel width mm
PB_H <- 220  # panel height mm

# --- Text sizes ---
txt_cell   <- scale_text(BASE_GENE, PB_W)
txt_count  <- scale_text(BASE_COUNT, PB_W)
txt_name   <- scale_text(BASE_GENE, PB_W)
txt_axis   <- scale_text(BASE_STAT, PB_W)
txt_brack  <- scale_text(BASE_GENE, PB_W) * 0.8

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
gap <- 0.15
xpos <- c(1, 2, 3,
           3 + gap + 1, 3 + gap + 2,
           3 + gap + 2 + gap + 1,
           3 + gap + 2 + gap + 1 + gap + 1,
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

# --- Module display names (biological labels) ---
mod_color_raw <- setNames(gsub("^ME", "", mod_order), mod_order)
mod_display <- mod_display_label(mod_color_raw)
mod_display_map <- setNames(mod_display, mod_order)

light_modules <- c("yellow", "greenyellow", "cyan", "tan", "pink", "salmon")

# =========================================================================
# Part 1: Module Gene Counts
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
            size = txt_count, fontface = "bold",
            color = ifelse(mod_counts$module_color %in% light_modules,
                           "grey30", "white")) +
  scale_x_continuous(expand = c(0, 0), limits = c(0, max(mod_counts$n_proteins) * 1.02)) +
  scale_y_discrete(labels = NULL) +
  labs(y = NULL, x = "Gene Counts") +
  FIG_THEME +
  theme(axis.text.y        = element_blank(),
        axis.ticks.y       = element_blank(),
        axis.text.x        = element_text(size = txt_axis),
        axis.title.x       = element_text(size = txt_axis),
        panel.grid.major.y = element_blank(),
        panel.grid.minor   = element_blank(),
        legend.position    = "none",
        plot.margin        = margin(2, 0, 2, 2))

# =========================================================================
# Part 2: Module names + Color sidebar
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
            size = txt_name, fontface = "bold", color = "black") +
  scale_x_continuous(limits = c(-1.2, 1.4), expand = c(0, 0)) +
  theme_void() +
  theme(plot.margin = margin(2, -5, 2, 0))

# =========================================================================
# Part 3: Bracket annotations above heatmap
# =========================================================================

xmin_all <- min(xpos) - 0.55
xmax_all <- max(xpos) + 0.55

brackets <- tribble(
  ~label,                       ~start,  ~end,
  "Study Design",               xpos[1], xpos[3],
  "Body Comp.\n(n=62)",         xpos[4], xpos[5],
  "Strength\n(n=50)",           xpos[6], xpos[6],
  "Fiber Morphology\n(n=44)",   xpos[7], xpos[8]
)
brackets <- brackets %>% mutate(mid = (start + end) / 2)

p_brackets <- ggplot() +
  geom_segment(data = brackets,
               aes(x = start - 0.4, xend = end + 0.4, y = 0.28, yend = 0.28),
               linewidth = 0.4, color = "grey30") +
  geom_segment(data = brackets,
               aes(x = start - 0.4, xend = start - 0.4, y = 0.28, yend = 0.14),
               linewidth = 0.4, color = "grey30") +
  geom_segment(data = brackets,
               aes(x = end + 0.4, xend = end + 0.4, y = 0.28, yend = 0.14),
               linewidth = 0.4, color = "grey30") +
  geom_text(data = brackets,
            aes(x = mid, y = 0.62, label = label),
            size = txt_brack, fontface = "bold", color = "grey25",
            lineheight = 0.85) +
  scale_x_continuous(limits = c(xmin_all, xmax_all), expand = c(0, 0)) +
  scale_y_continuous(limits = c(0, 1), expand = c(0, 0)) +
  theme_void() +
  theme(plot.margin = margin(2, 2, 0, 0))

# =========================================================================
# Part 4: Main heatmap
# =========================================================================

p_heat <- ggplot(heat_df, aes(x = xpos, y = module, fill = cor)) +
  geom_tile(color = "black", linewidth = 0.3) +
  geom_tile(data = heat_df %>% filter(pval < 0.001),
            color = "black", linewidth = 1.0, fill = NA) +
  geom_text(data = heat_df %>% filter(pval >= 0.05),
            aes(label = label), size = txt_cell, color = "black") +
  geom_text(data = heat_df %>% filter(pval < 0.05 & pval >= 0.001),
            aes(label = label), size = txt_cell, fontface = "bold", color = "black") +
  geom_text(data = heat_df %>% filter(pval < 0.001),
            aes(label = label), size = txt_cell * 1.1, fontface = "bold", color = "white") +
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
  FIG_THEME +
  theme(axis.text.x        = element_text(angle = 45, hjust = 1, size = txt_axis),
        axis.text.y        = element_blank(),
        axis.ticks.y       = element_blank(),
        axis.ticks.x       = element_blank(),
        panel.grid         = element_blank(),
        panel.border       = element_blank(),
        legend.position    = "none",
        plot.margin        = margin(0, 2, 2, -3))

# =========================================================================
# Combine using area-based design
# =========================================================================

design <- c(
  area(1, 4, 1, 12),   # 1: brackets (above heatmap only)
  area(2, 1, 10, 2),   # 2: gene counts
  area(2, 3, 10, 3),   # 3: sidebar
  area(2, 4, 10, 12)   # 4: heatmap
)

pB <- wrap_elements(p_brackets) + p_counts + p_sidebar + p_heat +
  plot_layout(design = design) +
  plot_annotation(
    title    = "Module-Trait Correlations",
    subtitle = "BH-corrected | * p<.05  ** p<.01  *** p<.001 | Phenotype correlations computed on available samples per trait",
    theme = theme(
      plot.title    = element_text(face = "bold", size = txt_axis * 2.5),
      plot.subtitle = element_text(size = txt_cell * 2, color = "grey30",
                                   face = "italic"),
      plot.margin   = margin(2, 2, 2, 2)
    )
  )

write_csv(heat_df, file.path(DAT, "02_panel_B_heatmap_data.csv"))

ggsave(file.path(RPT, "panel_B_heatmap.pdf"), pB,
       width = PB_W, height = PB_H, units = "mm",
       device = pdf_device, limitsize = FALSE)
ggsave(file.path(RPT, "panel_B_heatmap.png"), pB,
       width = PB_W, height = PB_H, units = "mm",
       dpi = 300, limitsize = FALSE)

# cor.test() CIs for every displayed cell (stat audit)
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

disp_mods <- rownames(cor_mat)

ci_df <- expand.grid(module = disp_mods, trait = trait_order,
                     stringsAsFactors = FALSE) %>%
  purrr::pmap_dfr(function(module, trait) {
    me_vec <- MEs[, module]; tr_vec <- audit_trait_mat[, trait]
    ok <- complete.cases(me_vec, tr_vec); n_pair <- sum(ok)
    if (n_pair < 4) return(tibble(module = module, trait = trait, n = n_pair,
                                   r = NA_real_, ci_lo = NA_real_, ci_hi = NA_real_, p_raw = NA_real_))
    ct <- cor.test(me_vec[ok], tr_vec[ok], method = "pearson")
    tibble(module = module, trait = trait, n = n_pair,
           r = round(ct$estimate, 4), ci_lo = round(ct$conf.int[1], 4),
           ci_hi = round(ct$conf.int[2], 4), p_raw = ct$p.value)
  })
ci_df$p_bh_global <- heat_df$pval[match(
  paste(ci_df$module, ci_df$trait),
  paste(heat_df$module, heat_df$trait)
)]
write_csv(ci_df, file.path(DAT, "02_panel_B_correlation_CIs.csv"))

message("  Panel B saved")
