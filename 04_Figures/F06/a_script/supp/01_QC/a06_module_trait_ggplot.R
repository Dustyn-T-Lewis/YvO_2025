# Supplementary: Module-Trait Correlation Heatmap (ggplot rewrite)
#
# Replaces the base-R corrplot module-trait heatmap from YvO_WGCNA_run.R
# with a FIG_THEME-styled ggplot version for publication consistency.
#
# Uses pre-computed correlation and BH-adjusted p-value CSVs from c_data/wgcna/.
# Displays biological module labels and cleaned trait names.
#
# Generates: 01_QC/SUPP_module_trait_heatmap.png

setwd(rprojroot::find_rstudio_root_file())
source("04_Figures/shared/style.R")

library(readr)
library(dplyr)
library(tidyr)
library(tibble)
library(stringr)

DAT <- "04_Figures/F06/c_data"
RPT <- "04_Figures/F06/b_reports/supp/01_QC"
dir.create(RPT, recursive = TRUE, showWarnings = FALSE)

pdf_device <- get_pdf_device()

# --- Load data
cor_df <- read_csv(file.path(DAT, "wgcna/wgcna_module_trait_correlations.csv"),
                    show_col_types = FALSE)
pval_df <- read_csv(file.path(DAT, "wgcna/wgcna_module_trait_pvalues_bh.csv"),
                     show_col_types = FALSE)
mod_bio <- read_csv(file.path(DAT, "mod_bio_labels.csv"), show_col_types = FALSE)

# --- Build display labels for modules
if (!"display_label" %in% colnames(mod_bio)) {
  mod_bio <- mod_bio %>%
    mutate(display_label = ifelse(is.na(bio_label),
                                  str_to_title(module_color),
                                  paste0(bio_label, " (", str_to_title(module_color), ")")))
}
mod_label_map <- setNames(mod_bio$display_label, paste0("ME", mod_bio$module_color))

# Order modules by protein count (largest first, grey last)
mod_order <- mod_bio %>%
  arrange(desc(n_proteins)) %>%
  pull(module_color)
me_order <- paste0("ME", c(mod_order, "grey"))

# --- Clean trait names for display
trait_labels <- c(
  age_num        = "Age",
  time_num       = "Time",
  interaction    = "Age \u00d7 Time",
  VL_thick_cm    = "VL Thickness",
  DXA_LBM_kg    = "Lean Body Mass",
  BMI            = "BMI",
  Type_I_fCSA    = "Type I fCSA",
  Type_II_fCSA   = "Type II fCSA",
  deadlift_1rm_kg = "Deadlift 1RM"
)

# --- Pivot to long format
cor_long <- cor_df %>%
  pivot_longer(-module, names_to = "trait", values_to = "r")

pval_long <- pval_df %>%
  pivot_longer(-module, names_to = "trait", values_to = "p_bh")

heat_df <- cor_long %>%
  left_join(pval_long, by = c("module", "trait")) %>%
  mutate(
    module_label = ifelse(module %in% names(mod_label_map),
                          mod_label_map[module],
                          gsub("^ME", "", module) %>% str_to_title()),
    trait_label  = ifelse(trait %in% names(trait_labels),
                          trait_labels[trait], trait),
    sig_label    = sig_stars(p_bh),
    # Text color: white on dark cells, grey20 on light cells
    text_col     = ifelse(abs(r) > 0.4, "white", "grey20"),
    # Cell text: correlation value + stars
    cell_text    = sprintf("%.2f", r),
    cell_full    = ifelse(sig_label == "", cell_text,
                          paste0(cell_text, "\n", sig_label))
  )

# Factor ordering
heat_df <- heat_df %>%
  mutate(
    module = factor(module, levels = rev(me_order)),
    module_label = factor(module_label,
                          levels = rev(c(mod_label_map[me_order[me_order %in% names(mod_label_map)]],
                                         "Grey"))),
    trait_label = factor(trait_label,
                         levels = trait_labels[names(trait_labels) %in% unique(heat_df$trait)])
  )

# --- Panel dimensions
PA_W <- 200
PA_H <- 160

txt_cell <- scale_text(BASE_GENE, PA_W) * 0.85
txt_star <- scale_text(BASE_GENE, PA_W) * 0.7
txt_axis <- scale_text(BASE_STAT, PA_W)

# --- Module color sidebar data
mod_color_df <- heat_df %>%
  distinct(module, module_label) %>%
  mutate(color_name = gsub("^ME", "", as.character(module)))

# --- Build heatmap
p <- ggplot(heat_df, aes(x = trait_label, y = module_label, fill = r)) +
  geom_tile(color = "white", linewidth = 0.4) +
  geom_text(aes(label = cell_text, color = text_col),
            size = txt_cell, fontface = "plain", vjust = 0.35) +
  geom_text(aes(label = sig_label, color = text_col),
            size = txt_star, fontface = "bold", vjust = 1.8) +
  scale_fill_gradient2(
    low = "#4393C3", mid = "white", high = "#D6604D",
    midpoint = 0, limits = c(-1, 1),
    name = "Pearson r",
    guide = guide_colorbar(barwidth = unit(4, "cm"), barheight = unit(0.35, "cm"),
                           title.position = "left",
                           title.theme = element_text(size = 7, face = "bold",
                                                       vjust = 0.8))
  ) +
  scale_color_identity() +
  labs(x = NULL, y = NULL,
       title = "Module\u2013Trait Correlations",
       subtitle = "Pearson correlation | * BH < 0.05  ** BH < 0.01  *** BH < 0.001") +
  FIG_THEME +
  theme(
    axis.text.x  = element_text(angle = 40, hjust = 1, size = txt_axis),
    axis.text.y  = element_text(size = txt_axis * 0.9),
    legend.position = "bottom",
    legend.margin = margin(t = -5),
    panel.grid    = element_blank()
  )

ggsave(file.path(RPT, "SUPP_module_trait_heatmap.png"), p,
       width = PA_W, height = PA_H, units = "mm",
       dpi = 300, limitsize = FALSE)

message("  Module-trait heatmap ggplot saved")
