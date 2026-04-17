# Supplementary — UpSet Plot: WGCNA Module x DEP List Overlap
# Shows how DEPs from each contrast distribute across WGCNA modules.
# Bridges the DEP-level analysis (F04/F05) with module-level analysis (F06).

setwd(rprojroot::find_rstudio_root_file())
source("04_Figures/shared/style.R")

library(readr)
library(dplyr)
library(tidyr)
library(ggplot2)
library(patchwork)


RPT <- "04_Figures/F06/b_reports/supp/02_analysis"
DAT <- "04_Figures/F06/c_data"
dir.create(RPT, recursive = TRUE, showWarnings = FALSE)

pdf_device <- get_pdf_device()

message("Supp e02: UpSet — DEP x module overlap...")

# --- Load DEP data from canonical combined results
dep_df <- read_csv("03_DEP/c_data/03_combined_results.csv", show_col_types = FALSE)

extract_dep <- function(df, contrast, set_prefix) {
  fdr_col <- paste0("adj.P.Val_", contrast)
  lfc_col <- paste0("logFC_", contrast)
  df %>%
    filter(!is.na(.data[[fdr_col]]) & .data[[fdr_col]] < 0.05) %>%
    transmute(gene,
              set = paste0(set_prefix, ifelse(.data[[lfc_col]] > 0, "Up", "Down")))
}

vol_aging <- extract_dep(dep_df, "Aging", "Aging_")
vol_ty    <- extract_dep(dep_df, "Training_Young", "Trn_Young_")
vol_to    <- extract_dep(dep_df, "Training_Old", "Trn_Old_")
vol_int   <- dep_df %>%
  filter(!is.na(adj.P.Val_Interaction) & adj.P.Val_Interaction < 0.05) %>%
  transmute(gene, set = "Interaction")

# --- Load module assignments
mod_assign <- read_csv(file.path(DAT, "wgcna", "wgcna_module_assignments.csv"),
                        show_col_types = FALSE)
mod_bio <- read_csv(file.path(DAT, "mod_bio_labels.csv"), show_col_types = FALSE)

# Combine DEPs with module membership
all_deps <- bind_rows(vol_aging, vol_ty, vol_to, vol_int)
dep_genes <- unique(all_deps$gene)

dep_modules <- mod_assign %>%
  filter(gene %in% dep_genes) %>%
  select(gene, module_color)

# Build binary membership matrix
dep_sets <- all_deps %>%
  distinct(gene, set) %>%
  mutate(present = 1L) %>%
  pivot_wider(names_from = set, values_from = present, values_fill = 0L) %>%
  left_join(dep_modules, by = "gene")

# --- Module color mapping
mod_hex <- c(
  blue = "#4393C3", turquoise = "#40E0D0", yellow = "#DAA520",
  green = "#2E8B57", brown = "#8B4513", red = "#D6604D",
  black = "#333333", pink = "#FF69B4", magenta = "#8B008B",
  grey = "#999999"
)

# --- Panel dimensions
PU_W <- 260
PU_H <- 220

txt_axis  <- scale_text(BASE_STAT, PU_W) * 1.0
txt_title <- scale_text(BASE_GENE, PU_W) * 1.3

# --- Build stacked bar: DEP counts per module, by contrast
dep_long <- all_deps %>%
  left_join(dep_modules, by = "gene") %>%
  filter(!is.na(module_color)) %>%
  count(module_color, set) %>%
  mutate(module_color = factor(module_color,
    levels = names(sort(tapply(n, module_color, sum), decreasing = TRUE))))

# Shorter display labels for DEP categories
set_labels <- c(
  Aging_Up       = "Aging Up",
  Aging_Down     = "Aging Dn",
  Trn_Young_Up   = "TY Up",
  Trn_Young_Down = "TY Dn",
  Trn_Old_Up     = "TO Up",
  Trn_Old_Down   = "TO Dn",
  Interaction    = "Int"
)

# Biological labels for modules
if (!is.null(mod_bio)) {
  mod_label_map <- setNames(mod_bio$display_label, mod_bio$module_color)
  dep_long <- dep_long %>%
    mutate(module_label = ifelse(module_color %in% names(mod_label_map),
                                 mod_label_map[as.character(module_color)],
                                 as.character(module_color)))
  dep_long$module_label <- factor(dep_long$module_label,
    levels = unique(dep_long$module_label[order(match(dep_long$module_color,
                                                       levels(dep_long$module_color)))]))
} else {
  dep_long$module_label <- dep_long$module_color
}

p_bar <- ggplot(dep_long, aes(x = module_label, y = n, fill = set)) +
  geom_col(position = "stack", width = 0.7, color = "black", linewidth = 0.2) +
  scale_fill_manual(
    labels = set_labels,
    values = c(
      Aging_Up        = scales::alpha("#D6604D", 0.8),
      Aging_Down      = scales::alpha("#4393C3", 0.8),
      Trn_Young_Up    = scales::alpha("#E05A4E", 0.5),
      Trn_Young_Down  = scales::alpha("#5DA5DA", 0.5),
      Trn_Old_Up      = scales::alpha("#E05A4E", 0.25),
      Trn_Old_Down    = scales::alpha("#5DA5DA", 0.25),
      Interaction     = "#9B7FBF"
    ),
    name = "DEP Category"
  ) +
  labs(x = "WGCNA Module", y = "Number of DEPs",
       title = "DEP Distribution Across WGCNA Modules") +
  FIG_THEME +
  theme(
    plot.title    = element_text(size = txt_title, face = "bold"),
    axis.title    = element_text(size = txt_axis),
    axis.text     = element_text(size = txt_axis * 0.85),
    axis.text.x   = element_text(angle = 45, hjust = 1),
    legend.position = "right",
    legend.text   = element_text(size = txt_axis * 0.75)
  )

# --- Build module proportion by contrast (heatmap-style)
dep_prop <- all_deps %>%
  left_join(dep_modules, by = "gene") %>%
  filter(!is.na(module_color)) %>%
  count(set, module_color) %>%
  group_by(set) %>%
  mutate(prop = n / sum(n)) %>%
  ungroup()

# Order modules by total DEP count
mod_order <- dep_prop %>%
  group_by(module_color) %>%
  summarise(total = sum(n), .groups = "drop") %>%
  arrange(desc(total)) %>%
  pull(module_color)

# Add biological labels to proportion data
if (!is.null(mod_bio)) {
  dep_prop <- dep_prop %>%
    mutate(module_label = ifelse(module_color %in% names(mod_label_map),
                                  mod_label_map[module_color],
                                  module_color))
} else {
  dep_prop$module_label <- dep_prop$module_color
}

dep_prop <- dep_prop %>%
  mutate(module_label = factor(module_label,
           levels = rev(unique(module_label[order(match(module_color, mod_order))]))),
         set = factor(set, levels = c("Aging_Up", "Aging_Down",
                                       "Trn_Young_Up", "Trn_Young_Down",
                                       "Trn_Old_Up", "Trn_Old_Down",
                                       "Interaction")))

p_heat <- ggplot(dep_prop, aes(x = set, y = module_label)) +
  geom_tile(aes(fill = prop), color = "white", linewidth = 0.5) +
  geom_text(aes(label = n), size = txt_axis * 0.85, color = "grey20") +
  scale_fill_gradient(low = "white", high = "#2166AC", name = "Proportion",
                      limits = c(0, NA)) +
  scale_x_discrete(labels = set_labels) +
  labs(x = NULL, y = "Module",
       title = "Module Enrichment per DEP Category") +
  FIG_THEME +
  theme(
    plot.title    = element_text(size = txt_title, face = "bold"),
    axis.title    = element_text(size = txt_axis),
    axis.text     = element_text(size = txt_axis * 0.85),
    axis.text.x   = element_text(angle = 40, hjust = 1),
    legend.position = "right"
  )

# --- Composite
composite <- p_bar / p_heat + plot_layout(heights = c(0.45, 0.55))

ggsave(file.path(RPT, "SUPP_upset_dep_modules.png"), composite,
       width = PU_W, height = PU_H, units = "mm",
       dpi = 300, limitsize = FALSE)

# --- Export data
write_csv(dep_long, file.path(DAT, "supp", "e02_dep_module_counts.csv"))
write_csv(dep_prop, file.path(DAT, "supp", "e02_dep_module_proportions.csv"))

message("  e02_upset_dep_modules saved")
