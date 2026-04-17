# Supplementary: WGCNA — 5-module eigengene exemplar strip
#
# Horizontal row of 5 small boxplots showing eigengene dynamics for the
# BH-significant modules, in biological-story order:
#   1. MEturquoise — Lipid Catabolism       — Aging ↑
#   2. MEbrown     — Redox / Glycolysis     — Aging ↓
#   3. MEmagenta   — Translation Initiation — Aging ↓
#   4. MEblue      — Cell Cycle / Proteost. — Training_Young ↑
#   5. MEgreen     — Oxidative Phosphoryl.  — ΔVL (within Young) ↓
#
# Each sub-panel: z-scored eigengene (y) vs Group × Time (x: Y-Pre, Y-Post,
# O-Pre, O-Post) as a boxplot with jittered subject points. Key finding tag
# is shown as a per-panel subtitle (coloured by module colour).
#
# Visual vocabulary matches the F06 panel_B_triptych.R middle sub-plot
# (boxplots + FIG_THEME + AGE_COLORS), and is used purely as a compact visual
# summary — no statistics are computed here; the significance tags are read
# from the existing pooled LMM audit.
#
# Reads:  F06/c_data/MEs.rds, meta.csv, mod_bio_labels.csv
# Writes: F06/b_reports/supp/02_analysis/SUPP_exemplar_strip.png

setwd(rprojroot::find_rstudio_root_file())
source("04_Figures/shared/style.R")

library(tidyverse)
library(patchwork)

DAT <- "04_Figures/F06/c_data"
RPT <- "04_Figures/F06/b_reports/supp/02_analysis"
dir.create(RPT, recursive = TRUE, showWarnings = FALSE)

pdf_device <- get_pdf_device()

# --- Data loading
MEs     <- readRDS(file.path(DAT, "MEs.rds"))
meta    <- read_csv(file.path(DAT, "meta.csv"),           show_col_types = FALSE)
mod_bio <- read_csv(file.path(DAT, "mod_bio_labels.csv"), show_col_types = FALSE)

meta <- meta %>%
  mutate(time = factor(time, levels = c("Pre", "Post")),
         age  = factor(age,  levels = c("Young", "Old")),
         group_time = factor(paste(age, time, sep = "-"),
                             levels = c("Young-Pre", "Young-Post",
                                        "Old-Pre",   "Old-Post")))

# --- Dynamic module selection: BH<0.05 in any LMM contrast OR |r|>0.4 anywhere
lmm_audit <- read_csv(file.path(DAT, "wgcna/wgcna_lmm_contrast_audit.csv"),
                       show_col_types = FALSE) %>%
  filter(module != "MEgrey")

# Gather all correlation/p-value files for |r|>0.4 check
pheno_files <- list(
  list(r = "wgcna_baseline_trait_correlations_young.csv",
       p = "wgcna_baseline_trait_pvalues_bh_young.csv"),
  list(r = "wgcna_baseline_trait_correlations_old.csv",
       p = "wgcna_baseline_trait_pvalues_bh_old.csv"),
  list(r = "wgcna_change_trait_correlations_young.csv",
       p = "wgcna_change_trait_pvalues_bh_young.csv"),
  list(r = "wgcna_change_trait_correlations_old.csv",
       p = "wgcna_change_trait_pvalues_bh_old.csv")
)

# Modules with BH<0.05 in any LMM contrast
sig_mods <- unique(lmm_audit$module[lmm_audit$p_bh < 0.05])
# Modules with |r|>0.4 in any LMM contrast
high_r_mods <- unique(lmm_audit$module[abs(lmm_audit$r_equiv) > 0.4])

# Also check phenotype correlations
for (pf in pheno_files) {
  rmat <- read_csv(file.path(DAT, "wgcna", pf$r), show_col_types = FALSE) %>%
    column_to_rownames("module") %>% as.matrix()
  pmat <- read_csv(file.path(DAT, "wgcna", pf$p), show_col_types = FALSE) %>%
    column_to_rownames("module") %>% as.matrix()
  for (mod in rownames(rmat)) {
    if (mod == "MEgrey") next
    if (any(pmat[mod, ] < 0.05, na.rm = TRUE)) sig_mods <- c(sig_mods, mod)
    if (any(abs(rmat[mod, ]) > 0.4, na.rm = TRUE)) high_r_mods <- c(high_r_mods, mod)
  }
}
qualifying_mods <- sort(unique(c(sig_mods, high_r_mods)))
qualifying_mods <- qualifying_mods[qualifying_mods != "MEgrey"]

# Order by biological story: aging first, then training, then phenotype
story_order <- c("MEturquoise", "MEbrown", "MEmagenta", "MEyellow",
                 "MEblue", "MEgreen", "MEred", "MEblack", "MEpink", "MEpurple")
exemplar_modules <- story_order[story_order %in% qualifying_mods]

message(sprintf("  Exemplar strip: %d qualifying modules (of 10)",
                length(exemplar_modules)))

# Auto-generate caption tags from strongest LMM association per module
caption_tags <- setNames(
  sapply(exemplar_modules, function(mod) {
    mod_lmm <- lmm_audit %>% filter(module == mod) %>%
      arrange(p_bh) %>% slice(1)
    if (nrow(mod_lmm) == 0) return("")
    stars <- sig_stars(mod_lmm$p_bh)
    ctr_label <- c(Aging = "Aging", Training_Young = "Trn_Y",
                   Training_Old = "Trn_O", Interaction = "Age\u00d7Trn")
    paste0(ctr_label[mod_lmm$contrast], " ", stars)
  }),
  exemplar_modules
)

# Biological labels — same convention as Panel A
if (!"display_label" %in% colnames(mod_bio)) {
  mod_bio <- mod_bio %>%
    mutate(display_label = paste0(bio_label, " (",
                                  stringr::str_to_title(module_color), ")"))
}
mod_display_vec <- setNames(mod_bio$display_label,
                            paste0("ME", mod_bio$module_color))

# --- Long-format ME data, filtered to the 5 exemplar modules
MEs_long <- as.data.frame(MEs) %>%
  tibble::rownames_to_column("sample_id") %>%
  pivot_longer(-sample_id, names_to = "module", values_to = "ME") %>%
  filter(module %in% exemplar_modules) %>%
  left_join(meta %>% dplyr::select(sample_id, age, time, group_time, subject),
            by = "sample_id") %>%
  group_by(module) %>%
  mutate(ME_z = as.numeric(scale(ME))) %>%
  ungroup()

# Per-module colour (WGCNA module colour for boxes)
module_color_map <- setNames(gsub("^ME", "", exemplar_modules),
                             exemplar_modules)

# --- Dimensions (dynamic: scale width to module count)
n_exemplar <- length(exemplar_modules)
STRIP_W <- max(260, n_exemplar * 46)  # 46mm per panel, min 260mm
STRIP_H <- 100  # mm tall

txt_cell    <- scale_text(BASE_GENE, STRIP_W) * 0.9
txt_axis    <- scale_text(BASE_STAT, STRIP_W) * 0.95
txt_title   <- scale_text(BASE_GENE, STRIP_W) * 1.0
txt_caption <- scale_text(BASE_STAT, STRIP_W) * 0.85

# Pre shown lighter, Post darker via alpha; single fill colour per age
build_panel <- function(mod) {
  dat <- MEs_long %>% filter(module == mod)
  mod_col <- module_color_map[[mod]]
  label   <- mod_display_vec[mod]
  if (is.na(label) || grepl("\\bNA\\b", label)) {
    col_name <- stringr::str_to_title(gsub("^ME", "", mod))
    mod_id   <- mod_bio$module_id[mod_bio$module_color == gsub("^ME", "", mod)]
    label <- if (length(mod_id) && !is.na(mod_id)) paste0(mod_id, ": ", col_name) else col_name
  }
  tag <- caption_tags[[mod]]

  ggplot(dat, aes(x = group_time, y = ME_z)) +
    geom_boxplot(aes(fill = age, alpha = time),
                 outlier.shape = NA,
                 colour = "black", linewidth = 0.35, width = 0.75) +
    geom_jitter(width = 0.09, size = 1.0,
                alpha = 0.50, colour = "grey25", na.rm = TRUE) +
    scale_fill_manual(values = c(Young = unname(AGE_COLORS["Young"]),
                                 Old   = unname(AGE_COLORS["Old"])),
                      guide = "none") +
    scale_alpha_manual(values = c(Pre = 0.55, Post = 0.95),
                       guide = "none") +
    scale_x_discrete(labels = c("Y-Pre", "Y-Post", "O-Pre", "O-Post")) +
    scale_y_continuous(expand = expansion(mult = c(0, 0))) +
    labs(x = NULL, y = NULL,
         title    = label,
         subtitle = tag) +
    FIG_THEME +
    theme(plot.title         = element_text(face = "bold",
                                             size = txt_title,
                                             colour = "grey15",
                                             lineheight = 0.9,
                                             margin = margin(b = 1)),
          plot.subtitle      = element_text(face = "bold",
                                             size = txt_caption * 1.15,
                                             colour = mod_col,
                                             margin = margin(b = 3)),
          axis.text.x        = element_text(size = txt_axis * 0.85,
                                             face = "bold",
                                             angle = 45,
                                             hjust = 1),
          axis.text.y        = element_text(size = txt_axis * 0.85),
          panel.grid.major.x = element_blank(),
          panel.grid.minor   = element_blank(),
          plot.margin        = margin(4, 4, 2, 4))
}

# --- Build panels and compose horizontally
panels <- lapply(exemplar_modules, build_panel)
n_row <- if (n_exemplar > 6) 2 else 1
n_col <- ceiling(n_exemplar / n_row)
if (n_row > 1) STRIP_H <- STRIP_H * 1.8

# Build subtitle from caption_tags
tag_summary <- paste(
  sapply(exemplar_modules, function(m) {
    col <- gsub("^ME", "", m)
    tag <- caption_tags[m]
    paste0(tag, " ", col)
  }),
  collapse = " | "
)

strip  <- wrap_plots(panels, ncol = n_col) +
  plot_annotation(
    title    = sprintf("Eigengene Exemplars \u2014 %d modules (BH<0.05 or |r|>0.4)",
                       n_exemplar),
    subtitle = paste0("Z-scored module eigengene by Group \u00d7 Time | ", tag_summary),
    theme = theme(
      plot.title    = element_text(face = "bold",
                                    size = txt_title * 1.25),
      plot.subtitle = element_text(size = txt_caption,
                                    colour = "grey30",
                                    lineheight = 0.95),
      plot.margin   = margin(4, 4, 4, 4)
    ))

ggsave(file.path(RPT, "SUPP_exemplar_strip.png"), strip,
       width = STRIP_W, height = STRIP_H, units = "mm", dpi = 300)

message(sprintf("  Saved: %s",
                file.path(RPT, "SUPP_exemplar_strip.png")))
