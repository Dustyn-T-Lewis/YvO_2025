# Figure 5 — Panel C: Per-Module Triptych
# Layout: Accent | Mini Heatmap | Sigmoid Sankey | PEA Bars | Eigengene Strip
# Also computes: module reversal (paired t + Wilcoxon), age gap closure
#
# STAT AUDIT: Eigengene reversal = paired t + Shapiro-Wilk + Wilcoxon sensitivity,
# BH across age modules. Age gap closure = Wilcoxon signed-rank + Hodges-Lehmann CI.
# Sankey = greedy 1:1 (non-inferential). GO bars = BH-corrected ORA.

setwd(rprojroot::find_rstudio_root_file())
source("04_Figures/shared/style.R")

suppressPackageStartupMessages({
  library(tidyverse)
  library(patchwork)
  library(ggrepel)
})

RPT <- "04_Figures/F5/b_reports"
DAT <- "04_Figures/F5/c_data"
dir.create(RPT, recursive = TRUE, showWarnings = FALSE)
dir.create(DAT, recursive = TRUE, showWarnings = FALSE)

pdf_device <- get_pdf_device()

meta       <- read_csv(file.path(DAT, "meta.csv"), show_col_types = FALSE)
meta$group <- factor(meta$group,
                     levels = c("Young_Pre", "Young_Post", "Old_Pre", "Old_Post"))

MEs        <- readRDS(file.path(DAT, "MEs.rds"))
kME_all    <- readRDS(file.path(DAT, "kME_all.rds"))
imp_mat    <- readRDS(file.path(DAT, "imp_mat.rds"))

module_df  <- read_csv(file.path(DAT, "wgcna/wgcna_module_assignments.csv"),
                       show_col_types = FALSE)
go_df      <- read_csv(file.path(DAT, "wgcna/wgcna_module_GO_enrichment.csv"),
                       show_col_types = FALSE)

# Trait correlations and BH p-values (data frames with 'module' column)
trait_cor  <- read_csv(file.path(DAT, "wgcna/wgcna_module_trait_correlations.csv"),
                       show_col_types = FALSE)
pval_bh    <- read_csv(file.path(DAT, "wgcna/wgcna_module_trait_pvalues_bh.csv"),
                       show_col_types = FALSE)

imp_annot  <- read_csv(file.path(DAT, "imp_annotations.csv"), show_col_types = FALSE)
imp_df     <- imp_annot  # alias used by kME matching below

mod_bio_labels_df  <- read_csv(file.path(DAT, "mod_bio_labels.csv"), show_col_types = FALSE)
mod_bio_labels_vec <- setNames(mod_bio_labels_df$bio_label, mod_bio_labels_df$module_color)
mod_display_label  <- function(color) {
  lbl <- mod_bio_labels_vec[color]
  lbl[is.na(lbl)] <- stringr::str_to_title(color[is.na(lbl)])
  paste0(lbl, " (", stringr::str_to_title(color), ")")
}

message("Panel C: per-module triptych...")

PC_W <- 580  # panel width mm
PC_H <- 340  # panel height mm

KEY_MODULES <- c("blue", "brown", "turquoise", "green", "black", "pink")
Z_CAP <- 2
X_CAP <- 12

GROUP_COLS <- c("Young_Pre", "Young_Post", "Old_Pre", "Old_Post")
GROUP_LABS <- c("Y-Pre", "Y-Post", "O-Pre", "O-Post")

MODULE_TRAIT <- c(
  blue = "age", brown = "training", turquoise = "age",
  green = "age", black = "age", pink = "age"
)

txt_accent <- scale_text(BASE_STAT, PC_W)
txt_annot  <- scale_text(BASE_GENE, PC_W)
txt_axis   <- scale_text(BASE_STAT, PC_W)
txt_tick   <- scale_text(BASE_GENE, PC_W)

group_samples <- list(
  Young_Pre  = meta$sample_id[meta$age == "Young" & meta$time == "Pre"],
  Young_Post = meta$sample_id[meta$age == "Young" & meta$time == "Post"],
  Old_Pre    = meta$sample_id[meta$age == "Old"   & meta$time == "Pre"],
  Old_Post   = meta$sample_id[meta$age == "Old"   & meta$time == "Post"]
)

trait_heatmap_df <- read_csv(file.path(DAT, "02_panel_B_heatmap_data.csv"),
                             show_col_types = FALSE)
KEY_TRAITS      <- c("age_num", "time_num", "VL_thick_cm", "Type_II_fCSA")
KEY_TRAIT_LABELS <- c("Age", "Training", "VL Thick.", "Type II fCSA")

module_sizes <- c()
heatmap_list   <- list()
enrich_list    <- list()
sankey_list    <- list()

for (mod in KEY_MODULES) {

  mod_proteins <- module_df %>% filter(module_color == mod)

  kme_col <- paste0("kME", mod)
  if (kme_col %in% colnames(kME_all)) {
    mod_proteins <- mod_proteins %>%
      left_join(
        tibble(uniprot_id = mod_proteins$uniprot_id) %>%
          left_join(
            tibble(uniprot_id = imp_df$uniprot_id,
                   kME = kME_all[match(imp_df$uniprot_id, rownames(imp_mat)), kme_col]),
            by = "uniprot_id"
          ),
        by = "uniprot_id"
      ) %>%
      arrange(desc(kME))
  }

  mod_ids   <- mod_proteins$uniprot_id
  mod_genes <- mod_proteins$gene
  n_mod     <- length(mod_ids)
  module_sizes[mod] <- n_mod

  mod_mat <- imp_mat[intersect(mod_ids, rownames(imp_mat)), , drop = FALSE]
  grp_means <- sapply(group_samples, function(samps) {
    rowMeans(mod_mat[, intersect(samps, colnames(mod_mat)), drop = FALSE], na.rm = TRUE)
  })
  grp_z <- t(scale(t(grp_means)))
  grp_z[grp_z >  Z_CAP] <-  Z_CAP
  grp_z[grp_z < -Z_CAP] <- -Z_CAP
  row_order <- intersect(mod_ids, rownames(grp_z))
  grp_z <- grp_z[row_order, , drop = FALSE]

  uid_to_gene <- setNames(mod_proteins$gene, mod_proteins$uniprot_id)

  heatmap_list[[mod]] <- as.data.frame(grp_z) %>%
    rownames_to_column("uniprot_id") %>%
    mutate(gene = uid_to_gene[uniprot_id]) %>%
    pivot_longer(cols = all_of(GROUP_COLS),
                 names_to = "group", values_to = "z") %>%
    mutate(module = mod)

  mod_enrich <- go_df %>%
    filter(module == mod) %>%
    mutate(database = paste0("GO:", ONTOLOGY)) %>%
    group_by(database) %>%
    slice_min(p.adjust, n = 3, with_ties = FALSE) %>%
    ungroup() %>%
    mutate(module = mod)

  enrich_list[[mod]] <- mod_enrich

  if (nrow(mod_enrich) > 0) {
    mod_enrich$gene_list <- strsplit(mod_enrich$geneID, "/")

    assigned <- tibble(gene = character(), pathway = character(),
                       database = character(), module = character())

    for (g in mod_genes) {
      best_idx <- NA
      for (j in seq_len(nrow(mod_enrich))) {
        if (g %in% mod_enrich$gene_list[[j]]) {
          best_idx <- j
          break
        }
      }
      if (!is.na(best_idx)) {
        assigned <- bind_rows(assigned, tibble(
          gene = g, pathway = mod_enrich$Description[best_idx],
          database = mod_enrich$database[best_idx], module = mod
        ))
      } else {
        assigned <- bind_rows(assigned, tibble(
          gene = g, pathway = "Unmapped",
          database = NA_character_, module = mod
        ))
      }
    }
    sankey_list[[mod]] <- assigned
  } else {
    sankey_list[[mod]] <- tibble(
      gene = mod_genes, pathway = "Unmapped",
      database = NA_character_, module = mod
    )
  }
}

me_long <- MEs %>%
  rownames_to_column("sample_id") %>%
  left_join(meta %>% dplyr::select(sample_id, group, subject, age, time),
            by = "sample_id")

eigen_data_list <- list()
for (mod in KEY_MODULES) {
  me_col <- paste0("ME", mod)
  if (!(me_col %in% colnames(MEs))) next
  eig <- me_long %>%
    dplyr::select(sample_id, group, subject, age, time, eigengene = all_of(me_col)) %>%
    mutate(module = mod)
  eigen_data_list[[mod]] <- eig
}
eigen_all <- bind_rows(eigen_data_list)

group_short <- c(Young_Pre = "Y-Pre", Young_Post = "Y-Post",
                 Old_Pre = "O-Pre", Old_Post = "O-Post")
eigen_all$group_label <- group_short[as.character(eigen_all$group)]
eigen_all$group_label <- factor(eigen_all$group_label,
                                 levels = c("Y-Pre", "Y-Post", "O-Pre", "O-Post"))

all_enrich <- bind_rows(enrich_list)
if (nrow(all_enrich) > 0) {
  max_neg_log10 <- max(-log10(all_enrich$p.adjust), na.rm = TRUE)
  shared_x_max <- min(X_CAP, max_neg_log10) + 2
} else {
  shared_x_max <- X_CAP + 2
}

build_triptych_row <- function(mod, show_xlab, shared_x_max) {

  mod_col    <- mod
  border_col <- darken_color(mod_col, 0.3)

  mod_proteins <- module_df %>% filter(module_color == mod)
  mod_ids <- mod_proteins$uniprot_id
  kme_col <- paste0("kME", mod)
  if (kme_col %in% colnames(kME_all)) {
    kme_vals <- kME_all[match(mod_ids, rownames(imp_mat)), kme_col]
    mod_ids <- mod_ids[order(-kme_vals)]
  }
  mod_genes <- mod_proteins$gene[match(mod_ids, mod_proteins$uniprot_id)]
  n_genes <- length(mod_ids)

  p_accent <- ggplot() +
    geom_rect(aes(xmin = 0, xmax = 1, ymin = 0, ymax = 1), fill = mod_col) +
    annotate("text", x = 0.5, y = 0.5,
             label = paste0(str_to_title(mod), "\nn=", n_genes),
             size = txt_accent, fontface = "bold",
             color = if (mod %in% c("yellow", "pink")) "grey20" else "white",
             lineheight = 0.85, angle = 90) +
    coord_cartesian(expand = FALSE) +
    theme_void() +
    theme(plot.margin = margin(0, 0, 0, 0))

  ht_long <- heatmap_list[[mod]] %>%
    mutate(
      uniprot_id = factor(uniprot_id, levels = rev(mod_ids)),
      group = factor(group, levels = GROUP_COLS)
    )

  ht_colors <- c("#2166AC", "#92C5DE", "white", "#F4A582", "#B2182B")

  p_ht <- ggplot(ht_long, aes(x = group, y = uniprot_id, fill = z)) +
    geom_raster() +
    scale_fill_gradientn(
      colours = ht_colors, limits = c(-Z_CAP, Z_CAP),
      oob = scales::squish, guide = "none"
    ) +
    scale_x_discrete(
      labels = if (show_xlab) GROUP_LABS else NULL,
      expand = expansion(0)
    ) +
    scale_y_discrete(expand = expansion(0)) +
    labs(x = NULL, y = NULL) +
    FIG_THEME +
    theme(
      axis.text.y  = element_blank(),
      axis.ticks.y = element_blank(),
      axis.text.x  = if (show_xlab) element_text(size = txt_tick, angle = 45,
                                                   hjust = 1, face = "bold")
                      else element_blank(),
      axis.ticks.x = if (show_xlab) element_line() else element_blank(),
      panel.border = element_rect(colour = border_col, linewidth = 0.8, fill = NA),
      legend.position = "none",
      plot.margin  = margin(t = 1, r = -5, b = if (show_xlab) 4 else 0, l = 0)
    )

  cl_links  <- sankey_list[[mod]]
  cl_enrich <- enrich_list[[mod]]
  mapped_links <- cl_links %>% filter(pathway != "Unmapped")
  mapped_pathways <- unique(mapped_links$pathway)
  n_pw <- length(mapped_pathways)

  if (n_pw == 0) {
    p_sankey <- ggplot() +
      annotate("text", x = 0.5, y = 0.5, label = "No mapped\npathways",
               size = txt_annot, color = "grey50") +
      theme_void() +
      theme(plot.margin = margin(t = 1, r = 0, b = if (show_xlab) 4 else 0, l = 0))
  } else {
    Y_SPAN <- n_genes
    X_GENE <- 0.0
    X_PW   <- 1.5
    BAR_W  <- 0.06

    gene_h   <- Y_SPAN / (n_genes * 1.15)
    gene_gap <- if (n_genes > 1) (Y_SPAN - n_genes * gene_h) / (n_genes - 1) else 0

    gene_bars <- tibble(
      gene = mod_genes, idx = seq_along(mod_genes),
      y_top = Y_SPAN - (idx - 1) * (gene_h + gene_gap),
      y_bot = y_top - gene_h,
      y_ctr = (y_top + y_bot) / 2,
      x_left = X_GENE - BAR_W / 2, x_right = X_GENE + BAR_W / 2,
      fill = mod_col
    )

    pw_h   <- Y_SPAN / (n_pw * 1.4)
    pw_gap <- if (n_pw > 1) (Y_SPAN - n_pw * pw_h) / (n_pw - 1) else 0

    pw_order <- cl_enrich %>%
      filter(Description %in% mapped_pathways) %>%
      arrange(p.adjust) %>%
      pull(Description)
    pw_order <- c(pw_order, setdiff(mapped_pathways, pw_order))

    pw_bars <- tibble(
      pathway = pw_order, idx = seq_along(pw_order),
      y_top = Y_SPAN - (idx - 1) * (pw_h + pw_gap),
      y_bot = y_top - pw_h,
      y_ctr = (y_top + y_bot) / 2,
      x_left = X_PW - BAR_W / 2, x_right = X_PW + BAR_W / 2
    ) %>%
      left_join(mapped_links %>% dplyr::select(pathway, database) %>% distinct(),
                by = "pathway") %>%
      mutate(fill = DB_COLORS[database])

    slot_data <- mapped_links %>%
      group_by(pathway) %>%
      mutate(slot_idx = row_number(), n_slots = n()) %>%
      ungroup() %>%
      left_join(pw_bars %>% dplyr::select(pathway, y_top, y_bot), by = "pathway") %>%
      mutate(
        slot_h   = (y_top - y_bot) / n_slots,
        slot_top = y_top - (slot_idx - 1) * slot_h,
        slot_bot = slot_top - slot_h,
        slot_ctr = (slot_top + slot_bot) / 2
      )

    ribbon_input <- slot_data %>%
      left_join(gene_bars %>% dplyr::select(gene, y_top_gene = y_top,
                                             y_bot_gene = y_bot, y_ctr_gene = y_ctr),
                by = "gene") %>%
      filter(!is.na(y_top_gene)) %>%
      mutate(
        ribbon_id  = paste0(mod, "_", gene, "_", pathway),
        cross_dist = abs(y_ctr_gene - slot_ctr)
      ) %>%
      arrange(desc(cross_dist))

    ribbon_polys <- pmap_dfr(ribbon_input, function(...) {
      args <- list(...)
      make_sigmoid_ribbon(
        x0 = X_GENE + BAR_W / 2, x1 = X_PW - BAR_W / 2,
        y0_top = args$y_top_gene, y0_bot = args$y_bot_gene,
        y1_top = args$slot_top, y1_bot = args$slot_bot,
        ribbon_id = args$ribbon_id
      ) %>% mutate(fill = DB_COLORS[args$database])
    })

    gene_rect <- gene_bars %>% dplyr::select(x_left, x_right, y_bot, y_top, fill)
    pw_rect   <- pw_bars %>% dplyr::select(x_left, x_right, y_bot, y_top, fill)

    pw_label_size <- if (n_pw <= 4) txt_annot * 1.1 else if (n_pw <= 7) txt_annot else txt_annot * 0.85

    p_sankey <- ggplot() +
      geom_polygon(data = ribbon_polys,
                   aes(x = x, y = y, group = ribbon_id, fill = fill),
                   alpha = 0.25, color = NA) +
      geom_rect(data = gene_rect,
                aes(xmin = x_left, xmax = x_right, ymin = y_bot, ymax = y_top, fill = fill),
                color = NA) +
      geom_rect(data = pw_rect,
                aes(xmin = x_left, xmax = x_right, ymin = y_bot, ymax = y_top, fill = fill),
                color = NA) +
      geom_text(data = pw_bars,
                aes(x = x_left - 0.04, y = y_ctr,
                    label = str_trunc(clean_pathway_name(pathway), 32, ellipsis = "...")),
                hjust = 1, size = pw_label_size, fontface = "bold", color = "grey20") +
      scale_fill_identity() +
      coord_cartesian(xlim = c(-0.04, 1.6), ylim = c(0, Y_SPAN), expand = FALSE) +
      theme_void() +
      theme(plot.margin = margin(t = 1, r = 0, b = if (show_xlab) 4 else 0, l = 0))
  }

  cl_enrich_vis <- cl_enrich %>% filter(Description %in% mapped_pathways)

  if (nrow(cl_enrich_vis) == 0) {
    p_bars <- ggplot() +
      annotate("text", x = 0.5, y = 0.5, label = "No significant\nenrichment",
               size = txt_annot, color = "grey50") +
      theme_void() +
      theme(plot.margin = margin(t = 1, r = 0, b = if (show_xlab) 4 else 0, l = 0))
  } else {

    bar_data <- cl_enrich_vis %>%
      mutate(
        neg_log10_p = -log10(p.adjust),
        clean_name  = clean_pathway_name(Description),
        gene_count  = Count,
        db_fill     = DB_COLORS[database],
        is_capped   = neg_log10_p > X_CAP,
        display_val = pmin(neg_log10_p, X_CAP)
      ) %>%
      left_join(pw_bars %>% dplyr::select(pathway, y_top, y_bot, y_ctr),
                by = c("Description" = "pathway"))

    bar_label_size <- if (nrow(bar_data) <= 4) txt_annot else if (nrow(bar_data) <= 7) txt_annot * 0.9 else txt_annot * 0.8

    p_bars <- ggplot(bar_data) +
      geom_rect(aes(xmin = 0, xmax = display_val, ymin = y_bot, ymax = y_top,
                     fill = db_fill),
                color = "black", linewidth = 0.25) +
      geom_text(aes(x = display_val, y = y_ctr, label = gene_count),
                hjust = -0.3, size = bar_label_size * 0.7, fontface = "bold") +
      geom_text(
        data = bar_data %>% filter(is_capped),
        aes(x = X_CAP - 0.2, y = y_ctr, label = "//"),
        hjust = 0.5, size = txt_annot * 0.8, fontface = "bold", color = "white"
      ) +
      geom_text(
        data = bar_data %>% filter(is_capped),
        aes(x = X_CAP + 0.2, y = y_ctr,
            label = sprintf("%.0f", neg_log10_p)),
        hjust = 0, size = txt_annot * 0.7, fontface = "italic", color = "grey30"
      ) +
      geom_vline(xintercept = -log10(0.05), linetype = "dashed",
                 color = "grey40", linewidth = 0.3) +
      scale_fill_identity() +
      scale_x_continuous(
        expand = expansion(mult = c(0, 0.05)),
        name   = if (show_xlab) expression(-log[10](p[adj])) else NULL
      ) +
      coord_cartesian(xlim = c(0, shared_x_max), ylim = c(0, Y_SPAN),
                      expand = FALSE) +
      labs(y = NULL) +
      FIG_THEME +
      theme(
        axis.text.y  = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.x  = if (show_xlab) element_text(size = txt_tick, face = "bold") else element_blank(),
        axis.ticks.x = if (show_xlab) element_line() else element_blank(),
        axis.title.x = if (show_xlab) element_text(size = txt_axis, face = "bold") else element_blank(),
        panel.border = element_rect(colour = "grey70", linewidth = 0.3, fill = NA),
        legend.position = "none",
        plot.margin  = margin(t = 1, r = 1, b = if (show_xlab) 4 else 0, l = 0)
      )
  }

  trait_mod <- trait_heatmap_df %>%
    filter(module == paste0("ME", mod), trait %in% KEY_TRAITS) %>%
    mutate(trait_label = factor(trait_label, levels = rev(KEY_TRAIT_LABELS)))

  p_trait <- ggplot(trait_mod, aes(x = 1, y = trait_label, fill = cor)) +
    geom_tile(color = "grey60", linewidth = 0.3) +
    geom_text(aes(label = stars), size = txt_annot, fontface = "bold", vjust = 0.75) +
    scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B",
                         midpoint = 0, limits = c(-1, 1), guide = "none") +
    scale_x_continuous(expand = expansion(0)) +
    labs(x = if (show_xlab) "Trait r" else NULL, y = NULL) +
    FIG_THEME +
    theme(
      axis.text.y    = element_text(size = txt_tick, face = "bold"),
      axis.text.x    = element_blank(),
      axis.ticks     = element_blank(),
      axis.title.x   = if (show_xlab) element_text(size = txt_tick, face = "bold.italic")
                        else element_blank(),
      panel.border   = element_rect(colour = "grey70", linewidth = 0.3, fill = NA),
      legend.position = "none",
      plot.margin    = margin(t = 1, r = 1, b = if (show_xlab) 4 else 0, l = 1)
    )

  eig_mod <- eigen_all %>% filter(module == mod)
  trait_type  <- MODULE_TRAIT[mod]
  trait_label <- if (trait_type == "age") "Age-Correlated Module" else "Training-Responsive Module"

  eig_slope <- eig_mod %>%
    mutate(age_fct  = factor(age, levels = c("Young", "Old")),
           time_fct = factor(time, levels = c("Pre", "Post")))

  eig_means <- eig_slope %>%
    group_by(age_fct, time_fct) %>%
    summarise(eigengene = mean(eigengene), .groups = "drop")

  p_eigen <- ggplot(eig_slope, aes(x = time_fct, y = eigengene)) +
    geom_line(aes(group = subject, color = age_fct),
              alpha = 0.25, linewidth = 0.3) +
    geom_point(aes(fill = age_fct), shape = 21, size = 1.0,
               alpha = 0.45, stroke = 0.2, color = "white") +
    geom_line(data = eig_means, aes(group = 1, color = age_fct),
              linewidth = 1.3) +
    geom_point(data = eig_means, aes(fill = age_fct), shape = 21,
               size = 2.8, stroke = 0.6, color = "white") +
    facet_wrap(~age_fct, nrow = 1) +
    scale_color_manual(values = AGE_COLORS, guide = "none") +
    scale_fill_manual(values = AGE_COLORS, guide = "none") +
    labs(x = NULL,
         y = if (show_xlab) "Eigengene" else NULL,
         title = trait_label) +
    FIG_THEME +
    theme(
      plot.title       = element_text(size = txt_axis, face = "bold.italic",
                                      color = "grey25", hjust = 0.5,
                                      margin = margin(0, 0, 2, 0)),
      strip.text       = element_text(size = txt_tick, face = "bold", color = "grey15"),
      strip.background = element_rect(fill = "grey93", colour = NA),
      axis.text.x  = if (show_xlab) element_text(size = txt_tick, face = "bold")
                      else element_blank(),
      axis.ticks.x = if (show_xlab) element_line() else element_blank(),
      axis.text.y  = element_text(size = txt_tick, face = "bold"),
      axis.title.y = if (show_xlab) element_text(size = txt_axis, face = "bold") else element_blank(),
      panel.border = element_rect(colour = border_col, linewidth = 0.5, fill = NA),
      panel.spacing = unit(3, "pt"),
      legend.position = "none",
      plot.margin  = margin(t = 1, r = 1, b = if (show_xlab) 4 else 0, l = 1)
    )

  (p_ht | p_sankey | p_bars | p_trait | p_eigen) +
    plot_layout(widths = c(0.08, 0.30, 0.26, 0.08, 0.18), nrow = 1)
}

panels <- lapply(seq_along(KEY_MODULES), function(i) {
  build_triptych_row(KEY_MODULES[i],
                     show_xlab = (i == length(KEY_MODULES)),
                     shared_x_max = shared_x_max)
})

raw_sizes <- module_sizes[KEY_MODULES]
row_heights <- sqrt(raw_sizes) / sum(sqrt(raw_sizes))
row_heights <- pmax(row_heights, 0.10)
row_heights <- row_heights / sum(row_heights)

panel_C <- wrap_plots(panels, ncol = 1, heights = row_heights) +
  plot_annotation(
    title    = "Module Functional Characterization",
    subtitle = "Z-scored group means  |  Gene-pathway Sankey  |  GO enrichment  |  Trait r  |  Eigengene slopes",
    theme = theme(
      plot.title    = element_text(face = "bold", size = txt_axis * 3),
      plot.subtitle = element_text(size = txt_axis * 2, color = "grey25",
                                   face = "bold.italic")
    )
  )

ggsave(file.path(RPT, "panel_C_triptych.pdf"), panel_C,
       width = PC_W, height = PC_H, units = "mm",
       device = pdf_device, limitsize = FALSE)
ggsave(file.path(RPT, "panel_C_triptych.png"), panel_C,
       width = PC_W, height = PC_H, units = "mm",
       dpi = 300, limitsize = FALSE)

write_csv(all_enrich, file.path(DAT, "03_panel_C_triptych_enrichment.csv"))
write_csv(bind_rows(sankey_list), file.path(DAT, "03_panel_C_sankey_links.csv"))
write_csv(eigen_all, file.path(DAT, "03_panel_C_eigengene_data.csv"))
write_csv(bind_rows(heatmap_list), file.path(DAT, "03_panel_C_heatmap_zscores.csv"))

# Module eigengene reversal (paired t-test + Shapiro-Wilk + Wilcoxon sensitivity)
age_modules <- pval_bh %>%
  dplyr::select(module, age_num) %>%
  filter(age_num < 0.05) %>%
  mutate(module_color = gsub("^ME", "", module)) %>%
  pull(module_color)

age_cors <- trait_cor %>%
  dplyr::select(module, cor = age_num) %>%
  mutate(module_color = gsub("^ME", "", module))

pre_samples  <- meta$sample_id[meta$time == "Pre"]
post_samples <- meta$sample_id[meta$time == "Post"]

me_pre_f5  <- MEs[pre_samples, , drop = FALSE]
me_post_f5 <- MEs[post_samples, , drop = FALSE]

pre_subjects  <- sub("_(Pre|Post)$", "", pre_samples)
post_subjects <- sub("_(Pre|Post)$", "", post_samples)
rownames(me_pre_f5)  <- pre_subjects
rownames(me_post_f5) <- post_subjects

old_subjects <- unique(meta$subject[meta$age == "Old"])
old_subjects <- intersect(old_subjects,
                           intersect(rownames(me_pre_f5), rownames(me_post_f5)))

young_pre_samps <- meta$sample_id[meta$age == "Young" & meta$time == "Pre"]

reversal_results <- list()

for (mod in age_modules) {
  me_col <- paste0("ME", mod)
  if (!(me_col %in% colnames(MEs))) next

  old_pre_me  <- me_pre_f5[old_subjects, me_col]
  old_post_me <- me_post_f5[old_subjects, me_col]

  t_res <- t.test(old_post_me, old_pre_me, paired = TRUE)

  diffs <- old_post_me - old_pre_me
  sw_test <- shapiro.test(diffs)
  wilcox_sensitivity <- if (sw_test$p.value < 0.05) {
    wilcox.test(old_post_me, old_pre_me, paired = TRUE, conf.int = TRUE)
  } else {
    NULL
  }

  young_pre_mean <- mean(MEs[young_pre_samps, me_col], na.rm = TRUE)
  old_pre_mean_val <- mean(old_pre_me)
  old_post_mean_val <- mean(old_post_me)

  toward_young <- sign(old_post_mean_val - old_pre_mean_val) ==
                  sign(young_pre_mean - old_pre_mean_val)

  age_gap <- abs(young_pre_mean - old_pre_mean_val)
  shift <- old_post_mean_val - old_pre_mean_val
  reversal_pct <- if (age_gap > 0) (shift / (young_pre_mean - old_pre_mean_val)) * 100 else NA

  age_r <- age_cors$cor[age_cors$module_color == mod]

  reversal_results <- c(reversal_results, list(tibble(
    module = mod,
    age_cor = if (length(age_r) > 0) round(age_r, 3) else NA_real_,
    mean_old_pre  = round(old_pre_mean_val, 4),
    mean_old_post = round(old_post_mean_val, 4),
    mean_young_pre = round(young_pre_mean, 4),
    t_stat = round(t_res$statistic, 3),
    p_value = t_res$p.value,
    ci_lo_mean_diff = round(t_res$conf.int[1], 4),
    ci_hi_mean_diff = round(t_res$conf.int[2], 4),
    shapiro_p = round(sw_test$p.value, 4),
    normality_ok = sw_test$p.value >= 0.05,
    wilcox_p_sensitivity = if (!is.null(wilcox_sensitivity)) wilcox_sensitivity$p.value else NA_real_,
    toward_young = toward_young,
    reversal_pct = round(reversal_pct, 1),
    n_old_subjects = length(old_subjects)
  )))
}

if (length(reversal_results) > 0) {
  rev_df <- bind_rows(reversal_results)
  rev_df$p_adj_bh <- p.adjust(rev_df$p_value, method = "BH")
  write_csv(rev_df, file.path(DAT, "07_module_reversal.csv"))
}

# Age gap closure (Wilcoxon signed-rank + Hodges-Lehmann CI)
all_modules <- unique(module_df$module_color)
all_modules <- all_modules[all_modules != "grey"]

young_pre_samps  <- meta$sample_id[meta$age == "Young" & meta$time == "Pre"]
young_post_samps <- meta$sample_id[meta$age == "Young" & meta$time == "Post"]
old_pre_samps    <- meta$sample_id[meta$age == "Old"   & meta$time == "Pre"]
old_post_samps   <- meta$sample_id[meta$age == "Old"   & meta$time == "Post"]

closure_results <- list()
for (mod in all_modules) {
  me_col <- paste0("ME", mod)
  if (!(me_col %in% colnames(MEs))) next

  gap_pre  <- abs(mean(MEs[young_pre_samps, me_col], na.rm = TRUE) -
                  mean(MEs[old_pre_samps, me_col], na.rm = TRUE))
  gap_post <- abs(mean(MEs[young_post_samps, me_col], na.rm = TRUE) -
                  mean(MEs[old_post_samps, me_col], na.rm = TRUE))
  closure  <- (gap_pre - gap_post) / gap_pre * 100

  closure_results <- c(closure_results, list(tibble(
    module = mod, gap_pre = round(gap_pre, 4), gap_post = round(gap_post, 4),
    closure_pct = round(closure, 1)
  )))
}

closure_df <- bind_rows(closure_results)

wilcox_res <- wilcox.test(closure_df$gap_pre, closure_df$gap_post,
                           paired = TRUE, conf.int = TRUE)

closure_df$gap_closure_wilcox_p <- wilcox_res$p.value
closure_df$gap_closure_pseudomedian <- wilcox_res$estimate
closure_df$gap_closure_ci_lo <- wilcox_res$conf.int[1]
closure_df$gap_closure_ci_hi <- wilcox_res$conf.int[2]
write_csv(closure_df, file.path(DAT, "08_age_gap_closure.csv"))

closure_long <- closure_df %>%
  pivot_longer(cols = c(gap_pre, gap_post), names_to = "timepoint",
               values_to = "gap") %>%
  mutate(timepoint = factor(ifelse(timepoint == "gap_pre", "Pre", "Post"),
                             levels = c("Pre", "Post")))

p_closure <- ggplot(closure_long, aes(x = timepoint, y = gap)) +
  geom_line(aes(group = module, color = module), alpha = 0.6, linewidth = 0.5) +
  geom_point(aes(fill = module), shape = 21, size = 2.5, stroke = 0.3,
             color = "white") +
  scale_color_identity() + scale_fill_identity() +
  annotate("text", x = 1.5, y = max(closure_long$gap) * 1.05,
           label = sprintf("Wilcoxon p = %.3g", wilcox_res$p.value),
           size = txt_annot, fontface = "bold", color = "grey25") +
  labs(title = "Age Gap Closure by Module",
       subtitle = "|Young_mean - Old_mean| eigengene: Pre vs Post training",
       x = NULL, y = "Age Gap (|ME_Young - ME_Old|)") +
  FIG_THEME +
  theme(legend.position = "none")

ggsave(file.path(RPT, "age_gap_closure.pdf"), p_closure,
       width = 140, height = 160, units = "mm", device = pdf_device)

message("  Panel C complete")
