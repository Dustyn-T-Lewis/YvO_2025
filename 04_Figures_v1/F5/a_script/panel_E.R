################################################################################
#   YvO Figure 5 — Panel E: Functionally Annotated Hub Protein Networks
#   6 key modules in a 2x3 grid, with GO-based hull annotations
#   + trait correlation sidebar heatmaps per module
################################################################################
#
# ── STAT AUDIT (2026-02-27) ──────────────────────────────────────────────────
#
# TOM THRESHOLD (75th percentile)
#   Edges below the 75th percentile of TOM values within the top-30 hub
#   sub-network are zeroed.  This is a visualization filter, not a
#   statistical threshold.  SENSITIVITY: at 50th percentile the graph is
#   dense and unreadable; at 90th it is over-sparse.  The 75th percentile
#   is a common heuristic (Langfelder & Horvath tutorials).
#   No formal sensitivity sweep is needed for a display-only threshold.
#
# kME (module membership)
#   Pearson r between each protein's expression and the module eigengene.
#   CIs are not displayed on the network nodes (impractical for n=30
#   nodes per module).  AUDIT ADDITION: we compute cor.test() CIs for
#   all hub proteins and export as supplementary CSV.
#
# GS (Gene Significance)
#   Pearson r between protein expression and the module's key phenotype.
#   Same CI note as kME above — exported in the supplementary CSV.
#
# FUNCTIONAL GROUP ASSIGNMENT
#   Greedy 1:1 gene-to-GO-term mapping (best p.adjust per gene, then
#   only terms with >= 3 genes form a hull).  This is descriptive
#   grouping for visualization, not a statistical test.  No correction
#   is applied or needed.
#
# AUDIT ADDITIONS:
#   1. cor.test() CIs for kME and GS of all hub proteins → CSV export
# ─────────────────────────────────────────────────────────────────────────────

source("04_Figures/F5/a_script/YvO_F5_setup.R")

suppressPackageStartupMessages({
  library(ggforce)
  library(concaveman)
  library(graphlayouts)
  library(tidygraph)
  library(shadowtext)
  library(ggnewscale)
})

# === Constants ================================================================
KEY_MODULES <- c("blue", "brown", "turquoise", "green", "black", "pink")
TOP_N       <- 30
EDGE_QUANT  <- 0.75
NET_POWER   <- 14
MAX_GROUPS  <- 4
MIN_GROUP_N <- 3

# Context ring (annulus of all remaining module proteins)
PERIPH_R_INNER_MULT <- 1.30   # inner annulus edge (× hub max radius)
PERIPH_R_OUTER_MULT <- 2.05   # outer annulus edge
PERIPH_ALPHA        <- 0.38
PERIPH_STROKE       <- 0.15

MODULE_GS_PHENO <- c(
  blue = "Type_II_fCSA", brown = "training", turquoise = "Type_II_fCSA",
  green = "VL_thick_cm", black = "age", pink = "age"
)
MODULE_GS_LABEL <- c(
  blue = "Type II fCSA", brown = "Training", turquoise = "Type II fCSA",
  green = "VL Thickness", black = "Age", pink = "Age"
)

HULL_PALETTE <- c(
  "#A6CEE3", "#B2DF8A", "#FDBF6F", "#CAB2D6",
  "#FB9A99", "#D9D9D9", "#FFFF99", "#B3B3B3"
)

# Human-readable trait labels for sidebar
TRAIT_LABELS <- c(
  age_num         = "Age",
  time_num        = "Training",
  interaction     = "Age\u00d7Train",
  VL_thick_cm     = "VL Thick.",
  DXA_LBM_kg      = "LBM",
  BMI             = "BMI",
  Type_I_fCSA     = "Type I",
  Type_II_fCSA    = "Type II",
  deadlift_1rm_kg = "DL 1RM"
)

# Trait display order (design factors first, then phenotypes)
TRAIT_ORDER <- c("age_num", "time_num", "interaction",
                 "VL_thick_cm", "Type_I_fCSA", "Type_II_fCSA",
                 "DXA_LBM_kg", "BMI", "deadlift_1rm_kg")

# === Prepare phenotype vectors ================================================
meta$age_binary  <- ifelse(meta$age == "Old", 1, 0)
meta$time_binary <- ifelse(meta$time == "Post", 1, 0)

# === Load module-trait correlations ===========================================
trait_cor_df <- read_csv("04_Figures/F5/c_data/wgcna/wgcna_module_trait_correlations.csv",
                         show_col_types = FALSE)
trait_pval_df <- read_csv("04_Figures/F5/c_data/wgcna/wgcna_module_trait_pvalues_bh.csv",
                          show_col_types = FALSE)

# === Load GO enrichment data ==================================================
go_df_full <- read_csv("04_Figures/F5/c_data/wgcna/wgcna_module_GO_enrichment.csv",
                        show_col_types = FALSE)

# === Trait sidebar heatmap builder ============================================
make_trait_sidebar <- function(mod) {
  me_name <- paste0("ME", mod)
  traits <- intersect(TRAIT_ORDER, names(trait_cor_df))

  cor_row <- trait_cor_df %>% filter(module == me_name)
  pval_row <- trait_pval_df %>% filter(module == me_name)

  if (nrow(cor_row) == 0) return(ggplot() + theme_void())

  sidebar_df <- tibble(trait = traits) %>%
    mutate(
      r     = as.numeric(cor_row[1, traits]),
      p_adj = as.numeric(pval_row[1, traits]),
      star  = case_when(
        p_adj < 0.001 ~ "***",
        p_adj < 0.01  ~ "**",
        p_adj < 0.05  ~ "*",
        TRUE          ~ ""
      ),
      label = sprintf("%+.2f%s", r, star),
      trait_label = TRAIT_LABELS[trait],
      trait_label = factor(trait_label, levels = rev(TRAIT_LABELS[traits]))
    )

  ggplot(sidebar_df, aes(x = 1, y = trait_label, fill = r)) +
    geom_tile(color = "white", linewidth = 0.4) +
    geom_text(aes(label = label), size = 1.9, color = "black", fontface = "bold") +
    scale_fill_gradient2(
      low = "#2166AC", mid = "white", high = "#B2182B",
      midpoint = 0, limits = c(-1, 1), guide = "none"
    ) +
    scale_x_continuous(expand = c(0, 0)) +
    labs(x = NULL, y = NULL) +
    theme_void() +
    theme(
      axis.text.y  = element_text(size = 5, hjust = 1, margin = margin(r = 1)),
      axis.text.x  = element_blank(),
      plot.margin  = margin(t = 18, b = 2, l = 0, r = 2)
    )
}

# === Functional group assignment ==============================================
assign_functional_groups <- function(mod, gene_names, max_groups = MAX_GROUPS) {
  mod_go <- go_df_full %>%
    filter(module == mod, ONTOLOGY == "BP") %>%
    arrange(p.adjust)

  if (nrow(mod_go) == 0) {
    mod_go <- go_df_full %>%
      filter(module == mod) %>%
      arrange(p.adjust)
  }
  if (nrow(mod_go) == 0) {
    return(setNames(rep("Other", length(gene_names)), gene_names))
  }

  gene_term_map <- mod_go %>%
    mutate(genes = strsplit(geneID, "/")) %>%
    unnest(genes) %>%
    filter(genes %in% gene_names) %>%
    group_by(genes) %>%
    slice_min(p.adjust, n = 1, with_ties = FALSE) %>%
    ungroup() %>%
    dplyr::select(gene = genes, go_term = Description, p_adj = p.adjust)

  term_counts <- gene_term_map %>% count(go_term, sort = TRUE)
  keep_terms <- term_counts %>%
    filter(n >= MIN_GROUP_N) %>%
    slice_head(n = max_groups) %>%
    pull(go_term)

  clean_term <- function(term) {
    term %>%
      str_replace_all("_", " ") %>%
      str_to_sentence() %>%
      str_replace("^Regulation of ", "Reg. ") %>%
      str_replace("^Positive regulation of ", "+Reg. ") %>%
      str_replace("^Negative regulation of ", "-Reg. ") %>%
      str_replace(" organization$", " org.") %>%
      str_replace(" organization ", " org. ") %>%
      str_replace(" biosynthetic process$", " biosynthesis") %>%
      str_replace(" catabolic process$", " catabolism") %>%
      str_replace(" metabolic process$", " metabolism") %>%
      str_replace("-based process$", " process") %>%
      str_replace("protein-containing complex", "complex") %>%
      str_replace(" development$", " dev.") %>%
      str_replace("Rna-templated dna", "RNA-templated DNA") %>%
      str_replace("Trna ", "tRNA ") %>%
      str_replace(" aminoacylation for protein translation", " charging") %>%
      str_trunc(42, ellipsis = "\u2026")
  }

  assignments <- setNames(rep("Other", length(gene_names)), gene_names)
  for (g in gene_names) {
    row <- gene_term_map %>% filter(gene == g)
    if (nrow(row) > 0 && row$go_term[1] %in% keep_terms) {
      assignments[g] <- clean_term(row$go_term[1])
    }
  }
  assignments
}

# === Per-module network builder ===============================================
build_hub_network <- function(mod, top_n = TOP_N, edge_quant = EDGE_QUANT) {

  mod_proteins <- module_df$uniprot_id[module_df$module_color == mod]
  mod_proteins <- intersect(mod_proteins, colnames(datExpr))
  if (length(mod_proteins) < 5) {
    warning(sprintf("Module '%s' has too few proteins (%d), skipping", mod, length(mod_proteins)))
    return(list(plot = NULL, node_data = NULL, periph_data = NULL))
  }

  mod_idx <- which(colnames(datExpr) %in% mod_proteins)
  adj_mod <- adjacency(datExpr[, mod_idx], power = NET_POWER, type = "signed")
  TOM_mod <- TOMsimilarity(adj_mod)
  colnames(TOM_mod) <- rownames(TOM_mod) <- colnames(datExpr)[mod_idx]

  kme_col <- paste0("kME", mod)
  if (!kme_col %in% colnames(kME_all)) {
    warning(sprintf("kME column '%s' not found", kme_col))
    return(list(plot = NULL, node_data = NULL, periph_data = NULL))
  }
  mod_kme <- setNames(kME_all[mod_proteins, kme_col], mod_proteins)
  mod_kme <- mod_kme[!is.na(mod_kme)]
  actual_n <- min(top_n, length(mod_kme))
  top_prots <- names(sort(mod_kme, decreasing = TRUE))[1:actual_n]

  TOM_sub <- TOM_mod[top_prots, top_prots]
  tom_vals <- TOM_sub[upper.tri(TOM_sub)]
  edge_thresh <- quantile(tom_vals, edge_quant)
  TOM_thresh <- TOM_sub
  TOM_thresh[TOM_thresh < edge_thresh] <- 0
  diag(TOM_thresh) <- 0

  g <- graph_from_adjacency_matrix(TOM_thresh, mode = "undirected", weighted = TRUE)

  pheno_key <- MODULE_GS_PHENO[[mod]]
  if (pheno_key == "age") {
    pheno_vec <- meta$age_binary[match(rownames(datExpr), meta$sample_id)]
  } else if (pheno_key == "training") {
    pheno_vec <- meta$time_binary[match(rownames(datExpr), meta$sample_id)]
  } else {
    pheno_vec <- meta[[pheno_key]][match(rownames(datExpr), meta$sample_id)]
  }

  # GS for ALL module proteins (hubs + periphery)
  all_mod_prots <- names(mod_kme)
  gs_all <- cor(datExpr[, all_mod_prots, drop = FALSE], pheno_vec,
                use = "pairwise.complete.obs")[, 1]
  gs_vals <- gs_all[top_prots]

  node_df <- tibble(
    uniprot_id = V(g)$name, module = mod,
    kME = mod_kme[V(g)$name], GS = gs_vals[V(g)$name],
    GS_pheno = MODULE_GS_LABEL[[mod]]
  ) %>%
    left_join(module_df %>% dplyr::select(uniprot_id, gene), by = "uniprot_id")

  func_groups <- assign_functional_groups(mod, node_df$gene)
  node_df$functional_group <- func_groups[node_df$gene]

  V(g)$kME  <- node_df$kME
  V(g)$GS   <- node_df$GS
  V(g)$gene <- node_df$gene
  V(g)$functional_group <- node_df$functional_group
  V(g)$label_show <- ifelse(
    rank(-node_df$kME) <= 8 | rank(-abs(node_df$GS)) <= 4,
    node_df$gene, NA
  )

  n_mod <- length(mod_proteins)

  set.seed(42)
  g_tidy <- as_tbl_graph(g)
  lay <- create_layout(g_tidy, layout = "fr")

  lay$functional_group <- V(g)$functional_group[match(lay$name, V(g)$name)]
  lay$kME        <- V(g)$kME[match(lay$name, V(g)$name)]
  lay$GS         <- V(g)$GS[match(lay$name, V(g)$name)]
  lay$label_show <- V(g)$label_show[match(lay$name, V(g)$name)]

  # --- Periphery ring: annulus of non-hub module proteins ---
  periph_prots <- setdiff(all_mod_prots, top_prots)
  periph_df <- NULL
  if (length(periph_prots) > 0) {
    hub_cx <- mean(lay$x)
    hub_cy <- mean(lay$y)
    hub_max_r <- max(sqrt((lay$x - hub_cx)^2 + (lay$y - hub_cy)^2))
    if (hub_max_r < 1e-6) hub_max_r <- 1  # safety for degenerate layouts

    r_inner <- hub_max_r * PERIPH_R_INNER_MULT
    r_outer <- hub_max_r * PERIPH_R_OUTER_MULT

    periph_kme <- mod_kme[periph_prots]
    periph_gs  <- gs_all[periph_prots]
    n_periph   <- length(periph_prots)

    # Angular positions: quasi-uniform with ±40% jitter
    base_angles <- seq(0, 2 * pi, length.out = n_periph + 1)[1:n_periph]
    set.seed(42 + n_periph)
    angle_jitter <- runif(n_periph, -0.4, 0.4) * (2 * pi / n_periph)
    angles <- base_angles + angle_jitter

    # Radial positions: higher kME → closer to inner edge, with ±8% jitter
    kme_range_p <- range(periph_kme, na.rm = TRUE)
    if (diff(kme_range_p) < 1e-8) {
      kme_frac <- rep(0.5, n_periph)
    } else {
      kme_frac <- (periph_kme - kme_range_p[1]) / diff(kme_range_p)
    }
    base_r <- r_outer - kme_frac * (r_outer - r_inner)
    r_jitter <- runif(n_periph, -0.08, 0.08) * (r_outer - r_inner)
    radii <- pmax(r_inner, pmin(r_outer, base_r + r_jitter))

    periph_x <- hub_cx + radii * cos(angles)
    periph_y <- hub_cy + radii * sin(angles)

    periph_gene <- module_df$gene[match(periph_prots, module_df$uniprot_id)]

    periph_df <- tibble(
      uniprot_id = periph_prots, module = mod,
      kME = periph_kme, GS = periph_gs,
      x = periph_x, y = periph_y,
      gene = periph_gene
    )
  }

  # --- Hull groups ---
  grp_counts <- table(lay$functional_group)
  hull_groups <- setdiff(names(grp_counts[grp_counts >= MIN_GROUP_N]), "Other")
  lay$hull_group <- ifelse(lay$functional_group %in% hull_groups,
                           lay$functional_group, NA)

  hull_color_map <- setNames(HULL_PALETTE[seq_along(hull_groups)], hull_groups)

  # --- Numbered centroid legend for hull labels ---
  hull_legend_text <- ""
  hull_centroids <- NULL
  if (length(hull_groups) > 0) {
    hull_data <- lay %>% filter(!is.na(hull_group))
    hull_centroids <- hull_data %>%
      group_by(hull_group) %>%
      summarise(x = mean(x), y = mean(y), .groups = "drop") %>%
      mutate(num = seq_len(n()))
    hull_legend_text <- paste(
      paste0(hull_centroids$num, "=", hull_centroids$hull_group),
      collapse = " | "
    )
  }

  # Build annotated ggraph plot
  p <- ggraph(lay) +
    annotate("rect", xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf,
             fill = mod, alpha = 0.04)

  if (length(hull_groups) > 0) {
    p <- p +
      geom_mark_hull(
        data = hull_data,
        aes(x = x, y = y, group = hull_group, fill = hull_group),
        concavity = 8, expand = unit(3, "mm"), radius = unit(3, "mm"),
        alpha = 0.10, linewidth = 0.3, colour = "grey50",
        show.legend = FALSE
      ) +
      scale_fill_manual(values = hull_color_map, guide = "none") +
      new_scale_fill()
  }

  p <- p +
    geom_edge_link0(aes(edge_width = weight), alpha = 0.10, color = "grey55") +
    scale_edge_width_continuous(range = c(0.15, 1.0), guide = "none")

  # Periphery ring layer (below hub nodes)
  if (!is.null(periph_df) && nrow(periph_df) > 0) {
    p <- p +
      geom_point(data = periph_df, aes(x = x, y = y, size = kME, fill = GS),
                 shape = 21, color = "grey60", stroke = PERIPH_STROKE,
                 alpha = PERIPH_ALPHA, inherit.aes = FALSE)
  }

  p <- p +
    geom_node_point(aes(size = kME, fill = GS), shape = 21,
                    color = "black", stroke = 0.3) +
    scale_size_continuous(range = c(1.2, 5.5), name = "kME") +
    scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B",
                         midpoint = 0, name = "GS")

  # Numbered centroid labels for hull groups
  if (!is.null(hull_centroids) && nrow(hull_centroids) > 0) {
    p <- p +
      geom_label(data = hull_centroids,
                 aes(x = x, y = y, label = num),
                 inherit.aes = FALSE,
                 size = 2.5, fontface = "bold",
                 fill = scales::alpha("white", 0.85),
                 linewidth = 0.2, label.padding = unit(1, "mm"))
  }

  p <- p +
    geom_node_text(aes(label = label_show), size = 2.2, repel = TRUE,
                   max.overlaps = 14, na.rm = TRUE,
                   fontface = "bold", color = "grey10",
                   bg.color = "white", bg.r = 0.12,
                   segment.size = 0.15, segment.color = "grey55",
                   box.padding = 0.2, point.padding = 0.15,
                   min.segment.length = 0.3)

  # Subtitle: GS phenotype + hull legend
  sub_lines <- sprintf("GS: %s", MODULE_GS_LABEL[[mod]])
  if (nchar(hull_legend_text) > 0) {
    sub_lines <- paste0(sub_lines, "\n", hull_legend_text)
  }

  p <- p +
    labs(title    = sprintf("%s (n=%d)", str_to_title(mod), n_mod),
         subtitle = sub_lines) +
    theme_void() +
    theme(
      plot.title    = element_text(face = "bold", size = 8),
      plot.subtitle = element_text(size = 5, color = "grey30", face = "italic"),
      legend.position = "none"
    )

  return(list(plot = p, node_data = node_df, periph_data = periph_df))
}

# === Build all 6 module networks ==============================================
cat("Building functionally annotated hub networks for 6 key modules...\n")
results <- lapply(KEY_MODULES, function(mod) {
  cat(sprintf("  Module: %s\n", mod))
  build_hub_network(mod, top_n = TOP_N, edge_quant = EDGE_QUANT)
})
names(results) <- KEY_MODULES

plots      <- lapply(results, `[[`, "plot")
node_data  <- lapply(results, `[[`, "node_data")
periph_data <- lapply(results, `[[`, "periph_data")

ok <- !sapply(plots, is.null)
plots      <- plots[ok]
node_data  <- node_data[ok]
periph_data <- periph_data[ok]

all_node_data  <- bind_rows(node_data)
all_periph_data <- bind_rows(periph_data)

# === Build trait sidebars =====================================================
cat("Building trait correlation sidebars...\n")
sidebars <- lapply(KEY_MODULES[ok], make_trait_sidebar)

# === Pair each network plot with its sidebar ==================================
paired_panels <- lapply(seq_along(plots), function(i) {
  plots[[i]] + sidebars[[i]] +
    plot_layout(widths = c(1, 0.14))
})

# === Shared legend (combined hub + periphery ranges) ==========================
combined_gs  <- c(all_node_data$GS, all_periph_data$GS)
combined_kme <- c(all_node_data$kME, all_periph_data$kME)
gs_range  <- range(combined_gs, na.rm = TRUE)
gs_limit  <- max(abs(gs_range)) * 1.05
kme_range <- range(combined_kme, na.rm = TRUE)

legend_df <- tibble(kME = combined_kme, GS = combined_gs)
legend_plot <- ggplot(legend_df, aes(x = kME, y = GS)) +
  geom_point(aes(size = kME, fill = GS), shape = 21, color = "black", stroke = 0.3) +
  scale_size_continuous(range = c(1.2, 5.5), name = "kME",
                        breaks = pretty(kme_range, n = 3)) +
  scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B",
                       midpoint = 0, limits = c(-gs_limit, gs_limit),
                       name = "Gene Significance") +
  guides(
    fill = guide_colorbar(barwidth = unit(30, "mm"), barheight = unit(2.5, "mm"),
                          title.position = "top", title.hjust = 0.5),
    size = guide_legend(title.position = "top", title.hjust = 0.5,
                        override.aes = list(fill = "grey70"))
  ) +
  theme_void() +
  theme(
    legend.position  = "bottom",
    legend.box       = "horizontal",
    legend.key.size  = unit(3, "mm"),
    legend.text      = element_text(size = 6),
    legend.title     = element_text(size = 7, face = "bold")
  )

shared_legend <- cowplot::get_legend(legend_plot)

# Sidebar legend: a mini colorbar for the trait heatmap
sidebar_legend <- ggplot(tibble(r = seq(-1, 1, 0.01)), aes(x = r, y = 1, fill = r)) +
  geom_tile() +
  scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B",
                       midpoint = 0, limits = c(-1, 1),
                       name = "Module-Trait r") +
  guides(fill = guide_colorbar(barwidth = unit(25, "mm"), barheight = unit(2.5, "mm"),
                                title.position = "top", title.hjust = 0.5)) +
  theme_void() +
  theme(legend.position = "bottom",
        legend.text = element_text(size = 6),
        legend.title = element_text(size = 7, face = "bold"))
sidebar_leg_grob <- cowplot::get_legend(sidebar_legend)

# Combine both legends
combined_legend <- cowplot::plot_grid(
  shared_legend, sidebar_leg_grob,
  nrow = 1, rel_widths = c(1, 0.7)
)

# === Grid assembly ============================================================
panel_grid <- wrap_plots(paired_panels, ncol = 3, nrow = 2)

panel_E <- panel_grid /
  wrap_elements(combined_legend) +
  plot_layout(heights = c(1, 0.07)) +
  plot_annotation(
    title    = "E  Hub Protein Networks",
    subtitle = paste0(
      "Top 30 by kME | node size = kME, color = gene significance | hulls = GO BP | sidebar = module-trait r\n",
      "Context ring: all remaining module proteins (size = kME, color = GS, alpha = 0.38) | radial position ~ kME"
    ),
    theme = theme(
      plot.title    = element_text(face = "bold", size = 12),
      plot.subtitle = element_text(size = 7.5, color = "grey30", face = "italic")
    )
  )

# ── STAT AUDIT ADDITION: cor.test() CIs for hub kME and GS ──────────────────
cat("Computing cor.test() CIs for hub kME and GS...\n")

hub_ci_list <- list()
for (i in seq_along(node_data)) {
  nd <- node_data[[i]]
  if (is.null(nd) || nrow(nd) == 0) next

  mod <- nd$module[1]
  me_col <- paste0("ME", mod)
  pheno_key <- MODULE_GS_PHENO[[mod]]

  # Phenotype vector for GS
  if (pheno_key == "age") {
    pheno_vec <- meta$age_binary[match(rownames(datExpr), meta$sample_id)]
  } else if (pheno_key == "training") {
    pheno_vec <- meta$time_binary[match(rownames(datExpr), meta$sample_id)]
  } else {
    pheno_vec <- meta[[pheno_key]][match(rownames(datExpr), meta$sample_id)]
  }

  for (j in seq_len(nrow(nd))) {
    uid <- nd$uniprot_id[j]
    if (!(uid %in% colnames(datExpr))) next

    prot_expr <- datExpr[, uid]

    # kME CI: cor(protein, module eigengene)
    me_vec <- MEs[, me_col]
    kme_ct <- cor.test(prot_expr, me_vec, method = "pearson")

    # GS CI: cor(protein, phenotype)
    complete_gs <- complete.cases(prot_expr, pheno_vec)
    gs_ct <- if (sum(complete_gs) >= 4) {
      cor.test(prot_expr[complete_gs], pheno_vec[complete_gs], method = "pearson")
    } else {
      NULL
    }

    hub_ci_list <- c(hub_ci_list, list(tibble(
      module = mod,
      uniprot_id = uid,
      gene = nd$gene[j],
      kME = round(kme_ct$estimate, 4),
      kME_ci_lo = round(kme_ct$conf.int[1], 4),
      kME_ci_hi = round(kme_ct$conf.int[2], 4),
      kME_p = kme_ct$p.value,
      GS_pheno = nd$GS_pheno[j],
      GS = if (!is.null(gs_ct)) round(gs_ct$estimate, 4) else NA_real_,
      GS_ci_lo = if (!is.null(gs_ct)) round(gs_ct$conf.int[1], 4) else NA_real_,
      GS_ci_hi = if (!is.null(gs_ct)) round(gs_ct$conf.int[2], 4) else NA_real_,
      GS_p = if (!is.null(gs_ct)) gs_ct$p.value else NA_real_,
      GS_n = sum(complete_gs)
    )))
  }
}

hub_ci_df <- bind_rows(hub_ci_list)
write_csv(hub_ci_df, file.path(DAT_DIR, "05_panel_E_hub_CIs.csv"))
cat(sprintf("  Exported: 05_panel_E_hub_CIs.csv (%d hub proteins with 95%% CIs)\n",
            nrow(hub_ci_df)))
# ─────────────────────────────────────────────────────────────────────────────

# === Save outputs =============================================================
write_csv(all_node_data, file.path(DAT_DIR, "05_panel_E_hub_network.csv"))
write_csv(all_periph_data, file.path(DAT_DIR, "05_panel_E_context_ring.csv"))

ggsave(file.path(RPT_DIR, "panel_E_hub_network.pdf"), panel_E,
       width = 450, height = 315, units = "mm",
       device = pdf, limitsize = FALSE)
ggsave(file.path(RPT_DIR, "panel_E_hub_network.png"), panel_E,
       width = 450, height = 315, units = "mm",
       dpi = 300, limitsize = FALSE)

cat("Panel E saved (6-module annotated hub network grid with context rings + trait sidebars)\n")
