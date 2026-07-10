# Sourced by 02_supp_panels.R — expects style.R + pathway_utils.R already loaded.

pacman::p_load(tidyverse, patchwork, ggrepel, WGCNA, igraph, ggraph, ggforce, concaveman, graphlayouts, tidygraph, ggnewscale)

allowWGCNAThreads()
set.seed(42)

BASE <- "04_Figures/F06"

RPT_SUPP_PNG <- file.path(BASE, "b_reports", "supp", "png", "modules")
RPT_SUPP_PDF <- file.path(BASE, "b_reports", "supp", "pdf", "modules")
DAT          <- file.path(BASE, "c_data")
dir.create(RPT_SUPP_PNG, recursive = TRUE, showWarnings = FALSE)
dir.create(RPT_SUPP_PDF, recursive = TRUE, showWarnings = FALSE)

pdf_device <- get_pdf_device()

stopifnot(
  "WGCNA gs_phenotype_choices.csv missing — run YvO_WGCNA_run.R first" =
    file.exists(file.path(DAT, "wgcna/gs_phenotype_choices.csv")),
  "meta.csv missing — run YvO_WGCNA_run.R first" =
    file.exists(file.path(DAT, "meta.csv")),
  "MEs.rds missing — run YvO_WGCNA_run.R first" =
    file.exists(file.path(DAT, "MEs.rds")),
  "kME_all.rds missing — run YvO_WGCNA_run.R first" =
    file.exists(file.path(DAT, "kME_all.rds")),
  "datExpr.rds missing — run YvO_WGCNA_run.R first" =
    file.exists(file.path(DAT, "datExpr.rds")),
  "WGCNA module assignments missing — run YvO_WGCNA_run.R first" =
    file.exists(file.path(DAT, "wgcna/wgcna_module_assignments.csv")),
  "WGCNA sft_summary missing — run YvO_WGCNA_run.R first" =
    file.exists(file.path(DAT, "wgcna/wgcna_sft_summary.csv"))
)

meta <- read_csv(file.path(DAT, "meta.csv"))
meta$group <- factor(meta$group,
                     levels = c("Young_Pre", "Young_Post", "Old_Pre", "Old_Post"))
MEs     <- readRDS(file.path(DAT, "MEs.rds"))
kME_all <- readRDS(file.path(DAT, "kME_all.rds"))
datExpr <- readRDS(file.path(DAT, "datExpr.rds"))
module_df <- read_csv(file.path(DAT, "wgcna/wgcna_module_assignments.csv"))
sft_csv   <- read.csv(file.path(DAT, "wgcna/wgcna_sft_summary.csv"))
NET_POWER <- sft_csv$selected_power[1]

KEY_MODULES <- module_df |>
  filter(module_color != "grey") |>
  count(module_color, sort = TRUE) |>
  pull(module_color)
bg_genes <- unique(module_df$gene)

mod_bio_labels_df  <- read_csv(file.path(DAT, "mod_bio_labels.csv"))
mod_bio_labels_vec <- setNames(mod_bio_labels_df$bio_label, mod_bio_labels_df$module_color)
display_label_vec  <- setNames(mod_bio_labels_df$display_label, mod_bio_labels_df$module_color)

meta$age_binary  <- ifelse(meta$age == "Old", 1, 0)
meta$age_num     <- meta$age_binary

gs_choices <- read_csv(file.path(DAT, "wgcna/gs_phenotype_choices.csv"))
MODULE_GS_PHENO <- setNames(gs_choices$gs_phenotype, gs_choices$module)
gs_label_map <- c(age_num = "Age (Young=0, Old=1)",
                  delta_VL = "\u0394VL Thickness (cm)",
                  delta_LBM = "\u0394Lean Body Mass (kg)",
                  VL_thick_cm = "VL Thickness (cm)",
                  DXA_LBM_kg = "Lean Body Mass (kg)",
                  BMI = "BMI (kg/m\u00b2)")
MODULE_GS_LABEL <- setNames(
  ifelse(gs_choices$gs_phenotype %in% names(gs_label_map),
         gs_label_map[gs_choices$gs_phenotype],
         tools::toTitleCase(gsub("_", " ", gs_choices$gs_phenotype))),
  gs_choices$module
)
for (m in KEY_MODULES) {
  if (!(m %in% names(MODULE_GS_PHENO))) {
    MODULE_GS_PHENO[[m]] <- "age_num"
    MODULE_GS_LABEL[[m]] <- "Age (Young=0, Old=1)"
  }
}

if ("delta_VL" %in% MODULE_GS_PHENO && !("delta_VL" %in% colnames(meta))) {
  vl_wide <- meta |>
    filter(!is.na(VL_thick_cm)) |>
    select(subject, time, VL_thick_cm) |>
    pivot_wider(names_from = time, values_from = VL_thick_cm, names_prefix = "VL_") |>
    mutate(delta_VL = VL_Post - VL_Pre)
  meta <- meta |> left_join(vl_wide |> select(subject, delta_VL), by = "subject")
  message(sprintf("  Computed delta_VL for %d subjects (non-NA: %d)",
                  nrow(vl_wide), sum(!is.na(vl_wide$delta_VL))))
}

uid2gene    <- setNames(module_df$gene, module_df$uniprot_id)

PD_W <- 170; PD_H <- 170
txt_gene  <- scale_text(BASE_GENE, PD_W)
txt_title <- scale_text(BASE_STAT, PD_W) * 1.8
txt_sub   <- scale_text(BASE_STAT, PD_W) * 1.3
txt_leg   <- scale_text(BASE_STAT, PD_W) * 1.4
txt_legt  <- scale_text(BASE_STAT, PD_W) * 1.2

HULL_PALETTE <- c("#1B9E77", "#D95F02", "#7570B3", "#E7298A",
                  "#66A61E", "#E6AB02", "#A6761D", "#666666")

message("Panel D: hub protein networks...")

pw_full <- build_pathway_collection(min_size = 15, max_size = 500, include_goslim = FALSE)

select_hubs_q90 <- function(mod) {
  mod_prots <- module_df$uniprot_id[module_df$module_color == mod]
  kme_col   <- paste0("kME", mod)
  matched   <- intersect(mod_prots, rownames(kME_all))
  mod_kme   <- setNames(kME_all[matched, kme_col], matched)
  mod_kme   <- mod_kme[!is.na(mod_kme)]
  q90       <- quantile(mod_kme, 0.90)
  names(mod_kme[mod_kme >= q90])
}

assign_groups_ora <- function(gene_names, max_groups = 4, min_group_n = 3) {
  clean_pw_name <- function(name) {
    name |>
      (\(x) gsub("^HALLMARK_|^GOSLIM_|^GOBP_|^REACTOME_|^KEGG_MEDICUS_", "", x))() |>
      (\(x) gsub("_", " ", x))() |>
      str_to_title() |>
      str_trunc(35)
  }

  ora_res <- tryCatch(
    run_ora_deduplicated(
      genes    = gene_names,
      universe = bg_genes,
      pathways = pw_full,
      jaccard_cutoff = 0.5,
      min_size = 10, max_size = 500,
      padj_cutoff = 0.05
    ),
    error = function(e) { message("  ORA error: ", e$message); NULL }
  )

  if (is.null(ora_res) || nrow(ora_res) == 0) {
    message("  No significant ORA results")
    return(setNames(rep("Other", length(gene_names)), gene_names))
  }

  ora_res <- ora_res[order(ora_res$padj), ]
  gene_map <- data.frame(gene = character(), pathway = character(), stringsAsFactors = FALSE)
  for (i in seq_len(nrow(ora_res))) {
    hits <- intersect(ora_res$overlapGenes[[i]], gene_names)
    if (length(hits))
      gene_map <- rbind(gene_map, data.frame(gene = hits, pathway = ora_res$pathway[i],
                                             stringsAsFactors = FALSE))
  }
  gene_map <- gene_map[!duplicated(gene_map$gene), ]

  term_counts <- table(gene_map$pathway)
  keep <- names(term_counts[term_counts >= min_group_n])
  keep <- head(keep[order(term_counts[keep], decreasing = TRUE)], max_groups)

  if (length(keep) < 2) {
    message("  Fewer than 2 qualifying groups")
    return(setNames(rep("Other", length(gene_names)), gene_names))
  }

  assignments <- setNames(rep("Other", length(gene_names)), gene_names)
  for (g in gene_names) {
    row <- gene_map[gene_map$gene == g, ]
    if (nrow(row) > 0 && row$pathway[1] %in% keep)
      assignments[g] <- clean_pw_name(row$pathway[1])
  }
  assignments
}

compute_gs_signed <- function(mod) {
  pheno_col <- MODULE_GS_PHENO[[mod]]
  pheno_vec <- meta[[pheno_col]]
  if (is.null(pheno_vec)) {
    warning(sprintf("  GS phenotype '%s' not found in meta for %s; using age_num", pheno_col, mod))
    pheno_vec <- meta[["age_num"]]
    if (is.null(pheno_vec)) pheno_vec <- meta[["age_binary"]]
  }
  names(pheno_vec) <- meta$sample_id
  valid_samps <- intersect(meta$sample_id[!is.na(pheno_vec)], rownames(datExpr))
  gs <- cor(datExpr[valid_samps, , drop = FALSE], pheno_vec[valid_samps],
            use = "pairwise.complete.obs")
  setNames(gs[, 1], rownames(gs))
}

build_network_hull <- function(mod) {

  message(sprintf("\n=== %s ===", toupper(mod)))

  hub_ids <- select_hubs_q90(mod)
  hub_genes <- uid2gene[hub_ids]; hub_genes <- hub_genes[!is.na(hub_genes)]
  hub_ids <- hub_ids[hub_ids %in% names(hub_genes)]
  n_mod <- sum(module_df$module_color == mod)
  n_hubs <- length(hub_ids)
  message(sprintf("  %d hubs (Q90 of %d)", n_hubs, n_mod))

  mod_prots <- intersect(module_df$uniprot_id[module_df$module_color == mod],
                         colnames(datExpr))
  adj_mod <- adjacency(datExpr[, mod_prots], power = NET_POWER, type = "signed hybrid")
  tom_mod <- TOMsimilarity(adj_mod, TOMType = "signed")
  colnames(tom_mod) <- rownames(tom_mod) <- mod_prots

  tom_sub <- tom_mod[hub_ids, hub_ids]
  tom_q90 <- quantile(tom_sub[upper.tri(tom_sub)], 0.90)

  g <- graph_from_adjacency_matrix(tom_sub, mode = "undirected", weighted = TRUE, diag = FALSE)
  E(g)$weight_orig <- E(g)$weight
  g <- delete_edges(g, which(E(g)$weight < tom_q90))
  iso <- which(degree(g) == 0)
  if (length(iso)) g <- delete_vertices(g, iso)

  node_uids  <- V(g)$name
  node_genes <- uid2gene[node_uids]
  node_kme   <- setNames(kME_all[node_uids, paste0("kME", mod)], node_uids)
  gs_signed  <- compute_gs_signed(mod)
  node_gs    <- gs_signed[node_uids]; node_gs[is.na(node_gs)] <- 0

  top_kme <- names(sort(node_kme, decreasing = TRUE))[seq_len(min(6, length(node_kme)))]
  top_gs  <- names(sort(abs(node_gs), decreasing = TRUE))[seq_len(min(3, length(node_gs)))]
  label_ids <- unique(c(top_kme, top_gs))

  groups <- assign_groups_ora(node_genes)
  node_groups <- groups[node_genes]

  V(g)$gene <- node_genes; V(g)$kME <- node_kme
  V(g)$gs <- node_gs; V(g)$func_grp <- node_groups

  set.seed(42)
  lay <- layout_with_stress(g)
  lay_df <- data.frame(x = lay[, 1], y = lay[, 2], name = V(g)$name)

  nd <- data.frame(x = lay_df$x, y = lay_df$y, name = node_uids,
                   gene = node_genes, kME = node_kme, GS = node_gs,
                   func_grp = node_groups, stringsAsFactors = FALSE)

  grp_counts <- table(nd$func_grp)
  nd$n_in_grp <- as.integer(grp_counts[nd$func_grp])

  grp_names <- setdiff(unique(nd$func_grp[nd$func_grp != "Other" & nd$n_in_grp >= 3]), NA)
  hull_colors <- setNames(HULL_PALETTE[seq_along(grp_names)], grp_names)

  tg <- as_tbl_graph(g)

  disp_label <- display_label_vec[[mod]]
  pheno_label <- MODULE_GS_LABEL[[mod]]

  hull_nd <- nd |> filter(func_grp != "Other", n_in_grp >= 3)

  p <- ggraph(tg, layout = "manual", x = lay_df$x, y = lay_df$y)

  if (nrow(hull_nd) > 0 && length(grp_names) > 0) {
    p <- p +
      geom_mark_hull(
        data = hull_nd,
        aes(x = x, y = y, group = func_grp, fill = func_grp),
        concavity = 2, expand = unit(2, "mm"), radius = unit(2, "mm"),
        alpha = 0.15, linewidth = 0.6, show.legend = TRUE,
        inherit.aes = FALSE
      ) +
      scale_fill_manual(values = hull_colors, name = "Pathway")
  }

  p <- p +
    geom_edge_link(aes(width = weight_orig), alpha = 0.55,
                   color = "grey30", show.legend = FALSE) +
    scale_edge_width_continuous(range = c(0.5, 2.0), name = "TOM",
                               guide = guide_legend(keywidth = unit(10, "mm")))

  p <- p +
    new_scale_fill() +
    geom_node_point(aes(size = kME, fill = gs), shape = 21,
                    color = "black", stroke = 0.5) +
    scale_size_continuous(range = c(2.0, 7.0), guide = "none") +
    scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B",
                         midpoint = 0, name = "GS",
                         guide = guide_colorbar(barwidth = unit(3, "mm"),
                                                barheight = unit(20, "mm")))

  p <- p +
    geom_label_repel(
      data = nd[nd$name %in% label_ids, ],
      aes(x = x, y = y, label = gene),
      size = txt_gene, fontface = "bold.italic",
      fill = alpha("white", 0.88), color = "grey10",
      linewidth = 0.15, label.padding = unit(1.0, "mm"),
      segment.size = 0.25, segment.color = "grey40",
      box.padding = 0.45, point.padding = 0.2,
      max.overlaps = 20, seed = 42, inherit.aes = FALSE
    )

  p <- p +
    labs(title    = disp_label,
         subtitle = sprintf("%d hubs | Node color: GS (%s)", n_hubs, pheno_label)) +
    theme_void() +
    theme(
      plot.title      = element_text(face = "bold", size = txt_title, hjust = 0.5),
      plot.subtitle   = element_text(size = txt_sub, hjust = 0.5, color = "grey40"),
      plot.background = element_rect(fill = "white", color = NA),
      plot.margin     = margin(4, 4, 4, 4),
      legend.position = "right",
      legend.title    = element_text(face = "bold", size = txt_leg),
      legend.text     = element_text(size = txt_legt)
    )

  list(plot = p, node_data = nd)
}

results <- lapply(KEY_MODULES, function(mod) build_network_hull(mod))
names(results) <- KEY_MODULES

hub_letters <- setNames(LETTERS[6:(5 + length(KEY_MODULES))], KEY_MODULES)
pathway_slug <- setNames(
  gsub("[/ ]+", "_", tolower(mod_bio_labels_vec)),
  names(mod_bio_labels_vec)
)

for (mod in KEY_MODULES) {
  res <- results[[mod]]
  slug <- pathway_slug[mod]
  if (is.na(slug)) slug <- mod
  fname <- sprintf("SUPP_hub_%s_%s", mod, slug)
  ggsave(file.path(RPT_SUPP_PNG, paste0(fname, ".png")), res$plot,
         width = PD_W, height = PD_H, units = "mm", dpi = 300)
  ggsave(file.path(RPT_SUPP_PDF, paste0(fname, ".pdf")), res$plot,
         width = PD_W, height = PD_H, units = "mm", device = pdf_device)
  message(sprintf("  Saved %s", fname))
}

node_data_list <- lapply(results, `[[`, "node_data")

all_node_df <- bind_rows(lapply(KEY_MODULES, function(mod) {
  nd <- node_data_list[[mod]]
  if (is.null(nd) || nrow(nd) == 0) return(NULL)
  tibble(uniprot_id = nd$name, module = mod, gene = nd$gene,
         kME = nd$kME, GS = nd$GS, functional_group = nd$func_grp,
         GS_pheno = MODULE_GS_LABEL[[mod]])
}))
write_csv(all_node_df, file.path(DAT, "04_panel_D_hub_network.csv"))

hub_ci_list <- list()
for (mod in KEY_MODULES) {
  nd <- node_data_list[[mod]]
  if (is.null(nd) || nrow(nd) == 0) next
  me_col <- paste0("ME", mod)
  pheno_col <- MODULE_GS_PHENO[[mod]]
  pheno_vec <- meta[[pheno_col]]
  if (is.null(pheno_vec)) pheno_vec <- meta[["age_num"]]
  names(pheno_vec) <- meta$sample_id
  pheno_vec <- pheno_vec[intersect(names(pheno_vec), rownames(datExpr))]

  for (j in seq_len(nrow(nd))) {
    uid <- nd$name[j]
    if (!(uid %in% colnames(datExpr))) next
    prot_expr <- datExpr[, uid]
    kme_ct <- cor.test(prot_expr, MEs[, me_col], method = "pearson")
    valid <- intersect(names(pheno_vec), names(prot_expr))
    gs_ct <- if (length(valid) >= 4) cor.test(prot_expr[valid], pheno_vec[valid], method = "pearson") else NULL
    hub_ci_list <- c(hub_ci_list, list(tibble(
      module = mod, uniprot_id = uid, gene = nd$gene[j],
      kME = round(kme_ct$estimate, 4),
      kME_ci_lo = round(kme_ct$conf.int[1], 4),
      kME_ci_hi = round(kme_ct$conf.int[2], 4),
      kME_p = kme_ct$p.value,
      GS_pheno = MODULE_GS_LABEL[[mod]],
      GS = if (!is.null(gs_ct)) round(gs_ct$estimate, 4) else NA_real_,
      GS_ci_lo = if (!is.null(gs_ct)) round(gs_ct$conf.int[1], 4) else NA_real_,
      GS_ci_hi = if (!is.null(gs_ct)) round(gs_ct$conf.int[2], 4) else NA_real_,
      GS_p = if (!is.null(gs_ct)) gs_ct$p.value else NA_real_,
      GS_n = length(valid)
    )))
  }
}
write_csv(bind_rows(hub_ci_list), file.path(DAT, "04_panel_D_hub_CIs.csv"))

message("  Panel D saved (individual modules)")

message("  Building Panel E composite...")

plots_bare <- lapply(KEY_MODULES, function(mod) {
  results[[mod]]$plot + theme(legend.position = "none",
                               plot.margin = margin(2, 2, 2, 2))
})

n_mods  <- length(plots_bare)
n_cols  <- min(3, n_mods)
n_rows  <- ceiling(n_mods / n_cols)

panel_E <- wrap_plots(plots_bare, ncol = n_cols) +
  plot_layout(widths = rep(1, n_cols))

PE_W <- n_cols * 170
PE_H <- n_rows * 180

ggsave(file.path(RPT_SUPP_PNG, "SUPP_networks_composite.png"), panel_E,
       width = PE_W, height = PE_H, units = "mm", dpi = 300, limitsize = FALSE)
ggsave(file.path(RPT_SUPP_PDF, "SUPP_networks_composite.pdf"), panel_E,
       width = PE_W, height = PE_H, units = "mm", device = pdf_device, limitsize = FALSE)
message("  Supp S networks composite saved")

message("  Panel D complete (individual + composite)")
