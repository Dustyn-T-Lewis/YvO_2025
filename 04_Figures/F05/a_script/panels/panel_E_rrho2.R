# F05 Panel E: RRHO2 Reversal Map + Per-Quadrant ORA
# Rebuilt to use RRHO2 package (matching F04 construction logic)
# Contrasts: Aging vs Training_Old
setwd(rprojroot::find_rstudio_root_file())
source("04_Figures/shared/style.R")
source("04_Figures/shared/pathway_utils.R")
library(tidyverse)
library(patchwork)
library(ggrepel)
library(msigdbr)
library(fgsea)
library(RRHO2)

PE_W <- 260

RPT     <- "04_Figures/F05/b_reports/panels"
RPT_SUP <- "04_Figures/F05/b_reports/supp/panels"
DAT <- "04_Figures/F05/c_data"
dir.create(RPT, recursive = TRUE, showWarnings = FALSE)
dir.create(RPT_SUP, recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(DAT, "panel_E"), recursive = TRUE, showWarnings = FALSE)

pdf_device <- get_pdf_device()

dep_df <- read_csv("03_DEP/c_data/03_combined_results.csv", show_col_types = FALSE)

rr_df <- dep_df %>%
  transmute(gene, t_aging = t_Aging, t_old = t_Training_Old) %>%
  filter(!is.na(t_aging) & !is.na(t_old)) %>%
  distinct(gene, .keep_all = TRUE)

n_shared <- nrow(rr_df)

# --- RRHO2 computation (stratified hypergeometric, Cahill et al. 2018)
# Both lists sorted descending internally by RRHO2_initialize:
#   list1 (Aging):        index 1 = most UP with aging
#   list2 (Training_Old): index 1 = most UP with training
list1 <- data.frame(gene = rr_df$gene, score = rr_df$t_aging, stringsAsFactors = FALSE)
list2 <- data.frame(gene = rr_df$gene, score = rr_df$t_old,   stringsAsFactors = FALSE)

rrho_obj <- RRHO2_initialize(
  list1, list2,
  labels          = c("Aging", "Training (Old)"),
  log10.ind       = TRUE,
  multipleTesting = "none",
  boundary        = 0.02,
  method          = "hyper",
  stepsize        = 20
)

hmat <- rrho_obj$hypermat
nr <- nrow(hmat); nc <- ncol(hmat)
message(sprintf("  RRHO2 matrix: %d x %d", nr, nc))

# Locate the NA boundary strip (zero-crossing separation)
na_rows <- which(apply(hmat, 1, function(r) all(is.na(r))))
na_cols <- which(apply(hmat, 2, function(c) all(is.na(c))))

if (length(na_rows) > 0 && length(na_cols) > 0) {
  row_before <- 1:(min(na_rows) - 1)
  row_after  <- (max(na_rows) + 1):nr
  col_before <- 1:(min(na_cols) - 1)
  col_after  <- (max(na_cols) + 1):nc
} else {
  mid <- floor(nr / 2)
  row_before <- 1:mid
  row_after  <- (mid + 1):nr
  col_before <- 1:mid
  col_after  <- (mid + 1):nc
}

# Both lists sorted descending (Up first):
# [row_before, col_before] = top of both lists = Aging Up + Training Up = Exacerbated Up
# [row_after,  col_after]  = bottom of both    = Aging Down + Training Down = Exacerbated Down
# [row_before, col_after]  = top Aging, bottom Training = Aging Up + Training Down = REVERSED
# [row_after,  col_before] = bottom Aging, top Training = Aging Down + Training Up = REVERSED
max_UU <- max(hmat[row_before, col_before], na.rm = TRUE)  # Exacerbated Up
max_DD <- max(hmat[row_after,  col_after],  na.rm = TRUE)  # Exacerbated Down
max_UD <- max(hmat[row_before, col_after],  na.rm = TRUE)  # Reversed (Aging Up / Tr Down)
max_DU <- max(hmat[row_after,  col_before], na.rm = TRUE)  # Reversed (Aging Down / Tr Up)

# Extract hotspot genes from RRHO2 quadrant gene lists
hotspot_genes <- list(
  UU = rrho_obj$genelist_uu$gene_list_overlap_uu,
  DD = rrho_obj$genelist_dd$gene_list_overlap_dd,
  UD = rrho_obj$genelist_ud$gene_list_overlap_ud,
  DU = rrho_obj$genelist_du$gene_list_overlap_du
)

n_UU <- length(hotspot_genes$UU)
n_DD <- length(hotspot_genes$DD)
n_UD <- length(hotspot_genes$UD)
n_DU <- length(hotspot_genes$DU)
message(sprintf("  Hotspot genes: UU(Exac Up)=%d, DD(Exac Dn)=%d, UD(Rev Up)=%d, DU(Rev Dn)=%d",
                n_UU, n_DD, n_UD, n_DU))

txt_quad <- scale_text(BASE_QUADRANT, PE_W)

# --- Jet colormap heatmap (Cahill et al. 2018 canonical appearance)
JET_COLORS <- c("#00007F", "blue", "#007FFF", "cyan", "#7FFF7F",
                "yellow", "#FF7F00", "red", "#7F0000")

hmat_df <- expand.grid(row = 1:nr, col = 1:nc) %>%
  mutate(neg_log10_p = as.vector(hmat))

max_val <- max(hmat_df$neg_log10_p, na.rm = TRUE)

# Corner anchors — labels with hjust/vjust = 0/1 sit flush against
# the outer corner of each quadrant (matching Panel A scatter style)
ann_x_left  <- min(row_before)
ann_x_right <- max(row_after)
ann_y_bot   <- min(col_before)
ann_y_top   <- max(col_after)

# Soft white background keeps labels readable over high-intensity heatmap regions
LABEL_FILL    <- scales::alpha("white", 0.85)
LABEL_PADDING <- unit(0.8, "mm")

pE_heat <- ggplot(hmat_df, aes(x = row, y = col, fill = neg_log10_p)) +
  geom_raster() +
  scale_fill_gradientn(
    colors   = JET_COLORS,
    limits   = c(0, max_val),
    na.value = "white",
    name     = expression(-log[10](P)),
    guide    = guide_colorbar(
      barwidth = unit(34, "mm"), barheight = unit(3.5, "mm"),
      title.position = "left", title.vjust = 0.5,
      title.theme = element_text(size = 10.2, face = "bold", color = "grey15"))
  ) +
  # Exacerbated Up (Aging Up + Training Up) — bottom-left
  annotate("label", x = ann_x_left, y = ann_y_bot,
           label = "Exacerbated Up",
           color = "grey15", fill = LABEL_FILL, label.size = 0,
           label.padding = LABEL_PADDING, label.r = unit(0.5, "mm"),
           fontface = "bold", size = txt_quad * 0.80,
           hjust = 0, vjust = 0) +
  # Exacerbated Down (Aging Down + Training Down) — top-right
  annotate("label", x = ann_x_right, y = ann_y_top,
           label = "Exacerbated Down",
           color = "grey15", fill = LABEL_FILL, label.size = 0,
           label.padding = LABEL_PADDING, label.r = unit(0.5, "mm"),
           fontface = "bold", size = txt_quad * 0.80,
           hjust = 1, vjust = 1) +
  # Reversed: Aging Up / Training Down — top-left
  annotate("label", x = ann_x_left, y = ann_y_top,
           label = "Reversed (Age\u2191 Tr\u2193)",
           color = "grey15", fill = LABEL_FILL, label.size = 0,
           label.padding = LABEL_PADDING, label.r = unit(0.5, "mm"),
           fontface = "bold", size = txt_quad * 0.80,
           hjust = 0, vjust = 1) +
  # Reversed: Aging Down / Training Up — bottom-right
  annotate("label", x = ann_x_right, y = ann_y_bot,
           label = "Reversed (Age\u2193 Tr\u2191)",
           color = "grey15", fill = LABEL_FILL, label.size = 0,
           label.padding = LABEL_PADDING, label.r = unit(0.5, "mm"),
           fontface = "bold", size = txt_quad * 0.80,
           hjust = 1, vjust = 0) +
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0)) +
  labs(
    title    = "Threshold-Free Reversal (RRHO2)",
    subtitle = sprintf("Stratified hypergeometric | %d shared genes | warm off-diagonal = training reverses aging | No MTC (Cahill et al. 2018)", n_shared),
    x = expression("Aging rank"~(Up %->% Down)),
    y = expression("Training (Old) rank"~(Up %->% Down))
  ) +
  FIG_THEME +
  theme(
    axis.text        = element_blank(),
    axis.title.x     = element_text(margin = margin(t = 2)),
    axis.title.y     = element_text(margin = margin(r = 2)),
    axis.ticks       = element_blank(),
    panel.border     = element_blank(),
    panel.grid.major = element_blank(),
    legend.position  = "bottom",
    legend.text      = element_text(size = 9, face = "bold"),
    legend.margin    = margin(-2, 6, 0, 0),
    plot.margin = margin(2, 2, 2, 2, "mm")
  ) +
  coord_fixed(ratio = 1)

# --- Export hotspot genes
hotspot_export <- bind_rows(
  tibble(quadrant = "Exacerbated_Up",              gene = hotspot_genes$UU),
  tibble(quadrant = "Exacerbated_Down",             gene = hotspot_genes$DD),
  tibble(quadrant = "Reversed_AgingUp_TrainingDown", gene = hotspot_genes$UD),
  tibble(quadrant = "Reversed_AgingDown_TrainingUp", gene = hotspot_genes$DU)
)
write_csv(hotspot_export, file.path(DAT, "panel_E", "rrho2_hotspot_genes.csv"))

# --- Per-quadrant ORA
pw_collection_E <- build_pathway_collection(min_size = 10, max_size = 500,
                                             include_goslim = FALSE,
                                             exclude_variants = TRUE)
all_genes_E <- rr_df$gene

run_quadrant_ora <- function(gene_set, quadrant_name) {
  if (length(gene_set) < 5) return(tibble())
  res <- tryCatch(
    run_ora_deduplicated(
      genes          = gene_set,
      universe       = all_genes_E,
      pathways       = pw_collection_E,
      jaccard_cutoff = 0.5,
      min_size       = 10,
      max_size       = 500,
      padj_cutoff    = 0.05
    ),
    error = function(e) { message("  ORA error: ", e$message); tibble() }
  )
  if (nrow(res) > 0) {
    res %>%
      mutate(quadrant = quadrant_name,
             pathway_label = clean_pathway_name(pathway),
             ID = pathway,
             p.adjust = padj,
             GeneRatio = paste0(overlap, "/", length(gene_set)),
             geneID = sapply(overlapGenes, paste, collapse = "/")) %>%
      arrange(padj, size)
  } else {
    tibble()
  }
}

ora_UU <- run_quadrant_ora(hotspot_genes$UU, "Exacerbated Up")
ora_DD <- run_quadrant_ora(hotspot_genes$DD, "Exacerbated Down")
ora_UD <- run_quadrant_ora(hotspot_genes$UD, "Reversed (Aging Up / Training Down)")
ora_DU <- run_quadrant_ora(hotspot_genes$DU, "Reversed (Aging Down / Training Up)")

ora_reversal     <- bind_rows(ora_UD, ora_DU)
ora_exacerbation <- bind_rows(ora_UU, ora_DD)

write_csv(ora_reversal,     file.path(DAT, "panel_E", "rrho2_ora_concordant.csv"))
write_csv(ora_exacerbation, file.path(DAT, "panel_E", "rrho2_ora_discordant.csv"))

if (nrow(ora_exacerbation) == 0) {
  write_csv(tibble(note = "No pathways enriched in exacerbation quadrants (padj<0.05)"),
            file.path(DAT, "panel_E", "rrho2_ora_discordant_note.csv"))
}

# Save RRHO heatmap (square, since coord_fixed)
ggsave(file.path(RPT, "MAIN_panel_E_rrho2.png"), pE_heat,
       width = PE_W, height = PE_W, units = "mm", dpi = 300)

# --- Supplementary: combined heatmap + ORA bar chart
txt_ora <- scale_text(BASE_STAT, PE_W)
ora_all <- bind_rows(ora_reversal, ora_exacerbation)

if (nrow(ora_all) > 0) {

  all_quadrant_names <- c("Reversed (Aging Up / Training Down)",
                          "Reversed (Aging Down / Training Up)",
                          "Exacerbated Up", "Exacerbated Down")
  quad_short <- c(
    "Reversed (Aging Up / Training Down)"  = "Reversed\n(Up \u2192 Down)",
    "Reversed (Aging Down / Training Up)"  = "Reversed\n(Down \u2192 Up)",
    "Exacerbated Up"                       = "Exacerbated Up",
    "Exacerbated Down"                     = "Exacerbated Down"
  )

  MAX_PER_QUAD <- 12
  bar_df <- ora_all %>%
    mutate(
      neg_log10_padj = -log10(p.adjust),
      pathway_label  = str_trunc(clean_pathway_name(pathway), 40),
      quadrant       = factor(quadrant, levels = all_quadrant_names)
    ) %>%
    filter(!is.na(quadrant)) %>%
    group_by(quadrant, .drop = FALSE) %>%
    arrange(desc(neg_log10_padj)) %>%
    slice_head(n = MAX_PER_QUAD) %>%
    ungroup() %>%
    filter(!is.na(neg_log10_padj)) %>%
    arrange(quadrant, neg_log10_padj) %>%
    mutate(uid = fct_inorder(paste0(pathway_label, "___", quadrant)))

  n_shown <- nrow(bar_df)
  n_total <- nrow(ora_all)

  pE_ora <- ggplot(bar_df, aes(x = neg_log10_padj, y = uid, fill = quadrant)) +
    geom_col(width = 0.75) +
    geom_text(aes(label = overlap), hjust = -0.3, size = txt_ora * 0.7,
              color = "grey30") +
    scale_y_discrete(labels = function(x) str_remove(x, "___.*$")) +
    scale_fill_manual(values = ORA_QUAD_COLORS_F3, guide = "none") +
    scale_x_continuous(expand = expansion(mult = c(0, 0))) +
    facet_grid(quadrant ~ ., scales = "free_y", space = "free_y",
               labeller = labeller(quadrant = quad_short)) +
    labs(title = "Enriched Pathways by Reversal Quadrant",
         subtitle = if (n_shown < n_total)
           sprintf("Top %d per quadrant (%d terms total)", MAX_PER_QUAD, n_total)
         else sprintf("%d terms total", n_total),
         x = expression(-log[10](p[adj])),
         y = NULL) +
    FIG_THEME +
    theme(
      plot.title       = element_text(size = 10, face = "bold", hjust = 0.5),
      plot.subtitle    = element_text(size = 8, hjust = 0.5),
      strip.text.y     = element_text(size = 7, face = "bold", angle = 0),
      strip.background = element_rect(fill = "grey95", color = NA),
      panel.grid.major.y = element_blank(),
      panel.grid.major.x = element_line(color = "grey92", linewidth = 0.3),
      panel.grid.minor = element_blank(),
      axis.text.y  = element_text(size = 7),
      plot.margin  = margin(2, 4, 2, 2, "mm")
    )

  pE_combined <- (pE_heat | pE_ora) + plot_layout(widths = c(1, 1.3))

  ggsave(file.path(RPT_SUP, "SUPP_panel_E_rrho2_ora.png"), pE_combined,
         width = 400, height = 220, units = "mm", dpi = 300)
}

# --- Summary CSV
rrho2_meta <- tibble(
  quadrant = c("Exacerbated_Up", "Exacerbated_Down",
               "Reversed_AgingUp_TrainingDown", "Reversed_AgingDown_TrainingUp"),
  max_neg_log10_pvalue = round(c(max_UU, max_DD, max_UD, max_DU), 2),
  n_hotspot_genes = c(n_UU, n_DD, n_UD, n_DU),
  n_ora_pathways = c(nrow(ora_UU), nrow(ora_DD), nrow(ora_UD), nrow(ora_DU)),
  matrix_rows = nr, matrix_cols = nc, n_shared_genes = n_shared
)
write_csv(rrho2_meta, file.path(DAT, "panel_E", "rrho2_summary.csv"))

message("F05 Panel E done")
