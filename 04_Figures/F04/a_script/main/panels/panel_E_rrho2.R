# F04 Panel E: RRHO2 Concordance Map + Per-Quadrant ORA
setwd(rprojroot::find_rstudio_root_file())
source("04_Figures/shared/style.R")
source("04_Figures/shared/print_scale_380.R")
source("04_Figures/shared/pathway_utils.R")
library(tidyverse)
library(ggrepel)
library(msigdbr)
library(fgsea)
library(RRHO2)

PE_W <- 89

RPT_PNG <- "04_Figures/F04/b_reports/main/png/panels"
RPT_PDF <- "04_Figures/F04/b_reports/main/pdf/panels"
DAT <- "04_Figures/F04/c_data"
dir.create(RPT_PNG, recursive = TRUE, showWarnings = FALSE)
dir.create(RPT_PDF, recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(DAT, "panel_E"), recursive = TRUE, showWarnings = FALSE)

pdf_device <- get_pdf_device()

dep_df <- read_csv("03_DEP/c_data/03_combined_results.csv", show_col_types = FALSE)

rr_df <- dep_df %>%
  transmute(gene, t_young = t_Training_Young, t_old = t_Training_Old) %>%
  filter(!is.na(t_young) & !is.na(t_old)) %>%
  distinct(gene, .keep_all = TRUE)

n_shared <- nrow(rr_df)

# --- RRHO2 computation (stratified hypergeometric, Cahill et al. 2018)
list1 <- data.frame(gene = rr_df$gene, score = rr_df$t_young, stringsAsFactors = FALSE)
list2 <- data.frame(gene = rr_df$gene, score = rr_df$t_old,   stringsAsFactors = FALSE)

rrho_obj <- RRHO2_initialize(
  list1, list2,
  labels          = c("Training (Young)", "Training (Old)"),
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
# [row_before, col_before] = top of both lists = Concordant Up (UU)
# [row_after,  col_after]  = bottom of both    = Concordant Down (DD)
# [row_before, col_after]  = top list1, bottom list2 = Discordant (Y Up / O Down)
# [row_after,  col_before] = bottom list1, top list2 = Discordant (Y Down / O Up)
max_UU <- max(hmat[row_before, col_before], na.rm = TRUE)
max_DD <- max(hmat[row_after,  col_after],  na.rm = TRUE)
max_UD <- max(hmat[row_before, col_after],  na.rm = TRUE)
max_DU <- max(hmat[row_after,  col_before], na.rm = TRUE)

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
message(sprintf("  Hotspot genes: UU=%d, DD=%d, UD=%d, DU=%d", n_UU, n_DD, n_UD, n_DU))

# Match panel A/B/C quadrant label size
txt_quad <- scale_text(BASE_QUADRANT, 146) * 1.15   # +1pt for RRHO2 quadrant readability

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
LABEL_PADDING <- unit(0.8 * PRINT_SCALE, "mm")

pE_heat <- ggplot(hmat_df, aes(x = row, y = col, fill = neg_log10_p)) +
  geom_raster() +
  scale_fill_gradientn(
    colors   = JET_COLORS,
    limits   = c(0, max_val),
    na.value = "white",
    name     = expression(-log[10](P)),
    guide    = guide_colorbar(
      barwidth = unit(15 * PRINT_SCALE, "mm"), barheight = unit(2 * PRINT_SCALE, "mm"),
      title.position = "left", title.vjust = 0.5,
      title.theme = element_text(size = FIG_LEGEND_TITLE, face = "bold", color = "grey15"))
  ) +
  annotate("label", x = ann_x_left, y = ann_y_bot,
           label = "Concordant Up",
           color = "grey15", fill = LABEL_FILL, label.size = 0,
           label.padding = LABEL_PADDING, label.r = unit(0.5, "mm"),
           fontface = "bold", size = txt_quad,
           hjust = 0, vjust = 0) +
  annotate("label", x = ann_x_right, y = ann_y_top,
           label = "Concordant Down",
           color = "grey15", fill = LABEL_FILL, label.size = 0,
           label.padding = LABEL_PADDING, label.r = unit(0.5, "mm"),
           fontface = "bold", size = txt_quad,
           hjust = 1, vjust = 1) +
  annotate("label", x = ann_x_left, y = ann_y_top,
           label = "Discordant (Y\u2191 O\u2193)",
           color = "grey15", fill = LABEL_FILL, label.size = 0,
           label.padding = LABEL_PADDING, label.r = unit(0.5, "mm"),
           fontface = "bold", size = txt_quad,
           hjust = 0, vjust = 1) +
  annotate("label", x = ann_x_right, y = ann_y_bot,
           label = "Discordant (Y\u2193 O\u2191)",
           color = "grey15", fill = LABEL_FILL, label.size = 0,
           label.padding = LABEL_PADDING, label.r = unit(0.5, "mm"),
           fontface = "bold", size = txt_quad,
           hjust = 1, vjust = 0) +
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0)) +
  labs(
    title    = "Threshold-Free Concordance (RRHO2)",
    subtitle = sprintf("Stratified hypergeometric | %d shared genes | warm corners = concordant gene regulation | No MTC (Cahill et al. 2018)", n_shared),
    x = expression("Training (Young) rank"~(Up %->% Down)),
    y = expression("Training (Old) rank"~(Up %->% Down))
  ) +
  FIG_THEME +
  theme(
    axis.text        = element_blank(),
    axis.title.x     = element_text(size = FIG_AXIS_TEXT, face = "bold", margin = margin(t = 2)),
    axis.title.y     = element_text(size = FIG_AXIS_TEXT, face = "bold", margin = margin(r = 2)),
    axis.ticks       = element_blank(),
    panel.border     = element_blank(),
    panel.grid.major = element_blank(),
    legend.position  = "bottom",
    legend.text      = element_text(size = FIG_LEGEND_TEXT, face = "bold"),
    legend.margin    = margin(2, 24, 0, 0, "mm"),
    plot.margin = margin(0, 0, 0, 0, "mm")
  ) +
  coord_fixed(ratio = 1)

hotspot_export <- bind_rows(
  tibble(quadrant = "UU", gene = hotspot_genes$UU),
  tibble(quadrant = "DD", gene = hotspot_genes$DD),
  tibble(quadrant = "UD", gene = hotspot_genes$UD),
  tibble(quadrant = "DU", gene = hotspot_genes$DU)
)
write_csv(hotspot_export, file.path(DAT, "panel_E", "rrho2_hotspot_genes.csv"))

pw_collection_E <- build_pathway_collection(min_size = 15, max_size = 500,
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
      min_size       = 15,
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

ora_UU <- run_quadrant_ora(hotspot_genes$UU, "Concordant Up")
ora_DD <- run_quadrant_ora(hotspot_genes$DD, "Concordant Down")
ora_UD <- run_quadrant_ora(hotspot_genes$UD, "Discordant (Y Up / O Down)")
ora_DU <- run_quadrant_ora(hotspot_genes$DU, "Discordant (Y Down / O Up)")

ora_concordant  <- bind_rows(ora_UU, ora_DD)
ora_discordant  <- bind_rows(ora_UD, ora_DU)

write_csv(ora_concordant,  file.path(DAT, "panel_E", "rrho2_ora_concordant.csv"))
write_csv(ora_discordant,  file.path(DAT, "panel_E", "rrho2_ora_discordant.csv"))

# Save RRHO heatmap (square, since coord_fixed)
ggsave(file.path(RPT_PNG, "MAIN_panel_E_rrho2.png"), pE_heat,
       width = PE_W, height = PE_W, units = "mm", dpi = 300)
ggsave(file.path(RPT_PDF, "MAIN_panel_E_rrho2.pdf"), pE_heat,
       width = PE_W, height = PE_W, units = "mm", device = pdf_device)

# --- Per-quadrant ORA summary
MAX_PER_QUAD <- 12

conv_keys <- character(0)
conv_path <- file.path(DAT, "panel_F_fry", "convergent_pathways.csv")
if (file.exists(conv_path)) {
  conv_df <- tryCatch(read_csv(conv_path, show_col_types = FALSE), error = function(e) tibble())
  if (!"note" %in% names(conv_df) && nrow(conv_df) > 0)
    conv_keys <- tolower(str_trim(conv_df$pathway))
}

quad_meta <- list(
  list(name = "Concordant Up",              slug = "concordant_up",     data = ora_UU, n_hot = n_UU),
  list(name = "Concordant Down",            slug = "concordant_down",   data = ora_DD, n_hot = n_DD),
  list(name = "Discordant (Y Up / O Down)", slug = "discordant_y_up",  data = ora_UD, n_hot = n_UD),
  list(name = "Discordant (Y Down / O Up)", slug = "discordant_y_down", data = ora_DU, n_hot = n_DU)
)

n_conv_total <- 0L
for (qm in quad_meta) {
  if (nrow(qm$data) == 0) next
  q_df <- qm$data %>%
    mutate(neg_log10_padj = -log10(padj),
           pathway_label  = clean_pathway_name(pathway)) %>%
    arrange(desc(neg_log10_padj)) %>%
    slice_head(n = MAX_PER_QUAD)
  if (length(conv_keys) > 0) {
    q_df <- q_df %>%
      mutate(is_conv = tolower(str_trim(pathway_label)) %in% conv_keys,
             pathway_label = ifelse(is_conv, paste0("\u2605 ", pathway_label), pathway_label))
    n_conv_total <- n_conv_total + sum(q_df$is_conv)
  }
  message(sprintf("  %s: %d ORA pathways (%d hotspot genes)", qm$slug, nrow(q_df), qm$n_hot))
}
if (n_conv_total > 0) message(sprintf("  %d ORA bars marked convergent with fry", n_conv_total))

rrho2_meta <- tibble(
  quadrant = c("Concordant_Up", "Concordant_Down",
               "Discordant_YoungUp_OldDown", "Discordant_YoungDown_OldUp"),
  max_neg_log10_pvalue = round(c(max_UU, max_DD, max_UD, max_DU), 2),
  n_hotspot_genes = c(n_UU, n_DD, n_UD, n_DU),
  n_ora_pathways = c(nrow(ora_UU), nrow(ora_DD), nrow(ora_UD), nrow(ora_DU)),
  matrix_rows = nr, matrix_cols = nc, n_shared_genes = n_shared
)
write_csv(rrho2_meta, file.path(DAT, "panel_E", "rrho2_summary.csv"))

# --- Export for composite ---
pE_title    <- "Threshold-Free Concordance (RRHO2)"
pE_subtitle <- sprintf("Stratified hypergeometric | %d shared genes | warm corners = concordant gene regulation | No MTC (Cahill et al. 2018)", n_shared)
pE_legend   <- NULL
pE_heat     <- strip_for_composite(pE_heat)

message("F04 Panel E done")
