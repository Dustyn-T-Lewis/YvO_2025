################################################################################
#   Figure 2 — Panel F: Expanded Significance Heatmap + Category-Pathway
#                        Sankey + Stacked Composition Bars
#
#   Generates:
#     b_reports/panel_F_interaction.pdf, panel_F_interaction.png
#     c_data/panel_F/*.csv
################################################################################

if (!exists("dep_df")) source("05_Figures/F2/a_script/YvO_F2_setup.R")

# ==============================================================================
# 1. EXPANDED 7-CATEGORY CLASSIFICATION
# ==============================================================================

expanded_df <- dep_df %>%
  filter(!is.na(logFC_Training_Young), !is.na(logFC_Training_Old)) %>%
  mutate(expanded_cat = classify_expanded(
    pi_score_Training_Young, pi_score_Training_Old, pi_score_Interaction,
    logFC_Training_Young, logFC_Training_Old
  )) %>%
  filter(!is.na(expanded_cat)) %>%
  left_join(imputation_df, by = "gene") %>%
  mutate(imputed = replace_na(imputed, FALSE))

# Soft cap: if >150 proteins, keep top N per category by |logFC|
MAX_GENES <- 150
if (nrow(expanded_df) > MAX_GENES) {
  expanded_df <- expanded_df %>%
    group_by(expanded_cat) %>%
    arrange(desc(abs(logFC_Training_Young) + abs(logFC_Training_Old))) %>%
    slice_head(n = ceiling(MAX_GENES / length(EXPANDED_ORDER))) %>%
    ungroup()
}

# ==============================================================================
# 2. GO SLIM SUPER-CATEGORY MAPPING
# ==============================================================================

suppressPackageStartupMessages(library(GO.db))

expanded_genes <- expanded_df$gene
all_genes      <- dep_df$gene

# --- 2a. GO:BP annotations for all genes (universe + foreground) ---
suppressMessages({
  all_entrez <- AnnotationDbi::mapIds(org.Hs.eg.db, keys = all_genes,
                    keytype = "SYMBOL", column = "ENTREZID", multiVals = "first")
  all_go <- AnnotationDbi::select(org.Hs.eg.db,
               keys = as.character(na.omit(all_entrez)),
               keytype = "ENTREZID", columns = c("SYMBOL", "GO", "ONTOLOGY"))
})
all_bp <- all_go %>% filter(ONTOLOGY == "BP", !is.na(GO)) %>% distinct(SYMBOL, GO)
fg_bp  <- all_bp %>% filter(SYMBOL %in% expanded_genes)

# --- 2b. GO Slim Generic BP terms (62 curated from goslim_generic.obo) ---
bp_slim <- c(
  "GO:0000278", "GO:0000910", "GO:0002181", "GO:0002376", "GO:0003012",
  "GO:0003013", "GO:0003014", "GO:0003016", "GO:0005975", "GO:0006091",
  "GO:0006260", "GO:0006281", "GO:0006310", "GO:0006325", "GO:0006351",
  "GO:0006355", "GO:0006399", "GO:0006457", "GO:0006520", "GO:0006629",
  "GO:0006766", "GO:0006886", "GO:0006913", "GO:0006914", "GO:0006954",
  "GO:0007005", "GO:0007010", "GO:0007018", "GO:0007031", "GO:0007059",
  "GO:0007126", "GO:0007155", "GO:0007163", "GO:0007586", "GO:0009100",
  "GO:0012501", "GO:0016071", "GO:0016192", "GO:0023052", "GO:0030154",
  "GO:0030163", "GO:0030198", "GO:0032200", "GO:0034330", "GO:0042060",
  "GO:0042180", "GO:0042254", "GO:0044782", "GO:0048856", "GO:0048870",
  "GO:0050877", "GO:0051604", "GO:0055085", "GO:0055086", "GO:0061024",
  "GO:0065003", "GO:0071941", "GO:0072659", "GO:0098542", "GO:0098754",
  "GO:0140014", "GO:1901135")

# 13 super-categories: slim GO ID -> biologically coherent group
SLIM_SUPER <- c(
  "GO:0003012" = "Muscle & Contractile", "GO:0003013" = "Circulatory System",
  "GO:0030198" = "ECM & Adhesion", "GO:0007155" = "ECM & Adhesion",
  "GO:0034330" = "ECM & Adhesion", "GO:0042060" = "ECM & Adhesion",
  "GO:0007010" = "Cytoskeleton & Motility", "GO:0048870" = "Cytoskeleton & Motility",
  "GO:0007018" = "Cytoskeleton & Motility", "GO:0007163" = "Cytoskeleton & Motility",
  "GO:0044782" = "Cytoskeleton & Motility",
  "GO:0023052" = "Signaling", "GO:0050877" = "Signaling",
  "GO:0002376" = "Immune & Inflammation", "GO:0006954" = "Immune & Inflammation",
  "GO:0098542" = "Immune & Inflammation",
  "GO:0006520" = "Metabolism", "GO:0005975" = "Metabolism",
  "GO:0006629" = "Metabolism", "GO:0006091" = "Metabolism",
  "GO:0042180" = "Metabolism", "GO:0055086" = "Metabolism",
  "GO:0006766" = "Metabolism", "GO:0071941" = "Metabolism",
  "GO:0098754" = "Metabolism", "GO:0007586" = "Metabolism",
  "GO:1901135" = "Metabolism",
  "GO:0007005" = "Mitochondria & Energy", "GO:0007031" = "Mitochondria & Energy",
  "GO:0006457" = "Protein Homeostasis", "GO:0030163" = "Protein Homeostasis",
  "GO:0006914" = "Protein Homeostasis", "GO:0051604" = "Protein Homeostasis",
  "GO:0065003" = "Protein Homeostasis", "GO:0009100" = "Protein Homeostasis",
  "GO:0055085" = "Transport", "GO:0016192" = "Transport",
  "GO:0006886" = "Transport", "GO:0006913" = "Transport",
  "GO:0072659" = "Transport", "GO:0061024" = "Transport",
  "GO:0006351" = "Gene Expression", "GO:0006355" = "Gene Expression",
  "GO:0016071" = "Gene Expression", "GO:0006399" = "Gene Expression",
  "GO:0002181" = "Gene Expression", "GO:0042254" = "Gene Expression",
  "GO:0006325" = "Gene Expression",
  "GO:0006281" = "DNA & Cell Cycle", "GO:0006260" = "DNA & Cell Cycle",
  "GO:0006310" = "DNA & Cell Cycle", "GO:0032200" = "DNA & Cell Cycle",
  "GO:0000278" = "DNA & Cell Cycle", "GO:0140014" = "DNA & Cell Cycle",
  "GO:0007059" = "DNA & Cell Cycle", "GO:0000910" = "DNA & Cell Cycle",
  "GO:0007126" = "DNA & Cell Cycle",
  "GO:0048856" = "Development", "GO:0030154" = "Development",
  "GO:0012501" = "Development", "GO:0003014" = "Development",
  "GO:0003016" = "Development")

# --- 2c. Map GO terms to slim terms via ancestor lookup ---
ancestors  <- as.list(GOBPANCESTOR)
all_go_ids <- unique(all_bp$GO)

go_to_slim <- setNames(
  lapply(all_go_ids, function(go_id) {
    hits <- character(0)
    if (go_id %in% bp_slim) hits <- go_id
    anc <- ancestors[[go_id]]
    if (!is.null(anc)) hits <- c(hits, intersect(anc, bp_slim))
    unique(hits)
  }),
  all_go_ids)

# Map ALL genes to slim -> super-category (for Fisher background)
all_gene_slim <- all_bp %>%
  mutate(slim_list = go_to_slim[GO]) %>%
  unnest(slim_list) %>%
  select(SYMBOL, slim = slim_list) %>%
  distinct()

all_gene_super <- all_gene_slim %>%
  mutate(super = SLIM_SUPER[slim]) %>%
  filter(!is.na(super)) %>%
  distinct(SYMBOL, super)

# --- 2d. Specificity-weighted 1:1 assignment for foreground genes ---
fg_gene_slim <- all_gene_slim %>% filter(SYMBOL %in% expanded_genes)
fg_gene_super <- fg_gene_slim %>%
  mutate(super = SLIM_SUPER[slim]) %>%
  filter(!is.na(super))

fg_term_counts <- fg_gene_slim %>% count(slim, name = "n_fg")

best_super <- fg_gene_super %>%
  left_join(fg_term_counts, by = "slim") %>%
  mutate(priority = ifelse(super == "Development", 2, 1)) %>%
  arrange(priority, n_fg) %>%
  group_by(SYMBOL) %>%
  slice_head(n = 1) %>%
  ungroup()

# Unmapped genes -> "Other"
unmapped <- setdiff(expanded_genes, best_super$SYMBOL)
if (length(unmapped) > 0) {
  best_super <- bind_rows(best_super,
    tibble(SYMBOL = unmapped, slim = "OTHER", super = "Other",
           n_fg = NA_integer_, priority = 3L))
}

# Merge categories with <2 genes into "Other"
small_cats <- best_super %>% count(super) %>%
  filter(n < 2, super != "Other") %>% pull(super)
if (length(small_cats) > 0) {
  best_super <- best_super %>%
    mutate(super = ifelse(super %in% small_cats, "Other", super))
}

# --- 2e. Build links_1to1 and Fisher enrichment ---
links_1to1 <- best_super %>%
  select(gene = SYMBOL, pathway = super) %>%
  left_join(expanded_df %>% select(gene, expanded_cat,
                                    logFC_Training_Young, logFC_Training_Old),
            by = "gene") %>%
  mutate(database = ifelse(pathway == "Other", "Other", "GO_Slim"))

# Fisher's exact test per super-category
active_supers <- unique(links_1to1$pathway[links_1to1$pathway != "Other"])
n_universe    <- length(all_genes)

fisher_results <- tibble(
  pathway_label = active_supers,
  pvalue = sapply(active_supers, function(s) {
    fg_in  <- sum(links_1to1$pathway == s)
    fg_out <- length(expanded_genes) - fg_in
    bg_in  <- n_distinct(all_gene_super$SYMBOL[all_gene_super$super == s])
    bg_out <- n_universe - bg_in
    fisher.test(matrix(c(fg_in, bg_in, fg_out, bg_out), 2, 2),
                alternative = "greater")$p.value
  })
) %>%
  mutate(p.adjust = p.adjust(pvalue, method = "BH"),
         ID = pathway_label, database = "GO_Slim")

if ("Other" %in% links_1to1$pathway) {
  fisher_results <- bind_rows(fisher_results, tibble(
    pathway_label = "Other", pvalue = 1, p.adjust = 1,
    ID = "OTHER", database = "Other"))
}
ora_top <- fisher_results

pw_pvals_map <- setNames(ora_top$pvalue, ora_top$pathway_label)
links_1to1 <- links_1to1 %>% mutate(pvalue = pw_pvals_map[pathway])

# ==============================================================================
# 3. ORDERING
# ==============================================================================

# Pathway ordering: least significant first, Other always last
pw_1to1_counts <- links_1to1 %>% count(pathway, name = "n_1to1")
pw_order_df <- pw_1to1_counts %>%
  left_join(ora_top %>% select(pathway_label, pvalue),
            by = c("pathway" = "pathway_label")) %>%
  arrange(desc(pvalue))
pw_order <- pw_order_df$pathway
if ("Other" %in% pw_order) pw_order <- c(setdiff(pw_order, "Other"), "Other")

# Gene ordering: by expanded_cat -> pathway -> logFC within group
gene_order <- links_1to1 %>%
  mutate(cat_idx = match(as.character(expanded_cat), EXPANDED_ORDER),
         pw_idx  = match(pathway, pw_order)) %>%
  arrange(cat_idx, pw_idx, logFC_Training_Young) %>%
  pull(gene)

# ==============================================================================
# 4. UNIFIED COORDINATE SYSTEM: HEATMAP + ANNOTATION + SANKEY + ENRICHMENT
# ==============================================================================

n_genes <- length(gene_order)
Y_SPAN  <- n_genes

# --- Gene row geometry (top-down, tiles touching within groups) ---
n_cat_breaks <- length(unique(links_1to1$expanded_cat)) - 1
gap_size     <- 1.2
total_gap    <- n_cat_breaks * gap_size
gene_h       <- (Y_SPAN - total_gap) / n_genes

gene_bar_df <- tibble(gene = gene_order, idx = seq_along(gene_order)) %>%
  mutate(ymax  = Y_SPAN - (idx - 1) * gene_h,
         ymin  = ymax - gene_h,
         y_ctr = (ymin + ymax) / 2) %>%
  left_join(links_1to1 %>% select(gene, expanded_cat) %>% distinct(), by = "gene")

# Category gaps between groups
cat_breaks <- gene_bar_df %>%
  mutate(cat_str = as.character(expanded_cat)) %>%
  mutate(is_break = cat_str != lag(cat_str, default = cat_str[1]))
cum_gap <- cumsum(cat_breaks$is_break) * gap_size
gene_bar_df <- gene_bar_df %>%
  mutate(ymax = ymax - cum_gap, ymin = ymin - cum_gap,
         y_ctr = (ymin + ymax) / 2)

# --- X-coordinate layout constants ---
HM_X_Y   <- 1        # Young column center
HM_X_O   <- 1.92     # Old column center
HM_HW    <- 0.45     # heatmap tile half-width
ANNO_X   <- HM_X_O + HM_HW + 0.01 + 0.12   # annotation strip center
ANNO_HW  <- 0.12     # annotation strip half-width
S_X_PW   <- 5.2      # pathway bar center x
S_PW_HW  <- 0.12     # pathway bar half-width
S_X_BAR  <- S_X_PW + S_PW_HW + 0.10   # enrichment bars start x
S_MAX_LEN <- 2.8     # max enrichment bar length

# --- Heatmap tile data ---
heatmap_df <- expanded_df %>%
  filter(gene %in% gene_order) %>%
  left_join(gene_bar_df %>% select(gene, y_ctr, ymin, ymax), by = "gene") %>%
  select(gene, y_ctr, ymin, ymax, logFC_Training_Young, logFC_Training_Old) %>%
  pivot_longer(cols = c(logFC_Training_Young, logFC_Training_Old),
               names_to = "contrast", values_to = "logFC") %>%
  mutate(x = ifelse(contrast == "logFC_Training_Young", HM_X_Y, HM_X_O))

fc_max <- max(abs(heatmap_df$logFC), na.rm = TRUE)

# --- Annotation strip (one rect per gene, colored by category) ---
anno_rects <- gene_bar_df %>%
  mutate(fill = EXPANDED_COLORS[as.character(expanded_cat)])

# --- Category group positions ---
cat_heatmap_pos <- gene_bar_df %>%
  group_by(expanded_cat) %>%
  summarise(y_bot = min(ymin), y_top = max(ymax), .groups = "drop") %>%
  mutate(y_ctr = (y_top + y_bot) / 2, bar_h = y_top - y_bot)

# Group outline rectangles
group_outlines <- cat_heatmap_pos %>%
  mutate(xmin = HM_X_Y - HM_HW, xmax = ANNO_X + ANNO_HW,
         outline_col = EXPANDED_COLORS[as.character(expanded_cat)])

# Shared y-limits
S_Y_MIN  <- min(gene_bar_df$ymin) - gene_h * 2
S_Y_MAX  <- max(gene_bar_df$ymax) + gene_h * 2
S_Y_SPAN <- S_Y_MAX - S_Y_MIN

# --- Sankey data ---
cat_pw_counts <- links_1to1 %>%
  filter(pathway != "Other") %>%
  count(expanded_cat, pathway, name = "n_proteins") %>%
  filter(n_proteins > 0)

cat_totals <- links_1to1 %>%
  filter(pathway != "Other") %>%
  count(expanded_cat, name = "total") %>%
  arrange(match(as.character(expanded_cat), EXPANDED_ORDER))

pw_totals <- links_1to1 %>%
  filter(pathway != "Other") %>%
  count(pathway, name = "total") %>%
  mutate(pw_idx = match(pathway, pw_order)) %>%
  filter(!is.na(pw_idx)) %>%
  arrange(pw_idx)

n_cats <- nrow(cat_totals)
n_pws  <- nrow(pw_totals)

# Category bars aligned to heatmap gene groups
cat_bars <- cat_totals %>%
  left_join(cat_heatmap_pos, by = "expanded_cat") %>%
  mutate(fill = EXPANDED_COLORS[as.character(expanded_cat)])

# Pathway bars (distributed within y-range)
pw_gap_frac <- 0.03
pw_usable   <- S_Y_SPAN * (1 - pw_gap_frac * max(n_pws - 1, 0) / max(n_pws, 1))
pw_gap_size <- if (n_pws > 1) (S_Y_SPAN - pw_usable) / (n_pws - 1) else 0

pw_bars <- pw_totals %>%
  mutate(bar_h = pw_usable / n_pws,
         y_top = S_Y_MAX - cumsum(c(0, head(bar_h + pw_gap_size, -1))),
         y_bot = y_top - bar_h,
         y_ctr = (y_top + y_bot) / 2)

# Pathway colors (pastel palette)
active_pw_order <- pw_totals$pathway
pw_base <- c("#F48FB1", "#FDAE91", "#E8E8A0", "#A8D8A8", "#8DD3C7",
             "#A2CEE5", "#B6C8E8", "#DEB4D4", "#C9A9A6", "#AED581",
             "#CE93D8", "#80DEEA", "#FFD54F", "#90A4AE", "#E57373",
             "#FFCC80", "#9FA8DA", "#B39DDB", "#80CBC4", "#81C784")
pw_palette <- if (length(active_pw_order) > length(pw_base)) {
  colorRampPalette(pw_base)(length(active_pw_order))
} else {
  pw_base[seq_along(active_pw_order)]
}
PW_COLORS <- setNames(pw_palette, active_pw_order)
pw_bars$fill <- PW_COLORS[pw_bars$pathway]

# Enrichment bar scaling (fall back to gene counts when p-values don't discriminate)
pw_enrichment <- ora_top %>%
  filter(pathway_label != "Other") %>%
  select(pathway_label, p.adjust) %>%
  mutate(neg_log10_padj = -log10(pmax(p.adjust, 1e-50)))

if (max(pw_enrichment$neg_log10_padj, na.rm = TRUE) < 0.01) {
  pw_gene_counts <- links_1to1 %>%
    filter(pathway != "Other") %>% count(pathway, name = "n_bar")
  pw_enrichment <- pw_enrichment %>%
    left_join(pw_gene_counts, by = c("pathway_label" = "pathway")) %>%
    mutate(neg_log10_padj = replace_na(n_bar, 0)) %>%
    select(-n_bar)
  ENRICH_AXIS_LABEL <- "Gene count"
} else {
  ENRICH_AXIS_LABEL <- expression(-log[10](p[adj]))
}

pw_bars <- pw_bars %>%
  left_join(pw_enrichment, by = c("pathway" = "pathway_label")) %>%
  mutate(neg_log10_padj = replace_na(neg_log10_padj, 0),
         pw_label = stringr::str_wrap(pathway, width = 18))

pw_label_size <- if (n_pws > 20) 4.8 else if (n_pws > 15) 5.5 else 6.0

# --- Sigmoid ribbons (from annotation strip to pathway bars) ---
make_sigmoid_ribbon <- function(x0, x1, y0_top, y0_bot, y1_top, y1_bot,
                                n_pts = 50, ribbon_id) {
  t <- seq(0, 1, length.out = n_pts)
  blend <- (1 - cos(pi * t)) / 2
  tibble(
    x = c(x0 + (x1 - x0) * t, rev(x0 + (x1 - x0) * t)),
    y = c(y0_top + (y1_top - y0_top) * blend,
          rev(y0_bot + (y1_bot - y0_bot) * blend)),
    ribbon_id = ribbon_id)
}

cat_cum <- setNames(rep(0, n_cats), as.character(cat_totals$expanded_cat))
pw_cum  <- setNames(rep(0, n_pws), pw_totals$pathway)
ribbon_list <- list()
ribbon_idx  <- 0

for (ct in as.character(cat_totals$expanded_cat)) {
  ct_contribs <- cat_pw_counts %>%
    filter(as.character(expanded_cat) == ct) %>%
    arrange(match(pathway, active_pw_order))
  for (r in seq_len(nrow(ct_contribs))) {
    pw <- ct_contribs$pathway[r]; n <- ct_contribs$n_proteins[r]
    if (n == 0) next
    ribbon_idx <- ribbon_idx + 1
    ct_row <- cat_bars %>% filter(as.character(expanded_cat) == ct)
    frac_ct <- n / ct_row$total
    y0_top <- ct_row$y_top - cat_cum[ct] * ct_row$bar_h
    y0_bot <- y0_top - frac_ct * ct_row$bar_h
    cat_cum[ct] <- cat_cum[ct] + frac_ct
    pw_row <- pw_bars %>% filter(pathway == pw)
    frac_pw <- n / pw_row$total
    y1_top <- pw_row$y_top - pw_cum[pw] * pw_row$bar_h
    y1_bot <- y1_top - frac_pw * pw_row$bar_h
    pw_cum[pw] <- pw_cum[pw] + frac_pw
    ribbon_list[[ribbon_idx]] <- make_sigmoid_ribbon(
      x0 = ANNO_X + ANNO_HW, x1 = S_X_PW - S_PW_HW,
      y0_top = y0_top, y0_bot = y0_bot,
      y1_top = y1_top, y1_bot = y1_bot,
      ribbon_id = paste0(ct, "->", pw)
    ) %>% mutate(expanded_cat = ct, pathway = pw, fill_col = EXPANDED_COLORS[ct])
  }
}
ribbons_df <- bind_rows(ribbon_list)

# --- Stacked enrichment bars (category-colored segments, length proportional to significance) ---
max_enrich <- max(pw_bars$neg_log10_padj, na.rm = TRUE)
S_SBAR_H   <- min(5.0, 0.90 * S_Y_SPAN / max(n_pws, 1))

stacked_rects <- list()
stacked_idx   <- 0
for (pw in active_pw_order) {
  pw_row <- pw_bars %>% filter(pathway == pw)
  if (is.na(pw_row$neg_log10_padj) || pw_row$neg_log10_padj == 0) next
  pw_contribs <- cat_pw_counts %>%
    filter(pathway == pw) %>%
    arrange(match(as.character(expanded_cat), EXPANDED_ORDER))
  bar_y_ctr <- pw_row$y_ctr
  bar_total_w <- (pw_row$neg_log10_padj / max_enrich) * S_MAX_LEN
  x_cursor <- S_X_BAR
  for (r in seq_len(nrow(pw_contribs))) {
    ct <- as.character(pw_contribs$expanded_cat[r])
    n  <- pw_contribs$n_proteins[r]
    if (n == 0) next
    stacked_idx <- stacked_idx + 1
    seg_w <- (n / pw_row$total) * bar_total_w
    stacked_rects[[stacked_idx]] <- tibble(
      xmin = x_cursor, xmax = x_cursor + seg_w,
      ymin = bar_y_ctr - S_SBAR_H / 2, ymax = bar_y_ctr + S_SBAR_H / 2,
      expanded_cat = ct, pathway = pw, fill = EXPANDED_COLORS[ct])
    x_cursor <- x_cursor + seg_w
  }
}
stacked_df <- bind_rows(stacked_rects)

# Count labels at bar ends
count_labels <- pw_bars %>%
  mutate(bar_end = S_X_BAR + (neg_log10_padj / max_enrich) * S_MAX_LEN,
         label = paste0("n=", total),
         x = pmax(bar_end, S_X_BAR) + 0.10)

# --- Category key (overlaid in enrichment area) ---
key_cats <- names(EXPANDED_COLORS)
key_x0   <- S_X_BAR + S_MAX_LEN * 0.45
key_y0   <- S_Y_MAX
key_df <- tibble(cat = key_cats, color = EXPANDED_COLORS[key_cats],
                 row = seq_along(key_cats)) %>%
  mutate(sx = key_x0, sy = key_y0 - (row - 1) * 2.0)

# Enrichment axis
enrich_ticks <- pretty(c(0, max_enrich), n = 4)
enrich_ticks <- enrich_ticks[enrich_ticks >= 0 & enrich_ticks <= max_enrich * 1.05]
if (!0 %in% enrich_ticks) enrich_ticks <- c(0, enrich_ticks)
grid_x  <- S_X_BAR + (enrich_ticks / max_enrich) * S_MAX_LEN
grid_df <- tibble(x = grid_x, label = as.character(round(enrich_ticks, 1)))

axis_y   <- S_Y_MIN + 1
tick_len <- 1.5

# Dynamic text sizes
gene_text_pt <- if (n_genes > 100) 7 else if (n_genes > 60) 8.5 else 10
header_y     <- S_Y_MAX + 1
sig_header_x <- S_X_BAR + S_MAX_LEN / 2

# ==============================================================================
# 5. BUILD UNIFIED PANEL F PLOT
# ==============================================================================

pF <- ggplot() +
  # Heatmap tiles (diverging logFC)
  geom_rect(data = heatmap_df,
            aes(xmin = x - HM_HW, xmax = x + HM_HW,
                ymin = ymin, ymax = ymax, fill = logFC),
            color = NA, linewidth = 0) +
  scale_fill_gradient2(
    low = "#2166AC", mid = "white", high = "#B2182B",
    midpoint = 0, limits = c(-fc_max, fc_max),
    name = expression(log[2]~FC),
    guide = guide_colorbar(direction = "horizontal", title.position = "top",
                           barwidth = unit(50, "mm"), barheight = unit(5, "mm"),
                           title.theme = element_text(size = 10, face = "bold"),
                           label.theme = element_text(size = 9))
  ) +
  # Switch to identity fill for remaining layers
  ggnewscale::new_scale_fill() +
  # Annotation strip
  geom_rect(data = anno_rects,
            aes(xmin = ANNO_X - ANNO_HW, xmax = ANNO_X + ANNO_HW,
                ymin = ymin, ymax = ymax, fill = fill),
            color = NA, linewidth = 0) +
  # Group outlines
  geom_rect(data = group_outlines,
            aes(xmin = xmin, xmax = xmax, ymin = y_bot, ymax = y_top),
            fill = NA, color = "black", linewidth = 0.5) +
  # Sigmoid Sankey ribbons
  geom_polygon(data = ribbons_df,
               aes(x = x, y = y, group = ribbon_id, fill = fill_col),
               alpha = 0.30, color = NA) +
  # Pathway bars
  geom_rect(data = pw_bars,
            aes(xmin = S_X_PW - S_PW_HW, xmax = S_X_PW + S_PW_HW,
                ymin = y_bot, ymax = y_top, fill = fill),
            color = "black", linewidth = 0.3) +
  # Pathway labels
  geom_text(data = pw_bars,
            aes(x = S_X_PW - S_PW_HW - 0.05, y = y_ctr, label = pw_label),
            hjust = 1, size = pw_label_size, fontface = "bold",
            color = "grey20", lineheight = 0.85) +
  # Enrichment bars (stacked, category-colored)
  geom_rect(data = stacked_df,
            aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, fill = fill),
            color = "black", linewidth = 0.2) +
  # Count labels
  geom_text(data = count_labels,
            aes(x = x, y = y_ctr, label = label),
            hjust = 0, size = 4.0, fontface = "bold", color = "grey30") +
  # Category key swatches
  geom_point(data = key_df, aes(x = sx, y = sy, color = color),
             shape = 15, size = pw_label_size, show.legend = FALSE) +
  scale_color_identity() +
  # Category key labels
  geom_text(data = key_df,
            aes(x = sx + 0.15, y = sy, label = cat),
            hjust = 0, size = pw_label_size, fontface = "bold", color = "grey20") +
  scale_fill_identity() +
  # Column headers
  annotate("text", x = HM_X_Y, y = header_y, label = "Young",
           size = 5, fontface = "bold") +
  annotate("text", x = HM_X_O, y = header_y, label = "Old",
           size = 5, fontface = "bold") +
  # Enrichment axis label
  annotate("text", x = sig_header_x, y = axis_y - tick_len - 3.5,
           label = ENRICH_AXIS_LABEL,
           size = 4.5, fontface = "bold", color = "grey20") +
  # Enrichment axis line
  annotate("segment", x = S_X_BAR, xend = S_X_BAR + S_MAX_LEN,
           y = axis_y, yend = axis_y, color = "black", linewidth = 0.5) +
  geom_segment(data = grid_df, aes(x = x, xend = x),
               y = axis_y, yend = axis_y - tick_len,
               color = "black", linewidth = 0.4) +
  geom_text(data = grid_df, aes(x = x, label = label),
            y = axis_y - tick_len - 1.5, size = 4.0, color = "grey30") +
  # Scales and theme
  scale_y_continuous(
    breaks = gene_bar_df$y_ctr, labels = gene_bar_df$gene,
    limits = c(axis_y - tick_len - 5, S_Y_MAX + 3),
    expand = c(0, 0)) +
  scale_x_continuous(
    limits = c(HM_X_Y - HM_HW - 0.05, S_X_BAR + S_MAX_LEN + 1.5),
    expand = c(0, 0)) +
  labs(x = NULL, y = NULL,
       title = "Significant Proteins: Multi-Category Heatmap & Pathway Enrichment",
       subtitle = sprintf("%d proteins across %d significance categories | GO Slim super-categories (Fisher, BH adj.)",
                          nrow(expanded_df), length(EXPANDED_ORDER))) +
  theme_minimal(base_size = 9) +
  theme(
    plot.title       = element_text(face = "bold", size = 12),
    plot.subtitle    = element_text(size = 9, color = "grey30", face = "italic"),
    axis.text.y      = element_text(size = gene_text_pt, face = "bold.italic", hjust = 1),
    axis.text.x      = element_blank(),
    axis.ticks.x     = element_blank(),
    axis.ticks.y     = element_blank(),
    panel.grid       = element_blank(),
    legend.position  = "inside",
    legend.position.inside = c(0.10, 0.01),
    legend.justification = c(0.5, 0),
    legend.text      = element_text(size = 9),
    legend.margin    = margin(t = 8),
    plot.margin      = margin(8, 5, 8, 2))

# ==============================================================================
# 6. EXPORT CSVs
# ==============================================================================

dir.create(file.path(DAT_DIR, "panel_F"), recursive = TRUE, showWarnings = FALSE)

expanded_df %>%
  transmute(gene, expanded_category = as.character(expanded_cat),
            logFC_Training_Young = round(logFC_Training_Young, 4),
            logFC_Training_Old   = round(logFC_Training_Old, 4),
            pi_score_Young       = round(pi_score_Training_Young, 6),
            pi_score_Old         = round(pi_score_Training_Old, 6),
            pi_score_Interaction = round(pi_score_Interaction, 6),
            imputed) %>%
  arrange(expanded_category, desc(abs(logFC_Training_Young) + abs(logFC_Training_Old))) %>%
  write_csv(file.path(DAT_DIR, "panel_F", "expanded_classification.csv"))

links_1to1 %>%
  transmute(gene, pathway, database, expanded_category = as.character(expanded_cat),
            logFC_Training_Young = round(logFC_Training_Young, 4),
            logFC_Training_Old = round(logFC_Training_Old, 4),
            pathway_pvalue = signif(pvalue, 4)) %>%
  arrange(expanded_category, pathway, gene) %>%
  write_csv(file.path(DAT_DIR, "panel_F", "expanded_sankey_links.csv"))

stacked_df %>%
  transmute(expanded_category = expanded_cat, pathway,
            xmin = round(xmin, 4), xmax = round(xmax, 4)) %>%
  write_csv(file.path(DAT_DIR, "panel_F", "stacked_bar_data.csv"))

# Legacy export (backward compat with supplementary workbook)
links_1to1 %>%
  transmute(gene, pathway, database,
            interaction_pattern = as.character(expanded_cat),
            logFC_Training_Young = round(logFC_Training_Young, 4),
            logFC_Training_Old = round(logFC_Training_Old, 4),
            pathway_pvalue = signif(pvalue, 4)) %>%
  arrange(interaction_pattern, pathway, gene) %>%
  write_csv(file.path(DAT_DIR, "panel_F", "sankey_links.csv"))

# ==============================================================================
# 7. SAVE PANEL F
# ==============================================================================

f_height <- max(220, nrow(links_1to1) * 3.8 + 80)

ggsave(file.path(RPT_DIR, "panel_F_interaction.pdf"), pF,
       width = 500, height = f_height, units = "mm", device = pdf, limitsize = FALSE)
ggsave(file.path(RPT_DIR, "panel_F_interaction.png"), pF,
       width = 500, height = f_height, units = "mm", dpi = 300, limitsize = FALSE)

message("Panel F complete: ", file.path(RPT_DIR, "panel_F_interaction.pdf"))
