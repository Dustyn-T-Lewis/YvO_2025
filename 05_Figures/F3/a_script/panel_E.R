################################################################################
#   Figure 3 — Panel E: Reversal Classification Heatmap + Category-Pathway
#                        Sankey + Stacked Composition Bars
#
#   4-pattern reversal classification (Reversed/Attenuated/Amplified/Concordant)
#   with Hallmark + GO:BP ORA, sigmoid Sankey ribbons, enrichment bars.
#
#   Layout mirrors F2 Panel F (unified coordinate system).
#
#   Generates:
#     b_reports/panel_E_classification.pdf, panel_E_classification.png
#     c_data/panel_E/classification.csv, sankey_links.csv, enrichment_bars.csv
################################################################################

if (!exists("dep_df")) source("05_Figures/F3/a_script/YvO_F3_setup.R")
suppressPackageStartupMessages(library(ggnewscale))

message("Panel E: reversal classification (heatmap + sankey + enrichment bars)...")

# ==============================================================================
# 1. PROTEIN SELECTION & 4-PATTERN REVERSAL CLASSIFICATION
# ==============================================================================

# Select aging-significant OR reversal-significant proteins
class_df <- dep_df %>%
  filter(!is.na(logFC_Aging), !is.na(logFC_Training_Old)) %>%
  filter(pi_score_Aging < 0.05 | pi_score_Reversal < 0.05) %>%
  left_join(imputation_df, by = "gene") %>%
  mutate(imputed = replace_na(imputed, FALSE))

# Soft cap at 150 genes (top by |logFC_Training_Old|)
MAX_GENES <- 150
if (nrow(class_df) > MAX_GENES) {
  class_df <- class_df %>%
    arrange(desc(abs(logFC_Training_Old))) %>%
    slice_head(n = MAX_GENES)
}

# Apply 4-pattern reversal classification
class_df <- class_df %>%
  mutate(pattern = classify_reversal_pattern(logFC_Aging, logFC_Training_Old,
                                             pi_score_Aging))

PATTERN_ORDER <- c("Reversed", "Attenuated", "Amplified", "Concordant")

n_class <- nrow(class_df)
cat_counts <- class_df %>% count(pattern) %>% deframe()
message(sprintf("  Classification: %d proteins across %d categories",
                n_class, length(cat_counts)))
for (ct in names(cat_counts))
  message(sprintf("    %s: %d", ct, cat_counts[ct]))

# ==============================================================================
# 2. ORA + PATHWAY MAPPING (Hallmark + GO:BP via enricher)
# ==============================================================================

message("  Running pooled ORA (Hallmark + GO:BP)...")

all_class_genes <- class_df$gene
all_genes       <- dep_df$gene

# ---- ORA: Hallmark + GO:BP ----
hallmark_t2g <- msigdbr(species = "Homo sapiens", collection = "H") %>%
  dplyr::select(term = gs_name, gene = gene_symbol) %>% distinct()
gobp_t2g <- msigdbr(species = "Homo sapiens", collection = "C5",
                     subcollection = "GO:BP") %>%
  dplyr::select(term = gs_name, gene = gene_symbol) %>% distinct()

run_ora <- function(t2g, db_name) {
  res <- tryCatch(
    enricher(gene = all_class_genes, universe = all_genes, TERM2GENE = t2g,
             pAdjustMethod = "BH", pvalueCutoff = 1, qvalueCutoff = 1,
             minGSSize = 3, maxGSSize = 500),
    error = function(e) NULL)
  if (!is.null(res) && nrow(as.data.frame(res)) > 0)
    as.data.frame(res) %>% mutate(database = db_name) else tibble()
}

ora_combined <- bind_rows(run_ora(hallmark_t2g, "Hallmark"),
                          run_ora(gobp_t2g, "GO:BP")) %>%
  mutate(pathway_label = clean_pathway_name(ID))

# ---- Select top pathways + greedy coverage ----
ora_sig <- ora_combined %>% filter(pvalue < 0.05) %>% arrange(pvalue)
ora_top <- bind_rows(
  ora_sig %>% filter(database == "Hallmark") %>% slice_head(n = 5),
  ora_sig %>% filter(database == "GO:BP")    %>% slice_head(n = 5)
) %>% group_by(pathway_label) %>% slice_min(pvalue, n = 1, with_ties = FALSE) %>%
  ungroup()
if (nrow(ora_top) < 4) ora_top <- ora_combined %>% arrange(pvalue) %>% slice_head(n = 8)

all_t2g    <- bind_rows(hallmark_t2g, gobp_t2g) %>% distinct()
mapped_now <- ora_top %>% dplyr::select(geneID) %>% separate_rows(geneID, sep = "/") %>%
  filter(geneID %in% all_class_genes) %>% pull(geneID) %>% unique()
orphans    <- setdiff(all_class_genes, mapped_now)
remaining  <- ora_combined %>% filter(!ID %in% ora_top$ID) %>% arrange(pvalue)

# Greedy rescue: add pathways for unmapped genes (cap total at 15)
max_total_pws <- 15
for (i in seq_len(15)) {
  if (length(orphans) == 0 || nrow(remaining) == 0) break
  if (nrow(ora_top) >= max_total_pws) break
  coverage <- remaining %>% rowwise() %>%
    mutate(hits = sum(str_split(geneID, "/")[[1]] %in% orphans)) %>%
    ungroup() %>% filter(hits >= 1) %>% arrange(desc(hits), pvalue)
  if (nrow(coverage) == 0) break
  best <- coverage %>% slice_head(n = 1)
  ora_top   <- bind_rows(ora_top, best)
  remaining <- remaining %>% filter(ID != best$ID)
  orphans   <- setdiff(orphans, str_split(best$geneID, "/")[[1]])
}
ora_top <- ora_top %>% group_by(pathway_label) %>%
  slice_min(pvalue, n = 1, with_ties = FALSE) %>% ungroup()

# ---- Build ALL gene-pathway links ----
sankey_links_all <- ora_top %>%
  dplyr::select(pathway = pathway_label, ID, geneID, pvalue, p.adjust, Count,
                database) %>%
  separate_rows(geneID, sep = "/") %>% rename(gene = geneID) %>%
  filter(gene %in% all_class_genes)

# Rescue orphans via membership in selected pathways
leftover <- setdiff(all_class_genes, unique(sankey_links_all$gene))
if (length(leftover) > 0) {
  rescue <- all_t2g %>% filter(gene %in% leftover, term %in% ora_top$ID) %>%
    left_join(ora_top %>% dplyr::select(ID, pathway = pathway_label, pvalue,
                                         p.adjust, Count, database),
              by = c("term" = "ID")) %>%
    group_by(gene) %>% slice_min(pvalue, n = 1, with_ties = FALSE) %>% ungroup() %>%
    mutate(ID = term) %>%
    dplyr::select(gene, pathway, ID, pvalue, p.adjust, Count, database)
  sankey_links_all <- bind_rows(sankey_links_all, rescue)
  leftover <- setdiff(leftover, rescue$gene)
}

# Final catch-all: orphans that matched no selected pathway -> "Other"
leftover <- setdiff(all_class_genes, unique(sankey_links_all$gene))
if (length(leftover) > 0) {
  other_links <- tibble(
    gene = leftover, pathway = "Other", ID = "OTHER",
    pvalue = 1, p.adjust = 1, Count = length(leftover), database = "Other"
  )
  sankey_links_all <- bind_rows(sankey_links_all, other_links)
  ora_top <- bind_rows(ora_top, tibble(
    ID = "OTHER", Description = "Other", GeneRatio = NA, BgRatio = NA,
    pvalue = 1, p.adjust = 1, qvalue = 1,
    geneID = paste(leftover, collapse = "/"),
    Count = length(leftover), database = "Other", pathway_label = "Other"
  ))
  message(sprintf("  %d orphan genes assigned to 'Other' pathway", length(leftover)))
}

sankey_links_all <- sankey_links_all %>%
  left_join(class_df %>% dplyr::select(gene, pattern, logFC_Aging,
                                        logFC_Training_Old),
            by = "gene")

# ---- Force 1:1 mapping: each gene -> best pathway ----
links_1to1 <- sankey_links_all %>%
  group_by(gene) %>%
  slice_min(pvalue, n = 1, with_ties = FALSE) %>%
  ungroup()

# Redistribution loop: move solo-pathway genes to multi-gene pathways (5 iter)
pw_counts <- links_1to1 %>% count(pathway, name = "n_1to1")
for (iter in seq_len(5)) {
  solo_pws <- pw_counts %>% filter(n_1to1 == 1) %>% pull(pathway)
  if (length(solo_pws) == 0) break
  solo_genes <- links_1to1 %>% filter(pathway %in% solo_pws) %>% pull(gene)
  multi_pws <- pw_counts %>% filter(n_1to1 >= 2) %>% pull(pathway)
  better <- sankey_links_all %>%
    filter(gene %in% solo_genes, pathway %in% multi_pws) %>%
    group_by(gene) %>% slice_min(pvalue, n = 1, with_ties = FALSE) %>% ungroup()
  if (nrow(better) == 0) break
  links_1to1 <- bind_rows(
    links_1to1 %>% filter(!gene %in% better$gene),
    better)
  pw_counts <- links_1to1 %>% count(pathway, name = "n_1to1")
}

# Drop pathways that lost all genes after 1:1 assignment
active_pws <- unique(links_1to1$pathway)
ora_top <- ora_top %>% filter(pathway_label %in% active_pws)

message(sprintf("  %d pathways, %d mapped (1:1)", nrow(ora_top), nrow(links_1to1)))

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

# Gene ordering: by pattern -> pathway (by p-value) -> logFC_Training_Old within group
gene_order <- links_1to1 %>%
  mutate(cat_idx = match(as.character(pattern), PATTERN_ORDER),
         pw_idx  = match(pathway, pw_order)) %>%
  arrange(cat_idx, pw_idx, logFC_Training_Old) %>%
  pull(gene)

# ==============================================================================
# 4. UNIFIED COORDINATE SYSTEM: HEATMAP + ANNOTATION + SANKEY + ENRICHMENT
# ==============================================================================

n_genes <- length(gene_order)
Y_SPAN  <- n_genes

# --- Gene row geometry (top-down, tiles touching within groups) ---
n_cat_breaks <- length(unique(links_1to1$pattern)) - 1
gap_size     <- 1.2
total_gap    <- n_cat_breaks * gap_size
gene_h       <- (Y_SPAN - total_gap) / n_genes

gene_bar_df <- tibble(gene = gene_order, idx = seq_along(gene_order)) %>%
  mutate(ymax  = Y_SPAN - (idx - 1) * gene_h,
         ymin  = ymax - gene_h,
         y_ctr = (ymin + ymax) / 2) %>%
  left_join(links_1to1 %>% select(gene, pattern) %>% distinct(), by = "gene")

# Category gaps between groups
cat_breaks <- gene_bar_df %>%
  mutate(cat_str = as.character(pattern)) %>%
  mutate(is_break = cat_str != lag(cat_str, default = cat_str[1]))
cum_gap <- cumsum(cat_breaks$is_break) * gap_size
gene_bar_df <- gene_bar_df %>%
  mutate(ymax = ymax - cum_gap, ymin = ymin - cum_gap,
         y_ctr = (ymin + ymax) / 2)

# --- X-coordinate layout constants ---
HM_X_Y   <- 1        # Aging column center
HM_X_O   <- 1.92     # Training (Old) column center
HM_HW    <- 0.45     # heatmap tile half-width
ANNO_X   <- HM_X_O + HM_HW + 0.01 + 0.12   # annotation strip center
ANNO_HW  <- 0.12     # annotation strip half-width
S_X_PW   <- 5.2      # pathway bar center x
S_PW_HW  <- 0.12     # pathway bar half-width
S_X_BAR  <- S_X_PW + S_PW_HW + 0.10   # enrichment bars start x
S_MAX_LEN <- 2.8     # max enrichment bar length

# --- Heatmap tile data ---
heatmap_df <- class_df %>%
  filter(gene %in% gene_order) %>%
  left_join(gene_bar_df %>% select(gene, y_ctr, ymin, ymax), by = "gene") %>%
  select(gene, y_ctr, ymin, ymax, logFC_Aging, logFC_Training_Old) %>%
  pivot_longer(cols = c(logFC_Aging, logFC_Training_Old),
               names_to = "contrast", values_to = "logFC") %>%
  mutate(x = ifelse(contrast == "logFC_Aging", HM_X_Y, HM_X_O))

fc_max <- max(abs(heatmap_df$logFC), na.rm = TRUE)

# --- Annotation strip (one rect per gene, colored by pattern) ---
anno_rects <- gene_bar_df %>%
  mutate(fill = PATTERN_COLORS[as.character(pattern)])

# --- Category group positions ---
cat_heatmap_pos <- gene_bar_df %>%
  group_by(pattern) %>%
  summarise(y_bot = min(ymin), y_top = max(ymax), .groups = "drop") %>%
  mutate(y_ctr = (y_top + y_bot) / 2, bar_h = y_top - y_bot)

# Group outline rectangles
group_outlines <- cat_heatmap_pos %>%
  mutate(xmin = HM_X_Y - HM_HW, xmax = ANNO_X + ANNO_HW,
         outline_col = PATTERN_COLORS[as.character(pattern)])

# Shared y-limits
S_Y_MIN  <- min(gene_bar_df$ymin) - gene_h * 2
S_Y_MAX  <- max(gene_bar_df$ymax) + gene_h * 2
S_Y_SPAN <- S_Y_MAX - S_Y_MIN

# --- Sankey data ---
cat_pw_counts <- links_1to1 %>%
  filter(pathway != "Other") %>%
  count(pattern, pathway, name = "n_proteins") %>%
  filter(n_proteins > 0)

cat_totals <- links_1to1 %>%
  filter(pathway != "Other") %>%
  count(pattern, name = "total") %>%
  arrange(match(as.character(pattern), PATTERN_ORDER))

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
  left_join(cat_heatmap_pos, by = "pattern") %>%
  mutate(fill = PATTERN_COLORS[as.character(pattern)])

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

cat_cum <- setNames(rep(0, n_cats), as.character(cat_totals$pattern))
pw_cum  <- setNames(rep(0, n_pws), pw_totals$pathway)
ribbon_list <- list()
ribbon_idx  <- 0

for (ct in as.character(cat_totals$pattern)) {
  ct_contribs <- cat_pw_counts %>%
    filter(as.character(pattern) == ct) %>%
    arrange(match(pathway, active_pw_order))
  for (r in seq_len(nrow(ct_contribs))) {
    pw <- ct_contribs$pathway[r]; n <- ct_contribs$n_proteins[r]
    if (n == 0) next
    ribbon_idx <- ribbon_idx + 1
    ct_row <- cat_bars %>% filter(as.character(pattern) == ct)
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
    ) %>% mutate(pattern = ct, pathway = pw, fill_col = PATTERN_COLORS[ct])
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
    arrange(match(as.character(pattern), PATTERN_ORDER))
  bar_y_ctr <- pw_row$y_ctr
  bar_total_w <- (pw_row$neg_log10_padj / max_enrich) * S_MAX_LEN
  x_cursor <- S_X_BAR
  for (r in seq_len(nrow(pw_contribs))) {
    ct <- as.character(pw_contribs$pattern[r])
    n  <- pw_contribs$n_proteins[r]
    if (n == 0) next
    stacked_idx <- stacked_idx + 1
    seg_w <- (n / pw_row$total) * bar_total_w
    stacked_rects[[stacked_idx]] <- tibble(
      xmin = x_cursor, xmax = x_cursor + seg_w,
      ymin = bar_y_ctr - S_SBAR_H / 2, ymax = bar_y_ctr + S_SBAR_H / 2,
      pattern = ct, pathway = pw, fill = PATTERN_COLORS[ct])
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
key_cats <- names(PATTERN_COLORS)
key_x0   <- S_X_BAR + S_MAX_LEN * 0.45
key_y0   <- S_Y_MAX
key_df <- tibble(cat = key_cats, color = PATTERN_COLORS[key_cats],
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
# 5. BUILD UNIFIED PANEL E PLOT
# ==============================================================================

pE <- ggplot() +
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
  annotate("text", x = HM_X_Y, y = header_y, label = "Aging",
           size = 5, fontface = "bold") +
  annotate("text", x = HM_X_O, y = header_y, label = "Training\n(Old)",
           size = 5, fontface = "bold", lineheight = 0.85) +
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
       title = "Reversal Classification: Heatmap & Pathway Enrichment",
       subtitle = sprintf("%d proteins | 4-pattern reversal classification | Hallmark + GO:BP ORA (BH adj.)",
                          nrow(class_df))) +
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

dir.create(file.path(DAT_DIR, "panel_E"), recursive = TRUE, showWarnings = FALSE)

class_df %>%
  filter(gene %in% gene_order) %>%
  transmute(gene, pattern = as.character(pattern),
            logFC_Aging = round(logFC_Aging, 4),
            logFC_Training_Old = round(logFC_Training_Old, 4),
            pi_score_Aging = round(pi_score_Aging, 6),
            pi_score_Reversal = round(pi_score_Reversal, 6),
            imputed) %>%
  arrange(pattern, desc(abs(logFC_Aging) + abs(logFC_Training_Old))) %>%
  write_csv(file.path(DAT_DIR, "panel_E", "classification.csv"))

links_1to1 %>%
  transmute(gene, pathway, database, pattern = as.character(pattern),
            logFC_Aging = round(logFC_Aging, 4),
            logFC_Training_Old = round(logFC_Training_Old, 4),
            pathway_pvalue = signif(pvalue, 4)) %>%
  arrange(pattern, pathway, gene) %>%
  write_csv(file.path(DAT_DIR, "panel_E", "sankey_links.csv"))

# Enrichment bars export
enrich_export <- ora_top %>%
  filter(pathway_label != "Other") %>%
  left_join(links_1to1 %>% count(pathway, name = "gene_count"),
            by = c("pathway_label" = "pathway")) %>%
  transmute(pathway = pathway_label, gene_count,
            gene_ratio = round(gene_count / n_class, 4),
            pvalue = signif(pvalue, 4),
            p.adjust = signif(p.adjust, 4),
            genes = geneID) %>%
  arrange(pvalue)
write_csv(enrich_export, file.path(DAT_DIR, "panel_E", "enrichment_bars.csv"))

# ==============================================================================
# 7. SAVE PANEL E
# ==============================================================================

f_height <- max(220, nrow(links_1to1) * 3.8 + 80)

ggsave(file.path(RPT_DIR, "panel_E_classification.pdf"), pE,
       width = 500, height = f_height, units = "mm", device = pdf, limitsize = FALSE)
ggsave(file.path(RPT_DIR, "panel_E_classification.png"), pE,
       width = 500, height = f_height, units = "mm", dpi = 300, limitsize = FALSE)

message("Panel E complete: ", file.path(RPT_DIR, "panel_E_classification.pdf"))
