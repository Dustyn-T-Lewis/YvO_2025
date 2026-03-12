################################################################################
#   Figure 3 — Panel F: Reversal Classification — Multi-Contrast Response &
#                        Pathway Enrichment (dumbbell | sankey | bar)
#
#   Uses the Reversal contrast (pi_score_Reversal < 0.05, ~180 DEPs) as the
#   input set, takes the top 30 by |logFC_Reversal| for readability, and
#   classifies into reversal response patterns based on Aging vs Training_Old
#   directionality.
#
#   Layout mirrors F2 Panel F (panel_F.R).
#
#   Generates:
#     b_reports/panel_F_classification.pdf, panel_F_classification.png
#     c_data/panel_F/classification.csv, sankey_links.csv, enrichment_bars.csv
################################################################################
#
# STAT AUDIT (2026-02-27)
# ---------------------------------------------------------------------------
# 1. Response pattern classification thresholds:
#    - "Reversed": sign(logFC_Aging) != sign(logFC_Training_Old)
#      Clear, biologically grounded — opposite-direction effects.        PASS
#    - "Attenuated": sig_A (pi_Aging < 0.05) AND
#      |logFC_Training_Old| < 0.5 * |logFC_Aging|
#      Threshold 0.5 is reasonable (training effect < half of aging).    PASS
#    - "Amplified": sig_A AND |logFC_Training_Old| > |logFC_Aging|
#      Training exceeds aging magnitude.                                 PASS
#    - "Concordant": fallback for remaining cases.                       PASS
#    - NOTE: Gap between Attenuated (< 0.5x) and Amplified (> 1.0x)
#      means proteins with 0.5x <= |ratio| <= 1.0x fall to "Concordant".
#      This is acceptable as a moderate-response catch-all.              PASS
#
# 2. ORA (Over-Representation Analysis):
#    - Universe: all_genes = dep_df$gene (full quantified proteome).     PASS
#      This is the correct universe — all proteins that could have been
#      detected, not just DEPs. Confirmed: dep_df contains all proteins
#      including non-significant ones.
#    - enricher() from clusterProfiler used with BH correction.          PASS
#    - Pathway selection: top 5 Hallmark + 5 GO:BP by p-value, then
#      greedy rescue for orphan genes (up to 15 total pathways).
#      Greedy approach ensures coverage; cap prevents over-fitting.      PASS
#    - pvalueCutoff = 1: all pathways tested (filtered post-hoc).        PASS
#
# 3. 1:1 gene-pathway assignment:
#    - Each gene assigned to best pathway (lowest p-value).              PASS
#    - Redistribution loop moves solo-pathway genes to multi-gene
#      pathways to reduce visual clutter. Deterministic (5 iterations).  PASS
#
# 4. No formal statistical test on classification proportions.
#    This is a descriptive visualization of response patterns.           PASS
# ---------------------------------------------------------------------------

if (!exists("dep_df")) source("04_Figures/F3/a_script/YvO_F3_setup.R")

message("Panel F: reversal classification (dumbbell + sankey + enrichment)...")

# ==============================================================================
# 1. PROTEIN SELECTION & CLASSIFICATION
# ==============================================================================

# Use Reversal DEPs (pi_score_Reversal < 0.05) — the core F3 contrast.
# Take top 30 by |logFC_Reversal| for readability (F2 Panel F caps at 30).
max_show <- 30

class_df <- dep_df %>%
  filter(pi_score_Reversal < 0.05) %>%
  filter(!is.na(logFC_Aging), !is.na(logFC_Training_Old),
         !is.na(logFC_Training_Young)) %>%
  arrange(desc(abs(logFC_Reversal))) %>%
  slice_head(n = max_show) %>%
  transmute(
    gene,
    logFC_Aging, logFC_Training_Young, logFC_Training_Old,
    logFC_Reversal,
    pi_Aging = pi_score_Aging,
    pi_TO    = pi_score_Training_Old,
    sig_A    = pi_score_Aging < 0.05
  ) %>%
  mutate(
    category = case_when(
      sign(logFC_Aging) != sign(logFC_Training_Old) ~ "Reversed",
      sig_A & abs(logFC_Training_Old) < abs(logFC_Aging) * 0.5 ~ "Attenuated",
      sig_A & abs(logFC_Training_Old) > abs(logFC_Aging)       ~ "Amplified",
      TRUE ~ "Concordant"
    )
  )

# 4-pattern response classification (colors matching F2 palette)
PATTERN_COLORS <- c(
  "Reversed"    = "#E05A4E",
  "Attenuated"  = "#5DA5DA",
  "Amplified"   = "#F4A460",
  "Concordant"  = "#78909C"
)
PATTERN_ORDER <- c("Reversed", "Attenuated", "Amplified", "Concordant")

class_df <- class_df %>%
  mutate(category = factor(category, levels = PATTERN_ORDER))

n_class <- nrow(class_df)
cat_counts <- class_df %>% count(category) %>% deframe()
message(sprintf("  Classification: %d proteins across %d categories",
                n_class, length(cat_counts)))
for (cat in names(cat_counts))
  message(sprintf("    %s: %d", cat, cat_counts[cat]))

# Named lookups
gene_pattern_lookup <- setNames(as.character(class_df$category), class_df$gene)

# ==============================================================================
# 2. ORA + PATHWAY MAPPING (pooled across all classified proteins)
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
  left_join(class_df %>% dplyr::select(gene, category, logFC_Aging,
                                        logFC_Training_Young, logFC_Training_Old),
            by = "gene")

# ---- Force 1:1 mapping: each gene -> best pathway ----
links_1to1 <- sankey_links_all %>%
  group_by(gene) %>%
  slice_min(pvalue, n = 1, with_ties = FALSE) %>%
  ungroup()

# Redistribute from over-represented pathways
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

# ---- Pathway and gene ordering ----
pw_1to1_counts <- links_1to1 %>% count(pathway, name = "n_1to1")
pw_pvals <- ora_top %>% dplyr::select(pathway_label, pvalue)
pw_order_df <- pw_1to1_counts %>%
  left_join(pw_pvals, by = c("pathway" = "pathway_label")) %>%
  arrange(desc(pvalue))
pw_order <- pw_order_df$pathway
if ("Other" %in% active_pws) {
  pw_order <- setdiff(pw_order, "Other")
  pw_order <- c(pw_order, "Other")
}

# Gene order: pattern-based, then pathway within group
gene_order <- links_1to1 %>%
  mutate(
    pattern     = gene_pattern_lookup[gene],
    pattern_idx = match(pattern, PATTERN_ORDER),
    pw_idx      = match(pathway, pw_order)
  ) %>%
  arrange(pattern_idx, pw_idx, logFC_Training_Old) %>%
  pull(gene)

message("  ORA + pathway mapping complete")

# ==============================================================================
# 3. COORDINATE CONSTRUCTION (Sankey geometry)
# ==============================================================================

message("  Building Sankey coordinates...")

# ---- Pathway colors: muted pastels (same palette as F2) ----
pw_palette <- c(
  "#F48FB1", "#FDAE91", "#E8E8A0", "#A8D8A8", "#8DD3C7",
  "#A2CEE5", "#B6C8E8", "#DEB4D4", "#C9A9A6", "#AED581",
  "#CE93D8", "#80DEEA", "#FFD54F", "#90A4AE", "#E57373",
  "#FFCC80", "#9FA8DA", "#B39DDB", "#80CBC4", "#81C784"
)
PW_COLORS <- setNames(pw_palette[seq_along(pw_order)], pw_order)
PW_COLORS["Other"] <- "#D0D0D0"

# ---- Coordinate system ----
n_genes <- nrow(links_1to1)
n_pw    <- length(pw_order)
Y_SPAN  <- n_genes

X_GENE <- 1.0;  X_PW <- 3.0;  BAR_W <- 0.06

gene_h   <- Y_SPAN / (n_genes * 1.15)
gene_gap <- (Y_SPAN - n_genes * gene_h) / max(n_genes - 1, 1)

pw_h   <- Y_SPAN / (n_pw * 1.4)
pw_gap <- (Y_SPAN - n_pw * pw_h) / max(n_pw - 1, 1)

# ---- Gene bar data frame ----
gene_bar_df <- tibble(
  gene = gene_order,
  idx  = seq_along(gene_order)
) %>% mutate(
  ymax = Y_SPAN - (idx - 1) * (gene_h + gene_gap),
  ymin = ymax - gene_h,
  y_ctr = (ymin + ymax) / 2,
  fill_col = PATTERN_COLORS[gene_pattern_lookup[gene]]
)

# ---- Pathway bar data frame ----
pw_bar_df <- tibble(
  pathway = pw_order,
  idx     = seq_along(pw_order)
) %>% mutate(
  ymax = Y_SPAN - (idx - 1) * (pw_h + pw_gap),
  ymin = ymax - pw_h,
  y_ctr = (ymin + ymax) / 2,
  fill_col = PW_COLORS[pathway]
)

# ---- Stacking within pathway bars ----
slot_df <- links_1to1 %>%
  mutate(pw_idx = match(pathway, pw_order)) %>%
  arrange(pw_idx, match(gene, gene_order)) %>%
  group_by(pathway) %>%
  mutate(k = n(), slot_idx = row_number()) %>%
  ungroup() %>%
  left_join(pw_bar_df %>% dplyr::select(pathway, pw_ymin = ymin, pw_ymax = ymax),
            by = "pathway") %>%
  mutate(
    slot_h    = (pw_ymax - pw_ymin) / k,
    slot_ymax = pw_ymax - (slot_idx - 1) * slot_h,
    slot_ymin = slot_ymax - slot_h
  )

# ---- Sigmoid ribbon polygons ----
make_sigmoid_ribbon <- function(x0, x1, y0_top, y0_bot, y1_top, y1_bot,
                                n_pts = 50, ribbon_id) {
  t <- seq(0, 1, length.out = n_pts)
  blend <- (1 - cos(pi * t)) / 2
  tibble(
    x = c(x0 + (x1 - x0) * t, rev(x0 + (x1 - x0) * t)),
    y = c(y0_top + (y1_top - y0_top) * blend,
          rev(y0_bot + (y1_bot - y0_bot) * blend)),
    ribbon_id = ribbon_id
  )
}

ribbon_input <- slot_df %>%
  left_join(gene_bar_df %>% dplyr::select(gene, gene_ymin = ymin, gene_ymax = ymax),
            by = "gene")

ribbon_polys <- pmap_dfr(ribbon_input, function(gene, pathway,
                                                  gene_ymin, gene_ymax,
                                                  slot_ymin, slot_ymax, ...) {
  make_sigmoid_ribbon(
    x0 = X_GENE + BAR_W / 2,  x1 = X_PW - BAR_W / 2,
    y0_top = gene_ymax, y0_bot = gene_ymin,
    y1_top = slot_ymax, y1_bot = slot_ymin,
    ribbon_id = paste(gene, pathway, sep = "->")
  ) %>% mutate(pathway = pathway, gene = gene)
})

# Ribbon fill: pathway color
ribbon_polys <- ribbon_polys %>%
  mutate(fill_col = PW_COLORS[pathway])

# Draw order: crossing ribbons behind shorter ones
ribbon_order <- ribbon_input %>%
  left_join(gene_bar_df %>% dplyr::select(gene, gene_yctr = y_ctr), by = "gene") %>%
  mutate(
    slot_yctr  = (slot_ymin + slot_ymax) / 2,
    crossing_d = abs(gene_yctr - slot_yctr),
    ribbon_id  = paste(gene, pathway, sep = "->")
  ) %>%
  arrange(desc(crossing_d)) %>%
  pull(ribbon_id)

ribbon_polys <- ribbon_polys %>%
  mutate(ribbon_id = factor(ribbon_id, levels = ribbon_order))

# ---- Enrichment bar data ----
srplot_summary <- links_1to1 %>%
  group_by(pathway) %>%
  summarise(Count_1to1 = n(), .groups = "drop")

dot_df <- ora_top %>%
  left_join(srplot_summary, by = c("pathway_label" = "pathway")) %>%
  left_join(pw_bar_df %>% transmute(pathway_label = as.character(pathway),
                                     dot_y = y_ctr),
            by = "pathway_label") %>%
  mutate(gene_ratio  = Count_1to1 / length(all_class_genes),
         neg_log10_p = -log10(pvalue))

# ==============================================================================
# 4. EXPORT CSVs
# ==============================================================================

message("  Exporting Panel F data CSVs...")

class_df %>%
  transmute(gene, logFC_Aging = round(logFC_Aging, 4),
            logFC_Training_Young = round(logFC_Training_Young, 4),
            logFC_Training_Old = round(logFC_Training_Old, 4),
            pi_score_Aging = round(pi_Aging, 6),
            pi_score_Training_Old = round(pi_TO, 6),
            category = as.character(category)) %>%
  arrange(category, desc(abs(logFC_Aging))) %>%
  write_csv(file.path(DAT_DIR, "panel_F", "classification.csv"))

links_1to1 %>%
  transmute(gene = as.character(gene), pathway = as.character(pathway), database,
            reversal_pattern = gene_pattern_lookup[gene],
            logFC_Aging = round(logFC_Aging, 4),
            logFC_Training_Young = round(logFC_Training_Young, 4),
            logFC_Training_Old = round(logFC_Training_Old, 4),
            pathway_pvalue = signif(pvalue, 4)) %>%
  arrange(reversal_pattern, pathway, gene) %>%
  write_csv(file.path(DAT_DIR, "panel_F", "sankey_links.csv"))

dot_df %>%
  transmute(pathway = pathway_label, database, gene_count = Count_1to1,
            gene_ratio = round(gene_ratio, 4), pvalue = signif(pvalue, 4),
            p_adjust = signif(p.adjust, 4),
            genes = ora_top$geneID[match(pathway_label, ora_top$pathway_label)]) %>%
  arrange(pvalue) %>%
  write_csv(file.path(DAT_DIR, "panel_F", "enrichment_bars.csv"))

message("  CSVs exported")

# ==============================================================================
# 5. DUMBBELL SUB-PANEL (left)
# ==============================================================================

message("  Building dumbbell sub-panel...")

dumbbell_df <- class_df %>%
  filter(gene %in% gene_order) %>%
  left_join(gene_bar_df %>% dplyr::select(gene, y_ctr), by = "gene") %>%
  mutate(pattern = as.character(category))

db_xrange <- range(c(dumbbell_df$logFC_Training_Young,
                      dumbbell_df$logFC_Training_Old,
                      dumbbell_df$logFC_Aging), na.rm = TRUE)
db_xpad <- diff(db_xrange) * 0.08

# Pattern group separators
pattern_genes <- split(gene_order, gene_pattern_lookup[gene_order])
sep_ys <- numeric(0)
for (i in seq_len(length(PATTERN_ORDER) - 1)) {
  grp_cur  <- pattern_genes[[PATTERN_ORDER[i]]]
  grp_next <- pattern_genes[[PATTERN_ORDER[i + 1]]]
  if (length(grp_cur) > 0 && length(grp_next) > 0) {
    y_last  <- gene_bar_df$ymin[gene_bar_df$gene == tail(grp_cur, 1)]
    y_first <- gene_bar_df$ymax[gene_bar_df$gene == head(grp_next, 1)]
    sep_ys  <- c(sep_ys, (y_last + y_first) / 2)
  }
}

# Pattern group midpoints for bracket labels
pattern_label_df <- tibble(pattern = PATTERN_ORDER) %>%
  rowwise() %>%
  mutate(
    genes_in = list(pattern_genes[[pattern]]),
    n_genes  = length(genes_in),
    mid_y    = if (n_genes > 0)
      mean(gene_bar_df$y_ctr[gene_bar_df$gene %in% genes_in]) else NA_real_
  ) %>%
  ungroup() %>%
  filter(n_genes > 0)

pF_dumbbell <- ggplot(dumbbell_df) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey60", linewidth = 0.3) +
  # Pattern group separators
  geom_hline(yintercept = sep_ys, linetype = "dotted", color = "grey50",
             linewidth = 0.4) +
  # Connecting segment: Aging to Training_Old (the reversal distance)
  geom_segment(aes(x = logFC_Aging, xend = logFC_Training_Old,
                    y = y_ctr, yend = y_ctr),
               color = "grey65", linewidth = 0.4) +
  # Aging: open circle (green)
  geom_point(aes(x = logFC_Aging, y = y_ctr), shape = 1, size = 1.8,
             stroke = 0.5, color = CONTRAST_COLORS["Aging"]) +
  # Training Old: filled (blue)
  geom_point(aes(x = logFC_Training_Old, y = y_ctr), shape = 16, size = 1.8,
             color = CONTRAST_COLORS["Training_Old"]) +
  # Training Young: filled (red)
  geom_point(aes(x = logFC_Training_Young, y = y_ctr), shape = 16, size = 1.8,
             color = CONTRAST_COLORS["Training_Young"]) +
  # Pattern group labels
  annotate("text",
           x = db_xrange[2] + db_xpad * 0.9,
           y = pattern_label_df$mid_y,
           label = pattern_label_df$pattern,
           hjust = 0, size = 2.2, fontface = "bold.italic",
           color = PATTERN_COLORS[pattern_label_df$pattern]) +
  scale_y_continuous(limits = c(0, Y_SPAN), expand = expansion(mult = 0.02),
                     breaks = gene_bar_df$y_ctr, labels = gene_bar_df$gene) +
  scale_x_continuous(position = "top") +
  coord_cartesian(xlim = c(db_xrange[1] - db_xpad,
                            db_xrange[2] + db_xpad * 3.5),
                  clip = "off") +
  labs(x = expression(log[2]~FC), y = NULL) +
  THEME_PUB +
  theme(
    axis.text.y        = element_text(size = 6.5, face = "italic", hjust = 1),
    axis.ticks.y       = element_blank(),
    panel.grid.major.y = element_blank(),
    panel.grid.minor   = element_blank(),
    panel.grid.major.x = element_line(color = "grey92", linewidth = 0.3),
    panel.border       = element_rect(color = "black", linewidth = 0.6),
    plot.margin        = margin(8, 0, 8, 2)
  )

# ==============================================================================
# 6. SANKEY SUB-PANEL (center)
# ==============================================================================

message("  Building Sankey sub-panel...")

ribbon_polys_f <- ribbon_polys %>% filter(pathway != "Other")
pw_bar_f       <- pw_bar_df %>% filter(pathway != "Other")

pF_sankey <- ggplot() +
  # Sigmoid ribbons (excluding "Other")
  geom_polygon(data = ribbon_polys_f %>% arrange(ribbon_id),
               aes(x = x, y = y, group = ribbon_id, fill = fill_col),
               alpha = 0.32, color = NA) +
  scale_fill_identity() +
  # Gene bars (pattern-colored, ALL genes including orphans)
  geom_rect(data = gene_bar_df,
            aes(xmin = X_GENE - BAR_W / 2, xmax = X_GENE + BAR_W / 2,
                ymin = ymin, ymax = ymax),
            fill = gene_bar_df$fill_col, color = NA) +
  # Pathway bars (excluding "Other")
  geom_rect(data = pw_bar_f,
            aes(xmin = X_PW - BAR_W / 2, xmax = X_PW + BAR_W / 2,
                ymin = ymin, ymax = ymax),
            fill = pw_bar_f$fill_col, color = NA) +
  # Pathway labels (left of bars, excluding "Other")
  geom_text(data = pw_bar_f,
            aes(x = X_PW - BAR_W / 2 - 0.05, y = y_ctr, label = pathway),
            hjust = 1, size = 2.8, fontface = "bold") +
  scale_y_continuous(limits = c(0, Y_SPAN), expand = expansion(mult = 0.02)) +
  coord_cartesian(xlim = c(X_GENE - 0.05, X_PW + 0.15), clip = "off") +
  theme_void() +
  theme(plot.margin = margin(8, 0, 8, 0))

# ==============================================================================
# 7. ENRICHMENT BAR SUB-PANEL (right)
# ==============================================================================

message("  Building enrichment bar sub-panel...")

dot_df_f <- dot_df %>% filter(pathway_label != "Other")

# Bar height matches pathway bar height for visual alignment
bar_half_h <- pw_h * 0.35

dot_df_f <- dot_df_f %>%
  mutate(
    bar_fill = PW_COLORS[pathway_label],
    ymin_bar = dot_y - bar_half_h,
    ymax_bar = dot_y + bar_half_h
  )

# Use the SAME y-axis as the Sankey (0 to Y_SPAN) so bars align 1:1
pF_dot <- ggplot(dot_df_f) +
  # Pathway-colored horizontal bars
  geom_rect(aes(xmin = 0, xmax = gene_ratio,
                ymin = ymin_bar, ymax = ymax_bar),
            fill = dot_df_f$bar_fill, color = "grey30", linewidth = 0.2) +
  # Count label at end of bar
  geom_text(aes(x = gene_ratio + 0.005, y = dot_y, label = Count_1to1),
            hjust = 0, size = 2.5, fontface = "bold") +
  scale_y_continuous(limits = c(0, Y_SPAN), expand = expansion(mult = 0.02)) +
  scale_x_continuous(expand = expansion(mult = c(0, 0.15))) +
  labs(x = "Gene Ratio", y = NULL) +
  THEME_PUB +
  theme(
    axis.text.y        = element_blank(),
    axis.ticks.y       = element_blank(),
    panel.grid.major.x = element_line(color = "grey92", linewidth = 0.3),
    panel.grid.major.y = element_blank(),
    panel.grid.minor   = element_blank(),
    panel.border       = element_rect(color = "black", linewidth = 0.6),
    legend.position    = "none",
    plot.margin        = margin(8, 5, 8, 0)
  )

# ==============================================================================
# 8. LEGEND + ASSEMBLY
# ==============================================================================

message("  Building legend and assembling Panel F...")

ks <- 0.10

contrast_key_f <- tibble(
  x     = 0,
  y     = c(3, 2, 1) * ks,
  label = c("Aging", "Training (Old)", "Training (Young)"),
  color = unname(CONTRAST_COLORS[c("Aging", "Training_Old", "Training_Young")]),
  shape = c(1, 16, 16)
)

# Only include patterns that actually appear
present_patterns <- PATTERN_ORDER[PATTERN_ORDER %in% unique(gene_pattern_lookup)]
pat_key_f <- tibble(
  x     = 3.0,
  y     = rev(seq_along(present_patterns)) * ks,
  label = present_patterns,
  color = unname(PATTERN_COLORS[present_patterns])
)
title_y <- (max(c(contrast_key_f$y, pat_key_f$y)) + 1.2 * ks)

pF_key <- ggplot() +
  # Contrast column
  annotate("text", x = 0, y = title_y, label = "Contrast",
           hjust = 0, size = 2.8, fontface = "bold", color = "grey30") +
  geom_point(data = contrast_key_f, aes(x = x, y = y),
             shape = contrast_key_f$shape, size = 3.0,
             color = contrast_key_f$color, stroke = 0.5) +
  geom_text(data = contrast_key_f, aes(x = x + 0.3, y = y, label = label),
            hjust = 0, size = 2.5, fontface = "bold",
            color = contrast_key_f$color) +
  # Response column
  annotate("text", x = 3.0, y = title_y, label = "Response",
           hjust = 0, size = 2.8, fontface = "bold", color = "grey30") +
  geom_point(data = pat_key_f, aes(x = x, y = y),
             shape = 15, size = 3.0, color = pat_key_f$color) +
  geom_text(data = pat_key_f, aes(x = x + 0.3, y = y, label = label),
            hjust = 0, size = 2.5, fontface = "bold",
            color = pat_key_f$color) +
  scale_x_continuous(limits = c(-0.3, 7)) +
  scale_y_continuous(limits = c(-0.3 * ks, title_y + 0.5 * ks)) +
  theme_void() +
  theme(plot.margin = margin(0, 5, 0, 0))

# ---- Patchwork area design (4 panels on 30x100 grid) ----
f_design <- c(
  patchwork::area(1,  1,  26, 28),   # dumbbell: full height, left 28%
  patchwork::area(1,  29, 26, 68),   # sankey:   full height, middle 40%
  patchwork::area(1,  69, 26, 100),  # bar:      full height, right 32%
  patchwork::area(27, 69, 30, 100)   # key:      below bar
)

pF <- pF_dumbbell + pF_sankey + pF_dot + pF_key +
  plot_layout(design = f_design) +
  plot_annotation(
    title    = "Reversal Classification: Multi-Contrast Response & Pathway Enrichment",
    subtitle = sprintf("Top %d Reversal DEPs (pi < 0.05) classified by Aging vs Training_Old response",
                       length(all_class_genes)),
    theme = theme(
      plot.title    = element_text(face = "bold", size = 10),
      plot.subtitle = element_text(size = 7, color = "grey30", face = "italic")
    )
  )

# Save
f_height <- max(200, nrow(links_1to1) * 7 + 60)

ggsave(file.path(RPT_DIR, "panel_F_classification.pdf"), pF,
       width = 350, height = f_height, units = "mm", device = pdf)
ggsave(file.path(RPT_DIR, "panel_F_classification.png"), pF,
       width = 350, height = f_height, units = "mm", dpi = 300)

message("  Panel F saved: ", file.path(RPT_DIR, "panel_F_classification.pdf"))
message("  Panel F complete")
