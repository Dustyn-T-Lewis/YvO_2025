################################################################################
#   Figure 2 — Panel F: Interaction DEPs — Multi-Contrast Response &
#                        Pathway Enrichment (dumbbell | sankey | bar)
#
#   Generates:
#     b_reports/panel_F_interaction.pdf, panel_F_interaction.png
#     c_data/panel_F/*.csv
################################################################################

if (!exists("dep_df")) source("04_Figures/F2/a_script/YvO_figure2_shared.R")

message("Panel F: interaction classification (dumbbell + sankey + enrichment)...")

# ==============================================================================
# 1. PROTEIN SELECTION & CLASSIFICATION
# ==============================================================================

CATEGORY_COLORS <- c("Down Young / Up Old" = "#E05A4E",
                     "Attenuated"          = "#5DA5DA",
                     "Up Young / Down Old" = "#9B7FBF")
CATEGORY_ORDER  <- c("Down Young / Up Old", "Attenuated", "Up Young / Down Old")

# 4-pattern response classification (biological gradient for merged panel)
PATTERN_COLORS <- c(
  "Concordant"  = "#78909C",
  "Attenuated"  = "#5DA5DA",
  "Blunted"     = "#F4A460",
  "Reversed"    = "#E05A4E"
)
PATTERN_ORDER <- c("Concordant", "Attenuated", "Blunted", "Reversed")

# Classify all interaction DEPs into categories
int_class_all <- dep_df %>%
  filter(sig_pi_Interaction == TRUE) %>%
  transmute(
    gene,
    logFC_Y = logFC_Training_Young,
    logFC_O = logFC_Training_Old,
    logFC_Aging = logFC_Aging,
    divergence = abs(logFC_Training_Young - logFC_Training_Old),
    category = case_when(
      abs(logFC_O) < 0.15 ~ "Attenuated",
      logFC_Y > 0 & logFC_O < 0 ~ "Up Young / Down Old",
      logFC_Y < 0 & logFC_O > 0 ~ "Down Young / Up Old",
      sign(logFC_Y) == sign(logFC_O) & abs(logFC_O) < abs(logFC_Y) * 0.5 ~ "Attenuated",
      sign(logFC_Y) == sign(logFC_O) & abs(logFC_Y) < abs(logFC_O) * 0.5 ~ "Attenuated",
      logFC_Y < 0 & logFC_O <= 0 ~ "Down Young / Up Old",
      logFC_Y > 0 & logFC_O >= 0 ~ "Up Young / Down Old",
      TRUE ~ "Attenuated"
    )
  ) %>% filter(!is.na(logFC_Y), !is.na(logFC_O)) %>%
  mutate(
    response_pattern = case_when(
      logFC_O < -0.20                        ~ "Concordant",
      logFC_O >= -0.20 & logFC_O <= 0.15     ~ "Attenuated",
      logFC_O >  0.15  & logFC_O <= 0.30     ~ "Blunted",
      logFC_O >  0.30                         ~ "Reversed",
      TRUE                                    ~ "Attenuated"
    ) %>% factor(levels = PATTERN_ORDER)
  )

# Named lookups
gene_category_lookup <- setNames(int_class_all$category, int_class_all$gene)
gene_pattern_lookup  <- setNames(as.character(int_class_all$response_pattern),
                                  int_class_all$gene)

# ==============================================================================
# 2. ORA + PATHWAY MAPPING (pooled across all classified proteins)
# ==============================================================================

message("  Running pooled ORA (Hallmark + GO:BP)...")

int_class <- int_class_all
cat_counts <- int_class %>% count(category) %>% deframe()
int_class <- int_class %>%
  mutate(facet_label = sprintf("%s (n = %d)", category, cat_counts[category]))

all_int_genes   <- int_class$gene
all_genes       <- dep_df$gene

# ---- ORA: Hallmark + GO:BP ----
hallmark_t2g <- msigdbr(species = "Homo sapiens", collection = "H") %>%
  select(term = gs_name, gene = gene_symbol) %>% distinct()
gobp_t2g <- msigdbr(species = "Homo sapiens", collection = "C5",
                     subcollection = "GO:BP") %>%
  select(term = gs_name, gene = gene_symbol) %>% distinct()

run_ora <- function(t2g, db_name) {
  res <- tryCatch(
    enricher(gene = all_int_genes, universe = all_genes, TERM2GENE = t2g,
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
mapped_now <- ora_top %>% select(geneID) %>% separate_rows(geneID, sep = "/") %>%
  filter(geneID %in% all_int_genes) %>% pull(geneID) %>% unique()
orphans    <- setdiff(all_int_genes, mapped_now)
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
  select(pathway = pathway_label, ID, geneID, pvalue, p.adjust, Count, database) %>%
  separate_rows(geneID, sep = "/") %>% rename(gene = geneID) %>%
  filter(gene %in% all_int_genes)

# Rescue orphans via membership in selected pathways
leftover <- setdiff(all_int_genes, unique(sankey_links_all$gene))
if (length(leftover) > 0) {
  rescue <- all_t2g %>% filter(gene %in% leftover, term %in% ora_top$ID) %>%
    left_join(ora_top %>% select(ID, pathway = pathway_label, pvalue, p.adjust,
                                  Count, database), by = c("term" = "ID")) %>%
    group_by(gene) %>% slice_min(pvalue, n = 1, with_ties = FALSE) %>% ungroup() %>%
    mutate(ID = term) %>% select(gene, pathway, ID, pvalue, p.adjust, Count, database)
  sankey_links_all <- bind_rows(sankey_links_all, rescue)
  leftover <- setdiff(leftover, rescue$gene)
}

# Final catch-all: orphans that matched no selected pathway -> "Other"
leftover <- setdiff(all_int_genes, unique(sankey_links_all$gene))
if (length(leftover) > 0) {
  other_links <- tibble(
    gene     = leftover,
    pathway  = "Other",
    ID       = "OTHER",
    pvalue   = 1,
    p.adjust = 1,
    Count    = length(leftover),
    database = "Other"
  )
  sankey_links_all <- bind_rows(sankey_links_all, other_links)

  # Add "Other" to ora_top so downstream code picks it up
  ora_top <- bind_rows(ora_top, tibble(
    ID = "OTHER", Description = "Other", GeneRatio = NA, BgRatio = NA,
    pvalue = 1, p.adjust = 1, qvalue = 1, geneID = paste(leftover, collapse = "/"),
    Count = length(leftover), database = "Other", pathway_label = "Other"
  ))
  message(sprintf("  %d orphan genes assigned to 'Other' pathway", length(leftover)))
}

sankey_links_all <- sankey_links_all %>%
  left_join(int_class %>% select(gene, category, logFC_Y, logFC_O, divergence),
            by = "gene")
mapped_class <- int_class %>% filter(gene %in% unique(sankey_links_all$gene))

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
mapped_class <- int_class %>% filter(gene %in% unique(links_1to1$gene))
leftover <- setdiff(all_int_genes, unique(links_1to1$gene))

message(sprintf("  %d pathways, %d mapped (1:1), %d excluded",
                nrow(ora_top), nrow(links_1to1), length(leftover)))

# ---- Pathway and gene ordering ----
pw_1to1_counts <- links_1to1 %>% count(pathway, name = "n_1to1")
pw_pvals <- ora_top %>% select(pathway_label, pvalue)
pw_order_df <- pw_1to1_counts %>%
  left_join(pw_pvals, by = c("pathway" = "pathway_label")) %>%
  arrange(desc(pvalue))  # least significant first -> level 1
pw_order <- pw_order_df$pathway
# Ensure "Other" is always last (remove from sorted position, re-append)
if ("Other" %in% unique(links_1to1$pathway)) {
  pw_order <- setdiff(pw_order, "Other")
  pw_order <- c(pw_order, "Other")
}

# Gene order: pattern-based (biological gradient), then pathway within group
gene_order <- links_1to1 %>%
  mutate(
    pattern     = gene_pattern_lookup[gene],
    pattern_idx = match(pattern, PATTERN_ORDER),
    pw_idx      = match(pathway, pw_order)
  ) %>%
  arrange(pattern_idx, pw_idx, logFC_O) %>%
  pull(gene)

message("  ORA + pathway mapping complete")

# ==============================================================================
# 2b. DATA PREPARATION (Interaction DEP multi-contrast dumbbell data)
# ==============================================================================

message("  Data prep: interaction DEP multi-contrast data...")

# Select interaction DEPs and gather logFC across all 3 contrasts
int_deps <- dep_df %>%
  filter(sig_pi_Interaction == TRUE) %>%
  transmute(
    gene,
    logFC_Aging          = logFC_Aging,
    logFC_Training_Young = logFC_Training_Young,
    logFC_Training_Old   = logFC_Training_Old,
    pi_score_Interaction = pi_score_Interaction,
    # Compute interaction magnitude for ordering
    interaction_magnitude = abs(logFC_Training_Young - logFC_Training_Old)
  ) %>%
  filter(!is.na(logFC_Aging), !is.na(logFC_Training_Young), !is.na(logFC_Training_Old)) %>%
  arrange(desc(interaction_magnitude))

n_int <- nrow(int_deps)
message(sprintf("  %d interaction DEPs", n_int))

# Truncate to top 30 if too many for readability
max_show <- 30
if (n_int > max_show) {
  int_deps <- int_deps %>% slice_head(n = max_show)
  message(sprintf("  Showing top %d by interaction magnitude", max_show))
}

# Shared ordering: use pathway-grouped gene_order from Sankey;
# append any Panel F proteins not in Sankey at the end (by magnitude)
e_genes <- as.character(int_deps$gene)
shared_order <- c(gene_order, setdiff(e_genes, gene_order))
int_deps <- int_deps %>%
  mutate(gene = factor(gene, levels = rev(shared_order)))

# Assign each gene its interaction category color
int_deps <- int_deps %>%
  mutate(cat = gene_category_lookup[as.character(gene)])

# Pivot to long format for plotting
int_long <- int_deps %>%
  pivot_longer(
    cols = c(logFC_Aging, logFC_Training_Young, logFC_Training_Old),
    names_to = "contrast",
    values_to = "logFC",
    names_prefix = "logFC_"
  ) %>%
  mutate(
    contrast_label = case_when(
      contrast == "Aging"          ~ "Aging",
      contrast == "Training_Young" ~ "Training (Young)",
      contrast == "Training_Old"   ~ "Training (Old)"
    ),
    contrast_label = factor(contrast_label,
                            levels = c("Aging", "Training (Young)", "Training (Old)"))
  )

message("  Data prep complete")

# ==============================================================================
# 3. COORDINATE CONSTRUCTION (Sankey geometry)
# ==============================================================================

message("  Building Sankey coordinates...")

# Save classification for reference
write_csv(int_class, file.path(DAT_DIR, "panel_F", "interaction_classification.csv"))

# ---- Pathway colors: muted pastels (same palette as F3) ----
pw_palette <- c(
  "#F48FB1", "#FDAE91", "#E8E8A0", "#A8D8A8", "#8DD3C7",
  "#A2CEE5", "#B6C8E8", "#DEB4D4", "#C9A9A6", "#AED581",
  "#CE93D8", "#80DEEA", "#FFD54F", "#90A4AE", "#E57373",
  "#FFCC80", "#9FA8DA", "#B39DDB", "#80CBC4", "#81C784"
)
PW_COLORS <- setNames(pw_palette[seq_along(pw_order)], pw_order)
PW_COLORS["Other"] <- "#D0D0D0"  # neutral grey for catch-all

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
  fill_col   = PW_COLORS[pathway]
)

# ---- Stacking within pathway bars ----
slot_df <- links_1to1 %>%
  mutate(pw_idx = match(pathway, pw_order)) %>%
  arrange(pw_idx, match(gene, gene_order)) %>%
  group_by(pathway) %>%
  mutate(
    k        = n(),
    slot_idx = row_number()
  ) %>%
  ungroup() %>%
  left_join(pw_bar_df %>% select(pathway, pw_ymin = ymin, pw_ymax = ymax),
            by = "pathway") %>%
  mutate(
    slot_h   = (pw_ymax - pw_ymin) / k,
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
  left_join(gene_bar_df %>% select(gene, gene_ymin = ymin, gene_ymax = ymax),
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
  left_join(gene_bar_df %>% select(gene, gene_yctr = y_ctr), by = "gene") %>%
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
  mutate(gene_ratio  = Count_1to1 / length(all_int_genes),
         neg_log10_p = -log10(pvalue))

# ==============================================================================
# 4. EXPORT CSVs
# ==============================================================================

message("  Exporting Panel F data CSVs...")

dot_df %>%
  transmute(pathway = pathway_label, database, gene_count = Count_1to1,
            gene_ratio = round(gene_ratio, 4), pvalue = signif(pvalue, 4),
            p_adjust = signif(p.adjust, 4),
            genes = ora_top$geneID[match(pathway_label, ora_top$pathway_label)]) %>%
  arrange(pvalue) %>% write_csv(file.path(DAT_DIR, "panel_F", "sankey_dot.csv"))

links_1to1 %>%
  transmute(gene = as.character(gene), pathway = as.character(pathway), database,
            interaction_pattern = category, logFC_Training_Young = round(logFC_Y, 4),
            logFC_Training_Old = round(logFC_O, 4),
            pathway_pvalue = signif(pvalue, 4)) %>%
  arrange(interaction_pattern, pathway, gene) %>%
  write_csv(file.path(DAT_DIR, "panel_F", "sankey_links.csv"))

mapped_class %>%
  transmute(gene, interaction_pattern = category,
            logFC_Training_Young = round(logFC_Y, 4),
            logFC_Training_Old = round(logFC_O, 4),
            logFC_Aging = round(logFC_Aging, 4),
            divergence = round(divergence, 4)) %>%
  arrange(interaction_pattern, desc(divergence)) %>%
  write_csv(file.path(DAT_DIR, "panel_F", "interaction_patterns.csv"))

# ---- SRplot-format export (bioinformatics.com.cn Sankey-Dot tool input) ----
srplot_genes <- links_1to1 %>%
  group_by(pathway) %>%
  summarise(geneID = paste(sort(as.character(gene)), collapse = "/"),
            Count  = n(), .groups = "drop")

srplot_df <- ora_top %>%
  filter(pathway_label %in% active_pws) %>%
  select(Description = pathway_label, pvalue) %>%
  left_join(srplot_genes, by = c("Description" = "pathway")) %>%
  mutate(GeneRatio = Count / length(all_int_genes)) %>%
  select(Description, GeneRatio, pvalue, geneID, Count) %>%
  filter(pvalue > 0) %>%   # SRplot requires no zero p-values
  arrange(pvalue)

write_csv(srplot_df, file.path(DAT_DIR, "panel_F", "srplot_input.csv"))
message(sprintf("  SRplot input: %d pathways exported", nrow(srplot_df)))

int_deps %>%
  transmute(
    gene         = as.character(gene),
    logFC_Aging          = round(logFC_Aging, 4),
    logFC_Training_Young = round(logFC_Training_Young, 4),
    logFC_Training_Old   = round(logFC_Training_Old, 4),
    pi_score_Interaction = round(pi_score_Interaction, 6),
    interaction_magnitude = round(interaction_magnitude, 4)
  ) %>%
  arrange(desc(interaction_magnitude)) %>%
  write_csv(file.path(DAT_DIR, "panel_F", "interaction_dot.csv"))

int_long %>%
  transmute(
    gene = as.character(gene),
    contrast = contrast_label,
    logFC = round(logFC, 4)
  ) %>%
  write_csv(file.path(DAT_DIR, "panel_F", "interaction_dot_long.csv"))

message("  CSVs exported")

# ==============================================================================
# 5. DUMBBELL SUB-PANEL (left)
# ==============================================================================

message("  Building dumbbell sub-panel...")

dumbbell_df <- int_class %>%
  filter(gene %in% gene_order) %>%
  left_join(gene_bar_df %>% select(gene, y_ctr), by = "gene")

db_xrange <- range(c(dumbbell_df$logFC_Y, dumbbell_df$logFC_O,
                      dumbbell_df$logFC_Aging), na.rm = TRUE)
db_xpad <- diff(db_xrange) * 0.08

# Pattern group separators
dumbbell_df <- dumbbell_df %>%
  mutate(pattern = gene_pattern_lookup[gene])

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
  geom_hline(yintercept = sep_ys, linetype = "dotted", color = "grey50", linewidth = 0.4) +
  # Connecting line: Training Young to Training Old
  geom_segment(aes(x = logFC_Y, xend = logFC_O, y = y_ctr, yend = y_ctr),
               color = "grey65", linewidth = 0.4) +
  # Aging: open circle
  geom_point(aes(x = logFC_Aging, y = y_ctr), shape = 1, size = 1.8, stroke = 0.5,
             color = CONTRAST_COLORS["Aging"]) +
  # Training Young: filled
  geom_point(aes(x = logFC_Y, y = y_ctr), shape = 16, size = 1.8,
             color = CONTRAST_COLORS["Training_Young"]) +
  # Training Old: filled
  geom_point(aes(x = logFC_O, y = y_ctr), shape = 16, size = 1.8,
             color = CONTRAST_COLORS["Training_Old"]) +
  # Pattern group labels on the right edge
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
  # Gene bars (category-colored, ALL genes including orphans)
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
  geom_text(aes(x = gene_ratio + 0.005, y = dot_y,
                label = Count_1to1),
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
pat_key_f <- tibble(
  x     = 3.0,
  y     = c(3, 2, 1, 0) * ks,
  label = PATTERN_ORDER,
  color = unname(PATTERN_COLORS[PATTERN_ORDER])
)
title_y <- 4 * ks  # headers one step above top items

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
    title    = "Interaction DEPs: Multi-Contrast Response & Pathway Enrichment",
    subtitle = sprintf("%d proteins with significant Age x Training interaction (pi < 0.05)",
                       length(all_int_genes)),
    theme = theme(
      plot.title    = element_text(face = "bold", size = 10),
      plot.subtitle = element_text(size = 7, color = "grey30", face = "italic")
    )
  )

f_height <- max(200, nrow(links_1to1) * 7 + 60)

ggsave(file.path(RPT_DIR, "panel_F_interaction.pdf"), pF,
       width = 350, height = f_height, units = "mm", device = pdf)
ggsave(file.path(RPT_DIR, "panel_F_interaction.png"), pF,
       width = 350, height = f_height, units = "mm", dpi = 300)

# Save
message("  Panel F saved: ", file.path(RPT_DIR, "panel_F_interaction.pdf"))
message("  Panel F complete")
