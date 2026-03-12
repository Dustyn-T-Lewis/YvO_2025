################################################################################
#   Figure 3 — Panel F: Aging & Reversal Classification — Heatmap +
#                        Category-Pathway Sankey + Stacked Composition Bars
#
#   All proteins significant in the Aging OR Reversal contrast (no top-N cap),
#   classified into 4 groups by direction. Two heatmap columns
#   (Aging, Training_Old), annotation strip, Sankey to GO Slim super-categories,
#   and stacked enrichment bars.
#
#   Mirrors F2 Panel F unified coordinate approach for cross-figure consistency.
#
#   Generates:
#     b_reports/panel_F_classification.pdf, panel_F_classification.png
#     c_data/panel_F/classification.csv, sankey_links.csv, enrichment_bars.csv
################################################################################
#
# STAT AUDIT (2026-03-03; gene-count bars + consolidated pathways 2026-03-03)
# ---------------------------------------------------------------------------
# 1. Classification (4 groups):
#    - "Aging Up": pi_score_Aging < 0.05 AND NOT Reversal, logFC_Aging > 0.
#    - "Aging Down": pi_score_Aging < 0.05 AND NOT Reversal, logFC_Aging < 0.
#    - "Rev. (Age Up)": pi_score_Reversal < 0.05, logFC_Aging > 0.
#      Aging pushed protein up; training reversed it down.
#    - "Rev. (Age Down)": pi_score_Reversal < 0.05, logFC_Aging < 0.
#      Aging pushed protein down; training reversed it up.
#    - Priority: Reversed > Aging (a protein sig in both is classified as
#      Reversed). Biologically sound — reversal is the stronger claim.     PASS
#
# 2. GO Slim Consolidated Pathway Mapping:
#    - Uses assign_go_slim_super() from shared/go_slim_categories.R.       PASS
#    - 12 consolidated pathways ("Signaling" removed — GO:0023052 too broad). PASS
#    - Consistent 1:1 assignment for both foreground AND background via the
#      same function (no apples-to-oranges comparison).                    PASS
#    - Fisher's exact test per consolidated pathway; BH correction.        PASS
#    - Bars show gene counts (not fold enrichment).                        PASS
#
# 3. Visualization is descriptive — no formal test on proportions.         PASS
# ---------------------------------------------------------------------------

if (!exists("dep_df")) source("04_Figures/F3/a_script/YvO_F3_setup.R")

message("Panel F: aging & reversal classification (heatmap + sankey + stacked bars)...")

# ==============================================================================
# 1. CLASSIFICATION: Aging Up / Aging Down / Rev. (Age Up) / Rev. (Age Down)
# ==============================================================================

F3_CLASS_ORDER <- c("Aging Up", "Aging Down", "Rev. (Age Up)", "Rev. (Age Down)")

F3_CLASS_COLORS <- c(
  "Aging Up"        = "#D6604D",  # warm red (aging-upregulated, not reversed)
  "Aging Down"      = "#4393C3",  # cool blue (aging-downregulated, not reversed)
  "Rev. (Age Up)"   = "#E65100",  # deep orange (was aging-up, training pushed down)
  "Rev. (Age Down)" = "#FFB300"   # amber (was aging-down, training pushed up)
)

class_df <- dep_df %>%
  filter(!is.na(logFC_Aging), !is.na(logFC_Training_Old)) %>%
  filter(pi_score_Aging < 0.05 | pi_score_Reversal < 0.05) %>%
  mutate(
    f3_class = case_when(
      pi_score_Reversal < 0.05 & logFC_Aging > 0  ~ "Rev. (Age Up)",
      pi_score_Reversal < 0.05 & logFC_Aging < 0  ~ "Rev. (Age Down)",
      pi_score_Reversal < 0.05 & logFC_Aging == 0 ~ "Rev. (Age Up)",
      pi_score_Aging < 0.05 & logFC_Aging > 0     ~ "Aging Up",
      pi_score_Aging < 0.05 & logFC_Aging < 0     ~ "Aging Down"
    )
  ) %>%
  filter(!is.na(f3_class)) %>%
  mutate(f3_class = factor(f3_class, levels = F3_CLASS_ORDER))

n_class <- nrow(class_df)
cat_counts <- class_df %>% count(f3_class) %>% deframe()
message(sprintf("  %d significant proteins across 4 classes", n_class))
for (ct in names(cat_counts))
  message(sprintf("    %s: %d", ct, cat_counts[ct]))

# ==============================================================================
# 2. GO SLIM CONSOLIDATED PATHWAY MAPPING (via shared utility)
# ==============================================================================

message("  Mapping to GO Slim consolidated pathways...")

sig_genes  <- class_df$gene
all_genes  <- dep_df$gene

# Foreground: 1:1 assignment for significant genes
best_super <- assign_go_slim_super(sig_genes, all_genes)

message(sprintf("  %d genes in %d consolidated pathways",
                nrow(best_super), n_distinct(best_super$super)))

# Build links_1to1 (gene -> pathway with class metadata)
links_1to1 <- best_super %>%
  select(gene, pathway = super) %>%
  mutate(pathway = as.character(pathway)) %>%
  left_join(class_df %>% select(gene, f3_class, logFC_Aging, logFC_Training_Old),
            by = "gene") %>%
  mutate(database = ifelse(pathway == "Other", "Other", "GO_Slim"))

# Fisher's exact test per consolidated pathway (for reference; bars now show gene counts)
active_supers <- unique(links_1to1$pathway[links_1to1$pathway != "Other"])

# Background: consistent 1:1 assignment via assign_go_slim_super() (same method as foreground)
all_super <- assign_go_slim_super(all_genes, all_genes)

fisher_results <- tibble(
  pathway_label = active_supers,
  pvalue = sapply(active_supers, function(s) {
    fg_in  <- sum(links_1to1$pathway == s)
    fg_out <- length(sig_genes) - fg_in
    bg_in  <- sum(all_super$super == s)
    bg_out <- nrow(all_super) - bg_in
    fisher.test(matrix(c(fg_in, bg_in, fg_out, bg_out), 2, 2),
                alternative = "greater")$p.value
  })
) %>%
  mutate(p.adjust = p.adjust(pvalue, method = "BH"),
         ID = pathway_label,
         database = "GO_Slim")

if ("Other" %in% links_1to1$pathway) {
  fisher_results <- bind_rows(fisher_results, tibble(
    pathway_label = "Other", pvalue = 1, p.adjust = 1,
    ID = "OTHER", database = "Other"
  ))
}

ora_top <- fisher_results

# Add per-gene pvalue (pathway's Fisher p) to links_1to1
pw_pvals_map <- setNames(ora_top$pvalue, ora_top$pathway_label)
links_1to1 <- links_1to1 %>% mutate(pvalue = pw_pvals_map[pathway])

message(sprintf("  %d pathways, %d mapped (1:1)", nrow(ora_top), nrow(links_1to1)))

# ==============================================================================
# 3. ORDERING
# ==============================================================================

# Pathway ordering by p-value (least significant first, Other always last)
pw_1to1_counts <- links_1to1 %>% count(pathway, name = "n_1to1")
pw_pvals <- ora_top %>% select(pathway_label, pvalue)
pw_order_df <- pw_1to1_counts %>%
  left_join(pw_pvals, by = c("pathway" = "pathway_label")) %>%
  arrange(desc(pvalue))
pw_order <- pw_order_df$pathway
if ("Other" %in% pw_order) {
  pw_order <- c(setdiff(pw_order, "Other"), "Other")
}

# Gene ordering: by class -> pathway -> logFC_Aging within group
gene_order <- links_1to1 %>%
  mutate(
    class_idx = match(as.character(f3_class), F3_CLASS_ORDER),
    pw_idx    = match(pathway, pw_order)
  ) %>%
  arrange(class_idx, pw_idx, logFC_Aging) %>%
  pull(gene)

# ==============================================================================
# 4-5. UNIFIED HEATMAP + ANNOTATION STRIP + SANKEY + ENRICHMENT BARS
#      Single ggplot for perfect y-alignment across all elements.
# ==============================================================================

message("  Building unified Panel F (single coordinate system)...")

n_genes <- length(gene_order)
Y_SPAN  <- n_genes

# --- Gene row geometry (top-down, tiles touching within groups) ---
n_cat_breaks <- length(unique(class_df$f3_class)) - 1
gap_size     <- 1.2        # narrow gap between class groups
total_gap    <- n_cat_breaks * gap_size
gene_h       <- (Y_SPAN - total_gap) / n_genes

gene_bar_df <- tibble(gene = gene_order, idx = seq_along(gene_order)) %>%
  mutate(
    ymax  = Y_SPAN - (idx - 1) * gene_h,
    ymin  = ymax - gene_h,
    y_ctr = (ymin + ymax) / 2
  ) %>%
  left_join(links_1to1 %>% select(gene, f3_class) %>% distinct(), by = "gene")

# Category gaps between groups only (tiles touch within groups)
cat_breaks <- gene_bar_df %>%
  mutate(cat_str = as.character(f3_class)) %>%
  mutate(is_break = cat_str != lag(cat_str, default = cat_str[1]))
cum_gap <- cumsum(cat_breaks$is_break) * gap_size
gene_bar_df <- gene_bar_df %>%
  mutate(
    ymax  = ymax - cum_gap,
    ymin  = ymin - cum_gap,
    y_ctr = (ymin + ymax) / 2
  )

# --- X-coordinate layout (all elements share one coordinate system) ---
HM_X_A   <- 1       # Aging column center
HM_X_TO  <- 1.92    # Training_Old column center
HM_HW    <- 0.45    # heatmap tile half-width
ANNO_X   <- HM_X_TO + HM_HW + 0.01 + 0.12   # annotation strip center
ANNO_HW  <- 0.12    # annotation strip half-width
S_X_PW   <- 5.2     # pathway bar center
S_PW_HW  <- 0.12    # pathway bar half-width
S_X_BAR  <- S_X_PW + S_PW_HW + 0.10   # enrichment bars start
S_MAX_LEN <- 3.4    # max bar length (wider bars for gene counts)

# --- Heatmap tile data ---
heatmap_df <- class_df %>%
  filter(gene %in% gene_order) %>%
  left_join(gene_bar_df %>% select(gene, y_ctr, ymin, ymax), by = "gene") %>%
  select(gene, y_ctr, ymin, ymax, logFC_Aging, logFC_Training_Old) %>%
  pivot_longer(cols = c(logFC_Aging, logFC_Training_Old),
               names_to = "contrast", values_to = "logFC") %>%
  mutate(x = ifelse(contrast == "logFC_Aging", HM_X_A, HM_X_TO))

fc_max <- max(abs(heatmap_df$logFC), na.rm = TRUE)

# --- Annotation strip data (one rect per gene, colored by class) ---
anno_rects <- gene_bar_df %>%
  mutate(fill = F3_CLASS_COLORS[as.character(f3_class)])

# --- Category group positions (for Sankey source — perfectly aligned) ---
cat_heatmap_pos <- gene_bar_df %>%
  group_by(f3_class) %>%
  summarise(y_bot = min(ymin), y_top = max(ymax), .groups = "drop") %>%
  mutate(y_ctr = (y_top + y_bot) / 2, bar_h = y_top - y_bot)

# Group outline rectangles (span heatmap + annotation strip)
group_outlines <- cat_heatmap_pos %>%
  mutate(
    xmin = HM_X_A - HM_HW,
    xmax = ANNO_X + ANNO_HW,
    outline_col = F3_CLASS_COLORS[as.character(f3_class)]
  )

# Shared y-limits
S_Y_MIN <- min(gene_bar_df$ymin) - gene_h * 2
S_Y_MAX <- max(gene_bar_df$ymax) + gene_h * 2
S_Y_SPAN <- S_Y_MAX - S_Y_MIN

# --- Sankey data ---
cat_pw_counts <- links_1to1 %>%
  filter(pathway != "Other") %>%
  count(f3_class, pathway, name = "n_proteins") %>%
  filter(n_proteins > 0)

cat_totals <- links_1to1 %>%
  filter(pathway != "Other") %>%
  count(f3_class, name = "total") %>%
  arrange(match(as.character(f3_class), F3_CLASS_ORDER))

pw_totals <- links_1to1 %>%
  filter(pathway != "Other") %>%
  count(pathway, name = "total") %>%
  mutate(pw_idx = match(pathway, pw_order)) %>%
  filter(!is.na(pw_idx)) %>%
  arrange(pw_idx)

n_cats <- nrow(cat_totals)
n_pws  <- nrow(pw_totals)

# Category bars (from heatmap gene group positions — exact alignment)
cat_bars <- cat_totals %>%
  left_join(cat_heatmap_pos, by = "f3_class") %>%
  mutate(fill = F3_CLASS_COLORS[as.character(f3_class)])

# Pathway bars (distributed within y-range, proportional to gene count)
pw_gap_frac <- 0.03
pw_usable   <- S_Y_SPAN * (1 - pw_gap_frac * max(n_pws - 1, 0) / max(n_pws, 1))
pw_gap_size <- if (n_pws > 1) (S_Y_SPAN - pw_usable) / (n_pws - 1) else 0

pw_bars <- pw_totals %>%
  mutate(
    bar_h = pw_usable / n_pws,
    y_top = S_Y_MAX - cumsum(c(0, head(bar_h + pw_gap_size, -1))),
    y_bot = y_top - bar_h,
    y_ctr = (y_top + y_bot) / 2
  )

# Pathway colors (from shared CONSOLIDATED_COLORS palette)
PW_COLORS <- CONSOLIDATED_COLORS[pw_totals$pathway]
PW_COLORS[is.na(PW_COLORS)] <- "#D0D0D0"
pw_bars$fill <- PW_COLORS

# Gene count per consolidated pathway (for bar scaling)
pw_gene_counts <- links_1to1 %>%
  filter(pathway != "Other") %>%
  count(pathway, name = "gene_count")

# Keep fold enrichment for reference/CSV export, but bars show gene counts
fg_total <- length(sig_genes)
bg_total <- nrow(all_super)

pw_enrichment <- ora_top %>%
  filter(pathway_label != "Other") %>%
  mutate(
    fg_in = sapply(pathway_label, function(s) sum(links_1to1$pathway == s)),
    bg_in = sapply(pathway_label, function(s) sum(all_super$super == s)),
    fold_enrichment = (fg_in / fg_total) / pmax(bg_in / bg_total, 1e-10)
  ) %>%
  select(pathway_label, p.adjust, fold_enrichment)

ENRICH_AXIS_LABEL <- "Gene count"

pw_bars <- pw_bars %>%
  left_join(pw_enrichment, by = c("pathway" = "pathway_label")) %>%
  mutate(fold_enrichment = replace_na(fold_enrichment, 0)) %>%
  left_join(pw_gene_counts, by = "pathway") %>%
  mutate(gene_count = replace_na(gene_count, 0L)) %>%
  mutate(pw_label = stringr::str_wrap(pathway, width = 18))

pw_label_size <- TXT_PF

# --- Sigmoid ribbons (from annotation strip to pathway bars) ---
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

cat_cum <- setNames(rep(0, n_cats), as.character(cat_totals$f3_class))
pw_cum  <- setNames(rep(0, n_pws), pw_totals$pathway)

ribbon_list <- list()
ribbon_idx  <- 0

for (ct in as.character(cat_totals$f3_class)) {
  ct_contribs <- cat_pw_counts %>%
    filter(as.character(f3_class) == ct) %>%
    arrange(match(pathway, pw_totals$pathway))

  for (r in seq_len(nrow(ct_contribs))) {
    pw <- ct_contribs$pathway[r]
    n  <- ct_contribs$n_proteins[r]
    if (n == 0) next

    ribbon_idx <- ribbon_idx + 1

    ct_row <- cat_bars %>% filter(as.character(f3_class) == ct)
    ct_n   <- ct_row$total
    ct_h   <- ct_row$bar_h
    frac_ct <- n / ct_n
    y0_top <- ct_row$y_top - cat_cum[ct] * ct_h
    y0_bot <- y0_top - frac_ct * ct_h
    cat_cum[ct] <- cat_cum[ct] + frac_ct

    pw_row <- pw_bars %>% filter(pathway == pw)
    pw_n   <- pw_row$total
    pw_h   <- pw_row$bar_h
    frac_pw <- n / pw_n
    y1_top <- pw_row$y_top - pw_cum[pw] * pw_h
    y1_bot <- y1_top - frac_pw * pw_h
    pw_cum[pw] <- pw_cum[pw] + frac_pw

    ribbon_poly <- make_sigmoid_ribbon(
      x0 = ANNO_X + ANNO_HW, x1 = S_X_PW - S_PW_HW,
      y0_top = y0_top, y0_bot = y0_bot,
      y1_top = y1_top, y1_bot = y1_bot,
      ribbon_id = paste0(ct, "->", pw)
    ) %>%
      mutate(f3_class = ct, pathway = pw,
             fill_col = F3_CLASS_COLORS[ct])

    ribbon_list[[ribbon_idx]] <- ribbon_poly
  }
}

ribbons_df <- bind_rows(ribbon_list)
message(sprintf("  Built %d Sankey ribbons", ribbon_idx))

# --- Stacked composition bars (length ∝ √gene count, fill = class composition) ---
# Sqrt scaling: compresses dominant categories (Metabolism) so smaller ones
# remain visually comparable. Axis shows real counts at sqrt-scaled positions.
max_count  <- max(pw_bars$gene_count, na.rm = TRUE)
sqrt_scale <- function(x) sqrt(x) / sqrt(max_count)   # 0-1 scaled
S_SBAR_H   <- min(pw_bars$bar_h[1] * 0.88, 0.90 * S_Y_SPAN / max(n_pws, 1))

# --- Grouped side-by-side bars (one sub-bar per category per pathway) ---
# Each pathway gets N thin bars (one per class present), lengths
# proportional to per-class gene count (sqrt scale), ordered top-to-bottom
# matching the heatmap class order.
max_cats_per_pw <- max(cat_pw_counts %>% count(pathway) %>% pull(n))
sub_h <- S_SBAR_H / max_cats_per_pw

grouped_rects <- list()
grouped_idx   <- 0

active_pw_order <- pw_totals$pathway

for (pw in active_pw_order) {
  pw_row <- pw_bars %>% filter(pathway == pw)
  if (is.na(pw_row$gene_count) || pw_row$gene_count == 0) next

  pw_contribs <- cat_pw_counts %>%
    filter(pathway == pw) %>%
    arrange(match(as.character(f3_class), F3_CLASS_ORDER))

  n_sub <- nrow(pw_contribs)
  total_sub_h <- n_sub * sub_h
  y_start <- pw_row$y_ctr + total_sub_h / 2

  for (r in seq_len(n_sub)) {
    ct <- as.character(pw_contribs$f3_class[r])
    n  <- pw_contribs$n_proteins[r]
    if (n == 0) next

    grouped_idx <- grouped_idx + 1
    bar_w <- sqrt_scale(n) * S_MAX_LEN

    y_top <- y_start - (r - 1) * sub_h
    y_bot <- y_top - sub_h

    grouped_rects[[grouped_idx]] <- tibble(
      xmin = S_X_BAR, xmax = S_X_BAR + bar_w,
      ymin = y_bot, ymax = y_top,
      f3_class = ct, pathway = pw,
      fill = F3_CLASS_COLORS[ct],
      n_genes = n
    )
  }
}

grouped_df <- bind_rows(grouped_rects)

count_labels <- grouped_df %>%
  group_by(pathway) %>%
  summarise(bar_end = max(xmax), n_total = sum(n_genes), .groups = "drop") %>%
  left_join(pw_bars %>% select(pathway, y_ctr), by = "pathway") %>%
  mutate(label = paste0("n=", n_total), x = bar_end + 0.10)

# --- Inline group labels (left of heatmap, replacing embedded key) ---
key_cat_counts <- links_1to1 %>%
  count(f3_class) %>%
  deframe()

inline_label_df <- cat_heatmap_pos %>%
  mutate(
    x     = HM_X_A - HM_HW - 0.15,
    label = sprintf("%s\n(n=%d)", as.character(f3_class),
                    ifelse(as.character(f3_class) %in% names(key_cat_counts),
                           key_cat_counts[as.character(f3_class)], 0L)),
    color = F3_CLASS_COLORS[as.character(f3_class)]
  )

# Gene count axis (sqrt-scaled positions, real count labels)
enrich_ticks <- pretty(c(0, max_count), n = 5)
enrich_ticks <- enrich_ticks[enrich_ticks >= 0 & enrich_ticks <= max_count * 1.05]
if (!0 %in% enrich_ticks) enrich_ticks <- c(0, enrich_ticks)
enrich_ticks <- unique(round(enrich_ticks))  # integer ticks
grid_x <- S_X_BAR + sqrt_scale(enrich_ticks) * S_MAX_LEN
grid_df <- tibble(x = grid_x, label = as.character(as.integer(enrich_ticks)))

axis_y <- S_Y_MIN + 1
tick_len <- 1.5

# Column headers (below heatmap)
header_y <- S_Y_MIN - gene_h * 0.5
sig_header_x <- S_X_BAR + S_MAX_LEN / 2

message(sprintf("  Built %d grouped bar segments", grouped_idx))

# --- Build unified Panel F plot ---
pF <- ggplot() +
  # ---- Heatmap tiles (first fill scale: diverging logFC) ----
  geom_rect(data = heatmap_df,
            aes(xmin = x - HM_HW, xmax = x + HM_HW,
                ymin = ymin, ymax = ymax, fill = logFC),
            color = NA, linewidth = 0) +
  scale_fill_gradient2(
    low = "#2166AC", mid = "white", high = "#B2182B",
    midpoint = 0, limits = c(-fc_max, fc_max),
    name = expression(log[2]~FC),
    guide = guide_colorbar(direction = "horizontal", title.position = "top",
                           barwidth = unit(40, "mm"), barheight = unit(6, "mm"),
                           title.theme = element_text(size = 11, face = "bold"),
                           label.theme = element_text(size = 9))
  ) +
  # ---- Switch to identity fill for all remaining layers ----
  ggnewscale::new_scale_fill() +
  # Annotation strip (class color per gene, directly adjacent to heatmap)
  geom_rect(data = anno_rects,
            aes(xmin = ANNO_X - ANNO_HW, xmax = ANNO_X + ANNO_HW,
                ymin = ymin, ymax = ymax, fill = fill),
            color = NA, linewidth = 0) +
  # Group outlines (span heatmap + annotation strip for visual distinction)
  geom_rect(data = group_outlines,
            aes(xmin = xmin, xmax = xmax, ymin = y_bot, ymax = y_top),
            fill = NA, color = "black", linewidth = 0.5) +
  # Ribbons (from annotation strip edge to pathway bars; lower alpha for ~250 genes)
  geom_polygon(data = ribbons_df,
               aes(x = x, y = y, group = ribbon_id, fill = fill_col),
               alpha = 0.22, color = NA) +
  # Pathway bars
  geom_rect(data = pw_bars,
            aes(xmin = S_X_PW - S_PW_HW, xmax = S_X_PW + S_PW_HW,
                ymin = y_bot, ymax = y_top, fill = fill),
            color = "black", linewidth = 0.3) +
  # Pathway labels (left of pathway bars)
  geom_text(data = pw_bars,
            aes(x = S_X_PW - S_PW_HW - 0.05, y = y_ctr, label = pw_label),
            hjust = 1, size = pw_label_size, fontface = "bold",
            color = "grey20", lineheight = 0.85) +
  # Enrichment bars (grouped by class, length = gene count)
  geom_rect(data = grouped_df,
            aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, fill = fill),
            color = "black", linewidth = 0.2) +
  # Count labels
  geom_text(data = count_labels,
            aes(x = x, y = y_ctr, label = label),
            hjust = 0, size = TXT_PF, fontface = "bold", color = "grey30") +
  # Inline group labels (left of heatmap, colored by class)
  geom_text(data = inline_label_df,
            aes(x = x, y = y_ctr, label = label, color = color),
            hjust = 1, size = TXT_PF * 0.85, fontface = "bold",
            lineheight = 0.85, show.legend = FALSE) +
  scale_color_identity() +
  scale_fill_identity() +
  # ---- Column headers ----
  annotate("text", x = HM_X_A, y = header_y, label = "Aging",
           size = TXT_PF, fontface = "bold") +
  annotate("text", x = HM_X_TO, y = header_y, label = "Trn. Old",
           size = TXT_PF, fontface = "bold") +
  # Gene count axis label at bottom of bars
  annotate("text", x = sig_header_x, y = axis_y - tick_len - 5.0,
           label = ENRICH_AXIS_LABEL,
           size = TXT_PF, fontface = "bold", color = "grey20") +
  # Gene count axis
  annotate("segment", x = S_X_BAR, xend = S_X_BAR + S_MAX_LEN,
           y = axis_y, yend = axis_y, color = "black", linewidth = 0.5) +
  geom_segment(data = grid_df, aes(x = x, xend = x),
               y = axis_y, yend = axis_y - tick_len,
               color = "black", linewidth = 0.4) +
  geom_text(data = grid_df, aes(x = x, label = label),
            y = axis_y - tick_len - 1.5, size = TXT_PF, color = "grey30") +
  # ---- Scales and theme ----
  scale_y_continuous(
    breaks = NULL,
    limits = c(axis_y - tick_len - 7.5, S_Y_MAX + 2),
    expand = c(0, 0)
  ) +
  scale_x_continuous(
    limits = c(HM_X_A - HM_HW - 3.0, S_X_BAR + S_MAX_LEN + 1.5),
    expand = c(0, 0)
  ) +
  labs(x = NULL, y = NULL,
       title = "Aging & Reversal Proteins: Functional Classification",
       subtitle = sprintf("%d proteins | 4 classes | 12 consolidated pathways | Gene count (sqrt scale)",
                          n_class)) +
  THEME_FIG +
  theme(
    axis.text.y      = element_blank(),
    axis.text.x      = element_blank(),
    axis.ticks.x     = element_blank(),
    axis.ticks.y     = element_blank(),
    panel.grid       = element_blank(),
    panel.border     = element_blank(),
    legend.position  = "none",
    plot.margin      = margin(8, 5, 8, 1)
  )

# ==============================================================================
# 6. EXPORT CSVs
# ==============================================================================

message("  Exporting Panel F data CSVs...")

# Classification
class_df %>%
  transmute(
    gene,
    f3_class = as.character(f3_class),
    logFC_Aging = round(logFC_Aging, 4),
    logFC_Training_Old = round(logFC_Training_Old, 4),
    pi_score_Aging = round(pi_score_Aging, 6),
    pi_score_Reversal = round(pi_score_Reversal, 6)
  ) %>%
  arrange(f3_class, desc(abs(logFC_Aging))) %>%
  write_csv(file.path(DAT_DIR, "panel_F", "classification.csv"))

# Sankey links
links_1to1 %>%
  transmute(gene, pathway, database,
            classification = as.character(f3_class),
            logFC_Aging = round(logFC_Aging, 4),
            logFC_Training_Old = round(logFC_Training_Old, 4),
            pathway_pvalue = signif(pvalue, 4)) %>%
  arrange(classification, pathway, gene) %>%
  write_csv(file.path(DAT_DIR, "panel_F", "sankey_links.csv"))

# Enrichment bars
pw_gene_list <- links_1to1 %>%
  group_by(pathway) %>%
  summarise(genes = paste(sort(gene), collapse = "/"), .groups = "drop")

ora_top %>%
  left_join(pw_gene_list, by = c("pathway_label" = "pathway")) %>%
  left_join(pw_1to1_counts, by = c("pathway_label" = "pathway")) %>%
  left_join(pw_enrichment %>% select(pathway_label, fold_enrichment),
            by = "pathway_label") %>%
  transmute(pathway = pathway_label, database, gene_count = n_1to1,
            fold_enrichment = round(fold_enrichment, 3),
            pvalue = signif(pvalue, 4),
            p_adjust = signif(p.adjust, 4), genes) %>%
  arrange(pvalue) %>%
  write_csv(file.path(DAT_DIR, "panel_F", "enrichment_bars.csv"))

message("  CSVs exported")

# ==============================================================================
# 7. SAVE PANEL F
# ==============================================================================

message("  Saving Panel F...")

f_height <- max(300, n_genes * 0.8 + 80)

ggsave(file.path(RPT_DIR, "panel_F_classification.pdf"), pF,
       width = 360, height = f_height, units = "mm", device = cairo_pdf, limitsize = FALSE)
ggsave(file.path(RPT_DIR, "panel_F_classification.png"), pF,
       width = 360, height = f_height, units = "mm", dpi = 300, limitsize = FALSE)

message("  Panel F saved: ", file.path(RPT_DIR, "panel_F_classification.pdf"))

# ==============================================================================
# 8. VERIFICATION SUMMARY
# ==============================================================================

message("\n  === Panel F Verification Summary ===")
message(sprintf("  Proteins per class:"))
for (ct in F3_CLASS_ORDER) {
  n_ct <- sum(links_1to1$f3_class == ct, na.rm = TRUE)
  if (n_ct > 0) message(sprintf("    %s: %d", ct, n_ct))
}
message(sprintf("  Consolidated pathways: %d (excl. Other)", length(active_pw_order)))
message(sprintf("  Pathway membership by class:"))
for (pw in active_pw_order) {
  pw_genes <- links_1to1 %>% filter(pathway == pw) %>% nrow()
  pw_padj <- pw_bars %>% filter(pathway == pw) %>% pull(p.adjust)
  message(sprintf("    %s: n=%d, p.adj=%.2e", pw, pw_genes,
                  ifelse(length(pw_padj) > 0, pw_padj, NA)))
}
sankey_cat_totals <- links_1to1 %>%
  filter(pathway != "Other") %>%
  count(f3_class) %>%
  pull(n) %>% sum()
message(sprintf("  Sankey source total: %d (should match mapped non-Other: %d) %s",
                sankey_cat_totals, sum(pw_totals$total),
                ifelse(sankey_cat_totals == sum(pw_totals$total), "OK", "MISMATCH")))
message("  Panel F complete")
