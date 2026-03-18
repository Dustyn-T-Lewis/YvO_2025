# F3 Panel F: Aging & Reversal Classification
setwd(rprojroot::find_rstudio_root_file())
source("04_Figures/shared/style.R")
source("04_Figures/shared/pathway_utils.R")
source("04_Figures/shared/go_slim_categories.R")

suppressPackageStartupMessages({
  library(tidyverse)
  library(ggnewscale)
  library(stringr)
})

PW <- 360; PH <- 300
RPT <- "04_Figures/F05/b_reports"
DAT <- "04_Figures/F05/c_data"
dir.create(RPT, recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(DAT, "panel_F"), recursive = TRUE, showWarnings = FALSE)

dep_df <- read_csv("03_DEP/c_data/03_combined_results.csv", show_col_types = FALSE)

pdf_device <- get_pdf_device()
PF_W <- 360

F3_CLASS_ORDER <- c("Aging Up", "Aging Down", "Rev. (Age Up)", "Rev. (Age Down)")

F3_CLASS_COLORS <- c(
  "Aging Up"        = "#D6604D",
  "Aging Down"      = "#4393C3",
  "Rev. (Age Up)"   = "#E65100",
  "Rev. (Age Down)" = "#FFB300"
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

sig_genes  <- class_df$gene
all_genes  <- dep_df$gene

# GO Slim-based classification: maps each gene to 1 of 15 biologically
# coherent categories via GO:BP annotation + GOBPANCESTOR hierarchy.
go_slim_result <- assign_go_slim_consolidated(sig_genes, all_genes)

links_1to1 <- go_slim_result %>%
  rename(pathway = consolidated) %>%
  mutate(pathway = as.character(pathway)) %>%
  left_join(class_df %>% select(gene, f3_class, logFC_Aging, logFC_Training_Old),
            by = "gene")

pw_1to1_counts <- links_1to1 %>% count(pathway, name = "n_1to1")
pw_order_df <- pw_1to1_counts %>% arrange(desc(n_1to1))
pw_order <- pw_order_df$pathway
if ("Other" %in% pw_order) {
  pw_order <- c(setdiff(pw_order, "Other"), "Other")
}

gene_order <- links_1to1 %>%
  mutate(
    class_idx = match(as.character(f3_class), F3_CLASS_ORDER),
    pw_idx    = match(pathway, pw_order)
  ) %>%
  arrange(class_idx, pw_idx, logFC_Aging) %>%
  pull(gene)

txt_pf <- scale_text(BASE_STAT, PF_W)

n_genes <- length(gene_order)
Y_SPAN  <- n_genes

n_cat_breaks <- length(unique(class_df$f3_class)) - 1
gap_size     <- 1.2
total_gap    <- n_cat_breaks * gap_size
gene_h       <- (Y_SPAN - total_gap) / n_genes

gene_bar_df <- tibble(gene = gene_order, idx = seq_along(gene_order)) %>%
  mutate(
    ymax  = Y_SPAN - (idx - 1) * gene_h,
    ymin  = ymax - gene_h,
    y_ctr = (ymin + ymax) / 2
  ) %>%
  left_join(links_1to1 %>% select(gene, f3_class) %>% distinct(), by = "gene")

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

HM_X_A   <- 1
HM_X_TO  <- 1.92
HM_HW    <- 0.45
ANNO_X   <- HM_X_TO + HM_HW + 0.01 + 0.12
ANNO_HW  <- 0.12
S_X_PW   <- 5.2
S_PW_HW  <- 0.12
S_X_BAR  <- S_X_PW + S_PW_HW + 0.10
S_MAX_LEN <- 3.4

heatmap_df <- class_df %>%
  filter(gene %in% gene_order) %>%
  left_join(gene_bar_df %>% select(gene, y_ctr, ymin, ymax), by = "gene") %>%
  select(gene, y_ctr, ymin, ymax, logFC_Aging, logFC_Training_Old) %>%
  pivot_longer(cols = c(logFC_Aging, logFC_Training_Old),
               names_to = "contrast", values_to = "logFC") %>%
  mutate(x = ifelse(contrast == "logFC_Aging", HM_X_A, HM_X_TO))

fc_max <- max(abs(heatmap_df$logFC), na.rm = TRUE)

anno_rects <- gene_bar_df %>%
  mutate(fill = F3_CLASS_COLORS[as.character(f3_class)])

cat_heatmap_pos <- gene_bar_df %>%
  group_by(f3_class) %>%
  summarise(y_bot = min(ymin), y_top = max(ymax), .groups = "drop") %>%
  mutate(y_ctr = (y_top + y_bot) / 2, bar_h = y_top - y_bot)

group_outlines <- cat_heatmap_pos %>%
  mutate(
    xmin = HM_X_A - HM_HW,
    xmax = ANNO_X + ANNO_HW,
    outline_col = F3_CLASS_COLORS[as.character(f3_class)]
  )

S_Y_MIN <- min(gene_bar_df$ymin) - gene_h * 2
S_Y_MAX <- max(gene_bar_df$ymax) + gene_h * 2
S_Y_SPAN <- S_Y_MAX - S_Y_MIN

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

cat_bars <- cat_totals %>%
  left_join(cat_heatmap_pos, by = "f3_class") %>%
  mutate(fill = F3_CLASS_COLORS[as.character(f3_class)])

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

pw_bars$fill <- CONSOLIDATED_COLORS[pw_bars$pathway]

pw_gene_counts <- links_1to1 %>%
  filter(pathway != "Other") %>%
  count(pathway, name = "gene_count")

pw_bars <- pw_bars %>%
  left_join(pw_gene_counts, by = "pathway") %>%
  mutate(gene_count = replace_na(gene_count, 0L)) %>%
  mutate(pw_label = stringr::str_wrap(pathway, width = 22))

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

max_count  <- max(pw_bars$gene_count, na.rm = TRUE)
sqrt_scale <- function(x) sqrt(x) / sqrt(max_count)
S_SBAR_H   <- min(pw_bars$bar_h[1] * 0.88, 0.90 * S_Y_SPAN / max(n_pws, 1))

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

enrich_ticks <- pretty(c(0, max_count), n = 5)
enrich_ticks <- enrich_ticks[enrich_ticks >= 0 & enrich_ticks <= max_count * 1.05]
if (!0 %in% enrich_ticks) enrich_ticks <- c(0, enrich_ticks)
enrich_ticks <- unique(round(enrich_ticks))
grid_x <- S_X_BAR + sqrt_scale(enrich_ticks) * S_MAX_LEN
grid_df <- tibble(x = grid_x, label = as.character(as.integer(enrich_ticks)))

axis_y <- S_Y_MIN + 1
tick_len <- 1.5

header_y <- S_Y_MIN - gene_h * 0.5
sig_header_x <- S_X_BAR + S_MAX_LEN / 2

message(sprintf("  Built %d grouped bar segments", grouped_idx))

pF <- ggplot() +
  geom_rect(data = heatmap_df,
            aes(xmin = x - HM_HW, xmax = x + HM_HW,
                ymin = ymin, ymax = ymax, fill = logFC),
            color = NA, linewidth = 0) +
  scale_fill_gradient2(
    low = "#2166AC", mid = "white", high = "#B2182B",
    midpoint = 0, limits = c(-fc_max, fc_max),
    name = expression(log[2]~FC)
  ) +
  ggnewscale::new_scale_fill() +
  geom_rect(data = anno_rects,
            aes(xmin = ANNO_X - ANNO_HW, xmax = ANNO_X + ANNO_HW,
                ymin = ymin, ymax = ymax, fill = fill),
            color = NA, linewidth = 0) +
  geom_rect(data = group_outlines,
            aes(xmin = xmin, xmax = xmax, ymin = y_bot, ymax = y_top),
            fill = NA, color = "black", linewidth = 0.5) +
  geom_polygon(data = ribbons_df,
               aes(x = x, y = y, group = ribbon_id, fill = fill_col),
               alpha = 0.22, color = NA) +
  geom_rect(data = pw_bars,
            aes(xmin = S_X_PW - S_PW_HW, xmax = S_X_PW + S_PW_HW,
                ymin = y_bot, ymax = y_top, fill = fill),
            color = "black", linewidth = 0.3) +
  geom_text(data = pw_bars,
            aes(x = S_X_PW - S_PW_HW - 0.05, y = y_ctr, label = pw_label),
            hjust = 1, size = txt_pf, fontface = "bold",
            color = "grey20", lineheight = 0.85) +
  geom_rect(data = grouped_df,
            aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, fill = fill),
            color = "black", linewidth = 0.2) +
  geom_text(data = count_labels,
            aes(x = x, y = y_ctr, label = label),
            hjust = 0, size = txt_pf, fontface = "bold", color = "grey30") +
  geom_text(data = inline_label_df,
            aes(x = x, y = y_ctr, label = label, color = color),
            hjust = 1, size = txt_pf * 0.85, fontface = "bold",
            lineheight = 0.85, show.legend = FALSE) +
  scale_color_identity() +
  scale_fill_identity() +
  annotate("text", x = HM_X_A, y = header_y, label = "Aging",
           size = txt_pf, fontface = "bold") +
  annotate("text", x = HM_X_TO, y = header_y, label = "Trn. Old",
           size = txt_pf, fontface = "bold") +
  annotate("text", x = sig_header_x, y = axis_y - tick_len - 5.0,
           label = "Gene count",
           size = txt_pf, fontface = "bold", color = "grey20") +
  annotate("segment", x = S_X_BAR, xend = S_X_BAR + S_MAX_LEN,
           y = axis_y, yend = axis_y, color = "black", linewidth = 0.5) +
  geom_segment(data = grid_df, aes(x = x, xend = x),
               y = axis_y, yend = axis_y - tick_len,
               color = "black", linewidth = 0.4) +
  geom_text(data = grid_df, aes(x = x, label = label),
            y = axis_y - tick_len - 1.5, size = txt_pf, color = "grey30") +
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
       subtitle = sprintf("%d proteins | 4 classes | %d GO Slim categories | Gene count (sqrt scale)",
                          n_class, length(active_pw_order))) +
  FIG_THEME +
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
  write_csv(file.path(DAT, "panel_F", "classification.csv"))

links_1to1 %>%
  transmute(gene, pathway, go_slim_id = slim,
            classification = as.character(f3_class),
            logFC_Aging = round(logFC_Aging, 4),
            logFC_Training_Old = round(logFC_Training_Old, 4)) %>%
  arrange(classification, pathway, gene) %>%
  write_csv(file.path(DAT, "panel_F", "sankey_links.csv"))

pw_gene_list <- links_1to1 %>%
  filter(pathway != "Other") %>%
  group_by(pathway) %>%
  summarise(genes = paste(sort(gene), collapse = "/"), .groups = "drop")

pw_1to1_counts %>%
  left_join(pw_gene_list, by = "pathway") %>%
  transmute(pathway, gene_count = n_1to1, genes) %>%
  arrange(desc(gene_count)) %>%
  write_csv(file.path(DAT, "panel_F", "category_counts.csv"))

f_height <- max(300, n_genes * 0.8 + 80)

ggsave(file.path(RPT, "panel_F_classification.pdf"), pF,
       width = PF_W, height = f_height, units = "mm", device = pdf_device,
       limitsize = FALSE)
ggsave(file.path(RPT, "panel_F_classification.png"), pF,
       width = PF_W, height = f_height, units = "mm", dpi = 300,
       limitsize = FALSE)

cat("F3 Panel F done\n")
