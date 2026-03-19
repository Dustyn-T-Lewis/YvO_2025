# Figure 3 — Panel B: UpSet Plot (Dual Bar + Dot Matrix)
# Sources panel_A for sig_sets. Outputs: pB_bars, pB_dots (ggplot objects)

setwd(rprojroot::find_rstudio_root_file())
source("04_Figures/F03/a_script/style.R")

suppressPackageStartupMessages({
  library(dplyr)
  library(tibble)
  library(tidyr)
  library(patchwork)
  library(ComplexHeatmap)
  library(purrr)
})

# Depends on Panel A for sig_sets, dir_map, all_genes, etc.
if (!exists("sig_sets")) source("04_Figures/F03/a_script/panel_A.R")

RPT <- "04_Figures/F03/b_reports"
DAT <- "04_Figures/F03/c_data"
dir.create(RPT, recursive = TRUE, showWarnings = FALSE)
dir.create(DAT, recursive = TRUE, showWarnings = FALSE)

pdf_device <- get_pdf_device()
PB_W <- 170

bin_mat <- sapply(sig_sets, function(s) as.integer(all_genes %in% s))
rownames(bin_mat) <- all_genes
colnames(bin_mat) <- SET_LABELS[colnames(bin_mat)]
label_to_contrast <- setNames(names(SET_LABELS), SET_LABELS)

cm <- make_comb_mat(bin_mat, mode = "intersect")
cs <- comb_size(cm)
keep <- cs > 0 & comb_degree(cm) > 0
cm_sub <- cm[keep]

comb_names_vec    <- comb_name(cm_sub)
n_comb            <- length(comb_names_vec)
set_names_ordered <- set_name(cm_sub)
up_counts    <- numeric(n_comb)
down_counts  <- numeric(n_comb)
mixed_counts <- numeric(n_comb)

for (i in seq_len(n_comb)) {
  members <- extract_comb(cm_sub, comb_names_vec[i])
  if (length(members) == 0) next
  bits <- as.logical(as.integer(strsplit(comb_names_vec[i], "")[[1]]))
  active_contrasts <- label_to_contrast[set_names_ordered[bits]]
  gene_dirs <- sapply(members, function(g) {
    dirs <- sapply(active_contrasts, function(ctr) {
      if (g %in% names(dir_map[[ctr]])) dir_map[[ctr]][g] else NA
    })
    dirs <- dirs[!is.na(dirs)]
    if (all(dirs == "Up")) "Up" else if (all(dirs == "Down")) "Down" else "Mixed"
  })
  up_counts[i]    <- sum(gene_dirs == "Up")
  down_counts[i]  <- sum(gene_dirs == "Down")
  mixed_counts[i] <- sum(gene_dirs == "Mixed")
}
n_mixed <- sum(mixed_counts)

# --- Pairwise overlap enrichment (Fisher's exact, one-sided) ---
overlap_tests <- list()
contrast_pairs <- combn(CONTRASTS, 2, simplify = FALSE)
n_bg <- length(all_genes)
for (pair in contrast_pairs) {
  a <- sig_sets[[pair[1]]]; b <- sig_sets[[pair[2]]]
  n_both <- length(intersect(a, b))
  n_a    <- length(a); n_b <- length(b)
  mat <- matrix(c(n_both, n_a - n_both, n_b - n_both,
                  n_bg - n_a - n_b + n_both), nrow = 2)
  ft <- fisher.test(mat, alternative = "greater")
  overlap_tests[[paste(pair, collapse = " & ")]] <- data.frame(
    set_A = pair[1], set_B = pair[2],
    n_A = n_a, n_B = n_b, overlap = n_both,
    expected = round(n_a * n_b / n_bg, 1),
    odds_ratio = round(ft$estimate, 2), p_value = ft$p.value
  )
}
overlap_df <- bind_rows(overlap_tests)
overlap_df$p_bh <- p.adjust(overlap_df$p_value, method = "BH")
n_sig_overlaps <- sum(overlap_df$p_bh < 0.05)
write.csv(overlap_df, file.path(DAT, "audit_panel_B_overlap_enrichment.csv"), row.names = FALSE)

display_total <- up_counts + down_counts
keep_display  <- display_total > 0
comb_ord <- which(keep_display)[order(-display_total[keep_display])]
up_ord   <- up_counts[comb_ord]
down_ord <- down_counts[comb_ord]

set_display_order <- c("Aging", "Tr.(Y)", "Tr.(O)", "Tr.(O)\u2013Tr.(Y)")
n_int          <- length(comb_ord)
comb_names_ord <- comb_names_vec[comb_ord]
set_order_ch   <- set_name(cm_sub)
set_y_levels   <- rev(set_display_order)

comb_deg_ord <- vapply(comb_names_ord, function(cn) {
  sum(as.integer(strsplit(cn, "")[[1]]))
}, integer(1))

bar_long <- tibble(
  x         = rep(seq_len(n_int), 2),
  direction = factor(rep(c("Up", "Down"), each = n_int), levels = c("Down", "Up")),
  count     = c(up_ord, down_ord),
  is_single = rep(comb_deg_ord == 1, 2)
)

dot_df <- expand_grid(x = seq_len(n_int), set = set_display_order)
dot_df$active <- vapply(seq_len(nrow(dot_df)), function(r) {
  bits <- as.integer(strsplit(comb_names_ord[dot_df$x[r]], "")[[1]])
  as.logical(bits[match(dot_df$set[r], set_order_ch)])
}, logical(1))
dot_df$set  <- factor(dot_df$set, levels = set_y_levels)
dot_df$ynum <- as.numeric(dot_df$set)

seg_list <- vector("list", n_int)
for (i in seq_len(n_int)) {
  bits <- as.logical(as.integer(strsplit(comb_names_ord[i], "")[[1]]))
  ypos <- match(set_order_ch[bits], set_y_levels)
  if (length(ypos) > 1)
    seg_list[[i]] <- tibble(x = i, ymin = min(ypos), ymax = max(ypos))
}
seg_df <- bind_rows(seg_list)

stripe_fills <- adjustcolor(unname(SET_DISPLAY_COLORS[set_y_levels]), alpha.f = 0.20)

bar_bg_list <- lapply(seq_len(n_int), function(i) {
  bits <- as.logical(as.integer(strsplit(comb_names_ord[i], "")[[1]]))
  if (sum(bits) != 1) return(NULL)
  tibble(xmin = i - 0.5, xmax = i + 0.5,
         fill = unname(SET_DISPLAY_COLORS[set_order_ch[bits]]))
})
bar_bg <- bind_rows(compact(bar_bg_list))

lbl_sz <- scale_text(BASE_COUNT, PB_W)

# Build overlap annotation text (top significant pairwise overlaps)
sig_overlaps <- overlap_df |>
  filter(p_bh < 0.05) |>
  arrange(p_bh) |>
  head(4) |>
  mutate(
    set_A_short = CTR_SHORT[set_A],
    set_B_short = CTR_SHORT[set_B],
    label = sprintf("%s x %s: OR=%.1f, %s",
                    set_A_short, set_B_short,
                    odds_ratio, sapply(p_bh, fmt_p))
  )

if (nrow(sig_overlaps) > 0) {
  overlap_annotation <- paste(sig_overlaps$label, collapse = "\n")
  top_or_txt <- sprintf("Strong %s<->%s convergence (OR=%.1f)",
                        sig_overlaps$set_A_short[1],
                        sig_overlaps$set_B_short[1],
                        sig_overlaps$odds_ratio[1])
} else {
  overlap_annotation <- "No significant pairwise overlaps"
  top_or_txt <- ""
}

pB_bars <- ggplot(bar_long, aes(x, count, fill = direction)) +
  geom_rect(data = bar_bg,
            aes(xmin = xmin, xmax = xmax, ymin = -Inf, ymax = Inf),
            fill = bar_bg$fill, alpha = 0.20,
            color = "grey70", linewidth = 0.2, inherit.aes = FALSE) +
  geom_col(position = position_dodge(width = 0.7), width = 0.6,
           color = "black", linewidth = 0.3) +
  geom_text(data = \(d) d |> filter(count > 0, is_single),
            aes(label = count, y = count / 2),
            position = position_dodge(width = 0.7), vjust = 0.5,
            size = lbl_sz - 0.8, color = "white", fontface = "bold") +
  geom_text(data = \(d) d |> filter(count > 0, !is_single),
            aes(label = count, y = count + 2),
            position = position_dodge(width = 0.7), vjust = 0,
            size = lbl_sz - 0.8, color = "black", fontface = "bold") +
  scale_fill_manual(values = c(Up = unname(DIR_COLORS["Up"]),
                                Down = unname(DIR_COLORS["Down"]))) +
  scale_x_continuous(expand = expansion(add = 0.3)) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.12))) +
  labs(y = "Intersection\nsize",
       title = "Contrast Overlap (UpSet)",
       subtitle = sprintf(
         "\u03A0 < 0.05 DEPs | %d/%d overlaps enriched | %s\n%s",
         n_sig_overlaps, nrow(overlap_df), top_or_txt,
         overlap_annotation
       ),
       tag = "B") +
  FIG_THEME +
  theme(axis.text.x  = element_blank(), axis.ticks.x = element_blank(),
        axis.title.x = element_blank(),
        panel.grid.major.x = element_blank(), panel.grid.minor = element_blank(),
        legend.position = "none",
        plot.margin = margin(2, 5.5, 0, 5.5))

pB_dots <- ggplot() +
  annotate("rect",
           xmin = 0.4, xmax = n_int + 0.6,
           ymin = seq_along(set_y_levels) - 0.40,
           ymax = seq_along(set_y_levels) + 0.40,
           fill = stripe_fills) +
  {if (nrow(seg_df) > 0)
    geom_segment(data = seg_df,
                 aes(x = x, xend = x, y = ymin, yend = ymax),
                 linewidth = 0.6, color = "grey25")} +
  geom_point(data = dot_df |> filter(!active),
             aes(x = x, y = ynum), color = "grey78", size = 1.6) +
  geom_point(data = dot_df |> filter(active),
             aes(x = x, y = ynum), color = "grey15", size = 1.6) +
  scale_x_continuous(expand = expansion(add = 0.3)) +
  scale_y_continuous(breaks = seq_along(set_y_levels), labels = set_y_levels,
                     expand = expansion(add = 0)) +
  coord_cartesian(ylim = c(0.55, 4.45)) +
  labs(x = NULL, y = NULL) +
  FIG_THEME +
  theme(axis.text.x  = element_blank(), axis.ticks = element_blank(),
        panel.grid   = element_blank(),
        panel.border = element_rect(color = "grey70", fill = NA, linewidth = 0.3),
        axis.text.y  = element_text(size = FIG_AXIS_TEXT - 1, face = "bold",
                                     margin = margin(r = 1)),
        plot.margin  = margin(0, 5.5, 0, 5.5))

pB_combined <- (pB_bars / pB_dots) + plot_layout(heights = c(0.80, 0.20))

ggsave(file.path(RPT, "panel_B_upset.pdf"), pB_combined,
       width = PB_W, height = 120, units = "mm", device = pdf_device)
ggsave(file.path(RPT, "panel_B_upset.png"), pB_combined,
       width = PB_W, height = 120, units = "mm", dpi = 300)
