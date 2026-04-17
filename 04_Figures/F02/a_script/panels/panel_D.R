# F02 — Panel D: DEP Contrast Overlap (UpSet Plot)
# Adapted from 04_Figures/F03/a_script/panel_C_upset.R
# Custom ggplot2 upset (dual bar + dot matrix) using ComplexHeatmap::make_comb_mat
# Pi-score significant DEPs across 4 contrasts, split by Up/Down direction
# Outputs: pD (ggdraw object), panel_D_upset_MAIN.pdf/.png

setwd(rprojroot::find_rstudio_root_file())
source("04_Figures/F02/a_script/style.R")

suppressPackageStartupMessages({
  library(dplyr)
  library(tibble)
  library(tidyr)
  library(readr)
  library(patchwork)
  library(cowplot)
  library(ComplexHeatmap)
  library(purrr)
})

DEP_FILE <- "03_DEP/c_data/03_combined_results.csv"
RPT      <- "04_Figures/F02/b_reports/panels"
DAT      <- "04_Figures/F02/c_data"
dir.create(RPT, recursive = TRUE, showWarnings = FALSE)
dir.create(DAT, recursive = TRUE, showWarnings = FALSE)

CONTRASTS <- c("Aging", "Training_Young", "Training_Old", "Interaction")
dep_df    <- read_csv(DEP_FILE, show_col_types = FALSE)
pdf_device <- get_pdf_device()

PC_W <- 165
PC_H <- 80

SET_LABELS <- c(Aging = "Aging", Training_Young = "Tr.(Y)",
                Training_Old = "Tr.(O)", Interaction = "Tr.(O)\u2013Tr.(Y)")

SET_DISPLAY_COLORS <- c("Aging"  = unname(CONTRAST_COLORS["Aging"]),
                        "Tr.(Y)" = unname(CONTRAST_COLORS["Training_Young"]),
                        "Tr.(O)" = unname(CONTRAST_COLORS["Training_Old"]),
                        "Tr.(O)\u2013Tr.(Y)" = unname(CONTRAST_COLORS["Interaction"]))

all_genes <- unique(dep_df$gene[!is.na(dep_df$gene)])

# --- Build Pi-score significant sets and direction maps ---
sig_sets <- list()
dir_map  <- list()

for (ctr in CONTRASTS) {
  sig_pi  <- dep_df[[paste0("sig_pi_", ctr)]]
  is_sig  <- !is.na(sig_pi) & sig_pi != 0
  sig_sets[[ctr]] <- dep_df$gene[is_sig]
  dir_map[[ctr]]  <- setNames(ifelse(sig_pi[is_sig] > 0, "Up", "Down"),
                               dep_df$gene[is_sig])
}

# --- Build combination matrix ---
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
write.csv(overlap_df, file.path(DAT, "panel_D_upset_overlap_enrichment.csv"), row.names = FALSE)

# --- Sort intersections by total size ---
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

lbl_sz <- scale_text(BASE_COUNT, PC_W)

# Gray vertical divider position: between bars with (up, down) counts
# (6, 4) and (2, 4), independent of direction orientation. If neither
# pattern is present, divider is not drawn.
divider_x <- {
  is_6_4 <- (up_ord == 6 & down_ord == 4) | (up_ord == 4 & down_ord == 6)
  is_2_4 <- (up_ord == 2 & down_ord == 4) | (up_ord == 4 & down_ord == 2)
  six_four_x <- which(is_6_4)
  two_four_x <- which(is_2_4)
  if (length(six_four_x) > 0 && length(two_four_x) > 0) {
    (max(six_four_x) + min(two_four_x)) / 2
  } else {
    NA_real_
  }
}

# Build overlap annotation text
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

pD_bars <- ggplot(bar_long, aes(x, count, fill = direction)) +
  geom_rect(data = bar_bg,
            aes(xmin = xmin, xmax = xmax, ymin = -Inf, ymax = Inf),
            fill = bar_bg$fill, alpha = 0.20,
            color = "grey70", linewidth = 0.2, inherit.aes = FALSE) +
  # Divider matches pD_bars' major-y gridline shade (grey80, linewidth 0.3);
  # only drawn in pD_bars so it stops at the x-axis and does not extend
  # into the dot matrix below.
  {if (!is.na(divider_x))
    geom_vline(xintercept = divider_x, color = "grey80",
               linewidth = 0.3, inherit.aes = FALSE)} +
  geom_col(position = position_dodge(width = 0.7), width = 0.6,
           color = "black", linewidth = 0.3) +
  geom_text(data = \(d) d |> filter(count > 0, is_single),
            aes(label = count, y = count / 2),
            position = position_dodge(width = 0.7), vjust = 0.5,
            size = lbl_sz - 0.8, color = "white", fontface = "bold") +
  geom_text(data = \(d) d |> filter(count > 0, !is_single),
            aes(label = count, y = count + 2),
            position = position_dodge(width = 0.7), vjust = 0,
            size = lbl_sz - 0.8, color = "black", fontface = "bold",
            check_overlap = TRUE) +
  scale_fill_manual(values = c(Up = unname(DIR_COLORS["Up"]),
                                Down = unname(DIR_COLORS["Down"]))) +
  scale_x_continuous(expand = expansion(add = 0)) +
  scale_y_continuous(expand = expansion(mult = c(0, 0))) +
  coord_cartesian(xlim = c(0.5, n_int + 0.5)) +
  labs(title = "Contrast Overlap (UpSet)",
       subtitle = sprintf("\u03A0 DEPs | %d/%d sig overlaps",
                          n_sig_overlaps, nrow(overlap_df)),
       y = NULL) +
  FIG_THEME +
  theme(plot.subtitle = element_text(size = FIG_SUBTITLE_SIZE,
                                     face = "bold.italic", color = "grey40"),
        axis.text.x  = element_blank(), axis.ticks.x = element_blank(),
        axis.title.x = element_blank(),
        panel.grid.major.x = element_blank(),
        # Darker y-grid so lines remain legible under 0.20-alpha contrast
        # rectangles AND in the unshaded overlap columns
        panel.grid.major.y = element_line(color = "grey80", linewidth = 0.3),
        panel.grid.minor   = element_blank(),
        legend.position = "none",
        # r=8 (inner) + r=2 (outer pD_pw) = 10pt total, matching pA r=10
        # for consistent visual space between the A/D and B/E columns.
        plot.margin = margin(2, 8, 0, 8))

pD_dots <- ggplot() +
  annotate("rect",
           xmin = 0.5, xmax = n_int + 0.5,
           ymin = seq_along(set_y_levels) - 0.45,
           ymax = seq_along(set_y_levels) + 0.45,
           fill = stripe_fills) +
  {if (nrow(seg_df) > 0)
    geom_segment(data = seg_df,
                 aes(x = x, xend = x, y = ymin, yend = ymax),
                 linewidth = 0.6, color = "grey25")} +
  geom_point(data = dot_df |> filter(!active),
             aes(x = x, y = ynum), color = "grey78", size = 1.6) +
  geom_point(data = dot_df |> filter(active),
             aes(x = x, y = ynum), color = "grey15", size = 1.6) +
  scale_x_continuous(expand = expansion(add = 0)) +
  scale_y_continuous(breaks = seq_along(set_y_levels), labels = set_y_levels,
                     expand = expansion(add = 0)) +
  coord_cartesian(xlim = c(0.5, n_int + 0.5), ylim = c(0.50, 4.50)) +
  labs(x = NULL, y = NULL) +
  FIG_THEME +
  theme(axis.text.x  = element_blank(), axis.ticks = element_blank(),
        panel.grid   = element_blank(),
        panel.border = element_rect(color = "grey70", fill = NA, linewidth = 0.3),
        axis.text.y  = element_text(size = FIG_AXIS_TEXT - 0.5, face = "bold",
                                     margin = margin(r = 1)),
        # Match pD_bars r=8 so the dot matrix and bar panel share identical
        # horizontal plot-region extents under the A/D column alignment.
        plot.margin  = margin(0, 8, 0, 8))

pD_pw <- (pD_bars / pD_dots) + plot_layout(heights = c(0.78, 0.22)) +
  plot_annotation(theme = theme(plot.margin = margin(t = 2, r = 2, b = 4, l = 0)))

# --- Direction key: EXACT mirror of panel E's make_key_plot() style
# (same square size = 2.8, same text size = 2.4, no background, theme_void) ---
dir_key_df_D <- tibble(
  label     = c("Direction", "Up", "Down"),
  y         = -cumsum(c(0, 0.004, 0.004)),
  fill      = c(NA, unname(DIR_COLORS["Up"]), unname(DIR_COLORS["Down"])),
  is_header = c(TRUE, FALSE, FALSE)
)

p_key_dir_D <- ggplot(dir_key_df_D |> dplyr::filter(!is_header)) +
  geom_point(aes(x = 0, y = y), shape = 22, size = 2.8,
             fill = dir_key_df_D$fill[!dir_key_df_D$is_header],
             color = "grey30", stroke = 0.4) +
  # Match Panel E's make_key_plot() text size (key_txt = 2.85) and weights
  geom_text(aes(x = 0.35, y = y, label = label),
            size = 2.85, color = "grey20", hjust = 0) +
  geom_text(data = dir_key_df_D |> dplyr::filter(is_header),
            aes(x = -0.1, y = y, label = label),
            size = 2.85, fontface = "bold", color = "grey20", hjust = 0) +
  scale_x_continuous(limits = c(-0.2, 2.0)) +
  # Mirror pE's shared_ylim exactly. pE uses c(min(db_df$y) - 0.005, 0.005)
  # where db_df has 6 rows (Database + 5 DB names), giving c(-0.025, 0.005).
  # Using the same ylim here gives pD's direction key an identical vertical
  # extent, entry spacing, and y-height to pE's direction key (entries sit
  # at the top with empty space below, matching pE's appearance).
  coord_cartesian(ylim = c(-0.025, 0.005), clip = "off") +
  theme_void() +
  theme(plot.margin = margin(0, 0, 0, 0))

# Place y-axis title + direction key manually via cowplot.
# Direction key: upper-right corner of the plot area (separator between
# (6,4) and (2,4) bars is drawn in-plot via geom_vline instead).
pD <- ggdraw(pD_pw) +
  draw_label("Intersection size", x = 0.05, y = 0.58, angle = 90,
             size = 10, fontface = "bold") +
  draw_plot(p_key_dir_D, x = 0.83, y = 0.70, width = 0.14, height = 0.22)

ggsave(file.path(RPT, "MAIN_panel_D_upset.png"), pD,
       width = PC_W, height = PC_H, units = "mm", dpi = 300)

message("F02 Panel D (upset, Pi-score) done")
