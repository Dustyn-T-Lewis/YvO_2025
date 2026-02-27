################################################################################
#   Figure 1 — Panel E: UpSet Plot (Dual Bar + Dot Matrix)
#
#   Requires from setup: dep_df, CONTRASTS, CONTRAST_COLORS, DIR_COLORS,
#                         THEME_PUB, RPT_DIR, KEY_TEXT, KEY_TITLE
#   Requires from panel_D: sig_sets, dir_map, all_genes, SET_LABELS,
#                           SET_DISPLAY_COLORS, pi_total, fdr_total
#   Outputs: pE_bars, pE_dots (ggplot objects), pKeys (legend key)
#
#   NOTE: This panel depends on panel_D.R for sig_sets, dir_map, etc.
#         Source panel_D.R before this script if running standalone.
################################################################################

if (!exists("meta")) source("04_Figures/F1/a_script/YvO_F1_setup.R")
if (!exists("sig_sets")) source("04_Figures/F1/a_script/panel_D.R")

# 1. KEYS =================================================================

# Unified key — Contrast Definitions (left) + Direction (right)
# Single ggplot ensures perfect alignment and consistent box sizing.

ctr_leg <- tibble(
  y     = c(1.35, 0.90, 0.45, 0),
  label = c("Aging:  Old_Pre \u2212 Young_Pre",
            "Tr. (Young):  Young_Post \u2212 Young_Pre",
            "Tr. (Old):  Old_Post \u2212 Old_Pre",
            "Interaction:  (Old_Post\u2212Old_Pre) \u2212 (Young_Post\u2212Young_Pre)"),
  fill  = unname(CONTRAST_COLORS[c("Aging", "Training_Young",
                                    "Training_Old", "Interaction")])
)

dir_leg <- tibble(
  y     = c(0.90, 0.45),
  label = c("Up", "Down"),
  fill  = unname(DIR_COLORS[c("Up", "Down")])
)

# X layout: contrast boxes 0..5.2, direction boxes 5.6..6.5
CTR_XMAX <- 5.2    # wide enough for the Interaction text
DIR_XMIN <- 5.6    # small gap between sections
DIR_XMAX <- 6.5    # direction boxes narrower (only "Up"/"Down")
BOX_HALF <- 0.18   # half-height of each box (same for both key types)

pKeys <- ggplot() +
  # -- Contrast definition boxes (left) --
  geom_rect(data = ctr_leg,
            aes(xmin = 0, xmax = CTR_XMAX, ymin = y - BOX_HALF, ymax = y + BOX_HALF),
            fill = alpha(ctr_leg$fill, 0.25), color = "grey70", linewidth = 0.2) +
  geom_text(data = ctr_leg,
            aes(x = 0.06, y = y, label = label),
            hjust = 0, size = 2.8, fontface = "bold", color = "grey15") +

  # -- Direction boxes (right, same row height) --
  geom_rect(data = dir_leg,
            aes(xmin = DIR_XMIN, xmax = DIR_XMAX, ymin = y - BOX_HALF, ymax = y + BOX_HALF),
            fill = dir_leg$fill, color = "grey70", linewidth = 0.2) +
  geom_text(data = dir_leg,
            aes(x = (DIR_XMIN + DIR_XMAX) / 2, y = y, label = label),
            color = "white", fontface = "bold", size = 2.8) +

  # -- Section titles (same y, same font) --
  annotate("text", x = 0, y = 1.62,
           label = "Contrast Definitions:",
           hjust = 0, size = KEY_TITLE, fontface = "bold", color = "grey25") +
  annotate("text", x = DIR_XMIN, y = 1.18,
           label = "Direction:",
           hjust = 0, size = KEY_TITLE, fontface = "bold", color = "grey25") +

  scale_x_continuous(limits = c(-0.05, 6.8), expand = c(0, 0)) +
  scale_y_continuous(limits = c(-0.25, 1.80), expand = c(0, 0)) +
  theme_void() +
  theme(plot.margin = margin(0, 5.5, 0, 5.5))

cat("Keys done\n")

# 2. PANEL E — UpSet (pi-score overlap) ===================================

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
up_counts   <- numeric(n_comb)
down_counts <- numeric(n_comb)

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
  up_counts[i]   <- sum(gene_dirs == "Up")
  down_counts[i] <- sum(gene_dirs == "Down")
}

# Count mixed-direction proteins excluded from bars
n_mixed <- sum(sapply(seq_len(n_comb), function(i) {
  members <- extract_comb(cm_sub, comb_names_vec[i])
  if (length(members) == 0) return(0L)
  bits <- as.logical(as.integer(strsplit(comb_names_vec[i], "")[[1]]))
  active_contrasts <- label_to_contrast[set_names_ordered[bits]]
  gene_dirs <- sapply(members, function(g) {
    dirs <- sapply(active_contrasts, function(ctr) {
      if (g %in% names(dir_map[[ctr]])) dir_map[[ctr]][g] else NA
    })
    dirs <- dirs[!is.na(dirs)]
    if (all(dirs == "Up")) "Up" else if (all(dirs == "Down")) "Down" else "Mixed"
  })
  sum(gene_dirs == "Mixed")
}))

display_total <- up_counts + down_counts
keep_display  <- display_total > 0
comb_ord <- which(keep_display)[order(-display_total[keep_display])]
up_ord   <- up_counts[comb_ord]
down_ord <- down_counts[comb_ord]

set_display_order <- c("Aging", "Training (Young)", "Training (Old)", "Interaction")
n_int          <- length(comb_ord)
comb_names_ord <- comb_names_vec[comb_ord]
set_order_ch   <- set_name(cm_sub)
set_y_levels   <- rev(set_display_order)

bar_long <- tibble(
  x         = rep(seq_len(n_int), 2),
  direction = factor(rep(c("Up", "Down"), each = n_int), levels = c("Down", "Up")),
  count     = c(up_ord, down_ord)
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

pE_bars <- ggplot(bar_long, aes(x, count, fill = direction)) +
  geom_rect(data = bar_bg,
            aes(xmin = xmin, xmax = xmax, ymin = -Inf, ymax = Inf),
            fill = bar_bg$fill, alpha = 0.20,
            color = "grey70", linewidth = 0.2, inherit.aes = FALSE) +
  geom_col(position = position_dodge(width = 0.7), width = 0.6,
           color = "black", linewidth = 0.3) +
  geom_text(aes(label = ifelse(count > 0, count, ""), y = count / 2),
            position = position_dodge(width = 0.7), vjust = 0.5,
            size = KEY_TEXT, color = "white", fontface = "bold") +
  scale_fill_manual(values = c(Up = unname(DIR_COLORS["Up"]),
                                Down = unname(DIR_COLORS["Down"]))) +
  scale_x_continuous(expand = expansion(add = 0.6)) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.05))) +
  labs(y = "Intersection\nsize",
       title = "E  Contrast Overlap (UpSet)",
       subtitle = sprintf("DEPs by pi < 0.05; bars split by direction | Pi: %d / FDR: %d total (%d mixed-direction proteins excluded)",
                          pi_total, fdr_total, n_mixed)) +
  THEME_PUB +
  theme(axis.text.x  = element_blank(), axis.ticks.x = element_blank(),
        axis.title.x = element_blank(),
        panel.grid.major.x = element_blank(), panel.grid.minor = element_blank(),
        legend.position = "none",
        plot.margin = margin(2, 5.5, 0, 5.5))

pE_dots <- ggplot() +
  annotate("rect",
           xmin = 0.4, xmax = n_int + 0.6,
           ymin = seq_along(set_y_levels) - 0.5,
           ymax = seq_along(set_y_levels) + 0.5,
           fill = stripe_fills) +
  {if (nrow(seg_df) > 0)
    geom_segment(data = seg_df,
                 aes(x = x, xend = x, y = ymin, yend = ymax),
                 linewidth = 0.6, color = "grey25")} +
  geom_point(data = dot_df |> filter(!active),
             aes(x = x, y = ynum), color = "grey78", size = 1.2) +
  geom_point(data = dot_df |> filter(active),
             aes(x = x, y = ynum), color = "grey15", size = 1.2) +
  scale_x_continuous(expand = expansion(add = 0.6)) +
  scale_y_continuous(breaks = seq_along(set_y_levels), labels = set_y_levels,
                     expand = expansion(add = 0.15)) +
  labs(x = NULL, y = NULL) +
  THEME_PUB +
  theme(axis.text.x  = element_blank(), axis.ticks = element_blank(),
        panel.grid   = element_blank(), panel.border = element_blank(),
        axis.text.y  = element_text(size = 6.5),
        plot.margin  = margin(0, 5.5, 0, 5.5))

cat("Panel E done\n")

pE_combined <- (pE_bars / pE_dots / pKeys) +
  plot_layout(heights = c(0.46, 0.18, 0.14))

ggsave(file.path(RPT_DIR, "panel_E_upset.pdf"), pE_combined,
       width = 170, height = 180, units = "mm", device = pdf)
ggsave(file.path(RPT_DIR, "panel_E_upset.png"), pE_combined,
       width = 170, height = 180, units = "mm", dpi = 300)
