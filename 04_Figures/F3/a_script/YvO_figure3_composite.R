################################################################################
#   Figure 3 — Composite Assembly
#   Sources all per-panel scripts, builds unified key strip, assembles with
#   patchwork into a single Figure_3.pdf/png.
#
#   Panels:
#     A — Volcano ring: Aging Effect
#     B — Volcano ring: Reversal Effect (Training - Aging)
#     C — Reversal scatter (logFC Aging vs logFC Training_Old)
#     D — RRHO2 threshold-free reversal map
#     E — fGSEA NES scatter (Aging vs Training_Old)
#     F — Reversal classification: dumbbell | sankey | enrichment bar
################################################################################

# === 0. Source all panel scripts (shared setup loaded automatically) =========

source("04_Figures/F3/a_script/YvO_figure3_shared.R")
source("04_Figures/F3/a_script/YvO_panel_AB.R")   # → pA, pB
source("04_Figures/F3/a_script/YvO_panel_C.R")     # → pC, pC_combined, pC_imp_key
source("04_Figures/F3/a_script/YvO_panel_D.R")     # → pD
source("04_Figures/F3/a_script/YvO_panel_E.R")     # → pE, pE_key, pE_combined
source("04_Figures/F3/a_script/YvO_panel_F.R")     # → pF

message("All F3 panel scripts sourced — assembling composite...")

# === 1. Build unified key strip =============================================

box_w <- 0.35
box_h <- KEY_BOX_HALF
item_gap <- 0.6

key_items <- list()
x_cursor <- 0

# Direction
key_items <- c(key_items, list(tibble(
  x = x_cursor, y = 0, label = "Direction:", type = "header", fill = NA_character_)))
x_cursor <- x_cursor + 1.8
for (nm in c("Up", "Down")) {
  key_items <- c(key_items, list(tibble(
    x = x_cursor, y = 0, label = NA_character_, type = "swatch",
    fill = DIR_COLORS[nm])))
  key_items <- c(key_items, list(tibble(
    x = x_cursor + box_w + 0.1, y = 0, label = nm, type = "item",
    fill = NA_character_)))
  x_cursor <- x_cursor + box_w + 0.1 + nchar(nm) * 0.22 + item_gap
}

# Database
x_cursor <- x_cursor + 0.4
key_items <- c(key_items, list(tibble(
  x = x_cursor, y = 0, label = "Database:", type = "header", fill = NA_character_)))
x_cursor <- x_cursor + 1.8
for (nm in c("Hallmark", "GO:BP")) {
  key_items <- c(key_items, list(tibble(
    x = x_cursor, y = 0, label = NA_character_, type = "swatch",
    fill = if (nm == "Hallmark") "#AA336A" else "#00796B")))
  key_items <- c(key_items, list(tibble(
    x = x_cursor + box_w + 0.1, y = 0, label = nm, type = "item",
    fill = NA_character_)))
  x_cursor <- x_cursor + box_w + 0.1 + nchar(nm) * 0.22 + item_gap
}

# Significance (from Panel E — F3 categories)
x_cursor <- x_cursor + 0.4
key_items <- c(key_items, list(tibble(
  x = x_cursor, y = 0, label = "Significance:", type = "header", fill = NA_character_)))
x_cursor <- x_cursor + 2.2
sig_nms <- c("Reversal", "Sig Both", "Sig Aging only", "Sig Training only")
for (nm in sig_nms) {
  key_items <- c(key_items, list(tibble(
    x = x_cursor, y = 0, label = NA_character_, type = "swatch",
    fill = SIG_COLORS[nm])))
  key_items <- c(key_items, list(tibble(
    x = x_cursor + box_w + 0.1, y = 0, label = nm, type = "item",
    fill = NA_character_)))
  x_cursor <- x_cursor + box_w + 0.1 + nchar(nm) * 0.22 + item_gap
}

key_df <- bind_rows(key_items)
sw <- key_df %>% filter(type == "swatch")
hd <- key_df %>% filter(type == "header")
it <- key_df %>% filter(type == "item")

pKey <- ggplot() +
  geom_rect(data = sw,
            aes(xmin = x, xmax = x + box_w, ymin = -box_h, ymax = box_h),
            fill = sw$fill, color = KEY_BORDER, linewidth = KEY_LW) +
  geom_text(data = hd, aes(x = x, y = 0, label = label),
            hjust = 0, size = KEY_TITLE, fontface = "bold", color = KEY_HDR_COL) +
  geom_text(data = it, aes(x = x, y = 0, label = label),
            hjust = 0, size = KEY_ITEM, fontface = "bold", color = KEY_ITEM_COL) +
  coord_cartesian(ylim = c(-1, 1), expand = FALSE) +
  theme_void() +
  theme(plot.margin = margin(t = 2, r = 4, b = 2, l = 4))

# === 2. Assemble composite ==================================================

# Row 1: A | B | C (volcanos + reversal scatter)
row1 <- (pA | pB | pC) + plot_layout(widths = c(0.33, 0.33, 0.34))

# Row 2: D | E | F (RRHO + NES scatter + reversal classification)
row2 <- (pD | pE | pF) + plot_layout(widths = c(0.25, 0.35, 0.40))

fig3 <- row1 / row2 / pKey +
  plot_layout(heights = c(0.44, 0.44, 0.12)) +
  plot_annotation(
    title = "Does Resistance Training Reverse the Aging Proteome?",
    subtitle = "Protein-level and pathway-level reversal analysis: Aging vs Training in older adults.",
    theme = theme(
      plot.title    = element_text(face = "bold", size = 11, hjust = 0.5),
      plot.subtitle = element_text(size = 8, color = "grey30", hjust = 0.5)
    )
  )

# === 3. Save =================================================================

ggsave(file.path(RPT_DIR, "Figure_3.pdf"), fig3,
       width = 380, height = 350, units = "mm", device = pdf, limitsize = FALSE)
ggsave(file.path(RPT_DIR, "Figure_3.png"), fig3,
       width = 380, height = 350, units = "mm", dpi = 300, limitsize = FALSE)

cat("Saved Figure_3.pdf and Figure_3.png\n")
