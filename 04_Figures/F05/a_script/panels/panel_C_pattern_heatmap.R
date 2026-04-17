# F05 Panel C: Aging Reversal Pattern Heatmap
# Layout: [legends] [sig|heatmap|quad] [sankey] [stacked bars in striped bg]
# 3 groups: Reversed Up, Reversed Down, Non-reversed
# Outputs: panel_C_pattern_heatmap_MAIN.pdf/.png

setwd(rprojroot::find_rstudio_root_file())
source("04_Figures/shared/style.R")
source("04_Figures/shared/go_slim_categories.R")

library(tidyverse)

RPT <- "04_Figures/F05/b_reports/panels"
DAT <- "04_Figures/F05/c_data"
dir.create(file.path(DAT, "panel_C_heatmap"), recursive = TRUE, showWarnings = FALSE)
dir.create(RPT, showWarnings = FALSE)

pdf_device <- get_pdf_device()

# =============================================================================
# 1. LOAD & CLASSIFY — 3 groups (Exacerbated merged to Non-reversed)
# =============================================================================
dep_df <- read_csv("03_DEP/c_data/03_combined_results.csv", show_col_types = FALSE)

sig_df <- dep_df %>%
  filter(!is.na(logFC_Aging), !is.na(logFC_Training_Old)) %>%
  filter(pi_score_Aging < 0.05 | pi_score_Training_Old < 0.05) %>%
  mutate(
    quadrant = case_when(
      logFC_Aging > 0 & logFC_Training_Old < 0 ~ "Reversed Up",
      logFC_Aging < 0 & logFC_Training_Old > 0 ~ "Reversed Down",
      TRUE ~ "Non-reversed"  # merged Exacerbated Up + Down
    ),
    sig_cat = case_when(
      pi_score_Aging < 0.05 & pi_score_Training_Old < 0.05 ~ "Both",
      pi_score_Aging < 0.05          ~ "Aging",
      pi_score_Training_Old < 0.05   ~ "Tr.(O)",
      TRUE ~ "NS"
    )
  )

QUAD_ORDER <- c("Reversed Up", "Reversed Down", "Non-reversed")
QUAD_COLORS <- c("Reversed Up" = "#D32F2F", "Reversed Down" = "#1976D2",
                 "Non-reversed" = "#388E3C")
QUAD_BG <- c("Reversed Up" = "#FFCDD2", "Reversed Down" = "#BBDEFB",
             "Non-reversed" = "#C8E6C9", "Tied" = "#EEEEEE")
ENDPOINT_COLORS <- c("Reversed Up" = "#8B0000", "Reversed Down" = "#0D47A1",
                     "Non-reversed" = "#1B5E20")

QUAD_LABELS <- c("Reversed Up"    = "Reversed:\nAging\u2191 Tr.(O)\u2193",
                 "Reversed Down"  = "Reversed:\nAging\u2193 Tr.(O)\u2191",
                 "Non-reversed"   = "Non-reversed")

SIG_COLORS <- c("Both" = "#2E7D32", "Aging" = "#4CAF50",
                "Tr.(O)" = "#5DA5DA", "NS" = "grey70")

# Display-label overrides for long GO Slim categories (rendering layer only).
# pathway_classification.csv export keeps the original category names.
PANEL_C_DISPLAY_LABELS <- c(
  "Carbohydrate & Energy Metabolism" = "Carb. & Energy Metab.",
  "Amino Acid & Cofactor Metabolism" = "AA & Cofactor Metab."
)

go_result <- assign_go_slim_consolidated(sig_df$gene, dep_df$gene)
sig_df <- sig_df %>%
  left_join(go_result %>% select(gene, consolidated), by = "gene") %>%
  mutate(pathway = ifelse(is.na(consolidated), "Other", as.character(consolidated)))

sig_df <- sig_df %>%
  mutate(quadrant = factor(quadrant, levels = QUAD_ORDER)) %>%
  arrange(quadrant, pathway, desc(logFC_Aging))

n_total <- nrow(sig_df)
message(sprintf("  %d significant proteins across %d quadrants", n_total,
                n_distinct(sig_df$quadrant)))

# =============================================================================
# 2. Y-COORDINATE LAYOUT
# =============================================================================
ROW_H <- 0.078

quad_counts <- sig_df %>% count(quadrant, .drop = FALSE) %>%
  mutate(quadrant = factor(quadrant, levels = QUAD_ORDER)) %>% arrange(quadrant)

y_pos <- numeric(n_total)
quad_starts <- numeric(nrow(quad_counts))
quad_ends   <- numeric(nrow(quad_counts))
idx <- 1; current_y <- 0

for (q in seq_len(nrow(quad_counts))) {
  nq <- quad_counts$n[q]
  quad_starts[q] <- current_y
  if (nq > 0) for (p in seq_len(nq)) {
    y_pos[idx] <- current_y + (p - 0.5) * ROW_H
    idx <- idx + 1
  }
  quad_ends[q] <- current_y + nq * ROW_H
  current_y <- current_y + nq * ROW_H
}
total_h <- current_y
sig_df$y <- y_pos
names(quad_starts) <- QUAD_ORDER
names(quad_ends)   <- QUAD_ORDER

# Bar region: 96% of heatmap vertical range (keys moved to bottom)
BAR_FRAC <- 0.86
BAR_YMIN <- 0                          # top of bar region = top of heatmap
BAR_YMAX <- total_h * BAR_FRAC         # bottom of bar region

# =============================================================================
# 3. PATHWAY LAYOUT
# =============================================================================
pw_counts <- sig_df %>%
  filter(pathway != "Other") %>%
  count(pathway, name = "n_prot") %>%
  arrange(desc(n_prot)) %>%
  filter(n_prot >= 2)
n_pw <- nrow(pw_counts)

row_height <- (BAR_YMAX - BAR_YMIN) / n_pw
pw_counts$y_center <- BAR_YMIN + row_height * (seq_len(n_pw) - 0.5)
pw_counts$y_top    <- BAR_YMIN + row_height * (seq_len(n_pw) - 1)
pw_counts$y_bot    <- BAR_YMIN + row_height * seq_len(n_pw)
BAR_H <- row_height * 0.78   # taller bars improve count-integer legibility

# Dominant quadrant: by count, break ties with sum(|logFC|)
dom_quad <- sig_df %>%
  filter(pathway %in% pw_counts$pathway) %>%
  group_by(pathway, quadrant) %>%
  summarise(n = n(),
            lfc_sum = sum(abs(logFC_Aging) + abs(logFC_Training_Old)),
            .groups = "drop") %>%
  group_by(pathway) %>%
  mutate(is_max = n == max(n), n_tied = sum(is_max)) %>%
  arrange(pathway, desc(n), desc(lfc_sum)) %>%
  slice_head(n = 1) %>%
  ungroup() %>%
  mutate(dom_quad = as.character(quadrant)) %>%
  select(pathway, dom_quad)

pw_counts <- pw_counts %>% left_join(dom_quad, by = "pathway")

# =============================================================================
# 4. X-COORDINATE LAYOUT — x-axis to ~30
# =============================================================================
STRIP_W <- 0.10; TILE_W <- 0.70

X_SIG   <- 0.8
X_COL1  <- X_SIG + STRIP_W/2 + TILE_W/2 + 0.01
X_COL2  <- X_COL1 + TILE_W + 0.01
X_QUAD  <- X_COL2 + TILE_W/2 + STRIP_W/2 + 0.01
HEAT_LEFT  <- X_SIG - STRIP_W/2
HEAT_RIGHT <- X_QUAD + STRIP_W/2

X_SANK_L  <- HEAT_RIGHT + 0.08
X_SANK_R  <- 3.2
X_BAR_L   <- 3.3
BAR_SCALE <- 0.20

count_max <- max(pw_counts$n_prot)
X_BAR_MAX <- max(X_BAR_L + 30 * BAR_SCALE, X_BAR_L + count_max * BAR_SCALE)

PW_OUT <- 250; PH_OUT <- 380

# =============================================================================
# 5. STACKED BAR DATA
# =============================================================================
bar_data <- sig_df %>%
  filter(pathway %in% pw_counts$pathway) %>%
  count(pathway, quadrant, name = "n_seg") %>%
  left_join(pw_counts %>% select(pathway, y_center, n_prot), by = "pathway") %>%
  group_by(pathway) %>%
  arrange(pathway, desc(n_seg)) %>%
  mutate(
    cum_n = cumsum(n_seg) - n_seg,
    xmin = X_BAR_L + cum_n * BAR_SCALE,
    xmax = X_BAR_L + (cum_n + n_seg) * BAR_SCALE,
    ymin = y_center - BAR_H / 2,
    ymax = y_center + BAR_H / 2
  ) %>% ungroup()

bg_stripes <- pw_counts %>%
  transmute(
    xmin = X_BAR_L - 0.05, xmax = X_BAR_MAX + 0.05,
    ymin = y_top, ymax = y_bot,
    fill = QUAD_BG[dom_quad]
  )

pw_labels <- pw_counts %>%
  mutate(display_pathway = ifelse(pathway %in% names(PANEL_C_DISPLAY_LABELS),
                                   PANEL_C_DISPLAY_LABELS[pathway], pathway)) %>%
  transmute(x = X_BAR_L + n_prot * BAR_SCALE + 0.08, y = y_center,
            label = ifelse(nchar(display_pathway) > 18,
                           stringr::str_wrap(display_pathway, width = 18),
                           display_pathway))

count_ticks <- tibble(
  val = pretty(c(0, count_max), n = 4),
  x = X_BAR_L + val * BAR_SCALE,
  y_tick_top = BAR_YMAX, y_tick_bot = BAR_YMAX + ROW_H * 0.7,
  y_label = BAR_YMAX + ROW_H * 1.6
) %>% filter(val >= 0, val <= count_max)

# =============================================================================
# 6. SANKEY — fractional y-bands
# =============================================================================
flow_df <- sig_df %>%
  filter(pathway %in% pw_counts$pathway) %>%
  count(quadrant, pathway, name = "n_flow") %>%
  filter(n_flow > 0)

source_bands <- flow_df %>%
  group_by(quadrant) %>%
  mutate(
    total_q = sum(n_flow), frac = n_flow / total_q,
    q_start = quad_starts[as.character(quadrant)],
    q_end   = quad_ends[as.character(quadrant)],
    q_height = q_end - q_start
  ) %>%
  arrange(quadrant, match(pathway, pw_counts$pathway)) %>%
  mutate(
    cum_frac = cumsum(frac) - frac,
    src_top = q_start + cum_frac * q_height,
    src_bot = q_start + (cum_frac + frac) * q_height
  ) %>% ungroup()

target_bands <- bar_data %>%
  group_by(pathway) %>%
  arrange(pathway, desc(n_seg)) %>%
  mutate(
    frac = n_seg / sum(n_seg),
    cum_frac = cumsum(frac) - frac,
    tgt_top = ymin + cum_frac * (ymax - ymin),
    tgt_bot = ymin + (cum_frac + frac) * (ymax - ymin)
  ) %>% ungroup() %>%
  select(pathway, quadrant, tgt_top, tgt_bot)

ribbon_df <- source_bands %>%
  select(quadrant, pathway, n_flow, src_top, src_bot) %>%
  left_join(target_bands, by = c("quadrant", "pathway"))

all_ribbons <- pmap_dfr(ribbon_df, function(quadrant, pathway, n_flow,
                                              src_top, src_bot, tgt_top, tgt_bot) {
  rid <- paste(quadrant, pathway, sep = "___")
  df <- make_sigmoid_ribbon(X_SANK_L, X_SANK_R, src_top, src_bot, tgt_top, tgt_bot,
                            n_pts = 60, ribbon_id = rid)
  df$quadrant <- quadrant; df$pathway <- pathway; df
})

endpoint_bars <- bar_data %>%
  transmute(xmin = X_SANK_R - 0.04, xmax = X_SANK_R + 0.04,
            ymin, ymax, quadrant = as.character(quadrant))

# =============================================================================
# 7. HEATMAP
# =============================================================================
fc_max <- max(abs(c(sig_df$logFC_Aging, sig_df$logFC_Training_Old)), na.rm = TRUE)

lfc_to_color <- function(v, fc_max) {
  v <- pmax(-fc_max, pmin(fc_max, v))
  ifelse(v >= 0,
         scales::seq_gradient_pal("#FFFFFF", "#B2182B")(v / fc_max),
         scales::seq_gradient_pal("#2166AC", "#FFFFFF")((v + fc_max) / fc_max))
}

heat_tiles <- bind_rows(
  sig_df %>% transmute(x = X_COL1, y, w = TILE_W, h = ROW_H,
                        fill = lfc_to_color(logFC_Aging, fc_max)),
  sig_df %>% transmute(x = X_COL2, y, w = TILE_W, h = ROW_H,
                        fill = lfc_to_color(logFC_Training_Old, fc_max))
)
sig_tiles <- sig_df %>%
  transmute(x = X_SIG, y, w = STRIP_W, h = ROW_H, fill = SIG_COLORS[sig_cat])
quad_tiles <- sig_df %>%
  transmute(x = X_QUAD, y, w = STRIP_W, h = ROW_H,
            fill = QUAD_COLORS[as.character(quadrant)])

divider_ys <- quad_ends[1:(length(QUAD_ORDER) - 1)]
divider_ys <- divider_ys[divider_ys > 0 & divider_ys < total_h]

col_headers <- tibble(x = c(X_COL1, X_COL2), y = total_h + ROW_H * 0.7,
                      label = c("Aging", "Tr.(O)"))

# =============================================================================
# 8. LEGENDS
# =============================================================================
n_g <- 50
HEAT_MID  <- (HEAT_LEFT + HEAT_RIGHT) / 2
GRAD_HALFW <- (HEAT_RIGHT - HEAT_LEFT) * 0.30   # narrower, centered under heatmap
GRAD_L    <- HEAT_MID - GRAD_HALFW
GRAD_R    <- HEAT_MID + GRAD_HALFW
grad_xs <- seq(GRAD_L, GRAD_R, length.out = n_g)
grad_h_legend <- tibble(
  xmin = grad_xs,
  xmax = lead(grad_xs, default = max(grad_xs) + diff(grad_xs)[1]),
  fv = seq(-fc_max, fc_max, length.out = n_g),
  fill = lfc_to_color(seq(-fc_max, fc_max, length.out = n_g), fc_max)
)
GRAD_Y <- total_h + ROW_H * 4.44

FONT_UNI <- 4.6  # unified size for count axis, gradient labels, key headers
FONT_BAR <- 3.2  # bar segment numbers — bumped from 2.8 for in-bar legibility
FONT_PW  <- FONT_UNI * 0.75 # pathway bar labels — reduced further for composite fit

sig_leg_cats <- c("Both", "Aging", "Tr.(O)")

# --- Keys at BOTTOM: sig centered under heatmap, quad centered under bars ---
KEY_SQ_SZ   <- 5.0
KEY_TXT_SZ  <- 3.8
KEY_HDR_SZ  <- 4.8
KEY_Y_BASE  <- BAR_YMAX + ROW_H * 15.5               # keys in whitespace below bar plot, above heatmap bottom
KEY_DY      <- ROW_H * 3.8                        # compensates for F05's taller y-data extent; ~matches F04's ~7.7 mm mm-spacing

# Both keys as vertical stacks, centered together under bar plot (sig on left, quad on right)
bar_mid     <- (X_BAR_L + X_BAR_MAX) / 2
KEY_X_SIG   <- bar_mid - 1.4
KEY_X_QUAD  <- bar_mid + 0.6

sig_cats_f05 <- c("Tr.(O)", "Aging", "Both")
sig_key_df <- tibble(
  x     = KEY_X_SIG,
  y     = KEY_Y_BASE + KEY_DY * seq_along(sig_cats_f05),
  label = c("Sig Training only", "Sig Aging only", "Sig Both"),
  fill  = SIG_COLORS[sig_cats_f05]
)

quad_key_df <- tibble(
  x     = KEY_X_QUAD,
  y     = KEY_Y_BASE + KEY_DY * seq_along(QUAD_ORDER),
  label = QUAD_ORDER,
  fill  = QUAD_COLORS[QUAD_ORDER]
)

# =============================================================================
# 9. RENDER
# =============================================================================
p <- ggplot() +
  geom_rect(data = bg_stripes,
            aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
            fill = bg_stripes$fill, color = "grey70", linewidth = 0.2) +
  geom_rect(data = heat_tiles,
            aes(xmin = x - w/2, xmax = x + w/2, ymin = y - h/2, ymax = y + h/2),
            fill = heat_tiles$fill, color = NA) +
  geom_rect(data = sig_tiles,
            aes(xmin = x - w/2, xmax = x + w/2, ymin = y - h/2, ymax = y + h/2),
            fill = sig_tiles$fill, color = NA) +
  geom_rect(data = quad_tiles,
            aes(xmin = x - w/2, xmax = x + w/2, ymin = y - h/2, ymax = y + h/2),
            fill = quad_tiles$fill, color = NA) +
  geom_segment(data = tibble(y = divider_ys),
               aes(x = X_SIG - STRIP_W/2, xend = X_QUAD + STRIP_W/2,
                   y = y, yend = y),
               color = "grey30", linewidth = 0.4) +
  geom_text(data = col_headers, aes(x = x, y = y, label = label),
            size = FONT_UNI - 1.0, fontface = "bold", color = "grey20") +
  geom_polygon(data = all_ribbons, aes(x = x, y = y, group = ribbon_id),
               fill = QUAD_COLORS[all_ribbons$quadrant], alpha = 0.40, color = NA) +
  geom_rect(data = endpoint_bars,
            aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
            fill = ENDPOINT_COLORS[endpoint_bars$quadrant], color = NA) +
  geom_rect(data = bar_data,
            aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
            fill = QUAD_COLORS[as.character(bar_data$quadrant)],
            color = "black", linewidth = 0.3) +
  geom_text(data = bar_data,
            aes(x = (xmin + xmax) / 2, y = (ymin + ymax) / 2, label = n_seg),
            size = FONT_BAR, fontface = "bold", color = "white") +
  geom_text(data = pw_labels, aes(x = x, y = y, label = label),
            size = FONT_PW, hjust = 0, fontface = "bold", color = "grey15",
            lineheight = 0.8) +
  annotate("segment", x = X_BAR_L, xend = X_BAR_MAX,
           y = BAR_YMAX, yend = BAR_YMAX, color = "grey20", linewidth = 0.5) +
  geom_segment(data = count_ticks,
               aes(x = x, xend = x, y = y_tick_top, yend = y_tick_bot),
               color = "grey20", linewidth = 0.3) +
  geom_text(data = count_ticks, aes(x = x, y = y_label, label = val),
            size = FONT_UNI, fontface = "bold", color = "grey20") +
  annotate("text", x = (X_BAR_L + X_BAR_MAX) / 2, y = BAR_YMAX + ROW_H * 4.5,
           label = "Protein count", size = FONT_UNI - 1.0, fontface = "bold", color = "grey20") +
  geom_rect(data = grad_h_legend,
            aes(xmin = xmin, xmax = xmax, ymin = GRAD_Y, ymax = GRAD_Y + ROW_H * 2.05),
            fill = grad_h_legend$fill, color = NA) +
  annotate("text", x = GRAD_L - 0.10, y = GRAD_Y + ROW_H * 3.4,
           label = sprintf("%.1f", -fc_max), size = 3.1, fontface = "bold",
           color = "grey20", hjust = 0) +
  annotate("text", x = GRAD_R + 0.10, y = GRAD_Y + ROW_H * 3.4,
           label = sprintf("+%.1f", fc_max), size = 3.1, fontface = "bold",
           color = "grey20", hjust = 1) +
  annotate("text", x = HEAT_MID, y = GRAD_Y + ROW_H * 3.4,
           label = "0", size = 3.1, fontface = "bold", color = "grey30") +
  annotate("text", x = GRAD_L - 0.08, y = GRAD_Y + ROW_H * 1.0,
           label = "logFC", size = 3.6, fontface = "bold", color = "grey15",
           hjust = 1, vjust = 0.5) +
  # --- Bottom keys: Significance (header centered above single row of swatches) ---
  annotate("text", x = KEY_X_SIG - 0.26, y = KEY_Y_BASE,
           label = "Significance", size = KEY_HDR_SZ, fontface = "bold",
           color = "grey25", hjust = 0) +
  geom_point(data = sig_key_df, aes(x = x - 0.26, y = y),
             shape = 22, size = KEY_SQ_SZ,
             fill = sig_key_df$fill, color = "grey30", stroke = 0.4) +
  geom_text(data = sig_key_df, aes(x = x - 0.02, y = y, label = label),
            size = KEY_TXT_SZ, color = "grey20", hjust = 0) +
  # --- Bottom keys: Quadrant (header centered above column of swatches) ---
  annotate("text", x = KEY_X_QUAD - 0.26, y = KEY_Y_BASE,
           label = "Quadrant", size = KEY_HDR_SZ, fontface = "bold",
           color = "grey25", hjust = 0) +
  geom_point(data = quad_key_df, aes(x = x - 0.26, y = y),
             shape = 22, size = KEY_SQ_SZ,
             fill = quad_key_df$fill, color = "grey30", stroke = 0.4) +
  geom_text(data = quad_key_df, aes(x = x - 0.02, y = y, label = label),
            size = KEY_TXT_SZ, color = "grey20", hjust = 0) +
  scale_y_reverse() +
  coord_cartesian(xlim = c(0.0, X_BAR_MAX + 2.0),
                  ylim = c(total_h + ROW_H * 11.0, -ROW_H * 0.05),
                  expand = FALSE) +
  labs(title = "Aging Reversal Patterns",
       subtitle = sprintf("%d proteins | GO Slim | %d pathways", n_total, n_pw)) +
  theme_void() +
  theme(plot.margin = margin(3, 3, 3, 3, "mm"),  # wide composite — harmonized with F04/F05 polish
        plot.title = element_text(face = "bold", size = 14, hjust = 0,
                                  margin = margin(l = 31.5, unit = "mm")),
        plot.subtitle = element_text(face = "italic", size = 10,
                                     hjust = 0, color = "grey40",
                                     margin = margin(l = 31.5, unit = "mm")),
        plot.title.position = "panel")

# =============================================================================
# 10. SAVE
# =============================================================================
ggsave(file.path(RPT, "MAIN_panel_C_pattern_heatmap.png"), p,
       width = PW_OUT, height = PH_OUT, units = "mm", dpi = 300)

# =============================================================================
# 11. DATA EXPORTS
# =============================================================================
sig_df %>%
  transmute(gene, quadrant = as.character(quadrant), sig_cat, pathway,
            logFC_Aging = round(logFC_Aging, 4),
            logFC_Training_Old = round(logFC_Training_Old, 4)) %>%
  write_csv(file.path(DAT, "panel_C_heatmap", "pattern_classification.csv"))
flow_df %>% write_csv(file.path(DAT, "panel_C_heatmap", "sankey_links.csv"))
bar_data %>%
  select(pathway, quadrant, n_seg, xmin, xmax) %>%
  write_csv(file.path(DAT, "panel_C_heatmap", "bar_data.csv"))

message("F05 Panel C (pattern heatmap) done")
