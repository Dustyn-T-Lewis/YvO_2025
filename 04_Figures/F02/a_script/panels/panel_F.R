# F02 — Panel F: DEP Rank Location (Barcode Plot)
# Adapted from 04_Figures/F03/a_script/panel_B.R
# Shows where Pi < 0.05 DEPs sit in the t-statistic-ranked proteome
# Dual density traces (Up/Down) + direction-colored barcode ticks
# Faceted rows (one per contrast) for composite alignment
# Outputs: pF (ggplot object), panel_F_barcode_MAIN.pdf/.png

setwd(rprojroot::find_rstudio_root_file())
source("04_Figures/F02/a_script/style.R")

library(dplyr)
library(readr)
library(tidyr)
library(ggplot2)


DEP_FILE <- "03_DEP/c_data/03_combined_results.csv"
RPT      <- "04_Figures/F02/b_reports/panels"
DAT      <- "04_Figures/F02/c_data"
dir.create(RPT, recursive = TRUE, showWarnings = FALSE)
dir.create(DAT, recursive = TRUE, showWarnings = FALSE)

CONTRASTS <- c("Aging", "Training_Young", "Training_Old", "Interaction")
dep_df    <- read_csv(DEP_FILE, show_col_types = FALSE)
pdf_device <- get_pdf_device()

PD_W <- 165
PD_H <- 120

# --- Build long-form data: rank position + DEP status per contrast ---
rank_list <- lapply(CONTRASTS, function(ctr) {
  t_col   <- paste0("t_", ctr)
  pi_col  <- paste0("pi_score_", ctr)
  lfc_col <- paste0("logFC_", ctr)

  dep_df |>
    filter(!is.na(.data[[t_col]])) |>
    arrange(.data[[t_col]]) |>
    mutate(
      rank_frac = seq_len(n()) / n(),
      is_dep    = !is.na(.data[[pi_col]]) & .data[[pi_col]] < 0.05,
      direction = case_when(
        !is_dep              ~ NA_character_,
        .data[[lfc_col]] > 0 ~ "Up",
        TRUE                 ~ "Down"
      ),
      contrast = ctr
    ) |>
    select(gene, contrast, rank_frac, is_dep, direction)
})
rank_df <- bind_rows(rank_list)
rank_df$contrast <- factor(rank_df$contrast, levels = CONTRASTS)

# DEP subset for ticks and density
dep_only <- rank_df |> filter(is_dep)
dep_only$direction <- factor(dep_only$direction, levels = c("Up", "Down"))

# DEP counts per contrast for annotation
dep_counts <- dep_only |>
  group_by(contrast) |>
  summarise(
    n_up    = sum(direction == "Up"),
    n_down  = sum(direction == "Down"),
    n_total = n(),
    .groups = "drop"
  ) |>
  mutate(label = sprintf("n = %d  (%d\u2191  %d\u2193)", n_total, n_up, n_down))

write.csv(dep_counts, file.path(DAT, "panel_F_barcode_enrichment.csv"), row.names = FALSE)

# --- Pre-compute density curves so we can normalize and control y-range ---
DENS_PAD <- 0.06
dens_list <- lapply(split(dep_only, dep_only$contrast), function(ctr_df) {
  lapply(split(ctr_df, ctr_df$direction, drop = TRUE), function(dir_df) {
    if (nrow(dir_df) < 2) return(NULL)
    d <- density(dir_df$rank_frac, adjust = 1.8,
                 from = -DENS_PAD, to = 1 + DENS_PAD, n = 512)
    tibble(x = d$x, y = d$y, direction = dir_df$direction[1],
           contrast = dir_df$contrast[1])
  }) |> bind_rows()
}) |> bind_rows()

# Normalize density per contrast: peak = 1.0 (makes panels comparable)
dens_list <- dens_list |>
  group_by(contrast) |>
  mutate(y_norm = y / max(y)) |>
  ungroup()
dens_list$direction <- factor(dens_list$direction, levels = c("Up", "Down"))
dens_list$contrast  <- factor(dens_list$contrast, levels = CONTRASTS)

TICK_DEPTH <- -0.25
ANNOT_SZ <- scale_text(BASE_STAT - 0.3, PD_W)
LABEL_NUDGE <- 0.06

# --- Compute peak positions for label placement ---
peak_pos <- dens_list |>
  group_by(contrast, direction) |>
  slice_max(y_norm, n = 1, with_ties = FALSE) |>
  ungroup() |>
  select(contrast, direction, peak_x = x, peak_y = y_norm)

# DEP counts per contrast x direction
n_down <- dep_only |> filter(direction == "Down") |> count(contrast) |> tibble::deframe()
n_up   <- dep_only |> filter(direction == "Up")   |> count(contrast) |> tibble::deframe()
DESC_DOWN <- c(Aging = "proteins lower in older vs young",
               Training_Young = "proteins dec. with training",
               Training_Old = "proteins dec. with training",
               Interaction = "proteins with greater Young response")
DESC_UP <- c(Aging = "proteins higher in older vs young",
             Training_Young = "proteins inc. with training",
             Training_Old = "proteins inc. with training",
             Interaction = "proteins with greater Old response")

# --- Build combined annotation data for faceted plot ---
bg_wash <- tibble(
  contrast = factor(CONTRASTS, levels = CONTRASTS),
  fill     = unname(CONTRAST_COLORS[CONTRASTS]),
  xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf
)

# Down direction annotations
ad_all <- peak_pos |>
  filter(direction == "Down") |>
  mutate(
    ctr      = as.character(contrast),
    label_x  = peak_x + LABEL_NUDGE,
    label_y  = peak_y * 0.85,
    label    = paste(n_down[ctr], DESC_DOWN[ctr])
  ) |>
  filter(!is.na(label))

# Up direction annotations
au_all <- peak_pos |>
  filter(direction == "Up") |>
  mutate(
    ctr      = as.character(contrast),
    label_x  = peak_x - LABEL_NUDGE,
    label_y  = if_else(ctr == "Interaction", 0.30, peak_y * 0.85),
    label    = paste(n_up[ctr], DESC_UP[ctr])
  ) |>
  filter(!is.na(label))

# Connector segments
cd_all <- ad_all |>
  mutate(x_start = peak_x, y_start = peak_y, x_end = label_x, y_end = label_y)
cu_all <- au_all |>
  mutate(x_start = peak_x, y_start = peak_y, x_end = label_x, y_end = label_y)

# --- Single faceted barcode plot ---
pF <- ggplot() +
  # Background contrast wash (per facet) — darkened to match C/D/E's
  # canonical 0.20 alpha (0.18 is a slight pull-back to avoid over-darkening)
  geom_rect(data = bg_wash,
            aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
            fill = bg_wash$fill, alpha = 0.18, inherit.aes = FALSE) +
  # Density ribbons
  geom_ribbon(data = dens_list,
              aes(x = x, ymin = 0, ymax = y_norm, fill = direction),
              alpha = 0.30, outline.type = "upper") +
  geom_line(data = dens_list,
            aes(x = x, y = y_norm, color = direction),
            linewidth = 0.5) +
  # Barcode ticks
  geom_segment(data = dep_only,
               aes(x = rank_frac, xend = rank_frac,
                   y = 0, yend = TICK_DEPTH, color = direction),
               linewidth = 0.35, alpha = 0.8) +
  # Zero line
  geom_hline(yintercept = 0, linewidth = 0.25, color = "grey50") +
  # Down connector segments
  {if (nrow(cd_all) > 0)
    geom_segment(data = cd_all,
                 aes(x = x_start, xend = x_end, y = y_start, yend = y_end),
                 linewidth = 0.3, color = unname(DIR_COLORS["Down"]),
                 alpha = 0.4, inherit.aes = FALSE)} +
  # Down labels
  {if (nrow(ad_all) > 0)
    geom_text(data = ad_all,
              aes(x = label_x, y = label_y, label = label),
              hjust = 0, vjust = 0.5, size = ANNOT_SZ,
              color = unname(DIR_COLORS["Down"]), fontface = "bold",
              inherit.aes = FALSE)} +
  # Up connector segments
  {if (nrow(cu_all) > 0)
    geom_segment(data = cu_all,
                 aes(x = x_start, xend = x_end, y = y_start, yend = y_end),
                 linewidth = 0.3, color = unname(DIR_COLORS["Up"]),
                 alpha = 0.4, inherit.aes = FALSE)} +
  # Up labels
  {if (nrow(au_all) > 0)
    geom_text(data = au_all,
              aes(x = label_x, y = label_y, label = label),
              hjust = 1, vjust = 0.5, size = ANNOT_SZ,
              color = unname(DIR_COLORS["Up"]), fontface = "bold",
              inherit.aes = FALSE)} +
  # Scales
  scale_fill_manual(values = c(Up = unname(DIR_COLORS["Up"]),
                                Down = unname(DIR_COLORS["Down"]))) +
  scale_color_manual(values = c(Up = unname(DIR_COLORS["Up"]),
                                 Down = unname(DIR_COLORS["Down"]))) +
  scale_x_continuous(labels = scales::percent_format(accuracy = 1)) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.02))) +
  coord_cartesian(xlim = c(-DENS_PAD, 1 + DENS_PAD),
                  ylim = c(TICK_DEPTH, 1.02)) +
  facet_grid(contrast ~ ., switch = "y",
             labeller = labeller(contrast = CTR_FACET)) +
  labs(title = "DEP Rank Location",
       subtitle = sprintf("%s proteins | %d \u03A0 DEPs (t-ranked)",
                          format(length(unique(rank_df$gene)), big.mark = ","),
                          sum(dep_counts$n_total)),
       x = "Rank position (by t-statistic)", y = NULL,
       tag = "F") +
  FIG_THEME +
  theme(
    plot.subtitle      = element_text(size = FIG_SUBTITLE_SIZE,
                                      face = "bold.italic", color = "grey40"),
    legend.position    = "none",
    axis.text.y        = element_blank(),
    axis.ticks.y       = element_blank(),
    strip.text.y.left  = element_text(face = "bold", size = FIG_AXIS_TEXT - 0.5,
                                       angle = 0, hjust = 1),
    strip.background   = element_blank(),
    strip.placement    = "outside",
    panel.grid.major.y = element_blank(),
    panel.grid.minor   = element_blank(),
    panel.spacing.y    = unit(2, "pt")
  )

ggsave(file.path(RPT, "MAIN_panel_F_barcode.png"), pF,
       width = PD_W, height = PD_H, units = "mm", dpi = 300)
