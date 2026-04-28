# Print-scale override for 380mm-wide composites
# J Physiol double-column = 178mm. At print, 380mm figures scale by 178/380 ≈ 0.47x.
# All text must be ~2.13x larger at source to print at spec (10-12pt tags, 7pt body).
# Source this AFTER style.R in any panel script destined for a 380mm composite.

PRINT_SCALE <- 380 / 178   # ≈ 2.13

# Override global text hierarchy (from style.R)
FIG_TITLE_SIZE    <- round(FIG_TITLE_SIZE * PRINT_SCALE, 1)     # 7 → 14.9
FIG_SUBTITLE_SIZE <- round(FIG_SUBTITLE_SIZE * PRINT_SCALE, 1)  # 4 → 8.5
FIG_STRIP_SIZE    <- round(FIG_STRIP_SIZE * PRINT_SCALE, 1)     # 5 → 10.7
FIG_AXIS_TEXT     <- round(FIG_AXIS_TEXT * PRINT_SCALE, 1)       # 5 → 10.7
FIG_LEGEND_TITLE  <- round(FIG_LEGEND_TITLE * PRINT_SCALE, 1)   # 5 → 10.7
FIG_LEGEND_TEXT   <- round(FIG_LEGEND_TEXT * PRINT_SCALE, 1)     # 4 → 8.5

# Override annotation base sizes (used by scale_text())
# In-plot annotation sizes: scale then drop 2pt to avoid overfilling panels
BASE_PATHWAY  <- round(BASE_PATHWAY * PRINT_SCALE, 1) - 2   # 6.0 → 4.0
BASE_GENE     <- round(BASE_GENE * PRINT_SCALE, 1) - 2      # 5.3 → 3.3
BASE_STAT     <- round(BASE_STAT * PRINT_SCALE, 1) - 2      # 5.3 → 3.3
BASE_QUADRANT <- round(BASE_QUADRANT * PRINT_SCALE, 1) - 2  # 6.0 → 4.0
BASE_COUNT    <- round(BASE_COUNT * PRINT_SCALE, 1) - 2     # 5.3 → 3.3

# Rebuild FIG_THEME with scaled sizes
FIG_THEME <- theme_bw(base_size = round(6 * PRINT_SCALE, 1), base_family = "Helvetica") +
  theme(
    plot.title         = element_text(face = "bold", size = FIG_TITLE_SIZE,
                                      margin = margin(b = 1)),
    plot.subtitle      = element_text(face = "bold.italic", size = FIG_SUBTITLE_SIZE,
                                      color = "grey30", margin = margin(t = 0, b = 2)),
    plot.tag           = element_text(face = "bold", size = FIG_TITLE_SIZE),
    strip.background   = element_blank(),
    strip.text         = element_text(face = "bold", size = FIG_STRIP_SIZE),
    axis.title.x       = element_text(face = "bold", size = FIG_AXIS_TEXT,
                                      margin = margin(t = 0)),
    axis.title.y       = element_text(face = "bold", size = FIG_AXIS_TEXT,
                                      margin = margin(r = -1)),
    axis.text          = element_text(size = FIG_AXIS_TEXT, color = "grey15"),
    legend.title       = element_text(face = "bold", size = FIG_LEGEND_TITLE,
                                      color = "grey20"),
    legend.text        = element_text(size = FIG_LEGEND_TEXT, color = "grey15"),
    legend.key.size    = unit(2.5 * PRINT_SCALE, "mm"),
    panel.grid.minor   = element_blank()
  )
