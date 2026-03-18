# F03 — Style
source("04_Figures/shared/style.R")

# --- F03-specific palettes ---
CONTRAST_COLORS <- c(
  Aging          = "#4CAF50",
  Training_Young = "#E05A4E",
  Training_Old   = "#5DA5DA",
  Interaction    = "#9B7FBF",
  Reversal       = "#FF8F00"
)

PCA_COLORS <- c(
  Young_Pre  = "#93C4DE",
  Young_Post = "#2166AC",
  Old_Pre    = "#F4A582",
  Old_Post   = "#B2182B"
)

PCA_SHAPES <- c(Young_Pre = 16, Young_Post = 17, Old_Pre = 16, Old_Post = 17)

# --- Contrast labels ---
CTR_SHORT <- c(
  Aging          = "Aging",
  Training_Young = "Tr.(Y)",
  Training_Old   = "Tr.(O)",
  Interaction    = "Tr.(O)\u2013Tr.(Y)",
  Reversal       = "Rev."
)

CTR_FACET <- c(
  Aging          = "Aging",
  Training_Young = "Tr.(Y)",
  Training_Old   = "Tr.(O)",
  Interaction    = "Tr.(O)\u2013Tr.(Y)"
)

CTR_AXIS <- c(
  Aging          = "Aging",
  Training_Young = "Tr.(Y)",
  Training_Old   = "Tr.(O)",
  Interaction    = "Tr.(O)\u2013Tr.(Y)"
)

# --- F03-specific sizing (overrides shared defaults) ---
BASE_COUNT    <- 4.0
BASE_GENE     <- 3.8
BASE_STAT     <- 4.0

FIG_TITLE_SIZE    <- 12
FIG_SUBTITLE_SIZE <- 9
FIG_STRIP_SIZE    <- 10
FIG_AXIS_TEXT     <- 9.5
FIG_LEGEND_TITLE  <- 10.5
FIG_LEGEND_TEXT   <- 9.5

# --- F03-specific theme (larger legend/axis text) ---
FIG_THEME <- theme_bw(base_size = 10) +
  theme(
    plot.title         = element_text(face = "bold", size = FIG_TITLE_SIZE),
    plot.subtitle      = element_text(face = "bold.italic", size = FIG_SUBTITLE_SIZE,
                                      color = "grey30"),
    plot.tag           = element_text(face = "bold", size = 15),
    strip.background   = element_blank(),
    strip.text         = element_text(face = "bold", size = FIG_STRIP_SIZE),
    axis.title         = element_text(face = "bold", size = 10),
    axis.text          = element_text(size = FIG_AXIS_TEXT, color = "grey15"),
    legend.title       = element_text(face = "bold", size = FIG_LEGEND_TITLE,
                                      color = "grey20"),
    legend.text        = element_text(size = FIG_LEGEND_TEXT, color = "grey15"),
    legend.key.size    = unit(3, "mm"),
    panel.grid.minor   = element_blank()
  )
