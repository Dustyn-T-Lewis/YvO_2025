# F7 Panel A — Proteomic Age Gap vs Training Response
# 1x2 scatter: Mahalanobis distance vs delta_VL / delta_LBM

setwd(rprojroot::find_rstudio_root_file())
source("04_Figures/F10/a_script/style.R")

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(patchwork)
})

DAT <- "04_Figures/F10/c_data"
RPT <- "04_Figures/F10/b_reports"

age_gap  <- read_csv(file.path(DAT, "01_age_gap.csv"), show_col_types = FALSE)
cors_df  <- read_csv(file.path(DAT, "01_age_gap_cors.csv"), show_col_types = FALSE) %>%
  filter(distance == "Mahalanobis")

outcome_labs <- c(delta_VL = "Delta VL (cm)", delta_LBM = "Delta LBM (kg)")

make_scatter <- function(df, xvar, yvar, cor_row) {
  ci_label <- sprintf("r = %.2f [%.2f, %.2f]\n%s, n = %d",
                       cor_row$r, cor_row$ci_lo, cor_row$ci_hi,
                       fmt_p(cor_row$p), cor_row$n)

  ggplot(df, aes(x = .data[[xvar]], y = .data[[yvar]])) +
    geom_smooth(method = "lm", se = TRUE, color = "grey40",
                fill = "grey85", linewidth = 0.7) +
    geom_point(aes(fill = age), shape = 21, size = 2.5, stroke = 0.4) +
    scale_fill_manual(values = AGE_COLORS, name = "Age Group") +
    annotate("text", x = Inf, y = Inf, label = ci_label,
             hjust = 1.1, vjust = 1.3, size = F7_ANNOT, color = "grey30") +
    labs(x = "Proteomic Age Gap\n(Mahalanobis distance from Young centroid)",
         y = outcome_labs[yvar]) +
    FIG_THEME +
    theme(legend.position = "bottom")
}

p_vl <- make_scatter(age_gap, "mahalanobis", "delta_VL",
                      cors_df %>% filter(outcome == "delta_VL"))
p_lbm <- make_scatter(age_gap, "mahalanobis", "delta_LBM",
                       cors_df %>% filter(outcome == "delta_LBM"))

panel_A <- (p_vl | p_lbm) +
  plot_layout(guides = "collect") +
  plot_annotation(
    title = "Proteomic Age Gap Predicts Training Response",
    subtitle = "Mahalanobis distance in ME space from Young centroid",
    tag_levels = list("A")
  ) &
  theme(legend.position = "bottom",
        plot.tag = element_text(face = "bold", size = F7_TAG))

dev <- get_pdf_device()
ggsave(file.path(RPT, "panel_A_age_gap.pdf"), panel_A,
       width = 200, height = 100, units = "mm", device = dev)
ggsave(file.path(RPT, "panel_A_age_gap.png"), panel_A,
       width = 200, height = 100, units = "mm", dpi = 300, bg = "white")

message("F7 Panel A done")
