# F7 Panel B — Within-Group ME Prediction
# 2x2 facet grid: rows = VL/LBM, cols = Young/Old

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

me_pre   <- readRDS(file.path(DAT, "me_pre.rds"))
pheno    <- readRDS(file.path(DAT, "pheno.rds"))
age_df   <- readRDS(file.path(DAT, "age.rds"))
top_me   <- read_csv(file.path(DAT, "02_top_me.csv"), show_col_types = FALSE)

outcome_labs <- c(delta_VL = "Delta VL (cm)", delta_LBM = "Delta LBM (kg)")

plots <- list()
for (i in seq_len(nrow(top_me))) {
  row <- top_me[i, ]
  grp <- row$group
  oc  <- row$outcome
  me_name <- row$me

  idx <- which(age_df$age == grp)
  plot_df <- tibble(
    x = me_pre[idx, me_name],
    y = pheno[[oc]][idx]
  )

  me_label <- gsub("^ME", "", me_name)
  ci_label <- sprintf("r = %.2f [%.2f, %.2f]\n%s, n = %d",
                       row$r, row$ci_lo, row$ci_hi, fmt_p(row$p), row$n)

  fill_col <- AGE_COLORS[grp]

  plots[[i]] <- ggplot(plot_df, aes(x = x, y = y)) +
    geom_smooth(method = "lm", se = TRUE, color = "grey40",
                fill = "grey85", linewidth = 0.7) +
    geom_point(fill = fill_col, shape = 21, size = 2.5, stroke = 0.4) +
    annotate("text", x = Inf, y = Inf, label = ci_label,
             hjust = 1.1, vjust = 1.3, size = F7_ANNOT, color = "grey30") +
    labs(title = sprintf("%s - %s", grp, outcome_labs[oc]),
         x = sprintf("Baseline %s eigengene", me_label),
         y = outcome_labs[oc]) +
    FIG_THEME +
    theme(plot.title = element_text(size = 10, face = "bold",
                                     color = fill_col))
}

panel_B <- wrap_plots(plots, ncol = 2) +
  plot_annotation(
    title = "Within-Group Module Eigengene Prediction",
    subtitle = "Best ME per group x outcome (uncorrected, exploratory)",
    tag_levels = list("B")
  ) &
  theme(plot.tag = element_text(face = "bold", size = F7_TAG))

dev <- get_pdf_device()
ggsave(file.path(RPT, "panel_B_within_me.pdf"), panel_B,
       width = 200, height = 180, units = "mm", device = dev)
ggsave(file.path(RPT, "panel_B_within_me.png"), panel_B,
       width = 200, height = 180, units = "mm", dpi = 300, bg = "white")

message("F7 Panel B done")
