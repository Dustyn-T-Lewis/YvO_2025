# F7 Panel D — Method Comparison Summary
# Horizontal dot plot of observed r across all piloted methods

setwd(rprojroot::find_rstudio_root_file())
source("04_Figures/F10/a_script/style.R")

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(ggplot2)
})

DAT <- "04_Figures/F10/c_data"
RPT <- "04_Figures/F10/b_reports"

methods_df <- read_csv(file.path(DAT, "04_method_comparison.csv"), show_col_types = FALSE)

# Add significance flag
methods_df <- methods_df %>%
  mutate(
    sig = ifelse(p < 0.05, "p < 0.05", "n.s."),
    method = factor(method, levels = rev(method))
  )

sig_colors <- c("p < 0.05" = "#2166AC", "n.s." = "grey60")

panel_D <- ggplot(methods_df, aes(x = r, y = method)) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey50") +
  geom_errorbar(aes(xmin = ci_lo, xmax = ci_hi),
                 width = 0.3, color = "grey50", linewidth = 0.4,
                 orientation = "y") +
  geom_point(aes(fill = sig), shape = 21, size = 3, stroke = 0.5) +
  scale_fill_manual(values = sig_colors, name = "Significance") +
  labs(title = "Baseline Prediction: Method Comparison",
       subtitle = "Observed Pearson r with 95% CI (all nominal, uncorrected)",
       x = "Pearson r (baseline predictor vs training response)",
       y = NULL) +
  FIG_THEME +
  theme(legend.position = "bottom",
        plot.tag = element_text(face = "bold", size = F7_TAG))

dev <- get_pdf_device()
ggsave(file.path(RPT, "panel_D_methods.pdf"), panel_D,
       width = 180, height = 120, units = "mm", device = dev)
ggsave(file.path(RPT, "panel_D_methods.png"), panel_D,
       width = 180, height = 120, units = "mm", dpi = 300, bg = "white")

message("F7 Panel D done")
