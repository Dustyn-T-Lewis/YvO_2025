# F02 — Panel B: logFC Density Histograms (Effect Size Distributions)
# Outputs: pB (ggplot object), MAIN_panel_B_logfc_density.{pdf,png}

setwd(rprojroot::find_rstudio_root_file())
source("04_Figures/F02/a_script/style.R")

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(readr)
  library(ggplot2)
})

PD_W <- 60; PD_H <- 120

RPT_PNG <- "04_Figures/F02/b_reports/main/png/panels"
RPT_PDF <- "04_Figures/F02/b_reports/main/pdf/panels"
DAT_DIR <- "04_Figures/F02/c_data"
dir.create(RPT_PNG, recursive = TRUE, showWarnings = FALSE)
dir.create(RPT_PDF, recursive = TRUE, showWarnings = FALSE)
dir.create(DAT_DIR, recursive = TRUE, showWarnings = FALSE)

dep_df <- read_csv("03_DEP/c_data/03_combined_results.csv",
                   show_col_types = FALSE)

pdf_device <- get_pdf_device()

# Long-format logFC for the visualization (3 displayed contrasts) AND for
# the supplement (all contrasts, including uniprot_id for traceability).
id_cols <- intersect(c("uniprot_id", "gene"), names(dep_df))
lfc_long_all <- dep_df |>
  dplyr::select(dplyr::all_of(id_cols), starts_with("logFC_")) |>
  pivot_longer(starts_with("logFC_"), names_to = "contrast", values_to = "logFC") |>
  mutate(contrast = str_remove(contrast, "logFC_")) |>
  filter(!is.na(logFC))

write_csv(lfc_long_all, file.path(DAT_DIR, "panel_B_logfc_long.csv"))

lfc_long <- lfc_long_all |>
  filter(contrast %in% c("Aging", "Training_Young", "Training_Old"))

lfc_long$contrast <- factor(lfc_long$contrast,
                            levels = c("Aging", "Training_Young", "Training_Old"))

set.seed(42)
boot_median_ci <- function(x, R = 2000, conf = 0.95) {
  meds <- replicate(R, median(sample(x, replace = TRUE)))
  qs   <- quantile(meds, c((1 - conf) / 2, (1 + conf) / 2))
  c(lower = unname(qs[1]), upper = unname(qs[2]))
}

lfc_stats <- lfc_long |>
  group_by(contrast) |>
  summarise(
    med_abs_lfc = median(abs(logFC)),
    ci_lo       = boot_median_ci(abs(logFC))[["lower"]],
    ci_hi       = boot_median_ci(abs(logFC))[["upper"]],
    n_above_05  = sum(abs(logFC) > 0.5),
    .groups = "drop"
  )

write_csv(lfc_stats, file.path(DAT_DIR, "panel_B_stats.csv"))

lfc_binwidth <- 4 / 50

lfc_stats$annotation <- sprintf(
  "Med.|logFC| = %.2f [%.2f, %.2f]\nn(>0.5) = %d",
  lfc_stats$med_abs_lfc, lfc_stats$ci_lo, lfc_stats$ci_hi, lfc_stats$n_above_05
)

blunt_file <- "03_DEP/c_data/06_blunting_diagnostics.csv"
power_file <- "03_DEP/c_data/08_power_analysis.csv"
dist_subtitle <- NULL

if (file.exists(blunt_file)) {
  blunt <- read_csv(blunt_file, show_col_types = FALSE)
  ks_row  <- blunt[blunt$test == "Kolmogorov-Smirnov", ]
  fl_row  <- blunt[blunt$test == "Fligner-Killeen", ]

  # Blunting ratio: median |logFC| Training_Old / Training_Young with bootstrap CI
  lfc_ty <- abs(lfc_long$logFC[lfc_long$contrast == "Training_Young"])
  lfc_to <- abs(lfc_long$logFC[lfc_long$contrast == "Training_Old"])
  obs_ratio <- median(lfc_to) / median(lfc_ty)

  set.seed(42)
  boot_ratios <- replicate(2000, {
    median(sample(lfc_to, replace = TRUE)) / median(sample(lfc_ty, replace = TRUE))
  })
  ratio_ci <- quantile(boot_ratios, c(0.025, 0.975))

  power_txt <- ""
  if (file.exists(power_file)) {
    pwr <- read_csv(power_file, show_col_types = FALSE)
    ty_min <- pwr$min_detectable_logFC[pwr$contrast == "Training_Young"]
    to_min <- pwr$min_detectable_logFC[pwr$contrast == "Training_Old"]
    power_txt <- sprintf(" | Min detect: Y=%.2f, O=%.2f", ty_min, to_min)
  }

  n_prot <- length(unique(lfc_long$gene))
  dist_subtitle <- sprintf(
    "n = %s | KS D = %.2f, %s | Blunting = %.2f [%.2f, %.2f]",
    format(n_prot, big.mark = ","),
    ks_row$statistic, fmt_p(ks_row$p_value),
    obs_ratio, ratio_ci[1], ratio_ci[2]
  )
} else {
  cat("  blunting_diagnostics.csv not found -- skipping annotation\n")
}

pB <- ggplot(lfc_long, aes(x = logFC, fill = contrast)) +
  geom_histogram(bins = 50, color = "black", linewidth = 0.2, alpha = 0.85) +
  geom_density(aes(y = after_stat(count) * lfc_binwidth),
               alpha = 0.15, linewidth = 0.5, color = "grey20") +
  geom_vline(xintercept = 0, linetype = "solid", color = "grey50", linewidth = 0.4) +
  geom_text(data = lfc_stats,
            aes(x = -1, y = Inf, label = annotation),
            inherit.aes = FALSE, hjust = -0.01, vjust = 1.15,
            size = scale_text(BASE_COUNT, PD_W) + 0.5,
            color = "grey20", fontface = "bold", lineheight = 0.9) +
  facet_wrap(~ contrast, ncol = 1, scales = "fixed",
             labeller = labeller(contrast = CTR_SHORT)) +
  coord_cartesian(xlim = c(-1, 1)) +
  scale_fill_manual(values = CONTRAST_COLORS[c("Aging", "Training_Young", "Training_Old")]) +
  labs(title = "Effect Size Distribution",
       subtitle = if (!is.null(dist_subtitle)) dist_subtitle else NULL,
       x = expression(bold(log[2]~FC)),
       # Single space reserves axis-title horizontal footprint (~4–5 mm
       # rotated) so B's plot panel left edge co-aligns with E's plot
       # panel (E has the visible "Significant pathways" y-axis title).
       y = " ",
       tag = "B") +
  FIG_THEME + theme(legend.position = "none",
                    strip.text = element_text(face = "bold", size = FIG_STRIP_SIZE, margin = margin(b = 2)),
                    panel.spacing.y = unit(4, "pt"),
                    plot.subtitle = element_text(size = FIG_SUBTITLE_SIZE,
                                                face = "bold.italic", color = "grey40"),
                    axis.text.y = element_text(size = FIG_AXIS_TEXT - 1.5, color = "grey40"),
                    axis.ticks.y = element_blank(),
                    # Match pE's effective horizontal margins so B/E column aligns.
                    # pE inner margin(0, 2, 0, 14) + outer margin(6, 2, 0, 2)
                    # → effective l = 16pt. Also gives "B" tag breathing room
                    # from pA and clearance from the title/subtitle.
                    plot.margin = margin(t = 6, r = 4, b = 0, l = 16))

ggsave(file.path(RPT_PNG, "MAIN_panel_B_logfc_density.png"), pB,
       width = PD_W, height = PD_H, units = "mm", dpi = 300)
ggsave(file.path(RPT_PDF, "MAIN_panel_B_logfc_density.pdf"), pB,
       width = PD_W, height = PD_H, units = "mm", device = pdf_device)
