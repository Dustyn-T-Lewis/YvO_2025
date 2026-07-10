#!/usr/bin/env Rscript
# F03 Supp — 5-panel diagnostic composite (2×2 + full-row)
# A: p-value histograms  B: Pi-score distributions  C: FDR distributions
# D: MA plots  E: outlier sensitivity (full bottom row)
#
# Panels A-C share the same histogram template (different column/color/label).

setwd(here::here())

pacman::p_load(dplyr, tidyr, readr, readxl, ggplot2, patchwork, cowplot)

source("04_Figures/shared/style.R")

BASE    <- "04_Figures/F03"
RPT_PNG <- file.path(BASE, "b_reports", "supp", "png", "panels")
RPT_PDF <- file.path(BASE, "b_reports", "supp", "pdf", "panels")
DAT     <- file.path(BASE, "c_data", "supp")
for (d in c(RPT_PNG, RPT_PDF, DAT)) dir.create(d, recursive = TRUE, showWarnings = FALSE)

pdf_dev <- get_pdf_device()
DEP_XLSX <- "03_DEP/c_data/03_DEP_results.xlsx"
CTRS <- c("Aging", "Training_Young", "Training_Old", "Interaction")

per_contrast <- lapply(CTRS, \(ctr) {
  as.data.frame(read_excel(DEP_XLSX, sheet = ctr))
})
names(per_contrast) <- CTRS

# Parametric histogram builder (collapses 3 near-identical panels)

make_dist_panel <- function(col, fill, xlab, vline = NULL, stat_fmt, title, tag) {
  hist_df <- bind_rows(lapply(CTRS, \(ctr) {
    tibble(contrast = ctr, value = per_contrast[[ctr]][[col]])
  })) |> filter(!is.na(value)) |>
    mutate(contrast = factor(contrast, levels = CTRS))

  n_sig <- hist_df |>
    summarise(n_sig = sum(value < 0.05), .by = contrast)

  n_bins <- 20
  p <- ggplot(hist_df, aes(value)) +
    geom_histogram(breaks = seq(0, 1, length.out = n_bins + 1),
                   fill = fill, color = "white", linewidth = 0.3)

  if (!is.null(vline))
    p <- p + geom_vline(xintercept = vline, linetype = "dashed",
                        color = "grey30", linewidth = 0.4)
  if (col == "P.Value") {
    uniform_ref <- hist_df |> summarise(h = n() / n_bins, .by = contrast)
    p <- p + geom_hline(data = uniform_ref, aes(yintercept = h),
                        linetype = "dashed", color = "grey40", linewidth = 0.4)
  }

  p <- p +
    geom_text(data = n_sig,
              aes(x = 0.5, y = Inf, label = CTR_SHORT[as.character(contrast)]),
              inherit.aes = FALSE, hjust = 0.5, vjust = 1.2,
              size = 2.0, fontface = "bold", color = "grey20") +
    geom_text(data = n_sig,
              aes(x = 0.5, y = Inf, label = sprintf(stat_fmt, n_sig)),
              inherit.aes = FALSE, hjust = 0.5, vjust = 3.8,
              size = 2.8, fontface = "bold", color = "grey40") +
    facet_wrap(~contrast, ncol = 2, scales = "free_y",
               labeller = labeller(contrast = CTR_SHORT)) +
    labs(title = title,
         subtitle = sprintf("%s proteins | 20 bins",
                            format(round(nrow(hist_df) / length(CTRS)), big.mark = ",")),
         x = xlab, y = "Proteins", tag = tag) +
    FIG_THEME + theme(strip.background = element_blank(), strip.text = element_blank())

  write_csv(hist_df, file.path(DAT, sprintf("panel_%s.csv", gsub("[. ]", "_", col))))
  ggsave(file.path(RPT_PNG, sprintf("SUPP_%s.png", gsub("[. ]", "_", title))), p,
         width = 89, height = 75, units = "mm", dpi = 300)
  p
}

pB <- make_dist_panel("P.Value", "#5DA5DA", "p_value", NULL,
                      "p < 0.05: %d", "Raw p-value distribution", "a")
pB_title <- "Raw p-value distribution"

pPi <- make_dist_panel("pi_score", "#E05A4E", "\u03A0-score", 0.05,
                       "\u03A0 < 0.05: %d", "\u03A0-score distribution", "b")
pPi_title <- "\u03A0-score distribution"

pFDR <- make_dist_panel("adj.P.Val", "#9B7FBF", "FDR (BH)", 0.05,
                        "FDR < 0.05: %d", "FDR distribution", "c")
pFDR_title <- "FDR distribution"

pB   <- strip_for_composite(pB)
pPi  <- strip_for_composite(pPi)
pFDR <- strip_for_composite(pFDR)

ma_df <- bind_rows(lapply(CTRS, \(ctr) {
  per_contrast[[ctr]] |>
    transmute(contrast = ctr, gene, average_intensity, logFC,
              sig_pi, direction = case_when(
                sig_pi ==  1 ~ "Up", sig_pi == -1 ~ "Down", TRUE ~ "NS"))
})) |>
  mutate(contrast  = factor(contrast, levels = CTRS),
         direction = factor(direction, levels = c("Up", "Down", "NS")))

n_dep <- ma_df |> filter(direction != "NS") |> count(contrast, name = "n_dep")

pC <- ggplot(ma_df, aes(average_intensity, logFC, color = direction)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey40", linewidth = 0.3) +
  geom_point(data = \(d) filter(d, direction == "NS"), alpha = 0.25, size = 0.6) +
  geom_point(data = \(d) filter(d, direction != "NS"), alpha = 0.85, size = 0.9) +
  geom_text(data = n_dep, aes(x = Inf, y = Inf, label = sprintf("\u03A0: %d", n_dep)),
            inherit.aes = FALSE, hjust = 1.05, vjust = 1.5, size = 2.2,
            fontface = "bold", color = "grey20") +
  scale_color_manual(values = DIR_COLORS, name = "\u03A0 < 0.05") +
  facet_wrap(~contrast, ncol = 2, labeller = labeller(contrast = CTR_SHORT)) +
  labs(title = "MA plots", x = "Mean log2 intensity", y = "logFC", tag = "d") +
  FIG_THEME + theme(strip.background = element_blank(), strip.text = element_blank(),
                    legend.position = "top", legend.key.size = unit(3, "mm"))
write_csv(ma_df, file.path(DAT, "panel_C_ma.csv"))
ggsave(file.path(RPT_PNG, "SUPP_panel_C_ma.png"), pC, width = 89, height = 75,
       units = "mm", dpi = 300)
pC_title <- "MA plots"
pC <- strip_for_composite(pC)

out_sens <- tryCatch(
  as.data.frame(read_excel(DEP_XLSX, sheet = "outlier_sensitivity")),
  error = \(e) NULL)

if (!is.null(out_sens) && nrow(out_sens) > 0) {
  long_df <- out_sens |>
    mutate(Contrast = factor(Contrast, levels = CTRS)) |>
    pivot_longer(c(FDR_full, FDR_reduced, Pi_full, Pi_reduced),
                 names_to = "metric_cohort", values_to = "n") |>
    mutate(metric = factor(sub("_.*", "", metric_cohort),
                           levels = c("FDR", "Pi"), labels = c("FDR < 0.05", "\u03A0 < 0.05")),
           cohort = factor(sub(".*_", "", metric_cohort), levels = c("full", "reduced")))

  pE <- ggplot(long_df, aes(cohort, n, fill = cohort)) +
    geom_col(width = 0.65) +
    geom_text(aes(label = n), vjust = -0.3, size = 2.8, fontface = "bold") +
    facet_grid(metric ~ Contrast, scales = "free_y", switch = "y",
               labeller = labeller(Contrast = CTR_SHORT)) +
    scale_fill_manual(values = c(full = "#2166AC", reduced = "#B2182B"),
                      labels = c("Full", "Outlier-removed"), name = NULL) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.22))) +
    labs(title = "DEP retention (outlier removal)", x = NULL, y = "DEPs", tag = "e") +
    FIG_THEME + theme(strip.text.x = element_blank(),
                      strip.text.y = element_text(face = "bold", size = FIG_STRIP_SIZE - 1),
                      legend.position = "top", legend.key.size = unit(3, "mm"))
} else {
  stop("outlier_sensitivity sheet missing from 03_DEP_results.xlsx — ",
       "run 03_DEP/a_script/02_generate_reports.R first")
}
pE_title <- "DEP retention (outlier removal)"
pE <- strip_for_composite(pE)

SUPP_PDF <- file.path(BASE, "b_reports", "supp", "pdf")
SUPP_PNG <- file.path(BASE, "b_reports", "supp", "png")
dir.create(SUPP_PDF, recursive = TRUE, showWarnings = FALSE)
dir.create(SUPP_PNG, recursive = TRUE, showWarnings = FALSE)

COMP_W <- 178; COMP_H <- 225
txt <- composite_text_sizes(COMP_H)
grid <- (pB | pPi) / (pFDR | pC) / pE &
  theme(plot.margin = margin(9, 4, 4, 4))

X_L <- 0.012; X_R <- 0.512; X_TTL <- 0.029; SUB_OFF <- 0.014
Y_R1 <- 0.974; Y_R2 <- 0.651; Y_R3 <- 0.321

all_titles <- list(pB_title, pPi_title, pFDR_title, pC_title, pE_title)
all_tags <- LETTERS[1:5]
all_xs <- c(X_L, X_R, X_L, X_R, X_L)
all_ys <- c(Y_R1, Y_R1, Y_R2, Y_R2, Y_R3)

composite <- ggdraw(grid)
for (i in seq_along(all_tags)) {
  composite <- composite +
    draw_label(all_tags[i], x = all_xs[i], y = all_ys[i],
               fontface = "bold", size = txt$tag, hjust = 0, vjust = 1) +
    draw_label(all_titles[[i]], x = all_xs[i] + X_TTL, y = all_ys[i],
               fontface = "bold", size = txt$title, hjust = 0, vjust = 1)
}

ggsave(file.path(SUPP_PDF, "SUPP_F03_composite.pdf"), composite,
       width = COMP_W, height = COMP_H, units = "mm", device = pdf_dev)
ggsave(file.path(SUPP_PNG, "SUPP_F03_composite.png"), composite,
       width = COMP_W, height = COMP_H, units = "mm", dpi = 300)

message("F03 supp composite done")
