# F00_QC Stage 03 — DEP detail composite
# 2x2 grid, 200x160 mm. Covers stage 03 diagnostics NOT already in F03 SUPP
# (which handles p-hist, MA, imputation sensitivity, outlier sensitivity).
# Panels:
#   A DEP count matrix heatmap          B Effect-size bootstrap forest
#   C Power analysis per contrast       D Training blunting scatter
# Output: 04_Figures/F00/b_reports/MAIN_F00_QC_stage03_DEP.{pdf,png}

setwd(rprojroot::find_rstudio_root_file())
source("04_Figures/shared/style.R")

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(readr)
  library(ggplot2)
  library(patchwork)
  library(scales)
})

RPT_DIR <- "04_Figures/F00/b_reports"
dir.create(RPT_DIR, recursive = TRUE, showWarnings = FALSE)
pdf_device <- get_pdf_device()

# --- Load ---
da_summ <- read_csv("03_DEP/c_data/02_DA_summary.csv", show_col_types = FALSE)
eff     <- read_csv("03_DEP/c_data/07_effect_size_bootstrap.csv", show_col_types = FALSE)
pwr     <- read_csv("03_DEP/c_data/08_power_analysis.csv", show_col_types = FALSE)
blunt   <- read_csv("03_DEP/c_data/06_blunting_diagnostics.csv", show_col_types = FALSE)
combo   <- read_csv("03_DEP/c_data/03_combined_results.csv", show_col_types = FALSE)

CONTRAST_ORDER <- c("Aging", "Training_Young", "Training_Old",
                    "Interaction", "Reversal")

# --- Panel A: DEP count matrix (5 contrasts x 4 thresholds) ---
# Sum up+down (exclude nonsig). Thresholds: p<0.05 (sig.PVal from nominal),
# FDR<0.10, FDR<0.05, Pi<0.05.
counts <- da_summ %>%
  filter(type %in% c("up", "down")) %>%
  group_by(contrast) %>%
  summarise(`p<0.05`    = sum(sig.PVal),
            `FDR<0.10`  = sum(sig.FDR.10),
            `FDR<0.05`  = sum(sig.FDR.05),
            `Pi<0.05`   = sum(sig.Pi),
            .groups = "drop") %>%
  pivot_longer(-contrast, names_to = "threshold", values_to = "n") %>%
  mutate(contrast  = factor(contrast,
                            levels = intersect(CONTRAST_ORDER,
                                               unique(contrast))),
         threshold = factor(threshold,
                            levels = c("p<0.05", "FDR<0.10",
                                       "FDR<0.05", "Pi<0.05")))

pA <- ggplot(counts, aes(threshold, contrast, fill = n)) +
  geom_tile(color = "white", linewidth = 0.8) +
  geom_text(aes(label = format(n, big.mark = ","),
                color = n > max(n) * 0.5),
            size = 2.8, fontface = "bold") +
  scale_color_manual(values = c("TRUE" = "white", "FALSE" = "grey20"),
                     guide = "none") +
  scale_fill_gradient(low = "#DEEBF7", high = "#08519C",
                      name = "DEPs") +
  labs(title = "DEP counts: contrast x threshold",
       subtitle = "up + down (nonsig excluded)",
       x = NULL, y = NULL, tag = "A") +
  FIG_THEME +
  theme(axis.text.x = element_text(angle = 15, hjust = 1, size = 7),
        axis.text.y = element_text(size = 7),
        legend.position = "right",
        legend.key.size = unit(3, "mm"),
        legend.text = element_text(size = 7))

# --- Panel B: effect-size bootstrap forest ---
eff2 <- eff %>%
  mutate(contrast = factor(contrast,
                           levels = rev(intersect(CONTRAST_ORDER,
                                                  unique(contrast)))))

pB <- ggplot(eff2, aes(median_absLFC, contrast, color = contrast)) +
  geom_errorbar(aes(xmin = ci_lower, xmax = ci_upper),
                width = 0.18, linewidth = 0.6, orientation = "y") +
  geom_point(size = 3) +
  geom_text(aes(label = sprintf("%.3f", median_absLFC)),
            vjust = -1.0, size = 2.6, color = "grey20") +
  scale_color_manual(values = CONTRAST_COLORS, guide = "none") +
  scale_x_continuous(expand = expansion(mult = c(0.05, 0.10))) +
  labs(title = "Effect-size bootstrap",
       subtitle = sprintf("median |logFC| with 95%% CI (n=%s)",
                          format(as.integer(median(eff$n_proteins, na.rm = TRUE)),
                                 big.mark = ",")),
       x = "median |logFC| (log2)", y = NULL, tag = "B") +
  FIG_THEME +
  theme(axis.text.y = element_text(size = 7))

# --- Panel C: power analysis per contrast ---
pwr2 <- pwr %>%
  mutate(contrast = factor(contrast,
                           levels = rev(intersect(CONTRAST_ORDER,
                                                  unique(contrast)))))

pC <- ggplot(pwr2, aes(min_detectable_d, contrast, fill = contrast)) +
  geom_col(width = 0.65) +
  geom_text(aes(label = sprintf("d=%.2f, power=%.2f, n=%d",
                                min_detectable_d, power, n_subjects)),
            hjust = -0.05, size = 2.5, color = "grey20") +
  scale_fill_manual(values = CONTRAST_COLORS, guide = "none") +
  scale_x_continuous(expand = expansion(mult = c(0, 0.38))) +
  labs(title = "Power analysis",
       subtitle = sprintf("min detectable Cohen d | alpha=%.2f | power=%.2f",
                          median(pwr$alpha, na.rm = TRUE),
                          median(pwr$power, na.rm = TRUE)),
       x = "min detectable d", y = NULL, tag = "C") +
  FIG_THEME +
  theme(axis.text.y = element_text(size = 7))

# --- Panel D: training blunting scatter ---
# Use combined results (wide). Columns: logFC_Training_Young / _Training_Old,
# FDR_Training_Young / _Training_Old.
need_cols <- c("logFC_Training_Young",   "logFC_Training_Old",
               "adj.P.Val_Training_Young", "adj.P.Val_Training_Old")
miss <- setdiff(need_cols, colnames(combo))
if (length(miss)) stop("Missing columns in 03_combined_results.csv: ",
                       paste(miss, collapse = ", "))

blunt_df <- combo %>%
  transmute(logFC_TY   = .data$logFC_Training_Young,
            logFC_TO   = .data$logFC_Training_Old,
            fdr_TY     = .data$adj.P.Val_Training_Young,
            fdr_TO     = .data$adj.P.Val_Training_Old) %>%
  filter(!is.na(logFC_TY), !is.na(logFC_TO)) %>%
  mutate(abs_TY = abs(logFC_TY),
         abs_TO = abs(logFC_TO),
         panel  = case_when(
           fdr_TY < 0.10 & fdr_TO < 0.10 ~ "Sig both",
           fdr_TY < 0.10                 ~ "Sig Young only",
           fdr_TO < 0.10                 ~ "Sig Old only",
           TRUE                          ~ "NS"))

max_lim <- max(c(blunt_df$abs_TY, blunt_df$abs_TO), na.rm = TRUE) * 1.05

ks_row  <- blunt %>% filter(grepl("Kolmogorov", test)) %>% slice_head(n = 1)
wil_row <- blunt %>% filter(grepl("Wilcoxon",  test)) %>% slice_head(n = 1)
cliff   <- blunt %>% filter(grepl("Cliff",     test)) %>% slice_head(n = 1)

annot_txt <- sprintf("KS p=%.1e | Wilcoxon p=%.1e\nCliff's delta = %.2f",
                     ks_row$p_value[1], wil_row$p_value[1],
                     cliff$statistic[1])

pD <- ggplot(blunt_df, aes(abs_TY, abs_TO, color = panel)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed",
              color = "grey55", linewidth = 0.4) +
  geom_point(alpha = 0.35, size = 0.85) +
  scale_color_manual(values = SIG_COLORS_F2, name = NULL) +
  scale_x_continuous(limits = c(0, max_lim), expand = c(0, 0)) +
  scale_y_continuous(limits = c(0, max_lim), expand = c(0, 0)) +
  coord_equal() +
  annotate("text", x = max_lim * 0.98, y = max_lim * 0.04,
           label = annot_txt, hjust = 1, vjust = 0,
           size = 2.5, color = "grey20") +
  labs(title = "Training blunting",
       subtitle = "|logFC| Young vs Old",
       x = "|logFC| Training_Young",
       y = "|logFC| Training_Old", tag = "D") +
  FIG_THEME +
  theme(legend.position = "top",
        legend.key.size = unit(3, "mm"),
        legend.text = element_text(size = 6.5))

# --- Composite ---
composite <- (pA | pB) / (pC | pD) +
  plot_annotation(
    title = "YvO Proteomics Pipeline - Stage 03: DEP detail",
    subtitle = sprintf(
      "DEP counts | effect size | power | training blunting | %s",
      format(Sys.Date(), "%Y-%m-%d")),
    theme = theme(plot.title    = element_text(face = "bold", size = 13),
                  plot.subtitle = element_text(face = "italic", size = 9,
                                               color = "grey30"))
  ) &
  theme(plot.tag = element_text(face = "bold", size = 14))

COMP_W <- 200
COMP_H <- 160

ggsave(file.path(RPT_DIR, "MAIN_F00_QC_stage03_DEP.pdf"),
       composite, width = COMP_W, height = COMP_H, units = "mm",
       device = pdf_device)
ggsave(file.path(RPT_DIR, "MAIN_F00_QC_stage03_DEP.png"),
       composite, width = COMP_W, height = COMP_H, units = "mm", dpi = 300)

message(sprintf("Stage 03 composite saved -> %s", RPT_DIR))
