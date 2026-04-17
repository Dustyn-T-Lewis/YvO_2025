# F00_QC Stage 01 — Normalization detail composite
# 2x2 grid, 200x160 mm. Goes deeper than the MAIN_F00_QC overview on the
# normalization stage specifically. Self-contained; reads only the stage 01
# report-intermediates RDS.
# Panels:
#   A Filter cascade (waterfall)       B Per-protein missingness distribution
#   C Eta-squared distribution         D Outlier consensus diagnostic
# Output: 04_Figures/F00/b_reports/MAIN_F00_QC_stage01_normalization.{pdf,png}

setwd(rprojroot::find_rstudio_root_file())
source("04_Figures/shared/style.R")

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(ggplot2)
  library(patchwork)
  library(scales)
})

RPT_DIR <- "04_Figures/F00/b_reports"
dir.create(RPT_DIR, recursive = TRUE, showWarnings = FALSE)
pdf_device <- get_pdf_device()

int01 <- readRDS("01_normalization/c_data/00_report_intermediates.rds")

# --- Panel A: filter cascade waterfall ---
fcasc <- int01$filter_log %>%
  mutate(step   = factor(step, levels = step),
         retained = n_after,
         removed  = ifelse(is.na(n_removed), 0, n_removed))

pA <- ggplot(fcasc, aes(step, retained)) +
  geom_col(fill = "#2166AC", width = 0.7) +
  geom_text(aes(label = format(retained, big.mark = ",")),
            vjust = -1.8, size = 2.8, fontface = "bold") +
  geom_text(aes(label = ifelse(removed > 0,
                               sprintf("-%s", format(removed, big.mark=",")),
                               "")),
            vjust = -0.4, size = 2.4, color = "#B2182B") +
  scale_y_continuous(expand = expansion(mult = c(0, 0.18)),
                     labels = comma) +
  labs(title = "Filter cascade",
       subtitle = sprintf("%s raw -> %s retained (%.1f%%)",
                          format(fcasc$retained[1], big.mark=","),
                          format(tail(fcasc$retained,1), big.mark=","),
                          tail(fcasc$pct_of_raw,1)),
       x = NULL, y = "Proteins retained", tag = "A") +
  FIG_THEME +
  theme(axis.text.x = element_text(angle = 22, hjust = 1, size = 7))

# --- Panel B: per-protein missingness distribution ---
# Compute % missing per protein from pre-normalization matrix.
prot_miss <- rowSums(is.na(int01$data_pre_outlier)) / ncol(int01$data_pre_outlier) * 100
miss_df <- tibble(pct_miss = prot_miss)

# Filter threshold: "Missingness (>=10 reps in >=1 group)" is the last filter step.
# Approx threshold line at 100*(1 - 10/16) = 37.5% (retained if any group has >=10 of 16).
# More accurate: show where the filter cut happens empirically (max value among retained).
thresh_pct <- quantile(prot_miss, 0.99, na.rm = TRUE)  # visual reference

pB <- ggplot(miss_df, aes(pct_miss)) +
  geom_histogram(binwidth = 2.5, fill = "#4393C3", color = "white", alpha = 0.85) +
  geom_vline(xintercept = thresh_pct, linetype = "dashed", color = "#B2182B", linewidth = 0.4) +
  annotate("text", x = thresh_pct, y = Inf,
           label = sprintf("99th pct = %.1f%%", thresh_pct),
           vjust = 1.6, hjust = 1.05, color = "#B2182B", size = 2.7) +
  scale_y_continuous(labels = comma) +
  labs(title = "Per-protein missingness (post-filter)",
       subtitle = sprintf("%s proteins | median = %.1f%%",
                          format(nrow(miss_df), big.mark=","),
                          median(prot_miss, na.rm=TRUE)),
       x = "% missing across samples", y = "Proteins", tag = "B") +
  FIG_THEME

# --- Panel C: eta-squared distribution ---
eta_df <- tibble(eta2 = as.numeric(int01$eta2_vals))
eta_med <- median(eta_df$eta2, na.rm = TRUE)
eta_90  <- quantile(eta_df$eta2, 0.90, na.rm = TRUE)

pC <- ggplot(eta_df, aes(eta2)) +
  geom_density(fill = "#4CAF50", alpha = 0.55, linewidth = 0.4) +
  geom_vline(xintercept = eta_med, linetype = "solid", color = "grey30",
             linewidth = 0.4) +
  geom_vline(xintercept = eta_90, linetype = "dashed", color = "#E05A4E",
             linewidth = 0.4) +
  annotate("text", x = eta_med, y = Inf, label = sprintf("median = %.2f", eta_med),
           vjust = 1.5, hjust = -0.1, size = 2.7, color = "grey20") +
  annotate("text", x = eta_90, y = Inf, label = sprintf("90th = %.2f", eta_90),
           vjust = 3.0, hjust = -0.1, size = 2.7, color = "#E05A4E") +
  labs(title = expression("Biological signal ("*eta^2*")"),
       subtitle = sprintf("%s proteins | primary factor explained variance",
                          format(nrow(eta_df), big.mark=",")),
       x = expression(eta^2), y = "Density", tag = "C") +
  FIG_THEME

# --- Panel D: outlier consensus diagnostic ---
od <- int01$outlier_diag %>%
  mutate(age = ifelse(grepl("^O", prefix), "Old", "Young"),
         age = factor(age, levels = c("Young", "Old")),
         consensus_lab = ifelse(consensus_outlier, "Consensus outlier",
                                ifelse(n_flags > 0, "Single flag", "Clean")))

pD <- ggplot(od, aes(mahal_dist, n_flags,
                     color = consensus_lab, shape = Timepoint)) +
  geom_jitter(width = 0.04, height = 0.12, size = 2.2, alpha = 0.85) +
  scale_color_manual(values = c("Clean" = "grey55",
                                "Single flag" = "#E6A100",
                                "Consensus outlier" = "#B2182B"),
                     name = NULL) +
  scale_shape_manual(values = SHAPE_TP, guide = "none") +
  scale_y_continuous(breaks = 0:4) +
  labs(title = "Outlier consensus diagnostic",
       subtitle = sprintf("%d samples | %d consensus outliers",
                          nrow(od), sum(od$consensus_outlier)),
       x = "Mahalanobis distance (PC1-3)",
       y = "Number of QC flags",
       tag = "D") +
  FIG_THEME +
  theme(legend.position = "top",
        legend.key.size = unit(3, "mm"),
        legend.text = element_text(size = 7))

# --- Composite ---
composite <- (pA | pB) / (pC | pD) +
  plot_annotation(
    title = "YvO Proteomics Pipeline - Stage 01: Normalization detail",
    subtitle = sprintf(
      "Filter cascade | per-protein missingness | biological signal | outlier consensus | %s",
      format(Sys.Date(), "%Y-%m-%d")),
    theme = theme(plot.title    = element_text(face = "bold", size = 13),
                  plot.subtitle = element_text(face = "italic", size = 9,
                                               color = "grey30"))
  ) &
  theme(plot.tag = element_text(face = "bold", size = 14))

COMP_W <- 200
COMP_H <- 160

ggsave(file.path(RPT_DIR, "MAIN_F00_QC_stage01_normalization.pdf"),
       composite, width = COMP_W, height = COMP_H, units = "mm",
       device = pdf_device)
ggsave(file.path(RPT_DIR, "MAIN_F00_QC_stage01_normalization.png"),
       composite, width = COMP_W, height = COMP_H, units = "mm", dpi = 300)

message(sprintf("Stage 01 composite saved -> %s", RPT_DIR))
