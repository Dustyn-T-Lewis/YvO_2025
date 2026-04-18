# Supplementary: Soft Threshold Selection (ggplot rewrite)
#
# Replaces the base-R soft threshold plot from YvO_WGCNA_run.R with a
# FIG_THEME-styled ggplot version for publication consistency.
#
# Re-runs pickSoftThreshold on the imputed matrix to recover full fitIndices
# (the pipeline only saves the selected power, not the sweep data).
#
# Generates: 01_QC/SUPP_soft_threshold.png

setwd(rprojroot::find_rstudio_root_file())
source("04_Figures/shared/style.R")

library(WGCNA)
library(readr)
library(dplyr)
library(tidyr)
library(patchwork)
library(ggrepel)

allowWGCNAThreads()

RPT_PNG <- "04_Figures/F06/b_reports/supp/01_QC/png/panels"
RPT_PDF <- "04_Figures/F06/b_reports/supp/01_QC/pdf/panels"
DAT     <- "04_Figures/F06/c_data"
dir.create(RPT_PNG, recursive = TRUE, showWarnings = FALSE)
dir.create(RPT_PDF, recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(DAT, "supp"), recursive = TRUE, showWarnings = FALSE)

pdf_device <- get_pdf_device()

# --- Load imputed data (same preprocessing as YvO_WGCNA_run.R)
df <- read_csv("02_Imputation/c_data/01_imputed.csv", show_col_types = FALSE)
ann_cols   <- c("uniprot_id", "protein", "gene", "description")
samp_names <- setdiff(names(df), ann_cols)
mat        <- as.matrix(df[, samp_names])
rownames(mat) <- df$uniprot_id
datExpr <- t(mat)

gsg <- goodSamplesGenes(datExpr, verbose = 0)
if (!gsg$allOK) datExpr <- datExpr[gsg$goodSamples, gsg$goodGenes]

cor <- WGCNA::cor

message(sprintf("Soft threshold ggplot: %d samples x %d proteins",
                nrow(datExpr), ncol(datExpr)))

# --- Run soft threshold sweep
powers <- 1:20
sft <- pickSoftThreshold(datExpr, powerVector = powers,
                          networkType = "signed", verbose = 2)

# Read selected power from pipeline output
sft_csv <- read_csv(file.path(DAT, "wgcna/wgcna_sft_summary.csv"),
                     show_col_types = FALSE)
soft_power <- sft_csv$selected_power[1]

# Build tidy data
fit_df <- tibble(
  power  = sft$fitIndices$Power,
  r2     = -sign(sft$fitIndices$slope) * sft$fitIndices$SFT.R.sq,
  mean_k = sft$fitIndices$mean.k.,
  slope  = sft$fitIndices$slope
) %>%
  mutate(selected = power == soft_power)

# --- Panel dimensions
PA_W <- 240
PA_H <- 110

txt_axis  <- scale_text(BASE_STAT, PA_W)
txt_title <- scale_text(BASE_GENE, PA_W) * 1.2
txt_label <- scale_text(BASE_GENE, PA_W) * 0.9

# --- Panel 1: Scale Independence (R² vs power)
p1 <- ggplot(fit_df, aes(power, r2)) +
  geom_hline(yintercept = 0.85, linetype = "dashed", color = "grey40",
             linewidth = 0.4) +
  geom_point(aes(fill = selected), shape = 21, size = 2.5, stroke = 0.4,
             color = "black") +
  geom_text_repel(aes(label = power), size = txt_label,
                  color = "grey25", fontface = "bold",
                  max.overlaps = 20, seed = 42,
                  min.segment.length = 0.3, segment.size = 0.25,
                  segment.color = "grey60") +
  scale_fill_manual(values = c("TRUE" = "#D6604D", "FALSE" = "grey70"),
                    guide = "none") +
  scale_x_continuous(breaks = seq(2, 20, 2)) +
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.2)) +
  annotate("text", x = 19, y = 0.88,
           label = expression(R^2 == 0.85~threshold),
           size = txt_label * 0.85, color = "grey40", hjust = 1) +
  labs(x = "Soft Threshold (power)",
       y = expression(Scale~Free~Topology~R^2),
       title = "Scale Independence") +
  FIG_THEME +
  theme(plot.title = element_text(size = txt_title, face = "bold"))

# --- Panel 2: Mean Connectivity
p2 <- ggplot(fit_df, aes(power, mean_k)) +
  geom_point(aes(fill = selected), shape = 21, size = 2.5, stroke = 0.4,
             color = "black") +
  geom_text_repel(aes(label = power), size = txt_label,
                  color = "grey25", fontface = "bold",
                  max.overlaps = 20, seed = 42,
                  min.segment.length = 0.3, segment.size = 0.25,
                  segment.color = "grey60") +
  scale_fill_manual(values = c("TRUE" = "#D6604D", "FALSE" = "grey70"),
                    guide = "none") +
  scale_x_continuous(breaks = seq(2, 20, 2)) +
  labs(x = "Soft Threshold (power)",
       y = "Mean Connectivity",
       title = "Mean Connectivity") +
  FIG_THEME +
  theme(plot.title = element_text(size = txt_title, face = "bold"))

# --- Composite
composite <- (p1 | p2) +
  plot_annotation(
    title = "Scale-Free Topology Fit",
    subtitle = sprintf("Signed network | selected power = %d (R\u00b2 = %.3f) | %s proteins",
                       soft_power, fit_df$r2[fit_df$power == soft_power],
                       format(ncol(datExpr), big.mark = ",")),
    theme = theme(
      plot.title    = element_text(face = "bold", size = FIG_TITLE_SIZE),
      plot.subtitle = element_text(face = "bold.italic", size = FIG_SUBTITLE_SIZE,
                                    color = "grey30"),
      plot.margin   = margin(4, 4, 4, 4)
    )
  )

ggsave(file.path(RPT_PNG, "SUPP_soft_threshold.png"), composite,
       width = PA_W, height = PA_H, units = "mm",
       dpi = 300, limitsize = FALSE)
ggsave(file.path(RPT_PDF, "SUPP_soft_threshold.pdf"), composite,
       width = PA_W, height = PA_H, units = "mm",
       device = pdf_device, limitsize = FALSE)

# Save fitIndices for future use (avoids re-running pickSoftThreshold)
write_csv(fit_df, file.path(DAT, "supp", "a05_sft_fit_indices.csv"))

message("  Soft threshold ggplot saved")
