# Sourced by 02_supp_panels.R — expects style.R already loaded.

library(readr)
library(dplyr)
library(tidyr)
library(patchwork)
library(ggrepel)

BASE <- here::here("04_Figures", "F06")

RPT_PNG <- file.path(BASE, "b_reports", "supp", "png", "panels")
RPT_PDF <- file.path(BASE, "b_reports", "supp", "pdf", "panels")
DAT     <- file.path(BASE, "c_data")
dir.create(RPT_PNG, recursive = TRUE, showWarnings = FALSE)
dir.create(RPT_PDF, recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(DAT, "supp"), recursive = TRUE, showWarnings = FALSE)

pdf_device <- get_pdf_device()

sft_fi <- readRDS(file.path(DAT, "wgcna/sft_fitIndices.rds"))

sft_csv <- read_csv(file.path(DAT, "wgcna/wgcna_sft_summary.csv"))
soft_power <- sft_csv$selected_power[1]
n_proteins <- sft_csv$n_proteins[1]

message(sprintf("Soft threshold ggplot: loading cached fitIndices (%d powers, %d proteins)",
                nrow(sft_fi), n_proteins))

fit_df <- tibble(
  power  = sft_fi$Power,
  r2     = -sign(sft_fi$slope) * sft_fi$SFT.R.sq,
  mean_k = sft_fi$mean.k.,
  slope  = sft_fi$slope
) |>
  mutate(selected = power == soft_power)

PA_W <- 240
PA_H <- 110

txt_axis  <- scale_text(BASE_STAT, PA_W)
txt_title <- scale_text(BASE_GENE, PA_W) * 1.6
txt_label <- scale_text(BASE_GENE, PA_W) * 0.9

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
  theme(plot.title = element_text(size = 10, face = "bold"))

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
  theme(plot.title = element_text(size = 10, face = "bold"))

composite <- (p1 | p2) +
  plot_annotation(
    title = "Scale-Free Topology Fit",
    subtitle = sprintf("Signed network | selected power = %d (R\u00b2 = %.3f) | %s proteins",
                       soft_power, fit_df$r2[fit_df$power == soft_power],
                       format(n_proteins, big.mark = ",")),
    theme = theme(
      plot.title    = element_text(face = "bold", size = 13),
      plot.subtitle = element_text(face = "bold.italic", size = 10,
                                    color = "grey30"),
      plot.margin   = margin(4, 4, 4, 4)
    )
  )

ggsave(file.path(RPT_PNG, "SUPP_soft_threshold.png"), composite,
       width = PA_W, height = PA_H, units = "mm", dpi = 300)
ggsave(file.path(RPT_PDF, "SUPP_soft_threshold.pdf"), composite,
       width = PA_W, height = PA_H, units = "mm", device = pdf_device)

write_csv(fit_df, file.path(DAT, "supp", "a05_sft_fit_indices.csv"))

message("  Soft threshold ggplot saved")
