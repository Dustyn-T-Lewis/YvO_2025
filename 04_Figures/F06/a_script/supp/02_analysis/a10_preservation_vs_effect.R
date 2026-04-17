# F06 SUPP — Preservation (Zsummary) vs strongest LMM effect size per module
#
# Integrates the two margins of Panel A (left = protein count, right = preservation)
# with the heatmap body (LMM r_equiv) into a single module-level scatter.
# Highlights modules that are BOTH robust (Zsummary > 10) and impactful (|r| large).
#
# Output: 04_Figures/F06/b_reports/supp/02_analysis/SUPP_preservation_vs_effect.{pdf,png}

setwd(rprojroot::find_rstudio_root_file())
source("04_Figures/shared/style.R")

library(readr)
library(dplyr)
library(ggplot2)
library(ggrepel)
library(stringr)

DAT <- "04_Figures/F06/c_data"
RPT <- "04_Figures/F06/b_reports/supp/02_analysis"
dir.create(RPT, recursive = TRUE, showWarnings = FALSE)
pdf_device <- get_pdf_device()

pres <- read_csv(file.path(DAT, "05_panel_E_preservation.csv"),
                 show_col_types = FALSE)
lmm  <- read_csv(file.path(DAT, "wgcna/wgcna_lmm_contrast_audit.csv"),
                 show_col_types = FALSE)
mod_bio <- read_csv(file.path(DAT, "mod_bio_labels.csv"),
                    show_col_types = FALSE)

# Strongest |r_equiv| per module across contrasts
eff <- lmm %>%
  filter(module != "MEgrey") %>%
  mutate(module_color = sub("^ME", "", module)) %>%
  group_by(module_color) %>%
  slice_max(abs(r_equiv), n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  dplyr::select(module_color, top_contrast = contrast,
                top_r = r_equiv, top_p_bh = p_bh)

df <- pres %>%
  filter(module != "grey") %>%
  left_join(eff, by = c("module" = "module_color")) %>%
  mutate(abs_r = abs(top_r),
         sig   = ifelse(!is.na(top_p_bh) & top_p_bh < 0.05, "BH < 0.05", "n.s."),
         display_lab = ifelse(is.na(bio_label) | bio_label == "" |
                                grepl("^M[0-9]+: NA$", bio_label),
                              str_to_title(module),
                              sub("^M[0-9]+:\\s*", "", bio_label)),
         label = sprintf("%s\n(%s: %+.2f)",
                         display_lab, top_contrast, top_r))

p <- ggplot(df, aes(x = Zsummary, y = abs_r)) +
  annotate("rect", xmin = 10, xmax = Inf, ymin = 0.4, ymax = Inf,
           fill = "#D6604D", alpha = 0.10) +
  annotate("text", x = Inf, y = Inf,
           label = "Robust &\nimpactful",
           hjust = 1.05, vjust = 1.3,
           size = 2.8, fontface = "bold.italic", color = "grey25") +
  geom_vline(xintercept = c(2, 10), linetype = "dashed",
             color = "grey50", linewidth = 0.3) +
  geom_hline(yintercept = 0.4, linetype = "dashed",
             color = "grey50", linewidth = 0.3) +
  geom_point(aes(fill = module, size = module_size, shape = sig),
             color = "black", stroke = 0.4) +
  geom_text_repel(aes(label = label), size = 2.3, lineheight = 0.85,
                  min.segment.length = 0, segment.size = 0.2,
                  max.overlaps = Inf, box.padding = 0.4) +
  scale_fill_identity() +
  scale_shape_manual(values = c("BH < 0.05" = 21, "n.s." = 24),
                     name = "Strongest LMM") +
  scale_size_continuous(range = c(3, 8), name = "Protein count") +
  scale_x_continuous(limits = c(0, max(df$Zsummary) * 1.1),
                     breaks = c(0, 2, 10, 20, 40)) +
  scale_y_continuous(limits = c(0, max(df$abs_r, 0.8) * 1.1)) +
  labs(title    = "Module Preservation vs. Strongest Effect Size",
       subtitle = "Zsummary (Langfelder 2011) | |r_equiv| from LMM contrast audit",
       x = "Preservation (Zsummary)",
       y = "|r_equiv| of strongest contrast") +
  FIG_THEME +
  theme(legend.position = "right",
        legend.key.size = unit(3, "mm"))

W <- 180; H <- 140
ggsave(file.path(RPT, "SUPP_preservation_vs_effect.png"), p,
       width = W, height = H, units = "mm", dpi = 300)
ggsave(file.path(RPT, "SUPP_preservation_vs_effect.pdf"), p,
       width = W, height = H, units = "mm", device = pdf_device)

message("  SUPP preservation-vs-effect saved")
