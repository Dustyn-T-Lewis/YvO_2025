# F03 SUPP — Panel C: MA plots per contrast
# One facet per contrast — logFC against mean log2 intensity. Points colored by
# Pi-score DEP flag. Provides mean-intensity dependence diagnostic complementary
# to the MAIN volcano rings.
#
# Reads:  03_DEP/c_data/04_per_contrast_results/{Aging,Training_Young,Training_Old,Interaction}.csv
# Writes: b_reports/supp/SUPP_panel_C_ma.png
#         c_data/supp/panel_C_ma.csv

setwd(rprojroot::find_rstudio_root_file())
source("04_Figures/shared/style.R")

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(ggplot2)
})

CTRS    <- c("Aging", "Training_Young", "Training_Old", "Interaction")
RPT_DIR <- "04_Figures/F03/b_reports/supp/panels"
DAT_DIR <- "04_Figures/F03/c_data/supp"
dir.create(RPT_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(DAT_DIR, recursive = TRUE, showWarnings = FALSE)

ma_df <- lapply(CTRS, function(ctr) {
  read_csv(sprintf("03_DEP/c_data/04_per_contrast_results/%s.csv", ctr),
           show_col_types = FALSE) %>%
    transmute(contrast = ctr,
              gene,
              average_intensity,
              logFC,
              adj.P.Val,
              pi_score,
              sig_pi,
              direction = dplyr::case_when(
                sig_pi ==  1 & logFC > 0 ~ "Up",
                sig_pi == -1 & logFC < 0 ~ "Down",
                TRUE                     ~ "NS"))
}) %>% bind_rows() %>%
  mutate(contrast  = factor(contrast, levels = CTRS),
         direction = factor(direction, levels = c("Up", "Down", "NS")))

write_csv(ma_df, file.path(DAT_DIR, "panel_C_ma.csv"))

n_dep <- ma_df %>%
  filter(direction != "NS") %>%
  count(contrast, name = "n_dep")

pC <- ggplot(ma_df, aes(average_intensity, logFC, color = direction)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey40", linewidth = 0.3) +
  geom_point(data = \(d) filter(d, direction == "NS"),
             alpha = 0.25, size = 0.6) +
  geom_point(data = \(d) filter(d, direction != "NS"),
             alpha = 0.85, size = 0.9) +
  geom_text(data = n_dep,
            aes(x = Inf, y = Inf,
                label = sprintf("Pi DEPs: %d", n_dep)),
            inherit.aes = FALSE,
            hjust = 1.05, vjust = 1.5,
            size = 3.0, fontface = "bold", color = "grey20") +
  scale_color_manual(values = DIR_COLORS, name = "Pi < 0.05") +
  facet_wrap(~ contrast, ncol = 2,
             labeller = labeller(contrast = CTR_SHORT)) +
  labs(title = "MA plots by contrast",
       subtitle = "logFC vs mean log2 intensity | red up, blue down (Pi-score DEPs)",
       x = "Mean log2 intensity",
       y = "log2 fold-change",
       tag = "B") +
  FIG_THEME +
  theme(strip.background = element_blank(),
        strip.text = element_text(face = "bold", size = FIG_STRIP_SIZE),
        legend.position = "top",
        legend.key.size = unit(3, "mm"))

PW <- 160; PH <- 120
ggsave(file.path(RPT_DIR, "SUPP_panel_C_ma.png"), pC,
       width = PW, height = PH, units = "mm", dpi = 300)

message("F03 SUPP panel C (MA) saved")
