# F00 Stage 03 Panel D: Training blunting scatter
setwd(rprojroot::find_rstudio_root_file())
source("04_Figures/shared/style.R")

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(ggplot2)
})

PW <- 100; PH <- 80
RPT_PNG <- "04_Figures/F00/b_reports/03_DEP/png/panels"
RPT_PDF <- "04_Figures/F00/b_reports/03_DEP/pdf/panels"
DAT     <- "04_Figures/F00/c_data/03_DEP"
dir.create(RPT_PNG, recursive = TRUE, showWarnings = FALSE)
dir.create(RPT_PDF, recursive = TRUE, showWarnings = FALSE)
dir.create(DAT, recursive = TRUE, showWarnings = FALSE)

if (!exists("combo"))
  combo <- read_csv("03_DEP/c_data/03_combined_results.csv", show_col_types = FALSE)
if (!exists("blunt"))
  blunt <- read_csv("03_DEP/c_data/06_blunting_diagnostics.csv", show_col_types = FALSE)

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

write.csv(blunt_df, file.path(DAT, "panel_D_blunting.csv"), row.names = FALSE)

pD_s03 <- ggplot(blunt_df, aes(abs_TY, abs_TO, color = panel)) +
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

ggsave(file.path(RPT_PNG, "panel_D_blunting.png"), pD_s03,
       width = PW, height = PH, units = "mm", dpi = 300)
ggsave(file.path(RPT_PDF, "panel_D_blunting.pdf"), pD_s03,
       width = PW, height = PH, units = "mm", device = get_pdf_device())
cat("F00 Stage03 Panel D done\n")
