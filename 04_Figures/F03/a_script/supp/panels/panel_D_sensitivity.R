# F03 SUPP — Panel D: DEP rank correlation under imputation perturbation
# Horizontal bars of Spearman rho(full, reduced-imputation) per contrast.
# Rho near 1 = DEP ranking is insensitive to imputation choices.
#
# Reads:  03_DEP/c_data/09_imputation_sensitivity.csv
# Writes: b_reports/supp/{png,pdf}/panels/SUPP_panel_D_sensitivity.{png,pdf}
#         c_data/supp/panel_D_sensitivity.csv

setwd(rprojroot::find_rstudio_root_file())
source("04_Figures/shared/style.R")

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(ggplot2)
})

RPT_PNG <- "04_Figures/F03/b_reports/supp/png/panels"
RPT_PDF <- "04_Figures/F03/b_reports/supp/pdf/panels"
DAT_DIR <- "04_Figures/F03/c_data/supp"
dir.create(RPT_PNG, recursive = TRUE, showWarnings = FALSE)
dir.create(RPT_PDF, recursive = TRUE, showWarnings = FALSE)
dir.create(DAT_DIR, recursive = TRUE, showWarnings = FALSE)

stopifnot(
  "03_DEP imputation sensitivity results missing" =
    file.exists("03_DEP/c_data/09_imputation_sensitivity.csv")
)

sens <- read_csv("03_DEP/c_data/09_imputation_sensitivity.csv",
                 show_col_types = FALSE) %>%
  mutate(contrast = factor(contrast,
                           levels = c("Aging", "Training_Young",
                                      "Training_Old", "Interaction")))

write_csv(sens, file.path(DAT_DIR, "panel_D_sensitivity.csv"))

# Order contrasts top to bottom by rho descending so the most-robust is on top
sens <- sens %>% arrange(desc(spearman_rho))
sens$contrast_f <- factor(sens$contrast, levels = rev(as.character(sens$contrast)))

pD <- ggplot(sens, aes(spearman_rho, contrast_f, fill = contrast)) +
  geom_col(width = 0.6) +
  geom_text(aes(label = sprintf("rho = %.4f (n = %s)",
                                spearman_rho,
                                format(n_proteins, big.mark = ","))),
            hjust = 1.05, color = "white", fontface = "bold", size = 3.2) +
  scale_fill_manual(values = CONTRAST_COLORS, guide = "none") +
  scale_x_continuous(limits = c(0, 1.0), breaks = seq(0, 1, 0.1),
                     expand = expansion(mult = c(0, 0.02))) +
  scale_y_discrete(labels = CTR_SHORT) +
  labs(title = "DEP rank stability under imputation perturbation",
       subtitle = "Spearman rho(all-methods sensitivity analysis) per contrast | 1.0 = perfect stability",
       x = "Spearman rho(full vs reduced imputation)",
       y = NULL,
       tag = "c") +
  FIG_THEME

PW <- 89; PH <- 75
ggsave(file.path(RPT_PNG, "SUPP_panel_D_sensitivity.png"), pD,
       width = PW, height = PH, units = "mm", dpi = 300)
ggsave(file.path(RPT_PDF, "SUPP_panel_D_sensitivity.pdf"), pD,
       width = PW, height = PH, units = "mm", device = get_pdf_device())

message("F03 SUPP panel D (imputation sensitivity) saved")

# --- Export for composite ---
pD_title    <- "DEP rank stability under imputation perturbation"
pD_subtitle <- "Spearman rho(all-methods sensitivity analysis) per contrast | 1.0 = perfect stability"
pD_legend   <- NULL
pD          <- strip_for_composite(pD)
