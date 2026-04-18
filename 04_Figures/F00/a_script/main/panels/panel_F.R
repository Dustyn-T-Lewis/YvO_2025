# F00 Panel F: Observed vs imputed density
setwd(rprojroot::find_rstudio_root_file())
source("04_Figures/shared/style.R")

suppressPackageStartupMessages({
  library(dplyr)
  library(tibble)
  library(ggplot2)
})

PW <- 100; PH <- 100
RPT_PNG <- "04_Figures/F00/b_reports/main/png/panels"
RPT_PDF <- "04_Figures/F00/b_reports/main/pdf/panels"
DAT     <- "04_Figures/F00/c_data"
dir.create(RPT_PNG, recursive = TRUE, showWarnings = FALSE)
dir.create(RPT_PDF, recursive = TRUE, showWarnings = FALSE)
dir.create(DAT, recursive = TRUE, showWarnings = FALSE)

if (!exists("int_imp"))
  int_imp <- readRDS("02_Imputation/c_data/00_report_intermediates.rds")

obs_vals <- as.numeric(int_imp$mat[!int_imp$was_na])
imp_vals <- as.numeric(int_imp$mat_imp[int_imp$was_na])

dens_df <- bind_rows(
  tibble(value = obs_vals, type = "Observed"),
  tibble(value = imp_vals, type = "Imputed (missForest)")
) %>%
  mutate(type = factor(type, levels = c("Observed", "Imputed (missForest)")))

oob <- int_imp$oob_error
oob_lab <- if (!is.null(oob) && !is.na(oob)) sprintf("OOB = %.3f", oob) else "OOB = NA"

# Summary stats for supplementary (full vectors are too large for Excel)
summary_df <- dens_df %>%
  group_by(type) %>%
  summarise(
    n = dplyr::n(),
    mean = mean(value, na.rm = TRUE),
    median = median(value, na.rm = TRUE),
    sd = sd(value, na.rm = TRUE),
    q25 = quantile(value, 0.25, na.rm = TRUE),
    q75 = quantile(value, 0.75, na.rm = TRUE),
    min = min(value, na.rm = TRUE),
    max = max(value, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(oob_error = ifelse(type == "Imputed (missForest)", oob, NA_real_))

write.csv(summary_df, file.path(DAT, "panel_F_imputation_summary.csv"), row.names = FALSE)

pF <- ggplot(dens_df, aes(value, fill = type, color = type)) +
  geom_density(alpha = 0.4, linewidth = 0.4) +
  scale_fill_manual(values = c(Observed = "#377EB8",
                               "Imputed (missForest)" = "#E41A1C"), name = NULL) +
  scale_color_manual(values = c(Observed = "#377EB8",
                                "Imputed (missForest)" = "#E41A1C"), name = NULL) +
  labs(title = "Imputation effect",
       subtitle = sprintf("%s observed | %s imputed | %s",
                          format(length(obs_vals), big.mark = ","),
                          format(length(imp_vals), big.mark = ","),
                          oob_lab),
       x = "log2 intensity", y = "Density", tag = "F") +
  FIG_THEME + theme(legend.position = "top",
                    legend.key.size = unit(3, "mm"),
                    legend.text = element_text(size = 7))

ggsave(file.path(RPT_PNG, "MAIN_panel_F_imputation_density.png"), pF,
       width = PW, height = PH, units = "mm", dpi = 300)
ggsave(file.path(RPT_PDF, "MAIN_panel_F_imputation_density.pdf"), pF,
       width = PW, height = PH, units = "mm", device = get_pdf_device())
cat("F00 Panel F done\n")
