# F00 Stage 02 Panel D: Sample-level mean intensity pre vs post imputation
setwd(rprojroot::find_rstudio_root_file())
source("04_Figures/shared/style.R")

suppressPackageStartupMessages({
  library(dplyr)
  library(tibble)
  library(ggplot2)
})

PW <- 100; PH <- 80
RPT_PNG <- "04_Figures/F00/b_reports/02_imputation/png/panels"
RPT_PDF <- "04_Figures/F00/b_reports/02_imputation/pdf/panels"
DAT     <- "04_Figures/F00/c_data/02_imputation"
dir.create(RPT_PNG, recursive = TRUE, showWarnings = FALSE)
dir.create(RPT_PDF, recursive = TRUE, showWarnings = FALSE)
dir.create(DAT, recursive = TRUE, showWarnings = FALSE)

if (!exists("int02"))
  int02 <- readRDS("02_Imputation/c_data/00_report_intermediates.rds")

mat     <- int02$mat
mat_imp <- int02$mat_imp
pre_mean  <- colMeans(mat,     na.rm = TRUE)
post_mean <- colMeans(mat_imp, na.rm = FALSE)

md <- tibble(Col_ID = colnames(mat), pre = pre_mean, post = post_mean)
if (!is.null(int02$meta) && "Col_ID" %in% colnames(int02$meta)) {
  md <- md %>% left_join(int02$meta %>% select(Col_ID, Group_Time) %>% distinct(),
                         by = "Col_ID")
} else {
  md$Group_Time <- "Sample"
}

r2 <- cor(md$pre, md$post)^2

write.csv(md, file.path(DAT, "panel_D_sample_intensity.csv"), row.names = FALSE)

pD_s02 <- ggplot(md, aes(pre, post, color = Group_Time)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed",
              color = "grey60", linewidth = 0.4) +
  geom_point(size = 2.2, alpha = 0.85) +
  scale_color_manual(values = GROUP_FILL, name = NULL,
                     na.value = "grey50") +
  labs(title = "Mean intensity: observed vs imputed",
       subtitle = sprintf("R^2 = %.3f | OOB = %.3f",
                          r2, int02$oob_error %||% NA_real_),
       x = "Per-sample mean (observed)",
       y = "Per-sample mean (imputed)", tag = "D") +
  FIG_THEME +
  theme(legend.position = "top",
        legend.key.size = unit(3, "mm"),
        legend.text = element_text(size = 7))

ggsave(file.path(RPT_PNG, "panel_D_sample_intensity.png"), pD_s02,
       width = PW, height = PH, units = "mm", dpi = 300)
ggsave(file.path(RPT_PDF, "panel_D_sample_intensity.pdf"), pD_s02,
       width = PW, height = PH, units = "mm", device = get_pdf_device())
cat("F00 Stage02 Panel D done\n")
