# F00_QC Stage 02 — Imputation detail composite
# 2x2 grid, 200x160 mm. Goes deeper than the MAIN_F00_QC overview on the
# missForest imputation stage. Self-contained; reads the stage 02 report
# intermediates RDS and the 23-method benchmark ranking CSV.
# Panels:
#   A Per-protein %missing x MAR/MNAR   B Benchmark top-5 methods
#   C MNAR audit shift distribution     D Sample-level mean intensity pre vs post
# Output: 04_Figures/F00/b_reports/MAIN_F00_QC_stage02_imputation.{pdf,png}

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

int02 <- readRDS("02_Imputation/c_data/00_report_intermediates.rds")
bench <- read_csv("02_Imputation/c_data/benchmark/04_composite_ranking.csv",
                   show_col_types = FALSE)

# --- Panel A: %missing x classification scatter ---
mc <- int02$miss_class %>%
  mutate(classification = factor(classification,
                                 levels = c("Complete", "MAR", "MNAR")))

n_total <- nrow(mc)
class_counts <- mc %>%
  count(classification) %>%
  mutate(lab = sprintf("%s: %d (%.1f%%)", classification, n, 100*n/n_total))

pA <- ggplot(mc, aes(mean_intensity, pct_miss, color = classification)) +
  geom_point(alpha = 0.45, size = 0.9) +
  scale_color_manual(values = PAL_CLASS, name = NULL,
                     labels = class_counts$lab) +
  scale_y_continuous(labels = function(x) paste0(x, "%")) +
  labs(title = "Per-protein missingness x class",
       subtitle = sprintf("%s proteins | k-means on mean intensity",
                          format(n_total, big.mark=",")),
       x = "Mean log2 intensity", y = "% missing", tag = "A") +
  FIG_THEME +
  theme(legend.position = "top",
        legend.key.size = unit(3, "mm"),
        legend.text = element_text(size = 7))

# --- Panel B: benchmark top-5 methods ---
topN <- bench %>%
  arrange(rank) %>%
  slice_head(n = 5) %>%
  mutate(method = factor(method, levels = rev(method)),
         is_winner = method == levels(method)[length(levels(method))])

pB <- ggplot(topN, aes(composite, method,
                        fill = ifelse(is_winner, "Winner", "Top 5"))) +
  geom_col(width = 0.7) +
  geom_text(aes(label = sprintf("%.3f", composite)),
            hjust = -0.1, size = 2.7, fontface = "bold") +
  scale_fill_manual(values = c(Winner = "#E41A1C", "Top 5" = "#377EB8"),
                    name = NULL) +
  scale_x_continuous(expand = expansion(mult = c(0, 0.15)),
                     limits = c(0, max(topN$composite) * 1.15)) +
  labs(title = "23-method benchmark (top 5)",
       subtitle = sprintf("Composite score | winner: %s (%.3f)",
                          int02$best %||% "missForest",
                          max(topN$composite)),
       x = "Composite score", y = NULL, tag = "B") +
  FIG_THEME +
  theme(legend.position = "top",
        legend.key.size = unit(3, "mm"),
        legend.text = element_text(size = 7),
        axis.text.y = element_text(size = 7))

# --- Panel C: MNAR audit shift distribution ---
mn_audit <- int02$mnar_audit
pC <- ggplot(mn_audit, aes(shift)) +
  geom_histogram(binwidth = 0.025, fill = "#E41A1C", alpha = 0.70,
                 color = "white", linewidth = 0.2) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey30",
             linewidth = 0.4) +
  geom_vline(xintercept = mean(mn_audit$shift, na.rm = TRUE),
             linetype = "solid", color = "#B2182B", linewidth = 0.4) +
  annotate("text", x = mean(mn_audit$shift, na.rm=TRUE), y = Inf,
           label = sprintf("mean = %+.3f", mean(mn_audit$shift, na.rm=TRUE)),
           vjust = 1.4, hjust = -0.1, size = 2.7, color = "#B2182B") +
  labs(title = "MNAR audit: imputation shift",
       subtitle = sprintf("%d MNAR proteins | pre-mean vs post-mean (log2)",
                          nrow(mn_audit)),
       x = "Shift (log2 intensity)", y = "MNAR proteins", tag = "C") +
  FIG_THEME

# --- Panel D: sample-level mean intensity pre vs post imputation ---
# Pre = mean of observed values only; post = mean of full imputed matrix.
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

pD <- ggplot(md, aes(pre, post, color = Group_Time)) +
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

# --- Composite ---
composite <- (pA | pB) / (pC | pD) +
  plot_annotation(
    title = "YvO Proteomics Pipeline - Stage 02: Imputation detail",
    subtitle = sprintf(
      "Missingness x class | method benchmark | MNAR audit | sample integrity | %s",
      format(Sys.Date(), "%Y-%m-%d")),
    theme = theme(plot.title    = element_text(face = "bold", size = 13),
                  plot.subtitle = element_text(face = "italic", size = 9,
                                               color = "grey30"))
  ) &
  theme(plot.tag = element_text(face = "bold", size = 14))

COMP_W <- 200
COMP_H <- 160

ggsave(file.path(RPT_DIR, "MAIN_F00_QC_stage02_imputation.pdf"),
       composite, width = COMP_W, height = COMP_H, units = "mm",
       device = pdf_device)
ggsave(file.path(RPT_DIR, "MAIN_F00_QC_stage02_imputation.png"),
       composite, width = COMP_W, height = COMP_H, units = "mm", dpi = 300)

message(sprintf("Stage 02 composite saved -> %s", RPT_DIR))
