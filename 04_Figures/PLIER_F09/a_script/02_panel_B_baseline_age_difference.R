# Panel B — Baseline age differences in pathway activity.
# Q1: What pathways distinguish Young from Old at baseline?
# Method: Welch t-test per LV (Pre-only samples), BH correction.

setwd(rprojroot::find_rstudio_root_file())
source("04_Figures/shared/style.R")

library(readr)
library(dplyr)
library(tidyr)
library(ggplot2)

DAT <- "04_Figures/PLIER_F09/c_data"
RPT <- "04_Figures/PLIER_F09/b_results"

# --- Load data ---
lv_scores <- read_csv(file.path(DAT, "02_lv_scores.csv"), show_col_types = FALSE)
dal <- readRDS("02_Imputation/c_data/01_DAList_imputed.rds")
meta <- as.data.frame(dal$metadata) |>
  mutate(subject = sub("_(Pre|Post)$", "", Col_ID))

# Merge LV scores with metadata
lv_meta <- lv_scores |>
  inner_join(meta, by = c("sample_id" = "Col_ID"))

# Pre-only
pre <- lv_meta |> filter(Timepoint == "Pre")
lv_cols <- grep("^LV\\d+$", names(pre), value = TRUE)

cat(sprintf("Baseline comparison: Young n=%d, Old n=%d\n",
            sum(pre$Group == "Young"), sum(pre$Group == "Old")))

# --- Welch t-test per LV ---
results <- lapply(lv_cols, function(lv) {
  y <- pre[[lv]][pre$Group == "Young"]
  o <- pre[[lv]][pre$Group == "Old"]
  tt <- t.test(y, o)
  d <- (mean(y) - mean(o)) / sqrt((sd(y)^2 + sd(o)^2) / 2)  # Cohen's d
  tibble(LV = lv,
         mean_young = mean(y), mean_old = mean(o),
         diff = mean(y) - mean(o),
         cohens_d = d,
         pvalue = tt$p.value)
}) |> bind_rows() |>
  mutate(padj = p.adjust(pvalue, method = "BH"),
         stars = sig_stars(padj))

# Add pathway labels from alignment summary
align_df <- read_csv(file.path(DAT, "05_lv_pathway_alignment_summary.csv"),
                     show_col_types = FALSE)
top_label <- align_df |>
  filter(FDR < 0.05) |>
  group_by(LV) |>
  slice_max(AUC, n = 1, with_ties = FALSE) |>
  ungroup() |>
  mutate(lv_name = paste0("LV", LV),
         pathway_label = clean_pathway_name(pathway))

results <- results |>
  left_join(top_label |> select(lv_name, pathway_label),
            by = c("LV" = "lv_name")) |>
  mutate(label = ifelse(is.na(pathway_label), LV, pathway_label)) |>
  arrange(pvalue)

write_csv(results, file.path(DAT, "12_panel_B_baseline_age_difference.csv"))

# --- Figure: top significant LVs ---
sig_res <- results |> filter(padj < 0.10) |> slice_head(n = 15)
if (nrow(sig_res) == 0) sig_res <- results |> slice_head(n = 10)

plot_df <- sig_res |>
  mutate(label = factor(label, levels = rev(label)))

pB <- ggplot(plot_df, aes(x = cohens_d, y = label, fill = cohens_d > 0)) +
  geom_col(width = 0.7, color = "black", linewidth = 0.3) +
  geom_text(aes(label = stars, x = cohens_d + sign(cohens_d) * 0.05),
            hjust = ifelse(plot_df$cohens_d > 0, 0, 1),
            size = 3, fontface = "bold") +
  scale_fill_manual(values = c("TRUE" = unname(AGE_COLORS["Young"]),
                                "FALSE" = unname(AGE_COLORS["Old"])),
                    labels = c("TRUE" = "Higher in Young",
                               "FALSE" = "Higher in Old"),
                    name = NULL) +
  geom_vline(xintercept = 0, linewidth = 0.4) +
  labs(title = "B  Baseline Age Differences in Pathway Activity",
       subtitle = sprintf("Welch t-test, BH-adjusted (Pre-only, n=%d)", nrow(pre)),
       x = "Cohen's d (Young - Old)", y = NULL) +
  FIG_THEME +
  theme(legend.position = c(0.85, 0.15),
        legend.background = element_rect(fill = alpha("white", 0.85), color = NA))

ggsave(file.path(RPT, "12_panel_B_baseline_age_difference.pdf"), pB,
       width = 180, height = 140, units = "mm")
ggsave(file.path(RPT, "12_panel_B_baseline_age_difference.png"), pB,
       width = 180, height = 140, units = "mm", dpi = 300)

n_sig <- sum(results$padj < 0.05)
cat(sprintf("Panel B: %d/%d LVs with padj < 0.05\n", n_sig, nrow(results)))
if (n_sig > 0) {
  top <- results |> filter(padj < 0.05) |> slice_head(n = 3)
  cat("  Top 3:\n")
  for (i in seq_len(nrow(top))) {
    cat(sprintf("    %s: d=%.2f, padj=%.3g\n", top$label[i], top$cohens_d[i], top$padj[i]))
  }
}
