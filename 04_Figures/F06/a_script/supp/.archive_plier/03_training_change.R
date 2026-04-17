# Supplementary: PLIER — Training-induced pathway changes within each age group.
# Q4: Which pathways change with training in Young? In Old?
# Method: Paired t-test per LV within each group, BH correction.

setwd(rprojroot::find_rstudio_root_file())
source("04_Figures/shared/style.R")

library(readr)
library(dplyr)
library(tidyr)
library(ggplot2)

DAT <- "04_Figures/F06/c_data/supp/plier"
RPT <- "04_Figures/F06/b_reports/supp/plier"
dir.create(DAT, recursive = TRUE, showWarnings = FALSE)
dir.create(RPT, recursive = TRUE, showWarnings = FALSE)

# --- Load data
lv_scores <- read_csv(file.path(DAT, "02_lv_scores.csv"), show_col_types = FALSE)
dal <- readRDS("02_Imputation/c_data/01_DAList_imputed.rds")
meta <- as.data.frame(dal$metadata) |>
  mutate(subject = sub("_(Pre|Post)$", "", Col_ID))

lv_meta <- lv_scores |>
  inner_join(meta, by = c("sample_id" = "Col_ID"))
lv_cols <- grep("^LV\\d+$", names(lv_meta), value = TRUE)

# --- Compute paired delta per subject
wide <- lv_meta |>
  select(subject, Group, Timepoint, all_of(lv_cols)) |>
  pivot_wider(names_from = Timepoint, values_from = all_of(lv_cols),
              names_sep = ".")

# --- Paired t-tests within each group
run_paired <- function(grp) {
  sub_df <- wide |> filter(Group == grp)
  n <- nrow(sub_df)

  lapply(lv_cols, function(lv) {
    pre_col <- paste0(lv, ".Pre")
    post_col <- paste0(lv, ".Post")
    if (!all(c(pre_col, post_col) %in% names(sub_df))) return(NULL)

    pre_vals <- sub_df[[pre_col]]
    post_vals <- sub_df[[post_col]]
    delta <- post_vals - pre_vals

    tt <- t.test(delta, mu = 0)
    tibble(
      LV = lv, group = grp, n = n,
      mean_pre = mean(pre_vals), mean_post = mean(post_vals),
      mean_delta = mean(delta), sd_delta = sd(delta),
      t_stat = tt$statistic, pvalue = tt$p.value
    )
  }) |> bind_rows()
}

results <- bind_rows(run_paired("Young"), run_paired("Old")) |>
  group_by(group) |>
  mutate(padj = p.adjust(pvalue, method = "BH")) |>
  ungroup() |>
  mutate(stars = sig_stars(padj))

# Add pathway labels
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
  mutate(label = ifelse(is.na(pathway_label), LV, pathway_label))

write_csv(results, file.path(DAT, "15_panel_E_training_change_within_group.csv"))

# --- Figure: side-by-side training effects
# Show LVs significant in at least one group
sig_lvs <- results |>
  filter(padj < 0.10) |>
  pull(LV) |>
  unique()

if (length(sig_lvs) == 0) sig_lvs <- results |>
  arrange(pvalue) |> slice_head(n = 10) |> pull(LV) |> unique()

plot_df <- results |>
  filter(LV %in% sig_lvs) |>
  mutate(group = factor(group, levels = c("Young", "Old")))

# Order by mean delta across groups
lv_order <- plot_df |>
  group_by(label) |>
  summarise(avg_delta = mean(mean_delta)) |>
  arrange(avg_delta) |>
  pull(label)
plot_df$label <- factor(plot_df$label, levels = lv_order)

pE <- ggplot(plot_df, aes(x = mean_delta, y = label, fill = group)) +
  geom_col(position = position_dodge(width = 0.7), width = 0.6,
           color = "black", linewidth = 0.3) +
  geom_errorbarh(aes(xmin = mean_delta - sd_delta / sqrt(n),
                     xmax = mean_delta + sd_delta / sqrt(n)),
                 position = position_dodge(width = 0.7),
                 height = 0.2, linewidth = 0.3) +
  geom_text(aes(label = stars,
                x = mean_delta + sign(mean_delta) * (sd_delta / sqrt(n) + 0.02)),
            position = position_dodge(width = 0.7),
            hjust = ifelse(plot_df$mean_delta > 0, 0, 1),
            size = 2.5, fontface = "bold", show.legend = FALSE) +
  geom_vline(xintercept = 0, linewidth = 0.4) +
  scale_fill_manual(values = AGE_COLORS, name = NULL) +
  labs(title = "E  Training-Induced Pathway Changes",
       subtitle = "Paired t-test (Post - Pre) within each age group, BH-adjusted",
       x = "Mean \u0394 LV Score (Post - Pre)", y = NULL) +
  FIG_THEME +
  theme(legend.position = "bottom")

ggsave(file.path(RPT, "15_panel_E_training_change_within_group.pdf"), pE,
       width = 200, height = 180, units = "mm")
ggsave(file.path(RPT, "15_panel_E_training_change_within_group.png"), pE,
       width = 200, height = 180, units = "mm", dpi = 300)

for (grp in c("Young", "Old")) {
  n_sig <- sum(results$group == grp & results$padj < 0.05)
  cat(sprintf("Panel E (%s): %d LVs with padj < 0.05\n", grp, n_sig))
}
