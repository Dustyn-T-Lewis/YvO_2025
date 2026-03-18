# Panel G — Age moderation of training response.
# Q6: Does age group moderate the pathway-level training response?
# Method: Welch t-test comparing delta LV (Young vs Old), BH correction.

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

lv_meta <- lv_scores |>
  inner_join(meta, by = c("sample_id" = "Col_ID"))
lv_cols <- grep("^LV\\d+$", names(lv_meta), value = TRUE)

# --- Compute subject-level delta LV ---
# Keep only subjects with both Pre and Post
paired_subjects <- lv_meta |>
  count(subject) |>
  filter(n == 2) |>
  pull(subject)

lv_wide <- lv_meta |>
  filter(subject %in% paired_subjects) |>
  select(subject, Group, Timepoint, all_of(lv_cols)) |>
  pivot_wider(names_from = Timepoint, values_from = all_of(lv_cols),
              names_sep = ".")

delta_df <- data.frame(subject = lv_wide$subject, Group = lv_wide$Group)
for (lv in lv_cols) {
  delta_df[[lv]] <- lv_wide[[paste0(lv, ".Post")]] - lv_wide[[paste0(lv, ".Pre")]]
}

n_young <- sum(delta_df$Group == "Young")
n_old <- sum(delta_df$Group == "Old")
cat(sprintf("Age moderation: Young n=%d, Old n=%d\n", n_young, n_old))

# --- Welch t-test: delta Young vs delta Old ---
results <- lapply(lv_cols, function(lv) {
  y <- delta_df[[lv]][delta_df$Group == "Young"]
  o <- delta_df[[lv]][delta_df$Group == "Old"]
  tt <- t.test(y, o)
  d <- (mean(y) - mean(o)) / sqrt((sd(y)^2 + sd(o)^2) / 2)
  tibble(
    LV = lv,
    mean_delta_young = mean(y), mean_delta_old = mean(o),
    diff = mean(y) - mean(o),
    cohens_d = d,
    pvalue = tt$p.value
  )
}) |> bind_rows() |>
  mutate(padj = p.adjust(pvalue, method = "BH"),
         stars = sig_stars(padj))

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
  mutate(label = ifelse(is.na(pathway_label), LV, pathway_label)) |>
  arrange(pvalue)

# Count pathway-aligned LVs with blunted response
aligned_lvs <- results |> filter(!is.na(pathway_label))
blunted <- aligned_lvs |>
  filter(padj < 0.05, abs(mean_delta_young) > abs(mean_delta_old))

write_csv(results, file.path(DAT, "17_panel_G_age_moderation.csv"))

# --- Figure ---
sig_res <- results |> filter(padj < 0.10) |> slice_head(n = 15)
if (nrow(sig_res) == 0) sig_res <- results |> slice_head(n = 10)

plot_df <- sig_res |>
  mutate(label = factor(label, levels = rev(label)))

pG <- ggplot(plot_df, aes(x = cohens_d, y = label, fill = cohens_d > 0)) +
  geom_col(width = 0.7, color = "black", linewidth = 0.3) +
  geom_text(aes(label = stars,
                x = cohens_d + sign(cohens_d) * 0.05),
            hjust = ifelse(plot_df$cohens_d > 0, 0, 1),
            size = 2.5, fontface = "bold", show.legend = FALSE) +
  geom_vline(xintercept = 0, linewidth = 0.4) +
  scale_fill_manual(values = c("TRUE" = unname(AGE_COLORS["Young"]),
                                "FALSE" = unname(AGE_COLORS["Old"])),
                    labels = c("TRUE" = "Larger delta in Young",
                               "FALSE" = "Larger delta in Old"),
                    name = NULL) +
  labs(title = "G  Age Moderation of Training Response",
       subtitle = sprintf("Welch t-test on \u0394 LV (Young vs Old), BH-adjusted (n=%d, %d)",
                           n_young, n_old),
       x = "Cohen's d (\u0394Young - \u0394Old)", y = NULL) +
  FIG_THEME +
  theme(legend.position = c(0.85, 0.15),
        legend.background = element_rect(fill = alpha("white", 0.85), color = NA))

ggsave(file.path(RPT, "17_panel_G_age_moderation.pdf"), pG,
       width = 180, height = 150, units = "mm")
ggsave(file.path(RPT, "17_panel_G_age_moderation.png"), pG,
       width = 180, height = 150, units = "mm", dpi = 300)

n_sig <- sum(results$padj < 0.05)
cat(sprintf("Panel G: %d/%d LVs with padj < 0.05\n", n_sig, nrow(results)))
cat(sprintf("  Blunted in Old (pathway-aligned): %d/%d\n",
            nrow(blunted), nrow(aligned_lvs)))
