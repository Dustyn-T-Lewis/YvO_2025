# Supplementary: PLIER — PLIER <-> WGCNA convergence.
# Correlate PLIER LV scores with WGCNA module eigengenes.

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

# --- Load PLIER LV scores
lv_scores <- read_csv(file.path(DAT, "02_lv_scores.csv"), show_col_types = FALSE)
lv_cols <- grep("^LV\\d+$", names(lv_scores), value = TRUE)

# --- Load WGCNA eigengenes (saved by YvO_WGCNA_run.R)
wgcna_dat <- "04_Figures/F06/c_data"
stopifnot(
  "WGCNA MEs.rds not found — run F06 first" =
    file.exists(file.path(wgcna_dat, "MEs.rds")),
  "WGCNA datExpr.rds not found — run F06 first" =
    file.exists(file.path(wgcna_dat, "datExpr.rds"))
)
MEs <- readRDS(file.path(wgcna_dat, "MEs.rds"))
datExpr <- readRDS(file.path(wgcna_dat, "datExpr.rds"))

# --- Match samples
me_df <- MEs |>
  as.data.frame() |>
  mutate(sample_id = rownames(datExpr))

merged <- lv_scores |>
  inner_join(me_df, by = "sample_id")

me_cols <- grep("^ME", names(merged), value = TRUE)
# Drop grey module
me_cols <- me_cols[me_cols != "MEgrey"]

cat(sprintf("Convergence analysis: %d samples, %d LVs, %d modules\n",
            nrow(merged), length(lv_cols), length(me_cols)))

# --- Pearson correlation matrix
lv_mat <- as.matrix(merged[, lv_cols])
me_mat <- as.matrix(merged[, me_cols])

cor_mat <- cor(lv_mat, me_mat, use = "pairwise.complete.obs")
n_obs <- nrow(merged)

# P-values
pval_mat <- matrix(NA, nrow = length(lv_cols), ncol = length(me_cols),
                   dimnames = list(lv_cols, me_cols))
for (i in seq_along(lv_cols)) {
  for (j in seq_along(me_cols)) {
    ct <- cor.test(lv_mat[, i], me_mat[, j])
    pval_mat[i, j] <- ct$p.value
  }
}

# --- Export full correlation matrix
full_cor <- as.data.frame(cor_mat) |>
  mutate(LV = rownames(cor_mat)) |>
  select(LV, everything())
write_csv(full_cor, file.path(DAT, "18_panel_H_wgcna_convergence_full.csv"))

# --- Filtered convergence (|r| > 0.3)
filt_long <- expand.grid(LV = lv_cols, ME = me_cols, stringsAsFactors = FALSE) |>
  mutate(r = as.vector(cor_mat),
         pvalue = as.vector(pval_mat)) |>
  filter(abs(r) > 0.3) |>
  mutate(padj = p.adjust(pvalue, method = "BH")) |>
  arrange(desc(abs(r)))

write_csv(filt_long, file.path(DAT, "18_panel_H_wgcna_convergence_filtered.csv"))

# --- Heatmap
# Add pathway labels to LVs
align_df <- read_csv(file.path(DAT, "05_lv_pathway_alignment_summary.csv"),
                     show_col_types = FALSE)
top_label <- align_df |>
  filter(FDR < 0.05) |>
  group_by(LV) |>
  slice_max(AUC, n = 1, with_ties = FALSE) |>
  ungroup() |>
  mutate(lv_name = paste0("LV", LV),
         pathway_label = clean_pathway_name(pathway, max_chars = 30))

# Only show LVs with at least one |r| > 0.3
lvs_to_show <- unique(filt_long$LV)
if (length(lvs_to_show) == 0) lvs_to_show <- lv_cols[1:min(10, length(lv_cols))]

heatmap_df <- expand.grid(LV = lvs_to_show, ME = me_cols, stringsAsFactors = FALSE) |>
  mutate(r = mapply(function(l, m) cor_mat[l, m], LV, ME)) |>
  left_join(top_label |> select(lv_name, pathway_label),
            by = c("LV" = "lv_name")) |>
  mutate(lv_label = ifelse(is.na(pathway_label), LV,
                           paste0(LV, ": ", pathway_label)),
         me_label = gsub("^ME", "", ME))

pH <- ggplot(heatmap_df, aes(x = me_label, y = lv_label, fill = r)) +
  geom_tile(color = "white", linewidth = 0.3) +
  geom_text(aes(label = ifelse(abs(r) > 0.3, sprintf("%.2f", r), "")),
            size = 2, color = "black") +
  scale_fill_gradient2(low = "#4393C3", mid = "white", high = "#D6604D",
                       midpoint = 0, limits = c(-1, 1), name = "r") +
  labs(title = "H  PLIER-WGCNA Convergence",
       subtitle = sprintf("Pearson r, n=%d samples (supplementary)", n_obs),
       x = "WGCNA Module", y = "PLIER Latent Variable") +
  FIG_THEME +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 7),
        axis.text.y = element_text(size = 6))

ggsave(file.path(RPT, "18_panel_H_wgcna_convergence.pdf"), pH,
       width = 220, height = 200, units = "mm")
ggsave(file.path(RPT, "18_panel_H_wgcna_convergence.png"), pH,
       width = 220, height = 200, units = "mm", dpi = 300)

cat(sprintf("Panel H: %d LV-module pairs with |r| > 0.3\n", nrow(filt_long)))
if (nrow(filt_long) > 0) {
  cat("  Strongest:\n")
  top3 <- filt_long |> slice_head(n = 3)
  for (i in seq_len(nrow(top3))) {
    cat(sprintf("    %s <-> %s: r=%.3f, padj=%.3g\n",
                top3$LV[i], top3$ME[i], top3$r[i], top3$padj[i]))
  }
}
