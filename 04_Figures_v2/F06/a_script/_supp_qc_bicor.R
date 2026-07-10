# Sourced by 02_supp_panels.R — expects style.R already loaded.
# Bicor sensitivity analysis: compares Pearson (main) vs biweight midcorrelation module overlap.

library(WGCNA)
library(tidyverse)

disableWGCNAThreads() # single-threaded so correlations/TOM are byte-reproducible across re-runs
set.seed(42)

BASE <- here::here("04_Figures_v2", "F06")

DATA_FILE <- here::here("02_imputation", "c_data", "01_imputed.csv")
PEARSON_MODS <- file.path(BASE, "c_data", "wgcna", "wgcna_module_assignments.csv")
RPT_PNG <- file.path(BASE, "b_reports", "supp", "png", "panels")
RPT_PDF <- file.path(BASE, "b_reports", "supp", "pdf", "panels")
DAT_OUT <- file.path(BASE, "c_data", "supp")

dir.create(RPT_PNG, recursive = TRUE, showWarnings = FALSE)
dir.create(RPT_PDF, recursive = TRUE, showWarnings = FALSE)
dir.create(DAT_OUT, recursive = TRUE, showWarnings = FALSE)

stopifnot(file.exists(DATA_FILE), file.exists(PEARSON_MODS))

df <- read_csv(DATA_FILE)
ann_cols <- c("uniprot_id", "protein", "gene", "description")
samp_names <- setdiff(names(df), ann_cols)
mat <- as.matrix(df[, samp_names])
rownames(mat) <- df$uniprot_id

datExpr <- t(mat)

gsg <- goodSamplesGenes(datExpr, verbose = 0)
if (!gsg$allOK) datExpr <- datExpr[gsg$goodSamples, gsg$goodGenes]

cor <- WGCNA::cor

message(sprintf("Bicor sensitivity: %d samples x %d proteins", nrow(datExpr), ncol(datExpr)))

powers <- seq_len(20)
sft_bicor <- pickSoftThreshold(datExpr,
  powerVector = powers,
  networkType = "signed",
  corFnc = "bicor",
  corOptions = list(maxPOutliers = 0.1),
  verbose = 2
)

r2_bicor <- -sign(sft_bicor$fitIndices$slope) * sft_bicor$fitIndices$SFT.R.sq
power_idx <- which(r2_bicor > 0.85)[1]
bicor_power <- if (!is.na(power_idx)) powers[power_idx] else 6

message(sprintf("  Bicor soft power: %d (R^2 = %.3f)", bicor_power, r2_bicor[bicor_power]))

net_bicor <- blockwiseModules(
  datExpr,
  power             = bicor_power,
  networkType       = "signed",
  TOMType           = "signed",
  corType           = "bicor",
  maxPOutliers      = 0.1,
  minModuleSize     = 30,
  mergeCutHeight    = 0.25,
  numericLabels     = TRUE,
  pamRespectsDendro = FALSE,
  saveTOMs          = FALSE,
  verbose           = 3
)

bicor_colors <- labels2colors(net_bicor$colors)
n_bicor <- length(unique(net_bicor$colors)) - (0 %in% net_bicor$colors)
message(sprintf("  Bicor modules: %d (+ grey)", n_bicor))

cor <- stats::cor # restore after WGCNA bicor/blockwise computations

pearson_df <- read_csv(PEARSON_MODS)

bicor_df <- tibble(
  uniprot_id   = colnames(datExpr),
  bicor_module = bicor_colors
)

compare <- pearson_df |>
  select(uniprot_id, pearson_module = module_color) |>
  inner_join(bicor_df, by = "uniprot_id")

pearson_mods <- setdiff(unique(compare$pearson_module), "grey")
bicor_mods <- setdiff(unique(compare$bicor_module), "grey")

jaccard_mat <- matrix(0,
  nrow = length(pearson_mods), ncol = length(bicor_mods),
  dimnames = list(pearson_mods, bicor_mods)
)

for (pm in pearson_mods) {
  p_set <- compare$uniprot_id[compare$pearson_module == pm]
  for (bm in bicor_mods) {
    b_set <- compare$uniprot_id[compare$bicor_module == bm]
    intersection <- length(intersect(p_set, b_set))
    union_size <- length(union(p_set, b_set))
    jaccard_mat[pm, bm] <- if (union_size > 0) intersection / union_size else 0
  }
}

jac_temp <- jaccard_mat
used_p <- character()
used_b <- character()
matched_rows <- list()

while (TRUE) {
  remaining <- jac_temp[!rownames(jac_temp) %in% used_p,
    !colnames(jac_temp) %in% used_b,
    drop = FALSE
  ]
  if (nrow(remaining) == 0 || ncol(remaining) == 0) break
  best_idx <- which(remaining == max(remaining), arr.ind = TRUE)[1, ]
  best_pm <- rownames(remaining)[best_idx[1]]
  best_bm <- colnames(remaining)[best_idx[2]]

  p_set <- compare$uniprot_id[compare$pearson_module == best_pm]
  b_set <- compare$uniprot_id[compare$bicor_module == best_bm]

  matched_rows <- c(matched_rows, list(tibble(
    pearson_module = best_pm,
    bicor_module   = best_bm,
    jaccard        = remaining[best_idx[1], best_idx[2]],
    n_shared       = length(intersect(p_set, b_set)),
    n_pearson      = length(p_set),
    n_bicor        = length(b_set)
  )))
  used_p <- c(used_p, best_pm)
  used_b <- c(used_b, best_bm)
}

matched <- bind_rows(matched_rows)

matched <- matched |> arrange(desc(jaccard))

write_csv(matched, file.path(DAT_OUT, "a02_bicor_sensitivity.csv"))
message(sprintf("\n  Matched modules: %d", nrow(matched)))
message(sprintf("  Mean Jaccard (matched): %.3f", mean(matched$jaccard)))
message(sprintf(
  "  Modules with Jaccard > 0.5: %d / %d",
  sum(matched$jaccard > 0.5), nrow(matched)
))

jac_long <- as_tibble(jaccard_mat, rownames = "pearson_module") |>
  pivot_longer(-pearson_module, names_to = "bicor_module", values_to = "jaccard")

pm_order <- c(
  matched$pearson_module,
  setdiff(pearson_mods, matched$pearson_module)
)
bm_order <- c(
  matched$bicor_module,
  setdiff(bicor_mods, matched$bicor_module)
)

jac_long <- jac_long |>
  mutate(
    pearson_module = factor(pearson_module, levels = rev(pm_order)),
    bicor_module = factor(bicor_module, levels = bm_order)
  )

p <- ggplot(jac_long, aes(bicor_module, pearson_module, fill = jaccard)) +
  geom_tile(colour = "white", linewidth = 0.3) +
  geom_text(aes(label = sprintf("%.2f", jaccard)),
    size = 2.8, colour = ifelse(jac_long$jaccard > 0.5, "white", "grey30")
  ) +
  scale_fill_gradient2(
    low = "white", mid = "#B2DFDB", high = "#00695C",
    midpoint = 0.35, limits = c(0, 1),
    name = "Jaccard\nIndex"
  ) +
  labs(
    x = "Bicor module", y = "Pearson module",
    title = "Pearson vs Bicor Module Overlap",
    subtitle = sprintf(
      "Mean matched Jaccard = %.3f  |  %d/%d modules > 0.5",
      mean(matched$jaccard),
      sum(matched$jaccard > 0.5), nrow(matched)
    )
  ) +
  FIG_THEME +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
    axis.text.y = element_text(size = 8),
    plot.title = element_text(size = 13, face = "bold"),
    plot.subtitle = element_text(size = 10)
  )

pdf_device <- get_pdf_device()

ggsave(file.path(RPT_PNG, "SUPP_bicor_sensitivity.png"), p,
  width = 180, height = 140, units = "mm", dpi = 300
)
ggsave(file.path(RPT_PDF, "SUPP_bicor_sensitivity.pdf"), p,
  width = 180, height = 140, units = "mm", device = pdf_device
)

message("  Bicor sensitivity analysis complete.")
