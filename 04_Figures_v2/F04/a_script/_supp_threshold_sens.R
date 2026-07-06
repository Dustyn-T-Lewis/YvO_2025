#!/usr/bin/env Rscript
# SUPP Panel C: Threshold Sensitivity
# Diagnostic for main Panel B — shows concordance pattern is stable across
# significance thresholds (Pi < 0.05, FDR < 0.05, FDR < 0.10, nominal p < 0.05).
# Sourced by 02_supp_panels.R — expects style.R already loaded.
# Exports: pS_thresh (ggplot)

BASE <- "04_Figures_v2/F04"
DAT <- file.path(BASE, "c_data", "panel_supp")
dir.create(DAT, recursive = TRUE, showWarnings = FALSE)

dep_df <- read_csv("03_DEP/c_data/03_combined_results.csv",
  show_col_types = FALSE
)

base_df <- dep_df |>
  transmute(gene,
    logFC_TY   = logFC_Training_Young,
    logFC_TO   = logFC_Training_Old,
    pi_TY      = pi_score_Training_Young,
    pi_TO      = pi_score_Training_Old,
    fdr_TY     = adj.P.Val_Training_Young,
    fdr_TO     = adj.P.Val_Training_Old,
    nom_TY     = P.Value_Training_Young,
    nom_TO     = P.Value_Training_Old
  ) |>
  filter(!is.na(logFC_TY), !is.na(logFC_TO))

# Quadrant assignment
assign_quad <- function(lfc_ty, lfc_to) {
  case_when(
    lfc_ty > 0 & lfc_to > 0 ~ "Concordant Up",
    lfc_ty < 0 & lfc_to < 0 ~ "Concordant Down",
    TRUE ~ "Discordant"
  )
}

# Threshold sweep
thresholds <- list(
  "Pi < 0.05" = function(d) d |> filter(pi_TY < 0.05 | pi_TO < 0.05),
  "FDR < 0.05" = function(d) d |> filter(fdr_TY < 0.05 | fdr_TO < 0.05),
  "FDR < 0.10" = function(d) d |> filter(fdr_TY < 0.10 | fdr_TO < 0.10),
  "Nom. p < 0.05" = function(d) d |> filter(nom_TY < 0.05 | nom_TO < 0.05)
)

results <- vector("list", length(thresholds))
for (i in seq_along(thresholds)) {
  thr_name <- names(thresholds)[i]
  sig_sub <- thresholds[[thr_name]](base_df)
  sig_sub <- sig_sub |> mutate(quadrant = assign_quad(logFC_TY, logFC_TO))

  counts <- sig_sub |> count(quadrant, name = "n")
  counts$threshold <- thr_name
  counts$n_total <- nrow(sig_sub)
  results[[i]] <- counts
}

sens_df <- bind_rows(results) |>
  mutate(
    threshold = factor(threshold, levels = names(thresholds)),
    quadrant = factor(quadrant, levels = c(
      "Concordant Up",
      "Concordant Down",
      "Discordant"
    ))
  )

write_csv(sens_df, file.path(DAT, "SUPP_threshold_sensitivity.csv"))
message("Threshold sensitivity:\n", paste(capture.output(print(sens_df)), collapse = "\n"))

QUAD_COLS <- c(
  "Concordant Up" = "#E57373",
  "Concordant Down" = "#64B5F6",
  "Discordant" = "#FFB74D"
)

pS_thresh <- ggplot(sens_df, aes(x = threshold, y = n, fill = quadrant)) +
  geom_col(
    position = "stack", width = 0.65,
    color = "white", linewidth = 0.3
  ) +
  geom_text(aes(label = n),
    position = position_stack(vjust = 0.5),
    size = BASE_COUNT, fontface = "bold", color = "white"
  ) +
  scale_fill_manual(values = QUAD_COLS, name = "Quadrant") +
  labs(
    title = "Concordance Pattern: Threshold Sensitivity",
    subtitle = "Protein counts per quadrant across significance criteria",
    x = NULL, y = "Significant proteins"
  ) +
  FIG_THEME

RPT_PNG <- file.path(BASE, "b_reports", "supp", "png", "panels")
RPT_PDF <- file.path(BASE, "b_reports", "supp", "pdf", "panels")
dir.create(RPT_PNG, recursive = TRUE, showWarnings = FALSE)
dir.create(RPT_PDF, recursive = TRUE, showWarnings = FALSE)

PW <- 89
PH <- 70
ggsave(file.path(RPT_PNG, "SUPP_threshold_sens.png"), pS_thresh,
  width = PW, height = PH, units = "mm", dpi = 300
)
ggsave(file.path(RPT_PDF, "SUPP_threshold_sens.pdf"), pS_thresh,
  width = PW, height = PH, units = "mm", device = pdf_device
)

pS_thresh <- strip_for_composite(pS_thresh)

message("SUPP Panel C (threshold sensitivity) done")
