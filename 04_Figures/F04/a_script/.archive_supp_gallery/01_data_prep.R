# Supplementary Enrichment Gallery — Data Preparation (F04: Blunting)
# Loads pre-computed enrichment CSV, pivots to wide, classifies patterns.
# Saves prep_blunting.rds for downstream viz scripts.
setwd(rprojroot::find_rstudio_root_file())
source("04_Figures/shared/style.R")
source("04_Figures/shared/pathway_utils.R")

library(tidyverse)

DAT <- "04_Figures/F04/c_data/supp"
dir.create(DAT, recursive = TRUE, showWarnings = FALSE)

# Training blunting (TY vs TO vs Interaction)
blunt_long <- read_csv("04_Figures/F04/c_data/panel_supp/enrichment_blunting.csv",
                        show_col_types = FALSE)

blunt_wide <- blunt_long %>%
  select(pathway, pathway_label, database, contrast, NES, padj, sig, size) %>%
  pivot_wider(
    id_cols     = c(pathway, pathway_label, database),
    names_from  = contrast,
    values_from = c(NES, padj, sig, size)
  ) %>%
  rename(
    NES_TY  = NES_Training_Young,
    NES_TO  = NES_Training_Old,
    NES_Int = NES_Interaction,
    padj_TY  = padj_Training_Young,
    padj_TO  = padj_Training_Old,
    padj_Int = padj_Interaction,
    sig_TY   = sig_Training_Young,
    sig_TO   = sig_Training_Old,
    sig_Int  = sig_Interaction
  ) %>%
  mutate(
    sig_TY  = replace_na(sig_TY,  FALSE),
    sig_TO  = replace_na(sig_TO,  FALSE),
    sig_Int = replace_na(sig_Int, FALSE),
    pattern = case_when(
      sig_TY & sig_TO & sign(NES_TY) == sign(NES_TO) ~ "Concordant",
      sig_TY & sig_TO & sign(NES_TY) != sign(NES_TO) ~ "Discordant",
      sig_TY & !sig_TO & !sig_Int ~ "Blunted",
      !sig_TY & sig_TO            ~ "Old-specific",
      sig_Int                     ~ "Interaction-only",
      TRUE                        ~ "Other"
    ),
    pattern = factor(pattern, levels = c("Blunted", "Concordant", "Old-specific",
                                          "Interaction-only", "Discordant", "Other")),
    bio_theme = classify_pathway_func(pathway),
    bio_theme = factor(bio_theme, levels = CONSOLIDATED_PATHWAY_ORDER),
    set_size  = coalesce(size_Training_Young, size_Training_Old, size_Interaction)
  )

cat(sprintf("F04 blunting: %d pathways\n", nrow(blunt_wide)))
cat("Pattern counts:\n")
print(table(blunt_wide$pattern))

saveRDS(blunt_wide, file.path(DAT, "prep_blunting.rds"))

cat("\nData prep complete. RDS file saved to", DAT, "\n")
