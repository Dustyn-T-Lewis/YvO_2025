# Supplementary Enrichment Gallery — Data Preparation (F05: Reversal)
# Loads pre-computed enrichment CSV, pivots to wide, classifies patterns.
# Saves prep_reversal.rds for downstream viz scripts.
setwd(rprojroot::find_rstudio_root_file())
source("04_Figures/shared/style.R")
source("04_Figures/shared/pathway_utils.R")

library(tidyverse)

DAT <- "04_Figures/F05/c_data/supp"
dir.create(DAT, recursive = TRUE, showWarnings = FALSE)

# Aging reversal (Aging vs Training_Old)
rev_long <- read_csv("04_Figures/F05/c_data/panel_supp/enrichment_reversal.csv",
                      show_col_types = FALSE)

rev_wide <- rev_long %>%
  select(pathway, pathway_label, database, contrast, NES, padj, sig, size) %>%
  pivot_wider(
    id_cols     = c(pathway, pathway_label, database),
    names_from  = contrast,
    values_from = c(NES, padj, sig, size)
  ) %>%
  rename(
    NES_Aging = NES_Aging,
    NES_TO    = NES_Training_Old,
    padj_Aging = padj_Aging,
    padj_TO    = padj_Training_Old,
    sig_Aging  = sig_Aging,
    sig_TO     = sig_Training_Old
  ) %>%
  mutate(
    sig_Aging = replace_na(sig_Aging, FALSE),
    sig_TO    = replace_na(sig_TO,    FALSE),
    reversal_ratio = ifelse(abs(NES_Aging) > 0.01, -NES_TO / NES_Aging, NA_real_),
    pattern = case_when(
      sig_Aging & sig_TO & sign(NES_Aging) != sign(NES_TO) ~ "Reversed",
      sig_Aging & sig_TO & sign(NES_Aging) == sign(NES_TO) ~ "Exacerbated",
      sig_Aging & !sig_TO ~ "Aging-specific",
      !sig_Aging & sig_TO ~ "Training-specific",
      TRUE                ~ "Other"
    ),
    pattern = factor(pattern, levels = c("Reversed", "Exacerbated", "Aging-specific",
                                          "Training-specific", "Other")),
    bio_theme = classify_pathway_func(pathway),
    bio_theme = factor(bio_theme, levels = CONSOLIDATED_PATHWAY_ORDER),
    set_size  = coalesce(size_Aging, size_Training_Old)
  )

cat(sprintf("F05 reversal: %d pathways\n", nrow(rev_wide)))
cat("Pattern counts:\n")
print(table(rev_wide$pattern))

saveRDS(rev_wide, file.path(DAT, "prep_reversal.rds"))

cat("\nData prep complete. RDS file saved to", DAT, "\n")
