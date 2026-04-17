# F05 Panel C ORA: Over-Representation Analysis on Reversal Pattern Groups
# Complement to panel_C heatmap — runs ORA on each of the 4 reversal classes:
#   Aging Up, Aging Down, Rev. (Age Up), Rev. (Age Down)
setwd(rprojroot::find_rstudio_root_file())
source("04_Figures/shared/style.R")
source("04_Figures/shared/pathway_utils.R")
library(readr)
library(dplyr)

DAT <- "04_Figures/F05/c_data"
dir.create(file.path(DAT, "panel_C"), recursive = TRUE, showWarnings = FALSE)

# --- Pattern classification (same as panel_C.R)
F3_CLASS_ORDER <- c("Aging Up", "Aging Down", "Rev. (Age Up)", "Rev. (Age Down)")

classify_f3 <- function(pi_aging, pi_training_old, logFC_aging, logFC_training_old) {
  dplyr::case_when(
    pi_aging < 0.05 & pi_training_old < 0.05 &
      sign(logFC_aging) != sign(logFC_training_old) &
      logFC_aging > 0  ~ "Rev. (Age Up)",
    pi_aging < 0.05 & pi_training_old < 0.05 &
      sign(logFC_aging) != sign(logFC_training_old) &
      logFC_aging <= 0 ~ "Rev. (Age Down)",
    pi_aging < 0.05 & logFC_aging > 0     ~ "Aging Up",
    pi_aging < 0.05 & logFC_aging <= 0    ~ "Aging Down",
    TRUE ~ NA_character_
  )
}

# --- Load and classify proteins
dep_df <- read_csv("03_DEP/c_data/03_combined_results.csv", show_col_types = FALSE)

pattern_df <- dep_df %>%
  filter(!is.na(logFC_Aging), !is.na(logFC_Training_Old)) %>%
  filter(pi_score_Aging < 0.05 | pi_score_Training_Old < 0.05) %>%
  mutate(pattern = classify_f3(pi_score_Aging, pi_score_Training_Old, logFC_Aging, logFC_Training_Old)) %>%
  filter(!is.na(pattern))

universe <- dep_df$gene[!is.na(dep_df$logFC_Aging)]
message(sprintf("  %d classified proteins | universe: %d", nrow(pattern_df), length(universe)))

# --- Build pathway collection
pw_collection <- build_pathway_collection(min_size = 15, max_size = 500,
                                           include_goslim = FALSE,
                                           exclude_variants = TRUE)

# --- ORA helper
MAX_BARS <- 12

run_set_ora <- function(genes, set_name) {
  if (length(genes) < 5) return(tibble())
  res <- tryCatch(
    run_ora_deduplicated(genes = genes, universe = universe,
                          pathways = pw_collection, jaccard_cutoff = 0.5,
                          min_size = 15, max_size = 500, padj_cutoff = 0.05),
    error = function(e) { message("  ORA error: ", e$message); tibble() }
  )
  if (nrow(res) > 0) {
    res %>% mutate(
      set = set_name,
      pathway_label = clean_pathway_name(pathway),
      neg_log10_padj = -log10(padj)
    ) %>% arrange(desc(neg_log10_padj)) %>% slice_head(n = MAX_BARS)
  } else tibble()
}

# --- Run ORA per pattern group
message("\n--- Pattern-group ORA ---")
pattern_results <- list()

for (pat in F3_CLASS_ORDER) {
  genes <- pattern_df$gene[pattern_df$pattern == pat]
  slug <- tolower(gsub("[- .()]", "_", pat))
  slug <- gsub("__+", "_", gsub("_$", "", slug))
  message(sprintf("  %s: %d genes", pat, length(genes)))
  ora <- run_set_ora(genes, pat)
  pattern_results[[slug]] <- ora
  if (nrow(ora) > 0) {
    message(sprintf("  %s: %d pathways enriched (%d genes)", slug, nrow(ora), length(genes)))
  } else {
    message(sprintf("  %s: no significant pathways (n=%d genes)", pat, length(genes)))
  }
}

# --- Save data
all_pattern_ora <- bind_rows(pattern_results)
if (nrow(all_pattern_ora) > 0)
  write_csv(all_pattern_ora, file.path(DAT, "panel_C", "ora_pattern_groups.csv"))

message("\nF05 Panel C ORA done")
