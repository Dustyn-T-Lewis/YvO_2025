# F04 Panel C ORA: Over-Representation Analysis on Pattern Classification Groups
# Complement to panel_C heatmap — runs ORA on each of the 4 response pattern groups:
#   Shared (sig TY + TO), Blunted (sig TY only), Age-resistant (sig TO only), Interaction only
setwd(rprojroot::find_rstudio_root_file())
source("04_Figures/shared/style.R")
source("04_Figures/shared/pathway_utils.R")
library(readr)
library(dplyr)

DAT <- "04_Figures/F04/c_data"
dir.create(file.path(DAT, "panel_C"), recursive = TRUE, showWarnings = FALSE)

# --- Pattern classification (same as panel_C.R)
PATTERN_ORDER <- c("Shared", "Blunted", "Age-resistant", "Interaction only")

classify_pattern_f2 <- function(pi_Y, pi_O, pi_int, threshold = 0.05) {
  dplyr::case_when(
    pi_Y < threshold & pi_O < threshold  ~ "Shared",
    pi_Y < threshold & pi_O >= threshold ~ "Blunted",
    pi_Y >= threshold & pi_O < threshold ~ "Age-resistant",
    pi_int < threshold                   ~ "Interaction only",
    TRUE                                 ~ NA_character_
  )
}

# --- Load and classify proteins
dep_df <- read_csv("03_DEP/c_data/03_combined_results.csv", show_col_types = FALSE)

pattern_df <- dep_df %>%
  filter(!is.na(logFC_Training_Young), !is.na(logFC_Training_Old)) %>%
  mutate(pattern = classify_pattern_f2(
    pi_score_Training_Young, pi_score_Training_Old, pi_score_Interaction
  )) %>%
  filter(!is.na(pattern))

universe <- dep_df$gene[!is.na(dep_df$logFC_Training_Young)]
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

for (pat in PATTERN_ORDER) {
  genes <- pattern_df$gene[pattern_df$pattern == pat]
  slug <- tolower(gsub("[- ]", "_", pat))
  message(sprintf("  %s: %d genes", pat, length(genes)))
  ora <- run_set_ora(genes, pat)
  pattern_results[[slug]] <- ora
  if (nrow(ora) > 0) {
    message(sprintf("  %s: %d significant pathways (%d genes)", slug, nrow(ora), length(genes)))
  } else {
    message(sprintf("  %s: no significant pathways (n=%d genes)", pat, length(genes)))
  }
}

# --- Save data
all_pattern_ora <- bind_rows(pattern_results)
if (nrow(all_pattern_ora) > 0)
  write_csv(all_pattern_ora, file.path(DAT, "panel_C", "ora_pattern_groups.csv"))

message("\nF04 Panel C ORA done")
