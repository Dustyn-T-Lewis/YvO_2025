#!/usr/bin/env Rscript
# SUPP Panel D: GO Slim Category Distribution
# Diagnostic for main Panel B — quantitative GO Slim category breakdown by
# concordance quadrant (Concordant Up / Concordant Down / Discordant).
# Sourced by 02_supp_panels.R — expects style.R already loaded.
# Exports: pS_goslim (ggplot)

source("04_Figures/shared/go_slim_categories.R")
library(readxl)

BASE <- "04_Figures/F04"
DAT  <- file.path(BASE, "c_data", "panel_supp")
dir.create(DAT, recursive = TRUE, showWarnings = FALSE)

xlsx_path <- file.path(BASE, "c_data", "F04_supplementary.xlsx")
pattern_df <- tryCatch({
  df <- read_excel(xlsx_path, sheet = "panel_B_pattern_class")
  message("Read pattern classification from xlsx")
  as.data.frame(df)
}, error = function(e) {
  message("xlsx read failed, recomputing from DEP results: ", e$message)
  NULL
})

if (is.null(pattern_df)) {
  dep_df <- read_csv("03_DEP/c_data/03_combined_results.csv",
                     show_col_types = FALSE)
  pattern_df <- dep_df |>
    filter(!is.na(logFC_Training_Young), !is.na(logFC_Training_Old)) |>
    filter(pi_score_Training_Young < 0.05 | pi_score_Training_Old < 0.05 |
           pi_score_Interaction < 0.05) |>
    mutate(quadrant = case_when(
      logFC_Training_Young > 0 & logFC_Training_Old > 0 ~ "Concordant Up",
      logFC_Training_Young < 0 & logFC_Training_Old < 0 ~ "Concordant Down",
      TRUE ~ "Discordant")) |>
    as.data.frame()
}

if (!"quadrant" %in% names(pattern_df)) {
  pattern_df <- pattern_df |>
    mutate(quadrant = case_when(
      logFC_Training_Young > 0 & logFC_Training_Old > 0 ~ "Concordant Up",
      logFC_Training_Young < 0 & logFC_Training_Old < 0 ~ "Concordant Down",
      TRUE ~ "Discordant"))
}

# GO Slim assignment
all_genes <- if ("gene" %in% names(pattern_df)) pattern_df$gene else pattern_df$SYMBOL
fg_genes  <- all_genes

slim_result <- assign_go_slim_consolidated(fg_genes = fg_genes,
                                            all_genes = fg_genes)

gene_col_nm <- if ("gene" %in% names(pattern_df)) "gene" else "SYMBOL"
pattern_quad <- transmute(pattern_df, gene = .data[[gene_col_nm]], quadrant)
slim_merged <- slim_result |>
  left_join(pattern_quad, by = "gene") |>
  filter(!is.na(quadrant), !is.na(consolidated))

# Summarise per category x quadrant
cat_quad <- slim_merged |>
  count(consolidated, quadrant, name = "n") |>
  mutate(quadrant = factor(quadrant,
                           levels = c("Concordant Up", "Concordant Down",
                                      "Discordant")))

cat_order <- cat_quad |>
  group_by(consolidated) |>
  summarise(total = sum(n), .groups = "drop") |>
  arrange(total) |>
  pull(consolidated)
cat_order <- c("Other", setdiff(cat_order, "Other"))

cat_quad <- cat_quad |>
  mutate(consolidated = factor(consolidated, levels = cat_order))

write_csv(cat_quad, file.path(DAT, "SUPP_goslim_distribution.csv"))
message("GO Slim distribution:\n", paste(capture.output(print(cat_quad)), collapse = "\n"))

QUAD_COLS <- c("Concordant Up" = "#E57373",
               "Concordant Down" = "#64B5F6",
               "Discordant" = "#FFB74D")

pS_goslim <- ggplot(cat_quad, aes(x = n, y = consolidated, fill = quadrant)) +
  geom_col(position = "stack", width = 0.7,
           color = "white", linewidth = 0.3) +
  scale_fill_manual(values = QUAD_COLS, name = "Quadrant") +
  labs(title = "GO Slim Category Distribution",
       subtitle = "Protein counts per GO Slim category by concordance quadrant",
       x = "Protein count", y = NULL) +
  FIG_THEME

RPT_PNG <- file.path(BASE, "b_reports", "supp", "png", "panels")
RPT_PDF <- file.path(BASE, "b_reports", "supp", "pdf", "panels")
dir.create(RPT_PNG, recursive = TRUE, showWarnings = FALSE)
dir.create(RPT_PDF, recursive = TRUE, showWarnings = FALSE)

PW <- 89; PH <- 70
ggsave(file.path(RPT_PNG, "SUPP_goslim_bars.png"), pS_goslim,
       width = PW, height = PH, units = "mm", dpi = 300)
ggsave(file.path(RPT_PDF, "SUPP_goslim_bars.pdf"), pS_goslim,
       width = PW, height = PH, units = "mm", device = pdf_device)

pS_goslim <- strip_for_composite(pS_goslim)

message("SUPP Panel D (GO Slim bars) done")
