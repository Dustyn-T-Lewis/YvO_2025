# Supplementary Figure: GO Slim Classification Documentation
# Output: shared/b_reports/Supp_GO_Slim_mapping.{pdf,png}
setwd(rprojroot::find_rstudio_root_file())
source("04_Figures/shared/style.R")
source("04_Figures/shared/go_slim_categories.R")

suppressPackageStartupMessages({
  library(tidyverse)
  library(gridExtra)
  library(grid)
  library(patchwork)
})

RPT <- "04_Figures/shared/b_reports"
dir.create(RPT, recursive = TRUE, showWarnings = FALSE)

pdf_device <- get_pdf_device()

dep_df <- read_csv("03_DEP/c_data/03_combined_results.csv", show_col_types = FALSE)
all_genes <- dep_df$gene

f2_genes <- dep_df %>%
  filter(!is.na(logFC_Training_Young), !is.na(logFC_Training_Old)) %>%
  mutate(
    sig_cat = classify_proteins_f2(
      pi_score_Training_Young, pi_score_Training_Old, pi_score_Interaction
    )
  ) %>%
  filter(sig_cat != "NS") %>%
  pull(gene)

f3_genes <- dep_df %>%
  filter(!is.na(logFC_Aging), !is.na(logFC_Training_Old)) %>%
  filter(pi_score_Aging < 0.05 | pi_score_Reversal < 0.05) %>%
  pull(gene)

message(sprintf("F2: %d significant genes, F3: %d significant genes",
                length(f2_genes), length(f3_genes)))

f2_result <- assign_go_slim_consolidated(f2_genes, all_genes)
f3_result <- assign_go_slim_consolidated(f3_genes, all_genes)

slim_terms <- AnnotationDbi::Term(bp_slim)

mapping_df <- tibble(
  go_id = names(SLIM_CONSOLIDATED),
  go_term = slim_terms[names(SLIM_CONSOLIDATED)],
  category = unname(SLIM_CONSOLIDATED)
) %>%
  arrange(match(category, CONSOLIDATED_PATHWAY_ORDER), go_id)

excluded_ids <- setdiff(bp_slim, names(SLIM_CONSOLIDATED))
excluded_df <- tibble(
  go_id = excluded_ids,
  go_term = slim_terms[excluded_ids],
  category = "— Excluded —"
)

full_mapping <- bind_rows(mapping_df, excluded_df)

table_df <- full_mapping %>%
  group_by(category) %>%
  summarise(
    go_terms = paste(go_term, collapse = "\n"),
    go_ids   = paste(go_id, collapse = "\n"),
    .groups  = "drop"
  ) %>%
  rename(Category = category, `GO Slim Terms` = go_terms, `GO IDs` = go_ids)

cat_colors <- c(CONSOLIDATED_COLORS, "— Excluded —" = "#FFCDD2")

tt <- gridExtra::ttheme_minimal(
  base_size = 7,
  core = list(
    fg_params = list(hjust = 0, x = 0.02, fontsize = 6.5),
    bg_params = list(fill = "white", col = "grey80", lwd = 0.5)
  ),
  colhead = list(
    fg_params = list(fontface = "bold", fontsize = 7.5, hjust = 0, x = 0.02),
    bg_params = list(fill = "grey90", col = "grey60", lwd = 0.5)
  )
)

tbl_grob <- gridExtra::tableGrob(table_df, rows = NULL, theme = tt)

for (i in seq_len(nrow(table_df))) {
  cat_name <- table_df$Category[i]
  if (cat_name %in% names(cat_colors)) {
    bg_col <- scales::alpha(cat_colors[cat_name], 0.35)
    # Row index in grob is i+1 (header is row 1)
    idx <- which(tbl_grob$layout$name == paste0("core-bg.", i, "-", 1))
    if (length(idx) > 0) {
      tbl_grob$grobs[[idx]]$gp$fill <- bg_col
    }
  }
}

pA <- ggplot() +
  annotation_custom(tbl_grob) +
  labs(title = "A. GO Slim Generic BP → 15 Consolidated Categories",
       subtitle = sprintf("62 curated GO Slim terms, 2 excluded (signaling, nervous system), %d mapped",
                          length(SLIM_CONSOLIDATED))) +
  theme_void() +
  theme(
    plot.title = element_text(face = "bold", size = 11, hjust = 0),
    plot.subtitle = element_text(size = 8, color = "grey40", hjust = 0),
    plot.margin = margin(5, 5, 5, 5)
  )

f2_mapped <- sum(f2_result$consolidated != "Other")
f3_mapped <- sum(f3_result$consolidated != "Other")

coverage_df <- tibble(
  Figure = c("F2 (Training)", "F3 (Aging/Reversal)"),
  Total  = c(nrow(f2_result), nrow(f3_result)),
  Mapped = c(f2_mapped, f3_mapped),
  Other  = c(nrow(f2_result) - f2_mapped, nrow(f3_result) - f3_mapped),
  Pct    = c(100 * f2_mapped / nrow(f2_result),
             100 * f3_mapped / nrow(f3_result))
)

f2_cats <- f2_result %>%
  count(consolidated, name = "F2") %>%
  rename(Category = consolidated)
f3_cats <- f3_result %>%
  count(consolidated, name = "F3") %>%
  rename(Category = consolidated)

cat_comparison <- full_join(f2_cats, f3_cats, by = "Category") %>%
  mutate(across(c(F2, F3), ~replace_na(.x, 0L))) %>%
  mutate(Category = factor(Category, levels = CONSOLIDATED_PATHWAY_ORDER)) %>%
  arrange(Category)

cat_long <- cat_comparison %>%
  pivot_longer(cols = c(F2, F3), names_to = "Figure", values_to = "n") %>%
  filter(n > 0)

pB <- ggplot(cat_long, aes(x = Figure, y = n, fill = Category)) +
  geom_col(position = "stack", width = 0.6, color = "white", linewidth = 0.3) +
  geom_text(aes(label = ifelse(n >= 3, n, "")),
            position = position_stack(vjust = 0.5),
            size = 2.5, fontface = "bold") +
  scale_fill_manual(values = CONSOLIDATED_COLORS, drop = FALSE) +
  labs(title = "B. Category Sizes per Figure",
       subtitle = sprintf("F2: %d/%d mapped (%.0f%%) | F3: %d/%d mapped (%.0f%%)",
                          coverage_df$Mapped[1], coverage_df$Total[1], coverage_df$Pct[1],
                          coverage_df$Mapped[2], coverage_df$Total[2], coverage_df$Pct[2]),
       x = NULL, y = "Number of proteins") +
  FIG_THEME +
  theme(
    legend.position = "right",
    legend.key.size = unit(3, "mm"),
    legend.text = element_text(size = 7),
    legend.title = element_blank(),
    plot.title = element_text(face = "bold", size = 11),
    plot.subtitle = element_text(size = 8, color = "grey40")
  ) +
  guides(fill = guide_legend(ncol = 1, reverse = TRUE))

composite <- pA / pB + plot_layout(heights = c(3, 2))

ggsave(file.path(RPT, "Supp_GO_Slim_mapping.pdf"), composite,
       width = 210, height = 297, units = "mm", device = pdf_device)
ggsave(file.path(RPT, "Supp_GO_Slim_mapping.png"), composite,
       width = 210, height = 297, units = "mm", dpi = 300)

write.csv(full_mapping, file.path(RPT, "go_slim_mapping_table.csv"), row.names = FALSE)
write.csv(coverage_df, file.path(RPT, "go_slim_coverage.csv"), row.names = FALSE)
write.csv(cat_comparison, file.path(RPT, "go_slim_category_comparison.csv"), row.names = FALSE)

cat("Supplementary GO Slim mapping figure done\n")
