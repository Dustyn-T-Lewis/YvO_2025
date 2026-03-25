# Figure 3 — Panel D: fGSEA Grouped Bar Chart (Pathway Enrichment)
# Hallmark + GO Slim + KEGG + Reactome + GO:BP curated collection
# DISPLAY + CACHE: per-database BH, no redundancy filter (raw results)
# Outputs: pD (ggplot object)

setwd(rprojroot::find_rstudio_root_file())
source("04_Figures/F03/a_script/style.R")
source("04_Figures/shared/pathway_utils.R")

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(tibble)
  library(tidyr)
  library(fgsea)
})

DEP_FILE <- "03_DEP/c_data/03_combined_results.csv"
RPT      <- "04_Figures/F03/b_reports"
DAT      <- "04_Figures/F03/c_data"
dir.create(RPT, recursive = TRUE, showWarnings = FALSE)
dir.create(DAT, recursive = TRUE, showWarnings = FALSE)

CONTRASTS <- c("Aging", "Training_Young", "Training_Old", "Interaction", "Reversal")
dep_df    <- read_csv(DEP_FILE, show_col_types = FALSE)
pdf_device <- get_pdf_device()
PC_W <- 110

# Build Hallmark + GO Slim + KEGG + Reactome curated collection
pw_collection <- build_pathway_collection(min_size = 10, max_size = 500)

set.seed(42)

# --- Per-database fGSEA with BH per database (NO Jaccard dedup) ---
# This gives us raw per-db BH results for display
fgsea_raw_all <- list()

for (ctr in CONTRASTS) {
  stats <- dep_df[[paste0("t_", ctr)]]
  names(stats) <- dep_df$gene
  stats <- sort(stats[!is.na(stats) & is.finite(stats)], decreasing = TRUE)

  pw_by_db <- split(names(pw_collection), classify_database(names(pw_collection)))

  ctr_results <- list()
  for (db_name in names(pw_by_db)) {
    db_pw <- pw_collection[pw_by_db[[db_name]]]
    if (length(db_pw) < 5) next
    raw <- fgseaMultilevel(
      pathways    = db_pw,
      stats       = stats,
      minSize     = 10,
      maxSize     = 500,
      nPermSimple = 10000,
      eps         = 0
    )
    raw <- as.data.frame(raw)
    raw$database <- db_name
    ctr_results[[db_name]] <- as_tibble(raw)
  }
  combined <- bind_rows(ctr_results)
  combined$contrast <- ctr
  fgsea_raw_all[[ctr]] <- combined
}

fgsea_raw <- bind_rows(fgsea_raw_all)

# --- Export raw per-database BH results (no Jaccard dedup) ---
fgsea_export <- fgsea_raw |>
  mutate(leadingEdge = sapply(leadingEdge, paste, collapse = ";")) |>
  arrange(database, contrast, padj)
write_csv(fgsea_export, file.path(DAT, "01_panel_D_fgsea_results.csv"))

# DEPENDENCY: F04 and F05 volcano panels read this cache.
# Run this script before F04/panel_A-C.R and F05/panel_A-B.R.
cache_dir_f2 <- "04_Figures/F04/c_data/shared"
cache_dir_f3 <- "04_Figures/F05/c_data/shared"
dir.create(cache_dir_f2, recursive = TRUE, showWarnings = FALSE)
dir.create(cache_dir_f3, recursive = TRUE, showWarnings = FALSE)
write_csv(fgsea_export, file.path(cache_dir_f2, "fgsea_tstat_all_v2.csv"))
write_csv(fgsea_export, file.path(cache_dir_f3, "fgsea_tstat_all_v2.csv"))

# --- Full GO:BP fGSEA (separate cache for multi-database volcano rings) ---
message("\n--- Running full GO:BP fGSEA ---")
fgsea_gobp_all <- list()
for (ctr in CONTRASTS) {
  stats <- dep_df[[paste0("t_", ctr)]]
  names(stats) <- dep_df$gene
  stats <- sort(stats[!is.na(stats) & is.finite(stats)], decreasing = TRUE)

  gobp_pw <- pw_collection[grepl("^GOBP_", names(pw_collection))]
  res <- fgseaMultilevel(
    pathways    = gobp_pw,
    stats       = stats,
    minSize     = 15,
    maxSize     = 500,
    nPermSimple = 10000,
    eps         = 0
  )
  res <- as_tibble(as.data.frame(res))
  res$database <- "GO:BP"
  res$contrast <- ctr
  fgsea_gobp_all[[ctr]] <- res
}

fgsea_gobp_combined <- bind_rows(fgsea_gobp_all)
fgsea_gobp_export <- fgsea_gobp_combined |>
  mutate(leadingEdge = sapply(leadingEdge, paste, collapse = ";")) |>
  arrange(contrast, padj)

write_csv(fgsea_gobp_export, file.path(DAT, "02_panel_D_fgsea_gobp.csv"))
# DEPENDENCY: F04 and F05 volcano panels read this cache.
# Run this script before F04/panel_A-C.R and F05/panel_A-B.R.
write_csv(fgsea_gobp_export, file.path(cache_dir_f2, "fgsea_gobp_v2.csv"))
write_csv(fgsea_gobp_export, file.path(cache_dir_f3, "fgsea_gobp_v2.csv"))

n_gobp_sig <- sum(!is.na(fgsea_gobp_combined$padj) & fgsea_gobp_combined$padj < 0.05)
message(sprintf("GO:BP cache: %d total rows, %d significant across %d contrasts",
                nrow(fgsea_gobp_combined), n_gobp_sig, length(CONTRASTS)))

# --- DISPLAY: use raw (pre-dedup) counts, exclude Reversal ---
DISPLAY_DBS       <- c("Hallmark", "GO Slim", "GO:BP", "KEGG", "Reactome")
DISPLAY_CONTRASTS <- c("Aging", "Training_Young", "Training_Old", "Interaction")

# Per-contrast, per-database totals (fgsea tests different counts per contrast
# because its internal size filter depends on gene overlap with ranked list)
db_ctr_totals <- fgsea_raw |>
  filter(database %in% DISPLAY_DBS, contrast %in% DISPLAY_CONTRASTS) |>
  group_by(contrast, database) |>
  summarise(n_total = n(), .groups = "drop")

count_df <- fgsea_raw |>
  filter(!is.na(padj), padj < 0.05,
         database %in% DISPLAY_DBS, contrast %in% DISPLAY_CONTRASTS) |>
  group_by(contrast, database) |>
  summarise(
    Up   = sum(NES > 0),
    Down = sum(NES < 0),
    .groups = "drop"
  ) |>
  tidyr::pivot_longer(cols = c(Up, Down), names_to = "direction",
                      values_to = "count") |>
  left_join(db_ctr_totals, by = c("contrast", "database")) |>
  mutate(fraction = count / n_total)

nonempty_dbs <- count_df |>
  group_by(database) |> filter(sum(count) > 0) |> pull(database) |> unique()
count_df <- count_df |> filter(database %in% nonempty_dbs)

# Facet labels: use median tested count per database across contrasts
db_label_n <- db_ctr_totals |>
  filter(database %in% nonempty_dbs) |>
  group_by(database) |>
  summarise(n_label = as.integer(median(n_total)), .groups = "drop")
db_labels <- setNames(
  sprintf("%s (n=%d testable)", db_label_n$database, db_label_n$n_label),
  db_label_n$database
)

count_df$contrast  <- factor(count_df$contrast, levels = DISPLAY_CONTRASTS)
count_df$database  <- factor(count_df$database, levels = intersect(DISPLAY_DBS, nonempty_dbs))
count_df$direction <- factor(count_df$direction, levels = c("Up", "Down"))

# --- Pathway-level blunting: Fisher's exact on sig/non-sig x Tr.(Y)/Tr.(O) ---
# Uses per-contrast totals (fgsea tests slightly different pathway counts per contrast)
sig_ty <- fgsea_raw |>
  filter(contrast == "Training_Young", database %in% DISPLAY_DBS) |>
  summarise(sig = sum(!is.na(padj) & padj < 0.05), total = n())
sig_to <- fgsea_raw |>
  filter(contrast == "Training_Old", database %in% DISPLAY_DBS) |>
  summarise(sig = sum(!is.na(padj) & padj < 0.05), total = n())

blunt_mat <- matrix(
  c(sig_ty$sig, sig_ty$total - sig_ty$sig,
    sig_to$sig, sig_to$total - sig_to$sig),
  nrow = 2, byrow = TRUE,
  dimnames = list(c("Tr.(Y)", "Tr.(O)"), c("sig", "non-sig"))
)
blunt_fisher <- fisher.test(blunt_mat, alternative = "greater")
blunt_ratio  <- sprintf("%.1fx", sig_ty$sig / max(sig_to$sig, 1))
blunt_p      <- fmt_p(blunt_fisher$p.value)

n_facets <- length(levels(count_df$database))
PC_H <- max(80, n_facets * 55)

lbl_sz <- scale_text(BASE_COUNT, PC_W)

pD <- ggplot(count_df, aes(x = contrast, y = fraction * 100, fill = direction)) +
  annotate("rect", xmin = 0.5, xmax = 1.5, ymin = -Inf, ymax = Inf,
           fill = CONTRAST_COLORS["Aging"], alpha = 0.20,
           color = "grey70", linewidth = 0.2) +
  annotate("rect", xmin = 1.5, xmax = 2.5, ymin = -Inf, ymax = Inf,
           fill = CONTRAST_COLORS["Training_Young"], alpha = 0.20,
           color = "grey70", linewidth = 0.2) +
  annotate("rect", xmin = 2.5, xmax = 3.5, ymin = -Inf, ymax = Inf,
           fill = CONTRAST_COLORS["Training_Old"], alpha = 0.20,
           color = "grey70", linewidth = 0.2) +
  annotate("rect", xmin = 3.5, xmax = 4.5, ymin = -Inf, ymax = Inf,
           fill = CONTRAST_COLORS["Interaction"], alpha = 0.20,
           color = "grey70", linewidth = 0.2) +
  geom_col(position = position_dodge(width = 0.7), width = 0.6,
           color = "black", linewidth = 0.3) +
  geom_text(aes(y = fraction * 100 / 2,
                label = ifelse(count > 0, count, "")),
            position = position_dodge(width = 0.7),
            vjust = 0.5, hjust = 0.5, size = lbl_sz,
            color = "white", fontface = "bold", show.legend = FALSE) +
  facet_wrap(~ database, ncol = 1, scales = "free_y", strip.position = "top",
             labeller = as_labeller(db_labels)) +
  scale_x_discrete(labels = CTR_SHORT) +
  scale_fill_manual(values = DIR_COLORS) +
  labs(title = "Fraction of Pathways Enriched Up/Down",
       subtitle = "fGSEA per-database BH corrected fraction",
       x = NULL, y = "% of database significant",
       tag = "D") +
  FIG_THEME +
  theme(axis.text.x    = element_text(angle = 35, hjust = 1,
                                     size = FIG_AXIS_TEXT - 0.5),
        legend.position = "none",
        strip.text     = element_text(size = FIG_STRIP_SIZE - 1, face = "bold",
                                      margin = margin(1, 0, 1, 0)),
        strip.placement = "outside")

nes_summary <- fgsea_raw |>
  filter(!is.na(padj), padj < 0.05) |>
  group_by(contrast, database) |>
  summarise(
    n_sig = n(), n_up = sum(NES > 0), n_down = sum(NES < 0),
    median_NES = median(NES), mean_NES = mean(NES), sd_NES = sd(NES),
    min_padj = min(padj), median_padj = median(padj),
    .groups = "drop"
  )
write.csv(nes_summary, file.path(DAT, "panel_D_nes_summary.csv"),
          row.names = FALSE)

ggsave(file.path(RPT, "panel_D_fgsea.pdf"), pD,
       width = PC_W, height = PC_H, units = "mm", device = pdf_device)
ggsave(file.path(RPT, "panel_D_fgsea.png"), pD,
       width = PC_W, height = PC_H, units = "mm", dpi = 300)
