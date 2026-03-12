# Figure 1 — Panel H: fGSEA Grouped Bar Chart (Pathway Enrichment)
# Full collection: Hallmark + C2:CP (all subcollections) + GO:BP
# Per-database BH correction. Exports fgsea cache for F2/F3.
# Outputs: pH (ggplot object)

setwd(rprojroot::find_rstudio_root_file())
source("04_Figures/shared/style.R")
source("04_Figures/shared/pathway_utils.R")

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(tibble)
  library(tidyr)
  library(fgsea)
  library(msigdbr)
})

DEP_FILE <- "03_DEP/c_data/03_combined_results.csv"
RPT      <- "04_Figures/F1/b_reports"
DAT      <- "04_Figures/F1/c_data"
dir.create(RPT, recursive = TRUE, showWarnings = FALSE)
dir.create(DAT, recursive = TRUE, showWarnings = FALSE)

CONTRASTS <- c("Aging", "Training_Young", "Training_Old", "Interaction")
dep_df    <- read_csv(DEP_FILE, show_col_types = FALSE)
pdf_device <- get_pdf_device()
PH_W <- 120

# Full C2:CP subcollections (KEGG, Reactome, WikiPathways, BioCarta, PID)

build_full_collection <- function(species = "Homo sapiens",
                                   min_size = 10, max_size = 500) {
  hallmark <- msigdbr::msigdbr(species = species, collection = "H")
  gobp     <- msigdbr::msigdbr(species = species, collection = "C5",
                                 subcollection = "GO:BP")

  cp_subs <- c("CP", "CP:KEGG_MEDICUS", "CP:REACTOME", "CP:WIKIPATHWAYS",
                "CP:BIOCARTA", "CP:PID")
  c2_list <- lapply(cp_subs, function(sub) {
    msigdbr::msigdbr(species = species, collection = "C2", subcollection = sub)
  })
  c2_all <- do.call(rbind, c2_list, quote = TRUE)

  disease_pat <- paste0("DISEASE|CANCER|TUMOR|CARCINOMA|LEUKEMIA|LYMPHOMA|",
                         "MELANOMA|GLIOMA|HEPATITIS|HIV|INFECTION|VIRAL|",
                         "BACTERIAL|PARASIT")
  c2_all <- c2_all[!grepl(disease_pat, c2_all$gs_name, ignore.case = TRUE), ]

  all_sets <- rbind(
    hallmark[, c("gs_name", "gene_symbol")],
    c2_all[, c("gs_name", "gene_symbol")],
    gobp[, c("gs_name", "gene_symbol")]
  )

  pw_list <- split(all_sets$gene_symbol, all_sets$gs_name)
  pw_list <- lapply(pw_list, unique)

  sizes <- vapply(pw_list, length, integer(1))
  pw_list <- pw_list[sizes >= min_size & sizes <= max_size]

  db <- classify_database(names(pw_list))
  db_tab <- table(db)
  message(sprintf("Full pathway collection: %d sets, size %d-%d",
                  length(pw_list), min_size, max_size))
  message(sprintf("  Databases: %s",
                  paste(sprintf("%s=%d", names(db_tab), db_tab), collapse = ", ")))
  pw_list
}

pw_collection <- build_full_collection(min_size = 10, max_size = 500)

set.seed(42)
pw_by_db <- split(names(pw_collection), classify_database(names(pw_collection)))

fgsea_perdb_all <- list()

for (ctr in CONTRASTS) {
  stats <- dep_df[[paste0("t_", ctr)]]
  names(stats) <- dep_df$gene
  stats <- sort(stats[!is.na(stats) & is.finite(stats)], decreasing = TRUE)

  for (db_name in names(pw_by_db)) {
    db_pw <- pw_collection[pw_by_db[[db_name]]]
    if (length(db_pw) < 5) next
    raw_db <- fgsea::fgseaMultilevel(
      pathways    = db_pw,
      stats       = stats,
      minSize     = 10,
      maxSize     = 500,
      nPermSimple = 10000,
      eps         = 0
    )
    raw_db <- as.data.frame(raw_db)
    raw_db$database <- db_name
    raw_db$contrast <- ctr
    fgsea_perdb_all[[paste0(ctr, "_", db_name)]] <- as_tibble(raw_db)
  }
}

fgsea_combined <- bind_rows(fgsea_perdb_all)

fgsea_export <- fgsea_combined |>
  mutate(leadingEdge = sapply(leadingEdge, paste, collapse = ";")) |>
  arrange(database, contrast, padj)
write_csv(fgsea_export, file.path(DAT, "06_panel_H_fgsea_results.csv"))

fgsea_hbp_cache <- fgsea_export %>% filter(database %in% c("Hallmark", "GO:BP"))
cache_dir <- "04_Figures/F2/c_data/shared"
dir.create(cache_dir, recursive = TRUE, showWarnings = FALSE)
write_csv(fgsea_hbp_cache, file.path(cache_dir, "fgsea_tstat_all_v2.csv"))

DISPLAY_DBS <- c("Hallmark", "KEGG", "Reactome", "WikiPathways", "GO:BP")

db_totals <- fgsea_combined |>
  filter(database %in% DISPLAY_DBS) |>
  distinct(pathway, database) |>
  count(database, name = "n_total")

count_df <- fgsea_combined |>
  filter(!is.na(padj), padj < 0.05, database %in% DISPLAY_DBS) |>
  group_by(contrast, database) |>
  summarise(
    Up   = sum(NES > 0),
    Down = sum(NES < 0),
    .groups = "drop"
  ) |>
  tidyr::pivot_longer(cols = c(Up, Down), names_to = "direction",
                      values_to = "count") |>
  left_join(db_totals, by = "database") |>
  mutate(fraction = count / n_total)

nonempty_dbs <- count_df |>
  group_by(database) |> filter(sum(count) > 0) |> pull(database) |> unique()
count_df <- count_df |> filter(database %in% nonempty_dbs)

db_labels <- db_totals |>
  filter(database %in% nonempty_dbs) |>
  mutate(label = sprintf("%s\n(n=%d)", database, n_total)) |>
  tibble::deframe()
db_labels <- setNames(
  sprintf("%s\n(n=%d)", db_totals$database, db_totals$n_total),
  db_totals$database
)[nonempty_dbs]

count_df$contrast  <- factor(count_df$contrast, levels = CONTRASTS)
count_df$database  <- factor(count_df$database, levels = intersect(DISPLAY_DBS, nonempty_dbs))
count_df$direction <- factor(count_df$direction, levels = c("Up", "Down"))

n_facets <- length(levels(count_df$database))
PH_H <- max(80, n_facets * 55)

lbl_sz <- scale_text(BASE_COUNT, PH_W)

pH <- ggplot(count_df, aes(x = contrast, y = fraction * 100, fill = direction)) +
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
  facet_grid(database ~ ., scales = "free_y",
             labeller = as_labeller(db_labels)) +
  scale_x_discrete(labels = CTR_SHORT) +
  scale_fill_manual(values = DIR_COLORS) +
  labs(title = "Pathway Enrichment",
       subtitle = "fGSEA (padj < 0.05)\nper-database BH correction",
       x = NULL, y = "% of database significant",
       tag = "H") +
  FIG_THEME +
  theme(axis.text.x    = element_text(angle = 35, hjust = 1,
                                     size = FIG_AXIS_TEXT - 1.5),
        legend.position = "none",
        strip.text.y   = element_text(size = FIG_STRIP_SIZE, face = "bold", angle = 0))

nes_summary <- fgsea_combined |>
  filter(!is.na(padj), padj < 0.05) |>
  group_by(contrast, database) |>
  summarise(
    n_sig = n(), n_up = sum(NES > 0), n_down = sum(NES < 0),
    median_NES = median(NES), mean_NES = mean(NES), sd_NES = sd(NES),
    min_padj = min(padj), median_padj = median(padj),
    .groups = "drop"
  )
write.csv(nes_summary, file.path(DAT, "audit_panel_H_nes_summary.csv"),
          row.names = FALSE)

ggsave(file.path(RPT, "panel_H_fgsea.pdf"), pH,
       width = PH_W, height = PH_H, units = "mm", device = pdf_device)
ggsave(file.path(RPT, "panel_H_fgsea.png"), pH,
       width = PH_W, height = PH_H, units = "mm", dpi = 300)
