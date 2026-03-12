################################################################################
#   Figure 1 — Panel H Options: Alternative Pathway Visualizations
#   Generates side-by-side comparison of size-bias-corrected approaches.
#   Outputs: PDF + PNG in b_reports/panel_H_options/
################################################################################

setwd(rprojroot::find_rstudio_root_file())
source("04_Figures_v2/shared/style.R")
source("04_Figures_v2/shared/pathway_utils.R")

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(tibble)
  library(tidyr)
  library(ggplot2)
  library(fgsea)
  library(msigdbr)
  library(stringr)
})

# ── Paths ────────────────────────────────────────────────────────────────────

DEP_FILE <- "03_DEP/c_data/03_combined_results.csv"
RPT      <- "04_Figures_v2/F1/b_reports/panel_H_options"
dir.create(RPT, recursive = TRUE, showWarnings = FALSE)

# ── Constants ────────────────────────────────────────────────────────────────

CONTRASTS <- c("Aging", "Training_Young", "Training_Old", "Interaction")
PI_THRESH <- 0.05

# ── Load data ────────────────────────────────────────────────────────────────

dep_df <- read_csv(DEP_FILE, show_col_types = FALSE)

# ── Build expanded pathway collection (all CP subcollections) ────────────────
# build_pathway_collection() only fetches C2:CP (19 generic sets), missing
# KEGG, Reactome, WikiPathways, BioCarta, PID. Build locally with all subs.

build_full_collection <- function(species = "Homo sapiens",
                                   min_size = 10, max_size = 500) {
  requireNamespace("msigdbr", quietly = TRUE)

  hallmark <- msigdbr::msigdbr(species = species, collection = "H")
  gobp     <- msigdbr::msigdbr(species = species, collection = "C5",
                                 subcollection = "GO:BP")

  # All C2:CP subcollections
  cp_subs <- c("CP", "CP:KEGG_MEDICUS", "CP:REACTOME", "CP:WIKIPATHWAYS",
                "CP:BIOCARTA", "CP:PID")
  c2_list <- lapply(cp_subs, function(sub) {
    msigdbr::msigdbr(species = species, collection = "C2", subcollection = sub)
  })
  c2_all <- do.call(rbind, c2_list, quote = TRUE)

  # Filter disease/cancer terms from C2
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
  message(sprintf("Full pathway collection: %d sets, size %d-%d", length(pw_list), min_size, max_size))
  message(sprintf("  Databases: %s", paste(sprintf("%s=%d", names(db_tab), db_tab), collapse = ", ")))
  pw_list
}

pw_collection <- build_full_collection(min_size = 10, max_size = 500)

# ── Run fGSEA ────────────────────────────────────────────────────────────────
# Two versions:
#   1. Deduped (cross-database Jaccard) — for Options A-D
#   2. Per-database (separate BH per database) — for heatmap (Options E/E2)
#      This ensures KEGG, Reactome, WikiPathways etc. aren't drowned by GO:BP.

set.seed(42)
fgsea_dedup_all <- list()
fgsea_perdb_all <- list()

# Split pathway collection by database
pw_by_db <- split(names(pw_collection), classify_database(names(pw_collection)))

for (ctr in CONTRASTS) {
  stats <- dep_df[[paste0("t_", ctr)]]
  names(stats) <- dep_df$gene
  stats <- sort(stats[!is.na(stats) & is.finite(stats)], decreasing = TRUE)

  # Deduped version (for Options A-D)
  dedup <- run_fgsea_deduplicated(
    ranks          = stats,
    pathways       = pw_collection,
    jaccard_cutoff = 0.5,
    nperm          = 10000,
    min_size       = 10,
    max_size       = 500
  )
  dedup$contrast <- ctr
  fgsea_dedup_all[[ctr]] <- dedup

  # Per-database fGSEA (separate BH per database)
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

  n_dedup <- sum(!is.na(dedup$padj) & dedup$padj < 0.05)
  cat(sprintf("  fgsea done: %s (dedup: %d sig)\n", ctr, n_dedup))
}

fgsea_combined <- bind_rows(fgsea_dedup_all)
fgsea_perdb    <- bind_rows(fgsea_perdb_all)

# Report per-database significance
perdb_sig <- fgsea_perdb |>
  filter(!is.na(padj), padj < 0.05) |>
  distinct(pathway, database) |>
  count(database, name = "n_sig")
cat("\nPer-database significant pathways (separate BH):\n")
print(as.data.frame(perdb_sig))

all_genes <- dep_df$gene

# ── PDF device ───────────────────────────────────────────────────────────────

pdf_device <- get_pdf_device()

# ── Contrast background shading helper ───────────────────────────────────────

contrast_bg <- function(n_contrasts = 4) {
  list(
    annotate("rect", xmin = 0.5, xmax = 1.5, ymin = -Inf, ymax = Inf,
             fill = CONTRAST_COLORS["Aging"], alpha = 0.20,
             color = "grey70", linewidth = 0.2),
    annotate("rect", xmin = 1.5, xmax = 2.5, ymin = -Inf, ymax = Inf,
             fill = CONTRAST_COLORS["Training_Young"], alpha = 0.20,
             color = "grey70", linewidth = 0.2),
    annotate("rect", xmin = 2.5, xmax = 3.5, ymin = -Inf, ymax = Inf,
             fill = CONTRAST_COLORS["Training_Old"], alpha = 0.20,
             color = "grey70", linewidth = 0.2),
    annotate("rect", xmin = 3.5, xmax = 4.5, ymin = -Inf, ymax = Inf,
             fill = CONTRAST_COLORS["Interaction"], alpha = 0.20,
             color = "grey70", linewidth = 0.2)
  )
}


# ── Expanded 20-category classifier (local, for Option D) ──────────────────
# Broader keyword rules to reduce "Other" bucket.

EXPANDED_CAT_ORDER <- c(
  "Muscle & Contractile", "Cytoskeleton & Motility", "ECM & Adhesion",
  "Lipid Metabolism", "Carbohydrate Metabolism", "Amino Acid Metabolism",
  "Mitochondria & Energy", "Redox & ROS",
  "Protein Homeostasis", "Protein Modification",
  "Vesicle & Transport", "Gene Expression & RNA",
  "Ribosome & Translation", "Signaling",
  "Immune & Inflammation", "Cell Death & Survival",
  "DNA & Cell Cycle", "Morphogenesis & Development",
  "Circulatory & Hemostasis", "Ion & Calcium",
  "Other"
)

EXPANDED_CAT_COLORS <- c(
  "Muscle & Contractile"      = "#E57373",
  "Cytoskeleton & Motility"   = "#FF8A65",
  "ECM & Adhesion"            = "#FFB74D",
  "Lipid Metabolism"           = "#FFD54F",
  "Carbohydrate Metabolism"    = "#DCE775",
  "Amino Acid Metabolism"      = "#AED581",
  "Mitochondria & Energy"      = "#4DB6AC",
  "Redox & ROS"                = "#26A69A",
  "Protein Homeostasis"        = "#4FC3F7",
  "Protein Modification"       = "#29B6F6",
  "Vesicle & Transport"        = "#7986CB",
  "Gene Expression & RNA"      = "#9575CD",
  "Ribosome & Translation"     = "#BA68C8",
  "Signaling"                  = "#CE93D8",
  "Immune & Inflammation"      = "#A1887F",
  "Cell Death & Survival"      = "#EF9A9A",
  "DNA & Cell Cycle"           = "#90A4AE",
  "Morphogenesis & Development" = "#B0BEC5",
  "Circulatory & Hemostasis"   = "#F48FB1",
  "Ion & Calcium"              = "#80CBC4",
  "Other"                      = "#D0D0D0"
)

classify_expanded <- function(ids) {
  rules <- list(
    "Muscle & Contractile"      = "MYOGEN|MYOFIBRIL|SARCOMERE|MUSCLE_|CONTRACTILE|ACTOMYOSIN|MYOSIN|I_BAND|STRIATED|SKELETAL_MUSCLE",
    "Cytoskeleton & Motility"   = "CYTOSKELET|ACTIN_BIND|ACTIN_FILAMENT|STRUCTURAL_MOLECULE|MOTIL|SUPRAMOLECUL|MICROTUB|INTERMEDIATE_FILAMENT|TUBULIN|SPINDLE",
    "ECM & Adhesion"            = "EXTRACELLULAR_MATRIX|COLLAGEN|BASEMENT_MEMBRANE|ADHESION|APICAL_JUNCTION|EMT|ENCAPSULATING|INTEGRIN|LAMININ|MESENCHYM|EPITHELIAL_MESENCHYMAL",
    "Lipid Metabolism"           = "LIPID|FATTY_ACID|ADIPOGEN|CHOLESTEROL|STEROID|SPHINGO|BILE_ACID|PEROXISOME|PHOSPHOLIPID|TRIGLYCERIDE",
    "Carbohydrate Metabolism"    = "GLYCOLY|GLUCONEO|GLYCOGEN|PENTOSE|PYRUVATE|GLUCOSE|FRUCTOSE|GALACTOSE|HEXOSE",
    "Amino Acid Metabolism"      = "AMINO_ACID|BRANCHED_CHAIN|UREA_CYCLE|GLUTATHIONE|TRYPTOPHAN|GLUTAMATE|ARGININE|METHIONINE|CYSTEINE",
    "Mitochondria & Energy"      = "MITOCHOND|OXIDATIVE_PHOSPH|ELECTRON_TRANSFER|RESPIRATORY|OXIDOREDUCT|TCA_CYCLE|CITRATE|NADH|ATP_SYNTH|PROTON_MOTIVE",
    "Redox & ROS"                = "REACTIVE_OXYGEN|PEROXIDE|ANTIOXIDANT|DETOXIF|XENOBIOT|REDOX|GLUTATHIONE_TRANSFER|OXIDANT",
    "Protein Homeostasis"        = "PROTEASOM|UBIQUITIN|AUTOPHAGY|MTORC1|PROTEIN_FOLD|UNFOLDED|CHAPERONE|HEAT_SHOCK|PROTEIN_STABIL|PROTEIN_CATAB|ENDOPLASMIC_RETICULUM.*STRESS",
    "Protein Modification"       = "PHOSPHORYL|GLYCOSYL|ACETYL|METHYLAT|SUMOYL|NEDDYL|UBIQUITIN_LIKE|PROTEIN_PROCESS|POST_TRANSLAT",
    "Vesicle & Transport"        = "VESICLE|ENDOCYT|EXOCYT|SECRETI|GOLGI|LYSOSOM|ENDOSOM|PHAGOCYT|INTRACELLULAR_TRANSPORT|PROTEIN_LOCALIZATION|MEMBRANE_FUSION",
    "Gene Expression & RNA"      = "TRANSCRIPT|SPLICEOSOM|MYC_TARGET|E2F_TARGET|RNA_PROCESS|RNA_LOCALI|MRNA|CHROMATIN|HISTONE|EPIGENET|GENE_SILENC|RNA_SPLICING|POLY_A",
    "Ribosome & Translation"     = "RIBOSOM|TRANSLAT(?!.*SIGNAL)|TRNA|PEPTIDE_BIOSYN",
    "Signaling"                  = "(?<!ELECTRO.{0,5})SIGNAL(?!.*RIBOSOM)|KINASE|ANDROGEN|ESTROGEN|NOTCH|WNT|HEDGEHOG|TGF_BETA|KRAS|PI3K|MTOR(?!C1)|RECEPTOR|GROWTH_FACTOR|MAPK|JAK_STAT|HORMONE|NF.KB|TNFA",
    "Immune & Inflammation"      = "IMMUN|INFLAMMA|INTERFERON|IL[0-9]|COMPLEMENT|ALLOGRAFT|LEUKOCYTE|LYMPHOCYTE|CYTOKINE|INNATE|ADAPTIVE",
    "Cell Death & Survival"      = "APOPTOSIS|PROGRAMMED_CELL_DEATH|NECROPTOSIS|FERROPTOSIS|CELL_DEATH|SURVIVAL|P53_PATHWAY|UV_RESPONSE",
    "DNA & Cell Cycle"           = "DNA_REPAIR|CELL_CYCLE|MITOTIC|G2M|CELL_DIVIS|CHECKPOINT|REPLICAT|CHROMOSOME|TELOMERE",
    "Morphogenesis & Development" = "MORPHOGEN|DIFFERENTIAT(?!.*MUSCLE)|EMBRYO|ORGANOGEN|ORGAN_DEVELOP|TISSUE_DEVELOP|LIMB|PATTERN|CARDIAC.*DEVELOP|ANGIOGEN|VASCULAT|BLOOD_VESSEL|NEUROGEN|AXON|SYNAP|DENDRIT",
    "Circulatory & Hemostasis"   = "COAGULAT|HEMOSTAS|PLATELET|FIBRIN|HEME|ERYTHROCYTE|BLOOD_COAG|SPERMATOGEN",
    "Ion & Calcium"              = "ION_TRANSPORT|CALCIUM|POTASSIUM|SODIUM|CHANNEL|CATION|ANION|METAL_ION"
  )
  vapply(toupper(ids), function(id) {
    for (cat in names(rules)) {
      if (grepl(rules[[cat]], id, perl = TRUE)) return(cat)
    }
    "Other"
  }, character(1), USE.NAMES = FALSE)
}


################################################################################
#   Option A: ORA Enrichment Bars (Field Standard)
################################################################################

cat("Building Option A: ORA Enrichment Bars\n")

ora_top_all <- list()
for (ctr in CONTRASTS) {
  pi_col <- paste0("pi_score_", ctr)
  sig_genes <- dep_df$gene[!is.na(dep_df[[pi_col]]) & dep_df[[pi_col]] < PI_THRESH]

  if (length(sig_genes) < 5) {
    cat(sprintf("  ORA skipped for %s: only %d DEPs\n", ctr, length(sig_genes)))
    next
  }

  ora_res <- run_ora_deduplicated(
    genes          = sig_genes,
    universe       = all_genes,
    pathways       = pw_collection,
    jaccard_cutoff = 0.5,
    min_size       = 10,
    max_size       = 500,
    padj_cutoff    = 0.05
  )

  if (nrow(ora_res) > 0) {
    ora_top <- ora_res |>
      arrange(padj) |>
      head(10) |>
      mutate(
        contrast = ctr,
        clean_name = clean_pathway_name(pathway, max_chars = 40),
        neg_log10_padj = -log10(padj)
      )
    ora_top_all[[ctr]] <- ora_top
  }
  cat(sprintf("  ORA done: %s (%d sig terms)\n", ctr, nrow(ora_res)))
}

ora_plot_df <- bind_rows(ora_top_all)

if (nrow(ora_plot_df) > 0) {
  ora_plot_df$contrast <- factor(ora_plot_df$contrast, levels = CONTRASTS)

  pH_A <- ggplot(ora_plot_df,
                 aes(x = neg_log10_padj,
                     y = reorder_within(clean_name, neg_log10_padj, contrast),
                     fill = database)) +
    geom_col(color = "black", linewidth = 0.2, width = 0.7) +
    facet_wrap(~ contrast, scales = "free_y", ncol = 1,
               labeller = as_labeller(CTR_FACET)) +
    scale_y_reordered() +
    scale_fill_manual(values = DB_COLORS, name = "Database") +
    labs(title = "Option A: ORA Enrichment",
         subtitle = "Top 10 enriched terms per contrast\n(Pi-score DEPs, Jaccard dedup)",
         x = expression(-log[10](p[adj])), y = NULL) +
    FIG_THEME +
    theme(axis.text.y     = element_text(size = 6.5),
          strip.text       = element_text(face = "bold", size = FIG_STRIP_SIZE),
          legend.position  = "bottom",
          legend.key.size  = unit(3, "mm"),
          legend.text      = element_text(size = 7))

  n_terms <- ora_plot_df |> group_by(contrast) |> summarise(n = n()) |> pull(n) |> sum()
  PHA_H <- max(120, n_terms * 4 + 40)

  ggsave(file.path(RPT, "option_A_ora_bars.pdf"), pH_A,
         width = 140, height = PHA_H, units = "mm", device = pdf_device)
  ggsave(file.path(RPT, "option_A_ora_bars.png"), pH_A,
         width = 140, height = PHA_H, units = "mm", dpi = 300)
  cat("Option A saved\n")
} else {
  cat("Option A: no significant ORA results\n")
}


################################################################################
#   Option B: Hallmark NES Heatmap
################################################################################

cat("Building Option B: Hallmark NES Heatmap\n")

hallmark_df <- fgsea_combined |>
  filter(database == "Hallmark")

hallmark_wide <- hallmark_df |>
  select(pathway, contrast, NES, padj) |>
  mutate(clean_name = clean_pathway_name(pathway, max_chars = 35))

heatmap_df <- hallmark_wide |>
  mutate(
    significant = !is.na(padj) & padj < 0.05,
    NES_display = ifelse(significant, NES, NA_real_),
    contrast = factor(contrast, levels = CONTRASTS)
  )

pathway_order <- heatmap_df |>
  group_by(clean_name) |>
  summarise(mean_abs_nes = mean(abs(NES), na.rm = TRUE), .groups = "drop") |>
  arrange(mean_abs_nes) |>
  pull(clean_name)

heatmap_df$clean_name <- factor(heatmap_df$clean_name, levels = pathway_order)
nes_lim <- max(abs(heatmap_df$NES), na.rm = TRUE)

pH_B <- ggplot(heatmap_df, aes(x = contrast, y = clean_name)) +
  geom_tile(data = heatmap_df |> filter(!significant),
            fill = "grey90", color = "white", linewidth = 0.3) +
  geom_tile(data = heatmap_df |> filter(significant),
            aes(fill = NES), color = "white", linewidth = 0.3) +
  geom_text(data = heatmap_df |> filter(significant),
            aes(label = sig_stars(padj)),
            size = 2.5, color = "black", vjust = 0.8) +
  scale_fill_gradient2(
    low = DIR_COLORS["Down"], mid = "white", high = DIR_COLORS["Up"],
    midpoint = 0, limits = c(-nes_lim, nes_lim),
    name = "NES", na.value = "grey90"
  ) +
  scale_x_discrete(labels = CTR_SHORT) +
  labs(title = "Option B: Hallmark NES Heatmap",
       subtitle = "50 MSigDB Hallmark gene sets\n(grey = padj >= 0.05)",
       x = NULL, y = NULL) +
  FIG_THEME +
  theme(axis.text.y      = element_text(size = 5.5),
        axis.text.x      = element_text(angle = 35, hjust = 1, size = 8),
        legend.position   = "right",
        legend.key.height = unit(12, "mm"),
        legend.key.width  = unit(3, "mm"),
        panel.grid        = element_blank())

n_hallmark <- length(unique(heatmap_df$clean_name))
PHB_H <- max(120, n_hallmark * 4 + 30)

ggsave(file.path(RPT, "option_B_hallmark_heatmap.pdf"), pH_B,
       width = 130, height = PHB_H, units = "mm", device = pdf_device)
ggsave(file.path(RPT, "option_B_hallmark_heatmap.png"), pH_B,
       width = 130, height = PHB_H, units = "mm", dpi = 300)
cat("Option B saved\n")


################################################################################
#   Option C: Hallmark-Only Bar Chart
################################################################################

cat("Building Option C: Hallmark-Only Bar Chart\n")

count_hallmark <- fgsea_combined |>
  filter(!is.na(padj), padj < 0.05, database == "Hallmark") |>
  group_by(contrast) |>
  summarise(
    Up   = sum(NES > 0),
    Down = sum(NES < 0),
    .groups = "drop"
  ) |>
  pivot_longer(cols = c(Up, Down), names_to = "direction", values_to = "count")

count_hallmark$contrast  <- factor(count_hallmark$contrast, levels = CONTRASTS)
count_hallmark$direction <- factor(count_hallmark$direction, levels = c("Up", "Down"))

PHC_W <- 100
PHC_H <- 80
lbl_sz_c <- scale_text(BASE_COUNT, PHC_W)

pH_C <- ggplot(count_hallmark, aes(x = contrast, y = count, fill = direction)) +
  contrast_bg() +
  geom_col(position = position_dodge(width = 0.7), width = 0.6,
           color = "black", linewidth = 0.3) +
  geom_text(data = count_hallmark |> filter(count > 0) |> mutate(mid = count / 2),
            aes(y = mid, label = count),
            position = position_dodge(width = 0.7),
            vjust = 0.5, hjust = 0.5, size = lbl_sz_c,
            color = "white", fontface = "bold") +
  scale_x_discrete(labels = CTR_SHORT) +
  scale_fill_manual(values = DIR_COLORS) +
  labs(title = "Option C: Hallmark-Only Bars",
       subtitle = "50 curated gene sets\n(Liberzon et al. 2015)",
       x = NULL, y = "Significant pathways") +
  FIG_THEME +
  theme(axis.text.x    = element_text(angle = 35, hjust = 1, size = 7),
        legend.position = "none")

ggsave(file.path(RPT, "option_C_hallmark_bars.pdf"), pH_C,
       width = PHC_W, height = PHC_H, units = "mm", device = pdf_device)
ggsave(file.path(RPT, "option_C_hallmark_bars.png"), pH_C,
       width = PHC_W, height = PHC_H, units = "mm", dpi = 300)
cat("Option C saved\n")


################################################################################
#   Option D: Expanded 20-Category Bars (fGSEA-dedup approach)
#   Uses significant fGSEA pathways (Jaccard-deduped) → expanded keyword
#   classifier with 20 biological categories to minimize "Other".
################################################################################

cat("Building Option D: Expanded Category Bars (fGSEA-dedup)\n")

consol_fgsea <- fgsea_combined |>
  filter(!is.na(padj), padj < 0.05) |>
  mutate(expanded_cat = classify_expanded(pathway))

# Report classification coverage
n_other <- sum(consol_fgsea$expanded_cat == "Other")
n_total <- nrow(consol_fgsea)
cat(sprintf("  Classification: %d/%d assigned (%.1f%%), %d Other\n",
            n_total - n_other, n_total, 100 * (1 - n_other / n_total), n_other))

# Count pathways per category per contrast
consol_count <- consol_fgsea |>
  count(contrast, expanded_cat, name = "n")

consol_count$contrast    <- factor(consol_count$contrast, levels = CONTRASTS)
consol_count$expanded_cat <- factor(consol_count$expanded_cat, levels = EXPANDED_CAT_ORDER)

# Drop empty categories
active_cats <- consol_count |>
  group_by(expanded_cat) |> summarise(total = sum(n)) |>
  filter(total > 0) |> pull(expanded_cat)
consol_count <- consol_count |> filter(expanded_cat %in% active_cats)

PHD_W <- 130
PHD_H <- 105

pH_D <- ggplot(consol_count, aes(x = contrast, y = n, fill = expanded_cat)) +
  contrast_bg() +
  geom_col(color = "black", linewidth = 0.2, width = 0.7) +
  scale_x_discrete(labels = CTR_SHORT) +
  scale_fill_manual(values = EXPANDED_CAT_COLORS, name = "Category",
                    drop = TRUE) +
  labs(title = "Option D: Biological Categories",
       subtitle = "fGSEA significant pathways\nclassified into 20 categories",
       x = NULL, y = "Significant pathways") +
  FIG_THEME +
  theme(axis.text.x    = element_text(angle = 35, hjust = 1, size = 7),
        legend.position = "right",
        legend.text     = element_text(size = 6),
        legend.title    = element_text(size = 8, face = "bold"),
        legend.key.size = unit(2.5, "mm"))

ggsave(file.path(RPT, "option_D_consolidated_bars.pdf"), pH_D,
       width = PHD_W, height = PHD_H, units = "mm", device = pdf_device)
ggsave(file.path(RPT, "option_D_consolidated_bars.png"), pH_D,
       width = PHD_W, height = PHD_H, units = "mm", dpi = 300)
cat("Option D saved\n")


################################################################################
#   Option E: Full NES Heatmap (all databases, per-database BH)
#   All pathways significant in >= 1 contrast, ordered by cross-contrast density.
#   Uses per-database fGSEA so KEGG/Reactome/WikiPathways aren't drowned by GO:BP.
################################################################################

cat("Building Option E: Full NES Heatmap (all databases)\n")

# Identify pathways significant in at least 1 contrast (per-database BH)
sig_pathways_perdb <- fgsea_perdb |>
  filter(!is.na(padj), padj < 0.05) |>
  pull(pathway) |> unique()

sig_db_counts <- fgsea_perdb |>
  filter(pathway %in% sig_pathways_perdb) |>
  distinct(pathway, database) |>
  count(database, name = "n_pathways")
cat("  Significant pathways by database (per-db BH):\n")
print(as.data.frame(sig_db_counts))

cat(sprintf("  %d unique pathways significant in >= 1 contrast\n",
            length(sig_pathways_perdb)))

# Build full matrix for those pathways
heat_full <- fgsea_perdb |>
  filter(pathway %in% sig_pathways_perdb) |>
  select(pathway, contrast, NES, padj, database) |>
  mutate(
    clean_name  = clean_pathway_name(pathway, max_chars = 50),
    significant = !is.na(padj) & padj < 0.05,
    contrast    = factor(contrast, levels = CONTRASTS)
  )

# Deduplicate clean_names: use pathway-level unique labels
pw_label_map <- heat_full |>
  distinct(pathway, clean_name, database) |>
  group_by(clean_name) |>
  mutate(n_dup = n()) |>
  ungroup() |>
  mutate(clean_label = ifelse(
    n_dup > 1,
    paste0(clean_name, " [", database, "]"),
    clean_name
  )) |>
  select(pathway, clean_label)

# Final pass: ensure truly unique (make.unique handles any remaining dupes)
pw_label_map$clean_label <- make.unique(pw_label_map$clean_label, sep = " ")

heat_full <- heat_full |>
  left_join(pw_label_map, by = "pathway") |>
  mutate(clean_name = clean_label) |>
  select(-clean_label)

# Ordering: n_sig (density) then mean |NES|
pathway_order_full <- heat_full |>
  group_by(pathway, clean_name, database) |>
  summarise(
    n_sig        = sum(significant),
    mean_abs_nes = mean(abs(NES), na.rm = TRUE),
    .groups = "drop"
  ) |>
  arrange(n_sig, mean_abs_nes)

heat_full$clean_name <- factor(heat_full$clean_name,
                                levels = pathway_order_full$clean_name)

# Database strip
db_strip <- pathway_order_full |>
  select(clean_name, database) |>
  mutate(clean_name = factor(clean_name, levels = pathway_order_full$clean_name))

nes_lim_full <- max(abs(heat_full$NES), na.rm = TRUE)

# Main heatmap
pH_E <- ggplot(heat_full, aes(x = contrast, y = clean_name)) +
  geom_tile(data = heat_full |> filter(!significant),
            fill = "grey93", color = "white", linewidth = 0.1) +
  geom_tile(data = heat_full |> filter(significant),
            aes(fill = NES), color = "white", linewidth = 0.1) +
  scale_fill_gradient2(
    low = DIR_COLORS["Down"], mid = "white", high = DIR_COLORS["Up"],
    midpoint = 0, limits = c(-nes_lim_full, nes_lim_full),
    name = "NES", na.value = "grey93"
  ) +
  scale_x_discrete(labels = CTR_SHORT) +
  labs(title = "Option E: Full Pathway NES Heatmap",
       subtitle = sprintf("%d pathways (sig in >= 1 contrast)\nordered by cross-contrast density",
                          length(sig_pathways_perdb)),
       x = NULL, y = NULL) +
  FIG_THEME +
  theme(axis.text.y       = element_text(size = 2.5),
        axis.text.x       = element_text(angle = 35, hjust = 1, size = 8),
        legend.position    = "right",
        legend.key.height  = unit(15, "mm"),
        legend.key.width   = unit(3, "mm"),
        panel.grid         = element_blank())

n_pw <- length(sig_pathways_perdb)
PHE_H <- max(200, n_pw * 2.5 + 40)
PHE_W <- 170

ggsave(file.path(RPT, "option_E_full_heatmap.pdf"), pH_E,
       width = PHE_W, height = PHE_H, units = "mm", device = pdf_device,
       limitsize = FALSE)
ggsave(file.path(RPT, "option_E_full_heatmap.png"), pH_E,
       width = PHE_W, height = PHE_H, units = "mm", dpi = 300,
       limitsize = FALSE)
cat("Option E saved\n")


################################################################################
#   Option E2: Full Heatmap with database color strip (patchwork)
################################################################################

cat("Building Option E2: Full NES Heatmap + database strip\n")

p_db_strip <- ggplot(db_strip, aes(x = 1, y = clean_name, fill = database)) +
  geom_tile(color = "white", linewidth = 0.1) +
  scale_fill_manual(values = DB_COLORS, name = "Database") +
  scale_x_continuous(breaks = NULL) +
  labs(x = NULL, y = NULL) +
  theme_void() +
  theme(axis.text.y  = element_blank(),
        legend.position = "bottom",
        legend.key.size = unit(3, "mm"),
        legend.text     = element_text(size = 6.5))

if (requireNamespace("patchwork", quietly = TRUE)) {
  library(patchwork)

  pH_E2 <- (p_db_strip + pH_E +
    plot_layout(widths = c(1, 20))) +
    plot_annotation(
      title    = "Option E2: Full Pathway NES Heatmap + Database",
      subtitle = sprintf("%d pathways ordered by cross-contrast significance density",
                         length(sig_pathways_perdb)),
      theme = theme(
        plot.title    = element_text(face = "bold", size = 12),
        plot.subtitle = element_text(face = "bold.italic", size = 9, color = "grey30")
      )
    )

  PHE2_H <- max(200, n_pw * 2.5 + 50)
  PHE2_W <- 190

  ggsave(file.path(RPT, "option_E2_full_heatmap_db.pdf"), pH_E2,
         width = PHE2_W, height = PHE2_H, units = "mm", device = pdf_device,
         limitsize = FALSE)
  ggsave(file.path(RPT, "option_E2_full_heatmap_db.png"), pH_E2,
         width = PHE2_W, height = PHE2_H, units = "mm", dpi = 300,
         limitsize = FALSE)
  cat("Option E2 saved (with database strip)\n")
} else {
  cat("Option E2 skipped: patchwork not installed\n")
}

################################################################################
#   Option F: Original Method with All Databases (Full Collection)
#   Same as current Panel H bar chart (up/down counts per database) but using
#   the full collection (all C2:CP subcollections, not just generic CP).
################################################################################

cat("Building Option F: Original Method (All Databases)\n")

DISPLAY_DBS <- c("Hallmark", "KEGG", "Reactome", "WikiPathways", "GO:BP",
                  "BioCarta", "PID")

count_df_f <- fgsea_combined |>
  filter(!is.na(padj), padj < 0.05, database %in% DISPLAY_DBS) |>
  group_by(contrast, database) |>
  summarise(
    Up   = sum(NES > 0),
    Down = sum(NES < 0),
    .groups = "drop"
  ) |>
  pivot_longer(cols = c(Up, Down), names_to = "direction", values_to = "count")

# Drop databases with no significant pathways
nonempty_dbs <- count_df_f |>
  group_by(database) |> filter(sum(count) > 0) |> pull(database) |> unique()
count_df_f <- count_df_f |> filter(database %in% nonempty_dbs)

count_df_f$contrast  <- factor(count_df_f$contrast, levels = CONTRASTS)
count_df_f$database  <- factor(count_df_f$database,
                                levels = intersect(DISPLAY_DBS, nonempty_dbs))
count_df_f$direction <- factor(count_df_f$direction, levels = c("Up", "Down"))

n_facets_f <- length(levels(count_df_f$database))
PHF_H <- max(80, n_facets_f * 55)
PHF_W <- 100
lbl_sz_f <- scale_text(BASE_COUNT, PHF_W)

pH_F <- ggplot(count_df_f, aes(x = contrast, y = count, fill = direction)) +
  contrast_bg() +
  geom_col(position = position_dodge(width = 0.7), width = 0.6,
           color = "black", linewidth = 0.3) +
  geom_text(data = \(d) d |> filter(count > 0) |> mutate(mid = count / 2),
            aes(y = mid, label = count),
            position = position_dodge(width = 0.7),
            vjust = 0.5, hjust = 0.5, size = lbl_sz_f,
            color = "white", fontface = "bold") +
  facet_grid(database ~ ., scales = "free_y") +
  scale_x_discrete(labels = CTR_SHORT) +
  scale_fill_manual(values = DIR_COLORS) +
  labs(title = "Option F: Original (All DBs)",
       subtitle = "fGSEA counts (padj < 0.05)\nJaccard dedup (0.5)",
       x = NULL, y = "Significant pathways") +
  FIG_THEME +
  theme(axis.text.x    = element_text(angle = 35, hjust = 1, size = 7),
        legend.position = "none",
        strip.text.y   = element_text(size = FIG_STRIP_SIZE, face = "bold",
                                       angle = 0))

ggsave(file.path(RPT, "option_F_original_all_dbs.pdf"), pH_F,
       width = PHF_W, height = PHF_H, units = "mm", device = pdf_device)
ggsave(file.path(RPT, "option_F_original_all_dbs.png"), pH_F,
       width = PHF_W, height = PHF_H, units = "mm", dpi = 300)
cat("Option F saved\n")


################################################################################
#   Option G: compareCluster Dotplot (clusterProfiler)
#   Rows = top enriched terms, Columns = contrasts
#   Dot size = GeneRatio, Color = p.adjust
################################################################################

cat("Building Option G: compareCluster Dotplot\n")

suppressPackageStartupMessages({
  library(clusterProfiler)
  library(enrichplot)
})

# Build gene lists from DEPs per contrast
gene_lists <- list()
for (ctr in CONTRASTS) {
  pi_col <- paste0("pi_score_", ctr)
  sig <- dep_df$gene[!is.na(dep_df[[pi_col]]) & dep_df[[pi_col]] < PI_THRESH]
  if (length(sig) > 0) gene_lists[[ctr]] <- sig
}

# Build term-to-gene mapping from our pathway collection
t2g <- data.frame(
  term = rep(names(pw_collection), vapply(pw_collection, length, integer(1))),
  gene = unlist(pw_collection, use.names = FALSE),
  stringsAsFactors = FALSE
)

# Add clean names
t2g$name <- clean_pathway_name(t2g$term, max_chars = 45)

# Run compareCluster
cc_res <- compareCluster(
  geneClusters = gene_lists,
  fun          = "enricher",
  TERM2GENE    = t2g[, c("term", "gene")],
  TERM2NAME    = t2g[, c("term", "name")] |> distinct(),
  universe     = all_genes,
  pvalueCutoff = 0.05,
  qvalueCutoff = 0.2
)

if (!is.null(cc_res) && nrow(as.data.frame(cc_res)) > 0) {
  pH_G <- dotplot(cc_res, showCategory = 15, font.size = 7) +
    scale_color_gradient(low = DIR_COLORS["Up"], high = "grey80",
                         name = "p.adjust") +
    labs(title = "Option G: compareCluster Dotplot",
         subtitle = "Top 15 enriched terms per contrast\n(ORA, Jaccard dedup collection)") +
    FIG_THEME +
    theme(axis.text.y    = element_text(size = 6),
          axis.text.x    = element_text(angle = 35, hjust = 1, size = 8),
          legend.position = "right")

  PHG_H <- 220
  PHG_W <- 200

  ggsave(file.path(RPT, "option_G_comparecluster_dotplot.pdf"), pH_G,
         width = PHG_W, height = PHG_H, units = "mm", device = pdf_device)
  ggsave(file.path(RPT, "option_G_comparecluster_dotplot.png"), pH_G,
         width = PHG_W, height = PHG_H, units = "mm", dpi = 300)
  cat("Option G saved\n")
} else {
  cat("Option G: no significant compareCluster results\n")
}


################################################################################
#   Option H: mitch Multi-Contrast Enrichment
#   Rank-MANOVA on all contrasts simultaneously (Kaspi & Ziemann 2020)
#   Built-in heatmap visualization
################################################################################

cat("Building Option H: mitch Multi-Contrast Enrichment\n")

suppressPackageStartupMessages(library(mitch))

# Build input: data.frame of t-statistics per contrast (genes as rows)
mitch_input <- dep_df |>
  select(gene, starts_with("t_")) |>
  column_to_rownames("gene")
colnames(mitch_input) <- sub("^t_", "", colnames(mitch_input))
mitch_input <- mitch_input[, CONTRASTS]
mitch_input <- mitch_input[complete.cases(mitch_input), ]

# Build gene set list for mitch (list of character vectors)
mitch_sets <- as.list(pw_collection)

# Run mitch
mitch_res <- mitch_calc(
  x         = mitch_input,
  genesets  = mitch_sets,
  priority  = "effect",
  resrows   = nrow(mitch_input),
  cores     = 1
)

# mitch built-in report (HTML) — skip, just extract results
mitch_df <- mitch_res$enrichment_result
mitch_sig <- mitch_df |>
  filter(p.adjustMANOVA < 0.05) |>
  arrange(p.adjustMANOVA)

cat(sprintf("  mitch: %d significant pathways (MANOVA padj < 0.05)\n",
            nrow(mitch_sig)))

# Build heatmap of top mitch pathways (s. columns = effect size per contrast)
n_top_mitch <- min(50, nrow(mitch_sig))
mitch_top <- mitch_sig |> head(n_top_mitch)

mitch_heat <- mitch_top |>
  select(set, starts_with("s.")) |>
  mutate(clean_name = clean_pathway_name(set, max_chars = 45)) |>
  mutate(clean_name = make.unique(clean_name, sep = " "))

# Pivot for ggplot heatmap
mitch_long <- mitch_heat |>
  pivot_longer(cols = starts_with("s."), names_to = "contrast",
               values_to = "effect") |>
  mutate(contrast = sub("^s\\.", "", contrast))

# Add significance from mitch results
mitch_padj <- mitch_top |>
  select(set, p.adjustMANOVA) |>
  rename(manova_padj = p.adjustMANOVA)
mitch_long <- mitch_long |> left_join(mitch_padj, by = "set")

# Order by mean absolute effect
pw_ord_mitch <- mitch_long |>
  group_by(clean_name) |>
  summarise(mean_eff = mean(abs(effect), na.rm = TRUE), .groups = "drop") |>
  arrange(mean_eff) |>
  pull(clean_name)

mitch_long$clean_name <- factor(mitch_long$clean_name, levels = pw_ord_mitch)
mitch_long$contrast   <- factor(mitch_long$contrast, levels = CONTRASTS)

eff_lim <- max(abs(mitch_long$effect), na.rm = TRUE)

pH_H <- ggplot(mitch_long, aes(x = contrast, y = clean_name, fill = effect)) +
  geom_tile(color = "white", linewidth = 0.3) +
  scale_fill_gradient2(
    low = DIR_COLORS["Down"], mid = "white", high = DIR_COLORS["Up"],
    midpoint = 0, limits = c(-eff_lim, eff_lim),
    name = "Effect"
  ) +
  scale_x_discrete(labels = CTR_SHORT) +
  labs(title = "Option H: mitch Multi-Contrast",
       subtitle = sprintf("Top %d MANOVA-significant pathways\n(Kaspi & Ziemann 2020)",
                          n_top_mitch),
       x = NULL, y = NULL) +
  FIG_THEME +
  theme(axis.text.y        = element_text(size = 5),
        axis.text.x        = element_text(angle = 35, hjust = 1, size = 8),
        legend.position     = "right",
        legend.key.height   = unit(12, "mm"),
        legend.key.width    = unit(3, "mm"),
        panel.grid          = element_blank())

PHH_H <- max(140, n_top_mitch * 4 + 30)

ggsave(file.path(RPT, "option_H_mitch_heatmap.pdf"), pH_H,
       width = 150, height = PHH_H, units = "mm", device = pdf_device)
ggsave(file.path(RPT, "option_H_mitch_heatmap.png"), pH_H,
       width = 150, height = PHH_H, units = "mm", dpi = 300)
cat("Option H saved\n")


################################################################################
#   Option I: NES Bubble Plot (custom ggplot2)
#   x = contrast, y = pathway, size = -log10(padj), color = NES
#   Top N from fGSEA deduped results
################################################################################

cat("Building Option I: NES Bubble Plot\n")

# Take top 15 pathways per contrast by padj
bubble_top <- fgsea_combined |>
  filter(!is.na(padj), padj < 0.05) |>
  group_by(contrast) |>
  arrange(padj) |>
  slice_head(n = 15) |>
  ungroup()

# All unique pathways across contrasts
bubble_pws <- unique(bubble_top$pathway)

# Get full data for those pathways across all contrasts
bubble_df <- fgsea_combined |>
  filter(pathway %in% bubble_pws) |>
  mutate(
    significant     = !is.na(padj) & padj < 0.05,
    neg_log10_padj  = -log10(pmax(padj, 1e-20)),
    clean_name      = clean_pathway_name(pathway, max_chars = 45),
    contrast        = factor(contrast, levels = CONTRASTS)
  )

bubble_df$clean_name <- make.unique(bubble_df$clean_name, sep = " ")

# Order by cross-contrast density then mean |NES|
pw_ord_bubble <- bubble_df |>
  group_by(clean_name) |>
  summarise(
    n_sig        = sum(significant),
    mean_abs_nes = mean(abs(NES), na.rm = TRUE),
    .groups = "drop"
  ) |>
  arrange(n_sig, mean_abs_nes) |>
  pull(clean_name)

bubble_df$clean_name <- factor(bubble_df$clean_name, levels = pw_ord_bubble)
nes_lim_bubble <- max(abs(bubble_df$NES), na.rm = TRUE)

pH_I <- ggplot(bubble_df |> filter(significant),
               aes(x = contrast, y = clean_name)) +
  geom_point(aes(size = neg_log10_padj, color = NES)) +
  # Non-significant as small grey dots
  geom_point(data = bubble_df |> filter(!significant),
             aes(x = contrast, y = clean_name),
             size = 0.5, color = "grey80") +
  scale_color_gradient2(
    low = DIR_COLORS["Down"], mid = "white", high = DIR_COLORS["Up"],
    midpoint = 0, limits = c(-nes_lim_bubble, nes_lim_bubble),
    name = "NES"
  ) +
  scale_size_continuous(range = c(1, 6), name = expression(-log[10](p[adj]))) +
  scale_x_discrete(labels = CTR_SHORT) +
  labs(title = "Option I: NES Bubble Plot",
       subtitle = sprintf("Top 15 per contrast (%d unique pathways)\nSize = significance, Color = direction",
                          length(bubble_pws)),
       x = NULL, y = NULL) +
  FIG_THEME +
  theme(axis.text.y      = element_text(size = 5.5),
        axis.text.x      = element_text(angle = 35, hjust = 1, size = 8),
        legend.position   = "right",
        panel.grid.major  = element_line(color = "grey95"),
        panel.grid.minor  = element_blank())

n_bubble_pw <- length(unique(bubble_df$clean_name))
PHI_H <- max(140, n_bubble_pw * 4 + 30)

ggsave(file.path(RPT, "option_I_nes_bubble.pdf"), pH_I,
       width = 170, height = PHI_H, units = "mm", device = pdf_device)
ggsave(file.path(RPT, "option_I_nes_bubble.png"), pH_I,
       width = 170, height = PHI_H, units = "mm", dpi = 300)
cat("Option I saved\n")


################################################################################
#   Option J: Ridgeline Plots (ggridges)
#   Shows distribution of t-statistics for genes in each pathway
#   One plot per contrast, top 12 pathways each
################################################################################

cat("Building Option J: Ridgeline Plots\n")

suppressPackageStartupMessages(library(ggridges))

ridge_plots <- list()

for (ctr in CONTRASTS) {
  stats_col <- paste0("t_", ctr)
  stats_vec <- dep_df[[stats_col]]
  names(stats_vec) <- dep_df$gene

  # Top 12 pathways by padj in this contrast
  top12 <- fgsea_combined |>
    filter(contrast == ctr, !is.na(padj), padj < 0.05) |>
    arrange(padj) |>
    head(12)

  if (nrow(top12) == 0) next

  # Build long data: gene + t-stat for genes in each pathway
  ridge_data <- lapply(seq_len(nrow(top12)), function(i) {
    pw_genes <- pw_collection[[top12$pathway[i]]]
    pw_genes <- pw_genes[pw_genes %in% names(stats_vec)]
    data.frame(
      pathway    = clean_pathway_name(top12$pathway[i], max_chars = 40),
      t_stat     = stats_vec[pw_genes],
      NES        = top12$NES[i],
      stringsAsFactors = FALSE
    )
  })
  ridge_df <- bind_rows(ridge_data)

  # Order by NES
  pw_order <- top12 |>
    mutate(clean = clean_pathway_name(pathway, max_chars = 40)) |>
    arrange(NES) |>
    pull(clean)
  ridge_df$pathway <- factor(ridge_df$pathway, levels = pw_order)

  p_ridge <- ggplot(ridge_df, aes(x = t_stat, y = pathway, fill = after_stat(x))) +
    geom_density_ridges_gradient(scale = 1.2, rel_min_height = 0.01,
                                  gradient_lwd = 0.3) +
    scale_fill_gradient2(
      low = DIR_COLORS["Down"], mid = "white", high = DIR_COLORS["Up"],
      midpoint = 0, name = "t-statistic"
    ) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "grey50",
               linewidth = 0.4) +
    labs(title = sprintf("Option J: Ridgeline — %s", CTR_SHORT[ctr]),
         subtitle = "Top 12 pathways by padj\nGene t-statistic distributions",
         x = "t-statistic", y = NULL) +
    FIG_THEME +
    theme(axis.text.y      = element_text(size = 6),
          legend.position   = "none")

  ridge_plots[[ctr]] <- p_ridge
}

# Save each individually and also combined
for (ctr in names(ridge_plots)) {
  ggsave(file.path(RPT, sprintf("option_J_ridgeline_%s.pdf", ctr)),
         ridge_plots[[ctr]],
         width = 140, height = 120, units = "mm", device = pdf_device)
  ggsave(file.path(RPT, sprintf("option_J_ridgeline_%s.png", ctr)),
         ridge_plots[[ctr]],
         width = 140, height = 120, units = "mm", dpi = 300)
}

# Combined 2x2 grid
if (length(ridge_plots) == 4 && requireNamespace("patchwork", quietly = TRUE)) {
  library(patchwork)
  pH_J <- (ridge_plots[[1]] + ridge_plots[[2]]) /
           (ridge_plots[[3]] + ridge_plots[[4]]) +
    plot_annotation(
      title = "Option J: Ridgeline Plots (all contrasts)",
      theme = theme(plot.title = element_text(face = "bold", size = 12))
    )
  ggsave(file.path(RPT, "option_J_ridgeline_combined.pdf"), pH_J,
         width = 280, height = 240, units = "mm", device = pdf_device)
  ggsave(file.path(RPT, "option_J_ridgeline_combined.png"), pH_J,
         width = 280, height = 240, units = "mm", dpi = 300)
}
cat("Option J saved\n")


################################################################################
#   Option K: mitch Scatterplot (2D projections)
#   Pairwise contrasts showing pathway effect sizes
################################################################################

cat("Building Option K: mitch Scatter (pairwise contrasts)\n")

# Use mitch results — s. columns are effect sizes per contrast
if (exists("mitch_sig") && nrow(mitch_sig) > 0) {
  mitch_scatter_df <- mitch_sig |>
    select(set, starts_with("s."), p.adjustMANOVA) |>
    mutate(
      clean_name  = clean_pathway_name(set, max_chars = 35),
      database    = classify_database(set),
      significant = p.adjustMANOVA < 0.05,
      neg_log10_p = -log10(pmax(p.adjustMANOVA, 1e-20))
    )

  # Key comparison: Aging vs Training_Young
  pH_K1 <- ggplot(mitch_scatter_df,
                   aes(x = s.Aging, y = s.Training_Young)) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "grey60") +
    geom_vline(xintercept = 0, linetype = "dashed", color = "grey60") +
    geom_point(aes(color = database, size = neg_log10_p), alpha = 0.7) +
    scale_color_manual(values = DB_COLORS, name = "Database") +
    scale_size_continuous(range = c(0.5, 4),
                          name = expression(-log[10](p[adj]^MANOVA))) +
    labs(title = "mitch: Aging vs Training (Young)",
         subtitle = "Each point = 1 pathway",
         x = "Effect (Aging)", y = "Effect (Training Young)") +
    FIG_THEME +
    theme(legend.position = "right")

  # Aging vs Training_Old
  pH_K2 <- ggplot(mitch_scatter_df,
                   aes(x = s.Aging, y = s.Training_Old)) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "grey60") +
    geom_vline(xintercept = 0, linetype = "dashed", color = "grey60") +
    geom_point(aes(color = database, size = neg_log10_p), alpha = 0.7) +
    scale_color_manual(values = DB_COLORS, name = "Database") +
    scale_size_continuous(range = c(0.5, 4),
                          name = expression(-log[10](p[adj]^MANOVA))) +
    labs(title = "mitch: Aging vs Training (Old)",
         x = "Effect (Aging)", y = "Effect (Training Old)") +
    FIG_THEME +
    theme(legend.position = "right")

  # Training_Young vs Training_Old
  pH_K3 <- ggplot(mitch_scatter_df,
                   aes(x = s.Training_Young, y = s.Training_Old)) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "grey60") +
    geom_vline(xintercept = 0, linetype = "dashed", color = "grey60") +
    geom_abline(slope = 1, intercept = 0, linetype = "dotted",
                color = "grey40") +
    geom_point(aes(color = database, size = neg_log10_p), alpha = 0.7) +
    scale_color_manual(values = DB_COLORS, name = "Database") +
    scale_size_continuous(range = c(0.5, 4),
                          name = expression(-log[10](p[adj]^MANOVA))) +
    labs(title = "mitch: Training Young vs Old",
         x = "Effect (Training Young)", y = "Effect (Training Old)") +
    FIG_THEME +
    theme(legend.position = "right")

  # Combined 1x3
  library(patchwork)
  pH_K <- pH_K1 + pH_K2 + pH_K3 +
    plot_layout(ncol = 3, guides = "collect") +
    plot_annotation(
      title    = "Option K: mitch Pairwise Scatter",
      subtitle = sprintf("%d MANOVA-significant pathways (Kaspi & Ziemann 2020)",
                         nrow(mitch_sig)),
      theme = theme(
        plot.title    = element_text(face = "bold", size = 12),
        plot.subtitle = element_text(face = "bold.italic", size = 9,
                                      color = "grey30")
      )
    ) &
    theme(legend.position = "bottom")

  ggsave(file.path(RPT, "option_K_mitch_scatter.pdf"), pH_K,
         width = 350, height = 130, units = "mm", device = pdf_device)
  ggsave(file.path(RPT, "option_K_mitch_scatter.png"), pH_K,
         width = 350, height = 130, units = "mm", dpi = 300)
  cat("Option K saved\n")
} else {
  cat("Option K skipped: no mitch results\n")
}


################################################################################
#   Option L: enrichplot cnetplot + emapplot (network-based)
#   Network where nodes = enriched terms, edges = shared genes
################################################################################

cat("Building Option L: Network Visualizations\n")

# Use compareCluster result from Option G if available
if (exists("cc_res") && !is.null(cc_res) && nrow(as.data.frame(cc_res)) > 0) {

  # emapplot — enrichment map (pathway similarity network)
  # Build term-presence ggplot heatmap (binary: significant or not per contrast)
  tryCatch({
    cc_df <- as.data.frame(cc_res)
    top_terms <- cc_df |>
      group_by(Description) |>
      summarise(min_p = min(p.adjust), .groups = "drop") |>
      arrange(min_p) |>
      head(30) |>
      pull(Description)

    presence_long <- cc_df |>
      filter(Description %in% top_terms) |>
      select(Description, Cluster, p.adjust) |>
      mutate(neg_log10_p = -log10(p.adjust))

    # Order terms by number of contrasts enriched, then by min padj
    term_ord <- presence_long |>
      group_by(Description) |>
      summarise(n_ctr = n(), min_p = min(p.adjust), .groups = "drop") |>
      arrange(n_ctr, -min_p) |>
      pull(Description)
    presence_long$Description <- factor(presence_long$Description, levels = term_ord)

    pH_L <- ggplot(presence_long,
                    aes(x = Cluster, y = Description, fill = neg_log10_p)) +
      geom_tile(color = "white", linewidth = 0.5) +
      scale_fill_gradient(low = "#FFCDD2", high = DIR_COLORS["Up"],
                           name = expression(-log[10](p[adj]))) +
      labs(title = "Option L: Term Presence Heatmap",
           subtitle = "Top 30 ORA terms across contrasts\n(color = significance strength)",
           x = NULL, y = NULL) +
      FIG_THEME +
      theme(axis.text.y       = element_text(size = 6),
            axis.text.x       = element_text(angle = 35, hjust = 1, size = 8),
            panel.grid        = element_blank(),
            legend.position    = "right",
            legend.key.height  = unit(12, "mm"),
            legend.key.width   = unit(3, "mm"))

    ggsave(file.path(RPT, "option_L_term_presence.pdf"), pH_L,
           width = 170, height = 180, units = "mm", device = pdf_device)
    ggsave(file.path(RPT, "option_L_term_presence.png"), pH_L,
           width = 170, height = 180, units = "mm", dpi = 300)
    cat("Option L saved\n")
  }, error = function(e) {
    cat(sprintf("Option L error: %s\n", conditionMessage(e)))
  })
} else {
  cat("Option L skipped: no compareCluster results\n")
}


################################################################################
#   Option M: Treemap of pathway categories per contrast
#   Area = number of significant pathways, color = mean NES
################################################################################

cat("Building Option M: Treemap\n")

if (requireNamespace("treemapify", quietly = TRUE)) {
  library(treemapify)

  tree_df <- consol_fgsea |>
    group_by(contrast, expanded_cat) |>
    summarise(
      n_paths  = n(),
      mean_NES = mean(NES),
      .groups  = "drop"
    ) |>
    filter(n_paths > 0)

  tree_df$contrast <- factor(tree_df$contrast, levels = CONTRASTS)

  pH_M <- ggplot(tree_df,
                  aes(area = n_paths, fill = mean_NES,
                      label = paste0(expanded_cat, "\n(", n_paths, ")"))) +
    geom_treemap(color = "white", linewidth = 0.5) +
    geom_treemap_text(color = "black", place = "centre", size = 6,
                       min.size = 3) +
    facet_wrap(~ contrast, labeller = as_labeller(CTR_FACET)) +
    scale_fill_gradient2(
      low = DIR_COLORS["Down"], mid = "white", high = DIR_COLORS["Up"],
      midpoint = 0, name = "Mean NES"
    ) +
    labs(title = "Option M: Pathway Treemap",
         subtitle = "Area = # significant pathways\nColor = mean NES direction") +
    theme(strip.text       = element_text(face = "bold", size = 9),
          plot.title        = element_text(face = "bold", size = 11),
          plot.subtitle     = element_text(size = 8, color = "grey30"),
          legend.position   = "bottom",
          legend.key.width  = unit(20, "mm"),
          legend.key.height = unit(3, "mm"))

  ggsave(file.path(RPT, "option_M_treemap.pdf"), pH_M,
         width = 220, height = 180, units = "mm", device = pdf_device)
  ggsave(file.path(RPT, "option_M_treemap.png"), pH_M,
         width = 220, height = 180, units = "mm", dpi = 300)
  cat("Option M saved\n")
} else {
  cat("Option M skipped: treemapify not installed\n")
}


cat("\nAll Panel H options complete.\n")
