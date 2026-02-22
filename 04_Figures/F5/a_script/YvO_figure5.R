################################################################################
#   Figure 5: WGCNA Co-Expression Architecture
#
#   Study: YvO — DIA-MS Proteomics of Skeletal Muscle
#   Design: 2x2 mixed (Age x Time) with repeated measures
#
#   Central question: What co-expression modules underlie the proteomic
#   landscape, how do they relate to experimental factors and phenotypic
#   outcomes, and which modules show age-dependent training responses?
#
#   Panels (to be added in later tasks):
#     A — Module dendrogram with color bar
#     B — Module-Trait correlation heatmap
#     C — Eigenprotein Pre/Post trajectories (key modules)
#     D — Pathway enrichment per module (GO:BP + Hallmark)
#     E — kME scatter (hub proteins vs module membership)
#     F — FCM x WGCNA cross-tabulation / overlap
#
#   Inputs:
#     - WGCNA network object (wgcna_network.rds)
#     - Module assignments (wgcna_module_assignments.csv)
#     - Hub proteins (wgcna_hub_proteins.csv)
#     - Module-trait correlations (wgcna_module_trait_correlations.csv)
#     - GO enrichment (wgcna_module_GO_enrichment.csv)
#     - Imputed expression matrix (01_imputed.csv)
#     - DAList with metadata (01_DAList_imputed.rds)
#     - DEP results (combined_results.csv)
################################################################################

# === 1. PACKAGES =============================================================

# Save stats::cor BEFORE loading WGCNA (which overwrites it)
.orig_cor <- stats::cor

suppressPackageStartupMessages({
  library(WGCNA)
  library(tidyverse)
  library(patchwork)
  library(clusterProfiler)
  library(msigdbr)
  library(org.Hs.eg.db)
  library(ggrepel)
  library(scales)
  library(grid)
})

# WGCNA::cor conflict management
# Override with WGCNA::cor for all WGCNA operations; restore on exit
cor <- WGCNA::cor
on.exit(assign("cor", .orig_cor, envir = .GlobalEnv), add = TRUE)

# Suppress WGCNA thread warnings
allowWGCNAThreads()

# === 2. SEED =================================================================

set.seed(42)

# === 3. PATH RESOLUTION ======================================================

.script_dir <- tryCatch(dirname(normalizePath(sys.frame(1)$ofile)),
                         error = function(e) {
                           args <- commandArgs(trailingOnly = FALSE)
                           f <- grep("^--file=", args, value = TRUE)
                           if (length(f)) dirname(normalizePath(sub("^--file=", "", f[1])))
                           else normalizePath(".")
                         })
BASE_DIR <- normalizePath(file.path(.script_dir, "..", "..", ".."))
FIG_DIR  <- normalizePath(file.path(.script_dir, ".."))
RPT_DIR  <- file.path(FIG_DIR, "b_reports")
DAT_DIR  <- file.path(FIG_DIR, "c_data")

# === 4. DIRECTORY CREATION ===================================================

dir.create(RPT_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(DAT_DIR, recursive = TRUE, showWarnings = FALSE)

# === 5. FILE PATH CONSTANTS ==================================================

IMPUTED_FILE <- file.path(BASE_DIR, "02_Imputation", "c_data", "01_imputed.csv")
DALIST_RDS   <- file.path(BASE_DIR, "02_Imputation", "c_data", "01_DAList_imputed.rds")
DEP_FILE     <- file.path(BASE_DIR, "03_DEP", "c_data", "combined_results.csv")
WGCNA_DIR    <- file.path(BASE_DIR, "05_WGCNA", "c_data")
MOD_ASSIGN   <- file.path(WGCNA_DIR, "wgcna_module_assignments.csv")
HUB_FILE     <- file.path(WGCNA_DIR, "wgcna_hub_proteins.csv")
TRAIT_COR    <- file.path(WGCNA_DIR, "wgcna_module_trait_correlations.csv")
GO_ENRICH    <- file.path(WGCNA_DIR, "wgcna_module_GO_enrichment.csv")
NET_RDS      <- file.path(WGCNA_DIR, "wgcna_network.rds")

# === 6. CANONICAL CONSTANTS ==================================================

AGE_COLORS <- c(Young = "#4393C3", Old = "#D6604D")
DB_COLORS  <- c(Hallmark = "#AA336A", "GO:BP" = "#00796B",
                "GO:CC" = "#26A69A", "GO:MF" = "#CD5C5C")
CLUSTER_COLORS <- c(C1 = "#E74C3C", C2 = "#3498DB", C3 = "#2ECC71",
                     C4 = "#F39C12", C5 = "#9B59B6", C6 = "#1ABC9C",
                     C7 = "#E67E22", C8 = "#34495E", C9 = "#D35400",
                     C10 = "#7F8C8D")
KEY_TEXT  <- 2.2
KEY_TITLE <- 2.3

GROUP_FILL <- c(
  Young_Pre  = scales::alpha("#4393C3", 0.5),
  Young_Post = "#4393C3",
  Old_Pre    = scales::alpha("#D6604D", 0.5),
  Old_Post   = "#D6604D"
)

THEME_PUB <- theme_bw(base_size = 8) +
  theme(plot.title       = element_text(face = "bold", size = 9),
        plot.subtitle    = element_text(size = 6.5, color = "grey30", face = "italic"),
        strip.background = element_blank(),
        strip.text       = element_text(face = "bold", size = 6.5),
        legend.key.size  = unit(3, "mm"))

# === 7. MODULE BIOLOGY LOOKUP ================================================

MODULE_BIOLOGY <- c(
  green       = "Mito. ETC / OxPhos",
  blue        = "Lipid/AA Catabolism",
  turquoise   = "RNA Processing/Splicing",
  pink        = "Ribosome (60S)",
  salmon      = "Ribosome (40S)",
  red         = "Sarcomeric/Contractile",
  magenta     = "ECM/Collagen",
  black       = "Complement/Serum",
  purple      = "ATP Synthase",
  cyan        = "Translation/Proteasome",
  yellow      = "Cytoskeletal/Transport",
  tan         = "Filamin Network",
  greenyellow = "Chaperonin/Folding",
  brown       = "Proteasome/Ubiquitin"
)

# === 8. HELPER FUNCTIONS =====================================================

clean_pathway_name <- function(name) {
  name |>
    str_remove("^HALLMARK_") |>
    str_remove("^GOBP_") |>
    str_remove("^GOCC_") |>
    str_remove("^GOMF_") |>
    str_replace_all("_", " ") |>
    str_to_title() |>
    str_replace("Mtorc1", "mTORC1") |>
    str_replace("Myc ", "MYC ") |>
    str_replace("E2f ", "E2F ") |>
    str_replace("Dna ", "DNA ") |>
    str_replace("Rna ", "RNA ") |>
    str_replace("Tnfa ", "TNFa ") |>
    str_replace("Uv ", "UV ") |>
    str_replace("G2m ", "G2M ") |>
    str_replace("Il6 ", "IL6 ") |>
    str_replace("Il2 ", "IL2 ") |>
    str_replace("Kras ", "KRAS ") |>
    str_replace("P53 ", "p53 ") |>
    str_replace("Tgf ", "TGF ") |>
    str_replace("Nf Kb", "NF-kB") |>
    str_replace("Atp ", "ATP ") |>
    str_replace("Nadh ", "NADH ") |>
    str_trunc(45, ellipsis = "...")
}

sig_stars <- function(padj) {
  dplyr::case_when(
    padj < 0.001 ~ "***",
    padj < 0.01  ~ "**",
    padj < 0.05  ~ "*",
    TRUE         ~ ""
  )
}

reorder_within <- function(x, by, within, fun = mean, sep = "___", ...) {
  new_x <- paste(x, within, sep = sep)
  stats::reorder(new_x, by, FUN = fun, ...)
}

scale_y_reordered <- function(..., sep = "___") {
  reg <- paste0(sep, ".+$")
  ggplot2::scale_y_discrete(labels = function(x) gsub(reg, "", x), ...)
}

# === 9. DATA LOADING =========================================================

cat("Loading data...\n")

# Module assignments (uniprot_id, module_color, module_num, gene)
mod_df <- read_csv(MOD_ASSIGN, show_col_types = FALSE)
cat(sprintf("  Module assignments: %d proteins, %d modules (incl. grey)\n",
            nrow(mod_df), length(unique(mod_df$module_color))))

# Hub proteins (uniprot_id, kME, module, gene)
hub_df <- read_csv(HUB_FILE, show_col_types = FALSE)
cat(sprintf("  Hub proteins: %d\n", nrow(hub_df)))

# Module-trait correlations (correlation values only, no p-values)
trait_cor_df <- read_csv(TRAIT_COR, show_col_types = FALSE)
cat(sprintf("  Trait correlations: %d modules x %d traits\n",
            nrow(trait_cor_df), ncol(trait_cor_df) - 1))

# GO enrichment results
go_enrich <- read_csv(GO_ENRICH, show_col_types = FALSE)
cat(sprintf("  GO enrichment: %d terms across %d modules\n",
            nrow(go_enrich), length(unique(go_enrich$module))))

# DEP results
dep_df <- read_csv(DEP_FILE, show_col_types = FALSE)
cat(sprintf("  DEP results: %d proteins\n", nrow(dep_df)))

# Imputed expression matrix
imp_df <- read_csv(IMPUTED_FILE, show_col_types = FALSE)
ann_cols <- c("uniprot_id", "protein", "gene", "description")
samp_names <- setdiff(names(imp_df), ann_cols)
mat <- as.matrix(imp_df[, samp_names])
rownames(mat) <- imp_df$uniprot_id
cat(sprintf("  Imputed matrix: %d proteins x %d samples\n", nrow(mat), ncol(mat)))

# DAList (for phenotype metadata)
dal <- readRDS(DALIST_RDS)

# === 10. BUILD datExpr (WGCNA CONVENTION: samples x proteins) ================

datExpr <- t(mat)
cat(sprintf("  datExpr: %d samples x %d proteins\n", nrow(datExpr), ncol(datExpr)))

# === 11. PARSE SAMPLE METADATA ===============================================

meta <- tibble(sample_id = samp_names) |>
  mutate(
    # Age group: O_ and OP_ prefixes = Old; Y_ and YP_ = Young
    age     = if_else(str_detect(sample_id, "^(O_|OP_)"), "Old", "Young"),
    # Timepoint: _Pre or _Post suffix
    time    = if_else(str_detect(sample_id, "_Post$"), "Post", "Pre"),
    # Subject ID: strip _Pre/_Post suffix
    subject = str_remove(sample_id, "_(Pre|Post)$"),
    # Group label
    group   = paste(age, time, sep = "_")
  )

cat(sprintf("  Samples: %d total (%d Young, %d Old)\n",
            nrow(meta), sum(meta$age == "Young"), sum(meta$age == "Old")))

# === 12. BUILD TRAIT MATRIX ==================================================

# Experimental design columns
traits <- meta |>
  mutate(
    age_num     = if_else(age == "Old", 1, 0),
    time_num    = if_else(time == "Post", 1, 0),
    interaction = age_num * time_num
  ) |>
  dplyr::select(sample_id, age_num, time_num, interaction)

# Phenotype columns from DAList metadata
if (!is.null(dal) && "metadata" %in% names(dal)) {
  dal_meta <- as.data.frame(dal$metadata)
  id_col <- if ("Col_ID" %in% names(dal_meta)) "Col_ID" else names(dal_meta)[1]

  pheno_cols <- c("VL_thick_cm", "DXA_LBM_kg", "BMI",
                   "Type_I_fCSA", "Type_II_fCSA", "deadlift_1rm_kg")

  # Convert character columns to numeric where needed
  for (pc in pheno_cols) {
    if (pc %in% names(dal_meta)) {
      vals <- dal_meta[[pc]]
      if (is.character(vals)) {
        vals[vals == "NA" | vals == "na" | vals == ""] <- NA
        vals <- as.numeric(vals)
      }
      dal_meta[[pc]] <- vals
    }
  }

  pheno_available <- intersect(pheno_cols, names(dal_meta))
  pheno_df <- dal_meta |>
    dplyr::select(all_of(c(id_col, pheno_available))) |>
    dplyr::rename(sample_id = !!id_col)

  traits <- traits |> left_join(pheno_df, by = "sample_id")
  cat(sprintf("  Phenotype columns joined: %s\n", paste(pheno_available, collapse = ", ")))
}

# Convert to numeric matrix aligned with datExpr row order
traits_mat <- traits |>
  column_to_rownames("sample_id") |>
  mutate(across(everything(), as.numeric))
traits_mat <- traits_mat[rownames(datExpr), , drop = FALSE]

cat(sprintf("  Trait matrix: %d samples x %d traits\n",
            nrow(traits_mat), ncol(traits_mat)))

# === 13. RECOMPUTE MODULE EIGENGENES =========================================

# Filter datExpr to only WGCNA-assigned proteins (those in module assignments)
wgcna_proteins <- mod_df$uniprot_id[mod_df$uniprot_id %in% colnames(datExpr)]
datExpr_wgcna  <- datExpr[, wgcna_proteins]
module_colors  <- setNames(mod_df$module_color, mod_df$uniprot_id)[wgcna_proteins]

cat(sprintf("  WGCNA subset: %d proteins\n", ncol(datExpr_wgcna)))

# moduleEigengenes returns ME{color} naming (MEturquoise, MEblue, etc.)
MEs <- moduleEigengenes(datExpr_wgcna, colors = module_colors)$eigengenes
MEs <- orderMEs(MEs)

cat(sprintf("  Module eigengenes: %d samples x %d modules (%s)\n",
            nrow(MEs), ncol(MEs), paste(colnames(MEs), collapse = ", ")))

# === 14. RECOMPUTE MODULE-TRAIT CORRELATIONS WITH P-VALUES ===================

# The stored trait_cor_df has correlation values only — recompute with p-values
module_trait_cor  <- WGCNA::cor(MEs, traits_mat, use = "p")
module_trait_pval <- corPvalueStudent(module_trait_cor, nrow(datExpr))

cat(sprintf("  Module-trait correlation matrix: %d x %d (with p-values)\n",
            nrow(module_trait_cor), ncol(module_trait_cor)))

# === 15. COMPUTE DELTA TRAITS PER SUBJECT ====================================

# Add phenotype data to meta for delta computation
meta_pheno <- meta |>
  left_join(traits |> dplyr::select(-age_num, -time_num, -interaction),
            by = "sample_id")

delta_traits <- meta_pheno |>
  group_by(subject) |>
  summarise(
    delta_VL  = VL_thick_cm[time == "Post"] - VL_thick_cm[time == "Pre"],
    delta_LBM = DXA_LBM_kg[time == "Post"] - DXA_LBM_kg[time == "Pre"],
    age       = dplyr::first(age),
    .groups   = "drop"
  )

cat(sprintf("  Delta traits: %d subjects (%d with non-NA delta_VL, %d with non-NA delta_LBM)\n",
            nrow(delta_traits),
            sum(!is.na(delta_traits$delta_VL)),
            sum(!is.na(delta_traits$delta_LBM))))

# === 16. COMPUTE kME (SIGNED MODULE MEMBERSHIP) ==============================

# signedKME returns columns named kME{color} (e.g., kMEturquoise, kMEblue)
kME_all <- signedKME(datExpr_wgcna, MEs)

cat(sprintf("  kME matrix: %d proteins x %d modules (%s ...)\n",
            nrow(kME_all), ncol(kME_all),
            paste(head(colnames(kME_all), 4), collapse = ", ")))

# === 17. LOAD WGCNA NETWORK OBJECT ===========================================

# Network object contains: colors (numeric 0-14), unmergedColors, MEs (ME0-ME14),
# dendrograms (list of 1), TOMFiles, blockGenes, blocks, MEsOK
net_obj <- readRDS(NET_RDS)

cat(sprintf("  Network object loaded: %d genes in block 1, %d dendrograms\n",
            length(net_obj$blockGenes[[1]]), length(net_obj$dendrograms)))

# === 18. MODULE SIZE COUNTS ==================================================

mod_sizes <- mod_df |>
  filter(module_color != "grey") |>
  count(module_color, name = "n_proteins")

cat(sprintf("  Module sizes (excluding grey): %d modules, range %d-%d proteins\n",
            nrow(mod_sizes), min(mod_sizes$n_proteins), max(mod_sizes$n_proteins)))

# === SETUP COMPLETE ===========================================================

cat("\n=== Setup complete ===\n")
cat(sprintf("  %d modules (excl. grey), %d proteins, %d samples\n",
            length(unique(mod_df$module_color[mod_df$module_color != "grey"])),
            nrow(mod_df), nrow(datExpr)))
cat(sprintf("  MEs recomputed: %s\n", paste(colnames(MEs), collapse = ", ")))
cat(sprintf("  kME columns: %s\n", paste(head(colnames(kME_all), 5), collapse = ", ")))
cat(sprintf("  Traits: %s\n", paste(colnames(traits_mat), collapse = ", ")))
cat("  Ready for panel construction.\n")
