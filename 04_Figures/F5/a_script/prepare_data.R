# F5 prepare_data.R: Load WGCNA outputs, compute MEs/kME, save intermediates.
# Run BEFORE panel scripts. STAT AUDIT: corPvalueStudent() uses n=62 for ALL
# traits — p-values slightly liberal for phenotype traits with missing data.

setwd(rprojroot::find_rstudio_root_file())
source("04_Figures/shared/style.R")

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(stringr)
  library(tibble)
  library(WGCNA)
})

allowWGCNAThreads()
set.seed(42)

DAT <- "04_Figures/F5/c_data"
RPT <- "04_Figures/F5/b_reports"
dir.create(DAT, recursive = TRUE, showWarnings = FALSE)
dir.create(RPT, recursive = TRUE, showWarnings = FALSE)

wgcna_dir <- file.path(DAT, "wgcna")
stopifnot(dir.exists(wgcna_dir))

net       <- readRDS(file.path(wgcna_dir, "wgcna_network.rds"))
module_df <- read_csv(file.path(wgcna_dir, "wgcna_module_assignments.csv"), show_col_types = FALSE)

module_colors <- module_df$module_color
n_modules <- length(unique(module_colors[module_colors != "grey"]))

imp_df <- read_csv("02_Imputation/c_data/01_imputed.csv", show_col_types = FALSE)
ann_cols <- c("uniprot_id", "protein", "gene", "description")
samp_names <- setdiff(names(imp_df), ann_cols)
imp_mat <- as.matrix(imp_df[, samp_names])
rownames(imp_mat) <- imp_df$uniprot_id

datExpr <- t(imp_mat)  # WGCNA convention: samples as rows

dal <- readRDS("02_Imputation/c_data/01_DAList_imputed.rds")
dal_meta <- as.data.frame(dal$metadata)
meta <- tibble(
  sample_id = dal_meta$Col_ID,
  subject   = sub("_(Pre|Post)$", "", dal_meta$Col_ID),
  age       = dal_meta$Group,
  time      = dal_meta$Timepoint,
  group     = dal_meta$Group_Time
)
meta$group <- factor(meta$group,
                     levels = c("Young_Pre", "Young_Post", "Old_Pre", "Old_Post"))

pheno_cols <- c("VL_thick_cm", "DXA_LBM_kg", "BMI",
                "Type_I_fCSA", "Type_II_fCSA", "deadlift_1rm_kg")
for (pc in pheno_cols) {
  if (pc %in% names(dal_meta)) {
    vals <- dal_meta[[pc]]
    if (!is.numeric(vals)) vals <- as.numeric(as.character(vals))
    meta[[pc]] <- vals[match(meta$sample_id, dal_meta$Col_ID)]
  }
}

MEs <- moduleEigengenes(datExpr, colors = module_colors)$eigengenes

kME_all <- signedKME(datExpr, MEs)

mod_bio_labels <- tibble(
  module_color = c("blue", "brown", "turquoise", "green", "black",
                   "pink", "yellow", "red", "magenta"),
  bio_label    = c("Catabolic Metabolism", "Cytoskeletal Remodeling",
                   "Protein Modification", "Muscle Structure",
                   "Translation Machinery", "Immune/Complement",
                   "Ubiquitin-Proteasome", "Red", "Magenta")
)

saveRDS(MEs,        file.path(DAT, "MEs.rds"))
saveRDS(kME_all,    file.path(DAT, "kME_all.rds"))
saveRDS(datExpr,    file.path(DAT, "datExpr.rds"))
saveRDS(imp_mat,    file.path(DAT, "imp_mat.rds"))
saveRDS(module_colors, file.path(DAT, "module_colors.rds"))
write_csv(meta,     file.path(DAT, "meta.csv"))
write_csv(mod_bio_labels, file.path(DAT, "mod_bio_labels.csv"))

imp_annot <- imp_df[, ann_cols]
write_csv(imp_annot, file.path(DAT, "imp_annotations.csv"))

message(sprintf("F5 prepare_data.R complete: %d modules, %d proteins, %d samples",
                n_modules, ncol(datExpr), nrow(datExpr)))
