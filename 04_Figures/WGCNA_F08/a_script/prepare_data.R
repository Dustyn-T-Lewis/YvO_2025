# F6 prepare_data.R: Eigengene splitting + prediction correlations
# Run BEFORE panel scripts. Saves intermediates to c_data/.
#
# STAT AUDIT (2026-02-27):
#   GS (gs_vl, gs_lbm): Pearson r of baseline expression vs delta phenotype.
#   kME computed on full datExpr (all timepoints) — standard WGCNA practice.

setwd(rprojroot::find_rstudio_root_file())
source("04_Figures/shared/style.R")

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(tibble)
  library(WGCNA)
})

allowWGCNAThreads()
set.seed(42)

DAT <- "04_Figures/WGCNA_F08/c_data"
RPT <- "04_Figures/WGCNA_F08/b_reports"
dir.create(DAT, recursive = TRUE, showWarnings = FALSE)
dir.create(RPT, recursive = TRUE, showWarnings = FALSE)

wgcna_dir <- if (dir.exists("04_Figures/WGCNA_F07/c_data/wgcna")) {
  "04_Figures/WGCNA_F07/c_data/wgcna"
} else if (dir.exists("04_Figures_v1/F5/c_data/wgcna")) {
  "04_Figures_v1/F5/c_data/wgcna"
} else {
  stop("WGCNA data not found in v2 or v1. Run Figure 5 first.")
}
message("Using WGCNA data from: ", wgcna_dir)

module_df <- read_csv(file.path(wgcna_dir, "wgcna_module_assignments.csv"),
                       show_col_types = FALSE)

imp_df <- read_csv("02_Imputation/c_data/01_imputed.csv", show_col_types = FALSE)
ann_cols <- c("uniprot_id", "protein", "gene", "description")
samp_names <- setdiff(names(imp_df), ann_cols)
imp_mat <- as.matrix(imp_df[, samp_names])
rownames(imp_mat) <- imp_df$uniprot_id
datExpr <- t(imp_mat)

dal <- readRDS("02_Imputation/c_data/01_DAList_imputed.rds")
dal_meta <- as.data.frame(dal$metadata)

meta <- tibble(
  sample_id   = dal_meta$Col_ID,
  Subject_ID  = dal_meta$Subject_ID,
  age         = dal_meta$Group,
  time        = dal_meta$Timepoint,
  group       = dal_meta$Group_Time,
  subject_key = sub("_(Pre|Post)$", "", dal_meta$Col_ID)
)

for (pc in c("VL_thick_cm", "DXA_LBM_kg")) {
  if (pc %in% names(dal_meta)) {
    meta[[pc]] <- as.numeric(as.character(dal_meta[[pc]]))
  }
}

module_colors <- module_df$module_color
MEs <- moduleEigengenes(datExpr, colors = module_colors)$eigengenes
pre_meta  <- meta %>% filter(time == "Pre")
post_meta <- meta %>% filter(time == "Post")

me_pre_raw  <- MEs[pre_meta$sample_id, , drop = FALSE]
me_post_raw <- MEs[post_meta$sample_id, , drop = FALSE]
rownames(me_pre_raw)  <- pre_meta$subject_key
rownames(me_post_raw) <- post_meta$subject_key

common_subj <- intersect(rownames(me_pre_raw), rownames(me_post_raw))
me_pre  <- me_pre_raw[common_subj, , drop = FALSE]
me_post <- me_post_raw[common_subj, , drop = FALSE]
delta_me <- me_post - me_pre

pheno_pre <- meta %>%
  filter(time == "Pre") %>%
  dplyr::select(subject_key, VL_thick_cm, DXA_LBM_kg) %>%
  rename(VL_Pre = VL_thick_cm, LBM_Pre = DXA_LBM_kg)

pheno_post <- meta %>%
  filter(time == "Post") %>%
  dplyr::select(subject_key, VL_thick_cm, DXA_LBM_kg) %>%
  rename(VL_Post = VL_thick_cm, LBM_Post = DXA_LBM_kg)

pheno_wide <- inner_join(pheno_pre, pheno_post, by = "subject_key") %>%
  mutate(delta_VL  = VL_Post  - VL_Pre,
         delta_LBM = LBM_Post - LBM_Pre) %>%
  filter(subject_key %in% common_subj)

subj_age <- meta %>%
  filter(time == "Pre") %>%
  dplyr::select(subject_key, age) %>%
  distinct()

pre_expr <- datExpr[pre_meta$sample_id, , drop = FALSE]
rownames(pre_expr) <- pre_meta$subject_key
pre_expr <- pre_expr[common_subj, , drop = FALSE]

delta_vl_vec  <- pheno_wide$delta_VL[match(common_subj, pheno_wide$subject_key)]
delta_lbm_vec <- pheno_wide$delta_LBM[match(common_subj, pheno_wide$subject_key)]

gs_vl <- cor(pre_expr, delta_vl_vec, use = "pairwise.complete.obs")
colnames(gs_vl) <- "GS_deltaVL"

gs_lbm <- cor(pre_expr, delta_lbm_vec, use = "pairwise.complete.obs")
colnames(gs_lbm) <- "GS_deltaLBM"

kME_all <- signedKME(datExpr, MEs)
pred_cor <- tibble(module = colnames(me_pre)) %>%
  rowwise() %>%
  mutate(
    r_vl  = cor(me_pre[common_subj, module],
                pheno_wide$delta_VL[match(common_subj, pheno_wide$subject_key)],
                use = "complete.obs"),
    r_lbm = cor(me_pre[common_subj, module],
                pheno_wide$delta_LBM[match(common_subj, pheno_wide$subject_key)],
                use = "complete.obs"),
    max_r = max(abs(r_vl), abs(r_lbm), na.rm = TRUE)
  ) %>% ungroup() %>%
  mutate(across(c(r_vl, r_lbm, max_r), as.numeric))

top3 <- pred_cor %>% arrange(desc(max_r)) %>% head(3) %>% pull(module)

# Biological module labels (14 entries: 7 curated + 7 identity fallbacks)
mod_bio_labels <- c(
  blue        = "Catabolic Metabolism",
  brown       = "Cytoskeletal Remodeling",
  turquoise   = "Protein Modification",
  green       = "Muscle Structure",
  black       = "Translation Machinery",
  pink        = "Immune/Complement",
  yellow      = "Ubiquitin-Proteasome",
  red         = "Red",
  magenta     = "Magenta",
  purple      = "Purple",
  cyan        = "Cyan",
  tan         = "Tan",
  salmon      = "Salmon",
  greenyellow = "Greenyellow"
)

outcome_labels <- c(delta_VL  = "Delta VL (cm)",
                    delta_LBM = "Delta LBM (kg)")

saveRDS(MEs,          file.path(DAT, "MEs.rds"))
saveRDS(me_pre,       file.path(DAT, "me_pre.rds"))
saveRDS(me_post,      file.path(DAT, "me_post.rds"))
saveRDS(delta_me,     file.path(DAT, "delta_me.rds"))
saveRDS(pre_expr,     file.path(DAT, "pre_expr.rds"))
saveRDS(datExpr,      file.path(DAT, "datExpr.rds"))
saveRDS(kME_all,      file.path(DAT, "kME_all.rds"))
saveRDS(gs_vl,        file.path(DAT, "gs_vl.rds"))
saveRDS(gs_lbm,       file.path(DAT, "gs_lbm.rds"))
saveRDS(module_colors, file.path(DAT, "module_colors.rds"))

write_csv(meta,       file.path(DAT, "meta.csv"))
write_csv(pheno_wide, file.path(DAT, "pheno_wide.csv"))
write_csv(subj_age,   file.path(DAT, "subj_age.csv"))
write_csv(pred_cor,   file.path(DAT, "pred_cor.csv"))
write_csv(module_df,  file.path(DAT, "module_df.csv"))

saveRDS(list(
  top3           = top3,
  common_subj    = common_subj,
  mod_bio_labels = mod_bio_labels,
  outcome_labels = outcome_labels,
  delta_vl_vec   = delta_vl_vec,
  delta_lbm_vec  = delta_lbm_vec
), file.path(DAT, "shared_objects.rds"))

message(sprintf("F6 prepare_data.R complete: %d subjects, %d modules, %d proteins",
                length(common_subj), ncol(MEs), ncol(datExpr)))
