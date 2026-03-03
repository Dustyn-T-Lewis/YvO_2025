################################################################################
#   YvO Figure 0 — Shared Setup
#   Sources: libraries, palettes, metadata, phenotype summaries, and statistics
#   Sourced by each panel_X.R script
################################################################################

# 0. LIBRARIES ================================================================

suppressPackageStartupMessages({
  library(readxl)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(patchwork)
  library(rstatix)
  library(ggsignif)
  library(scales)
})

setwd(rprojroot::find_rstudio_root_file())

# 1. PALETTES & THEME ========================================================

source("04_Figures/shared/palettes.R")

RPT_DIR <- "04_Figures/F0/b_reports"
DAT_DIR <- "04_Figures/F0/c_data"
dir.create(RPT_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(DAT_DIR, recursive = TRUE, showWarnings = FALSE)

# 2. LOAD METADATA ============================================================

meta <- read_excel("00_input/YvO_meta.xlsx")

# Coerce character columns to numeric
char_to_num <- c("BMI", "Type_I_fCSA", "Type_II_fCSA",
                 "deadlift_1rm_kg", "Total_Training_Volume_kg")
for (col in char_to_num) {
  if (col %in% names(meta) && is.character(meta[[col]])) {
    meta[[col]] <- suppressWarnings(as.numeric(meta[[col]]))
  }
}

# Create subject key and factor levels
meta <- meta %>%
  mutate(
    subject_key = sub("_(Pre|Post)$", "", Col_ID),
    Group       = factor(Group, levels = c("Young", "Old")),
    Timepoint   = factor(Timepoint, levels = c("Pre", "Post")),
    Group_Time  = factor(Group_Time,
                         levels = c("Young_Pre", "Young_Post",
                                    "Old_Pre", "Old_Post"))
  )

cat(sprintf("Loaded: %d samples, %d subjects\n",
            nrow(meta), length(unique(meta$subject_key))))

# 3. TRAINING VOLUME (Panel A) ===============================================

# Training volume is only in Post rows — extract one per subject
tv_df <- meta %>%
  filter(Timepoint == "Post") %>%
  select(subject_key, Group, Total_Training_Volume_kg) %>%
  rename(tv = Total_Training_Volume_kg) %>%
  filter(!is.na(tv))

stats_A <- t.test(tv ~ Group, data = tv_df)

tv_summary <- tv_df %>%
  group_by(Group) %>%
  summarise(
    n    = n(),
    mean = mean(tv),
    sd   = sd(tv),
    sem  = sd(tv) / sqrt(n()),
    .groups = "drop"
  )

# 4. DXA LEAN MASS & VL THICKNESS (Panels B, C) ==============================

# Wide format: one row per subject with Pre/Post/delta values
pheno_wide <- meta %>%
  select(subject_key, Group, Timepoint, DXA_LBM_kg, VL_thick_cm) %>%
  pivot_wider(
    names_from  = Timepoint,
    values_from = c(DXA_LBM_kg, VL_thick_cm),
    names_sep   = "_"
  ) %>%
  rename(
    DXA_Pre  = DXA_LBM_kg_Pre,
    DXA_Post = DXA_LBM_kg_Post,
    VL_Pre   = VL_thick_cm_Pre,
    VL_Post  = VL_thick_cm_Post
  ) %>%
  mutate(
    delta_DXA = DXA_Post - DXA_Pre,
    delta_VL  = VL_Post  - VL_Pre
  )

# 5. MIXED ANOVA — DXA (Panel B) =============================================

dxa_long <- meta %>%
  select(subject_key, Group, Timepoint, DXA_LBM_kg) %>%
  filter(!is.na(DXA_LBM_kg))

stats_B_anova <- rstatix::anova_test(
  data    = dxa_long,
  dv      = DXA_LBM_kg,
  wid     = subject_key,
  between = Group,
  within  = Timepoint
)

# 6. MIXED ANOVA — VL (Panel C) ==============================================

vl_long <- meta %>%
  select(subject_key, Group, Timepoint, VL_thick_cm) %>%
  filter(!is.na(VL_thick_cm))

stats_C_anova <- rstatix::anova_test(
  data    = vl_long,
  dv      = VL_thick_cm,
  wid     = subject_key,
  between = Group,
  within  = Timepoint
)

# 7. PAIRED T-TESTS WITHIN GROUP & BETWEEN-GROUP DELTA T-TESTS ===============

# DXA — paired t-tests within each group
dxa_young <- pheno_wide %>% filter(Group == "Young")
dxa_old   <- pheno_wide %>% filter(Group == "Old")

stats_B_paired_young <- t.test(dxa_young$DXA_Post, dxa_young$DXA_Pre, paired = TRUE)
stats_B_paired_old   <- t.test(dxa_old$DXA_Post,   dxa_old$DXA_Pre,   paired = TRUE)
stats_B_delta        <- t.test(delta_DXA ~ Group, data = pheno_wide)

# VL — paired t-tests within each group
vl_young <- pheno_wide %>% filter(Group == "Young")
vl_old   <- pheno_wide %>% filter(Group == "Old")

stats_C_paired_young <- t.test(vl_young$VL_Post, vl_young$VL_Pre, paired = TRUE)
stats_C_paired_old   <- t.test(vl_old$VL_Post,   vl_old$VL_Pre,   paired = TRUE)
stats_C_delta        <- t.test(delta_VL ~ Group, data = pheno_wide)

# 8. GROUP SUMMARY (means + SEM for bar plots) ===============================

group_summary <- meta %>%
  select(Group, Timepoint, Group_Time, DXA_LBM_kg, VL_thick_cm) %>%
  group_by(Group, Timepoint, Group_Time) %>%
  summarise(
    dxa_mean = mean(DXA_LBM_kg, na.rm = TRUE),
    dxa_sem  = sd(DXA_LBM_kg, na.rm = TRUE) / sqrt(sum(!is.na(DXA_LBM_kg))),
    vl_mean  = mean(VL_thick_cm, na.rm = TRUE),
    vl_sem   = sd(VL_thick_cm, na.rm = TRUE) / sqrt(sum(!is.na(VL_thick_cm))),
    .groups  = "drop"
  )

cat("F0 setup complete.\n")
