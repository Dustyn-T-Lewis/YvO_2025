################################################################################
#   Figure 0 — Composite Assembly + Supplementary Excel Workbook
#
#   Sources all per-panel scripts, assembles with patchwork into a single
#   Figure_0.pdf/png, and exports F0_supplementary.xlsx.
#
#   Panels:
#     A — Training volume bar chart (Young vs Old)
#     B — DXA lean body mass (Pre/Post + delta)
#     C — VL thickness (Pre/Post + delta)
################################################################################

# === 0. Source all panel scripts (setup loaded automatically by panel_A) =====

source("04_Figures_v1/F0/a_script/panel_A.R")   # -> pA
source("04_Figures_v1/F0/a_script/panel_B.R")   # -> pB
source("04_Figures_v1/F0/a_script/panel_C.R")   # -> pC

message("All F0 panel scripts sourced - assembling composite...")

# === 1. Composite figure (patchwork layout) ==================================

fig0 <- (pA / pB / pC) + plot_layout(heights = c(0.33, 0.33, 0.33))

ggsave(file.path(RPT_DIR, "Figure_0.pdf"), fig0,
       width = 180, height = 240, units = "mm", device = pdf)
ggsave(file.path(RPT_DIR, "Figure_0.png"), fig0,
       width = 180, height = 240, units = "mm", dpi = 300)

cat("Saved Figure_0.pdf and Figure_0.png\n")

# === 2. Supplementary Excel workbook ========================================

library(openxlsx)

# Sheet 1: A_training_volume
sheet_A <- tv_df %>%
  left_join(tv_summary %>% select(Group, mean, sem), by = "Group") %>%
  rename(group_mean = mean, group_sem = sem)

# Sheet 2: B_dxa_lean_mass (Pre/Post group summary)
sheet_B <- group_summary %>%
  select(Group, Timepoint, Group_Time, dxa_mean, dxa_sem)

# Sheet 3: B_dxa_deltas (per-subject deltas)
sheet_B_delta <- pheno_wide %>%
  select(subject_key, Group, DXA_Pre, DXA_Post, delta_DXA)

# Sheet 4: C_vl_thickness (Pre/Post group summary)
sheet_C <- group_summary %>%
  select(Group, Timepoint, Group_Time, vl_mean, vl_sem)

# Sheet 5: C_vl_deltas (per-subject deltas)
sheet_C_delta <- pheno_wide %>%
  select(subject_key, Group, VL_Pre, VL_Post, delta_VL)

# Sheet 6: statistics (all test results)
# Collect ANOVA results — strip rstatix class to avoid vctrs conflicts
anova_B <- tibble(measure = "DXA_LBM_kg",
                  Effect  = stats_B_anova$Effect,
                  DFn     = stats_B_anova$DFn,
                  DFd     = stats_B_anova$DFd,
                  F_stat  = stats_B_anova$F,
                  p       = stats_B_anova$p,
                  ges     = stats_B_anova$ges)
anova_C <- tibble(measure = "VL_thick_cm",
                  Effect  = stats_C_anova$Effect,
                  DFn     = stats_C_anova$DFn,
                  DFd     = stats_C_anova$DFd,
                  F_stat  = stats_C_anova$F,
                  p       = stats_C_anova$p,
                  ges     = stats_C_anova$ges)

# Collect t-test results
ttest_rows <- tibble(
  test = c("Training volume (Young vs Old)",
           "DXA paired (Young Pre vs Post)",
           "DXA paired (Old Pre vs Post)",
           "DXA delta (Young vs Old)",
           "VL paired (Young Pre vs Post)",
           "VL paired (Old Pre vs Post)",
           "VL delta (Young vs Old)"),
  statistic = c(stats_A$statistic,
                stats_B_paired_young$statistic,
                stats_B_paired_old$statistic,
                stats_B_delta$statistic,
                stats_C_paired_young$statistic,
                stats_C_paired_old$statistic,
                stats_C_delta$statistic),
  df = c(stats_A$parameter,
         stats_B_paired_young$parameter,
         stats_B_paired_old$parameter,
         stats_B_delta$parameter,
         stats_C_paired_young$parameter,
         stats_C_paired_old$parameter,
         stats_C_delta$parameter),
  p_value = c(stats_A$p.value,
              stats_B_paired_young$p.value,
              stats_B_paired_old$p.value,
              stats_B_delta$p.value,
              stats_C_paired_young$p.value,
              stats_C_paired_old$p.value,
              stats_C_delta$p.value),
  method = c(stats_A$method,
             stats_B_paired_young$method,
             stats_B_paired_old$method,
             stats_B_delta$method,
             stats_C_paired_young$method,
             stats_C_paired_old$method,
             stats_C_delta$method)
)

sheet_stats <- list(anova = bind_rows(anova_B, anova_C),
                    ttests = ttest_rows)

# Sheet 7: _metadata
sheet_meta <- tibble(
  item = c(
    "figure_title",
    "description",
    "groups",
    "statistical_methods",
    "sample_sizes",
    "software",
    "date_generated"
  ),
  value = c(
    "Figure 0: Participant Phenotype Characteristics",
    "Training volume (Young vs Old), DXA lean body mass (Pre/Post x Group), and vastus lateralis thickness (Pre/Post x Group).",
    "Young (n=17), Old (n=15); Pre and Post resistance training",
    "Welch t-test (training volume, deltas), mixed ANOVA (Group x Timepoint), paired t-tests (within-group Pre/Post)",
    sprintf("%d subjects (Young = %d, Old = %d), 64 observations",
            length(unique(meta$subject_key)),
            sum(tv_df$Group == "Young"), sum(tv_df$Group == "Old")),
    "R, rstatix, ggsignif, patchwork",
    format(Sys.Date(), "%Y-%m-%d")
  )
)

# Build workbook
wb <- createWorkbook()

addWorksheet(wb, "A_training_volume");   writeData(wb, "A_training_volume",   sheet_A)
addWorksheet(wb, "B_dxa_lean_mass");     writeData(wb, "B_dxa_lean_mass",     sheet_B)
addWorksheet(wb, "B_dxa_deltas");        writeData(wb, "B_dxa_deltas",        sheet_B_delta)
addWorksheet(wb, "C_vl_thickness");      writeData(wb, "C_vl_thickness",      sheet_C)
addWorksheet(wb, "C_vl_deltas");         writeData(wb, "C_vl_deltas",         sheet_C_delta)
addWorksheet(wb, "statistics_anova");    writeData(wb, "statistics_anova",    sheet_stats$anova)
addWorksheet(wb, "statistics_ttests");   writeData(wb, "statistics_ttests",   sheet_stats$ttests)
addWorksheet(wb, "_metadata");           writeData(wb, "_metadata",           sheet_meta)

# Bold headers on every sheet
hs <- createStyle(textDecoration = "bold")
for (sn in names(wb)) {
  addStyle(wb, sn, hs, rows = 1, cols = seq_len(20), gridExpand = TRUE)
  setColWidths(wb, sn, cols = seq_len(20), widths = "auto")
}

xlsx_path <- file.path(DAT_DIR, "F0_supplementary.xlsx")
saveWorkbook(wb, xlsx_path, overwrite = TRUE)
cat(sprintf("Saved %s (%d sheets)\n", xlsx_path, length(names(wb))))

cat("Figure 0 composite assembly complete.\n")
