# ── F0 Panel A: Training Volume ───────────────────────────────────────────────
setwd(rprojroot::find_rstudio_root_file())
source("04_Figures_v2/shared/style.R")

suppressPackageStartupMessages({
  library(readxl)
  library(dplyr)
  library(ggsignif)
})

# Panel dimensions
PW <- 70
PH <- 80

# Output paths
RPT <- "04_Figures_v2/F0/b_reports"
dir.create(RPT, recursive = TRUE, showWarnings = FALSE)

# ── Data ──────────────────────────────────────────────────────────────────────
meta <- read_excel("00_input/YvO_meta.xlsx")

char_to_num <- c("BMI", "Type_I_fCSA", "Type_II_fCSA",
                 "deadlift_1rm_kg", "Total_Training_Volume_kg")
for (col in char_to_num) {
  if (col %in% names(meta) && is.character(meta[[col]]))
    meta[[col]] <- suppressWarnings(as.numeric(meta[[col]]))
}

meta <- meta %>%
  mutate(
    subject_key = sub("_(Pre|Post)$", "", Col_ID),
    Group = factor(Group, levels = c("Young", "Old")),
    Timepoint = factor(Timepoint, levels = c("Pre", "Post"))
  ) %>%
  filter(subject_key != "Y_S05")

tv_df <- meta %>%
  filter(Timepoint == "Post") %>%
  select(subject_key, Group, Total_Training_Volume_kg) %>%
  rename(tv = Total_Training_Volume_kg) %>%
  filter(!is.na(tv))

stats_A <- t.test(tv ~ Group, data = tv_df)

# ── Plot ──────────────────────────────────────────────────────────────────────
bar_colors <- c(Young = unname(GROUP_FILL["Young_Post"]),
                Old   = unname(GROUP_FILL["Old_Post"]))

pA <- ggplot(tv_df, aes(x = Group, y = tv, fill = Group)) +
  geom_bar(stat = "summary", fun = mean, width = 0.6, color = "grey30", linewidth = 0.3) +
  geom_errorbar(stat = "summary", fun.data = mean_se, width = 0.2, linewidth = 0.4) +
  geom_jitter(width = 0.15, size = 1.5, alpha = 0.5, shape = 16, color = "grey30") +
  geom_signif(
    comparisons = list(c("Young", "Old")),
    annotations = fmt_p(stats_A$p.value),
    textsize = 3, tip_length = 0.02,
    y_position = max(tv_df$tv) * 1.15
  ) +
  scale_fill_manual(values = bar_colors) +
  scale_x_discrete(labels = c(
    Young = sprintf("Younger (n = %d)", sum(tv_df$Group == "Young")),
    Old   = sprintf("Older (n = %d)",   sum(tv_df$Group == "Old")))) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.15)),
                     labels = scales::label_comma()) +
  labs(title = "Training Volume", y = "Total training volume (kg)", x = NULL,
       tag = "A") +
  FIG_THEME + theme(legend.position = "none")

# ── Save ──────────────────────────────────────────────────────────────────────
ggsave(file.path(RPT, "panel_A_training_volume.pdf"), pA,
       width = PW, height = PH, units = "mm", device = get_pdf_device())
ggsave(file.path(RPT, "panel_A_training_volume.png"), pA,
       width = PW, height = PH, units = "mm", dpi = 300)
cat("F0 Panel A done\n")
