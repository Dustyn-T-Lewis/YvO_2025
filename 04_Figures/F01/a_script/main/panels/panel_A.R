# F01 Panel A: Training Volume
setwd(rprojroot::find_rstudio_root_file())
source("04_Figures/shared/style.R")

library(readxl)
library(dplyr)
library(ggsignif)

PW <- 90; PH <- 150
RPT_PNG <- "04_Figures/F01/b_reports/main/png/panels"
RPT_PDF <- "04_Figures/F01/b_reports/main/pdf/panels"
DAT <- "04_Figures/F01/c_data"
dir.create(RPT_PNG, recursive = TRUE, showWarnings = FALSE)
dir.create(RPT_PDF, recursive = TRUE, showWarnings = FALSE)
dir.create(DAT, recursive = TRUE, showWarnings = FALSE)

stopifnot("Input metadata missing: 00_input/YvO_meta.xlsx" = file.exists("00_input/YvO_meta.xlsx"))
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
  )

tv_df <- meta %>%
  filter(Timepoint == "Post") %>%
  select(subject_key, Group, Total_Training_Volume_kg) %>%
  rename(tv = Total_Training_Volume_kg) %>%
  filter(!is.na(tv))

stats_A <- t.test(tv ~ Group, data = tv_df)

sw_young <- shapiro.test(tv_df$tv[tv_df$Group == "Young"])
sw_old   <- shapiro.test(tv_df$tv[tv_df$Group == "Old"])
n_y <- sum(tv_df$Group == "Young"); n_o <- sum(tv_df$Group == "Old")
norm_sub <- sprintf("Welch t %s", fmt_p(stats_A$p.value))

# --- Audit CSV
audit_A <- tv_df %>%
  group_by(Group) %>%
  summarise(n = n(), mean = mean(tv), sd = sd(tv), sem = sd(tv) / sqrt(n()),
            .groups = "drop") %>%
  mutate(shapiro_p = c(sw_young$p.value, sw_old$p.value),
         t_test_p  = stats_A$p.value)
write.csv(audit_A, file.path(DAT, "panel_A_training_volume.csv"), row.names = FALSE)

bar_colors <- c(Young = unname(GROUP_FILL["Young_Post"]),
                Old   = unname(GROUP_FILL["Old_Post"]))

tv_df$tv_scaled <- tv_df$tv / 1e5

pA <- ggplot(tv_df, aes(x = Group, y = tv_scaled, fill = Group)) +
  annotate("rect", xmin = 0.5, xmax = 1.5, ymin = -Inf, ymax = Inf,
           fill = AGE_COLORS["Young"], alpha = 0.20, color = "grey85", linewidth = 0.15) +
  annotate("rect", xmin = 1.5, xmax = 2.5, ymin = -Inf, ymax = Inf,
           fill = AGE_COLORS["Old"], alpha = 0.20, color = "grey85", linewidth = 0.15) +
  geom_bar(stat = "summary", fun = mean, width = 0.6, color = "grey30", linewidth = 0.3) +
  geom_errorbar(stat = "summary", fun.data = mean_se, width = 0.2, linewidth = 0.4) +
  geom_jitter(width = 0.15, size = 1, alpha = 0.35, shape = 16, color = "grey30") +
  geom_signif(
    comparisons = list(c("Young", "Old")),
    annotations = fmt_p_plot(stats_A$p.value),
    parse = TRUE, textsize = 1.5, size = 0.3, tip_length = 0.02,
    y_position = max(tv_df$tv_scaled) * 1.15
  ) +
  scale_fill_manual(values = bar_colors) +
  scale_x_discrete(labels = c(Young = "Young", Old = "Old"),
                   expand = expansion(add = 0.3)) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.18))) +
  labs(title = "Training Volume", subtitle = norm_sub,
       y = expression(bold("Volume (" %*% 10^5 ~ "kg)")), x = NULL,
       tag = "a") +
  FIG_THEME + theme(legend.position = "none",
                    plot.margin = margin(2, 2, 2, 2))

ggsave(file.path(RPT_PNG, "MAIN_panel_A_training_volume.png"), pA,
       width = PW, height = PH, units = "mm", dpi = 300)
ggsave(file.path(RPT_PDF, "MAIN_panel_A_training_volume.pdf"), pA,
       width = PW, height = PH, units = "mm", device = get_pdf_device())
cat("F01 Panel A done\n")

# --- Export for composite ---
pA_title    <- "Training Volume"
pA_subtitle <- if (stats_A$p.value < 0.05) {
  paste0('bold(italic("', norm_sub, '"))')
} else {
  paste0('italic("', norm_sub, '")')
}
pA_legend   <- NULL
pA          <- strip_for_composite(pA)
