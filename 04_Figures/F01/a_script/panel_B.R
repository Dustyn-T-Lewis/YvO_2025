# F01 Panel B: DXA Lean Body Mass
setwd(rprojroot::find_rstudio_root_file())
source("04_Figures/F01/a_script/style.R")

suppressPackageStartupMessages({
  library(readxl)
  library(dplyr)
  library(tidyr)
  library(patchwork)
  library(ggsignif)
  library(rstatix)
})

PW <- 170; PH <- 80
RPT <- "04_Figures/F01/b_reports"
DAT <- "04_Figures/F01/c_data"
dir.create(RPT, recursive = TRUE, showWarnings = FALSE)
dir.create(DAT, recursive = TRUE, showWarnings = FALSE)

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
    Group       = factor(Group, levels = c("Young", "Old")),
    Timepoint   = factor(Timepoint, levels = c("Pre", "Post")),
    Group_Time  = factor(Group_Time,
                         levels = c("Young_Pre", "Young_Post", "Old_Pre", "Old_Post"))
  )

pheno_wide <- meta %>%
  select(subject_key, Group, Timepoint, DXA_LBM_kg) %>%
  pivot_wider(names_from = Timepoint, values_from = DXA_LBM_kg, names_sep = "_") %>%
  rename(DXA_Pre = Pre, DXA_Post = Post) %>%
  mutate(delta_DXA = DXA_Post - DXA_Pre)

dxa_long <- meta %>% select(subject_key, Group, Timepoint, DXA_LBM_kg) %>%
  filter(!is.na(DXA_LBM_kg))
stats_B_anova <- rstatix::anova_test(data = dxa_long, dv = DXA_LBM_kg,
                                      wid = subject_key,
                                      between = Group, within = Timepoint)

dxa_young <- pheno_wide %>% filter(Group == "Young")
dxa_old   <- pheno_wide %>% filter(Group == "Old")
stats_B_paired_young <- t.test(dxa_young$DXA_Post, dxa_young$DXA_Pre, paired = TRUE)
stats_B_paired_old   <- t.test(dxa_old$DXA_Post, dxa_old$DXA_Pre, paired = TRUE)
stats_B_delta        <- t.test(delta_DXA ~ Group, data = pheno_wide)

group_summary_B <- meta %>%
  filter(!is.na(DXA_LBM_kg)) %>%
  group_by(Group, Timepoint, Group_Time) %>%
  summarise(dxa_mean = mean(DXA_LBM_kg), dxa_sem = sd(DXA_LBM_kg) / sqrt(n()),
            .groups = "drop")

anova_tbl <- as.data.frame(stats_B_anova)
anova_sub <- sprintf("Age %s   Time %s   Interaction %s",
                     fmt_p(anova_tbl$p[anova_tbl$Effect == "Group"]),
                     fmt_p(anova_tbl$p[anova_tbl$Effect == "Timepoint"]),
                     fmt_p(anova_tbl$p[anova_tbl$Effect == "Group:Timepoint"]))

# --- Shapiro-Wilk normality ---
sw_yp <- shapiro.test(dxa_long$DXA_LBM_kg[dxa_long$Group == "Young" & dxa_long$Timepoint == "Pre"])
sw_yo <- shapiro.test(dxa_long$DXA_LBM_kg[dxa_long$Group == "Young" & dxa_long$Timepoint == "Post"])
sw_op <- shapiro.test(dxa_long$DXA_LBM_kg[dxa_long$Group == "Old" & dxa_long$Timepoint == "Pre"])
sw_oo <- shapiro.test(dxa_long$DXA_LBM_kg[dxa_long$Group == "Old" & dxa_long$Timepoint == "Post"])
sw_dy <- shapiro.test(dxa_young$delta_DXA)
sw_do <- shapiro.test(dxa_old$delta_DXA)
norm_sub <- sprintf("Shapiro-Wilk (delta): Young %s, Old %s",
                    fmt_p(sw_dy$p.value), fmt_p(sw_do$p.value))
full_sub <- paste0(anova_sub, "\n", norm_sub)

# --- Audit CSV ---
audit_B <- data.frame(
  test = c("paired_t_young", "paired_t_old", "unpaired_t_delta"),
  Group = c("Young", "Old", "Young vs Old"),
  statistic = c(stats_B_paired_young$statistic, stats_B_paired_old$statistic, stats_B_delta$statistic),
  p_value = c(stats_B_paired_young$p.value, stats_B_paired_old$p.value, stats_B_delta$p.value),
  df = c(stats_B_paired_young$parameter, stats_B_paired_old$parameter, stats_B_delta$parameter),
  mean_diff = c(stats_B_paired_young$estimate, stats_B_paired_old$estimate, diff(stats_B_delta$estimate)),
  ci_lo = c(stats_B_paired_young$conf.int[1], stats_B_paired_old$conf.int[1], stats_B_delta$conf.int[1]),
  ci_hi = c(stats_B_paired_young$conf.int[2], stats_B_paired_old$conf.int[2], stats_B_delta$conf.int[2]),
  shapiro_p = c(sw_dy$p.value, sw_do$p.value, NA)
)
write.csv(audit_B, file.path(DAT, "audit_panel_B.csv"), row.names = FALSE)

y_max_left <- max(meta$DXA_LBM_kg, na.rm = TRUE)

pB_left <- ggplot(meta, aes(x = Group_Time, y = DXA_LBM_kg, fill = Group_Time)) +
  geom_bar(stat = "summary", fun = mean, width = 0.65, color = "grey30", linewidth = 0.3) +
  geom_errorbar(stat = "summary", fun.data = mean_se, width = 0.2, linewidth = 0.4) +
  geom_jitter(width = 0.12, size = 1.2, alpha = 0.5, shape = 21, color = "black", stroke = 0.3) +
  geom_signif(comparisons = list(c("Young_Pre", "Young_Post")),
              annotations = fmt_p(stats_B_paired_young$p.value),
              y_position = y_max_left * 1.05, textsize = 2.5, tip_length = 0.01) +
  geom_signif(comparisons = list(c("Old_Pre", "Old_Post")),
              annotations = fmt_p(stats_B_paired_old$p.value),
              y_position = y_max_left * 1.05, textsize = 2.5, tip_length = 0.01) +
  annotate("text", x = 1.5, y = -Inf, label = "Young",
           vjust = 4.2, fontface = "bold", size = 3.2, color = "grey25") +
  annotate("text", x = 3.5, y = -Inf, label = "Old",
           vjust = 4.2, fontface = "bold", size = 3.2, color = "grey25") +
  scale_fill_manual(values = GROUP_FILL) +
  scale_x_discrete(labels = c(Young_Pre = "Pre", Young_Post = "Post",
                               Old_Pre = "Pre", Old_Post = "Post")) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.22))) +
  coord_cartesian(clip = "off") +
  labs(title = "DXA Lean Body Mass", subtitle = full_sub,
       y = "DXA LBM (kg)", x = NULL, tag = "B") +
  FIG_THEME +
  theme(plot.subtitle = element_text(size = 7, color = "grey40", face = "italic"),
        plot.margin = margin(5, 5, 20, 5), legend.position = "none")

delta_bar_colors <- c(Young = unname(GROUP_FILL["Young_Post"]),
                      Old   = unname(GROUP_FILL["Old_Post"]))

y_max_right <- max(pheno_wide$delta_DXA, na.rm = TRUE)

pB_right <- ggplot(pheno_wide, aes(x = Group, y = delta_DXA, fill = Group)) +
  geom_bar(stat = "summary", fun = mean, width = 0.55, color = "grey30", linewidth = 0.3) +
  geom_errorbar(stat = "summary", fun.data = mean_se, width = 0.15, linewidth = 0.4) +
  geom_jitter(width = 0.12, size = 1.2, alpha = 0.5, shape = 21, color = "black", stroke = 0.3) +
  geom_signif(comparisons = list(c("Young", "Old")),
              annotations = fmt_p(stats_B_delta$p.value),
              textsize = 2.5, tip_length = 0.02,
              y_position = y_max_right * 1.10) +
  scale_fill_manual(values = delta_bar_colors) +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.25))) +
  labs(y = "change in DXA LBM (kg)", x = NULL) +
  FIG_THEME + theme(legend.position = "none")

pB <- (pB_left | pB_right) + plot_layout(widths = c(0.65, 0.35))

ggsave(file.path(RPT, "panel_B_dxa_lbm.pdf"), pB,
       width = PW, height = PH, units = "mm", device = get_pdf_device())
ggsave(file.path(RPT, "panel_B_dxa_lbm.png"), pB,
       width = PW, height = PH, units = "mm", dpi = 300)
cat("F01 Panel B done\n")
