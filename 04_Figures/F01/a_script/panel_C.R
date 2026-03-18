# F01 Panel C: Vastus Lateralis Thickness
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
  select(subject_key, Group, Timepoint, VL_thick_cm) %>%
  pivot_wider(names_from = Timepoint, values_from = VL_thick_cm, names_sep = "_") %>%
  rename(VL_Pre = Pre, VL_Post = Post) %>%
  mutate(delta_VL = VL_Post - VL_Pre)

vl_long <- meta %>% select(subject_key, Group, Timepoint, VL_thick_cm) %>%
  filter(!is.na(VL_thick_cm))
stats_C_anova <- rstatix::anova_test(data = vl_long, dv = VL_thick_cm,
                                      wid = subject_key,
                                      between = Group, within = Timepoint)

vl_young <- pheno_wide %>% filter(Group == "Young")
vl_old   <- pheno_wide %>% filter(Group == "Old")
stats_C_paired_young <- t.test(vl_young$VL_Post, vl_young$VL_Pre, paired = TRUE)
stats_C_paired_old   <- t.test(vl_old$VL_Post, vl_old$VL_Pre, paired = TRUE)
stats_C_delta        <- t.test(delta_VL ~ Group, data = pheno_wide)

group_summary_C <- meta %>%
  filter(!is.na(VL_thick_cm)) %>%
  group_by(Group, Timepoint, Group_Time) %>%
  summarise(vl_mean = mean(VL_thick_cm), vl_sem = sd(VL_thick_cm) / sqrt(n()),
            .groups = "drop")

anova_tbl <- as.data.frame(stats_C_anova)
anova_sub <- sprintf("Age %s   Time %s   Interaction %s",
                     fmt_p(anova_tbl$p[anova_tbl$Effect == "Group"]),
                     fmt_p(anova_tbl$p[anova_tbl$Effect == "Timepoint"]),
                     fmt_p(anova_tbl$p[anova_tbl$Effect == "Group:Timepoint"]))

# --- Shapiro-Wilk normality ---
sw_yp <- shapiro.test(vl_long$VL_thick_cm[vl_long$Group == "Young" & vl_long$Timepoint == "Pre"])
sw_yo <- shapiro.test(vl_long$VL_thick_cm[vl_long$Group == "Young" & vl_long$Timepoint == "Post"])
sw_op <- shapiro.test(vl_long$VL_thick_cm[vl_long$Group == "Old" & vl_long$Timepoint == "Pre"])
sw_oo <- shapiro.test(vl_long$VL_thick_cm[vl_long$Group == "Old" & vl_long$Timepoint == "Post"])
sw_dy <- shapiro.test(vl_young$delta_VL)
sw_do <- shapiro.test(vl_old$delta_VL)
norm_sub <- sprintf("Shapiro-Wilk (delta): Young %s, Old %s",
                    fmt_p(sw_dy$p.value), fmt_p(sw_do$p.value))
full_sub <- paste0(anova_sub, "\n", norm_sub)

# --- Audit CSV ---
audit_C <- data.frame(
  test = c("paired_t_young", "paired_t_old", "unpaired_t_delta"),
  Group = c("Young", "Old", "Young vs Old"),
  statistic = c(stats_C_paired_young$statistic, stats_C_paired_old$statistic, stats_C_delta$statistic),
  p_value = c(stats_C_paired_young$p.value, stats_C_paired_old$p.value, stats_C_delta$p.value),
  df = c(stats_C_paired_young$parameter, stats_C_paired_old$parameter, stats_C_delta$parameter),
  mean_diff = c(stats_C_paired_young$estimate, stats_C_paired_old$estimate, diff(stats_C_delta$estimate)),
  ci_lo = c(stats_C_paired_young$conf.int[1], stats_C_paired_old$conf.int[1], stats_C_delta$conf.int[1]),
  ci_hi = c(stats_C_paired_young$conf.int[2], stats_C_paired_old$conf.int[2], stats_C_delta$conf.int[2]),
  shapiro_p = c(sw_dy$p.value, sw_do$p.value, NA)
)
write.csv(audit_C, file.path(DAT, "audit_panel_C.csv"), row.names = FALSE)

y_max_left <- max(meta$VL_thick_cm, na.rm = TRUE)

pC_left <- ggplot(meta, aes(x = Group_Time, y = VL_thick_cm, fill = Group_Time)) +
  geom_bar(stat = "summary", fun = mean, width = 0.65, color = "grey30", linewidth = 0.3) +
  geom_errorbar(stat = "summary", fun.data = mean_se, width = 0.2, linewidth = 0.4) +
  geom_jitter(width = 0.12, size = 1.2, alpha = 0.5, shape = 21, color = "black", stroke = 0.3) +
  geom_signif(comparisons = list(c("Young_Pre", "Young_Post")),
              annotations = fmt_p(stats_C_paired_young$p.value),
              y_position = y_max_left * 1.05, textsize = 2.5, tip_length = 0.01) +
  geom_signif(comparisons = list(c("Old_Pre", "Old_Post")),
              annotations = fmt_p(stats_C_paired_old$p.value),
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
  labs(title = "VL Thickness", subtitle = full_sub,
       y = "VL thickness (cm)", x = NULL, tag = "C") +
  FIG_THEME +
  theme(plot.subtitle = element_text(size = 7, color = "grey40", face = "italic"),
        plot.margin = margin(5, 5, 20, 5), legend.position = "none")

delta_bar_colors <- c(Young = unname(GROUP_FILL["Young_Post"]),
                      Old   = unname(GROUP_FILL["Old_Post"]))

y_max_right <- max(pheno_wide$delta_VL, na.rm = TRUE)

pC_right <- ggplot(pheno_wide, aes(x = Group, y = delta_VL, fill = Group)) +
  geom_bar(stat = "summary", fun = mean, width = 0.55, color = "grey30", linewidth = 0.3) +
  geom_errorbar(stat = "summary", fun.data = mean_se, width = 0.15, linewidth = 0.4) +
  geom_jitter(width = 0.12, size = 1.2, alpha = 0.5, shape = 21, color = "black", stroke = 0.3) +
  geom_signif(comparisons = list(c("Young", "Old")),
              annotations = fmt_p(stats_C_delta$p.value),
              textsize = 2.5, tip_length = 0.02,
              y_position = y_max_right * 1.10) +
  scale_fill_manual(values = delta_bar_colors) +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.25))) +
  labs(y = "change in VL thickness (cm)", x = NULL) +
  FIG_THEME + theme(legend.position = "none")

pC <- (pC_left | pC_right) + plot_layout(widths = c(0.65, 0.35))

ggsave(file.path(RPT, "panel_C_vl_thickness.pdf"), pC,
       width = PW, height = PH, units = "mm", device = get_pdf_device())
ggsave(file.path(RPT, "panel_C_vl_thickness.png"), pC,
       width = PW, height = PH, units = "mm", dpi = 300)
cat("F01 Panel C done\n")
