# F01 Panel B: DXA Lean Body Mass
setwd(rprojroot::find_rstudio_root_file())
source("04_Figures/F01/a_script/style.R")

library(readxl)
library(dplyr)
library(tidyr)
library(patchwork)
library(ggsignif)
library(rstatix)

PW <- 170; PH <- 80
RPT <- "04_Figures/F01/b_reports/panels"
DAT <- "04_Figures/F01/c_data"
dir.create(RPT, recursive = TRUE, showWarnings = FALSE)
dir.create(DAT, recursive = TRUE, showWarnings = FALSE)

meta <- read_excel("00_input/YvO_meta.xlsx") %>%
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

anova_tbl <- as.data.frame(stats_B_anova)
anova_sub <- sprintf("Age %s   Time %s   Interaction %s",
                     fmt_p(anova_tbl$p[anova_tbl$Effect == "Group"]),
                     fmt_p(anova_tbl$p[anova_tbl$Effect == "Timepoint"]),
                     fmt_p(anova_tbl$p[anova_tbl$Effect == "Group:Timepoint"]))

# --- Shapiro-Wilk normality (on deltas — the tested quantity)
sw_dy <- shapiro.test(dxa_young$delta_DXA)
sw_do <- shapiro.test(dxa_old$delta_DXA)
n_y <- nrow(dxa_young); n_o <- nrow(dxa_old)
norm_sub <- sprintf("n = %d (Y %d, O %d) | Shapiro-Wilk (delta): Y %s, O %s",
                    n_y + n_o, n_y, n_o,
                    fmt_p(sw_dy$p.value), fmt_p(sw_do$p.value))
full_sub <- paste0(anova_sub, "\n", norm_sub)

# --- Audit CSV
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
write.csv(audit_B, file.path(DAT, "panel_B_dxa_lbm.csv"), row.names = FALSE)

y_max_left <- max(meta$DXA_LBM_kg, na.rm = TRUE)

pB_left <- ggplot(meta, aes(x = Group_Time, y = DXA_LBM_kg, fill = Group_Time)) +
  # Background lane shading (F02/F03 style)
  annotate("rect", xmin = 0.5, xmax = 2.5, ymin = -Inf, ymax = Inf,
           fill = AGE_COLORS["Young"], alpha = 0.20, color = "grey85", linewidth = 0.15) +
  annotate("rect", xmin = 2.5, xmax = 4.5, ymin = -Inf, ymax = Inf,
           fill = AGE_COLORS["Old"], alpha = 0.20, color = "grey85", linewidth = 0.15) +
  geom_bar(stat = "summary", fun = mean, width = 0.65, color = "grey30", linewidth = 0.3) +
  geom_errorbar(stat = "summary", fun.data = mean_se, width = 0.2, linewidth = 0.4) +
  geom_jitter(width = 0.12, size = 1.2, alpha = 0.35, shape = 21, color = "black", stroke = 0.3) +
  geom_signif(comparisons = list(c("Young_Pre", "Young_Post")),
              annotations = fmt_p(stats_B_paired_young$p.value),
              y_position = y_max_left * 1.05, textsize = 2.5, tip_length = 0.01) +
  geom_signif(comparisons = list(c("Old_Pre", "Old_Post")),
              annotations = fmt_p(stats_B_paired_old$p.value),
              y_position = y_max_left * 1.05, textsize = 2.5, tip_length = 0.01) +
  annotate("text", x = 1.5, y = -Inf, label = "Young",
           vjust = 3, fontface = "bold", size = 3.2, color = "grey25") +
  annotate("text", x = 3.5, y = -Inf, label = "Old",
           vjust = 3, fontface = "bold", size = 3.2, color = "grey25") +
  scale_fill_manual(values = GROUP_FILL) +
  scale_x_discrete(labels = c(Young_Pre = "Pre", Young_Post = "Post",
                               Old_Pre = "Pre", Old_Post = "Post")) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.18))) +
  coord_cartesian(clip = "off") +
  labs(title = "DXA Lean Body Mass", subtitle = full_sub,
       y = "DXA LBM (kg)", x = NULL, tag = "B") +
  FIG_THEME +
  theme(plot.margin = margin(5, 5, 20, 5), legend.position = "none")

delta_bar_colors <- c(Young = unname(GROUP_FILL["Young_Post"]),
                      Old   = unname(GROUP_FILL["Old_Post"]))

y_max_right <- max(pheno_wide$delta_DXA, na.rm = TRUE)

pB_right <- ggplot(pheno_wide, aes(x = Group, y = delta_DXA, fill = Group)) +
  # Background lane shading (F02/F03 style)
  annotate("rect", xmin = 0.5, xmax = 1.5, ymin = -Inf, ymax = Inf,
           fill = AGE_COLORS["Young"], alpha = 0.20, color = "grey85", linewidth = 0.15) +
  annotate("rect", xmin = 1.5, xmax = 2.5, ymin = -Inf, ymax = Inf,
           fill = AGE_COLORS["Old"], alpha = 0.20, color = "grey85", linewidth = 0.15) +
  geom_bar(stat = "summary", fun = mean, width = 0.55, color = "grey30", linewidth = 0.3) +
  geom_errorbar(stat = "summary", fun.data = mean_se, width = 0.15, linewidth = 0.4) +
  geom_jitter(width = 0.12, size = 1.2, alpha = 0.35, shape = 21, color = "black", stroke = 0.3) +
  geom_signif(comparisons = list(c("Young", "Old")),
              annotations = fmt_p(stats_B_delta$p.value),
              textsize = 2.5, tip_length = 0.02,
              y_position = y_max_right * 1.10) +
  scale_fill_manual(values = delta_bar_colors) +
  scale_y_continuous(expand = expansion(mult = c(0.02, 0.22))) +
  labs(y = "change in DXA LBM (kg)", x = NULL) +
  FIG_THEME + theme(legend.position = "none",
                    plot.margin = margin(5, 5, 20, 5))

pB <- (pB_left | pB_right) + plot_layout(widths = c(0.65, 0.35))

ggsave(file.path(RPT, "MAIN_panel_B_dxa_lbm.png"), pB,
       width = PW, height = PH, units = "mm", dpi = 300)
cat("F01 Panel B done\n")
