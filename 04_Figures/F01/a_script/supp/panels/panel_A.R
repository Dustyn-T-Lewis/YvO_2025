# F01 Supplementary Panel A: Deadlift 1RM
setwd(rprojroot::find_rstudio_root_file())
source("04_Figures/shared/style.R")

library(readxl)
library(dplyr)
library(tidyr)
library(patchwork)
library(ggsignif)
library(rstatix)

PW <- 170; PH <- 80
RPT_PNG <- "04_Figures/F01/b_reports/supp/png/panels"
RPT_PDF <- "04_Figures/F01/b_reports/supp/pdf/panels"
DAT <- "04_Figures/F01/c_data/supp"
dir.create(RPT_PNG, recursive = TRUE, showWarnings = FALSE)
dir.create(RPT_PDF, recursive = TRUE, showWarnings = FALSE)
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

pheno_long <- meta %>%
  select(subject_key, Group, Timepoint, deadlift_1rm_kg) %>%
  filter(!is.na(deadlift_1rm_kg))

pheno_wide <- meta %>%
  select(subject_key, Group, Timepoint, deadlift_1rm_kg) %>%
  pivot_wider(names_from = Timepoint, values_from = deadlift_1rm_kg, names_sep = "_") %>%
  rename(DL_Pre = Pre, DL_Post = Post) %>%
  filter(!is.na(DL_Pre), !is.na(DL_Post)) %>%
  mutate(delta_DL = DL_Post - DL_Pre)

stats_anova <- rstatix::anova_test(data = pheno_long, dv = deadlift_1rm_kg,
                                    wid = subject_key,
                                    between = Group, within = Timepoint)

dl_young <- pheno_wide %>% filter(Group == "Young")
dl_old   <- pheno_wide %>% filter(Group == "Old")
stats_paired_young <- t.test(dl_young$DL_Post, dl_young$DL_Pre, paired = TRUE)
stats_paired_old   <- t.test(dl_old$DL_Post, dl_old$DL_Pre, paired = TRUE)
stats_delta        <- t.test(delta_DL ~ Group, data = pheno_wide)

anova_tbl <- as.data.frame(stats_anova)
anova_sub <- sprintf("Age %s   Time %s   Interaction %s",
                     fmt_p(anova_tbl$p[anova_tbl$Effect == "Group"]),
                     fmt_p(anova_tbl$p[anova_tbl$Effect == "Timepoint"]),
                     fmt_p(anova_tbl$p[anova_tbl$Effect == "Group:Timepoint"]))

sw_dy <- shapiro.test(dl_young$delta_DL)
sw_do <- shapiro.test(dl_old$delta_DL)
n_y <- sum(pheno_wide$Group == "Young")
n_o <- sum(pheno_wide$Group == "Old")
norm_sub <- sprintf("n = %d (Y %d, O %d) | Shapiro-Wilk (delta): Y %s, O %s",
                    n_y + n_o, n_y, n_o,
                    fmt_p(sw_dy$p.value), fmt_p(sw_do$p.value))
full_sub <- paste0(anova_sub, "\n", norm_sub)

audit_sA <- data.frame(
  test = c("paired_t_young", "paired_t_old", "unpaired_t_delta"),
  Group = c("Young", "Old", "Young vs Old"),
  statistic = c(stats_paired_young$statistic, stats_paired_old$statistic, stats_delta$statistic),
  p_value = c(stats_paired_young$p.value, stats_paired_old$p.value, stats_delta$p.value),
  df = c(stats_paired_young$parameter, stats_paired_old$parameter, stats_delta$parameter),
  mean_diff = c(stats_paired_young$estimate, stats_paired_old$estimate, diff(stats_delta$estimate)),
  ci_lo = c(stats_paired_young$conf.int[1], stats_paired_old$conf.int[1], stats_delta$conf.int[1]),
  ci_hi = c(stats_paired_young$conf.int[2], stats_paired_old$conf.int[2], stats_delta$conf.int[2]),
  shapiro_p = c(sw_dy$p.value, sw_do$p.value, NA)
)
write.csv(audit_sA, file.path(DAT, "panel_A_deadlift_1rm.csv"), row.names = FALSE)

# Filter long-form data to complete-pair subjects only
plot_long <- meta %>%
  filter(subject_key %in% pheno_wide$subject_key, !is.na(deadlift_1rm_kg))

y_max_left <- max(plot_long$deadlift_1rm_kg, na.rm = TRUE)

pSA_left <- ggplot(plot_long, aes(x = Group_Time, y = deadlift_1rm_kg, fill = Group_Time)) +
  annotate("rect", xmin = 0.5, xmax = 2.5, ymin = -Inf, ymax = Inf,
           fill = AGE_COLORS["Young"], alpha = 0.20, color = "grey85", linewidth = 0.15) +
  annotate("rect", xmin = 2.5, xmax = 4.5, ymin = -Inf, ymax = Inf,
           fill = AGE_COLORS["Old"], alpha = 0.20, color = "grey85", linewidth = 0.15) +
  geom_bar(stat = "summary", fun = mean, width = 0.65, color = "grey30", linewidth = 0.3) +
  geom_errorbar(stat = "summary", fun.data = mean_se, width = 0.2, linewidth = 0.4) +
  geom_jitter(width = 0.12, size = 1.2, alpha = 0.35, shape = 21, color = "black", stroke = 0.3) +
  geom_signif(comparisons = list(c("Young_Pre", "Young_Post")),
              annotations = fmt_p(stats_paired_young$p.value),
              y_position = y_max_left * 1.05, textsize = 2.5, tip_length = 0.01) +
  geom_signif(comparisons = list(c("Old_Pre", "Old_Post")),
              annotations = fmt_p(stats_paired_old$p.value),
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
  labs(title = "Deadlift 1RM", subtitle = full_sub,
       y = "Deadlift 1RM (kg)", x = NULL, tag = "A") +
  FIG_THEME +
  theme(plot.subtitle = element_text(size = 7, color = "grey40", face = "italic"),  # 7pt: supplementary compact panel
        plot.margin = margin(5, 5, 20, 5), legend.position = "none")

delta_bar_colors <- c(Young = unname(GROUP_FILL["Young_Post"]),
                      Old   = unname(GROUP_FILL["Old_Post"]))

y_max_right <- max(pheno_wide$delta_DL, na.rm = TRUE)

pSA_right <- ggplot(pheno_wide, aes(x = Group, y = delta_DL, fill = Group)) +
  annotate("rect", xmin = 0.5, xmax = 1.5, ymin = -Inf, ymax = Inf,
           fill = AGE_COLORS["Young"], alpha = 0.20, color = "grey85", linewidth = 0.15) +
  annotate("rect", xmin = 1.5, xmax = 2.5, ymin = -Inf, ymax = Inf,
           fill = AGE_COLORS["Old"], alpha = 0.20, color = "grey85", linewidth = 0.15) +
  geom_bar(stat = "summary", fun = mean, width = 0.55, color = "grey30", linewidth = 0.3) +
  geom_errorbar(stat = "summary", fun.data = mean_se, width = 0.15, linewidth = 0.4) +
  geom_jitter(width = 0.12, size = 1.2, alpha = 0.35, shape = 21, color = "black", stroke = 0.3) +
  geom_signif(comparisons = list(c("Young", "Old")),
              annotations = fmt_p(stats_delta$p.value),
              textsize = 2.5, tip_length = 0.02,
              y_position = y_max_right * 1.10) +
  scale_fill_manual(values = delta_bar_colors) +
  scale_y_continuous(expand = expansion(mult = c(0.02, 0.22))) +
  labs(y = "change in 1RM (kg)", x = NULL) +
  FIG_THEME + theme(legend.position = "none",
                    plot.margin = margin(5, 5, 20, 5))

pSA <- (pSA_left | pSA_right) + plot_layout(widths = c(0.65, 0.35))

ggsave(file.path(RPT_PNG, "SUPP_panel_A_deadlift_1rm.png"), pSA,
       width = PW, height = PH, units = "mm", dpi = 300)
ggsave(file.path(RPT_PDF, "SUPP_panel_A_deadlift_1rm.pdf"), pSA,
       width = PW, height = PH, units = "mm", device = get_pdf_device())
cat("F01 Supp Panel A done\n")
