# F01 pre/post panel template — sourced by scripts that define a `cfg` list.
#
# Required cfg: dv_col, y_label, delta_label, title, tag, output_prefix,
#   file_tag, audit_file, rpt_png, rpt_pdf, dat, file_prefix
# Optional cfg: y_breaks, y_labels, coerce_cols, filter_complete,
#   left_margin, right_margin, left_x_expand, use_plotmath_subtitle

stopifnot(exists("cfg"), is.list(cfg))

source(here::here("04_Figures", "shared", "style.R"))

library(readxl)
library(dplyr)
library(tidyr)
library(patchwork)
library(ggsignif)
library(rstatix)

if (is.null(cfg$coerce_cols))           cfg$coerce_cols   <- FALSE
if (is.null(cfg$filter_complete))       cfg$filter_complete <- FALSE
if (is.null(cfg$left_margin))           cfg$left_margin   <- margin(2, 2, 2, 2)
if (is.null(cfg$right_margin))          cfg$right_margin  <- margin(2, 2, 2, 2)
if (is.null(cfg$left_x_expand))         cfg$left_x_expand <- expansion(add = 0.3)
if (is.null(cfg$use_plotmath_subtitle)) cfg$use_plotmath_subtitle <- FALSE

PW <- 170; PH <- 80
for (d in c(cfg$rpt_png, cfg$rpt_pdf, cfg$dat))
  dir.create(d, recursive = TRUE, showWarnings = FALSE)

meta <- read_excel(here::here("00_input", "YvO_meta.xlsx"))

if (cfg$coerce_cols) {
  for (col in c("BMI", "Type_I_fCSA", "Type_II_fCSA",
                "deadlift_1rm_kg", "Total_Training_Volume_kg"))
    if (col %in% names(meta) && is.character(meta[[col]]))
      meta[[col]] <- suppressWarnings(as.numeric(meta[[col]]))
}

meta <- meta |>
  mutate(subject_key = sub("_(Pre|Post)$", "", Col_ID),
         Group      = factor(Group, levels = c("Young", "Old")),
         Timepoint  = factor(Timepoint, levels = c("Pre", "Post")),
         Group_Time = factor(Group_Time,
           levels = c("Young_Pre", "Young_Post", "Old_Pre", "Old_Post")))

DV <- cfg$dv_col

pheno_wide <- meta |>
  select(subject_key, Group, Timepoint, all_of(DV)) |>
  pivot_wider(names_from = Timepoint, values_from = all_of(DV)) |>
  rename(DV_Pre = Pre, DV_Post = Post) |>
  filter(!is.na(DV_Pre), !is.na(DV_Post)) |>
  mutate(delta_DV = DV_Post - DV_Pre)

pheno_long <- meta |>
  select(subject_key, Group, Timepoint, all_of(DV)) |>
  filter(!is.na(.data[[DV]]), subject_key %in% pheno_wide$subject_key)

stats_anova <- anova_test(pheno_long, dv = !!sym(DV),
                          wid = subject_key, between = Group, within = Timepoint)

grp_young <- pheno_wide |> filter(Group == "Young")
grp_old   <- pheno_wide |> filter(Group == "Old")
stats_paired_young <- t.test(grp_young$DV_Post, grp_young$DV_Pre, paired = TRUE)
stats_paired_old   <- t.test(grp_old$DV_Post,   grp_old$DV_Pre,   paired = TRUE)
stats_delta        <- t.test(delta_DV ~ Group, data = pheno_wide)

anova_tbl <- as.data.frame(stats_anova)
anova_sub <- sprintf("RM-ANOVA: Age %s   Time %s   Int. %s",
                     fmt_p(anova_tbl$p[anova_tbl$Effect == "Group"]),
                     fmt_p(anova_tbl$p[anova_tbl$Effect == "Timepoint"]),
                     fmt_p(anova_tbl$p[anova_tbl$Effect == "Group:Timepoint"]))

sw_dy <- shapiro.test(grp_young$delta_DV)
sw_do <- shapiro.test(grp_old$delta_DV)

audit_df <- data.frame(
  test      = c("paired_t_young", "paired_t_old", "unpaired_t_delta"),
  Group     = c("Young", "Old", "Young vs Old"),
  n         = c(nrow(grp_young), nrow(grp_old), nrow(pheno_wide)),
  statistic = c(stats_paired_young$statistic, stats_paired_old$statistic, stats_delta$statistic),
  p_value   = c(stats_paired_young$p.value,   stats_paired_old$p.value,   stats_delta$p.value),
  df        = c(stats_paired_young$parameter,  stats_paired_old$parameter,  stats_delta$parameter),
  mean_diff = c(stats_paired_young$estimate,   stats_paired_old$estimate,   diff(stats_delta$estimate)),
  ci_lo     = c(stats_paired_young$conf.int[1], stats_paired_old$conf.int[1], stats_delta$conf.int[1]),
  ci_hi     = c(stats_paired_young$conf.int[2], stats_paired_old$conf.int[2], stats_delta$conf.int[2]),
  shapiro_p = c(sw_dy$p.value, sw_do$p.value, NA))
write.csv(audit_df, file.path(cfg$dat, cfg$audit_file), row.names = FALSE)

plot_long <- if (cfg$filter_complete) {
  meta |> filter(subject_key %in% pheno_wide$subject_key, !is.na(.data[[DV]]))
} else meta

y_max_left <- max(plot_long[[DV]], na.rm = TRUE)

p_left <- ggplot(plot_long, aes(Group_Time, .data[[DV]], fill = Group_Time)) +
  annotate("rect", xmin = 0.5, xmax = 2.5, ymin = -Inf, ymax = Inf,
           fill = AGE_COLORS["Young"], alpha = 0.20, color = "grey85", linewidth = 0.15) +
  annotate("rect", xmin = 2.5, xmax = 4.5, ymin = -Inf, ymax = Inf,
           fill = AGE_COLORS["Old"], alpha = 0.20, color = "grey85", linewidth = 0.15) +
  geom_bar(stat = "summary", fun = mean, width = 0.65, color = "grey30", linewidth = 0.3) +
  geom_errorbar(stat = "summary", fun.data = mean_se, width = 0.2, linewidth = 0.4) +
  geom_jitter(width = 0.12, size = 0.8, alpha = 0.35, shape = 21, color = "black", stroke = 0.2) +
  geom_signif(comparisons = list(c("Young_Pre", "Young_Post")),
              annotations = fmt_p_plot(stats_paired_young$p.value),
              parse = TRUE, y_position = y_max_left * 1.05,
              textsize = 1.5, size = 0.3, tip_length = 0.01) +
  geom_signif(comparisons = list(c("Old_Pre", "Old_Post")),
              annotations = fmt_p_plot(stats_paired_old$p.value),
              parse = TRUE, y_position = y_max_left * 1.05,
              textsize = 1.5, size = 0.3, tip_length = 0.01) +
  scale_fill_manual(values = GROUP_FILL) +
  scale_x_discrete(labels = c(Young_Pre = "Pre", Young_Post = "Post",
                               Old_Pre = "Pre", Old_Post = "Post"),
                   expand = cfg$left_x_expand) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.18))) +
  labs(title = cfg$title, subtitle = anova_sub,
       y = cfg$y_label, x = NULL, tag = cfg$tag) +
  FIG_THEME + theme(plot.margin = cfg$left_margin, legend.position = "none")

delta_bar_colors <- c(Young = unname(GROUP_FILL["Young_Post"]),
                      Old   = unname(GROUP_FILL["Old_Post"]))
y_max_right <- max(pheno_wide$delta_DV, na.rm = TRUE)

p_right <- ggplot(pheno_wide, aes(Group, delta_DV, fill = Group)) +
  annotate("rect", xmin = 0.5, xmax = 1.5, ymin = -Inf, ymax = Inf,
           fill = AGE_COLORS["Young"], alpha = 0.20, color = "grey85", linewidth = 0.15) +
  annotate("rect", xmin = 1.5, xmax = 2.5, ymin = -Inf, ymax = Inf,
           fill = AGE_COLORS["Old"], alpha = 0.20, color = "grey85", linewidth = 0.15) +
  geom_bar(stat = "summary", fun = mean, width = 0.55, color = "grey30", linewidth = 0.3) +
  geom_errorbar(stat = "summary", fun.data = mean_se, width = 0.15, linewidth = 0.4) +
  geom_jitter(width = 0.12, size = 0.8, alpha = 0.35, shape = 21, color = "black", stroke = 0.2) +
  geom_signif(comparisons = list(c("Young", "Old")),
              annotations = fmt_p_plot(stats_delta$p.value),
              parse = TRUE, textsize = 1.5, size = 0.3, tip_length = 0.02,
              y_position = y_max_right * 1.10) +
  scale_fill_manual(values = delta_bar_colors) +
  scale_x_discrete(expand = expansion(add = 0.3)) +
  labs(y = cfg$delta_label, x = NULL) +
  FIG_THEME + theme(legend.position = "none",
                    axis.title.y = element_text(margin = margin(r = 1)),
                    plot.margin = cfg$right_margin)

if (!is.null(cfg$y_breaks)) {
  p_right <- p_right +
    scale_y_continuous(expand = expansion(mult = c(0.02, 0.22)),
                       breaks = cfg$y_breaks, labels = cfg$y_labels)
} else {
  p_right <- p_right + scale_y_continuous(expand = expansion(mult = c(0.02, 0.22)))
}

p_combo <- (p_left | p_right) + plot_layout(widths = c(0.65, 0.35))
ggsave(file.path(cfg$rpt_png, paste0(cfg$file_prefix, "_", cfg$file_tag, ".png")),
       p_combo, width = PW, height = PH, units = "mm", dpi = 300)
ggsave(file.path(cfg$rpt_pdf, paste0(cfg$file_prefix, "_", cfg$file_tag, ".pdf")),
       p_combo, width = PW, height = PH, units = "mm", device = get_pdf_device())
message(sprintf("F01 %s done", cfg$title))

pfx <- cfg$output_prefix
assign(paste0(pfx, "_title"), cfg$title, envir = .GlobalEnv)
assign(paste0(pfx, "_subtitle"),
       if (cfg$use_plotmath_subtitle) {
         fmt_anova_sub(anova_tbl$p[anova_tbl$Effect == "Group"],
                       anova_tbl$p[anova_tbl$Effect == "Timepoint"],
                       anova_tbl$p[anova_tbl$Effect == "Group:Timepoint"])
       } else anova_sub,
       envir = .GlobalEnv)
assign(paste0(pfx, "_left"),  strip_for_composite(p_left),  envir = .GlobalEnv)
assign(paste0(pfx, "_right"), strip_for_composite(p_right), envir = .GlobalEnv)
