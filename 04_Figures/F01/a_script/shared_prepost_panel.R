# ──────────────────────────────────────────────────────────────────────────────
# F01 shared pre/post panel template
# ──────────────────────────────────────────────────────────────────────────────
# Sourced by thin wrapper scripts that define a `cfg` list before sourcing.
#
# Required cfg fields:
#   dv_col          character  Column name in YvO_meta.xlsx (e.g. "DXA_LBM_kg")
#   y_label         expression or character  Y-axis label for left panel
#   delta_label     expression or character  Y-axis label for delta (right) panel
#   title           character  Plot title
#   tag             character  Panel letter tag (e.g. "b", "c")
#   output_prefix   character  Prefix for variable names exported to parent
#                              (e.g. "pB", "pC", "pSA", "pSB", "pSC")
#   file_tag        character  Filename label (e.g. "panel_B_dxa_lbm")
#   audit_file      character  Audit CSV filename (e.g. "panel_B_dxa_lbm.csv")
#   rpt_png         character  Output directory for PNG panels
#   rpt_pdf         character  Output directory for PDF panels
#   dat             character  Output directory for data CSVs
#   file_prefix     character  "MAIN" or "SUPP"
#
# Optional cfg fields:
#   y_breaks        numeric vector  Custom breaks for delta y-axis
#   y_labels        character vector  Custom labels for delta y-axis
#   coerce_cols     logical  Whether to coerce char_to_num columns (default FALSE)
#   filter_complete logical  Whether to filter long-form plot data to complete
#                            pairs only (default FALSE)
#   left_margin     margin() object  Custom plot.margin for left panel
#   right_margin    margin() object  Custom plot.margin for right panel
#   left_x_expand   expansion()  Custom x-axis expansion for left panel
#   use_plotmath_subtitle  logical  If TRUE, export fmt_anova_sub() for subtitle;
#                                   if FALSE, export plain text (default FALSE)
# ──────────────────────────────────────────────────────────────────────────────

stopifnot(exists("cfg"), is.list(cfg))

# --- Setup (load ggplot2 via style.R before defaults that use margin/expansion) ----
setwd(rprojroot::find_rstudio_root_file())
source("04_Figures/shared/style.R")

# --- Defaults for optional fields ----
if (is.null(cfg$coerce_cols))            cfg$coerce_cols   <- FALSE
if (is.null(cfg$filter_complete))        cfg$filter_complete <- FALSE
if (is.null(cfg$left_margin))            cfg$left_margin   <- margin(2, 2, 2, 2)
if (is.null(cfg$right_margin))           cfg$right_margin  <- margin(2, 2, 2, 2)
if (is.null(cfg$left_x_expand))          cfg$left_x_expand <- expansion(add = 0.3)
if (is.null(cfg$use_plotmath_subtitle))  cfg$use_plotmath_subtitle <- FALSE

suppressPackageStartupMessages({
  library(readxl)
  library(dplyr)
  library(tidyr)
  library(patchwork)
  library(ggsignif)
  library(rstatix)
})

PW <- 170; PH <- 80
dir.create(cfg$rpt_png, recursive = TRUE, showWarnings = FALSE)
dir.create(cfg$rpt_pdf, recursive = TRUE, showWarnings = FALSE)
dir.create(cfg$dat,     recursive = TRUE, showWarnings = FALSE)

# --- Load metadata ----
stopifnot("Input metadata missing: 00_input/YvO_meta.xlsx" = file.exists("00_input/YvO_meta.xlsx"))
meta <- read_excel("00_input/YvO_meta.xlsx")

if (cfg$coerce_cols) {
  char_to_num <- c("BMI", "Type_I_fCSA", "Type_II_fCSA",
                   "deadlift_1rm_kg", "Total_Training_Volume_kg")
  for (col in char_to_num) {
    if (col %in% names(meta) && is.character(meta[[col]]))
      meta[[col]] <- suppressWarnings(as.numeric(meta[[col]]))
  }
}

meta <- meta %>%
  mutate(
    subject_key = sub("_(Pre|Post)$", "", Col_ID),
    Group       = factor(Group, levels = c("Young", "Old")),
    Timepoint   = factor(Timepoint, levels = c("Pre", "Post")),
    Group_Time  = factor(Group_Time,
                         levels = c("Young_Pre", "Young_Post", "Old_Pre", "Old_Post"))
  )

DV <- cfg$dv_col

# --- Pivot wide & compute delta ----
pheno_wide <- meta %>%
  select(subject_key, Group, Timepoint, all_of(DV)) %>%
  pivot_wider(names_from = Timepoint, values_from = all_of(DV), names_sep = "_") %>%
  rename(DV_Pre = Pre, DV_Post = Post) %>%
  filter(!is.na(DV_Pre), !is.na(DV_Post)) %>%
  mutate(delta_DV = DV_Post - DV_Pre)

# --- Long-form for ANOVA (complete pairs only) ----
pheno_long <- meta %>%
  select(subject_key, Group, Timepoint, all_of(DV)) %>%
  filter(!is.na(.data[[DV]])) %>%
  filter(subject_key %in% pheno_wide$subject_key)

# --- RM-ANOVA ----
stats_anova <- rstatix::anova_test(
  data = pheno_long, dv = !!sym(DV),
  wid = subject_key, between = Group, within = Timepoint
)

# --- Paired t-tests within each group ----
grp_young <- pheno_wide %>% filter(Group == "Young")
grp_old   <- pheno_wide %>% filter(Group == "Old")
stats_paired_young <- t.test(grp_young$DV_Post, grp_young$DV_Pre, paired = TRUE)
stats_paired_old   <- t.test(grp_old$DV_Post,   grp_old$DV_Pre,   paired = TRUE)
stats_delta        <- t.test(delta_DV ~ Group, data = pheno_wide)

# --- ANOVA subtitle (plain text) ----
anova_tbl <- as.data.frame(stats_anova)
anova_sub <- sprintf("RM-ANOVA: Age %s   Time %s   Int. %s",
                     fmt_p(anova_tbl$p[anova_tbl$Effect == "Group"]),
                     fmt_p(anova_tbl$p[anova_tbl$Effect == "Timepoint"]),
                     fmt_p(anova_tbl$p[anova_tbl$Effect == "Group:Timepoint"]))
full_sub <- anova_sub

# --- Shapiro-Wilk normality on deltas ----
sw_dy <- shapiro.test(grp_young$delta_DV)
sw_do <- shapiro.test(grp_old$delta_DV)
n_y   <- nrow(grp_young)
n_o   <- nrow(grp_old)

# --- Audit CSV ----
audit_df <- data.frame(
  test      = c("paired_t_young", "paired_t_old", "unpaired_t_delta"),
  Group     = c("Young", "Old", "Young vs Old"),
  n         = c(n_y, n_o, nrow(pheno_wide)),
  statistic = c(stats_paired_young$statistic, stats_paired_old$statistic, stats_delta$statistic),
  p_value   = c(stats_paired_young$p.value,   stats_paired_old$p.value,   stats_delta$p.value),
  df        = c(stats_paired_young$parameter,  stats_paired_old$parameter,  stats_delta$parameter),
  mean_diff = c(stats_paired_young$estimate,   stats_paired_old$estimate,   diff(stats_delta$estimate)),
  ci_lo     = c(stats_paired_young$conf.int[1], stats_paired_old$conf.int[1], stats_delta$conf.int[1]),
  ci_hi     = c(stats_paired_young$conf.int[2], stats_paired_old$conf.int[2], stats_delta$conf.int[2]),
  shapiro_p = c(sw_dy$p.value, sw_do$p.value, NA)
)
write.csv(audit_df, file.path(cfg$dat, cfg$audit_file), row.names = FALSE)

# --- Plot data: optionally restrict to complete pairs ----
if (cfg$filter_complete) {
  plot_long <- meta %>%
    filter(subject_key %in% pheno_wide$subject_key, !is.na(.data[[DV]]))
} else {
  plot_long <- meta
}

y_max_left <- max(plot_long[[DV]], na.rm = TRUE)

# --- LEFT panel: pre/post grouped bar ----
p_left <- ggplot(plot_long, aes(x = Group_Time, y = .data[[DV]], fill = Group_Time)) +
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
  labs(title = cfg$title, subtitle = full_sub,
       y = cfg$y_label, x = NULL, tag = cfg$tag) +
  FIG_THEME +
  theme(plot.margin = cfg$left_margin, legend.position = "none")

# --- RIGHT panel: delta bar ----
delta_bar_colors <- c(Young = unname(GROUP_FILL["Young_Post"]),
                      Old   = unname(GROUP_FILL["Old_Post"]))

y_max_right <- max(pheno_wide$delta_DV, na.rm = TRUE)

p_right <- ggplot(pheno_wide, aes(x = Group, y = delta_DV, fill = Group)) +
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
  FIG_THEME +
  theme(legend.position = "none",
        axis.title.y = element_text(margin = margin(r = 1)),
        plot.margin = cfg$right_margin)

# Apply optional custom y-axis breaks on delta panel
if (!is.null(cfg$y_breaks)) {
  p_right <- p_right +
    scale_y_continuous(expand = expansion(mult = c(0.02, 0.22)),
                       breaks = cfg$y_breaks,
                       labels = cfg$y_labels)
} else {
  p_right <- p_right +
    scale_y_continuous(expand = expansion(mult = c(0.02, 0.22)))
}

# --- Combine & save standalone panel ----
p_combo <- (p_left | p_right) + plot_layout(widths = c(0.65, 0.35))

ggsave(file.path(cfg$rpt_png, paste0(cfg$file_prefix, "_", cfg$file_tag, ".png")),
       p_combo, width = PW, height = PH, units = "mm", dpi = 300)
ggsave(file.path(cfg$rpt_pdf, paste0(cfg$file_prefix, "_", cfg$file_tag, ".pdf")),
       p_combo, width = PW, height = PH, units = "mm", device = get_pdf_device())
cat(sprintf("F01 %s done\n", cfg$title))

# --- Export for composite (into calling environment) ----
pfx <- cfg$output_prefix

# Title (plain text)
assign(paste0(pfx, "_title"), cfg$title, envir = .GlobalEnv)

# Subtitle: plotmath expression string for main, plain text for supp
if (cfg$use_plotmath_subtitle) {
  assign(paste0(pfx, "_subtitle"),
         fmt_anova_sub(
           anova_tbl$p[anova_tbl$Effect == "Group"],
           anova_tbl$p[anova_tbl$Effect == "Timepoint"],
           anova_tbl$p[anova_tbl$Effect == "Group:Timepoint"]
         ),
         envir = .GlobalEnv)
} else {
  assign(paste0(pfx, "_subtitle"), full_sub, envir = .GlobalEnv)
}

assign(paste0(pfx, "_legend"), NULL, envir = .GlobalEnv)

# Stripped sub-panels for composite stitching
assign(paste0(pfx, "_left"),  strip_for_composite(p_left),  envir = .GlobalEnv)
assign(paste0(pfx, "_right"), strip_for_composite(p_right), envir = .GlobalEnv)
