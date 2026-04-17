# F02 Supp Panel D: CV% Violins (Inter-Individual Variability)
# Faceted by Age (Young | Old), Pre/Post on x-axis. Median labels w/ bootstrap CIs.
# Outputs: pD (ggplot object), panel_D_cv_SUPP.{pdf,png}

setwd(rprojroot::find_rstudio_root_file())
source("04_Figures/F02/a_script/style.R")

library(dplyr)
library(tidyr)
library(stringr)
library(readr)
library(ggplot2)
library(ggbeeswarm)

PB_W <- 110
PB_H <- 120

RPT     <- "04_Figures/F02/b_reports/supp/panels"
DAT_DIR <- "04_Figures/F02/c_data"
dir.create(RPT,     recursive = TRUE, showWarnings = FALSE)
dir.create(DAT_DIR, recursive = TRUE, showWarnings = FALSE)

norm_df <- read_csv("01_normalization/c_data/02_normalized.csv",
                    show_col_types = FALSE)

ann_cols   <- c("uniprot_id", "protein", "gene", "description")
samp_names <- setdiff(names(norm_df), ann_cols)

meta <- tibble(sample_id = samp_names) |>
  mutate(
    prefix   = str_extract(sample_id, "^[A-Z]+"),
    subj_num = str_extract(sample_id, "S\\d+"),
    time     = str_extract(sample_id, "(Pre|Post)$"),
    age      = if_else(str_detect(prefix, "^O"), "Old", "Young"),
    subject  = paste0(prefix, "_", subj_num),
    group    = paste(age, time, sep = "_")
  )
meta$age   <- factor(meta$age,  levels = c("Young", "Old"))
meta$time  <- factor(meta$time, levels = c("Pre", "Post"))
meta$group <- factor(meta$group,
                     levels = c("Young_Pre", "Young_Post", "Old_Pre", "Old_Post"))

pdf_device <- get_pdf_device()

# CV on linear scale (Brenes 2024)
lin_mat <- 2^as.matrix(norm_df[, samp_names])

cv_list <- lapply(levels(meta$group), function(g) {
  idx <- meta$sample_id[meta$group == g]
  sub <- lin_mat[, idx, drop = FALSE]
  cv_pct <- apply(sub, 1, function(x) {
    x <- x[!is.na(x)]
    if (length(x) < 2) return(NA_real_)
    sd(x) / mean(x) * 100
  })
  tibble(group = g, cv = cv_pct)
})
cv_df <- bind_rows(cv_list) |> filter(!is.na(cv))
cv_df$group <- factor(cv_df$group,
                      levels = c("Young_Pre", "Young_Post", "Old_Pre", "Old_Post"))

cv_df$age  <- factor(ifelse(grepl("Young", cv_df$group), "Young", "Old"),
                     levels = c("Young", "Old"))
cv_df$time <- factor(ifelse(grepl("Pre", cv_df$group), "Pre", "Post"),
                     levels = c("Pre", "Post"))

# Bootstrap 95% CI on median CV per group
set.seed(42)
boot_median_ci <- function(x, R = 2000, conf = 0.95) {
  meds <- replicate(R, median(sample(x, replace = TRUE)))
  qs   <- quantile(meds, c((1 - conf) / 2, (1 + conf) / 2))
  c(lower = unname(qs[1]), upper = unname(qs[2]))
}

# Pairwise Wilcoxon tests, BH corrected (audit only — not on figure)
bracket_comps <- list(c("Young_Pre", "Young_Post"), c("Old_Pre", "Old_Post"),
                      c("Young_Pre", "Old_Pre"),    c("Young_Post", "Old_Post"))
bracket_pvals_raw <- sapply(bracket_comps, function(pair)
  wilcox.test(cv_df$cv[cv_df$group == pair[1]],
              cv_df$cv[cv_df$group == pair[2]])$p.value)
bracket_pvals <- p.adjust(bracket_pvals_raw, method = "BH")

cliffs_delta <- function(x, y) {
  nx <- length(x); ny <- length(y)
  d  <- outer(x, y, function(a, b) sign(a - b))
  sum(d) / (nx * ny)
}
cliff_results <- data.frame(
  comparison = sapply(bracket_comps, paste, collapse = " vs "),
  p_raw      = bracket_pvals_raw,
  p_bh       = bracket_pvals,
  cliffs_d   = sapply(bracket_comps, function(pair)
    cliffs_delta(cv_df$cv[cv_df$group == pair[1]],
                 cv_df$cv[cv_df$group == pair[2]]))
)

cv_ci <- cv_df |>
  group_by(age, time, group) |>
  summarise(
    med    = median(cv),
    ci_lo  = boot_median_ci(cv)[["lower"]],
    ci_hi  = boot_median_ci(cv)[["upper"]],
    cv_max = max(cv),
    .groups = "drop"
  )

n_prot    <- nrow(norm_df)
grand_med <- median(cv_df$cv)
grand_ci  <- boot_median_ci(cv_df$cv)

# Delta median CV per age group (Post - Pre)
delta_cv <- cv_ci |>
  select(age, time, med) |>
  pivot_wider(names_from = time, values_from = med) |>
  mutate(delta = Post - Pre,
         arrow_label = sprintf("%+.1f%%", delta))

# Arrow annotation data: one arrow per age facet from Pre median to Post median
arrow_df <- delta_cv |>
  mutate(x = 1, xend = 2,
         y_mid = (Pre + Post) / 2)

sub_txt <- sprintf(
  paste0("Group-level CV%% per Age x Time | %s proteins\n",
         "median %.0f%% [%.0f-%.0f boot CI] | ",
         "Pre->Post: Y %+.1f%%, O %+.1f%%"),
  format(n_prot, big.mark = ","), grand_med, grand_ci[1], grand_ci[2],
  delta_cv$delta[delta_cv$age == "Young"],
  delta_cv$delta[delta_cv$age == "Old"]
)

pD <- ggplot(cv_df, aes(x = time, y = cv, fill = group)) +
  geom_violin(alpha = 0.5, linewidth = 0.3, color = "black", scale = "width") +
  geom_quasirandom(aes(color = group), alpha = 0.15, size = 0.5,
                   width = 0.25, groupOnX = TRUE, show.legend = FALSE) +
  geom_boxplot(width = 0.15, outlier.shape = NA, linewidth = 0.3,
               color = "black", fill = "white", coef = 0) +
  geom_hline(yintercept = 25, linetype = "dashed", color = "grey50",
             linewidth = 0.4) +
  geom_label(data = cv_ci,
             aes(x = time, y = cv_max + 3,
                 label = sprintf("%.0f%% [%.0f\u2013%.0f]", med, ci_lo, ci_hi)),
             size = scale_text(BASE_COUNT + 0.5, PB_W / 2),
             fontface = "bold", fill = alpha("white", 0.8),
             linewidth = 0.2, label.padding = unit(1.5, "pt"),
             hjust = 0, nudge_x = -0.4) +
  geom_segment(data = arrow_df,
               aes(x = x, xend = xend, y = Pre, yend = Post),
               inherit.aes = FALSE, color = "grey30",
               arrow = arrow(length = unit(1.5, "mm"), type = "closed"),
               linewidth = 0.6) +
  geom_label(data = arrow_df,
             aes(x = 1.5, y = y_mid, label = arrow_label),
             inherit.aes = FALSE, size = scale_text(BASE_COUNT, PB_W / 2),
             fontface = "bold.italic", fill = alpha("white", 0.85),
             label.padding = unit(1.5, "pt"), linewidth = 0.2,
             color = "grey30") +
  facet_wrap(~ age, nrow = 1) +
  scale_fill_manual(values = GROUP_FILL) +
  scale_color_manual(values = GROUP_FILL) +
  coord_cartesian(ylim = c(0, max(cv_ci$cv_max) + 12)) +
  labs(title = "Inter-Individual Variability (CV%)",
       subtitle = sub_txt,
       x = NULL, y = "CV (%)",
       tag = "D") +
  FIG_THEME +
  theme(legend.position = "none",
        panel.spacing   = unit(8, "mm"),
        plot.title      = element_text(margin = margin(b = 0)),
        plot.subtitle   = element_text(size = FIG_SUBTITLE_SIZE - 1.0,
                                       face = "bold.italic", color = "grey40",
                                       margin = margin(t = 0, b = 1)),
        strip.text      = element_text(face = "bold", size = FIG_STRIP_SIZE,
                                       margin = margin(t = 0, b = 1)),
        plot.margin     = margin(t = 0, r = 5.5, b = 5.5, l = 5.5))

write.csv(as.data.frame(cv_ci),
          file.path(DAT_DIR, "audit_panel_D_median_cv_ci.csv"), row.names = FALSE)
write.csv(cliff_results,
          file.path(DAT_DIR, "audit_panel_D_wilcoxon_effects.csv"), row.names = FALSE)

ggsave(file.path(RPT, "SUPP_panel_D_cv.png"), pD,
       width = PB_W, height = PB_H, units = "mm", dpi = 300)

cat("F02 Supp Panel D done\n")
