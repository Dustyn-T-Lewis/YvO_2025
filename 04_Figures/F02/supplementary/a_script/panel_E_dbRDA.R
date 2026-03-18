# Figure 2 — Panel E (db-RDA): Constrained Ordination
# Distance-based RDA with Condition(subject) to partial out repeated measures.
# Outputs: pE_dbrda (ggplot object), panel_E_dbRDA.pdf/.png

setwd(rprojroot::find_rstudio_root_file())
source("04_Figures/F02/a_script/style.R")

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(readr)
  library(ggplot2)
  library(vegan)
})

PE_RDA_W <- 145; PE_RDA_H <- 100

RPT_DIR <- "04_Figures/F02/b_reports"
DAT_DIR <- "04_Figures/F02/c_data"
dir.create(RPT_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(DAT_DIR, recursive = TRUE, showWarnings = FALSE)

imp_df <- read_csv("02_Imputation/c_data/01_imputed.csv",
                   show_col_types = FALSE)

ann_cols   <- c("uniprot_id", "protein", "gene", "description")
samp_names <- setdiff(names(imp_df), ann_cols)

meta <- tibble(sample_id = samp_names) |>
  mutate(
    prefix   = str_extract(sample_id, "^[A-Z]+"),
    subj_num = str_extract(sample_id, "S\\d+"),
    time     = str_extract(sample_id, "(Pre|Post)$"),
    age      = if_else(str_detect(prefix, "^O"), "Old", "Young"),
    subject  = paste0(prefix, "_", subj_num),
    group    = paste(age, time, sep = "_")
  )
meta$age     <- factor(meta$age,  levels = c("Young", "Old"))
meta$time    <- factor(meta$time, levels = c("Pre", "Post"))
meta$group   <- factor(meta$group,
                       levels = c("Young_Pre", "Young_Post", "Old_Pre", "Old_Post"))
meta$subject <- factor(meta$subject)

imp_mat <- as.matrix(imp_df[, samp_names])
rownames(imp_mat) <- imp_df$gene

pdf_device <- get_pdf_device()

# db-RDA: constrained ordination
# Two models: (1) full model without conditioning for between-subject (age) variance
#             (2) conditional on subject for within-subject (time) variance
dist_mat <- vegdist(scale(t(imp_mat)), method = "euclidean")

# Full model: age + time + age:time (no conditioning)
rda_full <- dbrda(dist_mat ~ age * time, data = meta)

# Within-subject model: time + Condition(subject)
# Note: Condition(subject) absorbs age since subject is nested within age
rda_within <- dbrda(dist_mat ~ time + Condition(subject), data = meta)

# Use full model for visualization (shows both between and within effects)
site_scores  <- as.data.frame(scores(rda_full, display = "sites", choices = 1:2))
site_scores$sample_id <- rownames(site_scores)
site_scores <- left_join(site_scores, meta, by = "sample_id")

# Variance from full model
eigenvals_full <- eigenvals(rda_full)
constrained_eig <- eigenvals_full[grepl("^dbRDA", names(eigenvals_full))]
total_inertia   <- rda_full$tot.chi
constrained_var <- rda_full$CCA$tot.chi / total_inertia * 100

ax_var <- 100 * constrained_eig[1:2] / total_inertia

# ANOVA: full model terms and within-subject model for time
set.seed(42)
anova_full   <- anova(rda_full, by = "terms", permutations = 999)
anova_within <- anova(rda_within, by = "terms", permutations = 999)

# Extract from data frames by row name matching
get_anova_val <- function(av, term, col) {
  rn <- rownames(av)
  idx <- grep(paste0("^", term, "$"), rn)
  if (length(idx) == 0) return(NA_real_)
  av[idx, col]
}

# Column name depends on vegan version: "Variance" or "SumOfSqs"
var_col <- intersect(c("Variance", "SumOfSqs"), colnames(anova_full))[1]

total_ss_full   <- sum(anova_full[[var_col]], na.rm = TRUE)
total_ss_within <- sum(anova_within[[var_col]], na.rm = TRUE)

age_r2  <- get_anova_val(anova_full, "age", var_col) / total_ss_full
age_pv  <- get_anova_val(anova_full, "age", "Pr(>F)")
time_r2 <- get_anova_val(anova_within, "time", var_col) / total_ss_within
time_pv <- get_anova_val(anova_within, "time", "Pr(>F)")
int_r2  <- get_anova_val(anova_full, "age:time", var_col) / total_ss_full
int_pv  <- get_anova_val(anova_full, "age:time", "Pr(>F)")

adj_r2 <- RsquareAdj(rda_full)$adj.r.squared

# Safe p-value formatter
safe_fmt_p <- function(p) {
  if (is.na(p)) return("NA")
  fmt_p(p)
}

# Build annotation
anova_label <- sprintf(
  "db-RDA (adj. R\u00b2 = %.3f)\nAge  R\u00b2 = %.3f,  %s\nTime  R\u00b2 = %.3f,  %s  (cond.)\nAge\u00d7Time  R\u00b2 = %.3f,  %s",
  adj_r2,
  age_r2,  safe_fmt_p(age_pv),
  time_r2, safe_fmt_p(time_pv),
  int_r2,  safe_fmt_p(int_pv)
)

subtitle_txt <- sprintf(
  "Constrained R\u00b2 = %.1f%% | Time tested with Condition(subject) | %s proteins",
  constrained_var, format(nrow(imp_df), big.mark = ",")
)

pE_dbrda <- ggplot(site_scores, aes(x = dbRDA1, y = dbRDA2,
                                     color = group, shape = group)) +
  stat_ellipse(aes(fill = group), geom = "polygon",
               alpha = 0.10, level = 0.80, show.legend = FALSE) +
  stat_ellipse(aes(group = group), level = 0.80, linewidth = 0.4,
               linetype = "dashed", show.legend = FALSE) +
  geom_point(size = 2.5, alpha = 0.85) +
  annotate("text", x = Inf, y = Inf, label = anova_label,
           hjust = 1.05, vjust = 1.15,
           size = scale_text(BASE_STAT - 1.2, PE_RDA_W), color = "grey30",
           fontface = "bold") +
  scale_color_manual(values = PCA_COLORS,
                     labels = c("Young Pre", "Young Post", "Old Pre", "Old Post"),
                     guide = guide_legend(override.aes = list(size = 2.5))) +
  scale_fill_manual(values = PCA_COLORS, guide = "none") +
  scale_shape_manual(values = PCA_SHAPES,
                     labels = c("Young Pre", "Young Post", "Old Pre", "Old Post")) +
  labs(title = "Constrained Ordination (db-RDA)",
       subtitle = subtitle_txt,
       x = sprintf("dbRDA1 (%.1f%%)", ax_var[1]),
       y = sprintf("dbRDA2 (%.1f%%)", ax_var[2]),
       tag = "E'") +
  FIG_THEME + theme(legend.position = c(0.88, 0.12),
                    legend.background = element_rect(fill = alpha("white", 0.8),
                                                     color = NA),
                    legend.title = element_blank(),
                    legend.text = element_text(size = FIG_LEGEND_TEXT),
                    legend.key.size = unit(3, "mm"),
                    legend.spacing.y = unit(0.5, "mm"))

audit_df <- data.frame(
  term         = c("age", "time_conditioned", "age:time"),
  R2           = as.numeric(c(age_r2, time_r2, int_r2)),
  p_value      = as.numeric(c(age_pv, time_pv, int_pv)),
  constrained_R2_pct = rep(constrained_var, 3),
  adj_R2       = rep(adj_r2, 3),
  dbRDA1_pct   = rep(as.numeric(ax_var[1]), 3),
  dbRDA2_pct   = rep(as.numeric(ax_var[2]), 3)
)
write.csv(audit_df, file.path(DAT_DIR, "audit_panel_E_dbRDA.csv"),
          row.names = FALSE)

ggsave(file.path(RPT_DIR, "panel_E_dbRDA.pdf"), pE_dbrda,
       width = PE_RDA_W, height = PE_RDA_H, units = "mm", device = pdf_device)
ggsave(file.path(RPT_DIR, "panel_E_dbRDA.png"), pE_dbrda,
       width = PE_RDA_W, height = PE_RDA_H, units = "mm", dpi = 300)
