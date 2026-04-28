#!/usr/bin/env Rscript
# F07 Supplementary — Per-module univariate ROC grid (standalone SUPP)
#
# 6 outcome rows x 8 modules = 48 univariate ROCs. Each cell asks:
#   "How well does this single module eigengene classify this outcome?"
#
# Sourced by 01_main_panels.R — expects style.R + figure_supplement_helpers.R
# already loaded.
# Exports: module_grid_summary.csv + module_grid_curves.csv -> c_data/module_grid/
#          SUPP_F07_module_grid.{png,pdf} -> b_reports/supp/{png,pdf}/panels/

suppressPackageStartupMessages({
  library(tidyverse); library(patchwork); library(pROC)
})

BASE    <- here::here("04_Figures", "F07")
DAT_OUT <- file.path(BASE, "c_data", "module_grid")
RPT_PNG <- file.path(BASE, "b_reports", "supp", "png", "panels")
RPT_PDF <- file.path(BASE, "b_reports", "supp", "pdf", "panels")
dir.create(DAT_OUT, recursive = TRUE, showWarnings = FALSE)
dir.create(RPT_PNG, recursive = TRUE, showWarnings = FALSE)
dir.create(RPT_PDF, recursive = TRUE, showWarnings = FALSE)

F06_SUPP <- here::here("04_Figures", "F06", "c_data", "F06_supplementary.xlsx")
stopifnot("F06 stitcher must run first: missing F06_supplementary.xlsx" =
  file.exists(F06_SUPP))

# ---- Data --------------------------------------------------------------------
MEs     <- read_matrix_sheet(F06_SUPP, "MEs",     "sample_id")
me_pre  <- read_matrix_sheet(F06_SUPP, "me_pre",  "subject_key")
me_post <- read_matrix_sheet(F06_SUPP, "me_post", "subject_key")
subj_age<- read_sheet_df(F06_SUPP, "metadata_subj_age")
pheno   <- read_sheet_df(F06_SUPP, "metadata_pheno_wide")
common_subj <- read_vector_sheet(F06_SUPP, "common_subj")

MODULES <- c("turquoise","blue","brown","yellow","green","red","black","pink")
ME_cols <- paste0("ME", MODULES)

# ---- Outcome vectors ---------------------------------------------------------
age_bin <- ifelse(subj_age$age[match(common_subj, subj_age$subject_key)] == "Old", 1, 0)

# Combined ME = (Pre + Post)/2 per subject
me_comb <- (as.matrix(me_pre[common_subj, ME_cols]) +
            as.matrix(me_post[common_subj, ME_cols])) / 2

# delta-VL responder: age-residualize, median split
pheno_s <- pheno[match(common_subj, pheno$subject_key), ]
ok_vl <- !is.na(pheno_s$delta_VL)
resid_vl <- rep(NA_real_, length(common_subj))
resid_vl[ok_vl] <- residuals(lm(pheno_s$delta_VL[ok_vl] ~ age_bin[ok_vl]))
vl_bin <- ifelse(resid_vl > median(resid_vl, na.rm = TRUE), 1, 0)

# Pre vs Post pooled (paired): stack Pre then Post for common_subj
me_pre_s  <- as.matrix(me_pre[common_subj, ME_cols])
me_post_s <- as.matrix(me_post[common_subj, ME_cols])
tp_bin <- c(rep(0, length(common_subj)), rep(1, length(common_subj)))  # 0=Pre, 1=Post
subj_id <- c(common_subj, common_subj)

# ---- AUC + permutation p -----------------------------------------------------
uni_auc <- function(y, x) {
  ok <- !is.na(y) & !is.na(x)
  if (length(unique(y[ok])) < 2) return(list(auc = NA, ci_lo = NA, ci_hi = NA, roc = NULL))
  r <- suppressMessages(roc(y[ok], x[ok], quiet = TRUE, direction = "auto"))
  # DeLong CI assumes independence — anti-conservative for paired Pre/Post.
  # Permutation p-value (below) is the primary inferential quantity.
  ci <- as.numeric(ci.auc(r))
  list(auc = as.numeric(auc(r)), ci_lo = ci[1], ci_hi = ci[3], roc = r)
}

perm_p_unpaired <- function(y, x, obs_auc, n_perm = 1000) {
  ok <- !is.na(y) & !is.na(x)
  y <- y[ok]; x <- x[ok]
  nulls <- replicate(n_perm, {
    ys <- sample(y)
    r <- suppressMessages(roc(ys, x, quiet = TRUE, direction = "auto"))
    max(as.numeric(auc(r)), 1 - as.numeric(auc(r)))
  })
  (sum(nulls >= max(obs_auc, 1 - obs_auc)) + 1) / (n_perm + 1)
}

perm_p_paired <- function(tp, subj, x, obs_auc, n_perm = 1000) {
  nulls <- replicate(n_perm, {
    tp_shuf <- tp
    for (s in unique(subj)) {
      idx <- which(subj == s)
      tp_shuf[idx] <- sample(tp_shuf[idx])
    }
    r <- suppressMessages(roc(tp_shuf, x, quiet = TRUE, direction = "auto"))
    max(as.numeric(auc(r)), 1 - as.numeric(auc(r)))
  })
  (sum(nulls >= max(obs_auc, 1 - obs_auc)) + 1) / (n_perm + 1)
}

set.seed(42)
rows <- list()
curves <- list()

message("Row 1: Age ~ Combined ME ...")
for (j in seq_along(MODULES)) {
  mod <- MODULES[j]; x <- me_comb[, ME_cols[j]]
  r <- uni_auc(age_bin, x)
  p <- perm_p_unpaired(age_bin, x, r$auc)
  rows[[length(rows)+1]] <- tibble(row="Age", module=mod, n=length(age_bin),
    n_pos=sum(age_bin==1), n_neg=sum(age_bin==0),
    auc=r$auc, ci_lo=r$ci_lo, ci_hi=r$ci_hi, perm_p=p)
  curves[[length(curves)+1]] <- tibble(row="Age", module=mod,
    fpr=1 - r$roc$specificities, tpr=r$roc$sensitivities)
}

message("Row 2: delta-VL responder ~ Pre ME ...")
for (j in seq_along(MODULES)) {
  mod <- MODULES[j]; x <- me_pre_s[, ME_cols[j]]
  ok <- !is.na(vl_bin)
  r <- uni_auc(vl_bin[ok], x[ok])
  p <- if (!is.na(r$auc)) perm_p_unpaired(vl_bin[ok], x[ok], r$auc) else NA
  rows[[length(rows)+1]] <- tibble(row="\u0394VL-responder", module=mod, n=sum(ok),
    n_pos=sum(vl_bin[ok]==1), n_neg=sum(vl_bin[ok]==0),
    auc=r$auc, ci_lo=r$ci_lo, ci_hi=r$ci_hi, perm_p=p)
  curves[[length(curves)+1]] <- tibble(row="\u0394VL-responder", module=mod,
    fpr=1 - r$roc$specificities, tpr=r$roc$sensitivities)
}

message("Row 3: Pre vs Post (pooled, paired) ~ raw ME ...")
me_pooled <- rbind(me_pre_s, me_post_s)
for (j in seq_along(MODULES)) {
  mod <- MODULES[j]; x <- me_pooled[, ME_cols[j]]
  r <- uni_auc(tp_bin, x)
  p <- perm_p_paired(tp_bin, subj_id, x, r$auc)
  rows[[length(rows)+1]] <- tibble(row="Pre vs Post", module=mod, n=length(tp_bin),
    n_pos=sum(tp_bin==1), n_neg=sum(tp_bin==0),
    auc=r$auc, ci_lo=r$ci_lo, ci_hi=r$ci_hi, perm_p=p)
  curves[[length(curves)+1]] <- tibble(row="Pre vs Post", module=mod,
    fpr=1 - r$roc$specificities, tpr=r$roc$sensitivities)
}

message("Row 4: Pre vs Post (YOUNG only) ~ raw ME ...")
young_subj <- common_subj[age_bin == 0]
y_pre  <- as.matrix(me_pre[young_subj,  ME_cols])
y_post <- as.matrix(me_post[young_subj, ME_cols])
y_pool <- rbind(y_pre, y_post)
y_tp   <- c(rep(0, length(young_subj)), rep(1, length(young_subj)))
y_sid  <- c(young_subj, young_subj)
for (j in seq_along(MODULES)) {
  mod <- MODULES[j]; x <- y_pool[, ME_cols[j]]
  r <- uni_auc(y_tp, x)
  p <- perm_p_paired(y_tp, y_sid, x, r$auc)
  rows[[length(rows)+1]] <- tibble(row="Pre vs Post (Young)", module=mod,
    n=length(y_tp), n_pos=sum(y_tp==1), n_neg=sum(y_tp==0),
    auc=r$auc, ci_lo=r$ci_lo, ci_hi=r$ci_hi, perm_p=p)
  curves[[length(curves)+1]] <- tibble(row="Pre vs Post (Young)", module=mod,
    fpr=1 - r$roc$specificities, tpr=r$roc$sensitivities)
}

message("Row 5: Pre vs Post (OLD only) ~ raw ME ...")
old_subj <- common_subj[age_bin == 1]
o_pre  <- as.matrix(me_pre[old_subj,  ME_cols])
o_post <- as.matrix(me_post[old_subj, ME_cols])
o_pool <- rbind(o_pre, o_post)
o_tp   <- c(rep(0, length(old_subj)), rep(1, length(old_subj)))
o_sid  <- c(old_subj, old_subj)
for (j in seq_along(MODULES)) {
  mod <- MODULES[j]; x <- o_pool[, ME_cols[j]]
  r <- uni_auc(o_tp, x)
  p <- perm_p_paired(o_tp, o_sid, x, r$auc)
  rows[[length(rows)+1]] <- tibble(row="Pre vs Post (Old)", module=mod,
    n=length(o_tp), n_pos=sum(o_tp==1), n_neg=sum(o_tp==0),
    auc=r$auc, ci_lo=r$ci_lo, ci_hi=r$ci_hi, perm_p=p)
  curves[[length(curves)+1]] <- tibble(row="Pre vs Post (Old)", module=mod,
    fpr=1 - r$roc$specificities, tpr=r$roc$sensitivities)
}

message("Row 6: High vs Low baseline LBM (age-residualized) ~ Pre ME ...")
ok_lbm <- !is.na(pheno_s$LBM_Pre)
resid_lbm <- rep(NA_real_, length(common_subj))
resid_lbm[ok_lbm] <- residuals(lm(pheno_s$LBM_Pre[ok_lbm] ~ age_bin[ok_lbm]))
lbm_bin <- ifelse(resid_lbm > median(resid_lbm, na.rm = TRUE), 1L, 0L)
for (j in seq_along(MODULES)) {
  mod <- MODULES[j]; x <- me_pre_s[, ME_cols[j]]
  ok <- !is.na(lbm_bin)
  r <- uni_auc(lbm_bin[ok], x[ok])
  p <- if (!is.na(r$auc)) perm_p_unpaired(lbm_bin[ok], x[ok], r$auc) else NA
  rows[[length(rows)+1]] <- tibble(row="LBM (High vs Low)", module=mod, n=sum(ok),
    n_pos=sum(lbm_bin[ok]==1), n_neg=sum(lbm_bin[ok]==0),
    auc=r$auc, ci_lo=r$ci_lo, ci_hi=r$ci_hi, perm_p=p)
  curves[[length(curves)+1]] <- tibble(row="LBM (High vs Low)", module=mod,
    fpr=1 - r$roc$specificities, tpr=r$roc$sensitivities)
}

summ <- bind_rows(rows) |>
  mutate(q_bh = p.adjust(perm_p, method = "BH"),
         sig  = case_when(
           is.na(q_bh)        ~ "ns",
           q_bh   < 0.05      ~ "q<.05",
           perm_p < 0.05      ~ "p<.05",
           TRUE               ~ "ns"),
         border_color = ifelse(sig == "ns", "grey80", "black"),
         border_lty   = ifelse(sig == "p<.05", "dashed", "solid"),
         border_lw    = case_when(sig == "q<.05" ~ 2.4,
                                  sig == "p<.05" ~ 2.4,
                                  TRUE           ~ 0.3))
curves_df <- bind_rows(curves)

write_csv(summ,      file.path(DAT_OUT, "module_grid_summary.csv"))
write_csv(curves_df, file.path(DAT_OUT, "module_grid_curves.csv"))

message("AUC ranges:")
print(summ |> group_by(row) |>
  summarise(min=min(auc,na.rm=TRUE), max=max(auc,na.rm=TRUE),
            n_qsig=sum(sig=="q<.05"), n_psig=sum(sig=="p<.05")))

# ---- Cell builder ------------------------------------------------------------
build_cell <- function(row_name, module) {
  rr <- summ |> filter(row == row_name, module == !!module)
  dd <- curves_df |> filter(row == row_name, module == !!module)
  col <- module
  if (is.na(rr$auc)) {
    auc_line <- "NA"; ci_line <- ""; p_line <- ""
  } else {
    auc_line <- sprintf("AUC %.2f", rr$auc)
    ci_line  <- sprintf("[%.2f, %.2f]", rr$ci_lo, rr$ci_hi)
    p_line   <- if (rr$perm_p < 0.001) "p < .001" else sprintf("p = %.3f", rr$perm_p)
  }
  txt_color <- "grey10"
  ggplot(dd, aes(fpr, tpr)) +
    geom_ribbon(aes(ymin = 0, ymax = tpr), fill = col, alpha = 0.32) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed",
                color = "grey70", linewidth = 0.3) +
    geom_line(color = col, linewidth = 1.15) +
    annotate("text", x = 0.97, y = 0.44, label = auc_line,
             hjust = 1, vjust = 0.5, size = 4.6, fontface = "bold",
             color = txt_color) +
    annotate("text", x = 0.97, y = 0.26, label = ci_line,
             hjust = 1, vjust = 0.5, size = 3.6, color = txt_color) +
    annotate("text", x = 0.97, y = 0.09, label = p_line,
             hjust = 1, vjust = 0.5, size = 4.0, fontface = "bold",
             color = txt_color) +
    scale_x_continuous(limits = c(0,1), breaks = c(0,1), expand = c(0,0)) +
    scale_y_continuous(limits = c(0,1), breaks = c(0,1), expand = c(0,0)) +
    coord_fixed() +
    theme_void() +
    theme(panel.border = element_rect(fill = NA,
                                      color = rr$border_color,
                                      linewidth = rr$border_lw,
                                      linetype = rr$border_lty),
          plot.margin = margin(2, 2, 2, 2))
}

# ---- Headers (transposed: outcomes = columns, modules = rows) ----------------
outcome_header <- function(top, bot="") {
  p <- ggplot() + theme_void() +
    annotate("text", x=0.5, y=0.70, label=top,
             fontface="bold", size=3.5, color="grey10") +
    coord_cartesian(xlim=c(0,1), ylim=c(0,1), expand=FALSE)
  if (nzchar(bot)) {
    p <- p + annotate("text", x=0.5, y=0.26, label=bot,
                      fontface="italic", size=2.8, color="grey40")
  }
  p + theme(plot.margin = margin(1, 3, 1, 3))
}

module_row_label <- function(mod) {
  ggplot() + theme_void() +
    annotate("rect", xmin=0, xmax=0.34, ymin=0.08, ymax=0.92, fill=mod) +
    annotate("text", x=0.68, y=0.5, label=stringr::str_to_title(mod),
             size=3.2, fontface="bold", color="grey10", hjust=0.5) +
    coord_cartesian(xlim=c(0,1), ylim=c(0,1), expand=FALSE) +
    theme(plot.margin = margin(2, 2, 2, 2))
}

# Outcome column definitions: internal row_name + header top/bot
OUTCOMES <- list(
  list(key="Age",                 top="Age (Y vs O)",      bot="Combined ME"),
  list(key="\u0394VL-responder",  top="\u0394 VL responder", bot="Pre ME"),
  list(key="Pre vs Post",         top="Pre vs Post (all)", bot="Raw ME"),
  list(key="Pre vs Post (Young)", top="Pre vs Post (Y)",   bot="Raw ME"),
  list(key="Pre vs Post (Old)",   top="Pre vs Post (O)",   bot="Raw ME"),
  list(key="LBM (High vs Low)",   top="Baseline LBM",      bot="Pre ME")
)

# ---- Assemble (transposed) ---------------------------------------------------
header_cells <- c(list(patchwork::plot_spacer()),
                  lapply(OUTCOMES, function(o) outcome_header(o$top, o$bot)))

row_cells <- function(mod) {
  c(list(module_row_label(mod)),
    lapply(OUTCOMES, function(o) build_cell(o$key, mod)))
}

module_rows <- lapply(MODULES, row_cells)
all_cells <- c(header_cells, unlist(module_rows, recursive = FALSE))

composite <- wrap_plots(all_cells, nrow = 9, ncol = 7,
                        widths  = c(1.05, rep(1, 6)),
                        heights = c(0.30, rep(1, 8)))

composite <- composite +
  plot_annotation(
    title = "F07 Supplementary \u2014 Per-module ROC grid: signal decomposition across outcomes",
    subtitle = paste0(
      "Each cell = univariate AUC from one module eigengene.\n",
      "Border: solid black = BH q<0.05; dashed = nominal p<0.05; grey = ns.  1000 permutations (Phipson\u2013Smyth adjusted); paired within subject for Pre vs Post rows."),
    caption = paste0(
      "Row 1 (Age): Combined ME = (Pre+Post)/2 per subject, n=30.  ",
      "Row 2 (\u0394VL responder): median-split of age-residualized \u0394VL, predictor = Pre ME.  ",
      "Rows 4-5 (age-specific Pre vs Post): paired within subject; Young n=15 \u00d7 2, Old n=15 \u00d7 2.  ",
      "Row 6 (LBM): median-split of age-residualized baseline LBM, predictor = Pre ME.  ",
      "Row 3 (Pre vs Post): raw ME values stacked across 30 subjects \u00d7 2 timepoints = 60 obs, paired permutation.  ",
      "AUC shown centered inside each ROC. Data: 04_Figures/F07/c_data/module_grid/."),
    theme = theme(
      plot.title    = element_text(face = "bold", size = FIG_TITLE_SIZE + 4, color = "grey10",
                                   margin = margin(t = 14, l = 52, b = 4, unit = "pt")),
      plot.subtitle = element_text(size = FIG_SUBTITLE_SIZE + 4, face = "italic", color = "grey30",
                                   margin = margin(l = 52, b = 12, unit = "pt"),
                                   lineheight = 1.25),
      plot.caption  = element_text(size = FIG_AXIS_TEXT, color = "grey45", hjust = 0,
                                   lineheight = 1.3)))

W_in <- 12; H_in <- 16
ggsave(file.path(RPT_PDF, "SUPP_F07_module_grid.pdf"), composite,
       width = W_in, height = H_in, units = "in", device = get_pdf_device())
ggsave(file.path(RPT_PNG, "SUPP_F07_module_grid.png"), composite,
       width = W_in, height = H_in, units = "in", dpi = 300)

message("Wrote: SUPP_F07_module_grid.{pdf,png}")
