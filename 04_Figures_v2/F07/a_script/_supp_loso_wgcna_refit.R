#!/usr/bin/env Rscript
# F07 Supplementary — Full WGCNA-refit LOSO sensitivity
#
# Strongest sensitivity for module-definition circularity. For each of
# 30 subjects, remove their 2 samples from datExpr, refit the WGCNA network
# on the remaining 58 samples (same parameters, soft power fixed at 12),
# match training modules to full-sample modules by Jaccard overlap,
# project held-out subject's MEs onto training PCs, compute held-out AUC.
# Outcome binarization (median split / age-residual) is fit on the training
# subjects only and applied to the held-out subject, so no label leaks in.
#
# Runtime: ~13-15 min (30 folds x ~25 s/fold blockwiseModules).
#
# Sourced by 01_main_panels.R — expects style.R + figure_supplement_helpers.R
# already loaded.

suppressPackageStartupMessages({
  library(tidyverse)
  library(WGCNA)
  library(pROC)
})

source("04_Figures_v2/F07/a_script/_f07_helpers.R")

disableWGCNAThreads() # single-threaded so correlations/TOM are byte-reproducible across re-runs
set.seed(42)
# blockwiseModules dispatches on cor(); point it at WGCNA's version for the
# refit loop, then restore base cor() below (top-level on.exit never fires).
cor <- WGCNA::cor

BASE <- "04_Figures_v2/F07"
DAT_OUT <- file.path(BASE, "c_data", "loso_auc")
RPT_PNG <- file.path(BASE, "b_reports", "supp", "png", "panels")
RPT_PDF <- file.path(BASE, "b_reports", "supp", "pdf", "panels")
dir.create(DAT_OUT, recursive = TRUE, showWarnings = FALSE)
dir.create(RPT_PNG, recursive = TRUE, showWarnings = FALSE)
dir.create(RPT_PDF, recursive = TRUE, showWarnings = FALSE)

# Inputs
F06_SUPP <- "04_Figures_v2/F06/c_data/F06_supplementary.xlsx"
stopifnot(
  "F06 must run first: missing F06_supplementary.xlsx" =
    file.exists(F06_SUPP)
)
datExpr <- readRDS("04_Figures_v2/F06/c_data/datExpr.rds")
mod_full <- readRDS("04_Figures_v2/F06/c_data/module_colors.rds")
soft_power <- read_csv("04_Figures_v2/F06/c_data/wgcna/wgcna_sft_summary.csv",
  show_col_types = FALSE
)$selected_power[1]
me_pre <- readRDS("04_Figures_v2/F06/c_data/me_pre.rds")
me_post <- readRDS("04_Figures_v2/F06/c_data/me_post.rds")
subj_age <- read_sheet_df(F06_SUPP, "metadata_subj_age")
pheno <- read_sheet_df(F06_SUPP, "metadata_pheno_wide")

in_sample_csv <- file.path(BASE, "c_data", "module_grid", "module_grid_summary.csv")
in_sample_xlsx <- file.path(BASE, "c_data", "F07_supplementary.xlsx")
in_sample <- if (file.exists(in_sample_csv)) {
  read_csv(in_sample_csv, show_col_types = FALSE)
} else {
  read_sheet_df(in_sample_xlsx, "module_grid_summary")
}

stopifnot(length(mod_full) == ncol(datExpr))
names(mod_full) <- colnames(datExpr)

common_subj <- intersect(rownames(me_pre), rownames(me_post))
MODULES <- c("turquoise", "blue", "brown", "yellow", "green", "red", "black", "pink")

# Top-12 pairs
top12 <- in_sample |>
  filter(!is.na(perm_p), !is.na(auc), module %in% MODULES) |>
  arrange(perm_p, desc(auc)) |>
  slice_head(n = 12)

# Outcome vectors. Age is a fixed label; the responder/LBM splits are fit on
# the n-1 training subjects per fold and applied to the held-out subject, so
# the binarization threshold never sees the subject it is scoring.
age_bin <- f07_age_bin(subj_age, common_subj)

pheno_s <- pheno[match(common_subj, pheno$subject_key), ]
vl_bin <- f07_loso_resid_median_split(pheno_s$delta_VL, age_bin)
lbm_bin <- f07_loso_resid_median_split(pheno_s$LBM_Pre, age_bin)

# Helpers
match_modules <- function(mod_train, mod_full) {
  train_lvls <- setdiff(unique(mod_train), "grey")
  full_lvls <- setdiff(unique(mod_full), "grey")
  jaccard <- matrix(0,
    nrow = length(full_lvls), ncol = length(train_lvls),
    dimnames = list(full_lvls, train_lvls)
  )
  for (f in full_lvls) {
    for (t in train_lvls) {
      a <- names(mod_full)[mod_full == f]
      b <- names(mod_train)[mod_train == t]
      jaccard[f, t] <- length(intersect(a, b)) / length(union(a, b))
    }
  }
  best <- apply(jaccard, 1, function(x) colnames(jaccard)[which.max(x)])
  best_j <- apply(jaccard, 1, max)
  tibble(full = rownames(jaccard), train = best, jaccard = best_j)
}

compute_me_projection <- function(datExpr_full, train_rows, mod_train_full_len,
                                  module_name) {
  prot_ids <- names(mod_train_full_len)[mod_train_full_len == module_name]
  if (length(prot_ids) < 2) {
    return(NULL)
  }
  fit <- fit_pc1(datExpr_full[train_rows, prot_ids, drop = FALSE])
  function(rows) project_pc1(fit, datExpr_full[rows, prot_ids, drop = FALSE])
}

refit_fold <- function(hold_subj, datExpr_full, mod_full_named, soft_power = 12L) {
  hold_ids <- paste0(hold_subj, c("_Pre", "_Post"))
  train_ids <- setdiff(rownames(datExpr_full), hold_ids)
  train_idx <- match(train_ids, rownames(datExpr_full))
  hold_idx <- match(hold_ids, rownames(datExpr_full))
  datExpr_train <- datExpr_full[train_idx, , drop = FALSE]

  net <- blockwiseModules(
    datExpr_train,
    power             = soft_power,
    networkType       = "signed",
    TOMType           = "signed",
    minModuleSize     = 30,
    mergeCutHeight    = 0.25,
    numericLabels     = TRUE,
    pamRespectsDendro = FALSE,
    saveTOMs          = FALSE,
    verbose           = 0
  )
  mod_train <- labels2colors(net$colors)
  names(mod_train) <- colnames(datExpr_train)
  mm <- match_modules(mod_train, mod_full_named)

  projs <- map(MODULES, function(mod_full_name) {
    row <- mm |> filter(full == mod_full_name)
    if (nrow(row) == 0 || row$jaccard < 0.5) {
      return(NULL)
    }
    fn <- compute_me_projection(datExpr_full, train_idx, mod_train, row$train)
    if (is.null(fn)) {
      return(NULL)
    }
    list(
      module = mod_full_name, match = row$train, jaccard = row$jaccard,
      me_hold = fn(hold_idx), me_hold_names = rownames(datExpr_full)[hold_idx]
    )
  }) |> set_names(MODULES)

  list(projs = projs, jaccard_summary = mm)
}

# Run LOSO over all subjects
n_subj <- length(common_subj)
per_fold <- vector("list", n_subj)
stability <- vector("list", n_subj)

t0 <- Sys.time()
for (i in seq_len(n_subj)) {
  s <- common_subj[i]
  message(sprintf(
    "[fold %d/%d] subject %s  (elapsed: %.1f min)",
    i, n_subj, s, as.numeric(Sys.time() - t0, units = "mins")
  ))
  out <- refit_fold(s, datExpr, mod_full, soft_power = soft_power)
  per_fold[[i]] <- out$projs
  stability[[i]] <- out$jaccard_summary |> mutate(fold_subject = s)
}
message(sprintf(
  "Done. Total elapsed: %.1f min",
  as.numeric(Sys.time() - t0, units = "mins")
))

cor <- stats::cor # restore base cor() now the WGCNA refit loop is done

# Aggregate per-pair LOSO AUCs
subj_pred <- function(module) {
  preds <- lapply(seq_len(n_subj), function(i) {
    p <- per_fold[[i]][[module]]
    if (is.null(p)) {
      return(NULL)
    }
    nm <- p$me_hold_names
    tibble(
      subject = common_subj[i],
      timepoint = sub(".*_(Pre|Post)$", "\\1", nm),
      me = p$me_hold
    )
  })
  bind_rows(preds)
}

score_pair <- function(i) {
  row_name <- top12$row[i]
  mod <- top12$module[i]
  pred <- subj_pred(mod)
  if (nrow(pred) == 0) {
    return(NULL)
  }

  wide <- pred |> pivot_wider(names_from = timepoint, values_from = me)
  if (!all(c("Pre", "Post") %in% names(wide))) {
    return(NULL)
  }
  wide <- wide |> mutate(Combined = (Pre + Post) / 2)

  if (row_name == "Age") {
    df <- tibble(subject = common_subj) |>
      left_join(wide, by = "subject") |>
      mutate(y = age_bin)
    df <- df |> filter(!is.na(Combined))
    if (length(unique(df$y)) < 2) {
      return(NULL)
    }
    r <- suppressMessages(roc(df$y, df$Combined, quiet = TRUE, direction = "auto"))
  } else if (row_name == "\u0394VL-responder") {
    df <- tibble(subject = common_subj) |>
      left_join(wide, by = "subject") |>
      mutate(y = vl_bin)
    df <- df |> filter(!is.na(Pre), !is.na(y))
    if (length(unique(df$y)) < 2) {
      return(NULL)
    }
    r <- suppressMessages(roc(df$y, df$Pre, quiet = TRUE, direction = "auto"))
  } else if (row_name == "LBM (High vs Low)") {
    df <- tibble(subject = common_subj) |>
      left_join(wide, by = "subject") |>
      mutate(y = lbm_bin)
    df <- df |> filter(!is.na(Pre), !is.na(y))
    if (length(unique(df$y)) < 2) {
      return(NULL)
    }
    r <- suppressMessages(roc(df$y, df$Pre, quiet = TRUE, direction = "auto"))
  } else if (row_name %in% c("Pre vs Post", "Pre vs Post (Young)", "Pre vs Post (Old)")) {
    df <- pred |> mutate(y = ifelse(timepoint == "Post", 1L, 0L))
    if (row_name == "Pre vs Post (Young)") {
      df <- df |> filter(subject %in% common_subj[age_bin == 0])
    } else if (row_name == "Pre vs Post (Old)") {
      df <- df |> filter(subject %in% common_subj[age_bin == 1])
    }
    df <- df |> filter(!is.na(me))
    if (length(unique(df$y)) < 2) {
      return(NULL)
    }
    r <- suppressMessages(roc(df$y, df$me, quiet = TRUE, direction = "auto"))
  } else {
    return(NULL)
  }

  ci <- as.numeric(ci.auc(r))
  tibble(
    row = row_name, module = mod,
    n_obs = length(r$cases) + length(r$controls),
    auc_insample = top12$auc[i],
    auc_loso_refit = as.numeric(auc(r)),
    ci_lo_refit = ci[1], ci_hi_refit = ci[3],
    drop = top12$auc[i] - as.numeric(auc(r))
  )
}

loso_refit <- bind_rows(lapply(seq_len(nrow(top12)), score_pair))
write_csv(loso_refit, file.path(DAT_OUT, "loso_wgcna_refit_summary.csv"))

stab_df <- bind_rows(stability) |>
  group_by(full) |>
  summarise(
    mean_jaccard = mean(jaccard),
    min_jaccard = min(jaccard),
    n_folds = n(),
    n_failed = sum(jaccard < 0.5)
  ) |>
  arrange(desc(mean_jaccard))
write_csv(stab_df, file.path(DAT_OUT, "loso_wgcna_refit_module_stability.csv"))

message("\n=== Full WGCNA-refit LOSO ===")
print(loso_refit |> dplyr::select(row, module, auc_insample, auc_loso_refit, drop))
message(sprintf("Median drop: %.3f", median(loso_refit$drop)))

message("\n=== Module stability (Jaccard of training vs full-sample assignments) ===")
print(stab_df)

# Plot
plot_df <- loso_refit |>
  mutate(
    pair = paste(module, row, sep = " | "),
    pair = factor(pair, levels = pair[order(-auc_insample)])
  )

p <- ggplot(plot_df) +
  geom_segment(
    aes(
      x = "in-sample", xend = "LOSO+refit",
      y = auc_insample, yend = auc_loso_refit, color = module
    ),
    linewidth = 1.2
  ) +
  geom_point(aes(x = "in-sample", y = auc_insample, color = module), size = 3.5) +
  geom_point(aes(x = "LOSO+refit", y = auc_loso_refit, color = module), size = 3.5) +
  geom_text(aes(x = "in-sample", y = auc_insample, label = sprintf("%.2f", auc_insample)),
    hjust = 1.2, size = 3.3
  ) +
  geom_text(aes(x = "LOSO+refit", y = auc_loso_refit, label = sprintf("%.2f", auc_loso_refit)),
    hjust = -0.2, size = 3.3
  ) +
  geom_hline(yintercept = 0.5, linetype = "dashed", color = "grey50") +
  scale_color_identity() +
  scale_y_continuous(limits = c(0.3, 1.0), breaks = seq(0.3, 1.0, 0.1)) +
  coord_cartesian(clip = "off") +
  facet_wrap(~pair, ncol = 3, scales = "free_y") +
  labs(
    title = "F07 Supp \u2014 Full WGCNA-refit LOSO sensitivity",
    subtitle = paste0(
      "Network refit on n-1 subjects per fold (power=12, same parameters as full-sample). ",
      "Training modules matched to full-sample modules by Jaccard overlap of ",
      "protein assignments; held-out subject's MEs computed via projection onto ",
      "training module PC1."
    ),
    x = NULL, y = "AUC",
    caption = sprintf(
      "Top-12 (row, module) pairs by in-sample perm_p. Median drop = %.3f.",
      median(loso_refit$drop)
    )
  ) +
  FIG_THEME +
  theme(
    plot.margin = margin(14, 18, 14, 18),
    strip.text = element_text(face = "bold", size = 9),
    panel.spacing = unit(8, "pt")
  )

W_in <- 12
H_in <- 10
ggsave(file.path(RPT_PDF, "SUPP_F07_loso_wgcna_refit.pdf"), p,
  width = W_in, height = H_in, units = "in", device = get_pdf_device()
)
ggsave(file.path(RPT_PNG, "SUPP_F07_loso_wgcna_refit.png"), p,
  width = W_in, height = H_in, units = "in", dpi = 300
)

message("\nWrote: SUPP_F07_loso_wgcna_refit.{pdf,png}")
message("Wrote: ", file.path(DAT_OUT, "loso_wgcna_refit_summary.csv"))
message("Wrote: ", file.path(DAT_OUT, "loso_wgcna_refit_module_stability.csv"))
