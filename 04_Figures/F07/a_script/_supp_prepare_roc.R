#!/usr/bin/env Rscript
# ROC pilot — screen 7 candidate classifiers using panel A's framework
# (LOOCV + top-k feature selection + logistic regression + permutation test).
# Output: classifier_pilot_summary.csv + classifier_pilot_curves.csv
#
# Sourced by 01_main_panels.R — expects figure_supplement_helpers.R already loaded.

pacman::p_load(tidyverse, pROC)

OUT <- "04_Figures/F07/c_data"
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)

# Load data
F06_SUPP <- "04_Figures/F06/c_data/F06_supplementary.xlsx"
stopifnot("F06 stitcher must run first: missing F06_supplementary.xlsx" =
  file.exists(F06_SUPP))
MEs     <- read_matrix_sheet(F06_SUPP, "MEs",     "sample_id")
me_pre  <- read_matrix_sheet(F06_SUPP, "me_pre",  "subject_key")
me_post <- read_matrix_sheet(F06_SUPP, "me_post", "subject_key")
subj_age<- read_sheet_df(F06_SUPP, "metadata_subj_age")
common_subj <- read_vector_sheet(F06_SUPP, "common_subj")
pheno   <- read_sheet_df(F06_SUPP, "metadata_pheno_wide")
mods    <- read_sheet_df(F06_SUPP, "WGCNA_module_assignments")

dep <- read_csv("03_DEP/c_data/03_combined_results.csv", show_col_types = FALSE)
imp <- read_csv("02_imputation/c_data/01_imputed.csv",   show_col_types = FALSE)

# Labels
true_age <- ifelse(subj_age$age[match(common_subj, subj_age$subject_key)] == "Old", 1, 0)
n_subj   <- length(common_subj)

# Core LOOCV engine (same as panel A)
run_topk_loocv <- function(labels, X, k_range = 2:5) {
  X <- as.matrix(X); n <- length(labels); probs <- numeric(n); ks <- integer(n)
  for (i in seq_len(n)) {
    tr_x <- X[-i,,drop=FALSE]; te_x <- X[i,,drop=FALSE]; tr_y <- labels[-i]
    r <- abs(cor(tr_x, tr_y)); r[is.na(r)] <- 0
    ranked <- names(sort(r[,1], decreasing = TRUE))
    # inner LOOCV
    inner <- setNames(numeric(length(k_range)), as.character(k_range))
    for (k in k_range) {
      ku <- min(k, length(ranked)); feats <- ranked[seq_len(ku)]
      dev <- 0
      for (j in seq_along(tr_y)) {
        fit <- tryCatch(suppressWarnings(glm(y~., binomial,
                 data = cbind(y=tr_y[-j], as.data.frame(tr_x[-j, feats, drop=FALSE])))),
                 error=function(e) NULL)
        if (is.null(fit)) { dev <- dev + log(2); next }
        p <- predict(fit, type="response", newdata=as.data.frame(tr_x[j,feats,drop=FALSE]))
        p <- pmin(pmax(p,1e-6),1-1e-6)
        dev <- dev - (tr_y[j]*log(p) + (1-tr_y[j])*log(1-p))
      }
      inner[as.character(k)] <- dev
    }
    best_k <- as.integer(names(which.min(inner))); ks[i] <- best_k
    feats <- ranked[seq_len(best_k)]
    fit <- tryCatch(suppressWarnings(glm(y~., binomial,
             data=cbind(y=tr_y, as.data.frame(tr_x[,feats,drop=FALSE])))),
             error=function(e) NULL)
    probs[i] <- if (!is.null(fit))
      predict(fit, type="response", newdata=as.data.frame(te_x[,feats,drop=FALSE])) else 0.5
  }
  list(probs=probs, ks=ks)
}

fast_loocv_auc <- function(labels, X, k) {
  X <- as.matrix(X); n <- length(labels); probs <- numeric(n)
  for (i in seq_len(n)) {
    tr_x <- X[-i,,drop=FALSE]; te_x <- X[i,,drop=FALSE]; tr_y <- labels[-i]
    r <- abs(cor(tr_x, tr_y)); r[is.na(r)] <- 0
    feats <- names(sort(r[,1], decreasing=TRUE))[seq_len(min(k, ncol(X)))]
    fit <- tryCatch(suppressWarnings(glm(y~., binomial,
             data=cbind(y=tr_y, as.data.frame(tr_x[,feats,drop=FALSE])))),
             error=function(e) NULL)
    probs[i] <- if (!is.null(fit))
      predict(fit, type="response", newdata=as.data.frame(te_x[,feats,drop=FALSE])) else 0.5
  }
  tryCatch(as.numeric(auc(roc(labels, probs, quiet=TRUE))), error=function(e) 0.5)
}

eval_clf <- function(name, labels, X, k_range, n_perm = 200) {
  message(sprintf("  [%s] n_feat=%d, n_obs=%d", name, ncol(X), nrow(X)))
  res <- run_topk_loocv(labels, X, k_range)
  roc_obj <- tryCatch(roc(labels, res$probs, quiet=TRUE), error=function(e) NULL)
  if (is.null(roc_obj)) return(NULL)
  auc_v <- as.numeric(auc(roc_obj)); ci <- as.numeric(ci.auc(roc_obj))
  k_use <- as.integer(median(res$ks))
  nulls <- replicate(n_perm, fast_loocv_auc(sample(labels), X, k_use))
  list(name=name, auc=auc_v, ci_lo=ci[1], ci_hi=ci[3],
       perm_p=(sum(nulls >= auc_v) + 1) / (n_perm + 1), probs=res$probs, labels=labels,
       n=length(labels), n_pos=sum(labels==1), n_neg=sum(labels==0),
       k_med=k_use, fpr=1-roc_obj$specificities, tpr=roc_obj$sensitivities)
}

# Build feature matrices
sample_cols <- setdiff(colnames(imp), c("protein","uniprot_id","gene","description"))
X_prot_all <- t(as.matrix(imp[, sample_cols]))
rownames(X_prot_all) <- sample_cols
colnames(X_prot_all) <- as.character(imp[[1]])

subj_of <- sub("_(Pre|Post)$", "", sample_cols)
X_prot_subj <- do.call(rbind, lapply(common_subj, function(s) {
  colMeans(X_prot_all[subj_of == s, , drop=FALSE], na.rm=TRUE)
}))
rownames(X_prot_subj) <- common_subj
sds <- apply(X_prot_subj, 2, sd, na.rm=TRUE)
X_prot_subj <- X_prot_subj[, sds > 0 & is.finite(sds)]

X_pheno <- pheno |> dplyr::filter(subject_key %in% common_subj) |>
  dplyr::select(subject_key, VL_Pre, LBM_Pre, BMI_Pre) |>
  column_to_rownames("subject_key")
X_pheno <- X_pheno[common_subj, ]
ok_pheno <- complete.cases(X_pheno)

me_avg <- (as.matrix(me_pre[common_subj,]) + as.matrix(me_post[common_subj,]))/2
X_mepheno <- cbind(me_avg[ok_pheno,], as.matrix(X_pheno[ok_pheno,]))

pheno_subj <- pheno |> filter(subject_key %in% common_subj) |>
  mutate(age_bin = ifelse(subject_key %in% subj_age$subject_key[subj_age$age=="Old"],1,0))
pheno_subj <- pheno_subj[match(common_subj, pheno_subj$subject_key), ]

make_responder <- function(y_raw, age_bin) {
  ok <- !is.na(y_raw)
  if (sum(ok) < 10) return(list(lab=rep(NA, length(y_raw)), ok=ok))
  res <- rep(NA_real_, length(y_raw))
  res[ok] <- residuals(lm(y_raw[ok] ~ age_bin[ok]))
  lab <- rep(NA_integer_, length(y_raw))
  lab[ok] <- as.integer(res[ok] > median(res[ok], na.rm=TRUE))
  list(lab=lab, ok=ok)
}
resp_vl  <- make_responder(pheno_subj$delta_VL,  pheno_subj$age_bin)
resp_lbm <- make_responder(pheno_subj$delta_LBM, pheno_subj$age_bin)

# Run classifiers
set.seed(42)
results <- list()

message("[1/7] Age ~ top-k proteins (avg Pre+Post)")
results$prot_age <- eval_clf("Age|proteins",
  true_age, X_prot_subj, k_range = c(5,10,20))

message("[2/7] Age ~ phenotype (VL/LBM/BMI pre)")
results$pheno_age <- eval_clf("Age|pheno",
  true_age[ok_pheno], as.matrix(X_pheno[ok_pheno,]), k_range = 2:3)

message("[3/7] Age ~ MEs + phenotype")
results$mepheno_age <- eval_clf("Age|ME+pheno",
  true_age[ok_pheno], X_mepheno, k_range = 2:5)

message("[4/7] delta-VL responder ~ baseline MEs")
ok4 <- !is.na(resp_vl$lab)
results$resp_vl <- eval_clf("\u0394VL-resp|ME_Pre",
  resp_vl$lab[ok4], as.matrix(me_pre[common_subj,])[ok4,], k_range = 2:5)

message("[5/7] delta-LBM responder ~ baseline MEs")
ok5 <- !is.na(resp_lbm$lab)
results$resp_lbm <- eval_clf("\u0394LBM-resp|ME_Pre",
  resp_lbm$lab[ok5], as.matrix(me_pre[common_subj,])[ok5,], k_range = 2:5)

# Protein-level classifiers (10-fold CV, n=2113)
message("[6/7] Aging DEP vs non-DEP (protein-intrinsic features)")
prot_feat <- dep |>
  transmute(uniprot_id,
            is_aging_dep = as.integer(sig_pi_Aging != 0),
            abs_lfc_TY   = abs(logFC_Training_Young)) |>
  left_join(mods |> dplyr::select(uniprot_id, module_color), by="uniprot_id")

abund_stats <- tibble(
  uniprot_id = as.character(imp[[1]]),
  mean_abund = rowMeans(as.matrix(imp[,sample_cols]), na.rm=TRUE),
  sd_abund   = apply(as.matrix(imp[,sample_cols]), 1, sd, na.rm=TRUE))
prot_feat <- prot_feat |> left_join(abund_stats, by="uniprot_id") |>
  filter(!is.na(is_aging_dep), !is.na(mean_abund), !is.na(sd_abund),
         !is.na(module_color))
mod_oh <- model.matrix(~ module_color - 1, data=prot_feat)
X_prot6 <- cbind(mean_abund=prot_feat$mean_abund,
                 sd_abund=prot_feat$sd_abund,
                 abs_lfc_TY=prot_feat$abs_lfc_TY,
                 mod_oh)
y6 <- prot_feat$is_aging_dep

# Stratified 10-fold logistic CV + permutation null. Used for classifiers 6 and 7.
fit_10fold_perm <- function(X, y, n_perm = 200, seed = 42) {
  set.seed(seed)
  folds <- integer(length(y))
  for (cls in unique(y)) {
    idx <- which(y == cls)
    folds[idx] <- sample(rep(1:10, length.out = length(idx)))
  }
  probs <- numeric(length(y))
  for (f in 1:10) {
    tr  <- folds != f
    fit <- tryCatch(suppressWarnings(glm(y~., binomial,
      data=cbind(y=y[tr], as.data.frame(X[tr,])))), error=function(e) NULL)
    if (!is.null(fit)) probs[!tr] <- predict(fit, type="response",
      newdata=as.data.frame(X[!tr,]))
  }
  roc_obj <- roc(y, probs, quiet=TRUE)
  null_aucs <- numeric(n_perm)
  for (p in seq_len(n_perm)) {
    ys <- sample(y); pr <- numeric(length(ys))
    for (f in 1:10) {
      tr  <- folds != f
      fit <- tryCatch(suppressWarnings(glm(y~., binomial,
        data=cbind(y=ys[tr], as.data.frame(X[tr,])))), error=function(e) NULL)
      if (!is.null(fit)) pr[!tr] <- predict(fit, type="response",
        newdata=as.data.frame(X[!tr,]))
    }
    null_aucs[p] <- tryCatch(as.numeric(auc(roc(ys, pr, quiet=TRUE))),
                              error=function(e) 0.5)
  }
  auc_val <- as.numeric(auc(roc_obj))
  ci_vals <- as.numeric(ci.auc(roc_obj))
  list(
    auc    = auc_val,
    ci_lo  = ci_vals[1], ci_hi = ci_vals[3],
    perm_p = (sum(null_aucs >= auc_val) + 1) / (n_perm + 1),
    probs  = probs, labels = y,
    n      = length(y), n_pos = sum(y == 1), n_neg = sum(y == 0),
    k_med  = ncol(X),
    fpr    = 1 - roc_obj$specificities, tpr = roc_obj$sensitivities)
}

cv6 <- fit_10fold_perm(X_prot6, y6, n_perm = 200)
results$prot_dep <- c(list(name = "AgingDEP|prot_features"), cv6)

message("[7/7] Reversed vs Exacerbated among aging DEPs")
rev_df <- dep |> filter(sig_pi_Aging != 0) |>
  transmute(uniprot_id,
            reversed = as.integer(sign(logFC_Aging) != sign(logFC_Training_Old)),
            abs_lfc_aging = abs(logFC_Aging),
            abs_lfc_TY   = abs(logFC_Training_Young)) |>
  left_join(mods |> dplyr::select(uniprot_id, module_color), by="uniprot_id") |>
  left_join(abund_stats, by="uniprot_id") |>
  (\(df) df[complete.cases(df), ])()

if (nrow(rev_df) > 30 && length(unique(rev_df$reversed)) == 2) {
  mod_oh7 <- model.matrix(~ module_color - 1, data=rev_df)
  X7 <- cbind(abs_lfc_aging=rev_df$abs_lfc_aging,
              abs_lfc_TY=rev_df$abs_lfc_TY,
              mean_abund=rev_df$mean_abund, sd_abund=rev_df$sd_abund, mod_oh7)
  y7 <- rev_df$reversed
  cv7 <- fit_10fold_perm(X7, y7, n_perm = 200)
  results$rev_vs_exa <- c(list(name = "Reversed|prot_features"), cv7)
}

# Summary table
summ <- map_dfr(results, function(r) {
  tibble(classifier=r$name, n=r$n, n_pos=r$n_pos, n_neg=r$n_neg,
         k_med=r$k_med, auc=r$auc, ci_lo=r$ci_lo, ci_hi=r$ci_hi,
         perm_p=r$perm_p)
})
print(summ |> arrange(desc(auc)))
write_csv(summ, file.path(OUT, "classifier_pilot_summary.csv"))

curves <- map_dfr(results, function(r)
  tibble(classifier=r$name, fpr=r$fpr, tpr=r$tpr))
write_csv(curves, file.path(OUT, "classifier_pilot_curves.csv"))

message("Done. See: ", OUT)
