# ROC pilot — screen 7 candidate classifiers using panel A's framework
# (LOOCV + top-k feature selection + logistic regression + permutation test).
# Output: roc_pilot_summary.csv + roc_pilot_curves.csv + roc_pilot_composite.png
#
# Methods references:
#   Ambroise & McLachlan 2002 PNAS — feature selection inside CV folds
#   Varma & Simon 2006 BMC Bioinformatics — nested CV for small-N
#   Saeys et al. 2007 Bioinformatics — filter + wrapper hybrid feature selection

suppressPackageStartupMessages({
  library(tidyverse); library(pROC); library(patchwork)
})
setwd(rprojroot::find_rstudio_root_file())
source("04_Figures/shared/figure_supplement_helpers.R")  # read_matrix_sheet / read_sheet_df / read_vector_sheet

OUT <- "04_Figures/F07/c_data"
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)

# ---- Load data --------------------------------------------------------------
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
imp <- read_csv("02_Imputation/c_data/01_imputed.csv",   show_col_types = FALSE)

# Labels
true_age <- ifelse(subj_age$age[match(common_subj, subj_age$subject_key)] == "Old", 1, 0)
n_subj   <- length(common_subj)

# ---- Core LOOCV engine (same as panel A) -------------------------------------
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
      for (j in seq_len(length(tr_y))) {
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

# ---- Build feature matrices --------------------------------------------------
# Subject-averaged protein abundance (Pre+Post)/2, matching panel A "Combined"
prot_ids <- imp$protein %||% imp[[1]]
sample_cols <- setdiff(colnames(imp), c("protein","uniprot_id","gene","description"))
X_prot_all <- t(as.matrix(imp[, sample_cols]))  # rows=samples, cols=proteins
rownames(X_prot_all) <- sample_cols
colnames(X_prot_all) <- as.character(imp[[1]])

# subject mapping: strip _Pre/_Post
subj_of <- sub("_(Pre|Post)$", "", sample_cols)
# average Pre+Post per subject, restrict to common_subj
X_prot_subj <- do.call(rbind, lapply(common_subj, function(s) {
  colMeans(X_prot_all[subj_of == s, , drop=FALSE], na.rm=TRUE)
}))
rownames(X_prot_subj) <- common_subj
# remove zero-variance proteins
sds <- apply(X_prot_subj, 2, sd, na.rm=TRUE)
X_prot_subj <- X_prot_subj[, sds > 0 & is.finite(sds)]

# Phenotype (baseline only) — VL_Pre, LBM_Pre, BMI_Pre
# Namespace select() because MASS::select gets masked over dplyr::select when
# this script is sourced after panel_A/panel_B in 90_stitch_figure.R.
X_pheno <- pheno %>% dplyr::filter(subject_key %in% common_subj) %>%
  dplyr::select(subject_key, VL_Pre, LBM_Pre, BMI_Pre) %>%
  column_to_rownames("subject_key")
X_pheno <- X_pheno[common_subj, ]
# drop rows with NA -> use complete cases for pheno models
ok_pheno <- complete.cases(X_pheno)

# MEs (combined = avg Pre+Post)
me_avg <- (as.matrix(me_pre[common_subj,]) + as.matrix(me_post[common_subj,]))/2

# ME + phenotype combined
X_mepheno <- cbind(me_avg[ok_pheno,], as.matrix(X_pheno[ok_pheno,]))

# Responder labels (median split of delta_VL, delta_LBM), age-residualized
pheno_subj <- pheno %>% filter(subject_key %in% common_subj) %>%
  mutate(age_bin = ifelse(subject_key %in% subj_age$subject_key[subj_age$age=="Old"],1,0))
pheno_subj <- pheno_subj[match(common_subj, pheno_subj$subject_key), ]

# delta_VL responder: residualize out age, then median split
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

# ---- Run classifiers ---------------------------------------------------------
set.seed(42)
results <- list()

message("[1/7] Age ~ top-k proteins (avg Pre+Post)")
results$prot_age <- eval_clf("Age|proteins",
  true_age, X_prot_subj, k_range = c(5,10,20))

message("[2/7] Age ~ phenotype (VL/LBM/BMI pre)")
# use only subjects with complete pheno
results$pheno_age <- eval_clf("Age|pheno",
  true_age[ok_pheno], as.matrix(X_pheno[ok_pheno,]), k_range = 2:3)

message("[3/7] Age ~ MEs + phenotype")
results$mepheno_age <- eval_clf("Age|ME+pheno",
  true_age[ok_pheno], X_mepheno, k_range = 2:5)

message("[4/7] Δ VL responder ~ baseline MEs")
ok4 <- !is.na(resp_vl$lab)
results$resp_vl <- eval_clf("ΔVL-resp|ME_Pre",
  resp_vl$lab[ok4], as.matrix(me_pre[common_subj,])[ok4,], k_range = 2:5)

message("[5/7] Δ LBM responder ~ baseline MEs")
ok5 <- !is.na(resp_lbm$lab)
results$resp_lbm <- eval_clf("ΔLBM-resp|ME_Pre",
  resp_lbm$lab[ok5], as.matrix(me_pre[common_subj,])[ok5,], k_range = 2:5)

# ---- Protein-level classifiers (10-fold CV, n=2113) --------------------------
message("[6/7] Aging DEP vs non-DEP (protein-intrinsic features)")
# Features: mean abundance, SD, missingness, module color (one-hot), |TY logFC|
prot_feat <- dep %>%
  transmute(uniprot_id,
            is_aging_dep = as.integer(sig_pi_Aging != 0),
            abs_lfc_TY   = abs(logFC_Training_Young)) %>%
  left_join(mods %>% dplyr::select(uniprot_id, module_color), by="uniprot_id")

# add abundance summary from imputed matrix
abund_stats <- tibble(
  uniprot_id = as.character(imp[[1]]),
  mean_abund = rowMeans(as.matrix(imp[,sample_cols]), na.rm=TRUE),
  sd_abund   = apply(as.matrix(imp[,sample_cols]), 1, sd, na.rm=TRUE))
prot_feat <- prot_feat %>% left_join(abund_stats, by="uniprot_id") %>%
  filter(!is.na(is_aging_dep), !is.na(mean_abund), !is.na(sd_abund),
         !is.na(module_color))
# module one-hot (drop grey as reference)
mod_oh <- model.matrix(~ module_color - 1, data=prot_feat)
X_prot6 <- cbind(mean_abund=prot_feat$mean_abund,
                 sd_abund=prot_feat$sd_abund,
                 abs_lfc_TY=prot_feat$abs_lfc_TY,
                 mod_oh)
y6 <- prot_feat$is_aging_dep

# 10-fold CV logistic (no feature selection needed; small feature set)
set.seed(42)
folds <- sample(rep(1:10, length.out=length(y6)))
probs6 <- numeric(length(y6))
for (f in 1:10) {
  tr <- folds != f
  fit <- suppressWarnings(glm(y~., binomial,
    data=cbind(y=y6[tr], as.data.frame(X_prot6[tr,]))))
  probs6[!tr] <- predict(fit, type="response",
    newdata=as.data.frame(X_prot6[!tr,]))
}
roc6 <- roc(y6, probs6, quiet=TRUE)
# permutation
n_perm6 <- 200
null6 <- numeric(n_perm6)
for (p in 1:n_perm6) {
  ys <- sample(y6); pr <- numeric(length(ys))
  for (f in 1:10) {
    tr <- folds != f
    fit <- tryCatch(suppressWarnings(glm(y~., binomial,
      data=cbind(y=ys[tr], as.data.frame(X_prot6[tr,])))), error=function(e) NULL)
    if (!is.null(fit)) pr[!tr] <- predict(fit, type="response",
      newdata=as.data.frame(X_prot6[!tr,]))
  }
  null6[p] <- tryCatch(as.numeric(auc(roc(ys, pr, quiet=TRUE))), error=function(e) 0.5)
}
results$prot_dep <- list(
  name="AgingDEP|prot_features",
  auc=as.numeric(auc(roc6)),
  ci_lo=as.numeric(ci.auc(roc6))[1], ci_hi=as.numeric(ci.auc(roc6))[3],
  perm_p=(sum(null6 >= as.numeric(auc(roc6))) + 1) / (n_perm6 + 1),
  probs=probs6, labels=y6,
  n=length(y6), n_pos=sum(y6==1), n_neg=sum(y6==0), k_med=ncol(X_prot6),
  fpr=1-roc6$specificities, tpr=roc6$sensitivities)

message("[7/7] Reversed vs Exacerbated among aging DEPs")
rev_df <- dep %>% filter(sig_pi_Aging != 0) %>%
  transmute(uniprot_id,
            reversed = as.integer(sign(logFC_Aging) != sign(logFC_Training_Old)),
            abs_lfc_aging = abs(logFC_Aging),
            abs_lfc_TY   = abs(logFC_Training_Young)) %>%
  left_join(mods %>% dplyr::select(uniprot_id, module_color), by="uniprot_id") %>%
  left_join(abund_stats, by="uniprot_id") %>%
  filter(complete.cases(.))

if (nrow(rev_df) > 30 && length(unique(rev_df$reversed)) == 2) {
  mod_oh7 <- model.matrix(~ module_color - 1, data=rev_df)
  X7 <- cbind(abs_lfc_aging=rev_df$abs_lfc_aging,
              abs_lfc_TY=rev_df$abs_lfc_TY,
              mean_abund=rev_df$mean_abund, sd_abund=rev_df$sd_abund, mod_oh7)
  y7 <- rev_df$reversed
  set.seed(42); folds7 <- sample(rep(1:10, length.out=length(y7)))
  probs7 <- numeric(length(y7))
  for (f in 1:10) {
    tr <- folds7 != f
    fit <- tryCatch(suppressWarnings(glm(y~., binomial,
      data=cbind(y=y7[tr], as.data.frame(X7[tr,])))), error=function(e) NULL)
    if (!is.null(fit)) probs7[!tr] <- predict(fit, type="response",
      newdata=as.data.frame(X7[!tr,]))
  }
  roc7 <- roc(y7, probs7, quiet=TRUE)
  null7 <- numeric(200)
  for (p in 1:200) {
    ys <- sample(y7); pr <- numeric(length(ys))
    for (f in 1:10) {
      tr <- folds7 != f
      fit <- tryCatch(suppressWarnings(glm(y~., binomial,
        data=cbind(y=ys[tr], as.data.frame(X7[tr,])))), error=function(e) NULL)
      if (!is.null(fit)) pr[!tr] <- predict(fit, type="response",
        newdata=as.data.frame(X7[!tr,]))
    }
    null7[p] <- tryCatch(as.numeric(auc(roc(ys, pr, quiet=TRUE))), error=function(e) 0.5)
  }
  results$rev_vs_exa <- list(
    name="Reversed|prot_features", auc=as.numeric(auc(roc7)),
    ci_lo=as.numeric(ci.auc(roc7))[1], ci_hi=as.numeric(ci.auc(roc7))[3],
    perm_p=(sum(null7 >= as.numeric(auc(roc7))) + 1) / (length(null7) + 1),
    probs=probs7, labels=y7,
    n=length(y7), n_pos=sum(y7==1), n_neg=sum(y7==0), k_med=ncol(X7),
    fpr=1-roc7$specificities, tpr=roc7$sensitivities)
}

# ---- Summary table -----------------------------------------------------------
summ <- map_dfr(results, function(r) {
  tibble(classifier=r$name, n=r$n, n_pos=r$n_pos, n_neg=r$n_neg,
         k_med=r$k_med, auc=r$auc, ci_lo=r$ci_lo, ci_hi=r$ci_hi,
         perm_p=r$perm_p)
})
print(summ %>% arrange(desc(auc)))
write_csv(summ, file.path(OUT, "roc_pilot_summary.csv"))

# ROC curves long-form
curves <- map_dfr(results, function(r)
  tibble(classifier=r$name, fpr=r$fpr, tpr=r$tpr))
write_csv(curves, file.path(OUT, "roc_pilot_curves.csv"))

# ---- Composite plot ----------------------------------------------------------
curves %>%
  ggplot(aes(fpr, tpr, color=classifier)) +
  geom_abline(slope=1, intercept=0, linetype="dashed", color="grey60") +
  geom_line(linewidth=1) +
  geom_text(data=summ, aes(x=0.6, y=0.1,
    label=sprintf("AUC=%.2f (p=%.3f)", auc, perm_p)), size=3, hjust=0,
    inherit.aes=FALSE, color="grey20") +
  facet_wrap(~classifier, ncol=3) +
  coord_fixed() + theme_bw() +
  labs(x="1 − Specificity", y="Sensitivity",
       title="ROC pilot — 7 candidate classifiers") +
  theme(legend.position="none", strip.text=element_text(face="bold", size=9))
# Debug composite omitted — pilot data consumed by xlsx only

message("Done. See: ", OUT)
