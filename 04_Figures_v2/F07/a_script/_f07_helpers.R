# F07 shared helpers
# Sourced by the F07 panel/supp scripts after their library() calls.
# Holds the pieces that were duplicated across scripts: PC1 projection,
# outcome-vector construction, LOOCV classifiers, and module pretty-naming.

# PC1 projection

# Fit the scaled-PC1 loading on a training matrix, oriented so the score
# tracks the row mean. Returns the loading plus the centring used to project
# new rows onto the same axis.
fit_pc1 <- function(train_mat) {
  Xtr <- scale(train_mat)
  center <- attr(Xtr, "scaled:center")
  sc <- attr(Xtr, "scaled:scale")
  sc[sc == 0 | !is.finite(sc)] <- 1
  sv <- svd(Xtr, nu = 0, nv = 1)
  v <- sv$v[, 1]
  if (stats::cor(as.numeric(Xtr %*% v), rowMeans(Xtr)) < 0) v <- -v
  list(v = v, center = center, scale = sc)
}

project_pc1 <- function(fit, newdata) {
  X <- sweep(newdata, 2, fit$center, "-")
  X <- sweep(X, 2, fit$scale, "/")
  as.numeric(X %*% fit$v)
}

# Outcome vectors

f07_age_bin <- function(subj_age, common_subj) {
  as.integer(subj_age$age[match(common_subj, subj_age$subject_key)] == "Old")
}

# Age-residualised median split on the full sample.
f07_resid_median_split <- function(y_raw, age_bin) {
  ok <- !is.na(y_raw)
  resid <- rep(NA_real_, length(y_raw))
  resid[ok] <- stats::residuals(stats::lm(y_raw[ok] ~ age_bin[ok]))
  as.integer(resid > stats::median(resid, na.rm = TRUE))
}

# Held-out version: for each subject the residualisation model and the split
# threshold are fit on the other subjects, then applied to that subject.
f07_loso_resid_median_split <- function(y_raw, age_bin) {
  vapply(seq_along(y_raw), function(i) {
    if (is.na(y_raw[i])) {
      return(NA_integer_)
    }
    train <- setdiff(which(!is.na(y_raw)), i)
    fit <- stats::lm(y_raw[train] ~ age_bin[train])
    cf <- stats::coef(fit)
    pred_i <- cf[1] + cf[2] * age_bin[i]
    resid_i <- y_raw[i] - pred_i
    thr <- stats::median(stats::residuals(fit), na.rm = TRUE)
    as.integer(resid_i > thr)
  }, integer(1))
}

# LOOCV classifiers

# LOOCV with nested top-k feature selection and logistic regression. The inner
# leave-one-out loop picks k by deviance, the outer loop scores the held-out
# subject. Returns out-of-fold probabilities plus selection bookkeeping.
run_topk_loocv <- function(labels, X, k_range = 2:5) {
  X <- as.matrix(X)
  n <- length(labels)
  probs <- numeric(n)
  best_k_log <- integer(n)
  selected <- vector("list", n)
  for (i in seq_len(n)) {
    train_x <- X[-i, , drop = FALSE]
    test_x <- X[i, , drop = FALSE]
    train_y <- labels[-i]
    n_train <- length(train_y)

    train_cors <- abs(stats::cor(train_x, train_y))
    train_cors[is.na(train_cors)] <- 0
    ranked <- names(sort(train_cors[, 1], decreasing = TRUE))

    inner_dev <- setNames(numeric(length(k_range)), as.character(k_range))
    for (k in k_range) {
      feats <- ranked[seq_len(min(k, length(ranked)))]
      dev <- 0
      for (j in seq_len(n_train)) {
        fit <- tryCatch(
          suppressWarnings(stats::glm(y ~ .,
            family = stats::binomial,
            data = cbind(
              y = train_y[-j],
              as.data.frame(train_x[-j, feats, drop = FALSE])
            )
          )),
          error = function(e) NULL
        )
        if (is.null(fit)) {
          dev <- dev + log(2)
          next
        }
        p_hat <- stats::predict(fit,
          type = "response",
          newdata = as.data.frame(train_x[j, feats, drop = FALSE])
        )
        p_hat <- pmin(pmax(p_hat, 1e-6), 1 - 1e-6)
        dev <- dev - (train_y[j] * log(p_hat) + (1 - train_y[j]) * log(1 - p_hat))
      }
      inner_dev[as.character(k)] <- dev
    }

    best_k <- as.integer(names(which.min(inner_dev)))
    best_k_log[i] <- best_k
    feats <- ranked[seq_len(best_k)]
    selected[[i]] <- feats

    fit <- tryCatch(
      suppressWarnings(stats::glm(y ~ .,
        family = stats::binomial,
        data = cbind(
          y = train_y,
          as.data.frame(train_x[, feats, drop = FALSE])
        )
      )),
      error = function(e) NULL
    )
    probs[i] <- if (!is.null(fit)) {
      stats::predict(fit,
        type = "response",
        newdata = as.data.frame(test_x[, feats, drop = FALSE])
      )
    } else {
      0.5
    }
  }

  all_features <- colnames(X)
  feature_freq <- setNames(integer(length(all_features)), all_features)
  for (s in selected) feature_freq[s] <- feature_freq[s] + 1L

  list(
    probs = probs, feature_freq = feature_freq,
    selected = selected, best_k_log = best_k_log, ks = best_k_log
  )
}

# Fixed-k LOOCV that returns only the AUC. Used inside permutation nulls.
fast_loocv_auc <- function(labels, X, k) {
  X <- as.matrix(X)
  n <- length(labels)
  probs <- numeric(n)
  for (i in seq_len(n)) {
    train_x <- X[-i, , drop = FALSE]
    test_x <- X[i, , drop = FALSE]
    train_y <- labels[-i]
    train_cors <- abs(stats::cor(train_x, train_y))
    train_cors[is.na(train_cors)] <- 0
    feats <- names(sort(train_cors[, 1], decreasing = TRUE))[seq_len(min(k, ncol(X)))]
    fit <- tryCatch(
      suppressWarnings(stats::glm(y ~ .,
        family = stats::binomial,
        data = cbind(
          y = train_y,
          as.data.frame(train_x[, feats, drop = FALSE])
        )
      )),
      error = function(e) NULL
    )
    probs[i] <- if (!is.null(fit)) {
      stats::predict(fit,
        type = "response",
        newdata = as.data.frame(test_x[, feats, drop = FALSE])
      )
    } else {
      0.5
    }
  }
  tryCatch(as.numeric(pROC::auc(pROC::roc(labels, probs, quiet = TRUE))),
    error = function(e) 0.5
  )
}

# Module naming

# Build a module -> biological-label translator bound to one label vector.
make_pretty_mod <- function(mod_bio_labels) {
  function(m) {
    short <- gsub("^ME", "", m)
    bio <- mod_bio_labels[short]
    if (is.na(bio)) stringr::str_to_title(short) else unname(bio)
  }
}
