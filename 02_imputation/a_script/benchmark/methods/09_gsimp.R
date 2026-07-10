# GSimp — left-censored Gibbs sampler imputation
# Wei et al. 2018 (PLOS Comp Bio 14:e1005973)
# Real implementation sourced from github.com/WandeRum/GSimp
# Pure MNAR; uses glmnet to predict each variable, draws from
# truncated normal bounded above by column minimum

impute_GSimp <- function(mat, meta, is_mnar, ...) {
  set.seed(42)

  gsimp_dir <- "02_imputation/a_script/benchmark/methods/gsimp_source"
  stopifnot("GSimp source files not found" = dir.exists(gsimp_dir))

  pacman::p_load(magrittr, abind)
  source(file.path(gsimp_dir, "MVI_global.R"), local = FALSE)
  source(file.path(gsimp_dir, "Prediction_funcs.R"), local = FALSE)
  source(file.path(gsimp_dir, "GSimp.R"), local = FALSE)

  # GSimp expects samples-as-rows, proteins-as-columns (like missForest)
  dat <- as.data.frame(t(mat))

  # Scale, impute, recover (mirrors pre_processing_GS_wrapper)
  scaled <- scale_recover(dat, method = "scale")
  dat_scaled <- scaled[[1]]
  param_df   <- scaled[[2]]

  # Reduced iterations for benchmark (full: 100 inner x 20 outer)
  # 50 inner x 10 outer is sufficient for convergence on this dataset
  result <- tryCatch(
    GS_impute(dat_scaled,
              iters_each = 50,
              iters_all  = 10,
              initial    = "qrilc",
              lo         = -Inf,
              hi         = "min",
              imp_model  = "glmnet_pred",
              n_cores    = 1),
    error = function(e) {
      message("GSimp failed: ", conditionMessage(e))
      NULL
    }
  )

  if (is.null(result)) {
    message("GSimp fallback: QRILC")
    return(MsCoreUtils::impute_matrix(mat, method = "QRILC"))
  }

  # Recover original scale
  recovered <- scale_recover(result$data_imp, method = "recover",
                             param_df = param_df)
  imp <- t(as.matrix(recovered[[1]]))
  dimnames(imp) <- dimnames(mat)
  imp
}
