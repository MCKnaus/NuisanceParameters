#' Ensemble learner training for short stacking
#'
#' @description
#' \code{\link{ensemble_short}} trains an ensemble learner on the given training data
#' in the spirit of a short-stacking procedure.
#'
#' @param ml List of methods built via \code{\link{create_method}} to be used in
#' ensemble model.
#' @param x Covariate matrix of training sample.
#' @param y Vector of outcomes of training sample.
#' @param cf_mat Logical matrix with k columns of indicators representing the different folds
#' (for example created by \code{\link{prep_cf_mat}}).
#' @param storeModels Specifies whether and where to store the models.
#' If smoother weights should be computed in a second step (using the weights method),
#' the models need to be stored somewhere (either in memory or on disk).
#' Value has to be an item from \code{c("No", "Memory", "Disk")} with "No" being
#' the default value.
#' @param path Optional path to save the list containing the ensemble models
#' from each cross-fitting fold to. Only considered, if \code{storeModels = "Disk"}.
#' @param quiet If FALSE, method that is currently computed is printed into the
#' console.
#'
#' @return List object containing:
#' \item{fit_cv}{matrix of dimension \code{nrow(x)} x \code{length(ml)} containing
#' the cross-fitted predictions of the machine learning methods from the ensemble}
#' \item{nnls_weights}{the weights that each machine learning method receives in
#' the ensemble}
#' \item{mse_cv}{cross-validated MSEs for each machine learning  method}
#' \item{ml}{list of fitted machine learning models from all cross-fitting folds
#' (if storeModels is not set to "No")}
#'
#' @export
#'
ensemble_short = function(ml,
                          x, y,
                          cf_mat,
                          storeModels = c("No", "Memory", "Disk"),
                          path = NULL,
                          quiet = TRUE) {

  # storage configuration
  storeModels = match.arg(storeModels)
  saveModels = (storeModels != "No")

  s = NULL
  if(saveModels) s = vector("list", length = ncol(cf_mat))

  # configure matrix to store predictions
  fit_cv = matrix(NA, nrow = nrow(x), ncol = length(ml))
  colnames(fit_cv) = sprintf("Method%s", seq(1:length(ml)))
  for (i in 1:length(ml)) {
    if (!is.null(ml[[i]]$name)) colnames(fit_cv)[i] = ml[[i]]$name
  }

  for (i in 1:ncol(cf_mat)) {

    if (isFALSE(quiet)) print(paste("Cross-fitting fold: ", toString(i)))

    fold = cf_mat[, i]
    x_tr = x[!fold, ]
    y_tr = y[!fold]
    x_te = x[fold, ]

    ml_fit = ensemble_core(ml, x_tr, y_tr, quiet = quiet)
    fit_cv[fold,] = predict.ensemble_core(ml_fit, ml, x_tr, y_tr, x_te, quiet = quiet)

    if(saveModels) s[[i]] = ml_fit

  }

  if(storeModels == "Disk") {
    if(is.null(path)) path = getwd()
    file_name = paste0(path, "/ensemble_short_", format(Sys.time(), "%Y%m%d%H%M%S"), ".rds")
    saveRDS(s, file_name)
    s = file_name
  }

  fit_cv[is.na(fit_cv)] = mean(y)
  mse_cv = colMeans((c(y) - fit_cv)^2)

  nnls_w = nnls_weights(X = fit_cv, y = y)

  names(mse_cv) = names(nnls_w) = colnames(fit_cv)

  output = list(
    "fit_cv" = fit_cv,
    "nnls_weights" = nnls_w,
    "mse_cv" = mse_cv,
    "ml" = s
  )

  class(output) = "ensemble_short"

  return(output)

}


#' Prediction of fitted values for short-stacked ensemble learner
#'
#' @description
#' Returns prediction of fitted values for a fitted and short-stacked ensemble
#' learner created by \code{\link{ensemble_short}}.
#'
#' @param object Short-stacked ensemble learner from \code{\link{ensemble_short}}.
#' @param ... Ignore unused arguments
#'
#' @return Vector containing ensemble predictions.
#'
#' @export
#'
#' @method predict ensemble_short
#'
predict.ensemble_short = function(object, ...) {
  np = as.vector(object$fit_cv %*% object$nnls_weights)
  return(np)
}


#' Extraction of smoother weights for short-stacked ensemble learner
#'
#' @description
#' Extract smoother weights from short-stacked ensemble model created by
#' \code{\link{ensemble_short}}.
#'
#' @param object Short-stacked ensemble learner from \code{\link{ensemble_short}}.
#' @param ml List of ML models by \code{\link{create_method}} that was used as
#' input for \code{\link{ensemble_short}}.
#' @param x Covariate matrix of training sample.
#' @param y Vector of outcomes of training sample.
#' @param cf_mat Logical matrix with k columns of indicators representing the different folds
#' (for example created by \code{\link{prep_cf_mat}}).
#' @param quiet If FALSE, method that is currently computed is printed into the
#' console.
#' @param ... Ignore unused arguments
#'
#' @return Matrix of dimension \code{nrow(x)} x \code{nrow(x)} containing
#' ensemble smoother weights.
#'
#' @export
#'
#' @method weights ensemble_short
#'
weights.ensemble_short = function(object,
                                  ml,
                                  x, y,
                                  cf_mat,
                                  quiet = TRUE,
                                  ...) {

  if(is.null(object$ml)) {
    stop("Ensemble models were not saved after training.")
  } else if (is.character(object$ml) & length(object$ml) == 1) {
    ml_fit = readRDS(object$ml)
  } else {
    ml_fit = object$ml
  }

  smoother_weights = matrix(0, nrow = nrow(x), ncol = nrow(x))

  for (i in 1:ncol(cf_mat)) {

    fold = cf_mat[, i]
    x_tr = x[!fold, ]
    y_tr = y[!fold]
    x_te = x[fold, ]

    w_array = weights.ensemble_core(object = ml_fit[[i]], ml = ml, x_tr = x_tr, y_tr = y_tr, x_te = x_te, quiet)
    smoother_weights[fold, !fold] = agg_array(w_array, object$nnls_weights)

  }

  return(smoother_weights)

}
