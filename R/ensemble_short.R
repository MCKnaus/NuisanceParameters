#' Ensemble learner training for short stacking
#'
#' @description
#' \code{\link{ensemble_short}} trains an ensemble learner on the given training data
#' in the spirit of a short-stacking procedure.
#'
#' @param method List of methods built via \code{\link{create_method}} to be used in
#'               ensemble model.
#' @param X Covariate matrix of training sample.
#' @param Y Vector of outcomes of training sample.
#' @param subset Logical vector indicating which observations to use for determining
#'               ensemble weights. If not provided, all observations are used.
#' @param cf_mat Logical matrix with k columns of indicators representing the different folds.
#' @param storeModels Specifies whether to store the models. Value has to be an item
#'                    from \code{c("No", "Memory")} with "No" being the default value.
#' @param quiet If FALSE, progress output is printed into console.
#' @param pb A progress bar object to track overall computation progress.
#' @param pb_np String indicating the current nuisance parameter being
#'              estimated/predicted (for progress bar updates).
#'
#' @return List object containing:
#' \item{fit_cv}{matrix of dimension \code{nrow(X)} X \code{length(method)} containing
#' the cross-fitted predictions of the machine learning methods from the ensemble}
#' \item{method}{list of fitted machine learning models from all cross-fitting folds
#' (if \code{storeModels} is not set to "No")}
#'
#' @keywords internal
#'
ensemble_short <- function(method,
                           X, Y,
                           cf_mat,
                           subset = NULL,
                           storeModels = c("No", "Memory"),
                           quiet = TRUE, pb = NULL, pb_np = NULL) {
  ## Checks
  if (is.null(subset)) subset <- rep(TRUE, nrow(X))
  storeModels <- match.arg(storeModels)
  saveModels <- (storeModels != "No")

  # Store models and cross-validated predictions
  s <- if (saveModels) vector("list", length = ncol(cf_mat)) else NULL
  fit_cv <- make_fit_cv(method = method, N = nrow(X), Y = Y)

  # Learn in a subset D=d (or Z=z) and predict into test data
  for (i in 1:ncol(cf_mat)) {
    fold <- cf_mat[, i]
    X_tr <- X[!fold & subset, ]
    Y_tr <- Y[!fold & subset]
    X_te <- X[fold, ]

    # On-the-fold hyperparameter tuning
    method_tuned <- tune_learners(type = "tune_on_fold", method = method, X = X_tr, Y = Y_tr)

    # Core functionality
    method_fit <- ensemble_core(method_tuned, X_tr, Y_tr, quiet = quiet, pb = pb, pb_cf = i, pb_cv = ".", pb_np = pb_np)
    preds <- predict.ensemble_core(method_fit, method_tuned, X_tr, Y_tr, X_te, quiet = quiet, pb = pb, pb_cf = i, pb_cv = ".", pb_np = pb_np)
    
    if (saveModels) s[[i]] <- method_fit
    if (length(dim(fit_cv)) == 3) { fit_cv[fold, , ] <- preds
    } else { fit_cv[fold, ] <- preds }
  }

  output <- list(
    "fit_cv" = fit_cv,
    "method" = s
  )

  class(output) <- "ensemble_short"
  return(output)
}


#' Prediction of fitted values for short-stacked ensemble learner
#'
#' @description
#' Returns prediction of fitted values for a fitted and short-stacked ensemble
#' learner created by \code{\link{ensemble_short}}.
#'
#' @param object Short-stacked ensemble learner from \code{\link{ensemble_short}}
#' @param w Ensemble weights to aggregate predictions from different learners (optional).
#'          Needs to be a vector of length \code{ncol(object$fit_cv)}.
#' @param ... Additional arguments.
#'
#' @return Vector containing ensemble predictions.
#'
#' @keywords internal
#'
#' @method predict ensemble_short
#'
predict.ensemble_short <- function(object, w = NULL, ...) {
  if (is.null(w)) w <- rep(1 / ncol(object$fit_cv), ncol(object$fit_cv))

  np <- as.vector(object$fit_cv %*% w)

  return(np)
}
