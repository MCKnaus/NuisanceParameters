#' Ensemble learner training for short stacking
#'
#' @description
#' \code{\link{ensemble_short}} trains an ensemble learner on the given training data
#' in the spirit of a short-stacking procedure.
#'
#' @param method List of methods built via \code{\link{create_method}} to be used in
#' ensemble model.
#' @param X Covariate matrix of training sample.
#' @param Y Vector of outcomes of training sample.
#' @param subset Logical vector indicating which observations to use for determining
#' ensemble weights. If not provided, all observations are used.
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
#' @param pb A progress bar object to track overall computation progress.
#' @param pb_np String indicating the current nuisance parameter being
#'   estimated/predicted (for progress bar updates).
#'
#' @return List object containing:
#' \item{fit_cv}{matrix of dimension \code{nrow(X)} X \code{length(method)} containing
#' the cross-fitted predictions of the machine learning methods from the ensemble}
#' \item{method}{list of fitted machine learning models from all cross-fitting folds
#' (if storeModels is not set to "No")}
#'
#' @export
#'
ensemble_short = function(method,
                          X, Y,
                          cf_mat,
                          subset = NULL,
                          storeModels = c("No", "Memory", "Disk"),
                          path = NULL,
                          quiet = TRUE, pb = NULL, pb_np = NULL
                          ) {


  ## Checks
  if(is.null(subset)) subset = rep(TRUE, nrow(X))
  storeModels = match.arg(storeModels)
  saveModels = (storeModels != "No")

  s = NULL
  if(saveModels) s = vector("list", length = ncol(cf_mat))

  # Matrix to store the cross-validated predictions
  fit_cv = make_fit_cv(method = method, N = nrow(X), Y = Y)

  for (i in 1:ncol(cf_mat)) {

    fold = cf_mat[, i]
    X_tr = X[!fold & subset, ] # learn in a subset D=d (or Z=z)
    Y_tr = Y[!fold & subset]
    X_te = X[fold, ]           # predict into test data
    
    ## Core functionality
    method_fit = ensemble_core(method, X_tr, Y_tr, quiet = quiet, pb = pb, pb_cf = i, pb_cv = ".", pb_np = pb_np)
    preds = predict.ensemble_core(method_fit, method, X_tr, Y_tr, X_te, quiet = quiet, pb = pb, pb_cf = i, pb_cv = ".", pb_np = pb_np)
    
    if (length(dim(fit_cv)) == 3) { fit_cv[fold, , ] = preds } else { fit_cv[fold, ] = preds }
    if(saveModels) s[[i]] = method_fit

  }

  output = list(
    "fit_cv" = fit_cv,
    "method" = s
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
#' @param object Short-stacked ensemble learner from \code{\link{ensemble_short}}
#' @param w Ensemble weights to aggregate predictions from different learners (optional).
#' Needs to be a vector of length \code{ncol(object$fit_cv)}.
#' @param ... Ignore unused arguments
#'
#' @return Vector containing ensemble predictions.
#'
#' @export
#'
#' @method predict ensemble_short
#'
predict.ensemble_short = function(object, w = NULL, ...) {

  if (is.null(w)) w = rep(1 / ncol(object$fit_cv), ncol(object$fit_cv))

  np = as.vector(object$fit_cv %*% w)

  return(np)

}