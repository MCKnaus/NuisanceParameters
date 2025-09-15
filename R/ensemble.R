#' Train ensemble learner with short stacking
#'
#' Train an ensemble learner on the given training data using a short-stacking procedure.
#'
#' @param methods List of methods created with \code{\link{create_method}} to be 
#'  used in the ensemble model.
#' @param X Covariate matrix of the training sample.
#' @param Y Vector of outcomes of the training sample.
#' @param subset Logical vector indicating which observations to use for 
#'  determining ensemble weights. If not provided, all observations are used.
#' @param cf_mat Logical matrix with \code{k} columns of indicators representing 
#'  the different folds.
#' @param storeModels Character string specifying whether to store the models. 
#'  Must be one of \code{c("No", "Memory", "Disk")} with "No" as the default.
#' @param quiet Logical. If \code{FALSE}, progress output is printed to the console.
#' @param pb A progress bar object to track overall computation progress.
#' @param pb_np String indicating the current nuisance parameter being 
#'  estimated or predicted (for progress bar updates).
#'
#' @return A list with class \code{"ensemble_short"} containing:
#' \describe{
#'   \item{cf_preds}{A matrix of dimension \code{nrow(X)} x \code{length(methods)} 
#'    containing the cross-fitted predictions of the machine learning methods from the ensemble.}
#'   \item{ens_models}{A list of fitted machine learning models from all 
#'    cross-fitting folds (if \code{storeModels} is not set to "No").}
#' }
#'
#' @keywords internal
ensemble_short <- function(methods,
                           X, Y,
                           cf_mat,
                           subset = NULL,
                           storeModels = c("No", "Memory", "Disk"),
                           quiet = TRUE, pb = NULL, pb_np = NULL) {
  # Checks
  if (is.null(subset)) subset <- rep(TRUE, nrow(X))
  storeModels <- match.arg(storeModels)

  cf_preds <- make_cf_preds(methods = methods, N = nrow(X), Y = Y)
  ens_models <- if (storeModels != "No") vector("list", length = ncol(cf_mat)) else NULL

  # Learn in a subset D=d (Z=z); predict into test data
  for (i in 1:ncol(cf_mat)) {
    fold <- cf_mat[, i]
    X_tr <- X[!fold & subset, ]
    Y_tr <- Y[!fold & subset]
    X_te <- X[fold, ]
    
    # On-the-fold hyperparameter tuning
    methods_tuned <- tune_learners(type = "tune_on_fold", methods = methods, X = X_tr, Y = Y_tr)
    
    fits <- ensemble_core(
      methods = methods_tuned, X_tr = X_tr, Y_tr = Y_tr, quiet = quiet, 
      pb = pb, pb_cf = i, pb_cv = ".", pb_np = pb_np
      )
    preds <- predict.ensemble_core(
      object = fits, methods = methods_tuned, X_tr = X_tr, Y_tr = Y_tr, X_te = X_te, 
      quiet = quiet, pb = pb, pb_cf = i, pb_cv = ".", pb_np = pb_np
      )
    
    if (storeModels != "No") ens_models[[i]] <- fits
    if (length(dim(cf_preds)) == 3) { cf_preds[fold, , ] <- preds } else { cf_preds[fold, ] <- preds }
  }
  output <- list(
    "cf_preds" = cf_preds,
    "ens_models" = ens_models
  )
  class(output) <- "ensemble_short"
  return(output)
}


#' Train ensemble learner with standard stacking
#'
#' Train an ensemble learner on the given training data using cross-validation for ensemble weights.
#'
#' @param methods List of methods created with \code{\link{create_method}} to be 
#'  used in the ensemble model.
#' @param X Covariate matrix of the training sample.
#' @param Y Vector of outcomes of the training sample.
#' @param nfolds Integer. Number of folds used in cross-validation of 
#'  ensemble weights. Default is 5.
#' @param quiet Logical. If \code{FALSE}, progress output is printed to the console.
#' @param pb A progress bar object to track overall computation progress.
#' @param pb_np String indicating the current nuisance parameter being 
#'  estimated or predicted (for progress bar updates).
#' @param pb_cf Integer indicating the current cross-fitting fold 
#'  (for progress bar updates).
#' @param pb_cv Integer indicating the current cross-validation fold 
#'  (for progress bar updates).
#'
#' @return A list with class \code{"ensemble"} containing:
#' \describe{
#'   \item{cf_preds}{A matrix of dimension \code{nrow(X)} x \code{length(methods)} 
#'    containing the cross-fitted predictions of the machine learning methods from the ensemble.}
#'   \item{nnls_w}{A numeric vector of weights that each machine learning 
#'    method receives in the ensemble.}
#'   \item{ens_models}{A list of fitted machine learning models from the full sample.}
#' }
#'
#' @keywords internal
ensemble <- function(methods,
                     X, Y,
                     nfolds = 5,
                     quiet = TRUE,
                     pb = NULL, pb_np = NULL, pb_cf = NULL, pb_cv = NULL) {
  
  cv_mat <- prep_cf_mat(length(Y), nfolds)
  cf_preds <- make_cf_preds(methods = methods, N = nrow(X), Y = Y)

  ## Multiple methods specified: cross-validation of ensemble weights
  if (length(methods) > 1) {
    for (i in 1:ncol(cv_mat)) {
      fold <- cv_mat[, i]
      X_tr <- X[!fold, ]
      Y_tr <- Y[!fold]
      X_te <- X[fold, ]

      fits <- ensemble_core(
        methods = methods, X_tr = X_tr, Y_tr = Y_tr, quiet = quiet, 
        pb = pb, pb_cf = pb_cf, pb_cv = i, pb_np = pb_np
        )
      preds <- predict.ensemble_core(
        object = fits, methods = methods, X_tr = X_tr, Y_tr = Y_tr, X_te = X_te, 
        quiet = quiet, pb = pb, pb_cf = pb_cf, pb_cv = i, pb_np = pb_np
        )
      
      if (length(dim(cf_preds)) == 3) { cf_preds[fold, , ] <- preds } else { cf_preds[fold, ] <- preds }
    }

    # Is multivalued propensity score being estimated?
    is_mult <- !is.null(methods[[1]]$multinomial)
    
    nnls_w <- nnls_weights(X = cf_preds, Y = Y, is_mult = is_mult, do_bfgs = FALSE)

    # Re-run all methods on the full sample (fs)
    fits_fs <- ensemble_core(
      methods = methods, X_tr = X, Y_tr = Y, quiet = quiet, 
      pb = pb, pb_cf = pb_cf, pb_cv = ".", pb_np = pb_np
      )
  }

  ## Single method specified: no weights needed
  else if (length(methods) == 1) {
    fits_fs <- ensemble_core(
      methods = methods, X_tr = X, Y_tr = Y, quiet = quiet, 
      pb = pb, pb_cf = pb_cf, pb_cv = ".", pb_np = pb_np
      )
    nnls_w <- stats::setNames(1, names(methods))
  }
  # "cf_preds" left unsaved due to standard stacking
  output <- list(
    "nnls_w" = nnls_w,
    "ens_models" = fits_fs
  )
  
  class(output) <- "ensemble"
  return(output)
}


#' Predict method for short-stacked ensemble
#'
#' Generate predictions from a fitted short-stacked ensemble learner created 
#' with \code{\link{ensemble_short}}.
#'
#' @param object An object of class \code{"ensemble_short"} returned by 
#'  \code{\link{ensemble_short}}.
#' @param w Optional numeric vector of ensemble weights used to aggregate 
#'  predictions from different learners. Must have length \code{ncol(object$cf_preds)}.
#' @param ... Ignored additional arguments.
#'
#' @return A numeric vector containing ensemble predictions.
#'
#' @keywords internal
#' @method predict ensemble_short
predict.ensemble_short <- function(object, w = NULL, ...) {
  if (is.null(w)) w <- rep(1 / ncol(object$cf_preds), ncol(object$cf_preds))
  np <- as.vector(object$cf_preds %*% w)
  
  return(np)
}


#' Predict method for standard-stacked ensemble
#'
#' Generate predictions from a trained ensemble standard-stacked learner 
#' created with \code{\link{ensemble}}.
#'
#' @param object An object of class \code{"ensemble"} returned by \code{\link{ensemble}}.
#' @param methods List of methods created with \code{\link{create_method}}.
#' @param X Covariate matrix of the training sample.
#' @param Y Vector of outcomes of the training sample.
#' @param Xnew Covariate matrix of the test sample.
#' @param quiet Logical. If \code{FALSE}, progress output is printed to the console.
#' @param pb A progress bar object to track overall computation progress.
#' @param pb_np String indicating the current nuisance parameter being 
#'  estimated or predicted (for progress bar updates).
#' @param pb_cf Integer indicating the current cross-fitting fold 
#'  (for progress bar updates).
#' @param ... Ignored additional arguments.
#'
#' @return A list containing:
#' \describe{
#'   \item{cf_preds}{A matrix of predictions for the test sample from the individual learners.}
#'   \item{np}{A numeric vector containing aggregated ensemble predictions.}
#' }
#'
#' @keywords internal
#' @method predict ensemble
predict.ensemble <- function(object,
                             methods,
                             X, Y, Xnew,
                             quiet = TRUE, pb = NULL, pb_cf = NULL, pb_np = NULL,
                             ...) {
  
  cf_preds <- predict.ensemble_core(
    object = object$ens_models, methods = methods, X_tr = X, Y_tr = Y, X_te = Xnew, 
    quiet = quiet, pb = pb, pb_cf = pb_cf, pb_cv = ".", pb_np = pb_np
    )
  
  # Is multiclass propensity score being estimated?
  is_mult <- !is.null(methods[[1]]$multinomial)
  
  if (is_mult) {
    np <- agg_array(a = cf_preds, w = object$nnls_w)
  } else {
    np <- as.vector(cf_preds %*% (if (length(object$ens_models) > 1) object$nnls_w else 1))
    }
  
  return(list("cf_preds" = cf_preds, "np" = np))
}


#' Ensemble core fitting
#'
#' Train the inputted machine learning methods for one given set of covariates 
#' and outcomes within the ensemble framework.
#'
#' @param methods List of methods created with \code{\link{create_method}} to be 
#'  used in the ensemble model.
#' @param X_tr Covariate matrix of the training sample.
#' @param Y_tr Vector of outcomes of the training sample.
#' @param quiet Logical. If \code{FALSE}, progress output is printed to the console.
#' @param pb A progress bar object to track overall computation progress.
#' @param pb_np String indicating the current nuisance parameter being estimated 
#'  or predicted (for progress bar updates).
#' @param pb_cf Integer indicating the current cross-fitting fold 
#'  (for progress bar updates).
#' @param pb_cv Integer indicating the current cross-validation fold 
#'  (for progress bar updates).
#'
#' @return An object of class \code{"ensemble_core"} containing the 
#'  trained models of the ensemble.
#'
#' @keywords internal
ensemble_core <- function(methods,
                          X_tr, Y_tr,
                          quiet = TRUE,
                          pb = NULL, pb_cf = NULL, pb_cv = NULL, pb_np = NULL) {
  fits <- list()
  
  # Loop over methods
  for (i in seq_along(methods)) {
    update_progress(pb = pb, pb_np = pb_np, pb_cf = pb_cf, pb_cv = pb_cv, task = "Fit", method = methods[[i]]$method)

    wrapper <- paste0(methods[[i]]$method, "_fit")
    m_name <- names(methods)[i]
    X_tr_sel <- if (!is.null(methods[[i]]$x_select)) X_tr[, methods[[i]]$x_select] else X_tr

    # Default case
    if (is.null(methods[[i]]$multinomial)) {
      m_fit <- do.call(wrapper, list(X = X_tr_sel, Y = Y_tr, arguments = methods[[i]]$arguments))
    } else {
      # Multiclass treatment case
      m_fit <- switch(methods[[i]]$multinomial,
        "one-vs-one" = {
          do.call(ovo_fit, list(X = X_tr_sel, Y = Y_tr, method = methods[[i]]$method, parallel = methods[[i]]$parallel))
        },
        "one-vs-rest" = {
          do.call(ovr_fit, list(X = X_tr_sel, Y = Y_tr, method = methods[[i]]$method))
        },
        "multiclass" = {
          do.call(wrapper, list(X = X_tr_sel, Y = Y_tr))
        }
      )
    }
    fits[[m_name]] <- m_fit
  }

  class(fits) <- "ensemble_core"
  return(fits)
}


#' Prediction of fitted values for ensemble
#'
#' Generate predictions for (new) covariates based on trained ensemble models 
#' from \code{\link{ensemble_core}}.
#'
#' @param object An object of class \code{"ensemble_core"} returned by \code{\link{ensemble_core}}.
#' @param methods List of raw methods created with \code{\link{create_method}}.
#' @param X_tr Covariate matrix of the training sample.
#' @param Y_tr Vector of outcomes of the training sample.
#' @param X_te Covariate matrix of the test sample.
#' @param quiet Logical. If \code{FALSE}, progress output is printed to the console.
#' @param pb A progress bar object to track overall computation progress.
#' @param pb_np String indicating the current nuisance parameter being 
#'  estimated or predicted (for progress bar updates).
#' @param pb_cf Integer indicating the current cross-fitting fold 
#'  (for progress bar updates).
#' @param pb_cv Integer indicating the current cross-validation fold 
#'  (for progress bar updates).
#' @param ... Ignored additional arguments.
#'
#' @return A matrix of dimension \code{nrow(X_te)} x \code{length(methods)} containing 
#' predictions from each method in the ensemble. If learners return multi-column predictions, 
#' an array of dimension \code{nrow(X_te)} x \code{ncol(m_pred)} x \code{length(methods)} is returned.
#'
#' @keywords internal
#' @method predict ensemble_core
predict.ensemble_core <- function(object, methods,
                                  X_tr, Y_tr, X_te,
                                  quiet = TRUE, pb = NULL, pb_cf = NULL, pb_cv = NULL, pb_np = NULL,
                                  ...) {
  preds <- NULL

  for (i in 1:length(object)) {
    update_progress(pb = pb, pb_np = pb_np, pb_cf = pb_cf, pb_cv = pb_cv, task = "Pred", method = methods[[i]]$method)

    wrapper <- paste0(methods[[i]]$method, "_fit")
    X_tr_sel <- if (!is.null(methods[[i]]$x_select)) X_tr[, methods[[i]]$x_select] else X_tr
    X_te_sel <- if (!is.null(methods[[i]]$x_select)) X_te[, methods[[i]]$x_select] else X_te

    # Default case
    if (is.null(methods[[i]]$multinomial)) {
      m_pred <- do.call(paste0("predict.", wrapper), list(object[[i]], X = X_tr_sel, Y = Y_tr, Xnew = X_te_sel))
    } else {

      # Multiclass treatment case
      m_pred <- switch(methods[[i]]$multinomial,
        "one-vs-one" = {
          do.call(predict.ovo_fit, list(object[[i]], X = X_tr_sel, Y = Y_tr, Xnew = X_te_sel, 
                                        method = methods[[i]]$method, parallel = methods[[i]]$parallel))
        },
        "one-vs-rest" = {
          do.call(predict.ovr_fit, list(object[[i]], X = X_tr_sel, Y = Y_tr, 
                                        Xnew = X_te_sel, method = methods[[i]]$method))
        },
        "multiclass" = {
          do.call(paste0("predict.", wrapper), list(object[[i]], X = X_tr_sel, Y = Y_tr, Xnew = X_te_sel))
        }
      )
    }

    # If prediction is vector-like or multi-column
    is_vec <- is.vector(m_pred) || ncol(as.matrix(m_pred)) == 1

    # Initialize
    if (is.null(preds)) {
      if (is_vec) {
        preds <- matrix(NA, nrow = nrow(X_te), ncol = length(object))
      } else {
        preds <- array(NA, dim = c(nrow(X_te), ncol(m_pred), length(object) + 1))
      }
    }
    
    # Assign
    if (is_vec) {
      preds[, i] <- m_pred
    } else {
      preds[, , i] <- as.matrix(m_pred)
    }
  }
  
  # Drop the dummy slice (robust when there is a single learner)
  if (!is_vec) { preds <- preds[, , 1:length(object), drop = FALSE] }

  return(preds)
}
