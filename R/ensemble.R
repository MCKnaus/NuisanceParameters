#' Ensemble learner training
#'
#' @description
#' \code{\link{ensemble}} trains an ensemble learner on the given training data.
#'
#' @param method List of methods built via \code{\link{create_method}} to be used in
#' ensemble model
#' @param X Covariate matrix of training sample
#' @param Y Vector of outcomes of training sample
#' @param nfolds Number of folds used in cross-validation of ensemble weights
#' (default \code{nfolds = 5})
#' @param cv_mat Optional list of logical matrices of indicators representing the different folds,
#'              possibly from already estimated \code{NuisanceParameters} object.
#' @param quiet If FALSE, method that is currently computed is printed into the
#' console
#' @param pb A progress bar object to track overall computation progress.
#' @param pb_np String indicating the current nuisance parameter being
#'   estimated/predicted (for progress bar updates).
#' @param pb_cf Integer indicating the current cross-fitting fold
#'   (for progress bar updates).
#' @param pb_cv Integer indicating the current cross-validation fold
#'   (for progress bar updates).
#'
#' @return List object containing:
#' \item{fit_cv}{matrix of dimension \code{nrow(X)} X \code{length(method)} containing
#' the cross-fitted predictions of the machine learning methods from the ensemble}
#' \item{nnls_weights}{the weights that each machine learning method receives in
#' the ensemble}
#' \item{method}{list of fitted machine learning models from all cross-fitting folds
#' (if storeModels is not set to "No")}
#' 
#'
#' @export
#'
ensemble = function(method,
                    X, Y,
                    nfolds = 5, cv_mat = NULL,
                    quiet = TRUE, 
                    pb = NULL, pb_np = NULL, pb_cf = NULL, pb_cv = NULL
                    ) {

  # matrix to store the cross-validated predictions
  fit_cv = make_fit_cv(method = method, N = nrow(X), Y = Y)

  # cross-validation matrix
  if (is.null(cv_mat)) {cv_mat = prep_cf_mat(length(Y), nfolds)}


  # -> multiple methods specified - cross-validation of ensemble weights

  if (length(method) > 1) {
    # loop over different folds for cross-validation of weights
    for (i in 1:ncol(cv_mat)) {

      fold = cv_mat[, i]
      X_tr = X[!fold, ]
      Y_tr = Y[!fold]
      X_te = X[fold, ]

      method_fit = ensemble_core(method, X_tr, Y_tr, quiet = TRUE, pb = pb, pb_cf = pb_cf, pb_cv = i, pb_np = pb_np)
      preds = predict.ensemble_core(method_fit, method, X_tr, Y_tr, X_te, quiet = TRUE, pb = pb, pb_cf = pb_cf, pb_cv = i, pb_np = pb_np)
      if (length(dim(fit_cv)) == 3) { fit_cv[fold, , ] = preds } else { fit_cv[fold, ] = preds }

    }
    
    ## ensemble weights estimation
    # check: is multivalued propensity score being estimated?
    is_multinomial = !is.null(method[[1]]$multinomial)

    if (is_multinomial) {
      # stack_result <- stack_propensity_scores(preds = fit_cv, treatment = Y)
      # nnls_w <- stack_result$weights
      
      Y_stack <- as.vector(one_hot(Y)) # stack outcomes into a K*N vector
      X_stack <- sapply(dimnames(fit_cv)[[3]], function(m) as.vector(fit_cv[, , m])) # flatten predictions (K*N × M matrix)
      
      nnls_w = nnls_weights(X = X_stack, Y = Y_stack)
      
    } else {
      nnls_w = nnls_weights(X = fit_cv, Y = Y)
    }
    
    # run all methods on the full sample (fs)
    method_fit_full = ensemble_core(method, X, Y, quiet = quiet, pb = pb, pb_cf = pb_cf, pb_cv = ".", pb_np = pb_np)
  }


  # -> single method specified - no weights needed

  else if (length(method) == 1) {
    method_fit_full = ensemble_core(method, X, Y, quiet = quiet, pb = pb, pb_cf = pb_cf, pb_cv = ".", pb_np = pb_np)
    nnls_w <- stats::setNames(1, names(method))
  }
  

  output = list(
    "fit_cv" = fit_cv,
    "cv_mat" = cv_mat,
    "nnls_weights" = nnls_w,
    "method_fit"= method_fit_full
  )

  class(output) = "ensemble"

  return(output)

}


#' Prediction of fitted values by cross-validated ensemble learner
#'
#' @description
#' \code{\link{predict.ensemble}} makes prediction of fitted values for a cross-
#' validated ensemble learner.
#'
#' @param object Trained ensemble object from \code{\link{ensemble}}
#' @param method List of raw methods built via \code{\link{create_method}}
#' @param X Covariate matrix of training sample
#' @param Y Vector of outcomes of training sample
#' @param Xnew Covariate matrix of test sample
#' @param quiet If FALSE, method that is currently computed is printed into the
#' console.
#' @param pb A progress bar object to track overall computation progress.
#' @param pb_np String indicating the current nuisance parameter being
#'   estimated/predicted (for progress bar updates).
#' @param pb_cf Integer indicating the current cross-fitting fold
#'   (for progress bar updates).
#' @param ... Additional arguments.
#'
#' @return Vector containing ensemble predictions.
#'
#' @export
#'
#' @method predict ensemble
#'
predict.ensemble = function(object,
                            method,
                            X, Y, Xnew,
                            quiet = TRUE, pb = NULL, pb_cf = NULL, pb_np = NULL,
                            ...) {

  if (length(object$method) > 1) {

    pred = predict.ensemble_core(object = object$method_fit, method = method, X_tr = X, Y_tr = Y, X_te = Xnew, quiet = quiet, pb = pb, pb_cf = pb_cf, pb_cv = ".", pb_np = pb_np)
    np = as.vector(pred %*% object$nnls_weights)

  } else if (length(object$method) == 1) {

    pred = predict.ensemble_core(object = object$method_fit, method = method, X_tr = X, Y_tr = Y, X_te = Xnew, quiet = quiet, pb = pb, pb_cf = pb_cf, pb_cv = ".", pb_np = pb_np)
    np = as.vector(pred)

  }

  return(list("fit_cv" = pred, "np" = np))

}


#' Ensemble core fitting
#'
#' @description
#' Core function of \code{\link{ensemble}} used for fitting (i.e. training) the
#' inputted (machine learning) methods for the ensemble model for one given set
#' of covariates and targets.
#'
#' @param method List of methods built via \code{\link{create_method}} to be used in
#' ensemble model
#' @param X_tr Covariate matrix of training sample
#' @param Y_tr Vector of outcomes of training sample
#' @param quiet If FALSE, method that is currently computed is printed into the
#' console
#' @param pb A progress bar object to track overall computation progress.
#' @param pb_np String indicating the current nuisance parameter being
#'   estimated/predicted (for progress bar updates).
#' @param pb_cf Integer indicating the current cross-fitting fold
#'   (for progress bar updates).
#' @param pb_cv Integer indicating the current cross-validation fold
#'   (for progress bar updates).
#'
#' @return Returns list object containing the trained models of the ensemble.
#'
#' @export
#'
ensemble_core = function(method,
                         X_tr, Y_tr,
                         quiet = TRUE, 
                         pb = NULL, pb_cf = NULL, pb_cv = NULL, pb_np = NULL
                         ) {

  method_fit = list()

  # Loop over methods
  for (i in seq_along(method)) {
    update_progress(pb = pb, pb_np = pb_np, pb_cf = pb_cf, pb_cv = pb_cv, task = "Fitting", method = method[[i]]$method)

    wrapper  <- paste0(method[[i]]$method, "_fit")
    m_name   <- names(method)[i]
    X_tr_sel <- if (!is.null(method[[i]]$x_select)) X_tr[, method[[i]]$x_select] else X_tr
    
    # Default case 
    if (is.null(method[[i]]$multinomial)) {
      fit <- do.call(wrapper, list(X = X_tr_sel, Y = Y_tr, arguments = method[[i]]$arguments))
      } else {
        # Multiclass treatment case
        fit <- switch(method[[i]]$multinomial,
                      "one-vs-one"  = { do.call(ovo_fit, list(X = X_tr_sel, Y = Y_tr, method = method[[i]]$method)) },
                      "one-vs-rest" = { do.call(ovr_fit, list(X = X_tr_sel, Y = Y_tr, method = method[[i]]$method)) },
                      "multiclass"  = { do.call(wrapper, list(X = X_tr_sel, Y = Y_tr)) }
        )
        }

    method_fit[[m_name]] <- fit
    
    }

  class(method_fit) <- "ensemble_core"
  
  return(method_fit)
}


#' Prediction of fitted values for ensemble
#'
#' @description
#' Predicts fitted values for some (new) covariate matrix based on the fully
#' trained models of the ensemble model from \code{\link{ensemble_core}}.
#'
#' @param object List of fitted methods (from \code{\link{ensemble_core}}
#' @param method List of raw methods built via \code{\link{create_method}}
#' @param X_tr Covariate matrix of training sample
#' @param Y_tr Vector of outcomes of training sample
#' @param X_te Covariate matrix of test sample
#' @param quiet If FALSE, method that is currently computed is printed into the
#' console
#' @param pb A progress bar object to track overall computation progress.
#' @param pb_np String indicating the current nuisance parameter being
#'   estimated/predicted (for progress bar updates).
#' @param pb_cf Integer indicating the current cross-fitting fold
#'   (for progress bar updates).
#' @param pb_cv Integer indicating the current cross-validation fold
#'   (for progress bar updates).
#' @param ... Additional arguments.
#'
#' @return Returns matrix of dimension \code{nrow(X_te)} X \code{length(method)}
#' that contains the predictions of fitted values for each method in the ensemble.
#'
#' @export
#'
#' @method predict ensemble_core
#'
predict.ensemble_core = function(object, method,
                                 X_tr, Y_tr, X_te,
                                 quiet = TRUE, pb = NULL, pb_cf = NULL, pb_cv = NULL, pb_np = NULL, 
                                 ...) {
  
  fit_mat = NULL
  
  for (i in 1:length(object)) {
    update_progress(pb = pb, pb_np = pb_np, pb_cf = pb_cf, pb_cv = pb_cv, task = "Predict", method = method[[i]]$method)
    
    wrapper  <- paste0(method[[i]]$method, "_fit")
    m_name   <- names(method)[i]
    X_tr_sel <- if (!is.null(method[[i]]$x_select)) X_tr[, method[[i]]$x_select] else X_tr
    X_te_sel <- if (!is.null(method[[i]]$x_select)) X_te[, method[[i]]$x_select] else X_te
    
    # Default case 
    if (is.null(method[[i]]$multinomial)) {
      pred <- do.call(paste0("predict.", wrapper), list(object[[i]], X = X_tr_sel, Y = Y_tr, Xnew = X_te_sel))
    } else {
      # Multiclass treatment case
      pred <- switch(method[[i]]$multinomial,
                     "one-vs-one"  = { do.call(predict.ovo_fit, list(object[[i]], X = X_tr_sel, Y = Y_tr, Xnew = X_te_sel, method = method[[i]]$method)) },
                     "one-vs-rest" = { do.call(predict.ovr_fit, list(object[[i]], X = X_tr_sel, Y = Y_tr, Xnew = X_te_sel, method = method[[i]]$method)) },
                     "multiclass"  = { do.call(paste0("predict.", wrapper), list(object[[i]], X = X_tr_sel, Y = Y_tr, Xnew = X_te_sel)) }
      )
      }
    
    # If prediction is vector-like or multi-column
    is_vec <- is.vector(pred) || ncol(as.matrix(pred)) == 1
    
    # Initialize
    if (is.null(fit_mat)) {
      if (is_vec) {
        fit_mat <- matrix(NA, nrow = nrow(X_te), ncol = length(object))
      } else {
        fit_mat <- array(NA, dim = c(nrow(X_te), ncol(pred), length(object)))
      }
    }
    
    # Assign
    if (is_vec) {
      fit_mat[, i] <- pred
    } else {
      fit_mat[, , i] <- as.matrix(pred)
    }
    }
    
    return(fit_mat)
    
  }