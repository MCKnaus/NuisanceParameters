#' Ensemble learner training
#'
#' @description
#' \code{\link{ensemble}} trains an ensemble learner on the given training data.
#'
#' @param ml List of methods built via \code{\link{create_method}} to be used in
#' ensemble model
#' @param x Covariate matrix of training sample
#' @param y Vector of outcomes of training sample
#' @param nfolds Number of folds used in cross-validation of ensemble weights
#' (default \code{nfolds = 5})
#' @param quiet If FALSE, method that is currently computed is printed into the
#' console
#'
#' @return List object containing:
#' \item{fit_cv}{matrix of dimension \code{nrow(x)} x \code{length(ml)} containing
#' the cross-fitted predictions of the machine learning methods from the ensemble}
#' \item{nnls_weights}{the weights that each machine learning method receives in
#' the ensemble}
#' \item{ml}{list of fitted machine learning models from all cross-fitting folds
#' (if storeModels is not set to "No")}
#'
#' @export
#'
ensemble = function(ml,
                    x, y,
                    nfolds = 5,
                    quiet = TRUE, 
                    pb = NULL, pb_np = NULL, pb_cf = NULL, pb_cv = NULL
                    ) {

  # matrix to store the cross-validated predictions
  fit_cv = make_fit_cv(ml = ml, n = nrow(x),  learner = "t")

  # cross-validation matrix
  cf_mat = prep_cf_mat(length(y), nfolds)


  ### multiple ml methods specified - cross-validation of ensemble weights ###

  if (length(ml) > 1) {
    # loop over different folds for cross-validation of weights
    for (i in 1:ncol(cf_mat)) {

      fold = cf_mat[, i]
      x_tr = x[!fold, ]
      y_tr = y[!fold]
      x_te = x[fold, ]

      ml_fit = ensemble_core(ml, x_tr, y_tr, quiet = TRUE, pb = pb, pb_cf = pb_cf, pb_cv = i, pb_np = pb_np)
      fit_cv[fold, ] = predict.ensemble_core(ml_fit, ml, x_tr, y_tr, x_te, quiet = TRUE, pb = pb, pb_cf = pb_cf, pb_cv = i, pb_np = pb_np)

    }

    # estimate ensemble weights
    nnls_w = nnls_weights(X = fit_cv, y = y)

    # run all methods on the full sample (fs)
    ml_fit_full = ensemble_core(ml, x, y, quiet = quiet, pb = pb, pb_cf = pb_cf, pb_cv = "fs", pb_np = pb_np)
  }


  ### only one ml method specified - no weights needed ###

  else if (length(ml) == 1) {
    ml_fit_full = ensemble_core(ml, x, y, quiet = quiet, pb = pb, pb_cf = pb_cf, pb_cv = ".", pb_np = pb_np)
    nnls_w = fit_cv = NULL
  }

  output = list(
    "fit_cv" = fit_cv,
    "nnls_weights" = nnls_w,
    "ml_fit"= ml_fit_full
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
#' @param ml List of raw methods built via \code{\link{create_method}}
#' @param x Covariate matrix of training sample
#' @param y Vector of outcomes of training sample
#' @param xnew Covariate matrix of test sample
#' @param quiet If FALSE, method that is currently computed is printed into the
#' console.
#' @param ... Ignore unused arguments
#'
#' @return Vector containing ensemble predictions.
#'
#' @export
#'
#' @method predict ensemble
#'
predict.ensemble = function(object,
                            ml,
                            x, y, xnew,
                            quiet = TRUE, pb = NULL, pb_cf = NULL, pb_np = NULL, 
                            ...) {

  if (length(object$ml) > 1) {

    pred = predict.ensemble_core(object = object$ml_fit, ml = ml, x_tr = x, y_tr = y, x_te = xnew, quiet = quiet, pb = pb, pb_cf = pb_cf, pb_cv = "op", pb_np = pb_np)
    np = as.vector(pred %*% object$nnls_weights)

  } else if (length(object$ml) == 1) {

    pred = predict.ensemble_core(object = object$ml_fit, ml = ml, x_tr = x, y_tr = y, x_te = xnew, quiet = quiet, pb = pb, pb_cf = pb_cf, pb_cv = "op", pb_np = pb_np)
    np = as.vector(pred)

  }

  return(list("fit_cv" = pred, "np" = np))

}


#' Smoother weights from cross-validated ensemble learner
#'
#' @description
#' \code{\link{weights.ensemble}} extract ensemble smoother weights for some
#' (new) covariate matrix based on a fitted and cross-validated ensemble model from
#' \code{\link{ensemble}}.
#'
#' @param object Trained ensemble object from \code{\link{ensemble}}
#' @param ml List of raw methods built via \code{\link{create_method}}
#' @param x Covariate matrix of training sample
#' @param y Vector of outcomes of training sample
#' @param xnew Covariate matrix of test sample
#' @param quiet If FALSE, method that is currently computed is printed into the
#' console
#' @param ... Ignore unused arguments
#'
#' @return Matrix of dimension \code{nrow(xnew)} x \code{nrow(x)} containing
#' ensemble smoother weights.
#'
#' @export
#'
#' @method weights ensemble
#'
weights.ensemble = function(object,
                            ml,
                            x, y, xnew,
                            quiet = TRUE,
                            ...) {

  if (length(object$ml) > 1) {

    w_array = weights.ensemble_core(object = object$ml_fit, ml = ml, x_tr = x, y_tr = y, x_te = xnew, quiet)
    smoother_weights = agg_array(w_array, object$nnls_weights)

  } else if (length(object$ml) == 1) {

    w_array = weights.ensemble_core(object = object$ml_fit, ml = ml, x_tr = x, y_tr = y, x_te = xnew, quiet)
    smoother_weights = array(w_array[, , 1], dim = dim(w_array)[-3])

  }

  return(smoother_weights)
}


#' Ensemble core fitting
#'
#' @description
#' Core function of \code{\link{ensemble}} used for fitting (i.e. training) the
#' inputted (machine learning) methods for the ensemble model for one given set
#' of covariates and targets.
#'
#' @param ml List of methods built via \code{\link{create_method}} to be used in
#' ensemble model
#' @param x_tr Covariate matrix of training sample
#' @param y_tr Vector of outcomes of training sample
#' @param quiet If FALSE, method that is currently computed is printed into the
#' console
#'
#' @return Returns list object containing the trained models of the ensemble.
#'
#' @export
#'
ensemble_core = function(ml,
                         x_tr, y_tr,
                         quiet = TRUE, pb = NULL, pb_cf = NULL, pb_cv = NULL, pb_np = NULL
                         ) {

  # initialize list object to be filled
  ml_fit = list()

  # loop over specified methods
  for (i in 1:length(ml)) {
    update_progress(pb = pb, pb_np = pb_np, pb_cf = pb_cf, pb_cv = pb_cv, task = "Fitting", method = ml[[i]]$method)

    wrapper = paste0(ml[[i]]$method, "_fit")

    # check whether subset of variables specified and run method accordingly
    if (is.null(ml[[i]]$x_select)) {
      fit = do.call(wrapper, list(x = x_tr, y = y_tr, arguments = ml[[i]]$arguments))
    } else {
      fit = do.call(wrapper, list(x = x_tr[, ml[[i]]$x_select], y = y_tr, arguments = ml[[i]]$arguments))
    }

    # add models to list object
    ml_fit[[wrapper]] = fit

  }

  class(ml_fit) = "ensemble_core"

  return(ml_fit)
}


#' Prediction of fitted values for ensemble
#'
#' @description
#' Predicts fitted values for some (new) covariate matrix based on the fully
#' trained models of the ensemble model from \code{\link{ensemble_core}}.
#'
#' @param object List of fitted methods (from \code{\link{ensemble_core}}
#' @param ml List of raw methods built via \code{\link{create_method}}
#' @param x_tr Covariate matrix of training sample
#' @param y_tr Vector of outcomes of training sample
#' @param x_te Covariate matrix of test sample
#' @param quiet If FALSE, method that is currently computed is printed into the
#' console
#' @param ... Ignore unused arguments
#'
#' @return Returns matrix of dimension \code{nrow(x_te)} x \code{length(ml)}
#' that contains the predictions of fitted values for each method in the ensemble.
#'
#' @export
#'
#' @method predict ensemble_core
#'
predict.ensemble_core = function(object, ml,
                                 x_tr, y_tr, x_te,
                                 quiet = TRUE, pb = NULL, pb_cf = NULL, pb_cv = NULL, pb_np = NULL,
                                 ...) {

  # initialize matrix to be filled
  fit_mat = matrix(NA, nrow(x_te), length(object))

  # loop over all trained models
  for (i in 1:length(object)) {
    update_progress(pb = pb, pb_np = pb_np, pb_cf = pb_cf, pb_cv = pb_cv, task = "Predict", method = ml[[i]]$method)

    # extract predictions
    if (is.null(ml[[i]]$x_select)) {
      temp = do.call(paste0("predict.", names(object)[i]), list(object[[i]], x = x_tr, y = y_tr, xnew = x_te))
    } else {
      temp = do.call(paste0("predict.", names(object)[i]), list(object[[i]], x = x_tr[, ml[[i]]$x_select], y = y_tr, xnew = x_te[, ml[[i]]$x_select]))
    }
    fit_mat[,i] = temp
  }

  return(fit_mat)

}


#' Ensemble smoother weights
#'
#' @description
#' Extract smoother weights for some (new) covariate matrix based on the fully
#' trained models of the ensemble model from \code{\link{ensemble_core}}.
#'
#' @param object List of fitted methods (from \code{\link{ensemble_core}}
#' @param ml List of raw methods built via \code{\link{create_method}}
#' @param x_tr Covariate matrix of training sample
#' @param y_tr Vector of outcomes of training sample
#' @param x_te Covariate matrix of test sample
#' @param quiet If FALSE, method that is currently computed is printed into the
#' console
#' @param ... Ignore unused arguments
#'
#' @return Returns a three-dimensional array of dimension
#' \code{nrow(xnew)} x \code{nrow(x)} x \code{length(ml)}
#' containing the smoother weights that deliver predictions of each method
#' where each row gives the weight that each training outcome received in the
#' prediction of the respective test set observation.
#'
#' @export
#'
#' @method weights ensemble_core
#'
weights.ensemble_core = function(object, ml,
                                 x_tr, y_tr, x_te,
                                 quiet = TRUE,
                                 ...) {

  # initialize 3D array to be filled
  w_array = array(data = 0, dim = c(nrow(x_te), nrow(x_tr), length(object)))

  # loop over all trained models
  for (i in 1:length(object)) {

    if (isFALSE(quiet)) print(paste0("Smoother weights: ", ml[[i]]$method))

    # extract smoother weights
    if (is.null(ml[[i]]$x_select)) {
      temp = do.call(paste0("weights.", names(object)[i]), list(object[[i]], x = x_tr, y = y_tr, xnew = x_te))
    } else {
      temp = do.call(paste0("weights.", names(object)[i]), list(object[[i]], x = x_tr[, ml[[i]]$x_select], y = y_tr, xnew = x_te[, ml[[i]]$x_select]))
    }

    w_array[, , i] = temp

  }

  return(w_array)

}


#' Method creation for ensemble
#'
#' @description
#' Creates the methods to be used in the subsequent \code{\link{ensemble}} model.
#'
#' @param method Choose method from \code{c("mean", "ols", "ridge", "plasso", "forest_grf", "lasso", "knn", "forest_drf")}.
#' To be continued.
#' @param x_select Optional logical vector of length equal to the number of
#' columns of the covariate matrix indicating which variables should be used by
#' this method. E.g. tree-based methods usually should not be provided with the
#' interactions that Lasso is using.
#' @param arguments Optional list containing the additional arguments that should be
#' passed to the underlying method.
#' @param name Optional string naming the method.
#'
#' @return List object that can be passed as input to \code{\link{ensemble}}
#'
#' @export
#'
#' @examples
#' # create list of ml methods for ensemble
#' ml = list(
#'  "ols" = create_method("ols"),
#'  "forest_grf" = create_method("forest_grf"),
#'  "knn" = create_method("knn", arguments = list("k" = 3))
#' )
#'
create_method = function(
    method = c("mean", "ols", "ridge", "plasso", "forest_grf", "lasso", "knn", "forest_drf"),
    x_select = NULL,
    arguments = list(),
    name = NULL) {

  # check if method is valid
  method = match.arg(method)
  # check if other inputs are valid
  if (!(is.null(arguments) | is.list(arguments))) stop("Provide either NULL or list for arguments.")
  if (!(is.null(x_select) | is.logical(x_select))) stop("Provide either NULL or logical for x_select.")
  if (!((is.character(name) & length(name) == 1) | is.null(name))) stop("Provide single string to name method.")

  list(method = method, arguments = arguments, x_select = x_select, name = name)
}
