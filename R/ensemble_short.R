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
#' @param subset Logical vector indicating which observations to use for determining
#' ensemble weights. If not provided, all observations are used.
#' @param cf_mat Logical matrix with k columns of indicators representing the different folds
#' (for example created by \code{\link{prep_cf_mat}}).
#' @param learner Vector of characters indicating whether to use S or T learner
#' or a combination of the two.
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
#' \item{ml}{list of fitted machine learning models from all cross-fitting folds
#' (if storeModels is not set to "No")}
#'
#' @export
#'
ensemble_short = function(ml,
                          x, y,
                          cf_mat,
                          subset = NULL,
                          learner = c("t", "s", "both"),
                          storeModels = c("No", "Memory", "Disk"),
                          path = NULL,
                          quiet = TRUE, pb = NULL, pb_np = NULL
                          ) {


  # check if subset vector is provided
  if(is.null(subset)) subset = rep(TRUE, nrow(x))

  # parameter configuration
  learner = match.arg(learner)
  if(all(subset == TRUE)) learner = "t"
  storeModels = match.arg(storeModels)
  saveModels = (storeModels != "No")

  s = NULL
  if(saveModels) s = vector("list", length = ncol(cf_mat))

  # matrix to store the cross-validated predictions
  fit_cv = make_fit_cv(ml = ml, n = nrow(x), learner = learner)

  for (i in 1:ncol(cf_mat)) {

    #if (isFALSE(quiet)) print(paste("Cross-fitting fold: ", toString(i)))

    fold = cf_mat[, i]
    x_tr = x[!fold, ]
    y_tr = y[!fold]
    x_te = x[fold, ]
    subset_tr = subset[!fold]


    if (learner %in% c("t", "both")) {

      ml_fit = ensemble_core(ml, x_tr, y_tr, quiet = quiet, pb = pb, pb_cf = i, pb_cv = ".", pb_np = pb_np)
      fit_cv[fold, grepl("^t", colnames(fit_cv))] = predict.ensemble_core(ml_fit, ml, x_tr, y_tr, x_te, quiet = quiet, pb = pb, pb_cf = i, pb_cv = ".", pb_np = pb_np)

      if(saveModels) s[[i]] = ml_fit

    }

    if (learner %in% c("s", "both")) {

      x_tr_s = cbind(x_tr, subset_tr*1)
      x_te_s = cbind(x_te, rep(1, nrow(x_te)))

      ml_fit = ensemble_core(ml, x_tr_s, y_tr, quiet = quiet, pb = pb, pb_cf = i, pb_cv = ".", pb_np = pb_np)
      fit_cv[fold, grepl("^s", colnames(fit_cv))] = predict.ensemble_core(ml_fit, ml, x_tr_s, y_tr, x_te_s, quiet = quiet, pb = pb, pb_cf = i, pb_cv = ".", pb_np = pb_np)

      if(saveModels) s[[ncol(cf_mat) + i]] = ml_fit

    }

  }

  # if(storeModels == "Disk") {
  #   if(is.null(path)) path = getwd()
  #   file_name = paste0(path, "/ensemble_short_", format(Sys.time(), "%Y%m%d%H%M%S"), ".rds")
  #   saveRDS(s, file_name)
  #   s = file_name
  # }

  output = list(
    "fit_cv" = fit_cv,
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
#' @param subset Logical vector indicating which observations to use for determining
#' ensemble weights. If not provided, all observations are used.
#' @param w Ensemble weights to aggregate predictions from different learners (optional).
#' Needs to be a vector of length \code{ncol(object$fit_cv)}.
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
                                  subset,
                                  w = NULL,
                                  cf_mat,
                                  quiet = TRUE,
                                  ...) {

  # load fitted ml object
  if(is.null(object$ml)) {
    stop("Ensemble models were not saved after training.")
  } else if (is.character(object$ml) & length(object$ml) == 1) {
    ml_fit = readRDS(object$ml)
  } else {
    ml_fit = object$ml
  }

  # estimate nnls weights if not provided
  if (is.null(w)) w = rep(1 / ncol(object$fit_cv), ncol(object$fit_cv))

  # create matrix to store smoother weights
  smoother_weights = matrix(0, nrow = nrow(x), ncol = nrow(x))

  # extract learner
  learner = unique((substr(colnames(object$fit_cv), 1, 1)))

  for (i in 1:ncol(cf_mat)) {

    fold = cf_mat[, i]
    x_tr = x[!fold, ]
    y_tr = y[!fold]
    x_te = x[fold, ]

    w_array = NULL
    if ("t" %in% learner) {
      w_array = weights.ensemble_core(object = ml_fit[[i]], ml = ml, x_tr = x_tr, y_tr = y_tr, x_te = x_te, quiet)
    }
    if ("s" %in% learner) {
      subset_tr = subset[!fold]
      x_tr_s = cbind(x_tr, subset_tr*1)
      x_te_s = cbind(x_te, rep(1, nrow(x_te)))
      w_array = abind::abind(w_array, weights.ensemble_core(object = ml_fit[[ncol(cf_mat) + i]], ml = ml, x_tr = x_tr_s, y_tr = y_tr, x_te = x_te_s, quiet))
    }

    smoother_weights[fold, !fold] = agg_array(w_array, w)

  }

  return(smoother_weights)

}
