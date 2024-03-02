#' Propensity score nuisance parameters prediction
#'
#' @description
#' \code{\link{nuisance_e}} makes a cross-fitted ensemble prediction of
#' propensity score nuisance parameters.
#'
#' @param ml List of methods to be used in \code{\link{ensemble}} estimation of propensity score.
#' Methods can be created by \code{\link{create_method}}.
#' @param w_mat Logical matrix of treatment indicators (n x T+1). For example created by \code{\link{prep_w_mat}}.
#' @param x Covariate matrix.
#' @param cf_mat Logical matrix with k columns of indicators representing the different folds
#' (for example created by \code{\link{prep_cf_mat}}).
#' @param cv Number of cross-validation when estimating ensemble (default 5)
#' @param path Optional path to save the \code{\link{ensemble}} objects for later inspection.
#' Saved as Ensemble_Wi where i is the number of the treatment in multiple treatment settings.
#' @param quiet If FALSE, method that is currently running is printed into console.
#'
#' @return Returns n x T+1 matrix with each column containing the propensity score for the treatment corresponding to w_mat.
#'
#' @export
#'
nuisance_e = function(ml,
                      w_mat, x, cf_mat,
                      cv = 5,
                      path = NULL,
                      quiet = TRUE) {

  if (isFALSE(quiet)) message("Propensity score")
  if (isFALSE(quiet)) which_stacking(cv)

  # define path if not provided
  if (is.null(path)) path = getwd()


  # initialize nuisance matrix
  e_mat = matrix(NA, nrow(w_mat), ncol(w_mat))
  colnames(e_mat) = colnames(w_mat)



  ### binary treatment case ###

  if (ncol(w_mat) == 2) {
    path_temp = paste0(path, "/Ensemble_W")
    e_mat[, 1] = nuisance_cf(ml = ml, y = w_mat[, 1], x = x, cf_mat = cf_mat, cv = cv, path = path_temp, quiet = quiet,
                             learner = "t", subset = NULL, weights = FALSE)
    e_mat[, 2] = 1 - e_mat[, 1]
  }



  ### multiple treatment case ###

  else if (ncol(w_mat) > 2) {
    path_temp = paste0(path, "/Ensemble_W", 1:ncol(w_mat))
    for (i in 1:ncol(w_mat)) {
      e_mat[, i] = nuisance_cf(ml = ml, y = w_mat[, i], x = x, cf_mat = cf_mat, cv = cv, path = path_temp[i], quiet = quiet,
                               learner = "t", subset = NULL, weights = FALSE)
    }
    e_mat = e_mat / rowSums(e_mat)
  }

  else {
    stop("Provide treatment indicator matrix with at least 2 columns")
  }

  return(e_mat)

}


#' Outcome model nuisance parameters prediction
#'
#' @description
#' \code{\link{nuisance_m}} makes a cross-fitted ensemble prediction of
#' the outcome regression nuisance parameters.
#'
#' @param ml List of methods to be used in \code{\link{ensemble}} estimation of
#' propensity score.
#' Methods can be created by \code{\link{create_method}}.
#' @param y Numerical vector containing the outcome variable.
#' @param w_mat Logical matrix of treatment indicators (n x T+1). For example
#' created by \code{\link{prep_w_mat}}.
#' @param x Covariate matrix.
#' @param cf_mat Logical matrix with k columns of indicators representing the different folds
#' (for example created by \code{\link{prep_cf_mat}})
#' @param learner Vector of characters indicating whether to use S or T learner
#' or a combination of the two.
#' @param cv Number of cross-validation when estimating ensemble (default 5).
#' @param weights If TRUE, prediction weights of the outcome nuisance extracted
#' and saved (requires to provide a path).
#' @param path Optional path to save the \code{\link{ensemble}} objects for later processing.
#' Saved as Ensemble_Yi where i is the number of the treatment in multiple treatment settings.
#' @param quiet If FALSE, method that is currently running is printed into console.
#'
#' @return Returns n x T+1 matrix with each column containing the predicted outcome for the treatment corresponding to w_mat.
#'
#' @export
#'
nuisance_m = function(ml, y, w_mat, x, cf_mat,
                      learner = c("t", "s", "both"),
                      cv = 5,
                      weights = FALSE,
                      path = NULL,
                      quiet = TRUE) {

  # learner configuration
  learner = match.arg(learner)

  if((learner %in% c("s", "both")) & (ncol(w_mat) > 2)) stop("S-Learner cannot be combined with more than two treatments.")

  if (isFALSE(quiet)) message("Outcome regression")
  if (isFALSE(quiet)) which_stacking(cv)

  # define path if not provided
  if (is.null(path)) path = getwd()
  path_temp = paste0(path, "/Ensemble_Y", 1:ncol(w_mat))

  # initialize nuisance matrix
  m_mat = matrix(NA, nrow(w_mat), ncol(w_mat))
  colnames(m_mat) = colnames(w_mat)

  # calculate outcome predictions
  for (i in 1:ncol(w_mat)) {

    m_mat[,i] = nuisance_cf(ml = ml, y = y, x = x, cf_mat = cf_mat, learner = learner, cv = cv,
                            subset = (w_mat[, i]), weights = weights, path = path_temp[i], quiet = quiet)

  }

  return(m_mat)
}


#' Cross-fitting of nuisance parameters
#'
#' @description
#' \code{\link{nuisance_cf}} makes a cross-fitted ensemble prediction using
#' \code{\link{ensemble}}.
#'
#' @param ml List of methods to be used in \code{\link{ensemble}} estimation.
#' Methods can be created by \code{\link{create_method}}.
#' @param y Vector of variable to be predicted.
#' @param x Matrix of covariates.
#' @param cf_mat Logical matrix with k columns of indicators representing the different folds
#' (for example created by \code{\link{prep_cf_mat}}).
#' @param learner Vector of characters indicating whether to use S or T learner
#' or a combination of the two.
#' @param cv Number of cross-validation folds when estimating ensemble model.
#' Default value is 1 which then evaluates to a short-stacking procedure which
#' is computationally less demanding than standard stacking.
#' @param subset Optional logical vector if only subset of data should be used for prediction.
#' @param weights If TRUE, prediction weights of the outcome nuisance extracted and saved.
#' @param path Path to save fit_cv matrix and non-negative least square weights to.
#' Optionally, weights as well.
#' @param quiet If FALSE, method that is currently running is printed into console.
#'
#' @return Returns a vector of length n containing nuisance parameters.
#'
#' @keywords internal
#'
nuisance_cf = function(ml, y, x, cf_mat,
                       learner = c("t", "s", "both"),
                       cv = 5,
                       subset = NULL,
                       weights = FALSE,
                       path,
                       quiet = TRUE) {


  ### Parameter Configuration ###

  learner = match.arg(learner)



  ### Checks ###

  if (is.numeric(cf_mat)) {
    if (!all(cf_mat %in% 0:1) | !is.matrix(cf_mat)) stop("Please provide cf_mat as binary indicator matrix. E.g. use function prep_cf_mat")
    if (nrow(cf_mat) != length(y)) stop("cf_mat indicator matrix nrows different from # of obs.")
    if (length(y) != sum(cf_mat)) stop("cf_mat indicator matrix does not sum to number of observations.")
  }
  if (isTRUE(weights) & is.null(path)) stop("Provide path if weights = TRUE to save ensemble objects with weights for later processing or set weights = FALSE.")

  if (is.null(subset)) subset = rep(TRUE, length(y))

  np = rep(NA, length(y))



  ### Short-Stacking ###

  if (cv == 1) {

    if(isTRUE(weights)) storeModels = "Memory" else storeModels = "No"

    ens = ensemble_short(ml = ml, x = x, y = y, subset = subset, cf_mat = cf_mat, learner = learner, storeModels = storeModels, path = path, quiet = quiet)
    nnls_w = nnls_weights(ens$fit_cv[subset, ], y[subset])
    np = predict(ens, w = nnls_w)

    fit_sub = list("fit_cv" = ens$fit_cv, "nnls_w" = nnls_w)
    class(fit_sub) = "ens.learner"

    saveRDS(fit_sub, paste0(path, ".rds"))

    if(isTRUE(weights)) {
      w = weights(ens, ml = ml, x = x, y = y, subset = subset, w = nnls_w, cf_mat = cf_mat, quiet = FALSE)
      saveRDS(w, paste0(path, "_Weights.rds")); rm(w);
    }

  ### Standard-Stacking ###

  } else if (cv > 1) {

    if(isTRUE(weights)) w = matrix(0, nrow = nrow(x), ncol = nrow(x))

    fit_sub = list()

    for (i in 1:ncol(cf_mat)) {

      if (isFALSE(quiet)) print(paste("Cross-fitting fold: ", toString(i)))

      fold = cf_mat[, i]
      x_tr = x[!fold & subset, ]
      y_tr = y[!fold & subset]
      x_te = x[fold, ]

      ens = ensemble(ml = ml, x = x_tr, y = y_tr, nfolds = cv, quiet = quiet)
      ens_pred = predict(ens, ml, x = x_tr, y = y_tr, xnew = x_te, quiet = quiet)
      np[fold] = ens_pred$np

      fit_sub[[i]] = list("fit_cv" = ens_pred$fit_cv, "nnls_w" = ens$nnls_weights)

      if(isTRUE(weights)) {
        w_array = weights(ens, ml, x = x_tr, y = y_tr, xnew = x_te, quiet = quiet)
        w[fold, !fold & subset] = apply(w_array, c(1, 2), function(x) sum(x * ens$nnls_weights))
      }

    }

    class(fit_sub) = "ens.learner"
    saveRDS(fit_sub, paste0(path, ".rds")); rm(fit_sub)

    if(isTRUE(weights)) {
      saveRDS(w, paste0(path, "_Weights.rds")); rm(w);
    }

  }

  return(np)

}
