#' Extraction of smoother weights for ensemble learner
#'
#' @description
#' Extract smoother weights from short- or standard-stacked ensemble model 
#' created by \code{\link{ensemble_short}} or \code{\link{ensemble}}.
#'
#' @param object Stacked ensemble learner from \code{\link{ensemble_short}} or \code{\link{ensemble}}.
#' @param ml List of ML models by \code{\link{create_method}} that was used as
#' input for \code{\link{ensemble_short}} or \code{\link{ensemble}}.
#' @param x Covariate matrix of training sample.
#' @param y Vector of outcomes of training sample.
#' @param subset Logical vector indicating which observations to use for determining
#' ensemble weights. If not provided, all observations are used.
#' @param nnls_w Ensemble weights to aggregate predictions from different learners (optional).
#' Needs to be a vector of length \code{ncol(object$fit_cv)}.
#' @param cv Number of cross-validation folds when estimating ensemble model.
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
get_smoother = function(object,
                               ml,
                               x, y,
                               subset= NULL,
                               nnls_w = NULL,
                               cv = cv,
                               cf_mat,
                               quiet = TRUE,
                               ...) {
  
  if (is.null(object)) stop("Ensemble models were not saved after training.")
  if (is.null(subset)) subset = rep(TRUE, length(y))
  
  
  ### Short-Stacking ###
  if (cv == 1) {
    w = weights(object[["ens_object"]], ml = ml, x = x, y = y, subset = subset, w = object[["nnls_w"]], cf_mat = cf_mat, quiet = quiet)
    
    
    ### Standard-Stacking ###
    
  } else if (cv > 1) {
    w = matrix(0, nrow = nrow(x), ncol = nrow(x))
    
    for (i in 1:ncol(cf_mat)) {
      fold = cf_mat[, i]
      x_tr = x[!fold & subset, ]
      y_tr = y[!fold & subset]
      x_te = x[fold, ]
      
      w_array = weights(object[[i]][["ens_object"]], ml, x = x_tr, y = y_tr, xnew = x_te, quiet = quiet)
      w[fold, !fold & subset] = apply(w_array, c(1, 2), function(x) sum(x * object[[i]][["ens_object"]]$nnls_weights))
    }
    
  }
  
  return(w)
  
}

