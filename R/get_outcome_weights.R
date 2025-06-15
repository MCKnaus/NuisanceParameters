#' Arithmetic mean smoother weights
#'
#' @description
#' Returns smoother weights for test sample based on
#' arithmetic mean fitting of the training sample.
#'
#' @param mean_fit Output of \code{\link{mean_fit}}
#' @param xnew Covariate matrix of test sample
#' @param ... Ignore unused arguments
#'
#' @return Matrix of smoother weights.
#'
#' @method weights mean_fit
#'
#' @keywords internal
#'
weights.mean_fit = function(mean_fit, xnew, ...) {
  
  w = matrix(1 / mean_fit$n, nrow = nrow(xnew), ncol = mean_fit$n)
  
  return(w)
  
}


#' Smoother weights from OLS prediction
#'
#' @description
#' Extract smoother weights for test sample from an OLS model.
#'
#' @param ols_fit Output of \code{\link{ols_fit}}
#' @param x Covariate matrix of training sample
#' @param xnew Covariate matrix of test sample
#' @param ... Ignore unused arguments
#'
#' @return Matrix of smoother weights.
#'
#' @method weights ols_fit
#'
#' @keywords internal
#'
weights.ols_fit = function(ols_fit, x, xnew, ...) {
  
  x = add_intercept(x)
  xnew = add_intercept(xnew)
  
  # remove variables that were dropped due to collinearity
  x = x[, !is.na(ols_fit)]
  xnew = xnew[, !is.na(ols_fit)]
  
  # calculate hat matrix
  hat_mat = xnew %*% solve(crossprod(x), tol = 2.225074e-308) %*% t(x)
  
  return(hat_mat)
  
}


#' Smoother weights from Ridge regression prediction
#'
#' @description
#' Extract smoother weights for test sample from a Ridge regression model.
#'
#' @param ridge_fit Output of \code{\link{ridge_fit}}
#' @param x Covariate matrix of training sample
#' @param y Vector of outcomes of training sample
#' @param xnew Covariate matrix of test sample.
#' If not provided, prediction is done for the training sample.
#' @param ... Ignore unused arguments
#'
#' @return Matrix of smoother weights.
#'
#' @method weights ridge_fit
#'
#' @keywords internal
#'
weights.ridge_fit = function(ridge_fit, x, y, xnew = NULL, ...) {
  
  if (is.null(xnew)) xnew = x
  n = nrow(x)
  
  x = scale(x, ridge_fit$x_means, ridge_fit$x_sds)
  x = add_intercept(x)
  xnew = scale(xnew, ridge_fit$x_means, ridge_fit$x_sds)
  xnew = add_intercept(xnew)
  
  p = ncol(x) - 1
  
  sd_y = sqrt(stats::var(y) * ((n - 1) / n))
  lambda = (1 / sd_y) * ridge_fit$lambda.min * n
  
  # calculate hat matrix
  # reference: https://stats.stackexchange.com/questions/129179/why-is-glmnet-ridge-regression-giving-me-a-different-answer-than-manual-calculat
  hat_mat = xnew %*% solve(crossprod(x) + lambda * diag(x = c(0, rep(1, p)))) %*% t(x)
  
  return(hat_mat)
}


#' Smoother weights from Post-Lasso prediction
#'
#' @description
#' Extract smoother weights for test sample from a Post-Lasso regression model.
#'
#' @param plasso_fit Output of \code{\link{plasso_fit}}
#' @param xnew Covariate matrix of test sample.
#' If not provided, prediction is done for the training sample.
#' @param ... Ignore unused arguments
#'
#' @return Matrix of smoother weights.
#'
#' @method weights plasso_fit
#'
#' @keywords internal
#'
weights.plasso_fit = function(plasso_fit, xnew = NULL, ...) {
  
  if (is.null(xnew)) xnew = plasso_fit$x
  
  x = add_intercept(plasso_fit$x)
  xnew = add_intercept(xnew)
  
  colnames(x)[1] = "(Intercept)"
  colnames(xnew) = colnames(x)
  
  xact = x[, plasso_fit$names_pl, drop = FALSE]
  xactnew = xnew[, plasso_fit$names_pl, drop = FALSE]
  
  hat_mat = xactnew %*% solve(crossprod(xact), tol = 2.225074e-308) %*% t(xact)
  
  return(hat_mat)
}


#' Smoother weights from Random Forest model
#'
#' @description
#' Extract smoother weights for test sample from a Random Forest model.
#'
#' @param forest_grf_fit Output of \code{\link{forest_grf_fit}}
#' @param xnew Covariate matrix of test sample.
#' If not provided, prediction is done for the training sample.
#' @param ... Ignore unused arguments
#'
#' @return Matrix of smoother weights.
#'
#' @method weights forest_grf_fit
#'
#' @keywords internal
#'
weights.forest_grf_fit = function(forest_grf_fit, xnew = NULL, ...) {
  
  if(is.null(xnew)) xnew = forest_grf_fit$X.orig
  
  if (utils::packageVersion("grf") < "2.0.0") {
    w = grf::get_sample_weights(forest_grf_fit, newdata = xnew)
  } else {
    w = grf::get_forest_weights(forest_grf_fit, newdata = xnew)
  }
  w = as.matrix(w)
  
  return(w)
}


#' Smoother weights from k-Nearest-Neighbor algorithm
#'
#' @description
#' Extract smoother weights for k-Nearest Neighbor algorithm.
#' This comes quite naturally here as the weight will be \code{1 / k} for all
#' 'neighbors' and 0 for all 'non-neighbors' (for a given test set observation).
#'
#' @param arguments Output of \code{\link{knn_fit}}
#' @param x Covariate matrix of training sample
#' @param xnew Covariate matrix of test sample
#' @param ... Ignore unused arguments
#'
#' @return Matrix of smoother weights.
#'
#' @method weights knn_fit
#'
#' @keywords internal
#'
weights.knn_fit = function(arguments, x, xnew = NULL, ...) {
  
  if (is.null(xnew)) xnew = x
  if (is.null(arguments$k)) {
    k = 10
  } else if ( all.equal(arguments$k, as.integer(arguments$k))) {
    k = arguments$k
  } else {
    k = 10
  }
  
  distance = as.matrix(FastKNN::Distance_for_KNN_test(xnew, x))
  
  get_binary_vector = function(row) {
    min_indices = order(row)[1:k]
    binary_vector = rep(0, length(row))
    binary_vector[min_indices] = 1
    
    return(binary_vector)
  }
  
  w = apply(distance, 1, get_binary_vector)
  w = t(w) / k
  
  return(w)
}


#' Smoother weights from Distributional Random Forest model
#'
#' @description
#' Extract smoother (or adaptive nearest neighbor) weights for test sample from
#' Distributional Random Forest model.
#'
#' @param ridge_fit Output of \code{\link{forest_drf_fit}}
#' @param xnew Covariate matrix of test sample.
#' If not provided, prediction is done for the training sample.
#' @param ... Ignore unused arguments
#'
#' @return Matrix of smoother weights.
#'
#' @method weights forest_drf_fit
#'
#' @keywords internal
#'
weights.forest_drf_fit = function(forest_drf_fit, xnew = NULL, ...) {
  
  if(is.null(xnew)) xnew = forest_drf_fit$X.orig
  
  class(forest_drf_fit) = "drf"
  w = as.matrix(predict(forest_drf_fit, newdata = xnew)$weights)
  
  return(w)
}
