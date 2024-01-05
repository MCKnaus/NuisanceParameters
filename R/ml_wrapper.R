#' Arithmetic mean
#'
#' @description
#' \code{\link{mean_fit}} acts as prediction function by calculating the
#' arithmetic mean of the training outcomes.
#'
#' @param x Covariate matrix of training sample
#' @param y Vector of outcomes of training sample
#' @param ... Ignore unused arguments
#'
#' @return Returns list containing mean and number of observations.
#'
#' @keywords internal
#'
mean_fit = function(x, y, ...) {
  mean_value = mean(y)
  output = list(
    "mean" = mean_value,
    "n" = nrow(x)
  )
  class(output) = "mean_fit"
  return(output)
}


#' Predict after arithmetic mean fitting
#'
#' @description
#' Predicts arithmetic mean.
#'
#' @param mean_fit Output of \code{\link{mean_fit}}
#' @param xnew Covariate matrix of test sample. If no test sample provided,
#' prediction is done for training sample.
#' @param ... Ignore unused arguments
#'
#' @return Vector of fitted values.
#'
#' @method predict mean_fit
#'
#' @keywords internal
#'
predict.mean_fit = function(mean_fit, xnew = NULL, ...) {

  if (is.null(xnew)) {
    fit = rep(mean_fit$mean, mean_fit$n)
  } else {
    fit = rep(mean_fit$mean, nrow(xnew))
  }

  return(fit)

}


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


#' Fits OLS
#'
#' @description
#' \code{\link{ols_fit}} fits ordinary least squares (OLS) to the given data.
#'
#' @param x Covariate matrix of training sample
#' @param y Vector of outcomes of training sample
#' @param ... Ignore unused arguments
#'
#' @return Returns OLS coefficients
#'
#' @importFrom stats lm.fit
#'
#' @keywords internal
#'
ols_fit = function(x, y, ...) {
  x = add_intercept(x)
  ols_coef = stats::lm.fit(x, y)$coefficients
  class(ols_coef) = "ols_fit"
  return(ols_coef)
}


#' Predictions based on OLS
#'
#' @description
#' Prediction based on fitted ordinary least squares models. The method also
#' provides prediction weights if required.
#'
#' @param ols_fit Output of \code{\link{ols_fit}}
#' @param x Covariate matrix of training sample
#' @param y Vector of outcomes of training sample
#' @param xnew Covariate matrix of test sample. If no test sample provided,
#' prediction is done for training sample.
#' @param ... Ignore unused arguments
#'
#' @return Vector of fitted values.
#'
#' @method predict ols_fit
#'
#' @keywords internal
#'
predict.ols_fit = function(ols_fit, x, y, xnew = NULL, ...) {

  if (is.null(xnew)) xnew = x

  xnew = add_intercept(xnew)
  xnew = xnew[, !is.na(ols_fit)]

  fit = xnew %*% matrix(ols_fit, ncol = 1)

  return(fit)

}

#' Smoother weights from OLS prediction
#'
#' @description
#' Extract smoother weights for test sample from an OLS model.
#'
#' @param ols_fit Output of \code{\link{ols_fit}}
#' @param x Covariate matrix of training sample
#' @param y Vector of outcomes of training sample
#' @param xnew Covariate matrix of test sample
#' @param ... Ignore unused arguments
#'
#' @return Matrix of smoother weights.
#'
#' @method weights ols_fit
#'
#' @keywords internal
#'
weights.ols_fit = function(ols_fit, x, y, xnew, ...) {

  x = add_intercept(x)
  xnew = add_intercept(xnew)

  # remove variables that were dropped due to collinearity
  x = x[, !is.na(ols_fit)]
  xnew = xnew[, !is.na(ols_fit)]

  # calculate hat matrix
  hat_mat = xnew %*% solve(crossprod(x), tol = 2.225074e-308) %*% t(x)

  return(hat_mat)

}


#' Fits Ridge regression
#'
#' @description
#' \code{\link{ridge_fit}} estimates cross-validated Ridge regression based on
#' the \code{\link{glmnet}} package.
#'
#' @param x Matrix of covariates
#' @param y vector of outcomes
#' @param args List of arguments passed to \code{\link[glmnet]{glmnet}}
#'
#' @return An object with S3 class \code{\link[glmnet]{glmnet}}
#'
#' @importFrom glmnet cv.glmnet
#'
#' @keywords internal
#'
ridge_fit = function(x, y, args = list()) {
  ridge = do.call(glmnet::cv.glmnet, c(list(x = x, y = y, alpha = 0), args))
  class(ridge) = "ridge_fit"
  return(ridge)
}


#' Predictions based on Ridge regression
#'
#' @description
#' Prediction based on fitted Ridge regression model.
#'
#' @param ridge_fit Output of \code{\link{ridge_fit}}
#' @param xnew Covariate matrix of test sample.
#' If not provided, prediction is done for the training sample.
#' @param ... Ignore unused arguments
#'
#' @return Vector of predictions for xnew.
#'
#' @importFrom glmnet cv.glmnet
#'
#' @method predict ridge_fit
#'
#' @keywords internal
#'
predict.ridge_fit = function(ridge_fit, xnew = NULL, ...) {

  class(ridge_fit) = "cv.glmnet"
  if (is.null(xnew)) xnew = ridge_fit$x

  fit = predict(ridge_fit, newx = xnew, type = "response")

  return(fit)
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
weights.ridge_fit = function(ridge_fit, x, y, xnew, ...) {
  n = nrow(x)
  p = ncol(x)
  if (is.null(xnew)) xnew = x

  x = scale(x)
  x = add_intercept(x)
  xnew = scale(xnew)
  xnew = add_intercept(xnew)

  # calculate hat matrix
  # reference: https://stats.stackexchange.com/questions/129179/why-is-glmnet-ridge-regression-giving-me-a-different-answer-than-manual-calculat
  hat_mat = xnew %*% solve(crossprod(x) + ridge_fit$lambda.min * n / stats::sd(y) * diag(x = c(0, rep(1, p)))) %*% t(x)

  return(hat_mat)
}


#' Fits Post-Lasso regression
#'
#' @description
#' \code{\link{plasso_fit}} estimates cross-validated Post-Lasso regression.
#'
#' @param x Matrix of covariates (number of observations times number of covariates matrix)
#' @param y vector of outcomes
#' @param args List of arguments passed to \code{\link[plasso]{cv.plasso}}
#'
#' @importFrom plasso cv.plasso
#'
#' @return An object with S3 class \code{\link[plasso]{cv.plasso}}
#'
#' @keywords internal
#'
plasso_fit = function(x, y, args = list()) {
  p = do.call(plasso::cv.plasso, c(list(x = x, y = y), args))
  class(p) = "plasso_fit"
  return(p)
}


#' Predictions based on Post-Lasso regression
#'
#' @description
#' Prediction of fitted values (for a potentially new set of covaraites xnew)
#' based on a trained Post-Lasso model.
#'
#' @param plasso_fit Output of \code{\link{plasso_fit}}
#' @param xnew Covariate matrix of test sample.
#' If not provided, prediction is done for the training sample.
#' @param ... Ignore unused arguments
#'
#' @return Vector of fitted values.
#'
#' @importFrom plasso cv.plasso
#'
#' @method predict plasso_fit
#'
#' @keywords internal
#'
predict.plasso_fit = function(plasso_fit, xnew = NULL, ...) {

  class(plasso_fit) = "cv.plasso"
  if (is.null(xnew)) xnew = plasso_fit$x

  fit = predict(plasso_fit, newx = xnew, type = "response", s = "optimal", se_rule = 0)

  return(fit)
}


#' Smoother weights from Post-Lasso prediction
#'
#' @description
#' Extract smoother weights for test sample from a Post-Lasso regression model.
#'
#' @param ridge_fit Output of \code{\link{plasso_fit}}
#' @param x Covariate matrix of training sample
#' @param y Vector of outcomes of training sample
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
weights.plasso_fit = function(plasso_fit, x, y, xnew = NULL, ...) {

  if (is.null(xnew)) xnew = x
  x = add_intercept(x)
  xnew = add_intercept(xnew)

  nm_act = names(coef(plasso_fit$lasso_full)[, plasso_fit$ind_min_pl])[which(coef(plasso_fit$lasso_full)[, plasso_fit$ind_min_pl] != 0)]

  xact = x[, nm_act, drop = FALSE]
  xactnew = xnew[, nm_act, drop = FALSE]

  coef = stats::lm.fit(xact, y)$coefficients
  xact = xact[, !is.na(coef)]
  xactnew = xactnew[, !is.na(coef)]

  hat_mat = xactnew %*% solve(crossprod(xact), tol = 2.225074e-308) %*% t(xact)

  return(hat_mat)
}


#' Fits Random Forest
#'
#' @description
#' \code{\link{forest_grf_fit}} fits a random forest using the
#' \code{\link{grf}} package.
#'
#' @param x Matrix of covariates
#' @param y vector of outcomes
#' @param args List of arguments passed to \code{\link[grf]{regression_forest}}
#'
#' @importFrom grf regression_forest
#'
#' @return An object with S3 class \code{\link[grf]{regression_forest}}
#'
#' @keywords internal
#'
forest_grf_fit = function(x, y, args = list()) {
  rf = do.call(grf::regression_forest, c(list(X = x, Y = y), args))
  class(rf) = "forest_grf_fit"
  return(rf)
}


#' Predictions based on Random Forest
#'
#' @description
#' Prediction of fitted values (for a potentially new set of covaraites xnew)
#' from a Random Forest.
#'
#' @param plasso_fit Output of \code{\link{forest_grf_fit}}
#' @param xnew Covariate matrix of test sample.
#' If not provided, prediction is done for the training sample.
#' @param ... Ignore unused arguments
#'
#' @return Vector of fitted values.
#'
#' @import grf
#'
#' @method predict forest_grf_fit
#'
#' @keywords internal
#'
predict.forest_grf_fit = function(forest_grf_fit, xnew = NULL, ...) {

  if (is.null(xnew)) xnew = x

  class(forest_grf_fit) = "regression_forest"
  fit = predict(forest_grf_fit, newdata = xnew)$prediction

  return(fit)
}


#' Smoother weights from Random Forest model
#'
#' @description
#' Extract smoother weights for test sample from a Random Forest model.
#'
#' @param ridge_fit Output of \code{\link{forest_grf_fit}}
#' @param xnew Covariate matrix of test sample.
#' If not provided, prediction is done for the training sample.
#' @param ... Ignore unused arguments
#'
#' @return Matrix of smoother weights.
#'
#' @import grf
#'
#' @method weights forest_grf_fit
#'
#' @keywords internal
#'
weights.forest_grf_fit = function(forest_grf_fit, xnew, ...) {

  if (utils::packageVersion("grf") < "2.0.0") {
    w = grf::get_sample_weights(forest_grf_fit, newdata = xnew)
  } else {
    w = grf::get_forest_weights(forest_grf_fit, newdata = xnew)
  }

  return(w)
}


#' Fits Lasso regression
#'
#' @description
#' \code{\link{lasso_fit}} estimates cross-validated Lasso regression based on
#' the \code{\link{glmnet}} package.
#'
#' @param x Matrix of covariates (number of observations times number of covariates matrix)
#' @param y vector of outcomes
#' @param args List of arguments passed to \code{\link[glmnet]{cv.glmnet}}
#'
#' @return An object with S3 class \code{\link[glmnet]{cv.glmnet}}
#'
#' @importFrom glmnet cv.glmnet
#'
#' @keywords internal
#'
lasso_fit = function(x, y, args = list()) {
  lasso = do.call(glmnet::cv.glmnet, c(list(x = x, y = y), args))
  class(lasso) = "lasso_fit"
  return(lasso)
}


#' Lasso prediction
#'
#' @description
#' Prediction of fitted values based on fine-tuned Lasso regression model.
#'
#' @param lasso_fit Output of \code{\link{lasso_fit}}
#' @param xnew Covariate matrix of test sample.
#' If not provided, prediction is done for the training sample.
#' @param ... Ignore unused arguments
#'
#' @return Vector of fitted values.
#'
#' @importFrom glmnet cv.glmnet
#'
#' @method predict lasso_fit
#'
#' @keywords internal
#'
predict.lasso_fit = function(lasso_fit, xnew = NULL, ...) {

  class(lasso_fit) = "cv.glmnet"
  if (is.null(xnew)) xnew = lasso_fit$x

  fit = predict(lasso_fit, newx = xnew, type = "response")

  return(fit)
}


#' Pseudo fitting of k-Nearest-Neighbour to training data
#'
#' @description
#' \code{\link{knn_fit}} is a pseudo fitting function where only the kNN object
#' is initialized but not fitted as the fitting process comes naturally with
#' the prediction part.
#'
#' @param args List of arguments passed to \code{\link[FastKNN]{k.nearest.neighbors}}
#' @param ... Ignore unused arguments
#'
#' @return A list of possible arguments for \code{\link[FastKNN]{k.nearest.neighbors}}
#'
#'
#' @keywords internal
#'
knn_fit = function(args = list(), ...) {
  class(args) = "knn_fit"
  return(args)
}


#' Prediction of fitted values based on the k-Nearest-Neighbor algorithm
#'
#' @description
#' Predictions by k-Nearest-Neighbor algorithm.
#' Note that if no \code{k} is explicitly specified, this evaluates to
#' \code{k = 10} as default value.
#'
#' @param args Output of \code{\link{knn_fit}}
#' @param x Covariate matrix of training sample
#' @param y Vector of outcomes of training sample
#' @param xnew Covariate matrix of test sample
#' @param ... Ignore unused arguments
#'
#' @return Vector of fitted values.
#'
#' @importFrom FastKNN Distance_for_KNN_test k.nearest.neighbors
#'
#' @method predict knn_fit
#'
#' @keywords internal
#'
predict.knn_fit = function(args, x, y, xnew = NULL, ...) {

  # get smoother weights
  w = weights(args, x = x, xnew = xnew)
  # multiply with training outcome vector
  fit = (w %*% y) / k

  return(fit)
}


#' Smoother weights from k-Nearest-Neighbor algorithm
#'
#' @description
#' Extract smoother weights for k-Nearest Neighbor algorithm.
#' This comes quite naturally here as the weight will be \code{1 / k} for all
#' 'neighbors' and 0 for all 'non-neighbors' (for a given test set observation).
#'
#' @param args Output of \code{\link{knn_fit}}
#' @param x Covariate matrix of training sample
#' @param xnew Covariate matrix of test sample
#' @param ... Ignore unused arguments
#'
#' @return Matrix of smoother weights.
#'
#' @importFrom FastKNN Distance_for_KNN_test k.nearest.neighbors
#'
#' @method weights knn_fit
#'
#' @keywords internal
#'
weights.knn_fit = function(args, x, xnew = NULL, ...) {

  if (is.null(xnew)) xnew = x
  if (is.null(args[["k"]])) {
    k = 10
  } else if ( all.equal(args[["k"]], as.integer(args[["k"]])) ) {
    k = args[["k"]]
  } else {
    k = 10
  }

  distance = as.matrix(FastKNN::Distance_for_KNN_test(xnew, x))

  knn_search = function(row_index, distance_matrix, k) {
    binary_vector = rep(0, ncol(distance_matrix))
    indices = FastKNN::k.nearest.neighbors(row_index, distance_matrix, k = k)
    binary_vector[indices] = 1
    return(binary_vector)
  }

  w = t(sapply(1:nrow(distance), FUN = knn_search, distance_matrix = distance, k = k))

  return(w)
}
