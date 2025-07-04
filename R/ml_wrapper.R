#' Arithmetic mean
#'
#' @description
#' \code{\link{mean_fit}} acts as prediction function by calculating the
#' arithmetic mean of the training outcomes.
#'
#' @param X Covariate matrix of training sample
#' @param Y Vector of outcomes of training sample
#' @param ... Ignore unused arguments
#'
#' @return Returns list containing mean and number of observations.
#'
#' @keywords internal
#'
mean_fit = function(X, Y, ...) {
  mean_value = mean(Y)
  output = list(
    "mean" = mean_value,
    "n" = nrow(X)
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
#' @param Xnew Covariate matrix of test sample. If no test sample provided,
#' prediction is done for training sample.
#' @param ... Ignore unused arguments
#'
#' @return Vector of fitted values.
#'
#' @method predict mean_fit
#'
#' @keywords internal
#'
predict.mean_fit = function(mean_fit, Xnew = NULL, ...) {

  if (is.null(Xnew)) {
    fit = rep(mean_fit$mean, mean_fit$n)
  } else {
    fit = rep(mean_fit$mean, nrow(Xnew))
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
#' @param Xnew Covariate matrix of test sample
#' @param ... Ignore unused arguments
#'
#' @return Matrix of smoother weights.
#'
#' @method weights mean_fit
#'
#' @keywords internal
#'
weights.mean_fit = function(mean_fit, Xnew, ...) {

  w = matrix(1 / mean_fit$n, nrow = nrow(Xnew), ncol = mean_fit$n)

  return(w)

}


#' Fits OLS
#'
#' @description
#' \code{\link{ols_fit}} fits ordinary least squares (OLS) to the given data.
#'
#' @param X Covariate matrix of training sample
#' @param Y Vector of outcomes of training sample
#' @param ... Ignore unused arguments
#'
#' @return Returns OLS coefficients
#'
#' @keywords internal
#'
ols_fit = function(X, Y, ...) {
  X = add_intercept(X)
  ols_coef = stats::lm.fit(X, Y)$coefficients
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
#' @param Xnew Covariate matrix of test sample.
#' @param ... Ignore unused arguments
#'
#' @return Vector of fitted values.
#'
#' @method predict ols_fit
#'
#' @keywords internal
#'
predict.ols_fit = function(ols_fit, Xnew = NULL, ...) {

  Xnew = add_intercept(Xnew)
  Xnew = Xnew[, !is.na(ols_fit)]

  fit = as.vector(Xnew %*% matrix(ols_fit, ncol = 1))

  return(fit)

}

#' Smoother weights from OLS prediction
#'
#' @description
#' Extract smoother weights for test sample from an OLS model.
#'
#' @param ols_fit Output of \code{\link{ols_fit}}
#' @param X Covariate matrix of training sample
#' @param Xnew Covariate matrix of test sample
#' @param ... Ignore unused arguments
#'
#' @return Matrix of smoother weights.
#'
#' @method weights ols_fit
#'
#' @keywords internal
#'
weights.ols_fit = function(ols_fit, X, Xnew, ...) {

  X = add_intercept(X)
  Xnew = add_intercept(Xnew)

  # remove variables that were dropped due to collinearity
  X = X[, !is.na(ols_fit)]
  Xnew = Xnew[, !is.na(ols_fit)]

  # calculate hat matrix
  hat_mat = Xnew %*% solve(crossprod(X), tol = 2.225074e-308) %*% t(X)

  return(hat_mat)

}


#' Fits Ridge regression
#'
#' @description
#' \code{\link{ridge_fit}} estimates cross-validated Ridge regression based on
#' the \code{\link{glmnet}} package.
#'
#' @param X Matrix of covariates
#' @param Y Vector of outcomes
#' @param arguments List of arguments passed to \code{\link[glmnet]{glmnet}}
#'
#' @return An object with S3 class \code{\link[glmnet]{glmnet}}
#'
#' @keywords internal
#'
ridge_fit = function(X, Y, arguments = list()) {
  x_means = matrixStats::colMeans2(X)
  x_sds = matrixStats::colSds(X)
  x_std = scale(X, x_means, x_sds)
  ridge = do.call(glmnet::cv.glmnet, c(list(x = x_std, y = Y, alpha = 0, standardize = FALSE, intercept = TRUE), arguments))
  ridge[["x_means"]] = x_means
  ridge[["x_sds"]] = x_sds
  class(ridge) = "ridge_fit"
  return(ridge)
}


#' Predictions based on Ridge regression
#'
#' @description
#' Prediction based on fitted Ridge regression model.
#'
#' @param ridge_fit Output of \code{\link{ridge_fit}}
#' @param Xnew Covariate matrix of test sample.
#' @param ... Ignore unused arguments
#'
#' @return Vector of predictions for Xnew.
#'
#' @method predict ridge_fit
#'
#' @keywords internal
#'
predict.ridge_fit = function(ridge_fit, Xnew, ...) {

  Xnew = scale(Xnew, ridge_fit$x_means, ridge_fit$x_sds)

  class(ridge_fit) = "cv.glmnet"

  fit = as.vector(predict(ridge_fit, newx = Xnew, s = "lambda.min"))

  return(fit)
}


#' Smoother weights from Ridge regression prediction
#'
#' @description
#' Extract smoother weights for test sample from a Ridge regression model.
#'
#' @param ridge_fit Output of \code{\link{ridge_fit}}
#' @param X Covariate matrix of training sample
#' @param Y Vector of outcomes of training sample
#' @param Xnew Covariate matrix of test sample.
#' If not provided, prediction is done for the training sample.
#' @param ... Ignore unused arguments
#'
#' @return Matrix of smoother weights.
#'
#' @method weights ridge_fit
#'
#' @keywords internal
#'
weights.ridge_fit = function(ridge_fit, X, Y, Xnew = NULL, ...) {

  if (is.null(Xnew)) Xnew = X
  n = nrow(X)

  X = scale(X, ridge_fit$x_means, ridge_fit$x_sds)
  X = add_intercept(X)
  Xnew = scale(Xnew, ridge_fit$x_means, ridge_fit$x_sds)
  Xnew = add_intercept(Xnew)

  p = ncol(X) - 1

  sd_y = sqrt(stats::var(Y) * ((n - 1) / n))
  lambda = (1 / sd_y) * ridge_fit$lambda.min * n

  # calculate hat matrix
  # reference: https://stats.stackexchange.com/questions/129179/why-is-glmnet-ridge-regression-giving-me-a-different-answer-than-manual-calculat
  hat_mat = Xnew %*% solve(crossprod(X) + lambda * diag(X = c(0, rep(1, p)))) %*% t(X)

  return(hat_mat)
}


#' Fits Post-Lasso regression
#'
#' @description
#' \code{\link{plasso_fit}} estimates cross-validated Post-Lasso regression.
#'
#' @param X Matrix of covariates (number of observations times number of covariates matrix)
#' @param Y vector of outcomes
#' @param arguments List of arguments passed to \code{\link[plasso]{cv.plasso}}
#'
#' @return An object with S3 class \code{\link[plasso]{cv.plasso}}
#'
#' @keywords internal
#'
plasso_fit = function(X, Y, arguments = list()) {
  p = do.call(plasso::cv.plasso, c(list(x = X, y = Y), arguments))
  class(p) = "plasso_fit"
  return(p)
}


#' Predictions based on Post-Lasso regression
#'
#' @description
#' Prediction of fitted values (for a potentially new set of covaraites Xnew)
#' based on a trained Post-Lasso model.
#'
#' @param plasso_fit Output of \code{\link{plasso_fit}}
#' @param Xnew Covariate matrix of test sample.
#' If not provided, prediction is done for the training sample.
#' @param ... Ignore unused arguments
#'
#' @return Vector of fitted values.
#'
#' @method predict plasso_fit
#'
#' @keywords internal
#'
predict.plasso_fit = function(plasso_fit, Xnew = NULL, ...) {

  class(plasso_fit) = "cv.plasso"
  if (is.null(Xnew)) Xnew = plasso_fit$X

  fit = as.vector(predict(plasso_fit, newx = Xnew, type = "response", s = "optimal", se_rule = 0)$plasso)

  return(fit)
}


#' Smoother weights from Post-Lasso prediction
#'
#' @description
#' Extract smoother weights for test sample from a Post-Lasso regression model.
#'
#' @param plasso_fit Output of \code{\link{plasso_fit}}
#' @param Xnew Covariate matrix of test sample.
#' If not provided, prediction is done for the training sample.
#' @param ... Ignore unused arguments
#'
#' @return Matrix of smoother weights.
#'
#' @method weights plasso_fit
#'
#' @keywords internal
#'
weights.plasso_fit = function(plasso_fit, Xnew = NULL, ...) {

  if (is.null(Xnew)) Xnew = plasso_fit$X

  X = add_intercept(plasso_fit$X)
  Xnew = add_intercept(Xnew)

  colnames(X)[1] = "(Intercept)"
  colnames(Xnew) = colnames(X)

  xact = X[, plasso_fit$names_pl, drop = FALSE]
  xactnew = Xnew[, plasso_fit$names_pl, drop = FALSE]

  hat_mat = xactnew %*% solve(crossprod(xact), tol = 2.225074e-308) %*% t(xact)

  return(hat_mat)
}


#' Fits Random Forest
#'
#' @description
#' \code{\link{forest_grf_fit}} fits a random forest using the
#' \code{\link{grf}} package.
#'
#' @param X Matrix of covariates
#' @param Y vector of outcomes
#' @param arguments List of arguments passed to \code{\link[grf]{regression_forest}}
#'
#' @return An object with S3 class \code{\link[grf]{regression_forest}}
#'
#' @keywords internal
#'
forest_grf_fit = function(X, Y, arguments = list()) {
  rf = do.call(grf::regression_forest, c(list(X = X, Y = Y), arguments))
  class(rf) = "forest_grf_fit"
  return(rf)
}


#' Predictions based on Random Forest
#'
#' @description
#' Prediction of fitted values (for a potentially new set of covariates Xnew)
#' from a Random Forest.
#'
#' @param forest_grf_fit Output of \code{\link{forest_grf_fit}}
#' @param Xnew Covariate matrix of test sample.
#' If not provided, prediction is done for the training sample.
#' @param ... Ignore unused arguments
#'
#' @return Vector of fitted values.
#'
#' @method predict forest_grf_fit
#'
#' @keywords internal
#'
predict.forest_grf_fit = function(forest_grf_fit, Xnew = NULL, ...) {

  if(is.null(Xnew)) Xnew = forest_grf_fit$X.orig

  class(forest_grf_fit) = "regression_forest"
  fit = predict(forest_grf_fit, newdata = Xnew)$prediction

  return(fit)
}


#' Smoother weights from Random Forest model
#'
#' @description
#' Extract smoother weights for test sample from a Random Forest model.
#'
#' @param forest_grf_fit Output of \code{\link{forest_grf_fit}}
#' @param Xnew Covariate matrix of test sample.
#' If not provided, prediction is done for the training sample.
#' @param ... Ignore unused arguments
#'
#' @return Matrix of smoother weights.
#'
#' @method weights forest_grf_fit
#'
#' @keywords internal
#'
weights.forest_grf_fit = function(forest_grf_fit, Xnew = NULL, ...) {

  if(is.null(Xnew)) Xnew = forest_grf_fit$X.orig

  if (utils::packageVersion("grf") < "2.0.0") {
    w = grf::get_sample_weights(forest_grf_fit, newdata = Xnew)
  } else {
    w = grf::get_forest_weights(forest_grf_fit, newdata = Xnew)
  }
  w = as.matrix(w)

  return(w)
}


#' Fits Lasso regression
#'
#' @description
#' \code{\link{lasso_fit}} estimates cross-validated Lasso regression based on
#' the \code{\link{glmnet}} package.
#'
#' @param X Matrix of covariates (number of observations times number of covariates matrix)
#' @param Y vector of outcomes
#' @param arguments List of arguments passed to \code{\link[glmnet]{cv.glmnet}}
#'
#' @return An object with S3 class \code{\link[glmnet]{cv.glmnet}}
#'
#' @keywords internal
#'
lasso_fit = function(X, Y, arguments = list()) {
  lasso = do.call(glmnet::cv.glmnet, c(list(x = X, y = Y), arguments))
  class(lasso) = "lasso_fit"
  return(lasso)
}


#' Lasso prediction
#'
#' @description
#' Prediction of fitted values based on fine-tuned Lasso regression model.
#'
#' @param lasso_fit Output of \code{\link{lasso_fit}}
#' @param Xnew Covariate matrix of test sample.
#' If not provided, prediction is done for the training sample.
#' @param ... Ignore unused arguments
#'
#' @return Vector of fitted values.
#'
#' @method predict lasso_fit
#'
#' @keywords internal
#'
predict.lasso_fit = function(lasso_fit, Xnew = NULL, ...) {

  class(lasso_fit) = "cv.glmnet"
  if (is.null(Xnew)) Xnew = lasso_fit$X

  fit = predict(lasso_fit, newx = Xnew, type = "response")

  return(fit)
}


#' Pseudo fitting of k-Nearest-Neighbour to training data
#'
#' @description
#' \code{\link{knn_fit}} is a pseudo fitting function where only the kNN object
#' is initialized but not fitted as the fitting process comes naturally with
#' the prediction part.
#'
#' @param arguments arguments of arguments passed to \code{\link[FastKNN]{k.nearest.neighbors}}
#' @param ... Ignore unused arguments
#'
#' @return A list of possible arguments for \code{\link[FastKNN]{k.nearest.neighbors}}.
#'
#' @keywords internal
#'
knn_fit = function(arguments = list(), ...) {
  class(arguments) = "knn_fit"
  return(arguments)
}


#' Smoother weights from k-Nearest-Neighbor algorithm
#'
#' @description
#' Extract smoother weights for k-Nearest Neighbor algorithm.
#' This comes quite naturally here as the weight will be \code{1 / k} for all
#' 'neighbors' and 0 for all 'non-neighbors' (for a given test set observation).
#'
#' @param arguments Output of \code{\link{knn_fit}}
#' @param X Covariate matrix of training sample
#' @param Xnew Covariate matrix of test sample
#' @param ... Ignore unused arguments
#'
#' @return Matrix of smoother weights.
#'
#' @method weights knn_fit
#'
#' @keywords internal
#'
weights.knn_fit = function(arguments, X, Xnew = NULL, ...) {

  if (is.null(Xnew)) Xnew = X
  if (is.null(arguments$k)) {
    k = 10
  } else if ( all.equal(arguments$k, as.integer(arguments$k))) {
    k = arguments$k
  } else {
    k = 10
  }

  distance = as.matrix(FastKNN::Distance_for_KNN_test(Xnew, X))

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


#' Prediction of fitted values based on the k-Nearest-Neighbor algorithm
#'
#' @description
#' Predictions by k-Nearest-Neighbor algorithm.
#' Note that if no \code{k} is explicitly specified, this evaluates to
#' \code{k = 10} as default value.
#'
#' @param arguments Output of \code{\link{knn_fit}}
#' @param X Covariate matrix of training sample
#' @param Y Vector of outcomes of training sample
#' @param Xnew Covariate matrix of test sample
#' @param ... Ignore unused arguments
#'
#' @return Vector of fitted values.
#'
#' @method predict knn_fit
#'
#' @keywords internal
#'
predict.knn_fit = function(arguments, X, Y, Xnew = NULL, ...) {

  # get smoother weights
  w = weights.knn_fit(arguments, X = X, Xnew = Xnew)
  # multiply with training outcome vector
  fit = as.vector(w %*% Y)

  return(fit)
}


#' Fits Distributional Random Forest
#'
#' @description
#' \code{\link{forest_drf_fit}} fits a distributional random forest using the
#' \code{\link{drf}} package.
#'
#' @param X Matrix of covariates
#' @param Y vector of outcomes
#' @param arguments List of arguments passed to \code{\link[drf]{drf}}
#'
#' @return An object with S3 class \code{\link[drf]{drf}}
#'
#' @keywords internal
#'
forest_drf_fit = function(X, Y, arguments = list()) {
  rf = do.call(drf::drf, c(list(X = X, Y = Y, splitting.rule = "FourierMMD"), arguments))
  class(rf) = "forest_drf_fit"
  return(rf)
}


#' Predictions based on Distributional Random Forest
#'
#' @description
#' Prediction of fitted values (for a new set of covariates)
#' from a Random Forest.
#'
#' @param forest_drf_fit Output of \code{\link{forest_drf_fit}}
#' @param Xnew Covariate matrix of test sample.
#' If not provided, prediction is done for the training sample.
#' @param ... Ignore unused arguments
#'
#' @return Vector of fitted values.
#'
#' @method predict forest_drf_fit
#'
#' @keywords internal
#'
predict.forest_drf_fit = function(forest_drf_fit, Xnew = NULL, functional = "mean", ...) {

  if(is.null(Xnew)) Xnew = forest_drf_fit$X.orig

  class(forest_drf_fit) = "drf"
  pred = predict(forest_drf_fit, newdata = Xnew, functional = functional)
  fit = as.vector(pred[[functional]])

  return(fit)
}


#' Smoother weights from Distributional Random Forest model
#'
#' @description
#' Extract smoother (or adaptive nearest neighbor) weights for test sample from
#' Distributional Random Forest model.
#'
#' @param ridge_fit Output of \code{\link{forest_drf_fit}}
#' @param Xnew Covariate matrix of test sample.
#' If not provided, prediction is done for the training sample.
#' @param ... Ignore unused arguments
#'
#' @return Matrix of smoother weights.
#'
#' @method weights forest_drf_fit
#'
#' @keywords internal
#'
weights.forest_drf_fit = function(forest_drf_fit, Xnew = NULL, ...) {

  if(is.null(Xnew)) Xnew = forest_drf_fit$X.orig

  class(forest_drf_fit) = "drf"
  w = as.matrix(predict(forest_drf_fit, newdata = Xnew)$weights)

  return(w)
}
