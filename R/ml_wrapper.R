#' Arithmetic mean
#'
#' @description
#' \code{\link{mean_fit}} acts as prediction function by calculating the
#' arithmetic mean of the training outcomes.
#'
#' @param x Covariate matrix of training sample
#' @param y Vector of outcomes of training sample
#'
#' @return Returns list containing mean and number of observations.
#'
#' @keywords internal
#'
mean_fit = function(x, y) {
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
#' Predicts arithmetic mean and provides prediction weights if required.
#'
#' @param mean_fit Output of \code{\link{mean_fit}}
#' @param x Covariate matrix of training sample
#' @param y Vector of outcomes of training sample
#' @param xnew Covariate matrix of test sample
#' @param weights If TRUE, weights underlying the prediction for xnew calculated
#'
#' @return List element containing:
#' \item{prediction}{vector of predictions for xnew}
#' \item{weights}{if \code{weights=TRUE} prediction weights of dimension \code{nrow(xnew)} x \code{nrow(x)}
#' containing the weights that deliver predictions where each row gives the weight that each training
#' outcome received in the prediction for xnew.}
#'
#' @method predict mean_fit
#'
#' @keywords internal
#'
predict.mean_fit = function(mean_fit, x, y, xnew = NULL, weights = FALSE) {
  if (is.null(xnew)) fit = rep(mean_fit$mean, nrow(x))
  else fit = rep(mean_fit$mean, nrow(xnew))

  if (isTRUE(weights)) w = matrix(1 / length(y), nrow(xnew), nrow(x))
  else w = NULL

  return(
    list("prediction" = fit, "weights" = w)
  )
}


#' Fits OLS
#'
#' @description
#' \code{\link{ols_fit}} fits ordinary least squares (OLS) to the given data.
#'
#' @param x Covariate matrix of training sample
#' @param y Vector of outcomes of training sample
#'
#' @return Returns OLS coefficients
#'
#' @importFrom stats lm.fit
#'
#' @keywords internal
#'
ols_fit = function(x, y) {
  x = cbind(rep(1, nrow(x)), x)
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
#' @param xnew Covariate matrix of test sample
#' @param weights If TRUE, weights underlying the prediction for xnew calculated
#'
#' @return Returns list containing:
#' \item{prediction}{vector of predictions for xnew}
#' \item{weights}{If \code{weights=TRUE} prediction weights of dimension \code{nrow(xnew)} x \code{nrow(x)}
#' containing the weights that deliver predictions where each row gives the weight that each training
#' outcome received in the prediction for xnew.}
#'
#' @method predict ols_fit
#'
#' @keywords internal
#'
predict.ols_fit = function(ols_fit, x, y, xnew = NULL, weights = FALSE) {
  if (is.null(xnew)) xnew = x

  x = cbind(rep(1, nrow(x)), x)
  xnew = cbind(rep(1, nrow(xnew)), xnew)

  # remove variables that were dropped due to collinearity
  x = x[, !is.na(ols_fit)]
  xnew = xnew[, !is.na(ols_fit)]

  # calculate hat matrix
  hat_mat = xnew %*% solve(crossprod(x), tol = 2.225074e-308) %*% t(x)
  fit = hat_mat %*% y

  if (weights == FALSE) hat_mat = NULL

  return(
    list("prediction" = fit, "weights" = hat_mat)
  )
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
#' The method also provides prediction weights if required.
#'
#' @param ridge_fit Output of \code{\link{ridge_fit}}
#' @param x Covariate matrix of training sample
#' @param y Vector of outcomes of training sample
#' @param xnew Covariate matrix of test sample
#' @param weights If TRUE, weights underlying the prediction for xnew calculated
#'
#' @return Returns list containing:
#' \item{prediction}{vector of predictions for xnew}
#' \item{weights}{If \code{weights=TRUE} prediction weights of dimension \code{nrow(xnew)} x \code{nrow(x)}
#' containing the weights that deliver predictions where each row gives the weight that each training
#' outcome received in the prediction for xnew.}
#'
#' @importFrom stats sd
#' @importFrom glmnet predict.glmnet
#'
#' @method predict ridge_fit
#'
#' @keywords internal
#'
predict.ridge_fit = function(ridge_fit, x, y, xnew = NULL, weights = FALSE) {
  if (is.null(xnew)) xnew = x

  fit = glmnet::predict.glmnet(ridge_fit, newx = xnew, type = "response")

  if (weights == FALSE) hat_mat = NULL
  else {
    # Get covariate matrices
    n = nrow(x)
    p = ncol(x)
    x = scale(x)
    x = cbind(rep(1, nrow(x)), x)
    xnew = scale(xnew)
    xnew = cbind(rep(1, nrow(xnew)), xnew)

    # Calculate hat matrix, see also (https://stats.stackexchange.com/questions/129179/why-is-glmnet-ridge-regression-giving-me-a-different-answer-than-manual-calculat)
    hat_mat = xnew %*% solve(crossprod(x) + ridge_fit$lambda.min  * n / stats::sd(y) * diag(x = c(0, rep(1, p)))) %*% t(x)
    fit = hat_mat %*% y
  }

  return(
    list("prediction" = fit, "weights" = hat_mat)
  )
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
#' Prediction based on Post-Lasso and provides prediction weights if required.
#'
#' @param plasso_fit Output of \code{\link{plasso_fit}}
#' @param x Covariate matrix of training sample
#' @param y Vector of outcomes of training sample
#' @param xnew Covariate matrix of test sample
#' @param weights If TRUE, weights underlying the prediction for xnew calculated
#'
#' @return Returns list containing:
#' \item{prediction}{vector of predictions for xnew}
#' \item{weights}{If \code{weights=TRUE} prediction weights of dimension \code{nrow(xnew)} x \code{nrow(x)}
#' containing the weights that deliver predictions where each row gives the weight that each training
#' outcome received in the prediction for xnew.}
#'
#' @importFrom stats lm.fit
#'
#' @method predict plasso_fit
#'
#' @keywords internal
#'
predict.plasso_fit = function(plasso_fit, x, y, xnew = NULL, weights = FALSE) {
  if (is.null(xnew)) xnew = x
  x = add_intercept(x)
  xnew = add_intercept(xnew)

  # fitted values for post lasso
  nm_act = names(coef(plasso_fit$lasso_full)[, plasso_fit$ind_min_pl])[which(coef(plasso_fit$lasso_full)[, plasso_fit$ind_min_pl] != 0)]

  xact = x[, nm_act, drop = FALSE]
  xactnew = xnew[, nm_act, drop = FALSE]

  # remove potentially collinear variables
  coef = stats::lm.fit(xact, y)$coefficients
  xact = xact[, !is.na(coef)]
  xactnew = xactnew[, !is.na(coef)]

  hat_mat = xactnew %*% solve(crossprod(xact), tol = 2.225074e-308) %*% t(xact)
  fit_plasso = hat_mat %*% y
  if (weights == FALSE) hat_mat = NULL

  return(
    list("prediction" = fit_plasso, "weights" = hat_mat)
  )
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
#' Prediction based on Random Forest and provides prediction weights if
#' required.
#'
#' @param forest_grf_fit Output of \code{\link{forest_grf_fit}}
#' @param x Covariate matrix of training sample
#' @param y Vector of outcomes of training sample
#' @param xnew Covariate matrix of test sample
#' @param weights If TRUE, weights underlying the prediction for xnew calculated
#'
#' @return Returns list containing:
#' \item{prediction}{vector of predictions for xnew}
#' \item{weights}{If \code{weights=TRUE} prediction weights of dimension \code{nrow(xnew)} x \code{nrow(x)}
#' containing the weights that deliver predictions where each row gives the weight that each training
#' outcome received in the prediction for xnew.}
#'
#' @import grf
#' @importFrom utils packageVersion
#' @importFrom stats predict
#'
#' @method predict forest_grf_fit
#'
#' @keywords internal
#'
predict.forest_grf_fit = function(forest_grf_fit, x, y, xnew = NULL, weights = FALSE) {
  if (is.null(xnew)) xnew = x

  class(forest_grf_fit) = "regression_forest"
  fit = predict(forest_grf_fit, newdata = xnew)$prediction

  if (weights == TRUE) {
    if (utils::packageVersion("grf") < "2.0.0") w = grf::get_sample_weights(forest_grf_fit, newdata = xnew)
    else  w = grf::get_forest_weights(forest_grf_fit, newdata = xnew)
  }
  else w = NULL

  return(
    list("prediction" = fit, "weights" = w)
  )
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


#' Predictions based on Lasso regression
#'
#' @description
#' Prediction based on fitted Lasso regression model.
#' The method also provides prediction weights if required.
#'
#' @param lasso_fit Output of \code{\link{lasso_fit}}
#' @param x Covariate matrix of training sample
#' @param y Vector of outcomes of training sample
#' @param xnew Covariate matrix of test sample
#' @param weights Always FALSE as
#'
#' @return Returns list containing:
#' \item{prediction}{vector of predictions for xnew}
#' \item{weights}{Not available for Lasso, only for Post-Lasso}
#'
#' @importFrom glmnet predict.glmnet
#'
#' @method predict lasso_fit
#'
#' @keywords internal
#'
predict.lasso_fit = function(lasso_fit, x, y, xnew = NULL, weights = FALSE) {

  if (isTRUE(weights)) stop("No weighted representation of Lasso available.", call. = FALSE)
  if (is.null(xnew)) xnew = x

  class(lasso_fit) = "cv.glmnet"
  fit = predict(lasso_fit, newx = xnew, type = "response", s = "lambda.min")

  return(
    list("prediction" = fit, "weights" = NULL)
  )
}
