#' Fit arithmetic mean model
#'
#' Calculates the arithmetic mean of the training outcomes for use as a constant prediction model.
#'
#' @param X Covariate matrix (ignored, for compatibility).
#' @param Y Numeric vector of training outcomes.
#' @param ... Ignored additional arguments.
#'
#' @return An object of class `mean_fit` containing:
#'   \item{mean}{The arithmetic mean of `Y`}
#'   \item{N}{Number of observations}
#'
#' @keywords internal
mean_fit <- function(X, Y, ...) {
  mean_value <- mean(Y)
  output <- list("mean" = mean_value, "N" = nrow(X))
  
  class(output) <- "mean_fit"
  return(output)
}


#' Predict method for arithmetic mean model
#'
#' Returns a constant prediction for all observations, equal to the mean of the training sample.
#' Covariates `Xnew` are ignored.
#'
#' @param mean_fit An object of class `mean_fit` from \code{\link{mean_fit}}.
#' @param Xnew Covariate matrix for predictions. Only `nrow(Xnew)` is used.
#' @param ... Ignored additional arguments.
#'
#' @return A numeric vector of predictions, all equal to the trained mean.
#'
#' @method predict mean_fit
#' @keywords internal
predict.mean_fit <- function(mean_fit, Xnew = NULL, ...) {
  if (is.null(Xnew)) {
    fit <- rep(mean_fit$mean, mean_fit$N)
  } else {
    fit <- rep(mean_fit$mean, nrow(Xnew))
  }

  return(fit)
}


#' Fit ordinary least squares model
#'
#' Fits an ordinary least squares (OLS) regression to the given data using
#' \code{stats::lm.fit()}.
#'
#' @param X Numeric matrix of covariates for training sample.
#' @param Y Numeric vector of outcomes for training sample.
#' @param ... Ignored additional arguments.
#'
#' @return Numeric vector of OLS coefficients with class \code{ols_fit}.
#'
#' @keywords internal
ols_fit <- function(X, Y, ...) {
  X <- add_intercept(X)
  ols_coef <- stats::lm.fit(X, Y)$coefficients
  
  class(ols_coef) <- "ols_fit"
  return(ols_coef)
}


#' Predict method for OLS models
#'
#' Generates predictions from a fitted ordinary least squares model.
#'
#' @param ols_fit An object of class \code{ols_fit} from \code{\link{ols_fit}}.
#' @param Xnew Numeric matrix of covariates for test sample. If \code{NULL},
#'   predictions are generated for training data.
#' @param ... Ignored additional arguments.
#'
#' @return Numeric vector of predicted values.
#'
#' @method predict ols_fit
#' @keywords internal
predict.ols_fit <- function(ols_fit, Xnew = NULL, ...) {
  Xnew <- add_intercept(Xnew)
  Xnew <- Xnew[, !is.na(ols_fit)]

  fit <- as.vector(Xnew %*% matrix(ols_fit, ncol = 1))

  return(fit)
}


#' Fit ridge regression model
#'
#' Estimates ridge regression using the \code{glmnet} package.
#' Standardizes covariates and selects optimal lambda via
#' cross-validation.
#'
#' @param X Numeric matrix of covariates for training sample.
#' @param Y Numeric vector of outcomes for training sample.
#' @param arguments List of arguments passed to \code{\link[glmnet]{cv.glmnet}}.
#'
#' @return An object of class \code{ridge_fit} containing the fitted
#'   \code{\link[glmnet]{cv.glmnet}} model.
#'
#' @keywords internal
ridge_fit <- function(X, Y, arguments = list()) {
  if (!requireNamespace("glmnet", quietly = TRUE)) {
    stop(
      "The 'glmnet' package is not installed. 'Ridge' learner is optional (in Suggests).\n",
      "Install it with: install.packages('glmnet')", 
      call. = FALSE
    )
  }

  x_means <- colMeans(X, na.rm = TRUE)
  x_sds <- apply(X, 2, stats::sd, na.rm = TRUE)
  x_std <- scale(X, x_means, x_sds)
  ridge <- do.call(glmnet::cv.glmnet, c(list(x = x_std, y = Y, alpha = 0, standardize = FALSE, intercept = TRUE), arguments))
  ridge[["x_means"]] <- x_means
  ridge[["x_sds"]] <- x_sds

  class(ridge) <- c(class(ridge), "ridge_fit")
  return(ridge)
}


#' Predict method for ridge regression models
#'
#' Generates predictions from a fitted ridge regression model.
#' Applies the same standardization used during training to new data.
#'
#' @param ridge_fit An object of class \code{ridge_fit} from \code{\link{ridge_fit}}.
#' @param Xnew Numeric matrix of covariates for test sample.
#' @param ... Ignored additional arguments.
#'
#' @return Numeric vector of predictions for \code{Xnew}.
#'
#' @method predict ridge_fit
#' @keywords internal
predict.ridge_fit <- function(ridge_fit, Xnew, ...) {
  Xnew <- scale(Xnew, ridge_fit$x_means, ridge_fit$x_sds)
  fit <- as.vector(predict(ridge_fit, newx = Xnew, s = "lambda.min"))

  return(fit)
}


#' Fit post-lasso regression model
#'
#' Estimates cross-validated post-lasso regression using the \code{plasso} package.
#' Performs lasso selection followed by OLS on selected variables.
#'
#' @param X Numeric matrix of covariates for training sample.
#' @param Y Numeric vector of outcomes for training sample.
#' @param arguments List of arguments passed to \code{\link[plasso]{cv.plasso}}.
#'
#' @return An object of class \code{plasso_fit} containing the fitted
#'   \code{\link[plasso]{cv.plasso}} model.
#'
#' @keywords internal
plasso_fit <- function(X, Y, arguments = list()) {
  if (!requireNamespace("plasso", quietly = TRUE)) {
    stop(
      "The 'plasso' package is not installed. This learner is optional (in Suggests).\n",
      "Install it with: install.packages('plasso')", 
      call. = FALSE
    )
  }

  plasso <- do.call(plasso::cv.plasso, c(list(x = X, y = Y), arguments))
  class(plasso) <- c(class(plasso), "plasso_fit")

  return(plasso)
}


#' Predict method for post-lasso regression models
#'
#' Generates predictions from a fitted post-lasso regression model.
#' Returns predictions for either new data or the training sample.
#'
#' @param plasso_fit An object of class \code{plasso_fit} from \code{\link{plasso_fit}}.
#' @param Xnew Numeric matrix of covariates for test sample. If \code{NULL},
#'   predictions are generated for the training data.
#' @param ... Ignored additional arguments.
#'
#' @return Numeric vector of predicted values.
#'
#' @method predict plasso_fit
#' @keywords internal
predict.plasso_fit <- function(plasso_fit, Xnew = NULL, ...) {
  if (is.null(Xnew)) Xnew <- plasso_fit$x
  fit <- as.vector(predict(plasso_fit, newx = Xnew, type = "response", s = "optimal", se_rule = 0)$plasso)

  return(fit)
}


#' Fit post-lasso regression using hdm package
#'
#' Estimates post-lasso regression under homoscedastic and heteroscedastic
#' non-Gaussian disturbances using the \code{hdm} package.
#'
#' @param X Numeric matrix of covariates for training sample.
#' @param Y Numeric vector of outcomes for training sample.
#' @param arguments List of arguments passed to \code{\link[hdm]{rlasso}}.
#'
#' @return An object of class \code{rlasso_fit} containing the fitted
#'   \code{\link[hdm]{rlasso}} model.
#'
#' @keywords internal
rlasso_fit <- function(X, Y, arguments = list()) {
  if (!requireNamespace("hdm", quietly = TRUE)) {
    stop(
      "The 'hdm' package is not installed. 'rlasso' learner is optional (in Suggests).\n",
      "Install it with: install.packages('hdm')",
      call. = FALSE
    )
  }

  rlasso <- do.call(hdm::rlasso, c(list(x = X, y = Y), arguments))
  class(rlasso) <- c(class(rlasso), "rlasso_fit")

  return(rlasso)
}


#' Predict method for hdm rlasso models
#'
#' Generates predictions from a fitted post-lasso regression model
#' estimated using the \code{hdm} package.
#'
#' @param rlasso_fit An object of class \code{rlasso_fit} from \code{\link{rlasso_fit}}.
#' @param Xnew Numeric matrix of covariates for test sample. If \code{NULL},
#'   predictions are generated for the training data.
#' @param ... Ignored additional arguments.
#'
#' @return Numeric vector of predicted values.
#'
#' @method predict rlasso_fit
#' @keywords internal
predict.rlasso_fit <- function(rlasso_fit, Xnew = NULL, ...) {
  fit <- as.vector(predict(rlasso_fit, newdata = Xnew))

  return(fit)
}


#' Fit random forest using grf package
#'
#' Fits a random forest for regression using the \code{grf} package.
#'
#' @param X Numeric matrix of covariates for training sample.
#' @param Y Numeric vector of outcomes for training sample.
#' @param arguments List of arguments passed to \code{\link[grf]{regression_forest}}.
#'
#' @return An object of class \code{forest_grf_fit} containing the fitted
#'   \code{\link[grf]{regression_forest}} model.
#'
#' @keywords internal
forest_grf_fit <- function(X, Y, arguments = list()) {
  if (!requireNamespace("grf", quietly = TRUE)) {
    stop(
      "The 'grf' package is not installed. 'regression_forest' learner is optional (in Suggests).\n",
      "Install it with: install.packages('grf')", 
      call. = FALSE
    )
  }
  
  rf <- do.call(grf::regression_forest, c(list(X = X, Y = Y), arguments))
  
  class(rf) <- c(class(rf), "forest_grf_fit")
  return(rf)
}


#' Predict method for grf random forest models
#'
#' Generates predictions from a fitted random forest model
#' estimated using the \code{grf} package.
#'
#' @param forest_grf_fit An object of class \code{forest_grf_fit} from 
#'   \code{\link{forest_grf_fit}}.
#' @param Xnew Numeric matrix of covariates for test sample. If \code{NULL},
#'   predictions are generated for the training data.
#' @param ... Ignored additional arguments.
#'
#' @return Numeric vector of predicted values.
#'
#' @method predict forest_grf_fit
#' @keywords internal
predict.forest_grf_fit <- function(forest_grf_fit, Xnew = NULL, ...) {
  if (is.null(Xnew)) Xnew <- forest_grf_fit$X.orig
  fit <- predict(forest_grf_fit, newdata = Xnew)$prediction

  return(fit)
}


#' Fit XGBoost model
#'
#' Fits an XGBoost model using the \code{xgboost} package.
#'
#' @param X Numeric matrix of covariates for training sample.
#' @param Y Numeric vector of outcomes for training sample.
#' @param arguments List of arguments passed to \code{\link[xgboost]{xgb.train}}.
#'
#' @details
#' Parameter restrictions for the smoother matrix multiplied by the outcome 
#' vector to be (roughly) equivalent to the outcome nuisance vector 
#' \itemize{
#'   \item \strong{Restricted}: \code{reg_alpha = 0}, \code{subsample = 1}, 
#'         \code{max_delta_step = 0}, \code{base_score = 0}, 
#'         \code{objective = "reg:squarederror"}, \code{booster = "gbtree"}
#'   \item \strong{Free}: \code{nrounds}, \code{eta}, \code{gamma}, \code{max_depth},
#'         \code{min_child_weight}, \code{sampling_method}, \code{colsample_bytree},
#'         \code{colsample_bylevel}, \code{colsample_bynode}, \code{lambda},
#'         \code{tree_method}, \code{refresh_leaf}, \code{grow_policy},
#'         \code{max_leaves}, \code{max_bin}
#'   \item \strong{Grey zone}: All other parameters (equivalence not guaranteed)
#' }
#'
#' @return An object of class \code{xgboost_fit} containing the fitted
#'   \code{\link[xgboost]{xgb.Booster}} model.
#'
#' @keywords internal
xgboost_fit <- function(X, Y, arguments = list()) {
  if (!requireNamespace("xgboost", quietly = TRUE)) {
    stop(
      "The 'xgboost' package is not installed. This learner is optional (in Suggests).\n",
      "Install it with: install.packages('xgboost')",
      call. = FALSE
    )
  }
  
  free_params <- c(
    "nrounds", "eta", "gamma", "max_depth", "min_child_weight", "sampling_method",
    "colsample_bytree", "colsample_bylevel", "colsample_bynode", "lambda", 
    "tree_method", "refresh_leaf", "grow_policy", "max_leaves", "max_bin"
  )
  
  fixed <- list(
    alpha          = 0,
    subsample      = 1,
    max_delta_step = 0,
    base_score     = 0,
    objective      = "reg:squarederror",
    booster        = "gbtree"
  )
  
  # sanitize restricted params
  for (nm in names(fixed)) {
    if (!is.null(arguments[[nm]]) && !identical(arguments[[nm]], fixed[[nm]])) {
      message("xgboost_fit: resetting ", nm, " to ", fixed[[nm]],
              " (user-specified value ignored).")
    }
    arguments[[nm]] <- fixed[[nm]]
  }
  
  # grey zone check
  grey <- setdiff(names(arguments), c(free_params, names(fixed)))
  if (length(grey) > 0) {
    message("xgboost_fit: detected grey-zone parameter(s): ",
            paste(grey, collapse = ", "), ". S*Y = Y.hat is not guaranteed.")
  }
  
  if (is.null(arguments$eta)) arguments$eta <- 0.3
  if (is.null(arguments$lambda)) arguments$lambda <- 1

  nrounds <- if (is.null(arguments[["nrounds"]])) 100 else nrounds
  arguments[["nrounds"]] <- NULL
  
  dtrain <- xgboost::xgb.DMatrix(data = as.matrix(X), label = Y)
  xgb <- do.call(xgboost::xgb.train, c(list(data = dtrain, nrounds = nrounds, params = arguments)))
  
  class(xgb) = c(class(xgb), "xgboost_fit")
  return(xgb)
}


#' Predict method for XGBoost models
#'
#' Generates predictions from a fitted XGBoost model.
#'
#' @param xgboost_fit An object of class \code{xgboost_fit} from \code{\link{xgboost_fit}}.
#' @param Xnew Numeric matrix of covariates for test sample. If \code{NULL},
#'   predictions are generated for the training data.
#' @param ... Ignored additional arguments.
#'
#' @return Numeric vector of predicted values.
#'
#' @method predict xgboost_fit
#' @keywords internal
predict.xgboost_fit <- function(xgboost_fit, Xnew = NULL, ...) {
  dtest <- xgboost::xgb.DMatrix(data = as.matrix(Xnew))
  fit <- predict(xgboost_fit, newdata = dtest)

  return(fit)
}


#' Fit lasso regression model
#'
#' Estimates lasso regression using the \code{glmnet} package.
#' Selects optimal lambda penalty parameter via cross-validation.
#'
#' @param X Numeric matrix of covariates for training sample.
#' @param Y Numeric vector of outcomes for training sample.
#' @param arguments List of arguments passed to \code{\link[glmnet]{cv.glmnet}}.
#'
#' @return An object of class \code{lasso_fit} containing the fitted
#'   \code{\link[glmnet]{cv.glmnet}} model.
#'
#' @keywords internal
lasso_fit <- function(X, Y, arguments = list()) {
  if (!requireNamespace("glmnet", quietly = TRUE)) {
    stop(
      "The 'glmnet' package is not installed. 'Lasso' learner is optional (in Suggests).\n",
      "Install it with: install.packages('glmnet')",
      call. = FALSE
    )
  }

  lasso <- do.call(glmnet::cv.glmnet, c(list(x = X, y = Y), arguments))
  
  class(lasso) = c(class(lasso), "lasso_fit")
  return(lasso)
}


#' Predict method for lasso regression models
#'
#' Generates predictions from a fitted lasso regression model
#' using the optimal lambda selected during cross-validation.
#'
#' @param lasso_fit An object of class \code{lasso_fit} from \code{\link{lasso_fit}}.
#' @param Xnew Numeric matrix of covariates for test sample. If \code{NULL},
#'   predictions are generated for the training data.
#' @param ... Ignored additional arguments.
#'
#' @return Numeric vector of predicted values.
#'
#' @method predict lasso_fit
#' @keywords internal
predict.lasso_fit <- function(lasso_fit, Xnew = NULL, ...) {
  if (is.null(Xnew)) Xnew <- lasso_fit$X
  fit <- predict(lasso_fit, newx = Xnew, type = "response", s = "lambda.min")

  return(fit)
}


#' Fit distributional random forest model
#'
#' Fits a distributional random forest using the \code{drf} package with
#' Fourier MMD splitting rule.
#'
#' @param X Numeric matrix of covariates for training sample.
#' @param Y Numeric vector of outcomes for training sample.
#' @param arguments List of arguments passed to \code{\link[drf]{drf}}.
#'
#' @return An object of class \code{forest_drf_fit} containing the fitted
#'   \code{\link[drf]{drf}} model.
#'
#' @keywords internal
forest_drf_fit <- function(X, Y, arguments = list()) {
  if (!requireNamespace("drf", quietly = TRUE)) {
    stop(
      "The 'drf' package is not installed. This learner is optional (in Suggests).\n",
      "Install it with: install.packages('drf')",
      call. = FALSE
    )
  }

  rf <- do.call(drf::drf, c(list(X = X, Y = Y, splitting.rule = "FourierMMD"), arguments))
  
  class(rf) <- c(class(rf), "forest_drf_fit")
  return(rf)
}


#' Predict method for distributional random forest models
#'
#' Generates predictions from a fitted distributional random forest model.
#' By default, returns predictions for the mean functional, but supports
#' other distributional characteristics.
#'
#' @param forest_drf_fit An object of class \code{forest_drf_fit} from 
#'   \code{\link{forest_drf_fit}}.
#' @param Xnew Numeric matrix of covariates for test sample. If \code{NULL},
#'   predictions are generated for the training data.
#' @param functional Character specifying the functional to predict. Defaults
#'   to \code{"mean"} for mean predictions.
#' @param ... Ignored additional arguments.
#'
#' @return Numeric vector of predicted values.
#'
#' @method predict forest_drf_fit
#' @keywords internal
predict.forest_drf_fit <- function(forest_drf_fit, Xnew = NULL, functional = "mean", ...) {
  if (is.null(Xnew)) Xnew <- forest_drf_fit$X.orig

  pred <- predict(forest_drf_fit, newdata = Xnew, functional = functional)
  fit <- as.vector(pred[[functional]])

  return(fit)
}


#' Pseudo-fitting for k-Nearest Neighbors
#'
#' Creates a k-NN model object without actual fitting, as k-NN fitting occurs
#' naturally during prediction.
#'
#' @param arguments List of arguments to be used in k-NN prediction. Typically
#'   includes \code{k} (number of neighbors).
#' @param ... Ignored additional arguments.
#'
#' @return A list of class \code{knn_fit} containing the arguments for k-NN.
#'
#' @keywords internal
knn_fit <- function(arguments = list(), ...) {
  class(arguments) <- "knn_fit"
  return(arguments)
}


#' Predict method for k-Nearest Neighbors
#'
#' Generates predictions using the k-Nearest Neighbors algorithm with internal
#' Euclidean distance calculation, following the same approach as the
#' \code{Distance_for_KNN_test()} function of the \code{FastKNN} package.
#'
#' @param arguments Output of \code{\link{knn_fit}} of class \code{knn_fit}. 
#'  List of arguments to be used in k-NN prediction. Typically includes \code{k} 
#'  (number of neighbors). If no \code{k} is specified, defaults to \code{k = 10}.
#' @param X Numeric matrix of covariates for training sample.
#' @param Y Numeric vector of outcomes for training sample.
#' @param Xnew Numeric matrix of covariates for test sample. If \code{NULL},
#'   predictions are generated for the training data.
#' @param ... Ignored additional arguments.
#'
#' @return Numeric vector of predicted values.
#'
#' @method predict knn_fit
#' @keywords internal
predict.knn_fit <- function(arguments, X, Y, Xnew = NULL, ...) {
  if (is.null(Xnew)) Xnew = X
  k <- if (is.null(arguments$k)) 10 else arguments$k
  
  # distance = as.matrix(FastKNN::Distance_for_KNN_test(Xnew, X))
  euclidean_dist <- t(apply(Xnew, 1, function(x_i) {
    sqrt(rowSums(sweep(X, 2, x_i)^2))
    }))
    
  get_binary_vec = function(row) {
    min_indices = order(row)[1:k]
    binary_vec = rep(0, length(row))
    binary_vec[min_indices] = 1
    return(binary_vec)
    }
    
    w = apply(euclidean_dist, 1, get_binary_vec)
    w = t(w) / k
    
    fit <- as.vector(w %*% Y)
  
  return(fit)
}