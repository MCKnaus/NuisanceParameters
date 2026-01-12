#' Fit arithmetic mean model
#'
#' Calculates the arithmetic mean of the training outcomes for use as a
#' constant prediction model.
#'
#' @param X Covariate matrix (ignored, for compatibility).
#' @param Y Numeric vector of training outcomes.
#' @param ... Ignored additional arguments.
#'
#' @return An object of class \code{mean_fit} containing:
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
#' Returns a constant prediction for all observations, equal to the mean of the
#' training sample. Covariates \code{Xnew} are ignored.
#'
#' @param object An object of class \code{mean_fit} from \code{\link{mean_fit}}.
#' @param Xnew Covariate matrix for predictions. Only \code{nrow(Xnew)} is used.
#' @param ... Ignored additional arguments.
#'
#' @return A numeric vector of predictions, all equal to the trained mean.
#'
#' @method predict mean_fit
#' @keywords internal
#' @exportS3Method
predict.mean_fit <- function(object, Xnew = NULL, ...) {
  if (is.null(Xnew)) {
    fit <- rep(object$mean, object$N)
  } else {
    fit <- rep(object$mean, nrow(Xnew))
  }

  return(fit)
}


#' Fit ordinary least squares model
#'
#' Fits an ordinary least squares (OLS) regression to the given data using
#' \code{stats::lm.fit()}. An intercept is always added to the covariate matrix.
#'
#' @param X Numeric matrix of covariates for training sample.
#' @param Y Numeric vector of outcomes for training sample.
#' @param ... Ignored additional arguments.
#'
#' @return An object of class \code{ols_fit} containing estimated coefficients
#'         and the training covariate matrix.
#'
#' @keywords internal
ols_fit <- function(X, Y, ...) {
  X <- add_intercept(X)
  ols_coef <- stats::lm.fit(X, Y)$coefficients

  output <- list(ols_coef = ols_coef, X = X)

  class(output) <- "ols_fit"
  return(output)
}


#' Predict method for OLS models
#'
#' Generates predictions from a fitted ordinary least squares model.
#'
#' @param object An object of class \code{ols_fit} from \code{\link{ols_fit}}.
#' @param Xnew Numeric matrix of covariates for test sample. If \code{NULL},
#'             predictions are generated for training data.
#' @param ... Ignored additional arguments.
#'
#' @return Numeric vector of predicted values.
#'
#' @method predict ols_fit
#' @keywords internal
#' @exportS3Method
predict.ols_fit <- function(object, Xnew = NULL, ...) {
  if (is.null(Xnew)) Xnew <- object$X

  Xnew <- add_intercept(Xnew)
  Xnew <- Xnew[, !is.na(object)]

  fit <- as.vector(Xnew %*% matrix(object$ols_coef, ncol = 1))

  return(fit)
}


#' Fit ridge regression model
#'
#' Estimates ridge regression using the \code{glmnet} package.
#' Standardizes covariates and selects optimal lambda via cross-validation.
#'
#' @param X Numeric matrix of covariates for training sample.
#' @param Y Numeric vector of outcomes for training sample.
#' @param arguments List of arguments passed to \code{cv.glmnet}.
#'
#' @return An object of class \code{ridge_fit} containing the fitted
#'   \code{cv.glmnet} model.
#'
#' @keywords internal
ridge_fit <- function(X, Y, arguments = list()) {
  if (!requireNamespace("glmnet", quietly = TRUE)) {
    stop("The 'glmnet' package is not installed. 'Ridge' learner is optional (in Suggests).\n",
      "Install it with: install.packages('glmnet')",
      call. = FALSE
    )
  }

  x_means <- colMeans(X)
  x_sds <- sqrt(apply(X, 2, stats::var))
  x_std <- scale(X, x_means, x_sds)

  ridge <- do.call(glmnet::cv.glmnet, c(list(
    x = x_std, y = Y, alpha = 0,
    standardize = FALSE, intercept = TRUE,
    thresh = 1e-20
  ), arguments))
  ridge$x_means <- x_means
  ridge$x_sds <- x_sds

  class(ridge) <- c(class(ridge), "ridge_fit")
  return(ridge)
}


#' Predict method for ridge regression models
#'
#' Generates predictions from a fitted ridge regression model.
#' Applies the same standardization used during training to new data.
#'
#' @param object An object of class \code{ridge_fit} from \code{\link{ridge_fit}}.
#' @param Xnew Numeric matrix of covariates for test sample.
#' @param ... Ignored additional arguments.
#'
#' @return Numeric vector of predictions for \code{Xnew}.
#'
#' @method predict ridge_fit
#' @keywords internal
#' @exportS3Method
predict.ridge_fit <- function(object, Xnew = NULL, ...) {
  X <- object[["call"]][["x"]]
  Y <- object[["call"]][["y"]]
  
  if (is.null(Xnew)) {
    Xnew <- X
    Xnew <- add_intercept(Xnew)
  } else {
    Xnew <- scale(Xnew, object$x_means, object$x_sds)
    Xnew <- add_intercept(Xnew)
  }

  X <- add_intercept(X)
  N <- nrow(X)
  p <- ncol(X) - 1

  sd_y <- sqrt(stats::var(Y) * ((N - 1) / N))
  lambda <- (1 / sd_y) * object$lambda.min * N

  # Reference: https://stats.stackexchange.com/questions/129179/why-is-glmnet-ridge-regression-giving-me-a-different-answer-than-manual-calculat
  hat_mat <- Xnew %*% solve(crossprod(X) + lambda * diag(c(0, rep(1, p)))) %*% t(X)
  fit <- as.vector(hat_mat %*% Y)

  return(fit)
}


#' Fit post-lasso regression model
#'
#' Estimates cross-validated post-lasso regression using the \code{plasso} package.
#' Performs lasso selection followed by OLS on selected variables.
#'
#' @param X Numeric matrix of covariates for training sample.
#' @param Y Numeric vector of outcomes for training sample.
#' @param arguments List of arguments passed to \code{cv.plasso}.
#'
#' @return An object of class \code{plasso_fit} containing the fitted
#'   \code{cv.plasso} model.
#'
#' @keywords internal
plasso_fit <- function(X, Y, arguments = list()) {
  if (!requireNamespace("plasso", quietly = TRUE)) {
    stop("The 'plasso' package is not installed. This learner is optional (in Suggests).\n",
      "Install it with: devtools::install_github('stefan-1997/plasso')",
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
#' @param object An object of class \code{plasso_fit} from \code{\link{plasso_fit}}.
#' @param Xnew Numeric matrix of covariates for test sample. If \code{NULL},
#'             predictions are generated for the training data.
#' @param ... Ignored additional arguments.
#'
#' @return Numeric vector of predicted values.
#'
#' @method predict plasso_fit
#' @keywords internal
#' @exportS3Method
predict.plasso_fit <- function(object, Xnew = NULL, ...) {
  if (is.null(Xnew)) Xnew <- object$x

  fit <- as.vector(predict(object, newx = Xnew, type = "response", s = "optimal", se_rule = 0)$plasso)

  return(fit)
}


#' Fit post-lasso regression using hdm package
#'
#' Estimates post-lasso regression under homoscedastic and heteroscedastic
#' non-Gaussian disturbances using the \code{hdm} package.
#'
#' @param X Numeric matrix of covariates for training sample.
#' @param Y Numeric vector of outcomes for training sample.
#' @param arguments List of arguments passed to \code{rlasso}.
#'
#' @return An object of class \code{rlasso_fit} containing the fitted
#'   \code{rlasso} model.
#'
#' @keywords internal
rlasso_fit <- function(X, Y, arguments = list()) {
  if (!requireNamespace("hdm", quietly = TRUE)) {
    stop("The 'hdm' package is not installed. 'rlasso' learner is optional (in Suggests).\n",
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
#' @param object An object of class \code{rlasso_fit} from \code{\link{rlasso_fit}}.
#' @param Xnew Numeric matrix of covariates for test sample. If \code{NULL},
#'             predictions are generated for the training data.
#' @param ... Ignored additional arguments.
#'
#' @return Numeric vector of predicted values.
#'
#' @method predict rlasso_fit
#' @keywords internal
#' @exportS3Method
predict.rlasso_fit <- function(object, Xnew = NULL, ...) {
  if (is.null(Xnew)) Xnew <- object$model

  fit <- as.vector(predict(object, newdata = Xnew))

  return(fit)
}


#' Fit random forest using grf package
#'
#' Fits a random forest for regression using the \code{grf} package.
#'
#' @param X Numeric matrix of covariates for training sample.
#' @param Y Numeric vector of outcomes for training sample.
#' @param arguments List of arguments passed to \code{regression_forest}.
#'
#' @return An object of class \code{forest_grf_fit} containing the fitted
#'   \code{regression_forest} model.
#'
#' @keywords internal
forest_grf_fit <- function(X, Y, arguments = list()) {
  if (!requireNamespace("grf", quietly = TRUE)) {
    stop("The 'grf' package is not installed. 'regression_forest' learner is optional (in Suggests).\n",
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
#' @param object An object of class \code{forest_grf_fit} from
#'   \code{\link{forest_grf_fit}}.
#' @param Xnew Numeric matrix of covariates for test sample. If \code{NULL},
#'             predictions are generated for the training data.
#' @param ... Ignored additional arguments.
#'
#' @return Numeric vector of predicted values.
#'
#' @method predict forest_grf_fit
#' @keywords internal
#' @exportS3Method
predict.forest_grf_fit <- function(object, Xnew = NULL, ...) {
  if (is.null(Xnew)) Xnew <- object$X.orig

  fit <- predict(object, newdata = Xnew)$prediction

  return(fit)
}


#' Fit XGBoost model
#'
#' Fits an XGBoost model using the \code{xgboost} package.
#'
#' @param X Numeric matrix of covariates for training sample.
#' @param Y Numeric vector of outcomes for training sample.
#' @param arguments List of arguments passed to \code{xgb.train}.
#'
#' @details
#' For successful smoother extraction in \pkg{OutcomeWeights}, certain parameter
#' restrictions are recommended. If not specified by the user, the following
#' values are suggested:
#' \itemize{
#'   \item \code{max_delta_step = 0}
#'   \item \code{base_score = mean(Y)} 
#'   \item \code{subsample = 1} 
#'   \item \code{alpha = 0}
#' }
#' User-specified values in \code{arguments} take precedence over suggestions.
#'
#' @return An object of class \code{xgboost_fit} containing the fitted
#'   \code{xgb.Booster} model.
#'
#' @keywords internal
#' 
xgboost_fit <- function(X, Y, arguments = list()) {
  if (!requireNamespace("xgboost", quietly = TRUE)) {
    stop("The 'xgboost' package is not installed. This learner is optional (in Suggests).\n",
      "Install it with: install.packages('xgboost')",
      call. = FALSE
    )
  }
  
  # Set the defaults (for potential smoother to be extracted)
  if (is.null(arguments$max_delta_step)) arguments$max_delta_step <- 0
  if (is.null(arguments$subsample)) arguments$subsample <- 1
  if (is.null(arguments$alpha)) arguments$alpha <- 0
  
  nrounds <- if (is.null(arguments$nrounds)) 100 else arguments$nrounds
  arguments$nrounds <- NULL

  dtrain <- xgboost::xgb.DMatrix(data = as.matrix(X), label = Y)
  
  # [Smoother] Fix for floating-point precision mismatch (XGBoost -> float32, mean(Y) -> float64)
  if (is.null(arguments$base_score)) arguments$base_score <- mean(xgboost::getinfo(dtrain, "label"))
  
  booster <- do.call(xgboost::xgb.train, c(list(data = dtrain, nrounds = nrounds, params = arguments)))

  xgb <- list(
    booster = booster,
    X = X,
    Y = Y
  )

  class(xgb) <- c("xgboost_fit")
  return(xgb)
}


#' Predict method for XGBoost models
#'
#' Generates predictions from a fitted XGBoost model.
#'
#' @param object An object of class \code{xgboost_fit} from \code{\link{xgboost_fit}}.
#' @param Xnew Numeric matrix of covariates for test sample. If \code{NULL},
#'             predictions are generated for the training data.
#' @param ... Ignored additional arguments.
#'
#' @return Numeric vector of predicted values.
#'
#' @method predict xgboost_fit
#' @keywords internal
#' @exportS3Method
predict.xgboost_fit <- function(object, Xnew = NULL, ...) {
  if (is.null(Xnew)) Xnew <- object$X

  xgb <- object$booster
  dtest <- xgboost::xgb.DMatrix(data = as.matrix(Xnew))
  
  fit <- predict(object = xgb, newdata = dtest)

  return(fit)
}


#' Fit lasso regression model
#'
#' Estimates lasso regression using the \code{glmnet} package.
#' Selects optimal lambda penalty parameter via cross-validation.
#'
#' @param X Numeric matrix of covariates for training sample.
#' @param Y Numeric vector of outcomes for training sample.
#' @param arguments List of arguments passed to \code{cv.glmnet}.
#'
#' @return An object of class \code{lasso_fit} containing the fitted
#'   \code{cv.glmnet} model.
#'
#' @keywords internal
lasso_fit <- function(X, Y, arguments = list()) {
  if (!requireNamespace("glmnet", quietly = TRUE)) {
    stop("The 'glmnet' package is not installed. 'Lasso' learner is optional (in Suggests).\n",
      "Install it with: install.packages('glmnet')",
      call. = FALSE
    )
  }

  lasso <- do.call(glmnet::cv.glmnet, c(list(x = X, y = Y), arguments))

  class(lasso) <- c(class(lasso), "lasso_fit")
  return(lasso)
}


#' Predict method for lasso regression models
#'
#' Generates predictions from a fitted lasso regression model
#' using the optimal lambda selected during cross-validation.
#'
#' @param object An object of class \code{lasso_fit} from \code{\link{lasso_fit}}.
#' @param Xnew Numeric matrix of covariates for test sample. If \code{NULL},
#'   predictions are generated for the training data.
#' @param ... Ignored additional arguments.
#'
#' @return Numeric vector of predicted values.
#'
#' @method predict lasso_fit
#' @keywords internal
#' @exportS3Method
predict.lasso_fit <- function(object, Xnew = NULL, ...) {
  if (is.null(Xnew)) Xnew <- object[["call"]][["x"]]

  fit <- as.vector(predict(object, newx = Xnew, type = "response", s = "lambda.min"))

  return(fit)
}


#' Fit distributional random forest model
#'
#' Fits a distributional random forest using the \code{drf} package with
#' Fourier MMD splitting rule.
#'
#' @param X Numeric matrix of covariates for training sample.
#' @param Y Numeric vector of outcomes for training sample.
#' @param arguments List of arguments passed to \code{drf}.
#'
#' @return An object of class \code{forest_drf_fit} containing the fitted
#'   \code{drf} model.
#'
#' @keywords internal
forest_drf_fit <- function(X, Y, arguments = list()) {
  if (!requireNamespace("drf", quietly = TRUE)) {
    stop("The 'drf' package is not installed. This learner is optional (in Suggests).\n",
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
#' @param object An object of class \code{forest_drf_fit} from
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
#' @exportS3Method
predict.forest_drf_fit <- function(object, Xnew = NULL, functional = "mean", ...) {
  if (is.null(Xnew)) Xnew <- object$X.orig

  fit <- as.vector(predict(object, newdata = Xnew, functional = functional))

  return(fit)
}


#' Fit Random Forest for continuous outcomes using \code{ranger}
#'
#' Fits a Random Forest model using the \code{ranger} package for continuous outcomes.
#'
#' @param X Covariate matrix of training sample.
#' @param Y Vector of continuous outcomes of training sample.
#' @param arguments List of arguments passed to \code{ranger::ranger()}.
#'
#' @details
#' For successful smoother extraction in the \pkg{OutcomeWeights} package, 
#' \code{keep.inbag = TRUE} is set.
#'
#' @return An object of class \code{ranger} containing the fitted model.
#'
#' @keywords internal
ranger_fit <- function(X, Y, arguments = list()) {
  if (!requireNamespace("ranger", quietly = TRUE)) {
    stop("The 'ranger' package is not installed. Random Forest is optional (in Suggests).\n",
         "Install it with: install.packages('ranger')",
         call. = FALSE
    )
  }

  model <- do.call(ranger::ranger, c(list(x = X, y = Y, keep.inbag = TRUE), arguments))
  
  return(model)
}


#' Predict method for Random Forest continuous outcome fits
#'
#' Generates predicted values from fitted Random Forest models for continuous outcomes.
#'
#' @param object Output of \code{\link{ranger_fit}}.
#' @param Xnew Covariate matrix of test sample. If \code{NULL}, uses training data.
#' @param ... Ignored additional arguments.
#'
#' @return A numeric vector of predicted continuous outcomes.
#'
#' @method predict ranger_fit
#' @keywords internal
#' @exportS3Method
predict.ranger_fit <- function(object, Xnew = NULL, ...) {
  if (is.null(Xnew)) Xnew <- object$call$x
  
  fit <- predict(object, data = Xnew)$predictions
  
  return(fit)
}


#' Pseudo-fitting for k-Nearest Neighbors
#'
#' Creates a k-NN model object without actual fitting, as k-NN fitting occurs
#' during prediction.
#'
#' @param X Numeric matrix of covariates for training sample.
#' @param Y Numeric vector of outcomes for training sample.
#' @param arguments List of arguments to be used in k-NN prediction. Typically
#'   includes \code{k} (number of neighbors).
#'
#' @return A list of class \code{knn_fit} containing:
#'   \item{X}{Training covariate matrix}
#'   \item{Y}{Training outcome vector}
#'   \item{arguments}{List of arguments for k-NN prediction}
#'
#' @keywords internal
knn_fit <- function(X, Y, arguments = list()) {
  knn <- list(
    X = X,
    Y = Y,
    arguments = arguments
  )
  class(knn) <- "knn_fit"
  return(knn)
}


#' Predict method for k-Nearest Neighbors
#'
#' Generates predictions using the k-Nearest Neighbors algorithm with internal
#' Euclidean distance calculation, following the same approach as the
#' \code{Distance_for_KNN_test()} function of the \code{FastKNN} package.
#'
#' @param object Output of \code{\link{knn_fit}} of class \code{knn_fit}.
#'               List of arguments to be used in k-NN prediction. Typically
#'               includes \code{k}. (number of neighbors). If no \code{k} is
#'               specified, defaults to \code{k = 10}.
#' @param Xnew Numeric matrix of covariates for test sample. If \code{NULL},
#'             predictions are generated for the training data.
#' @param ... Ignored additional arguments.
#'
#' @return Numeric vector of predicted values.
#'
#' @method predict knn_fit
#' @keywords internal
#' @exportS3Method
predict.knn_fit <- function(object, Xnew = NULL, ...) {
  X <- object$X
  Y <- object$Y

  if (is.null(Xnew)) Xnew <- X
  k <- if (is.null(object$arguments$k)) 10 else object$arguments$k

  # Credit: FastKNN package
  euclidean_dist <- t(apply(Xnew, 1, function(x_i) {
    sqrt(rowSums(sweep(X, 2, x_i)^2))
  }))

  get_binary_vec <- function(row) {
    min_indices <- order(row)[1:k]
    binary_vec <- rep(0, length(row))
    binary_vec[min_indices] <- 1
    return(binary_vec)
  }

  w <- apply(euclidean_dist, 1, get_binary_vec)
  w <- t(w) / k

  fit <- as.vector(w %*% Y)

  return(fit)
}


#' Fit (multinomial) logistic regression using \code{glmnet}
#'
#' Uses \code{family = "binomial"} for binary outcomes and \code{family = "multinomial"}
#' for multiclass outcomes.
#'
#' @param X Covariate matrix of training sample
#' @param Y Vector of outcomes of training sample
#' @param arguments List of arguments passed to \code{glmnet::glmnet()}
#'
#' @return An object of class \code{glmnet} containing the fitted model.
#'
#' @keywords internal
logit_fit <- function(X, Y, arguments = list()) {
  if (!requireNamespace("glmnet", quietly = TRUE)) {
    stop("The 'glmnet' package is not installed. Logistic regression is optional (in Suggests).\n",
      "Install it with: install.packages('glmnet')",
      call. = FALSE
    )
  }

  if (length(unique(Y)) == 2) {
    logit <- do.call(glmnet::glmnet, c(list(
      x = X, y = Y,
      family = "binomial",
      type.measure = "class"
    ), arguments))
  } else {
    logit <- do.call(glmnet::glmnet, c(list(
      x = X, y = Y,
      family = "multinomial",
      type.measure = "class"
    ), arguments))
  }

  return(logit)
}


#' Predict method for \code{glmnet} logistic regression fits
#'
#' Generates predictions from fitted logistic regression models.
#'
#' @param object Output of \code{\link{logit_fit}}.
#' @param Xnew Covariate matrix of test sample. If \code{NULL}, uses training data.
#' @param ... Ignored additional arguments.
#'
#' @return A data frame of predicted class probabilities.
#'   \describe{
#'     \item{Binary outcome}{A two-column data frame with probabilities for each class.}
#'     \item{Multiclass outcome}{A data frame with one column per class, containing the predicted probabilities.}
#'   }
#'
#' @keywords internal
predict.logit_fit <- function(object, Xnew = NULL, ...) {
  if (is.null(Xnew)) Xnew <- object[["call"]][["x"]]
  lambda <- min(object$lambda)

  fit <- predict(object, newx = as.matrix(Xnew), s = lambda, type = "response")

  if (ncol(fit) == 1) fit <- as.vector(fit)

  return(fit)
}


#' Fit binary logistic or probit regression using \code{glm}
#'
#' Fits binary regression models using \code{glm}. Supports both logistic
#'  (logit) and probit link functions. For binary outcomes only.
#'
#' @param X Covariate matrix of training sample.
#' @param Y Vector of binary outcomes.
#' @param arguments List of arguments passed to \code{stats::glm}. To use probit,
#'   specify \code{family = binomial(link = "probit")}. Default is logit.
#'
#' @return An object of class \code{glm} containing the fitted model.
#'
#' @keywords internal
glm_fit <- function(X, Y, arguments = list()) {
  if (length(unique(Y)) != 2) {
    stop(
      "glm_fit only supports binary outcomes (exactly 2 levels). Found: ",
      length(unique(Y)), " levels"
    )
  }

  data <- data.frame(Y = as.factor(Y), X)

  if (!"family" %in% names(arguments)) {
    arguments$family <- stats::binomial(link = "logit")
  }

  model <- do.call(stats::glm, c(list(formula = Y ~ ., data = data), arguments))

  return(model)
}


#' Predict method for binary GLM fits
#'
#' Generates predicted class probabilities from fitted logistic or probit regression models.
#'
#' @param object Output of \code{\link{glm_fit}}.
#' @param Xnew Covariate matrix of test sample. If \code{NULL}, uses training data.
#' @param ... Ignored additional arguments.
#'
#' @return A numeric vector of probabilities for the positive class.
#'
#' @keywords internal
predict.glm_fit <- function(object, Xnew = NULL, ...) {
  if (is.null(Xnew)) Xnew <- as.matrix(object$data[, -1, drop = FALSE])

  fit <- predict(object, newdata = as.data.frame(Xnew), type = "response")

  return(fit)
}


#' Fit (multinomial) logistic regression using \code{nnet}
#'
#' Fits a logistic regression model using the \code{nnet} package.
#' Uses \code{multinom()} for binary and multiclass outcomes.
#'
#' @param X Covariate matrix of training sample.
#' @param Y Vector of outcomes of training sample.
#' @param arguments List of arguments passed to \code{nnet::multinom()}.
#'
#' @return An object of class \code{nnet} containing the fitted model.
#'
#' @keywords internal
logit_nnet_fit <- function(X, Y, arguments = list()) {
  if (!requireNamespace("nnet", quietly = TRUE)) {
    stop("The 'nnet' package is not installed. Multinomial logit is optional (in Suggests).\n",
      "Install it with: install.packages('nnet')",
      call. = FALSE
    )
  }

  data <- data.frame(Y = as.factor(Y), X)
  arguments <- utils::modifyList(list(trace = FALSE), arguments)

  model <- do.call(nnet::multinom, c(list(data = data, formula = Y ~ ., model = TRUE), arguments))

  return(model)
}


#' Predict method for \code{nnet} logistic regression fits
#'
#' Generates predicted class probabilities from fitted logistic regression models.
#'
#' @param object Output of \code{\link{logit_nnet_fit}}.
#' @param Xnew Covariate matrix of test sample. If \code{NULL}, uses training data.
#' @param ... Ignored additional arguments.
#'
#' @return Predicted probabilities:
#'   \describe{
#'     \item{Binary outcome}{A numeric vector of probabilities for the positive class.}
#'     \item{Multiclass outcome}{A numeric matrix with one column per class.}
#'   }
#'
#' @keywords internal
predict.logit_nnet_fit <- function(object, Xnew = NULL, ...) {
  if (is.null(Xnew)) Xnew <- as.matrix(object$model[, -1, drop = FALSE])

  fit <- predict(object, newdata = Xnew, type = "probs")

  return(fit)
}


#' Fit Gaussian Naive Bayes using \code{naivebayes}
#'
#' Fits a Gaussian Naive Bayes model using the \code{naivebayes} package.
#'
#' @param X Covariate matrix of training sample.
#' @param Y Vector of outcomes of training sample.
#' @param arguments List of arguments passed to \code{naivebayes::gaussian_naive_bayes()}.
#'
#' @return An object of class \code{gaussian_naive_bayes} containing the fitted model.
#'
#' @keywords internal
nb_gaussian_fit <- function(X, Y, arguments = list()) {
  if (!requireNamespace("naivebayes", quietly = TRUE)) {
    stop("The 'naivebayes' package is not installed. Gaussian Naive Bayes is optional (in Suggests).\n",
      "Install it with: install.packages('naivebayes')",
      call. = FALSE
    )
  }

  Y <- as.factor(Y)

  model <- do.call(naivebayes::gaussian_naive_bayes, c(list(x = X, y = Y), arguments))

  return(model)
}


#' Predict method for Gaussian Naive Bayes fits
#'
#' Generates predicted class probabilities from fitted Gaussian Naive Bayes models.
#'
#' @param object Output of \code{\link{nb_gaussian_fit}}.
#' @param Xnew Covariate matrix of test sample. If \code{NULL}, uses training data.
#' @param ... Ignored additional arguments.
#'
#' @return Predicted probabilities:
#'   \describe{
#'     \item{Binary outcome}{A numeric vector of probabilities for the positive class.}
#'     \item{Multiclass outcome}{A numeric matrix with one column per class.}
#'   }
#'
#' @keywords internal
predict.nb_gaussian_fit <- function(object, Xnew = NULL, ...) {
  if (is.null(Xnew)) Xnew <- object[["call"]][["x"]]

  fit <- predict(object, newdata = Xnew, type = "prob")

  if (ncol(fit) == 2) fit <- fit[, 2, drop = TRUE]

  return(fit)
}


#' Fit Bernoulli Naive Bayes using \code{naivebayes}
#'
#' Fits a Bernoulli Naive Bayes model using the \code{naivebayes} package.
#'
#' @param X Covariate matrix of training sample.
#' @param Y Vector of outcomes of training sample.
#' @param arguments List of arguments passed to \code{naivebayes::naive_bayes()}.
#'
#' @return An object of class \code{bernoulli_naive_bayes} containing the fitted model.
#'
#' @keywords internal
nb_bernoulli_fit <- function(X, Y, arguments = list()) {
  if (!requireNamespace("naivebayes", quietly = TRUE)) {
    stop("The 'naivebayes' package is not installed. Bernoulli Naive Bayes is optional (in Suggests).\n",
      "Install it with: install.packages('naivebayes')",
      call. = FALSE
    )
  }

  Y <- as.factor(Y)

  model <- do.call(naivebayes::naive_bayes, c(list(x = X, y = Y), arguments))

  return(model)
}


#' Predict method for Bernoulli Naive Bayes fits
#'
#' Generates predicted class probabilities from fitted Bernoulli Naive Bayes models.
#'
#' @param object Output of \code{\link{nb_bernoulli_fit}}.
#' @param Xnew Covariate matrix of test sample. If \code{NULL}, uses training data.
#' @param ... Ignored additional arguments.
#'
#' @return Predicted probabilities:
#'   \describe{
#'     \item{Binary outcome}{A numeric vector of probabilities for the positive class.}
#'     \item{Multiclass outcome}{A numeric matrix with one column per class.}
#'   }
#'
#' @keywords internal
predict.nb_bernoulli_fit <- function(object, Xnew = NULL, ...) {
  if (is.null(Xnew)) Xnew <- object[["call"]][["x"]]

  fit <- predict(object, newdata = Xnew, type = "prob")

  if (ncol(fit) == 2) fit <- fit[, 2, drop = TRUE]

  return(fit)
}


#' Fit gradient boosting model using \code{xgboost}
#'
#' Fits a gradient boosting model using the \code{xgboost} package.
#' Uses \code{objective = "binary:logistic"} for binary outcomes and
#' \code{objective = "multi:softprob"} for multiclass outcomes.
#'
#' @param X Covariate matrix of training sample.
#' @param Y Vector of outcomes of training sample.
#' @param arguments List of arguments passed to \code{xgboost::xgboost()}.
#'
#' @return An object of class \code{xgb.Booster} containing the fitted model.
#'
#' @keywords internal
xgboost_prop_fit <- function(X, Y, arguments = list()) {
  if (!requireNamespace("xgboost", quietly = TRUE)) {
    stop("The 'xgboost' package is not installed. This learner is optional (in Suggests).\n",
      "Install it with: install.packages('xgboost')",
      call. = FALSE
    )
  }

  K <- length(unique(Y))
  dtrain <- xgboost::xgb.DMatrix(data = X, label = Y)
  
  nrounds <- if (is.null(arguments$nrounds)) 100 else arguments$nrounds
  arguments$nrounds <- NULL
  
  if (K == 2) {
    # Binary classification
    params <- list(
      objective = "binary:logistic",
      eval_metric = "logloss"
    )
  } else {
    if (min(Y) == 1) Y <- Y - 1
    # Multi-class classification
    params <- list(
      objective = "multi:softprob",
      num_class = K,
      eval_metric = "mlogloss"
    )
  }
  
  # Defaults override user arguments
  params <- utils::modifyList(params, arguments)
  
  booster <- xgboost::xgb.train(
    params = params,
    data = dtrain,
    nrounds = nrounds
  )
  
  model <- list(
    booster = booster,
    X = X,
    Y = Y
  )

  return(model)
}


#' Predict method for \code{xgboost} fits
#'
#' Generates predicted class probabilities from fitted gradient boosting models.
#'
#' @param object Output of \code{\link{xgboost_prop_fit}}.
#' @param Xnew Covariate matrix of test sample. If \code{NULL}, uses training data.
#' @param ... Ignored additional arguments.
#'
#' @return Predicted probabilities:
#'   \describe{
#'     \item{Binary outcome}{A numeric vector of probabilities for the positive class.}
#'     \item{Multiclass outcome}{A numeric matrix with one column per class.}
#'   }
#'
#' @keywords internal
predict.xgboost_prop_fit <- function(object, Xnew = NULL, ...) {
  if (is.null(Xnew)) Xnew <- object$X
  
  Y <- object$Y
  xgb <- object$booster

  fit <- predict(object = xgb, newdata = Xnew, type = "prob")

  if (length(unique(Y)) > 2) {
    fit <- matrix(fit, nrow = nrow(Xnew), ncol = length(unique(Y)), byrow = TRUE)
    colnames(fit) <- sort(unique(Y))
  }

  return(fit)
}


#' Fit Support Vector Machine using \code{e1071}
#'
#' Fits a Support Vector Machine (SVM) model using the \code{e1071} package.
#' Uses a radial kernel and \code{C-classification} by default, with probability estimates enabled.
#'
#' @param X Covariate matrix of training sample.
#' @param Y Vector of outcomes of training sample.
#' @param arguments List of arguments passed to \code{e1071::svm()}.
#'
#' @return An object of class \code{svm} containing the fitted model.
#'
#' @keywords internal
svm_fit <- function(X, Y, arguments = list()) {
  if (!requireNamespace("e1071", quietly = TRUE)) {
    stop("The 'e1071' package is not installed. SVM learner is optional (in Suggests).\n",
      "Install it with: install.packages('e1071')",
      call. = FALSE
    )
  }

  model <- do.call(e1071::svm, c(list(
    y = Y, x = X, probability = TRUE, kernel = "radial", type = "C-classification"
  ), arguments))

  return(model)
}


#' Predict method for Support Vector Machine fits
#'
#' Generates predicted class probabilities from fitted Support Vector Machine models.
#'
#' @param object Output of \code{\link{svm_fit}}.
#' @param Xnew Covariate matrix of test sample. If \code{NULL}, uses training data.
#' @param ... Ignored additional arguments.
#'
#' @return Predicted probabilities:
#'   \describe{
#'     \item{Binary outcome}{A numeric vector of probabilities for the positive class.}
#'     \item{Multiclass outcome}{A numeric matrix with one column per class.}
#'   }
#'
#' @keywords internal
predict.svm_fit <- function(object, Xnew = NULL, ...) {
  if (is.null(Xnew)) Xnew <- object[["call"]][["x"]]

  pred <- predict(object, newdata = Xnew, probability = TRUE)
  fit <- attr(pred, "probabilities")

  if (ncol(fit) == 2) fit <- fit[, 2, drop = TRUE]

  return(fit)
}


#' Fit Probability Forest using \code{grf}
#'
#' Fits a Probability Forest model using the \code{grf} package.
#'
#' @param X Covariate matrix of training sample.
#' @param Y Vector of outcomes of training sample.
#' @param arguments List of arguments passed to \code{grf::probability_forest()}.
#'
#' @return An object of class \code{probability_forest} containing the fitted model.
#'
#' @keywords internal
prob_forest_fit <- function(X, Y, arguments = list()) {
  if (!requireNamespace("grf", quietly = TRUE)) {
    stop("The 'grf' package is not installed. Probability Forest is optional (in Suggests).\n",
      "Install it with: install.packages('grf')",
      call. = FALSE
    )
  }

  model <- do.call(grf::probability_forest, c(list(X = X, Y = as.factor(Y)), arguments))

  return(model)
}


#' Predict method for Probability Forest fits
#'
#' Generates predicted class probabilities from fitted Probability Forest models.
#'
#' @param object Output of \code{\link{prob_forest_fit}}.
#' @param Xnew Covariate matrix of test sample. If \code{NULL}, uses training data.
#' @param ... Ignored additional arguments.
#'
#' @return Predicted probabilities:
#'   \describe{
#'     \item{Binary outcome}{A numeric vector of probabilities for the positive class.}
#'     \item{Multiclass outcome}{A numeric matrix with one column per class.}
#'   }
#'
#' @keywords internal
predict.prob_forest_fit <- function(object, Xnew = NULL, ...) {
  if (is.null(Xnew)) Xnew <- object[["X.orig"]]

  fit <- predict(object, newdata = Xnew)$predictions

  if (ncol(fit) == 2) fit <- fit[, 2, drop = TRUE]

  return(fit)
}


#' Fit k-Nearest Neighbor using \code{kknn}
#'
#' Fits a k-Nearest Neighbor (kNN) model using the \code{kknn} package.
#'
#' @param X Covariate matrix of training sample.
#' @param Y Vector of outcomes of training sample.
#' @param arguments List of arguments passed to \code{kknn::train.kknn()}.
#'
#' @return An object of class \code{kknn} containing the fitted model.
#'
#' @keywords internal
knn_prop_fit <- function(X, Y, arguments = list()) {
  if (!requireNamespace("kknn", quietly = TRUE)) {
    stop("The 'kknn' package is not installed. k-Nearest Neighbor is optional (in Suggests).\n",
      "Install it with: install.packages('kknn')",
      call. = FALSE
    )
  }

  data <- data.frame(Y = as.factor(Y), X)

  model <- do.call(kknn::train.kknn, c(list(formula = Y ~ ., data = data), arguments))

  return(model)
}


#' Predict method for k-Nearest Neighbor fits
#'
#' Generates predicted class probabilities from fitted k-Nearest Neighbor models.
#'
#' @param object Output of \code{\link{knn_prop_fit}}.
#' @param Xnew Covariate matrix of test sample. If \code{NULL}, uses training data.
#' @param ... Ignored additional arguments.
#'
#' @return Predicted probabilities:
#'   \describe{
#'     \item{Binary outcome}{A numeric vector of probabilities for the positive class.}
#'     \item{Multiclass outcome}{A numeric matrix with one column per class.}
#'   }
#'
#' @keywords internal
predict.knn_prop_fit <- function(object, Xnew = NULL, ...) {
  if (is.null(Xnew)) Xnew <- object$data[-1]

  fit <- predict(object, newdata = as.data.frame(Xnew), type = "prob")

  if (ncol(fit) == 2) fit <- fit[, 2, drop = TRUE]

  return(fit)
}


#' Fit Probability Forest for propensity scores using \code{ranger}
#'
#' Fits a Probability Forest model using the \code{ranger} package for propensity score estimation.
#'
#' @param X Covariate matrix of training sample.
#' @param Y Vector of outcomes of training sample.
#' @param arguments List of arguments passed to \code{ranger::ranger()}.
#'
#' @return An object of class \code{ranger} containing the fitted model.
#'
#' @keywords internal
ranger_prop_fit <- function(X, Y, arguments = list()) {
  if (!requireNamespace("ranger", quietly = TRUE)) {
    stop("The 'ranger' package is not installed. Probability Forest is optional (in Suggests).\n",
         "Install it with: install.packages('ranger')",
         call. = FALSE
    )
  }
  
  data <- data.frame(Y = as.factor(Y), X)
  
  model <- do.call(ranger::ranger, c(list(data = data, formula = Y ~ ., probability = TRUE), arguments))
  
  return(model)
}


#' Predict method for Probability Forest propensity score fits
#'
#' Generates predicted class probabilities from fitted Probability Forest models for propensity scores.
#'
#' @param object Output of \code{\link{ranger_prop_fit}}.
#' @param Xnew Covariate matrix of test sample. If \code{NULL}, uses training data.
#' @param ... Ignored additional arguments.
#'
#' @return Predicted probabilities:
#'   \describe{
#'     \item{Binary outcome}{A numeric vector of probabilities for the positive class.}
#'     \item{Multiclass outcome}{A numeric matrix with one column per class.}
#'   }
#'
#' @keywords internal
predict.ranger_prop_fit <- function(object, Xnew = NULL, ...) {
  if (is.null(Xnew)) Xnew <- object$call$data
  
  data <- as.data.frame(Xnew)
  data$Y <- as.factor(0)
  
  fit <- predict(object, data = data)$predictions
  
  if (ncol(fit) == 2) fit <- fit[, 2, drop = TRUE]
  
  return(fit)
}
