#' Smoother weights for outcome models
#'
#' Extract smoother matrices for nuisance parameter outcome models from a
#' \code{NuisanceParameters} object created by \code{\link{nuisance_parameters}}.
#'
#' @param object A \code{NuisanceParameters} object created by
#'   \code{\link{nuisance_parameters}}, or a path to an RDS file containing one.
#' @param NuPa Character vector specifying which nuisance parameters to extract 
#'  smoothers for. Available outcome nuisance parameters 
#'  are \code{c("Y.hat", "Y.hat.d", "Y.hat.z")}.
#' @param quiet If \code{FALSE}, the method currently computed is printed to the console.
#' @param subset Logical vector indicating which observations to use when determining
#'   ensemble weights. If not provided, all observations are used.
#' @param ... Additional arguments.
#'
#' @return A list of matrices of dimension \code{nrow(X)} × \code{nrow(X)} containing ensemble smoother
#'   matrices
#'
#' @export
get_outcome_weights <- function(object,
                                NuPa,
                                quiet = TRUE,
                                subset = NULL,
                                ...) {
  # Input check
  NuPa <- match.arg(NuPa, choices = c("Y.hat", "Y.hat.d", "Y.hat.z"), several.ok = TRUE)
  
  if (inherits(object, "NuisanceParameters")) {
    object <- object
  } else if (is.character(object) && length(object) == 1) {
    loaded <- readRDS(object)
    if (!is.list(loaded) || is.null(loaded$models)) {
      stop("Loaded object is not a valid NuisanceParameters list.")
    }
    object <- list(models = loaded$models, numbers = loaded$numbers)
  } else {
    stop("np_object must be a NuisanceParameters object or a valid path to an RDS file.")
  }
  
  X <- object$numbers$X
  Y <- object$numbers$Y
  cv <- object$numbers$cv
  d_mat <- object$numbers$d_mat
  z_mat <- object$numbers$z_mat
  cf_mat <- object$numbers$cf_mat
  methods <- object$numbers$methods
  models <- object$models

  # Lasso check
  if (any(sapply(methods, function(met) {
    any(sapply(met, function(x) identical(x$method, "lasso")))
  }))) {stop("Smoother matrices cannot be extracted if method = 'lasso' was used.")}

  # Models check
  m_names <- paste0(NuPa, "_m")
  outcome_m <- m_names[m_names %in% names(models) & !vapply(models[m_names], is.character, logical(1))]
  missing_m <- setdiff(m_names, outcome_m)

  if (length(missing_m) > 0) {
    message(
      "Missing models: ", paste(missing_m, collapse = ", "),
      "\nContinuing with outcome models: ", paste(outcome_m, collapse = ", "))
  }
  if (length(outcome_m) == 0) stop("No valid (ensemble) models found in object$models.")
  if (is.null(subset)) subset <- rep(TRUE, length(Y))

  results <- list()

  # Process each model
  for (m_name in outcome_m) {
    current_model <- models[[m_name]]; nupa <- sub("_m$", "", m_name)
    nupa_smoother <- list()

    ## Short stacking
    if (cv == 1) {
      # Case 1: Single ensemble model (i.e. Y.hat)
      if ("ens_object" %in% names(current_model)) {
        nupa_smoother <- get_smoother(
          object = current_model,
          methods = methods[[nupa]],
          X = X,
          Y = Y,
          cv = cv,
          cf_mat = cf_mat,
          quiet = quiet
        )
      }
      # Case 2: List of models (i.e. Y.hat.d, Y.hat.z)
      else {
        for (i in seq_along(current_model)) {
          sub_element <- current_model[[i]]
          subset <- get_subset(sub_element, models, d_mat, z_mat)
          
          nupa_smoother[[i]] <- get_smoother(
            object = sub_element,
            methods = methods[[nupa]],
            X = X,
            Y = Y,
            subset = subset,
            cv = cv,
            cf_mat = cf_mat,
            quiet = quiet
            )
        }
      }
    }

    ## Standard stacking
    else {
      # Case 1: List of models
      if (!inherits(current_model, "ens.learner")) {
        for (i in seq_along(current_model)) {
          sub_element <- current_model[[i]]
          
          if (is.list(sub_element)) {
            subset <- get_subset(sub_element, models, d_mat, z_mat)
            
            nupa_smoother[[i]] <- get_smoother(
              object = sub_element,
              methods = methods[[nupa]],
              X = X,
              Y = Y,
              subset = subset,
              cv = cv,
              cf_mat = cf_mat,
              quiet = quiet
            )
          }
        }
      }
      # Case 2: Single ensemble model
      else {
        nupa_smoother <- get_smoother(
          object = current_model,
          methods = methods[[nupa]],
          X = X,
          Y = Y,
          subset = NULL,
          cv = cv,
          cf_mat = cf_mat,
          quiet = quiet
        )
      }
    }
    results[[nupa]] <- nupa_smoother
  }

  return(results)
}


#' Smoother weights from ensemble learner
#'
#' Extract smoother weights from a short- or standard-stacked ensemble model
#' created by \code{\link{ensemble_short}} or \code{\link{ensemble}}.
#'
#' @param object Stacked ensemble learner from \code{\link{ensemble_short}} or \code{\link{ensemble}}.
#' @param methods List of methods built via \code{\link{create_method}}.
#' @param X Covariate matrix of training sample.
#' @param Y Vector of outcomes of training sample.
#' @param subset Logical vector indicating which observations to use when determining
#'   ensemble weights. If not provided, all observations are used.
#' @param cf_mat Logical matrix with \code{k} columns of indicators representing folds,
#'   typically created by \code{\link{prep_cf_mat}}.
#' @param cv Integer specifying the number of cross-validation folds used in estimating
#'   the ensemble model.
#' @param quiet If \code{FALSE}, the method currently computed is printed to the console.
#'
#' @return A matrix of dimension \code{nrow(X)} × \code{nrow(X)} containing ensemble smoother weights.
#'
#' @keywords internal
get_smoother <- function(object,
                         methods,
                         X, Y,
                         subset = NULL,
                         cf_mat, cv,
                         quiet = TRUE) {
  if (is.null(subset)) subset <- rep(TRUE, length(Y))

  ## Short-Stacking
  if (cv == 1) {
    smoother_w <- weights(
      object = object$ens_object, methods = methods, X = X, Y = Y, subset = subset, 
      nnls_w = object$nnls_w, cf_mat = cf_mat, quiet = quiet
      )

    ## Standard-Stacking
  } else if (cv > 1) {
    smoother_w <- matrix(0, nrow = nrow(X), ncol = nrow(X))

    for (i in 1:ncol(cf_mat)) {
      fold <- cf_mat[, i]
      X_tr <- X[!fold & subset, ]
      Y_tr <- Y[!fold & subset]
      X_te <- X[fold, ]

      smoother_w_fold <- weights(
        object = object[[i]]$ens_object, methods = methods, nnls_w = object[[i]]$nnls_w,
        X = X_tr, Y = Y_tr, Xnew = X_te, quiet = quiet
        )
      
      smoother_w[fold, !fold & subset] <- smoother_w_fold
    }
  }

  return(smoother_w)
}


#' Smoother weights from short-stacked ensemble learner
#'
#' Extract smoother weights from a short-stacked ensemble model created by
#' \code{\link{ensemble_short}}.
#'
#' @param object Short-stacked ensemble learner from \code{\link{ensemble_short}}.
#' @param methods List of method models from \code{\link{create_method}} used 
#'  as input for \code{\link{ensemble_short}}.
#' @param X Covariate matrix of training sample.
#' @param Y Vector of outcomes of training sample.
#' @param subset Logical vector indicating which observations to use for 
#'  determining ensemble weights. If not provided, all observations are used.
#' @param nnls_w Ensemble weights to aggregate predictions from different learners. 
#'  Must be a vector of length \code{ncol(object$cf_preds)}. Optional.
#' @param cf_mat Logical matrix with \code{k} columns of indicators representing 
#'  folds (for example created by \code{\link{prep_cf_mat}}).
#' @param quiet If \code{FALSE}, the method currently computed is printed to the console.
#' @param ... Ignore unused arguments.
#'
#' @return A matrix of dimension \code{nrow(X)} × \code{nrow(X)} containing ensemble smoother weights.
#'
#' @method weights ensemble_short
#' @keywords internal
weights.ensemble_short <- function(object,
                                   methods,
                                   X, Y,
                                   cf_mat,
                                   subset = NULL,
                                   nnls_w = NULL,
                                   quiet = TRUE,
                                   ...) {
  
  if (is.null(nnls_w)) nnls_w <- rep(1 / ncol(object$cf_preds), ncol(object$cf_preds))
  if (is.null(subset)) subset <- rep(TRUE, length(Y))
  
  ens_models <- object$ens_models
  smoother_w <- matrix(0, nrow = nrow(X), ncol = nrow(X))
  
  ## Smoother extraction
  for (i in 1:ncol(cf_mat)) {
    fold <- cf_mat[, i]
    
    # Subset is responsible for T-learner logic
    test_idx <- which(fold)
    train_idx <- which(!fold)
    train_idx_sub <- which(!fold & subset)
    idx_mapping <- match(train_idx_sub, train_idx)
    
    X_tr <- X[train_idx_sub, ]
    Y_tr <- Y[train_idx_sub]
    X_te <- X[test_idx, ]
    
    # Get T-learner smoother weights (test × train&subset x methods)
    w_array <- weights.ensemble_core(
      object = ens_models[[i]], methods = methods, 
      X_tr = X_tr, Y_tr = Y_tr, X_te = X_te, quiet = quiet)
    
    # Pad (to test × full train x methods) and insert raw weights
    w_padded <- array(0, dim = c(length(test_idx), length(train_idx), length(methods)))
    w_padded[, idx_mapping, ] <- w_array
    
    smoother_w[fold, !fold] <- agg_array(a = w_padded, w = nnls_w)
  }
  
  return(smoother_w)
}


#' Smoother weights from cross-validated ensemble learner
#'
#' Extract ensemble smoother weights for a covariate matrix based on a
#' cross-validated ensemble model created by \code{\link{ensemble}}.
#'
#' @param object Trained ensemble object from \code{\link{ensemble}}.
#' @param methods List of methods built via \code{\link{create_method}}.
#' @param X Covariate matrix of training sample.
#' @param Y Vector of outcomes of training sample.
#' @param Xnew Covariate matrix of test sample.
#' @param nnls_w Ensemble weights to aggregate predictions from different learners. 
#'  Must be a vector of length \code{ncol(object$cf_preds)}. Optional.
#' @param quiet If \code{FALSE}, the method currently computed is printed to the console.
#' @param ... Ignore unused arguments.
#'
#' @return A matrix of dimension \code{nrow(Xnew)} × \code{nrow(X)} containing ensemble smoother weights.
#'
#' @method weights ensemble
#' @keywords internal
weights.ensemble <- function(object,
                             methods,
                             X, Y, Xnew,
                             nnls_w = NULL,
                             quiet = TRUE,
                             ...) {
  
  w_array <- weights.ensemble_core(
    object = object$ens_models, methods = methods, 
    X_tr = X, Y_tr = Y, X_te = Xnew, quiet = quiet)
  
  if (length(object$ens_models) > 1) {
    smoother_w <- agg_array(a = w_array, w = nnls_w)
  } else {
    smoother_w <- array(w_array[, , 1], dim = dim(w_array)[-3])
  }
  
  return(smoother_w)
}


#' Ensemble smoother weights
#'
#' Extract smoother weights for a (new) covariate matrix based on the fully trained
#' models of an ensemble created by \code{\link{ensemble_core}}.
#'
#' @param object List of trained models from \code{\link{ensemble_core}}.
#' @param methods List of methods built via \code{\link{create_method}}.
#' @param X_tr Covariate matrix of training sample.
#' @param Y_tr Vector of outcomes of training sample.
#' @param X_te Covariate matrix of test sample.
#' @param quiet If \code{FALSE}, the method currently computed is printed to the console.
#' @param ... Ignore unused arguments.
#'
#' @return A three-dimensional array of dimension
#'   \code{nrow(X_te)} × \code{nrow(X_tr)} × \code{length(methods)} containing smoother weights.
#'   Each slice corresponds to one method, where each row gives the weights
#'   assigned to training outcomes in predicting the respective test observation.
#'
#' @method weights ensemble_core
#' @keywords internal
weights.ensemble_core <- function(object, 
                                  methods,
                                  X_tr, Y_tr, X_te,
                                  quiet = TRUE,
                                  ...) {

  w_array <- array(data = 0, dim = c(nrow(X_te), nrow(X_tr), length(object)))
  
  for (i in 1:length(object)) {
    if (isFALSE(quiet)) print(paste0("Smoother weights: ", methods[[i]]$method))
    wrapper <- paste0(methods[[i]]$method, "_fit")
    
    x_select <- methods[[i]]$x_select
    X_tr_sel <- if (is.null(x_select)) X_tr else X_tr[, x_select]
    X_te_sel <- if (is.null(x_select)) X_te else X_te[, x_select]
    
    ## Smoother weights extraction
    w_array[, , i] <- do.call(paste0("weights.", wrapper), list(object[[i]], X = X_tr_sel, Y = Y_tr, Xnew = X_te_sel))
  }
  
  return(w_array)
}


#' Compute smoother weights from arithmetic mean
#'
#' Generate smoother weights for a test sample using arithmetic mean fitting of the training sample.
#'
#' @param mean_fit An object returned by \code{\link{mean_fit}}.
#' @param Xnew Covariate matrix of the test sample.
#' @param ... Ignored arguments.
#'
#' @return A matrix of smoother weights.
#'
#' @keywords internal
#' @method weights mean_fit
weights.mean_fit <- function(mean_fit, Xnew, ...) {
  w <- matrix(1 / mean_fit$N, nrow = nrow(Xnew), ncol = mean_fit$N)

  return(w)
}


#' Compute smoother weights from OLS
#'
#' Generate smoother weights for a test sample from an OLS model.
#'
#' @param ols_fit An object returned by \code{\link{ols_fit}}.
#' @param X Covariate matrix of the training sample.
#' @param Xnew Covariate matrix of the test sample.
#' @param ... Ignored arguments.
#'
#' @return A matrix of smoother weights.
#'
#' @keywords internal
#' @method weights ols_fit
weights.ols_fit <- function(ols_fit, X, Xnew, ...) {
  X <- add_intercept(X)
  Xnew <- add_intercept(Xnew)

  # Remove variables dropped due to multi-collinearity
  X <- X[, !is.na(ols_fit)]
  Xnew <- Xnew[, !is.na(ols_fit)]

  # Hat matrix
  hat_mat <- Xnew %*% solve(crossprod(X), tol = 2.225074e-308) %*% t(X)

  return(hat_mat)
}


#' Compute smoother weights from ridge regression
#'
#' Generate smoother weights for a test sample from a ridge regression model.
#'
#' @param ridge_fit An object returned by \code{\link{ridge_fit}}.
#' @param X Covariate matrix of the training sample.
#' @param Y Vector of outcomes of the training sample.
#' @param Xnew Optional covariate matrix of the test sample. If not provided, predictions are computed for the training sample.
#' @param ... Ignored arguments.
#'
#' @return A matrix of smoother weights.
#'
#' @keywords internal
#' @method weights ridge_fit
weights.ridge_fit <- function(ridge_fit, X, Y, Xnew = NULL, ...) {
  if (is.null(Xnew)) Xnew <- X
  N <- nrow(X)

  X <- scale(X, ridge_fit$x_means, ridge_fit$x_sds)
  X <- add_intercept(X)
  p <- ncol(X) - 1
  Xnew <- scale(Xnew, ridge_fit$x_means, ridge_fit$x_sds)
  Xnew <- add_intercept(Xnew)

  sd_y <- sqrt(stats::var(Y) * ((N - 1) / N))
  lambda <- (1 / sd_y) * ridge_fit$lambda.min * N

  # Hat matrix. Reference: https://stats.stackexchange.com/questions/129179/why-is-glmnet-ridge-regression-giving-me-a-different-answer-than-manual-calculat
  hat_mat <- Xnew %*% solve(crossprod(X) + lambda * diag(x = c(0, rep(1, p)))) %*% t(X)

  return(hat_mat)
}


#' Compute smoother weights from post-Lasso
#'
#' Generate smoother weights for a test sample from a post-Lasso regression model.
#'
#' @param plasso_fit An object returned by \code{\link{plasso_fit}}.
#' @param Xnew Optional covariate matrix of the test sample. 
#'  If not provided, predictions are computed for the training sample.
#' @param ... Ignored arguments.
#'
#' @return A matrix of smoother weights.
#'
#' @keywords internal
#' @method weights plasso_fit
weights.plasso_fit <- function(plasso_fit, Xnew = NULL, ...) {
  if (is.null(Xnew)) Xnew <- plasso_fit$X

  X <- add_intercept(plasso_fit$x)
  Xnew <- add_intercept(Xnew)

  colnames(X)[1] <- "(Intercept)"
  colnames(Xnew) <- colnames(X)

  Xact <- X[, plasso_fit$names_pl, drop = FALSE]
  Xactnew <- Xnew[, plasso_fit$names_pl, drop = FALSE]

  hat_mat <- Xactnew %*% solve(crossprod(Xact), tol = 2.225074e-308) %*% t(Xact)

  return(hat_mat)
}

#' Smoother weights from \code{hdm::rlasso} prediction
#'
#' Extract smoother weights for test sample from a post-lasso regression model.
#'
#' @param rlasso_fit Output of \code{\link{rlasso_fit}}.
#' @param Xnew Covariate matrix of test sample. If not provided, prediction is done for the training sample.
#' @param ... Ignore unused arguments.
#'
#' @return A matrix of smoother weights.
#'
#' @method weights rlasso_fit
#' @keywords internal
weights.rlasso_fit <- function(rlasso_fit, Xnew = NULL, ...) {
  if (!rlasso_fit$options$post) stop("Smoother matrix requires post-lasso specification.")
  if (is.null(Xnew)) Xnew <- rlasso_fit$model

  X <- add_intercept(rlasso_fit$model)
  Xnew <- add_intercept(Xnew)

  colnames(X)[1] <- "(Intercept)"
  colnames(Xnew) <- colnames(X)
  
  active_vars <- which(rlasso_fit$index)
  active_vars <- c("(Intercept)", names(rlasso_fit$beta)[active_vars])

  Xact <- X[, active_vars, drop = FALSE]
  Xactnew <- Xnew[, active_vars, drop = FALSE]

  hat_mat <- Xactnew %*% solve(crossprod(Xact), tol = 2.225074e-308) %*% t(Xact)

  return(hat_mat)
}


#' Smoother weights from random forest model
#'
#' Extract smoother weights for test sample from a random forest model.
#'
#' @param forest_grf_fit Output of \code{\link{forest_grf_fit}}.
#' @param Xnew Covariate matrix of test sample. If not provided, 
#'  prediction is done for the training sample.
#' @param ... Ignore unused arguments.
#'
#' @return A matrix of smoother weights.
#'
#' @method weights forest_grf_fit
#' @keywords internal
weights.forest_grf_fit <- function(forest_grf_fit, Xnew = NULL, ...) {
  if (is.null(Xnew)) Xnew <- forest_grf_fit$X.orig

  w <- grf::get_forest_weights(forest_grf_fit, newdata = Xnew)
  w <- as.matrix(w)

  return(w)
}


#' Smoother weights from XGBoost model
#'
#' Extract smoother weights for test sample from an XGBoost model.
#'
#' @param xgboost_fit Output of \code{\link{xgboost_fit}}.
#' @param X Covariate matrix of training sample.
#' @param Xnew Covariate matrix of test sample. 
#'  If not provided, prediction is done for the training sample.
#' @param ... Ignore unused arguments.
#'
#' @return A matrix of smoother weights.
#'
#' @method weights xgboost_fit
#' @keywords internal
weights.xgboost_fit <- function(xgboost_fit, X, Xnew = NULL, ...) {
  dtrain <- xgboost::xgb.DMatrix(data = as.matrix(X))
  dtest <- xgboost::xgb.DMatrix(data = as.matrix(Xnew))

  w <- get_xgboost_weights(xgboost_fit, dtrain = dtrain, dtest = dtest)
  w <- as.matrix(w$S_test)

  return(w)
}


#' Smoother weights from k-nearest neighbor algorithm
#'
#' Extract smoother weights from a k-nearest neighbor algorithm.
#' Each weight equals \code{1 / k} for neighbors and 0 for non-neighbors.
#'
#' @param arguments Output of \code{\link{knn_fit}}.
#' @param X Covariate matrix of training sample.
#' @param Xnew Covariate matrix of test sample.
#' @param ... Ignore unused arguments.
#'
#' @return A matrix of smoother weights.
#'
#' @method weights knn_fit
#' @keywords internal
weights.knn_fit <- function(arguments, X, Xnew = NULL, ...) {
  if (is.null(Xnew)) Xnew <- X
  k <- if (is.null(arguments$k)) 10 else arguments$k

  # distance = as.matrix(FastKNN::Distance_for_KNN_test(Xnew, X))
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

  return(w)
}


#' Smoother weights from distributional random forest model
#'
#' Extract smoother (adaptive nearest neighbor) weights for test sample
#' from a distributional random forest model.
#'
#' @param forest_drf_fit Output of \code{\link{forest_drf_fit}}.
#' @param Xnew Covariate matrix of test sample. 
#'  If not provided, prediction is done for the training sample.
#' @param ... Ignore unused arguments.
#'
#' @return A matrix of smoother weights.
#'
#' @method weights forest_drf_fit
#' @keywords internal
weights.forest_drf_fit <- function(forest_drf_fit, Xnew = NULL, ...) {
  if (is.null(Xnew)) Xnew <- forest_drf_fit$X.orig

  w <- as.matrix(predict(forest_drf_fit, newdata = Xnew)$weights)

  return(w)
}


#' Determine the Corresponding Data Subset for a Model Element
#'
#' Helper function to identify which column in \code{d_mat} or \code{z_mat} corresponds to
#' a given model element from \code{np_models$Y.hat.d_m} or \code{np_models$Y.hat.z_m}.
#'
#' @param sub_element A model element (e.g., from \code{np_models$Y.hat.d_m} or \code{np_models$Y.hat.z_m}).
#' @param np_models A list containing ensemble models with components \code{Y.hat.d_m} and \code{Y.hat.z_m}.
#' @param d_mat Logical matrix of treatment indicators.
#' @param z_mat Logical matrix of instrument indicators.
#'
#' @return The corresponding column from \code{d_mat} or \code{z_mat} if a match is found;
#'   otherwise, \code{NULL}.
#'
#' @keywords internal
get_subset <- function(sub_element, np_models, d_mat, z_mat) {
  # Check for matches in Y.hat.d_m
  for (i in seq_along(np_models$Y.hat.d_m)) {
    if (identical(sub_element, np_models$Y.hat.d_m[[i]])) {
      return(d_mat[, i])
    }
  }
  # Check for matches in Y.hat.z_m
  for (i in seq_along(np_models$Y.hat.z_m)) {
    if (identical(sub_element, np_models$Y.hat.z_m[[i]])) {
      return(z_mat[, i])
    }
  }
  NULL
}
