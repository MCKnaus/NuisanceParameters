#' Extraction of smoother weights for ensemble learner
#'
#' @description
#' Extract smoother weights from short- or standard-stacked ensemble model 
#' created by \code{\link{ensemble_short}} or \code{\link{ensemble}}.
#'
#' @param NuPa String vector specifying which nuisance parameters to extract smoother matrices for.
#' @param np_object Stacked ensemble learner from \code{\link{ensemble_short}} or \code{\link{ensemble}}.
#' @param method List of models by \code{\link{create_method}} that was used as
#' iobjectut for \code{\link{ensemble_short}} or \code{\link{ensemble}}.
#' @param X Covariate matrix of training sample.
#' @param Y Vector of outcomes of training sample.
#' @param d_mat Logical matrix of treatment indicators (N X T+1). For example created by \code{\link{prep_w_mat}}.
#' @param z_mat Logical matrix of instrument indicators (N X T+1). For example created by \code{\link{prep_w_mat}}.
#' @param subset Logical vector indicating which observations to use for determining
#' ensemble weights. If not provided, all observations are used.
#' @param cv Number of cross-validation folds when estimating ensemble model.
#' @param cf_mat Logical matrix with k columns of indicators representing the different folds
#' (for example created by \code{\link{prep_cf_mat}}).
#' @param quiet If FALSE, method that is currently computed is printed into the
#' console.
#' @param ... Ignore unused arguments
#'
#' @return Matrix of dimension \code{nrow(X)} X \code{nrow(X)} containing
#' ensemble smoother weights, or list of such matrices if multiple models are requested.
#'
#' @export
get_outcome_weights = function(np_object, 
                               NuPa, 
                               quiet = TRUE,
                               subset = NULL,
                               ...) {
  
  ## Input checks 
  
  if (is.null(np_object)) {
    stop("Ensemble models were not saved after training.")
    
  } else if (is.character(np_object) && length(np_object) == 1) {
    loaded_object <- readRDS(np_object)
    np_object <- list(models = loaded_object$models, 
                      numbers = loaded_object$numbers)
    
  } else if (!is.list(np_object)) {
    stop("np_object must be either a list or a valid path to an RDS file.")
  }
  
  model_names <- paste0(NuPa, "_m")
  existing_models <- model_names[model_names %in% names(np_object[["models"]])]
  missing_models <- setdiff(model_names, existing_models)
  
  if (length(missing_models) > 0) {
    message("The following ensemble models were not found: ", 
            paste(missing_models, collapse = ", "), "\n",
            "Continuing with available models: ", 
            paste(existing_models, collapse = ", "))
  }
  
  if (length(existing_models) == 0) {
    stop("No valid (ensemble) models found in np_object[['models']].")
  }
  
  if (is.null(subset)) subset = rep(TRUE, length(Y))
  
  # Extract numbers 
  
  cv = np_object$numbers$cv
  d_mat = np_object$numbers$d_mat
  z_mat = np_object$numbers$z_mat
  cf_mat = np_object$numbers$cf_mat
  method = np_object$numbers$method
  X = np_object$numbers$X
  Y = np_object$numbers$Y
  
  # Helper function to determine subset
  
  get_subset <- function(sub_element, np_models) {
    if (identical(sub_element, np_models[["Y.hat.d_m"]][[1]])) return(d_mat[, 1])
    if (identical(sub_element, np_models[["Y.hat.d_m"]][[2]])) return(d_mat[, 2])
    if (identical(sub_element, np_models[["Y.hat.z_m"]][[1]])) return(z_mat[, 1])
    if (identical(sub_element, np_models[["Y.hat.z_m"]][[2]])) return(z_mat[, 2])
    NULL
  }
  
  # Process each model
  
  results <- list()
  results <- lapply(existing_models, function(model) {
    current_model <- np_object[["models"]][[model]]
    
    if (cv == 1) {
      if ("ens_object" %in% names(current_model)) {
        return(get_smoother(current_model, method = method, X = X, Y = Y, 
                            nnls_w = NULL, cv = cv, cf_mat = cf_mat, quiet = quiet))
      } else {
        return(lapply(seq_along(current_model), function(i) {
          sub_element <- current_model[[i]]
          if (is.list(sub_element) && "ens_object" %in% names(sub_element)) {
            subset <- get_subset(sub_element, np_object[["models"]])
            get_smoother(sub_element, method = method, X = X, Y = Y, 
                         subset = subset, nnls_w = NULL, cv = cv, 
                         cf_mat = cf_mat, quiet = quiet)
          }
        }))
      }
    } else { # cv > 1
      if (!inherits(current_model, "ens.learner")) {
        return(lapply(seq_along(current_model), function(i) {
          sub_element <- current_model[[i]]
          if (is.list(sub_element)) {
            subset <- get_subset(sub_element, np_object[["models"]])
            get_smoother(sub_element, method = method, X = X, Y = Y,
                         subset = subset, nnls_w = NULL, cv = cv,
                         cf_mat = cf_mat, quiet = quiet)
          }
        }))
      } else {
        return(get_smoother(current_model, method = method, X = X, Y = Y,
                            subset = NULL, nnls_w = NULL, cv = cv,
                            cf_mat = cf_mat, quiet = quiet))
      }
    }
  })
  
  names(results) <- existing_models
  return(results)
}


#' Extraction of smoother weights for ensemble learner
#'
#' @description
#' Extract smoother weights from short- or standard-stacked ensemble model 
#' created by \code{\link{ensemble_short}} or \code{\link{ensemble}}.
#'
#' @param object Stacked ensemble learner from \code{\link{ensemble_short}} or \code{\link{ensemble}}.
#' @param method List of models by \code{\link{create_method}} that was used as
#' input for \code{\link{ensemble_short}} or \code{\link{ensemble}}.
#' @param X Covariate matrix of training sample.
#' @param Y Vector of outcomes of training sample.
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
#' @return Matrix of dimension \code{nrow(X)} X \code{nrow(X)} containing
#' ensemble smoother weights.
#'
#' @export
get_smoother = function(object,
                        method,
                        X, Y,
                        subset = NULL,
                        nnls_w = NULL,
                        cv = cv,
                        cf_mat,
                        quiet = TRUE,
                        ...) {
  
  if (is.null(object)) stop("Ensemble models were not saved after training.")
  if (is.null(subset)) subset = rep(TRUE, length(Y))
  
  
  ### Short-Stacking ###
  if (cv == 1) {
    w = weights(object[["ens_object"]], method = method, X = X, Y = Y, subset = subset, w = object[["nnls_w"]], cf_mat = cf_mat, quiet = quiet)
    
    
    ### Standard-Stacking ###
    
  } else if (cv > 1) {
    w = matrix(0, nrow = nrow(X), ncol = nrow(X))
    
    for (i in 1:ncol(cf_mat)) {
      fold = cf_mat[, i]
      X_tr = X[!fold & subset, ]
      Y_tr = Y[!fold & subset]
      X_te = X[fold, ]
      
      w_array = weights(object[[i]][["ens_object"]], method, X = X_tr, Y = Y_tr, Xnew = X_te, quiet = quiet)
      w[fold, !fold & subset] = apply(w_array, c(1, 2), function(x) sum(x * object[[i]][["ens_object"]]$nnls_weights))
    }
    
  }
  
  return(w)
  
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
  
  w = matrix(1 / mean_fit$N, nrow = nrow(Xnew), ncol = mean_fit$N)
  
  return(w)
  
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
  N = nrow(X)
  
  X = scale(X, ridge_fit$x_means, ridge_fit$x_sds)
  X = add_intercept(X)
  Xnew = scale(Xnew, ridge_fit$x_means, ridge_fit$x_sds)
  Xnew = add_intercept(Xnew)
  
  p = ncol(X) - 1
  
  sd_y = sqrt(stats::var(Y) * ((N - 1) / N))
  lambda = (1 / sd_y) * ridge_fit$lambda.min * N
  
  # calculate hat matrix
  # reference: https://stats.stackexchange.com/questions/129179/why-is-glmnet-ridge-regression-giving-me-a-different-answer-than-manual-calculat
  hat_mat = Xnew %*% solve(crossprod(X) + lambda * diag(X = c(0, rep(1, p)))) %*% t(X)
  
  return(hat_mat)
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
  
  Xact = X[, plasso_fit$names_pl, drop = FALSE]
  Xactnew = Xnew[, plasso_fit$names_pl, drop = FALSE]
  
  hat_mat = Xactnew %*% solve(crossprod(Xact), tol = 2.225074e-308) %*% t(Xact)
  
  return(hat_mat)
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


#' Smoother weights from XGBoost model
#'
#' @description
#' Extract smoother weights for test sample from a XGBoost model.
#'
#' @param xgboost_fit Output of \code{\link{xgboost_fit}}
#' @param Xnew Covariate matrix of test sample.
#' If not provided, prediction is done for the training sample.
#' @param ... Ignore unused arguments
#'
#' @return Matrix of smoother weights.
#'
#' @method weights xgboost_fit
#'
#' @keywords internal
#'
weights.xgboost_fit = function(xgboost_fit, X, Xnew = NULL, ...) {
  
  dtrain = xgboost::xgb.DMatrix(data = as.matrix(X))
  dtest = xgboost::xgb.DMatrix(data = as.matrix(Xnew))
  
  w = get_xgboost_weights(xgboost_fit, dtrain=dtrain, dtest=dtest)
  w = as.matrix(w$S_test)
  
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
