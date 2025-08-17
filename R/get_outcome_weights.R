#' Extraction of smoother weights for ensemble learner
#'
#' @description
#' Extract smoother weights from short- or standard-stacked ensemble model 
#' created by \code{\link{ensemble_short}} or \code{\link{ensemble}}.
#'
#' @param np_object A \code{NuisanceParameters} object created by 
#'   \code{\link{nuisance_parameters}} function
#' @param NuPa String vector specifying which nuisance parameters to extract smoother matrices for.
#' @param subset Logical vector indicating which observations to use for determining
#' ensemble weights. If not provided, all observations are used.
#' @param quiet If FALSE, method that is currently computed is printed into the
#' console.
#' @param ... Additional arguments
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
  
  ## Input check
  load_np_object <- function(obj) {
    if (inherits(obj, "NuisanceParameters")) return(obj)
    if (is.character(obj) && length(obj) == 1) {
      loaded <- readRDS(obj)
      if (!is.list(loaded) || is.null(loaded$models)) stop("Loaded object is not a valid NuisanceParameters list.")
      return(list(models = loaded$models, numbers = loaded$numbers))
    }
    stop("np_object must be a NuisanceParameters object or a valid path to an RDS file.")
  }
  
  np_object <- load_np_object(np_object)
  
  
  ## Lasso check
  if (any(sapply(np_object[["numbers"]][["method"]], function(method_group) {
    any(sapply(method_group, function(x) identical(x[["method"]], "lasso")))
  }))) { stop("Smoother matrices cannot be extracted if method = 'lasso' was used.") }
  
  
  ## Models check
  model_names <- paste0(NuPa, "_m")
  existing_models <- model_names[model_names %in% names(np_object[["models"]])]
  missing_models  <- setdiff(model_names, existing_models)
  
  if (length(missing_models) > 0) {
    message("Missing models: ", paste(missing_models, collapse = ", "), 
            "\nContinuing with available models: ", paste(existing_models, collapse = ", "))}
  if (length(existing_models) == 0) stop("No valid (ensemble) models found in np_object[['models']].")
  
  
  ## Extract numbers 
  cv = np_object$numbers$cv
  d_mat = np_object$numbers$d_mat
  z_mat = np_object$numbers$z_mat
  cf_mat = np_object$numbers$cf_mat
  method = np_object$numbers$method
  X = np_object$numbers$X
  Y = np_object$numbers$Y
  
  if (is.null(subset)) subset = rep(TRUE, length(Y))
  results = list()
  
  ## Process each model
  for (model_name in existing_models) {
    current_model <- np_object[["models"]][[model_name]]
    method_name <- sub("_m$", "", model_name)
    model_results <- list()
    
    ## Short stacking
    if (cv == 1) {
      # Case 1: Single ensemble model (i.e. Y.hat)
      if ("ens_object" %in% names(current_model)) {
        model_results <- get_smoother(current_model, 
                                      method = method[[method_name]], 
                                      X = X, 
                                      Y = Y,
                                      cv = cv, 
                                      cf_mat = cf_mat, 
                                      quiet = quiet)
      } 
      # Case 2: List of models (i.e. Y.hat.d)
      else {
        for (i in seq_along(current_model)) {
          sub_element <- current_model[[i]]
          if (is.list(sub_element) && "ens_object" %in% names(sub_element)) {
            subset <- get_subset(sub_element, np_object[["models"]], d_mat, z_mat)
            model_results[[i]] <- get_smoother(sub_element,
                                               method = method[[method_name]],
                                               X = X,
                                               Y = Y,
                                               subset = subset,
                                               cv = cv,
                                               cf_mat = cf_mat,
                                               quiet = quiet)
          }
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
            subset <- get_subset(sub_element, np_object[["models"]], d_mat, z_mat)
            model_results[[i]] <- get_smoother(sub_element,
                                               method = method[[method_name]],
                                               X = X,
                                               Y = Y,
                                               subset = subset,
                                               cv = cv,
                                               cf_mat = cf_mat,
                                               quiet = quiet)
          }
        }
      }
      # Case 2: Single ensemble model
      else {
        model_results <- get_smoother(current_model,
                                      method = method[[method_name]],
                                      X = X,
                                      Y = Y,
                                      subset = NULL,
                                      cv = cv,
                                      cf_mat = cf_mat,
                                      quiet = quiet)
      }
    }
    
    results[[NuPa]] <- model_results
  }
  
  names(results) <- NuPa
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
#' @param cv Number of cross-validation folds when estimating ensemble model.
#' @param cf_mat Logical matrix with k columns of indicators representing the different folds
#' (for example created by \code{\link{prep_cf_mat}}).
#' @param quiet If FALSE, method that is currently computed is printed into the
#' console.
#'
#' @return Matrix of dimension \code{nrow(X)} X \code{nrow(X)} containing 
#' ensemble smoother weights.
#'
#' @export
get_smoother = function(object,
                        method,
                        X, Y,
                        subset = NULL,
                        cf_mat, cv,
                        quiet = TRUE
                        ) {
  
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
  hat_mat = Xnew %*% solve(crossprod(X) + lambda * diag(x = c(0, rep(1, p)))) %*% t(X)
  
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
  
  X = add_intercept(plasso_fit$x)
  Xnew = add_intercept(Xnew)
  
  colnames(X)[1] = "(Intercept)"
  colnames(Xnew) = colnames(X)
  
  Xact = X[, plasso_fit$names_pl, drop = FALSE]
  Xactnew = Xnew[, plasso_fit$names_pl, drop = FALSE]
  
  hat_mat = Xactnew %*% solve(crossprod(Xact), tol = 2.225074e-308) %*% t(Xact)
  
  return(hat_mat)
}


#' Smoother weights from hdm::rlasso prediction
#'
#' @description
#' Extract smoother weights for test sample from a Lasso regression model.
#' Only works if post-lasso estimation was performed.
#'
#' @param rlasso_fit Output of \code{\link{rlasso_fit}}
#' @param Xnew Covariate matrix of test sample.
#' If not provided, prediction is done for the training sample.
#' @param ... Ignore unused arguments
#'
#' @return Matrix of smoother weights.
#'
#' @method weights rlasso_fit
#'
#' @keywords internal
#'
weights.rlasso_fit = function(rlasso_fit, Xnew = NULL, ...) {
  
  # Check if post-lasso was run
  if (!isTRUE(rlasso_fit[["options"]][["post"]])) {
    stop("Smoother matrix can be extracted only if post-lasso specification is run")
  }
  
  if (is.null(Xnew)) Xnew <- rlasso_fit$model
  
  X = add_intercept(rlasso_fit$model)
  Xnew = add_intercept(Xnew)
  
  colnames(X)[1] = "(Intercept)"
  colnames(Xnew) = colnames(X)
  
  active_vars <- which(rlasso_fit$index)
  active_vars <- c("(Intercept)", names(rlasso_fit$beta)[active_vars])
  
  Xact = X[, active_vars, drop = FALSE]
  Xactnew = Xnew[, active_vars, drop = FALSE]
  
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
  
  w = grf::get_forest_weights(forest_grf_fit, newdata = Xnew)
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
#' @param forest_drf_fit Output of \code{\link{forest_drf_fit}}
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


#' Extraction of smoother weights for short-stacked ensemble learner
#'
#' @description
#' Extract smoother weights from short-stacked ensemble model created by
#' \code{\link{ensemble_short}}.
#'
#' @param object Short-stacked ensemble learner from \code{\link{ensemble_short}}.
#' @param method List of method models by \code{\link{create_method}} that was used as
#' input for \code{\link{ensemble_short}}.
#' @param X Covariate matrix of training sample.
#' @param Y Vector of outcomes of training sample.
#' @param subset Logical vector indicating which observations to use for determining
#' ensemble weights. If not provided, all observations are used.
#' @param w Ensemble weights to aggregate predictions from different learners (optional).
#' Needs to be a vector of length \code{ncol(object$fit_cv)}.
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
#'
#' @method weights ensemble_short
#'
weights.ensemble_short = function(object,
                                  method,
                                  X, Y,
                                  subset = NULL,
                                  w = NULL,
                                  cf_mat,
                                  quiet = TRUE,
                                  ...) {
  
  # load fitted method object
  if(is.null(object$method)) {
    stop("Ensemble models were not saved after training.")
  } else if (is.character(object$method) & length(object$method) == 1) {
    method_fit = readRDS(object$method)
  } else {
    method_fit = object$method
  }
  
  if (is.null(w)) w = rep(1 / ncol(object$fit_cv), ncol(object$fit_cv))
  if (is.null(subset)) subset = rep(TRUE, length(Y))
  
  # create matrix to store smoother weights
  smoother_weights = matrix(0, nrow = nrow(X), ncol = nrow(X))
  
  
  for (i in 1:ncol(cf_mat)) {
    
    fold = cf_mat[, i]
    
    test_idx = which(fold)
    train_idx = which(!fold)
    train_idx_t = which(!fold & subset)
    
    # Create mapping from train_idx_t to positions in train_idx
    train_pos_in_full = match(train_idx_t, train_idx)
    
    X_tr = X[train_idx_t, ]
    Y_tr = Y[train_idx_t]
    X_te = X[test_idx, ]
    
    # Get T-learner smoother weights (test × train&subset)
    w_array_raw = weights.ensemble_core(object = method_fit[[i]], method = method, X_tr = X_tr, Y_tr = Y_tr, X_te = X_te, quiet)
    
    # Pad T-learner weights (to test × full train) and insert raw weights
    padded_w = array(0, dim = c(length(test_idx), length(train_idx), length(method)))
    padded_w[, train_pos_in_full, ] = w_array_raw      
    
    smoother_weights[fold, !fold] = agg_array(padded_w, w)
    
  }
  
  return(smoother_weights)
  
}


#' Smoother weights from cross-validated ensemble learner
#'
#' @description
#' \code{\link{weights.ensemble}} extract ensemble smoother weights for some
#' (new) covariate matrix based on a fitted and cross-validated ensemble model from
#' \code{\link{ensemble}}.
#'
#' @param object Trained ensemble object from \code{\link{ensemble}}
#' @param method List of raw methods built via \code{\link{create_method}}
#' @param X Covariate matrix of training sample
#' @param Y Vector of outcomes of training sample
#' @param Xnew Covariate matrix of test sample
#' @param quiet If FALSE, method that is currently computed is printed into the
#' console
#' @param ... Ignore unused arguments
#'
#' @return Matrix of dimension \code{nrow(Xnew)} X \code{nrow(X)} containing
#' ensemble smoother weights.
#'
#' @export
#'
#' @method weights ensemble
#'
weights.ensemble = function(object,
                            method,
                            X, Y, Xnew,
                            quiet = TRUE,
                            ...) {
  
  if (length(object$method_fit) > 1) {
    
    w_array = weights.ensemble_core(object = object$method_fit, method = method, X_tr = X, Y_tr = Y, X_te = Xnew, quiet)
    smoother_weights = agg_array(w_array, object$nnls_weights)
    
  } else if (length(object$method_fit) == 1) {
    
    w_array = weights.ensemble_core(object = object$method_fit, method = method, X_tr = X, Y_tr = Y, X_te = Xnew, quiet)
    smoother_weights = array(w_array[, , 1], dim = dim(w_array)[-3])
    
  }
  
  return(smoother_weights)
}


#' Ensemble smoother weights
#'
#' @description
#' Extract smoother weights for some (new) covariate matrix based on the fully
#' trained models of the ensemble model from \code{\link{ensemble_core}}.
#'
#' @param object List of fitted methods (from \code{\link{ensemble_core}}
#' @param method List of raw methods built via \code{\link{create_method}}
#' @param X_tr Covariate matrix of training sample
#' @param Y_tr Vector of outcomes of training sample
#' @param X_te Covariate matrix of test sample
#' @param quiet If FALSE, method that is currently computed is printed into the
#' console
#' @param ... Ignore unused arguments
#'
#' @return Returns a three-dimensional array of dimension
#' \code{nrow(Xnew)} X \code{nrow(X)} X \code{length(method)}
#' containing the smoother weights that deliver predictions of each method
#' where each row gives the weight that each training outcome received in the
#' prediction of the respective test set observation.
#'
#' @export
#'
#' @method weights ensemble_core
#'
weights.ensemble_core = function(object, method,
                                 X_tr, Y_tr, X_te,
                                 quiet = TRUE,
                                 ...) {
  
  # initialize 3D array to be filled
  w_array = array(data = 0, dim = c(nrow(X_te), nrow(X_tr), length(object)))
  
  # loop over all trained models
  for (i in 1:length(object)) {
    
    if (isFALSE(quiet)) print(paste0("Smoother weights: ", method[[i]]$method))
    
    wrapper  <- paste0(method[[i]]$method, "_fit")
    
    # browser()
    
    # extract smoother weights
    if (is.null(method[[i]]$x_select)) {
      temp = do.call(paste0("weights.", wrapper), list(object[[i]], X = X_tr, Y = Y_tr, Xnew = X_te))
    } else {
      temp = do.call(paste0("weights.", wrapper), list(object[[i]], X = X_tr[, method[[i]]$x_select], Y = Y_tr, Xnew = X_te[, method[[i]]$x_select]))
    }
    
    w_array[, , i] = temp
    
  }
  
  return(w_array)
  
}


#' Determine the Corresponding Data Subset for a Model Element
#'
#' @description
#' Helper function to identify which column in `d_mat` or `z_mat` corresponds to 
#' a given model element from `np_models$Y.hat.d_m` or `np_models$Y.hat.z_m`.
#'
#' @param sub_element A model element (e.g., from `np_models$Y.hat.d_m` or `np_models$Y.hat.z_m`).
#' @param np_models A list containing ensemble models with components `Y.hat.d_m` and `Y.hat.z_m`.
#' @param d_mat Logical matrix of treatment indicators.
#' @param z_mat Logical matrix of instrument indicators.
#' @param x Covariate matrix.
#' @param z_mat Matrix of instrumental variables (columns correspond to `Y.hat.z_m` models).
#'
#' @return The corresponding column from `d_mat` or `z_mat` if a match is found; 
#'   otherwise, `NULL`.
#'
#' @noRd
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