#' Tune learner hyperparameters
#'
#' Tunes hyperparameters for specified learners in a \code{methods} list.
#'
#' @param type Tuning type, either "\code{full_sample}" or "\code{fold}".
#' @param methods List of methods for ensemble estimation.
#' @param X Covariate matrix.
#' @param Y Numeric vector containing the outcome variable.
#'
#' @return Modified methods list with tuned parameters.
#'
#' @keywords internal
#'
tune_learners <- function(type = c("full_sample", "fold"),
                          methods,
                          X, Y) {
  type <- match.arg(type)
  
  to_tune <- which(vapply(methods, function(m)
    isTRUE(m$tuning == type), logical(1L)))
  
  if (!length(to_tune)) return(methods)
  
  for (i in to_tune) {
    mtd <- methods[[i]]
    X_sub <- if (is.null(mtd$x_select)) X else X[, mtd$x_select, drop = FALSE]
    
    fun_name <- paste0(mtd$method, "_tune")
    wrapper  <- get0(fun_name, mode = "function", inherits = TRUE)
    
    if (is.null(wrapper)) {
      warning("No tuner found for method '", mtd$method, "'. Skipping.")
      next
    }
    
    tuned_params <- wrapper(X = X_sub, Y = Y)
    
    args <- methods[[i]][["arguments"]]
    if (is.null(args)) args <- list()
    
    methods[[i]][["arguments"]] <- utils::modifyList(args, tuned_params)
  }
  
  return(methods)
}


#' Tune \code{ranger} Random Forest Hyperparameters for Regression
#'
#' Tunes \code{mtry}, \code{min.node.size}, and \code{sample.fraction} 
#' parameters for ranger random forests for regression tasks. 
#' Requires \pkg{tuneRanger} and \pkg{mlr} packages.
#'
#' @param X Covariate matrix or data frame.
#' @param Y Numeric vector containing the outcome variable.
#'
#' @return List containing tuned hyperparameters: mtry, min.node.size, and sample.fraction.
#' 
#' @details 
#' Tuning uses 1000 trees, 2 threads, and 70 iterations for tuning.
#'
#' @keywords internal
#' 
ranger_tune <- function(X, Y) {
  if (!requireNamespace("tuneRanger", quietly = TRUE) || 
      !requireNamespace("mlr", quietly = TRUE)) {
    stop("To enable ranger tuning, please install both tuneRanger and mlr packages.\n",
         "Install them with: install.packages(c('tuneRanger', 'mlr'))",
         call. = FALSE)
  }
  
  data <- data.frame(Y = Y, X)
  
  # Create an mlr regression task
  task = mlr::makeRegrTask(data = data, target = "Y")
  
  tuned_object = tuneRanger::tuneRanger(
    task, 
    measure = list(mlr::rmse), 
    num.trees = 1000, 
    num.threads = 2, 
    iters = 70, 
    show.info = FALSE
  )
  
  # Keep the tuned hyperparameters: mtry, min.node.size & sample.fraction
  tuned_params <- as.list(tuned_object[["recommended.pars"]][1, 1:3])
  return(tuned_params)
}


#' Tune \code{ranger} Random Forest Hyperparameters for Classification
#'
#' Tunes \code{mtry}, \code{min.node.size}, and \code{sample.fraction} 
#' parameters for ranger random forests for classification tasks. Uses 
#' \code{loglos} for binary classification and \code{multiclass.brier} 
#' for multi-class classification. Requires \pkg{tuneRanger} and \pkg{mlr} packages.
#'
#' @param X Covariate matrix or data frame
#' @param Y Factor vector containing the outcome variable
#'
#' @return List containing tuned hyperparameters: mtry, min.node.size, and sample.fraction
#' 
#' @details 
#' Tuning uses 1000 trees, 2 threads, and 70 iterations for tuning.
#'
#' @keywords internal
#' 
ranger_prop_tune <- function(X, Y) {
  if (!requireNamespace("tuneRanger", quietly = TRUE) || 
      !requireNamespace("mlr", quietly = TRUE)) {
    stop("To enable ranger tuning, please install both tuneRanger and mlr packages.\n",
         "Install them with: install.packages(c('tuneRanger', 'mlr'))",
         call. = FALSE)
  }
  
  data <- data.frame(Y = Y, X)
  
  # Binary or multiclass classification?
  if (length(unique(Y)) == 2) {
    measure <- list(mlr::logloss)
  } else {
    measure <- list(mlr::multiclass.brier)
  }
  
  # Create an mlr classification task
  task = mlr::makeClassifTask(data = data, target = "Y")

  tuned_object = tuneRanger::tuneRanger(
    task, 
    measure = measure, 
    num.trees = 1000, 
    num.threads = 2, 
    iters = 70,
    show.info = FALSE
  )
  
  # Keep the tuned hyperparameters: mtry, min.node.size & sample.fraction
  tuned_params <- as.list(tuned_object[["recommended.pars"]][1, 1:3])
  return(tuned_params)
}


#' Tune \code{grf} Random Forest Hyperparameters
#'
#' Performs default \code{grf} hyperparameter tuning.
#'
#' @param X Covariate matrix or data frame.
#' @param Y Numeric vector containing the outcome variable.
#'
#' @return List containing tuned hyperparameters
#'
#' @keywords internal
#' 
forest_grf_tune <- function(X, Y) {
  if (!requireNamespace("grf", quietly = TRUE)) {
    stop("To enable grf tuning, please install the grf package.\n",
         "Install it with: install.packages('grf')",
         call. = FALSE)
  }
  
  tuned_forest <- grf::regression_forest(
    X = X,
    Y = Y,
    tune.parameters = "all",
    num.trees = 50
  )
  
  tuned_params <- tuned_forest[["tunable.params"]]
  return(tuned_params)
}


#' Hyperparameter Tuning for XGBoost via Hyperband
#'
#' Implements the Hyperband multi-armed bandit for hyperparameter optimization.
#' Progressively allocates resources (boosting rounds)
#' to the most promising configurations across successive rungs of evaluation.
#'
#' @param X Covariate matrix.
#' @param Y Numeric vector containing the outcome variable.
#' @param max_rounds Maximum number of boosting rounds for final evaluation. Default is 100.
#' @param n_configs Initial number of random hyperparameter configurations to sample. Default is 50.
#' @param eta_downfactor Downsampling factor between rungs (η in Hyperband terminology). Default is 3.
#' @param nfold Number of cross-validation folds. Default is 5.
#' @param metrics Evaluation metric for early stopping. Options include "rmse", "mae",
#'                "logloss", "error", etc. Default is "rmse".
#' @param seed Random seed. Default is 123.
#'
#' @return A list with two components:
#' \itemize{
#'   \item{\code{params}}: Named list of optimal hyperparameter values
#'   \item{\code{best_score}}: The best cross-validation score achieved (minimum test error)
#' }
#'
#' @details
#' The function implements the Hyperband algorithm with three resource levels (rungs):
#' \enumerate{
#'   \item \strong{First rung:} 10 boosting rounds - evaluate all `n_configs` configurations
#'   \item \strong{Second rung:} 30 boosting rounds - keep top 1/η configurations
#'   \item \strong{Final rung:} `max_rounds` boosting rounds - keep top configuration
#' }
#'
#' Tunable hyperparameters include the learning rate (\code{eta}), tree depth (\code{max_depth}),
#' minimum child weight (\code{min_child_weight}), \code{gamma}, column subsampling fractions
#' (\code{colsample_bytree}, \code{colsample_bylevel}, \code{colsample_bynode}),
#' L2 regularization (\code{lambda}), \code{grow_policy}, \code{max_leaves}, and \code{max_bin}.
#' The search space uses log-uniform sampling for \code{eta} and \code{lambda},
#' and uniform or discrete sampling for the others. Parameters kept fixed are:
#' \code{subsample = 1}, \code{alpha = 0}, \code{max_delta_step = 0},
#' \code{tree_method = "hist"}, \code{objective = "reg:squarederror"},
#' and \code{base_score = mean(Y)}.
#'
#' @keywords internal
#' 
xgboost_tune <- function(X, Y,
                         max_rounds = 100, 
                         n_configs = 50,
                         eta_downfactor = 3, 
                         nfold = 5,
                         metrics = "rmse", 
                         seed = 123) {
  set.seed(seed)
  
  sample_params <- function() {
    list(
      eta = 10^stats::runif(1, log10(0.01), log10(0.3)),
      max_depth = sample(2:10, 1),
      min_child_weight = stats::runif(1, 0, 10),
      gamma = stats::runif(1, 0, 10),
      colsample_bytree = stats::runif(1, 0.5, 1),
      colsample_bylevel = stats::runif(1, 0.5, 1),
      colsample_bynode = stats::runif(1, 0.5, 1),
      lambda = 10^stats::runif(1, log10(1e-3), log10(10)),
      max_leaves = sample(0:10, 1),
      max_bin = sample(10:256, 1),
      tree_method = "hist",
      ## Must be fixed:
      objective = "reg:squarederror",
      subsample = 1,
      alpha = 0,
      max_delta_step = 0,
      base_score = mean(Y)
    )
  }
  
  configs <- replicate(n_configs, sample_params(), simplify = FALSE)
  budgets <- c(15, 35, max_rounds)
  
  best_score <- Inf
  best_params <- NULL
  
  dtrain <- xgboost::xgb.DMatrix(
    data = as.matrix(X),
    label = Y
  )
  
  for (budget in budgets) {
    scores <- numeric(length(configs))
    
    for (i in seq_along(configs)) {
      params <- configs[[i]]
      cv <- tryCatch(
        xgboost::xgb.cv(
          data = dtrain,
          params = params,
          nfold = nfold,
          nrounds = budget,
          metrics = metrics,
          verbose = FALSE,
          early_stopping_rounds = 10
        ),
        error = function(e) NULL
      )
      scores[i] <- if (!is.null(cv)) 
        min(cv$evaluation_log[[paste0("test_", metrics, "_mean")]]) 
      else NA
    }
    
    # Remove configs with NA scores
    valid_idx <- which(!is.na(scores))
    if (length(valid_idx) == 0) break
    scores <- scores[valid_idx]
    configs <- configs[valid_idx]
    
    # Update best params from this rung
    rung_best_idx <- which.min(scores)
    if (scores[rung_best_idx] < best_score) {
      best_score <- scores[rung_best_idx]
      best_params <- configs[[rung_best_idx]]
    }
    
    # Keep top configs for next rung
    keep_n <- max(1, floor(length(configs) / eta_downfactor))
    top_idx <- order(scores)[1:keep_n]
    configs <- configs[top_idx]
  }
  
  best_params$base_score <- NULL
  
  return(best_params)
}