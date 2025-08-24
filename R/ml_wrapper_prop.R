#' This function calculates the (multinomial) logistic regression based on the \code{\link{glmnet}} package 
#'
#' @param X Covariate matrix of training sample
#' @param Y Vector of outcomes of training sample
#' @param arguments List of arguments passed to  \code{\link{glmnet}}
#' @import glmnet
#'
#' @return An object with S3 class "glmnet"
#'
logit_fit = function(X,Y,arguments=list()){
  if (length(unique(Y))==2){
    logit = do.call(glmnet::glmnet, c(list(x=X,y=Y, family = "binomial", type.measure = "class"),arguments))
  }else{
    logit = do.call(glmnet::glmnet, c(list(x=X,y=Y, family = "multinomial", type.measure = "class"),arguments))
  }
  logit
}


#' Prediction based on (multinomial) logistic regression.
#' @param logit_fit Output of \code{\link{glmnet}} or \code{\link{logit_fit}}
#' @param X Covariate matrix of training sample
#' @param Y Vector of outcomes of training sample
#' @param Xnew Covariate matrix of test sample
#' @param weights Always FALSE as
#'
#' @return Returns list containing:
#' \item{prediction}{vector of predictions for Xnew}
#' \item{weights}{Not supported for propensity score estimation}
#'
#' @keywords internal
#'
predict.logit_fit = function(logit_fit,X,Y,Xnew=NULL,weights=FALSE){
  if (is.null(Xnew)) Xnew = X
  
  if (weights==TRUE) {
    warning("Weights are not supported for propensity score estimation.")
  }
  
  lambda = min(logit_fit$lambda)
  fit <- as.data.frame(predict(logit_fit, newx = as.matrix(Xnew), s = lambda, type = "response"))
  
  if (length(logit_fit$classnames)==2){
    fit[,2] = 1 - fit[,1]
    colnames(fit) = rev(logit_fit$glmnet.fit$classnames)
      } else {
    colnames(fit) = logit_fit$classnames
  }
  
  return(fit)
}


#' This function calculates the (multinomial) logistic regression based on the \code{\link{nnet}} package 
#'
#' @param X Covariate matrix of training sample
#' @param Y Vector of outcomes of training sample
#' @param arguments List of arguments passed to  \code{\link{nnet}}
#' @import nnet 
#' @importFrom utils modifyList
#'
#' @return An object with S3 method for class 'nnet'
#'
logit_nnet_fit = function(X,Y,arguments=list()){
  data <- data.frame(Y = as.factor(Y), X)
  arguments <- utils::modifyList(list(trace = FALSE), arguments)
  
  model = do.call(nnet::multinom, c(list(data=data, formula = Y ~., model=TRUE), arguments))
  model
}

#' Prediction based on (multinomial) logistic regression.
#' @param multinom_fit Output of \code{\link{nnet}}
#' @param X Covariate matrix of training sample
#' @param Y Vector of outcomes of training sample
#' @param Xnew Covariate matrix of test sample
#' @param weights Always FALSE as
#'
#' @return Returns list containing:
#' \item{prediction}{matrix of predictions for Xnew}
#' \item{weights}{Not supported for propensity score estimation}
#'
#' @keywords internal
#'
predict.logit_nnet_fit = function(multinom_fit, X, Y, Xnew=NULL, weights=FALSE){
  if (is.null(Xnew)) Xnew = X
  if (weights==TRUE) {
    warning("Weights are not supported for propensity score estimation.")
  }
  
  data = as.data.frame(Xnew)
  data$Y <- as.factor(0)
  
  fit <- as.data.frame(predict(multinom_fit, newdata = Xnew, type="probs"))
  
  if (length(unique(Y))==2){
    fit[,2] = 1-fit[,1]
    colnames(fit) = rev(sort(unique(Y)))
  }
  return(fit)
}


#' This function calculates the Gaussian Naive Bayes model based on the \code{\link{naivebayes}} package 
#'
#' @param X Covariate matrix of training sample
#' @param Y Vector of outcomes of training sample
#' @param arguments List of arguments passed to  \code{\link{naivebayes}}
#' @import naivebayes
#'
#' @return An object with S3 method for class 'gaussian_naive_bayes'
#'
nb_gaussian_fit = function(X,Y,arguments=list()){
  Y = as.factor(Y)
  model = do.call(naivebayes::gaussian_naive_bayes, c(list(x=X,y=Y),arguments))
  model
}

#' Prediction based on Gaussian Naive Bayes model.
#' @param nb_gaussian_fit Output of \code{\link{gaussian_naive_bayes}} or \code{\link{nb_gaussian_fit}}
#' @param X Covariate matrix of training sample
#' @param Y Vector of outcomes of training sample
#' @param Xnew Covariate matrix of test sample
#' @param weights Always FALSE as
#'
#' @return Returns list containing:
#' \item{prediction}{vector of predictions for Xnew}
#' \item{weights}{Not supported for propensity score estimation}
#'
#' @keywords internal
#'
predict.nb_gaussian_fit = function(nb_gaussian_fit,X,Y,Xnew=NULL,weights=FALSE){
  if (is.null(Xnew)) Xnew = X
  if (weights==TRUE) {
    warning("Weights are not supported for propensity score estimation.")
  }
  
  fit = predict(nb_gaussian_fit, newdata=Xnew, type = "prob")
  return(fit)
}

#' This function calculates the bernoulli Naive Bayes model based on the \code{\link{naivebayes}} package 
#'
#' @param X Covariate matrix of training sample
#' @param Y Vector of outcomes of training sample
#' @param arguments List of arguments passed to  \code{\link{naivebayes}}
#' @import naivebayes
#'
#' @return An object with S3 method for class 'bernoulli_naive_bayes'
#'
nb_bernoulli_fit = function(X,Y,arguments=list()){
  Y = as.factor(Y)
  model = do.call(naivebayes::naive_bayes, c(list(x=X,y=Y),arguments))
  model
}

#' Prediction based on bernoulli Naive Bayes model.
#' @param nb_bernoulli_fit Output of \code{\link{naive_bayes}} or \code{\link{nb_bernoulli_fit}}
#' @param X Covariate matrix of training sample
#' @param Y Vector of outcomes of training sample
#' @param Xnew Covariate matrix of test sample
#' @param weights Always FALSE as
#'
#' @return Returns list containing:
#' \item{prediction}{vector of predictions for Xnew}
#' \item{weights}{Not supported for propensity score estimation}
#'
#' @keywords internal
#'
predict.nb_bernoulli_fit = function(nb_bernoulli_fit,X,Y,Xnew=NULL,weights=FALSE){
  if (is.null(Xnew)) Xnew = X
  if (weights==TRUE) {
    warning("Weights are not supported for propensity score estimation.")
  }
  
  fit = predict(nb_bernoulli_fit, newdata=Xnew, type = "prob")
  return(fit)
}


#' This function calculates the xgboost model based on the \code{\link{xgboost}} package 
#'
#' @param X Covariate matrix of training sample
#' @param Y Vector of outcomes of training sample
#' @param arguments List of arguments passed to  \code{\link{xgboost}}
#' @import xgboost
#'
#' @return An object with S3 method for class 'xgb.Booster'
#'

# @Maren: seems not to work with the do.call() function
xgboost_prop_fit = function(X,Y,arguments=list(nrounds=40, eta=0.01)){
  K = length(unique(Y))
  if (K==2){
     model = do.call(xgboost::xgboost, c(list(data=X,label=Y, verbose=0,
                                     objective = "binary:logistic", eval_metric = "logloss"),arguments))
  }else{
    if (min(Y)==1) Y <- Y - 1
    model = do.call(xgboost::xgboost, c(list(data=X,label=Y, num_class = K, verbose=0,
                                    objective = "multi:softprob", eval_metric = "mlogloss"),arguments))
  }
  
  model
}


#' Prediction based on xgboost model.
#' @param xgboost_fit Output of \code{\link{xgboost}} or \code{\link{xgboost_fit}}
#' @param X Covariate matrix of training sample
#' @param Y Vector of outcomes of training sample
#' @param Xnew Covariate matrix of test sample
#' @param weights Always FALSE as
#'
#' @return Returns list containing:
#' \item{prediction}{vector of predictions for Xnew}
#' \item{weights}{Not supported for propensity score estimation}
#'
#' @keywords internal
#'
predict.xgboost_prop_fit = function(xgboost_fit,X,Y,Xnew=NULL,weights=FALSE){
  if (is.null(Xnew)) Xnew = X
  if (weights==TRUE) {
    warning("Weights are not supported for propensity score estimation.")
  }
  
  fit = predict(xgboost_fit, newdata=Xnew, type = "prob") 
  if (length(unique(Y))==2){
    fit = as.data.frame(fit)
    fit[,2] = 1 - fit[,1]
    colnames(fit) = rev(sort(unique(Y)))
  }else{
    fit = matrix(fit, nrow = nrow(Xnew), ncol = length(unique(Y)), byrow = TRUE)
    colnames(fit) = sort(unique(Y))
    
  }
  return(fit)
}


#' This function calculates the Support Vector Machines model based on the \code{e1071} package 
#'
#' @param X Covariate matrix of training sample
#' @param Y Vector of outcomes of training sample
#' @param arguments List of arguments passed to  \code{\link{svm}}
#' @import e1071
#'
#' @return An object with S3 method for class 'svm'
#'

svm_fit = function(X,Y,arguments=list(gamma = 0.1)){
  model = do.call(e1071::svm, c(list(y = Y, x = X, probability = TRUE, kernel = "radial", 
                              type = 'C-classification'), arguments))
  return(model)
}

#' Prediction based on Support Vector Machines model.
#' @param svm_fit Output of \code{\link{svm_fit}}
#' @param X Covariate matrix of training sample
#' @param Y Vector of outcomes of training sample
#' @param Xnew Covariate matrix of test sample
#' @param weights Always FALSE as
#'
#' @return Returns list containing:
#' \item{prediction}{vector of predictions for Xnew}
#' \item{weights}{Not supported for propensity score estimation}
#'
#' @keywords internal
#'
predict.svm_fit <- function(svm_fit, X, Y, Xnew = NULL, weights = FALSE) {
  if (is.null(Xnew)) Xnew = X
  if (weights==TRUE) {
    warning("Weights are not supported for propensity score estimation.")
  }
  
  # Perform predictions using the trained SVM model
  pred <- predict(svm_fit, newdata = Xnew, probability = TRUE)
  fit <- attr(pred, "probabilities")
  
  return(fit)
}


#' Calculates Probability Forest fit using the \code{\link{grf}} package
#'
#' @param X Matrix of covariates
#' @param Y vector of outcomes
#' @param arguments List of arguments passed to  \code{\link{probability_forest}}
#' @import grf
#'
#' @return An object with S3 class "probability_forest"
#'
#' @keywords internal
#'
prob_forest_fit = function(X,Y,arguments=list(num.trees=1000)){
  model = do.call(grf::probability_forest, c(list(X = X, Y = as.factor(Y)), arguments))
  model
}

#' Prediction based on Probability Forest fit.
#' @param prob_forest_fit Output of \code{\link{prob_forest_fit}}
#' @param X Covariate matrix of training sample
#' @param Y Vector of outcomes of training sample
#' @param Xnew Covariate matrix of test sample
#' @param weights Always FALSE
#'
#' @return Returns list containing:
#' \item{prediction}{vector of predictions for Xnew}
#' \item{weights}{Not supported for propensity score estimation}
#'
#' @keywords internal
#'
predict.prob_forest_fit = function(prob_forest_fit, X,Y,Xnew=NULL,weights=FALSE){
  if (is.null(Xnew)) Xnew = X
  
  if (weights) { w <- grf::get_forest_weights(forest_grf_fit, newdata = Xnew) } else { w <- NULL }
  
  fit = predict(prob_forest_fit, newdata=Xnew)$predictions
  return(fit)
}


#' Calculates k-Nearest Neighbor model using the \code{\link{kknn}} package
#'
#' @param X Matrix of covariates
#' @param Y vector of outcomes
#' @param arguments List of arguments passed to  \code{\link{kknn}}
#' @import kknn
#'
#' @return An object with S3 class "kknn"
#'
#' @keywords internal
#'
knn_prop_fit = function(X,Y,arguments=list(kmax=floor(0.05*nrow(X)))){ 
  data <- data.frame(Y = as.factor(Y), X)
  model = do.call(kknn::train.kknn ,c(list(formula = Y ~ ., data = data), arguments))
  model
}

#' Prediction based on k-Nearest Neighbor model.
#' @param knn_fit Output of \code{\link{kknn}} or \code{\link{knn_fit}}
#' @param X Covariate matrix of training sample
#' @param Y Vector of outcomes of training sample
#' @param Xnew Covariate matrix of test sample
#' @param weights Always FALSE
#'
#' @return Returns list containing:
#' \item{prediction}{vector of predictions for Xnew}
#' \item{weights}{Not supported for propensity score estimation}
#'
#' @keywords internal
#'
predict.knn_prop_fit = function(knn_fit, X,Y,Xnew=NULL,weights=FALSE){
  if (is.null(Xnew)) Xnew = X
  if (weights==TRUE) {
    warning("Weights are not supported for propensity score estimation.")
  }
  Xnew = as.data.frame(Xnew)
  fit = predict(knn_fit, newdata = Xnew, type='prob')
  return(fit)
}


#' Calculates Probability Forest fit using the \code{\link{ranger}} package
#'
#' @param X Matrix of covariates
#' @param Y vector of outcomes
#' @param arguments List of arguments passed to  \code{\link{ranger}}
#' @import ranger
#'
#' @return An object with S3 class "ranger"
#'
#' @keywords internal
#'
ranger_fit = function(X,Y,arguments=list(num.trees = 1000,
                                    min.node.size=0.1*nrow(X), 
                                    mtry=ceiling(sqrt(ncol(X))))){ #arguments=list()
  data <- data.frame(Y = as.factor(Y), X)

  model = do.call(ranger::ranger, c(list(data=data, formula= Y~., probability=TRUE), arguments))
  model
}


#' Prediction based on the \code{\link{ranger}} Probability Forest fit.
#' @param ranger_fit Output of \code{\link{ranger}} or \code{\link{prob_forest_fit}}
#' @param X Covariate matrix of training sample
#' @param Y Vector of outcomes of training sample
#' @param Xnew Covariate matrix of test sample
#' @param weights Always FALSE
#'
#' @return Returns list containing:
#' \item{prediction}{vector of predictions for Xnew}
#' \item{weights}{Not supported for propensity score estimation}
#'
#' @keywords internal
#'
predict.ranger_fit = function(ranger_fit, X,Y,Xnew=NULL,weights=FALSE){
  if (is.null(Xnew)) Xnew = X
  if (weights==TRUE) {
    warning("Weights are not supported for propensity score estimation.")
  }
  
  data = as.data.frame(Xnew)
  data$Y <- as.factor(0)
  
  fit <- predict(ranger_fit, data = data)$predictions 
  return(fit)
}


#####
#' Fits One-Vs-One (OvO) binary classifiers for multi-class classification
#' with optional parallelization if foreach/doParallel are available
#'
#' @param X Covariate matrix of training sample
#' @param Y Vector of outcomes of training sample
#' @param method Method to use for fitting binary classifiers (e.g., "logit" or "xgboost")
#' @param parallel Logical. Try to run in parallel if packages are available. Defaults to TRUE.
#' @param verbose Logical. If TRUE, messages about parallel/sequential choice are printed.
#' 
#' @return List of fitted OvO binary classifiers
#'
ovo_fit <- function(X, Y, method = "logit", parallel = FALSE, verbose = TRUE) {
  
  if (min(Y) == 0) Y <- Y + 1
  
  class_labels <- sort(unique(Y))
  num_classes <- length(class_labels)
  
  # ---- Decide if parallel backend is available ----
  use_parallel <- FALSE
  if (parallel && requireNamespace("foreach", quietly = TRUE)) {
    
    # if no parallel backend registered, try to auto-register
    if (foreach::getDoParWorkers() <= 1 && requireNamespace("doParallel", quietly = TRUE)) {
      n_cores <- max(1, parallel::detectCores() - 1)
      cl <- parallel::makeCluster(n_cores)
      doParallel::registerDoParallel(cl)
      if (verbose) message("Registered parallel backend with ", n_cores, " cores.")
    }
    
    # re-check if backend > 1 worker
    use_parallel <- foreach::getDoParWorkers() > 1
  }
  
  if (use_parallel) {
    if (verbose) {
      message("[OvO Par]")
    }
    
    ovo_classifiers <- foreach::foreach(i = 1:(num_classes - 1), .combine = "c", .export = c(paste0(method, "_fit"))) %dopar% {
      classifiers <- list()
      for (j in (i + 1):num_classes) {
        class_i <- class_labels[i]
        class_j <- class_labels[j]
        
        subset_x <- X[Y %in% c(class_i, class_j), ]
        subset_y <- Y[Y %in% c(class_i, class_j)]
        
        if (method == "xgboost") {
          subset_y <- ifelse(subset_y == i, 1, 0)
        }
        
        classifiers[[paste(class_i, class_j, sep = "_")]] <- do.call(
          paste0(method, "_fit"),
          list(Y = subset_y, X = subset_x)
        )
      }
      classifiers
    }
  } else {
    if (parallel && verbose) {
      message("Parallel backend not available. Falling back to sequential OvO fitting.\n",
              "Install 'foreach' and 'doParallel' to enable parallelization.")
    }
    
    ovo_classifiers <- list()
    for (i in 1:(num_classes - 1)) {
      for (j in (i + 1):num_classes) {
        class_i <- class_labels[i]
        class_j <- class_labels[j]
        
        subset_x <- X[Y %in% c(class_i, class_j), ]
        subset_y <- Y[Y %in% c(class_i, class_j)]
        
        if (method == "xgboost") {
          subset_y <- ifelse(subset_y == i, 1, 0)
        }
        
        ovo_classifiers[[paste(class_i, class_j, sep = "_")]] <- do.call(
          paste0(method, "_fit"),
          list(Y = subset_y, X = subset_x)
          
        )
      }
    }
  }
  
  return(ovo_classifiers)
}


#' Predicts class probabilities using One-Vs-One (OvO) fitted classifiers 
#' with optional parallelization (auto-registers a backend if available)
#'
#' @param ovo_fit Output list of \code{\link{ovo_fit}}
#' @param X Covariate matrix of training sample
#' @param Y Vector of outcomes of training sample
#' @param Xnew Covariate matrix of test sample
#' @param weights Flag indicating whether weights are used (unsupported for propensity score estimation)
#' @param method Method used for fitting the binary classifiers (e.g., "logit" or "xgboost")
#' @param parallel Logical. Try to run in parallel if packages are available. Defaults to TRUE.
#' @param verbose Logical. If TRUE, messages about parallel/sequential choice are printed.
#'
#' @return Matrix of class probability predictions for \code{Xnew}
#'
#' @keywords internal
#'
predict.ovo_fit <- function(ovo_fit, X, Y, Xnew = NULL, 
                            weights = FALSE, method = "logit", 
                            parallel = FALSE, verbose = TRUE) {
  
  if (is.null(Xnew)) Xnew <- X
  if (weights) {
    warning("Weights are not supported for propensity score estimation.")
  }
  
  n_classifiers <- length(ovo_fit)
  n_samples <- nrow(Xnew)
  n_classes <- length(unique(unlist(strsplit(names(ovo_fit), "_"))))
  
  if (min(Y) == 0) Y <- Y + 1
  
  # ---- Build q-matrix tensor ----
  q_matrix_tensor <- array(NA, dim = c(n_samples, n_classes, n_classes))
  
  for (i in seq_len(n_classifiers)) {
    class_names <- unlist(strsplit(names(ovo_fit)[[i]], "_"))
    class_i <- class_names[1]
    class_j <- class_names[2]
    
    subset_y <- Y[Y %in% c(class_i, class_j)]
    
    # Predict probabilities using the i-th binary classifier
    fit_raw <- do.call(paste0("predict.", method, "_fit"), 
                       list(ovo_fit[[i]], X = X, Y = subset_y, Xnew = Xnew))
    
    if (is.list(fit_raw) && !is.data.frame(fit_raw)) fit_raw <- fit_raw$prediction
    
    fit_raw <- as.data.frame(fit_raw, check.names = FALSE)
    
    if (ncol(fit_raw) == 1) {
      fit_raw[, 2] <- 1 - fit_raw[, 1]
    }
    if (!(all(class_names %in% colnames(fit_raw)))) {
      colnames(fit_raw) <- class_names
    }
    
    q_matrix_tensor[, as.numeric(class_i), as.numeric(class_j)] <- fit_raw[[class_i]]
    q_matrix_tensor[, as.numeric(class_j), as.numeric(class_i)] <- fit_raw[[class_j]]
  }
  
  # ---- Parallel decision ----
  use_parallel <- FALSE
  if (parallel && requireNamespace("foreach", quietly = TRUE)) {
    if (foreach::getDoParWorkers() <= 1 && requireNamespace("doParallel", quietly = TRUE)) {
      n_cores <- max(1, parallel::detectCores() - 1)
      cl <- parallel::makeCluster(n_cores)
      doParallel::registerDoParallel(cl)
      if (verbose) message("Registered parallel backend with ", n_cores, " cores.")
    }
    use_parallel <- foreach::getDoParWorkers() > 1
  }
  
  # ---- Optimization ----
  if (use_parallel) {
    if (verbose) message("[OvO Par]")
    opt_results <- foreach::foreach(row = 1:n_samples, .combine = "rbind", .export = "kl_convergence") %dopar% {
                                      opt_result <- optim(rep(1/n_classes, n_classes), 
                                                          kl_convergence, 
                                                          q_matrix = q_matrix_tensor[row, , ],
                                                          method = "L-BFGS-B", 
                                                          lower = rep(0, n_classes), upper = rep(1, n_classes))
                                      opt_result$par / sum(opt_result$par)
                                    }
    fit <- opt_results
  } else {
    if (parallel && verbose) {
      message("Parallel backend not available. Falling back to sequential prediction.")
    }
    opt_results <- apply(q_matrix_tensor, MARGIN = 1, function(q_matrix_row) {
      opt_result <- optim(rep(1/n_classes, n_classes), 
                          kl_convergence, 
                          q_matrix = matrix(q_matrix_row, nrow = n_classes),
                          method = "L-BFGS-B", 
                          lower = rep(0, n_classes), upper = rep(1, n_classes))
      opt_result$par / sum(opt_result$par)
    })
    fit <- t(simplify2array(opt_results))
  }
  
  colnames(fit) <- unique(unlist(strsplit(names(ovo_fit), "_")))
  return(fit)
}


#' Fits One-Vs-Rest (OvR) binary classifiers for multi-class classification 
#' @param X Covariate matrix of training sample
#' @param Y Vector of outcomes of training sample
#' @param method Method to use for fitting binary classifiers (e.g., "logit" or "xgboost")
#' 
#' @return List of fitted OvR binary classifiers
#'
ovr_fit <- function(X, Y, method = "logit") {
  class_labels <- sort(unique(Y))
  num_classes <- length(class_labels)
  ovr_classifiers <- list()
  
  if (min(Y)==0) Y = Y+1
  
  # Create binary classifiers for each pair of classes
  for (i in 1:(num_classes)) {
      class_i <- class_labels[i]      
      binarized_y <- ifelse(Y == i, 1, 0)
  
      ovr_classifiers[[paste(class_i)]] <- do.call(
        paste0(method, "_fit"),
        list(Y = binarized_y, X = X)
      )
  }
  
  return(ovr_classifiers)
}


#' Predicts class probabilities using One-Vs-One (OvR) fitted classifiers 
#' @param ovr_fit Output list of \code{\link{ovr_fit}} 
#' @param X Covariate matrix of training sample
#' @param Y Vector of outcomes of training sample
#' @param Xnew Covariate matrix of test sample
#' @param weights Flag indicating whether weights are used (unsupported for propensity score estimation)
#' @param method Method used for fitting the binary classifiers (e.g., "logit" or "xgboost")
#'
#' @return Returns list containing:
#' \item{prediction}{matrix of predictions for Xnew}
#' \item{weights}{Not supported for propensity score estimation}
#'
#' @keywords internal
#'
predict.ovr_fit <- function(ovr_fit,X,Y, Xnew=NULL,weights=FALSE, method = "logit") {
  if (is.null(Xnew)) Xnew = X
  if (weights==TRUE) {
    warning("Weights are not supported for propensity score estimation.")
  }
  
  n_classifiers <- length(ovr_fit)
  n_samples <- nrow(Xnew)
  n_classes <- n_classifiers
  
  if (min(Y)==0) Y = Y+1
  
  # Initialize an empty matrix to store the probabilities
  fit <- as.data.frame(matrix(0, nrow = n_samples, ncol = n_classes))
  colnames(fit) <- seq.int(1, n_classes)
  
  for (i in seq_len(n_classifiers)) {
    binarized_y <- ifelse(Y == i, 1, 0)
    
    # Predict probabilities using the i-th binary classifier
    fit_raw <- do.call(
      paste0("predict.", method, "_fit"), 
      list(ovr_fit[[i]], X = X, Y = binarized_y, Xnew = Xnew)
    )
    
    # If fit_raw is a list, extract prediction if it's not already a data frame
    if (is.list(fit_raw) && !is.data.frame(fit_raw)) {fit_raw <- fit_raw$prediction}
    
    fit_raw <- as.data.frame(fit_raw)
    
    # If two columns, keep only the column named "1"
    if (ncol(fit_raw) == 2) {fit_raw <- fit_raw[, 1, drop = FALSE]}
    
    # Update the corresponding column in e_hat
    fit[, i] <- fit[, i] + fit_raw[, 1]
  }
  
  # Normalize the probabilities to sum up to 1 for each sample
  fit <- fit / rowSums(fit)
  
  return(fit)
}


#' Computes KL divergence-based convergence measure
#' @param p Probability vector representing a probability distribution
#' @param q_matrix Square matrix representing transition probabilities
#' 
#' @return Returns the KL divergence-based convergence measure
#'
kl_convergence <- function(p, q_matrix) {
  epsilon <- 1e-7
  p <- p / sum(p) # Normalize p to sum to 1
  K <- length(p)
  total_divergence <- 0
  
  for (i in 1:(K-1)) {
    for (j in (i+1):K) {
      u_ij <- p[i] / (p[i] + p[j] + epsilon)
      r_ij <- q_matrix[i, j]
      
      if (!is.na(r_ij)) {
        total_divergence <- total_divergence + 
          r_ij * log((r_ij + epsilon) / (u_ij + epsilon)) + 
          (1 - r_ij) * log(((1 - r_ij) + epsilon) / ((1 - u_ij) + epsilon)) 
      }
    }
  }
  
  return(total_divergence)
}