#' Fits One-Vs-One (OvO) binary classifiers 
#' 
#' Fits One-Vs-One (OvO) binary classifiers for multi-class classification
#' with optional parallelization if \code{foreach} and \code{doParallel} packages are available.
#'
#' @param X Covariate matrix of training sample.
#' @param Y Vector of outcomes of training sample.
#' @param method Method to use for fitting binary classifiers (e.g., "logit" or "xgboost").
#' @param parallel Logical. Tries to run in parallel if \code{foreach} and 
#'  \code{doParallel} are available. Defaults to FALSE.
#' @param quiet Logical. If FALSE, messages about parallel/sequential choice are printed.
#'
#' @return List of fitted OvO binary classifiers.
#'
#' @keywords internal
ovo_fit <- function(X, Y, method, parallel = FALSE, quiet = TRUE) {
  if (min(Y) == 0) Y <- Y + 1

  class_labels <- sort(unique(Y))
  num_classes  <- length(class_labels)
  use_parallel <- FALSE
  
  # Is parallel backend available? 
  if (parallel && requireNamespace("foreach", quietly = TRUE)) {
    # if no parallel backend registered, try to auto-register
    if (foreach::getDoParWorkers() <= 1 && requireNamespace("doParallel", quietly = TRUE)) {
      n_cores <- max(1, parallel::detectCores() - 1)
      cl <- parallel::makeCluster(n_cores)
      doParallel::registerDoParallel(cl)
      message("Registered parallel backend with ", n_cores, " cores.")
    }
    use_parallel <- foreach::getDoParWorkers() > 1
  }

  if (use_parallel) {
    if (!quiet) message("[OvO Par]")
    `%dopar%` <- foreach::`%dopar%`

    ovo_classifiers <- foreach::foreach(i = 1:(num_classes - 1), .combine = "c", 
                                        .export = c(paste0(method, "_fit"))) %dopar% {
      classifiers <- list()
      for (j in (i + 1):num_classes) {
        class_i <- class_labels[i]; class_j <- class_labels[j]
        idx <- Y %in% c(class_i, class_j)
        subset_x <- X[idx, ]
        subset_y <- Y[idx]

        if (method == "xgboost_prop") {
          subset_y <- ifelse(subset_y == i, 1, 0)
        }

        classifiers[[paste(class_i, class_j, sep = "_")]] <- do.call(
          paste0(method, "_fit"), list(Y = subset_y, X = subset_x))
      }
      classifiers
    }
  } else {
    if (parallel) message("Parallel backend not available. Falling back to sequential OvO fitting.")

    ovo_classifiers <- list()
    for (i in 1:(num_classes - 1)) {
      for (j in (i + 1):num_classes) {
        class_i <- class_labels[i]; class_j <- class_labels[j]
        idx <- Y %in% c(class_i, class_j)
        subset_x <- X[idx, ]
        subset_y <- Y[idx]

        if (method == "xgboost_prop") {
          subset_y <- ifelse(subset_y == i, 1, 0)
        }

        ovo_classifiers[[paste(class_i, class_j, sep = "_")]] <- do.call(
          paste0(method, "_fit"), list(Y = subset_y, X = subset_x))
      }
    }
  }

  return(ovo_classifiers)
}


#' Predicts class probabilities using One-Vs-One (OvO) routine
#' 
#' Predicts class probabilities using One-Vs-One (OvO) fitted classifiers
#' with optional parallelization (auto-registers a backend if \code{foreach} 
#' and \code{doParallel} packages are available).
#'
#' @param object Output list of \code{\link{ovo_fit}}.
#' @param X Covariate matrix of training sample.
#' @param Y Vector of outcomes of training sample.
#' @param Xnew Covariate matrix of test sample.
#' @param method Method used for fitting the binary classifiers (e.g., "logit").
#' @param parallel Logical. Tries to run in parallel if \code{foreach} and 
#'  \code{doParallel} packages are available. Defaults to FALSE.
#' @param quiet Logical. If FALSE, messages about parallel/sequential choice are printed.
#' @param ... Ignored additional arguments.
#'
#' @return A matrix of predicted class probabilities for \code{Xnew} with one 
#' column per class, or a vector of positive-class probabilities for binary classification. 
#'
#' @keywords internal
predict.ovo_fit <- function(object, X, Y, Xnew = NULL, method,
                            parallel = FALSE, quiet = TRUE, ...) {
  if (is.null(Xnew)) Xnew <- X
  use_parallel <- FALSE

  n_classifiers <- length(object)
  n_samples <- nrow(Xnew)
  n_classes <- length(unique(unlist(strsplit(names(object), "_"))))

  if (min(Y) == 0) Y <- Y + 1

  # Build q-matrix tensor
  q_matrix_tensor <- array(NA, dim = c(n_samples, n_classes, n_classes))

  for (i in seq_len(n_classifiers)) {
    class_names <- unlist(strsplit(names(object)[[i]], "_"))
    class_i <- class_names[1]
    class_j <- class_names[2]
    subset_y <- Y[Y %in% c(class_i, class_j)]

    # Predict probabilities using the i-th binary classifier
    fit_raw <- do.call(
      paste0("predict.", method, "_fit"),
      list(object[[i]], X = X, Y = subset_y, Xnew = Xnew)
    )

    # Ensure the proper data format (need both classes)
    fit_raw <- data.frame(1 - fit_raw, fit_raw, check.names = FALSE)
    if (!(all(class_names %in% colnames(fit_raw)))) colnames(fit_raw) <- class_names

    q_matrix_tensor[, as.numeric(class_i), as.numeric(class_j)] <- fit_raw[[class_i]]
    q_matrix_tensor[, as.numeric(class_j), as.numeric(class_i)] <- fit_raw[[class_j]]
  }

  # Parallel option
  if (parallel && requireNamespace("foreach", quietly = TRUE)) {
    if (foreach::getDoParWorkers() <= 1 && requireNamespace("doParallel", quietly = TRUE)) {
      n_cores <- max(1, parallel::detectCores() - 1)
      cl <- parallel::makeCluster(n_cores)
      doParallel::registerDoParallel(cl)
      message("Registered parallel backend with ", n_cores, " cores.")
    }
    use_parallel <- foreach::getDoParWorkers() > 1
  }

  # ---- Optimization ----
  if (use_parallel) {
    if (!quiet) message("[OvO Par]")
    `%dopar%` <- foreach::`%dopar%`
    
    opt_results <- foreach::foreach(row = 1:n_samples, .combine = "rbind", .export = "kl_convergence") %dopar% {
      opt_result <- stats::optim(rep(1 / n_classes, n_classes),
        kl_convergence,
        q_matrix = q_matrix_tensor[row, , ],
        method = "L-BFGS-B",
        lower = rep(0, n_classes), upper = rep(1, n_classes)
      )
      opt_result$par / sum(opt_result$par)
    }
    fit <- opt_results
  } else {
    if (parallel) message("Parallel backend not available. Falling back to sequential prediction.")

    opt_results <- apply(q_matrix_tensor, MARGIN = 1, function(q_matrix_row) {
      opt_result <- stats::optim(rep(1 / n_classes, n_classes),
        kl_convergence,
        q_matrix = matrix(q_matrix_row, nrow = n_classes),
        method = "L-BFGS-B",
        lower = rep(0, n_classes), upper = rep(1, n_classes)
      )
      opt_result$par / sum(opt_result$par)
    })
    fit <- t(simplify2array(opt_results))
  }
  
  colnames(fit) <- unique(unlist(strsplit(names(object), "_")))
  if (ncol(fit) == 2) fit <- fit[, 2]

  return(fit)
}


#' Fit One-vs-Rest (OvR) classifiers
#' 
#' Fits a collection of binary classifiers in a one-vs-rest scheme for 
#' multi-class classification. Each class is modeled against all others.  
#'
#' @param X Covariate matrix of training sample.
#' @param Y Vector of outcomes of training sample.
#' @param method Character string specifying the binary classifier to use 
#'   (e.g., \code{"logit"}, \code{"xgboost"}).
#'
#' @return A named list of fitted binary classifiers, one for each class.
#'
#' @keywords internal
ovr_fit <- function(X, Y, method) {
  class_labels <- sort(unique(Y))
  num_classes <- length(class_labels)
  ovr_classifiers <- list()

  if (min(Y) == 0) Y <- Y + 1

  # Create binary classifiers for each pair of classes
  for (i in 1:(num_classes)) {
    class_i <- class_labels[i]
    binarized_y <- ifelse(Y == i, 1, 0)

    ovr_classifiers[[paste(class_i)]] <- do.call(
      paste0(method, "_fit"), list(Y = binarized_y, X = X))
  }

  return(ovr_classifiers)
}


#' Predict method for One-vs-Rest (OvR) classifiers
#' 
#' Generates predicted class probabilities using a set of fitted OvR binary classifiers.  
#'
#' @param object Output list from \code{\link{ovr_fit}}.
#' @param X Covariate matrix of training sample.
#' @param Y Vector of outcomes of training sample.
#' @param Xnew Covariate matrix of test sample. If \code{NULL}, uses training data.
#' @param method Character string specifying the binary classifier used for fitting.
#' @param ... Ignored additional arguments.
#'
#' @return A matrix of predicted class probabilities with one column per class, 
#'   or a vector of positive-class probabilities for binary classification.  
#'
#' @keywords internal
predict.ovr_fit <- function(object, X, Y, Xnew = NULL, method, ...) {
  if (is.null(Xnew)) Xnew <- X
  if (min(Y) == 0) Y <- Y + 1
  
  n_classifiers <- length(object)
  n_samples <- nrow(Xnew)
  n_classes <- n_classifiers
  
  fit <- matrix(0, nrow = n_samples, ncol = n_classes)
  colnames(fit) <- seq.int(1, n_classes)
  
  for (i in seq_len(n_classifiers)) {
    binarized_y <- ifelse(Y == i, 1, 0)
    
    # Predict probabilities using the i-th binary classifier
    fit_raw <- do.call(
      paste0("predict.", method, "_fit"),
      list(object[[i]], X = X, Y = binarized_y, Xnew = Xnew)
    )
    # Update the corresponding column in fit
    fit[, as.character(i)] <- fit_raw
  }
  fit <- fit / rowSums(fit)
  if (ncol(fit) == 2) fit <- fit[, 2]
  
  return(fit)
}


#' KL divergence-based convergence measure
#' 
#' Computes a convergence measure between an empirical distribution and 
#' a transition probability matrix using KL divergence.  
#'
#' @param p Numeric probability vector representing a distribution.
#' @param q_matrix Square matrix of transition probabilities.
#'
#' @return A numeric value representing the KL divergence-based convergence measure.  
#'
#' @keywords internal
kl_convergence <- function(p, q_matrix) {
  epsilon <- 1e-7
  p <- p / sum(p)
  K <- length(p)
  total_divergence <- 0

  for (i in 1:(K - 1)) {
    for (j in (i + 1):K) {
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
