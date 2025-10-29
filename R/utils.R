#' One-hot encoding of treatment/instrument vector
#'
#' @description
#' \code{\link{prep_indicator_mat}} creates a logical matrix of binary indicators
#' (N x C) where each column represents a unique category from the input vector.
#' Useful for encoding treatments/instruments.
#'
#' @param x Treatment (instrument) vector.
#'          Provide as factor to control category ordering, otherwise arranged
#'          in ascending order or alphabetically.
#'
#' @return Logical matrix (N x C) where C is the number of unique categories in `x`.
#'         Each column is a binary indicator for one category.
#'
#' @examples
#' \dontrun{
#' D <- factor(c("A", "B", "A", "C", "B", "C"))
#' d_mat <- prep_indicator_mat(D)
#' head(d_mat)
#' }
#'
#' @keywords internal
prep_indicator_mat <- function(x) {
  if (length(unique(x)) <= 1) {
    stop("Input vector must contain at least two unique values.")
  }
  if (!is.factor(x)) {
    x <- factor(x, levels = sort(unique(x)))
  }

  # One-hot encoding
  mat <- stats::model.matrix(~ 0 + x)
  colnames(mat) <- gsub("x", "", colnames(mat))
  return(mat == 1L)
}


#' Compatibility check of cluster vector
#'
#' @description
#' \code{\link{check_cluster_compatibility}} checks if the cross-fitting procedure
#' is feasible given the cluster vector and desired number of cross-fitting folds.
#'
#' @param cluster A vector representing the cluster assignments
#' @param cf The desired number of folds for cross-fitting
#'
#' @return Does not return anything; throws an error if checks fail.
#'
#' @keywords internal
#'
check_cluster_compatibility <- function(cluster, cf) {
  num_clusters <- length(unique(cluster))

  if (num_clusters < cf) {
    stop(
      "The number of clusters is less than the desired number of folds.
      Either choose a smaller number of folds or do not specify a cluster vector."
    )
  }

  cluster_shares <- table(cluster) / length(cluster)
  max_cluster_share <- max(cluster_shares)

  if (max_cluster_share > (1 / cf) * 0.9) {
    stop(
      "There is a high imbalance in the cluster sizes. This poses a problem for the cross-fitting procedure.
      Either choose a smaller number of folds or do not specify a cluster vector."
    )
  }
}


#' Cross-fitting fold indicators
#'
#' Creates a matrix of binary cross-fitting fold indicators (\eqn{N \times cf})
#'
#' @param N Number of observations
#' @param cf Number of cross-fitting folds
#' @param cluster Optional vector of cluster variable if cross-fitting should account
#'                for clusters within the data.
#' @param d_mat Optional logical matrix of treatment indicators (\eqn{N x T+1} with \eqn{T}
#'              being the number of treatments). For example created by \code{\link{prep_indicator_mat}}.
#'              If specified, cross-fitting folds will preserve the treatment ratios from full sample.
#'              However, if cluster vector is provided, \code{d_mat} is ignored due to computational
#'              constraints and randomization as well as feasibility issues.
#'
#' @return Logical matrix of cross-fitting folds (N x # folds).
#'
#' @examples
#' \dontrun{
#' N <- 1000
#' d_mat <- prep_indicator_mat(x = sample(3, N, replace = TRUE))
#' cf_mat <- prep_cf_matrix(N = N, cf = 5, d_mat = d_mat)
#' head(cf_mat)
#' }
#'
#' @keywords internal
prep_cf_matrix <- function(N,
                        cf,
                        cluster = NULL,
                        d_mat = NULL) {
  if (!is.null(cluster) & !is.null(d_mat)) {
    warning("Both cluster vector and treatment matrix were provided. 
            Only the cluster vector will be used.")
    d_mat <- NULL
  }
  
  # Only one fold (i.e. no cross-fitting)
  if (cf == 1) {
    cf_mat <- matrix(rep(TRUE, N), ncol = 1)
    
    # Neither treatment matrix nor cluster vector provided
  } else if (is.null(d_mat) & is.null(cluster)) {
    rnd_id <- sample(1:N, N)
    fold <- factor(as.numeric(cut(rnd_id, breaks = stats::quantile(rnd_id, probs = seq(0, 1, length = cf + 1)), include.lowest = TRUE)))
    cf_mat <- (stats::model.matrix(~ 0 + fold) == 1)
    
    # Treatment matrix but no cluster vector provided
  } else if (!is.null(d_mat) & is.null(cluster)) {
    cf_mat <- matrix(NA, nrow = N, ncol = cf)
    nw <- colSums(d_mat)
    
    for (i in 1:ncol(d_mat)) {
      cf_mat_w <- matrix(FALSE, nrow = nw[i], ncol = cf)
      rnd_id <- sample(1:nw[i], nw[i])
      fold <- as.numeric(cut(rnd_id, breaks = stats::quantile(rnd_id, probs = seq(0, 1, length = cf + 1)), include.lowest = TRUE))
      
      for (j in 1:cf) {
        cf_mat_w[fold == j, j] <- TRUE
      }
      
      cf_mat[d_mat[, i], ] <- cf_mat_w
    }
    
    # No treatment matrix but cluster vector provided
  } else if (is.null(d_mat) & !is.null(cluster)) {
    check_cluster_compatibility(cluster, cf)
    
    rnd_id <- sample(1:length(unique(cluster)), length(unique(cluster)))
    fold <- as.numeric(cut(rnd_id, breaks = stats::quantile(rnd_id, probs = seq(0, 1, length = cf + 1)), include.lowest = TRUE))
    fold <- factor(fold[match(cluster, unique(cluster))])
    cf_mat <- (stats::model.matrix(~ 0 + fold) == 1)
    
    cf_mat_balance <- colSums(cf_mat) / nrow(cf_mat)
    
    if (any(cf_mat_balance < (1 / cf) * 0.5)) {
      stop("High cluster size imbalance detected. Choose fewer folds or omit cluster vector.")
    }
  }
  
  colnames(cf_mat) <- sprintf("CF %d", 1:cf)
  
  return(cf_mat)
}


#' Estimate ensemble weights using various methods
#'
#' Estimates ensemble weights using non-negative least squares (NNLS), BFGS 
#' optimization, single best learner selection, ordinary least squares (OLS), 
#' and simple averaging. Supports binary/continuous and multinomial outcomes. 
#' Weights are normalized to sum to 1 for all methods.
#'
#' @param X For continuous/binary outcomes: a numeric matrix of predictions
#'   with rows as observations and columns as base learners.
#'   For multinomial outcomes: a 3D array of predictions with dimensions
#'   \eqn{N × K × M}, where \eqn{N} is the number of observations,
#'   \eqn{K} the number of classes, and \eqn{M} the number of learners.
#' @param Y A numeric vector of observed outcomes. For multinomial outcomes,
#'   a factor or integer vector of class labels.
#' @param subset Optional logical vector if only subset of data should be used.
#' @param is_mult Logical. Indicates whether the outcome is multinomial
#'   (\code{TRUE}) or binary/continuous (\code{FALSE}).
#' @param ensemble_type Method for calculating ensemble weights:
#'   \describe{
#'     \item{\code{"nnls"}}{Non-negative least squares; weights sum to 1 (default)}
#'     \item{\code{"bfgs"}}{BFGS optimization (for multivalued treatments only; falls back to \code{nnls} otherwise)}
#'     \item{\code{"singlebest"}}{Weight of 1 on the learner with the lowest RMSE}
#'     \item{\code{"ols"}}{Ordinary least squares regression weights (falls back to \code{nnls} for multivalued treatments)}
#'     \item{\code{"average"}}{Equal weights for all learners}
#'   }
#'
#' @return A numeric vector of ensemble weights.
#'
#' @examples
#' \dontrun{
#' X <- matrix(c(0.2, 0.3, 0.4, 0.5, 0.6, 0.7), ncol = 2)
#' Y <- c(0.25, 0.45, 0.65)
#' ens_w <- ens_weights_maker(X = X, Y = Y)
#' }
#' @keywords internal
ens_weights_maker <- function(X, Y, 
                              ensemble_type = "nnls",
                              subset = NULL, 
                              is_mult = FALSE
                              ) {
  # Checks
  if (is.null(subset)) subset <- rep(TRUE, length(Y))
  n_learners <- if (is_mult) dim(X)[3] else ncol(X)
  
  # Method validation
  if (ensemble_type == "bfgs" && !is_mult) {
    #message("BFGS only available for multinomial outcomes. Using 'nnls' instead.")
    ensemble_type <- "nnls"
  }
  if (ensemble_type == "ols" && is_mult) {
    #message("OLS not implemented for multinomial outcomes. Using 'nnls' instead.")
    ensemble_type <- "nnls"
  }
  
  # RMSE helper function
  get_rmse <- function(Y.hat, Y) {
    sqrt(mean((Y - Y.hat)^2))
  }
  
  if (ensemble_type == "singlebest") {
    # Find learner with lowest RMSE
    if (is_mult) {
      # 3D array case (multiclass)
      rmses <- apply(X[subset, , , drop = FALSE], 3, function(preds) {
        get_rmse(as.vector(preds), Y[subset])
      })
    } else {
      # 2D matrix case
      rmses <- apply(X[subset, , drop = FALSE], 2, function(preds) {
        get_rmse(preds, Y[subset])
      })
    }
    best_idx <- which.min(rmses)
    ens_w <- numeric(n_learners)
    ens_w[best_idx] <- 1
    
  } else if (ensemble_type == "average") {
    # Equal weights for all learners
    ens_w <- rep(1 / n_learners, n_learners)
    
  } else if (ensemble_type == "ols") {
    # OLS regression weights
    if (length(dim(X)) == 3) {
      Y_stack <- as.vector(one_hot(Y[subset]))
      X_stack <- apply(X[subset, , , drop = FALSE], 3, as.vector)
      ens_w <- qr.solve(a = X_stack, b = Y_stack)
    } else {
      ens_w <- qr.solve(a = as.matrix(X[subset, ]), b = Y[subset])
    }

  } else if (ensemble_type == "bfgs") {
    # Multinomial - BFGS optimization
    softmax <- function(z) { z <- z - max(z); exp(z) / sum(exp(z)) }
    one_hot_Y <- one_hot(Y = Y)
    
    obj_fun <- function(raw_w) {
      weights <- softmax(raw_w)
      ens_pred <- agg_array(a = X, w = weights)
      mean((one_hot_Y - ens_pred)^2)
    }
    
    bfgs <- stats::optim(
      par = rep(0, dim(X)[3]),
      fn = obj_fun,
      method = "BFGS",
      control = list(reltol = 1e-8)
    )
    ens_w <- softmax(bfgs$par)

  } else if (ensemble_type == "nnls") {
    # NNLS
    if (is_mult) {
      # Multinomial - stack Y (K*N vector) and flatten X (K*N × M matrix)
      Y_stack <- as.vector(one_hot(Y[subset]))
      X_stack <- apply(X[subset, , , drop = FALSE], 3, as.vector)
      ens_w <- nnls::nnls(A = as.matrix(X_stack), b = Y_stack)$x
    } else {
      # Binary/continuous outcome
      ens_w <- nnls::nnls(A = as.matrix(X[subset, ]), b = Y[subset])$x
    }
    # Handle degenerate all-zero weights
    if (sum(ens_w) == 0) ens_w <- rep(1 / n_learners, n_learners)
    ens_w <- ens_w / sum(ens_w)
  }
  
  names(ens_w) <- if (is_mult) dimnames(X)[[3]] else colnames(X)
  return(ens_w)
}


#' Adds an intercept to a matrix
#'
#' \code{\link{add_intercept}} adds an intercept to a matrix.
#'
#' @param mat Any matrix.
#'
#' @return Matrix with intercept.
#'
#' @keywords internal
add_intercept <- function(mat) {
  if (is.null(dim(mat))) {
    mat <- as.matrix(mat, ncol = 1)
  }
  if (all(mat[, 1] == 1)) {
    return(mat)
  } else {
    mat <- cbind(rep(1, nrow(mat)), mat)
    return(mat)
  }
}


#' Aggregate along an array dimension with weights
#'
#' Computes the weighted sum of an array along a specified dimension.
#' By default, aggregates along the second dimension.
#'
#' @param a A numeric array of at least 2 dimensions.
#' @param w A numeric vector of weights. Its length must match the size
#'   of the dimension being aggregated.
#' @param dim Integer. The dimension over which to aggregate
#'   (default is \code{2}).
#'
#' @return A numeric array with the same dimensions as \code{a},
#'   except with the specified dimension removed. For example, if
#'   \code{a} is \eqn{N \times M} and \code{dim = 2}, the result is
#'   a length-\eqn{N} vector. If \code{a} is \eqn{N \times M \times K},
#'   the result is an \eqn{N \times K} matrix.
#'
#' @keywords internal
agg_array <- function(a, w, dim = 2) {
  apply(a, c(1, dim), function(x) sum(x * w))
}


#' Prepare \code{cf_preds} matrix for cross-fitted ensemble predictions
#'
#' Creates a matrix or array to store cross-fitted ensemble predictions.
#'
#' @param methods List of methods built via \code{\link{create_method}} to be used in ensemble
#' @param N Number of observations
#' @param Y Outcome vector (used to determine number of classes for multinomial outcomes)
#'
#' @return A matrix of NAs with dimensions \code{N} x \code{length(methods)} for binary outcomes,
#'         or a 3D array with dimensions \code{N} x \code{K} x \code{length(methods)} for multinomial
#'         outcomes (where K is number of unique classes in Y)
#'
#' @keywords internal
#'
make_cf_preds <- function(methods, N, Y) {
  is_mult <- !is.null(methods[[1]]$multinomial)
  n_methods <- length(methods)
  method_names <- names(methods)

  if (is_mult) {
    n_classes <- length(unique(Y))
    cf_preds <- array(NA, dim = c(N, n_classes, n_methods))
    class_names <- if (!is.null(levels(Y))) levels(Y) else seq_len(n_classes) - 1
    dimnames(cf_preds) <- list(NULL, class_names, method_names)
  } else {
    cf_preds <- matrix(NA, nrow = N, ncol = n_methods)
    colnames(cf_preds) <- if (!is.null(method_names)) method_names else sprintf("methods%s", seq_len(n_methods))
  }

  return(cf_preds)
}


#' Short- or Standard-stacking message
#'
#' Prints a message into the console whether short- or standard-stacking is used.
#'
#' @param cv Cross-validation value for ensemble weights.
#'
#' @return Message in console. No object is returned.
#'
#' @keywords internal
which_stacking <- function(cv = 1) {
  if (cv == 1) {
    message("Short-stacking is used.")
  } else if (cv > 1) {
    message("Standard-stacking is used.")
  }
}


#' Update progress bar
#'
#' Updates progress bar with standardized formatting.
#'
#' @param pb Progress bar object from \code{progress} package
#' @param pb_np Current nuisance parameter (character)
#' @param pb_cf Current cross-fitting fold number (integer)
#' @param pb_cv Current cross-validation fold number (integer)
#' @param task Task description (character)
#' @param method Method name (character)
#'
#' @keywords internal
update_progress <- function(pb, pb_np, pb_cf, pb_cv, task, method) {
  if (is.null(pb)) {
    return(invisible(NULL))
  }

  # Center-align with padding
  format_center <- function(x, width) {
    x <- substr(as.character(x), 1, width)
    pad <- width - nchar(x)
    paste0(strrep(" ", ceiling(pad / 2)), x, strrep(" ", floor(pad / 2)))
  }

  # Tokens length
  pb$tick(tokens = list(
    nuisance = format_center(pb_np, 8),
    pb_cf    = format_center(pb_cf, 2),
    pb_cv    = format_center(pb_cv, 2),
    task     = format_center(task, 4),
    model    = format_center(method, 10)
  ))
}


#' Extract ensemble weights from the output of \code{nuisance_cf}
#'
#' @param np_cf A cross-fitted nuisance parameters object containing ensemble models.
#'   The object can have two possible structures:
#'   \enumerate{
#'     \item List of sublists, each containing \code{ens_w} elements
#'     \item Top-level \code{ens_w} element directly accessible
#'   }
#'
#' @return Ensemble weights as named numeric vectors. Returns a list of weights
#'   if multiple sublists are present, or a single named vector if \code{ens_w}
#'   is directly accessible.
#'
#' @keywords internal
extract_w <- function(np_cf) {
  if ("ens_w" %in% names(np_cf$models_nupa)) {
    # Structure 2: ens_w is top-level element
    np_cf$models_nupa$ens_w
  } else {
    # Structure 1: ens_w is within each sublist
    lapply(np_cf$models_nupa, function(x) x$ens_w)
  }
}


#' Format ensemble weights for short- or standard-stacking
#'
#' @param ens_weights List of ensemble weights to be formatted
#' @param cv Integer indicating number of cross-validation folds
#'
#' @return A data.frame (\code{cv = 1}) or list of data.frames (\code{cv > 1}) with formatted weights
#'
#' @keywords internal
format_weights <- function(ens_weights, cv) {
  # Short stacking: Single matrix with all learners
  if (cv == 1) {
    all_methods <- unique(unlist(lapply(ens_weights, names)))
    ens_weights_mat <- matrix(NA,
      nrow = length(all_methods),
      ncol = length(ens_weights),
      dimnames = list(all_methods, names(ens_weights))
    )

    for (param in names(ens_weights)) {
      ens_weights_mat[names(ens_weights[[param]]), param] <- ens_weights[[param]]
    }

    out <- as.data.frame(ens_weights_mat)
    out[is.na(out)] <- 0
    class(out) <- c("ens_weights_short", "data.frame")

    # Standard stacking: List of fold-specific weights
  } else {
    out <- list()

    for (NuPa in names(ens_weights)) {
      all_methods <- names(ens_weights[[NuPa]][[1]])

      fold_weights <- matrix(0,
        nrow = length(all_methods),
        ncol = length(ens_weights[[NuPa]]),
        dimnames = list(all_methods, NULL)
      )

      # Fill in weights for each fold
      for (i in seq_along(ens_weights[[NuPa]])) {
        fold_weights[, i] <- ens_weights[[NuPa]][[i]]
      }

      colnames(fold_weights) <- paste0("fold", seq_len(ncol(fold_weights)))
      out[[NuPa]] <- as.data.frame(fold_weights)
    }

    class(out) <- c("ens_weights_stand", "list")
  }

  return(out)
}


#' Generate a color palette for ensemble weights visuals
#'
#' Creates a consistent color palette for visualizing ensemble method weights.
#' Returns a predefined set of colors for up to 12 methods, and generates
#' additional colors using \code{colorRampPalette} when more methods are needed.
#'
#' @param n Number of colors to return (number of distinct methods)
#'
#' @return A character vector of hex color codes of length `n`
#'
#' @keywords internal
get_palette <- function(n) {
  palette <- c(
    "#FB8072", "#80B1D3", "#FFED6F", "#BEBADA",
    "#8DD3C7", "#FDB462", "#B3DE69", "#BC80BD",
    "#FCCDE5", "#D9D9D9", "#FFFFB3", "#8DA0CB"
  )
  if (n > length(palette)) {
    palette <- grDevices::colorRampPalette(palette)(n)
  }
  return(palette[1:n])
}


#' Visualize standard-stacked ensemble weights
#'
#' Creates a stacked bar plot showing the distribution of ensemble weights across
#' different learners and cross-validation folds for standard stacking.
#'
#' @param x An object of class `ens_weights_stand` containing cross-validated
#'   ensemble weights from standard stacking
#' @param ncols Number of columns for the facet grid (default: 3)
#' @param base_size Base font size for plot elements (default: 12)
#' @param ... Additional arguments passed to ggplot2 functions
#'
#' @return A ggplot object showing ensemble weights by fold and learner.
#'
#' @examples
#' \dontrun{
#' weights <- list(
#'   "Y.hat" = data.frame(
#'     fold1 = c(0.65, 0.10, 0.20, 0.05),
#'     fold2 = c(0.72, 0.08, 0.15, 0.05),
#'     row.names = c("ols", "forest_grf", "xgboost", "rlasso")
#'   ),
#'   "D.hat" = data.frame(
#'     fold1 = c(0.38, 0.10, 0.42, 0.10),
#'     fold2 = c(0.61, 0.10, 0.19, 0.10),
#'     row.names = c("ols", "forest_grf", "xgboost", "rlasso")
#'   )
#' )
#' class(weights) <- c("ens_weights_stand", "list")
#' plot(weights, ncols = 2)
#' }
#'
#' @export
#' @method plot ens_weights_stand
plot.ens_weights_stand <- function(x,
                                   ncols = 3,
                                   base_size = 12,
                                   ...) {
  # Input validation
  if (!inherits(x, "ens_weights_stand")) {
    stop("x must be an object of class 'ens_weights_stand'")
  }

  prepare_data <- function(model_name, model_data) {
    # Transpose model data to have folds as rows and methods as columns
    df <- as.data.frame(t(model_data))
    df$fold <- rownames(df)

    # Melt data from wide to long format for ggplot
    melted <- data.frame(
      fold = rep(df$fold, times = ncol(df) - 1),
      variable = rep(names(df)[-ncol(df)], each = nrow(df)),
      value = unlist(df[-ncol(df)]),
      stringsAsFactors = FALSE
    )

    melted$model <- model_name
    return(melted)
  }

  plot_data <- do.call(rbind, mapply(prepare_data, names(x), x, SIMPLIFY = FALSE))

  all_methods <- unique(plot_data$variable)
  method_palette <- stats::setNames(get_palette(length(all_methods)), all_methods)

  p <- ggplot2::ggplot(plot_data, ggplot2::aes(x = fold, y = value, fill = variable)) +
    ggplot2::geom_col(position = "stack", width = 0.7) +
    ggplot2::scale_fill_manual(values = method_palette, name = "Method") +
    # ggplot2::scale_y_continuous(limits = c(0, 1.01), expand = c(0, 0)) +
    ggplot2::facet_wrap(~model, ncol = ncols, scales = "free_x") +
    ggplot2::labs(x = NULL, y = "Ensemble Weight") +
    ggplot2::theme_minimal(base_size = base_size) +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(angle = 45, hjust = 0.75, size = base_size * 0.6),
      axis.text.y = ggplot2::element_text(size = base_size * 0.7),
      legend.position = "bottom",
      legend.text = ggplot2::element_text(size = base_size * 0.6),
      legend.title = ggplot2::element_text(size = base_size * 0.7)
    )

  return(p)
}


#' Visualize short-stacked ensemble weights
#'
#' Plots single set of ensemble weights for short stacking,
#' showing the weight assigned to each learner.
#'
#' @param x An object of class \code{ens_weights_short} containing ensemble weights
#' @param base_size Base font size for plot elements (default: 12)
#' @param ... Additional arguments passed to ggplot
#'
#' @return A ggplot object showing learner weights as a stacked bar plot
#'
#' @examples
#' \dontrun{
#' weights <- data.frame(
#'   "Y.hat" = c(0.48, 0.10, 0.42),
#'   "D.hat" = c(0.59, 0.10, 0.31),
#'   row.names = c("ols", "forest_grf", "xgboost")
#' )
#' class(weights) <- c("ens_weights_short", "data.frame")
#' plot(weights)
#' }
#'
#' @export
#' @method plot ens_weights_short
plot.ens_weights_short <- function(x,
                                   base_size = 12,
                                   ...) {
  # Input validation
  if (!inherits(x, "ens_weights_short")) {
    stop("x must be an object of class 'ens_weights_short'")
  }

  df <- as.data.frame(x)
  df$method <- rownames(df)

  df_long <- data.frame(
    method = rep(df$method, times = ncol(df) - 1),
    variable = rep(names(df)[-ncol(df)], each = nrow(df)),
    value = unlist(df[-ncol(df)]),
    stringsAsFactors = FALSE
  )

  all_methods <- unique(df$method)
  method_palette <- stats::setNames(get_palette(length(all_methods)), all_methods)

  ggplot2::ggplot(df_long, ggplot2::aes(x = variable, y = value, fill = method)) +
    ggplot2::geom_col(position = "stack", width = 0.7) +
    ggplot2::scale_fill_manual(values = method_palette, name = "Method") +
    ggplot2::labs(x = "Nuisance Parameter", y = "Weight", fill = "Method") +
    ggplot2::theme_minimal(base_size = base_size) +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))
}


#' Setup a Progress Bar from the \code{progress} package
#'
#' Creates a customized progress bar to track progress.
#'
#' @param NuPa Character vector specifying the nuisance parameters to estimate.
#'             Currently supported options:
#'             \code{c("Y.hat", "Y.hat.d", "Y.hat.z", "D.hat", "D.hat.z", "Z.hat")}
#' @param n_d Integer, number of unique treatment values.
#' @param n_z Integer, number of unique instrument values.
#' @param cf_folds Integer, number of cross-fitting folds.
#' @param cv_folds Integer, number of cross-validation folds.
#' @param methods List of methods to use for \code{\link{ensemble}} estimation.
#'
#' @return A progress bar object from the \code{progress} package, configured with:
#' \itemize{
#'   \item Custom format showing percentage, counts, nuisance parameter status
#'   \item Dynamic tokens for cross-fitting (pb_cf) and cross-validation (pb_cv) folds
#'   \item Task and model information
#' }
#'
#' @seealso \code{\link[progress]{progress_bar}} for the underlying progress bar implementation.
#'
#' @keywords internal
setup_pb <- function(NuPa, 
                     n_d, n_z, 
                     cf_folds, cv_folds, 
                     methods
                     ) {
  count_nested_lists <- function(methods, n_d, n_z) {
    multipliers <- list("Y.hat.d" = n_d, "Y.hat.z" = n_z, "D.hat.z" = n_z)

    # Calculate weighted counts
    sum(sapply(names(methods), function(key) {
      base_count <- length(methods[[key]])
      multiplier <- ifelse(key %in% names(multipliers), multipliers[[key]], 1)
      base_count * multiplier
    }))
  }

  # Final ticks: cf folds × methods × 2 (i.e. fitting+prediction) × cv folds
  total_ticks <- cf_folds * count_nested_lists(methods, n_d, n_z) * 2 * if (cv_folds > 1) (cv_folds + 1) else 1 # && length(methods) > 1

  pb <- progress::progress_bar$new(
    format = "[:bar] :percent | :current/:total | :nuisance | cf =:pb_cf, cv =:pb_cv | :task :model",
    total = total_ticks, clear = TRUE, width = 80, force = FALSE
  )

  return(pb)
}


#' Method creation for ensemble
#'
#' Creates methods to be used in the subsequent \code{\link{ensemble}} model.
#'
#' @param method Choose method from available options. See details for supported methods.
#' @param x_select Optional logical vector (length = number of covariates) indicating which
#'  variables to use. For example, tree-based methods typically exclude the
#'  interactions used by Lasso.
#' @param arguments Optional list of additional arguments passed to the underlying method.
#' @param tuning Hyperparameter tuning options for GRF's regression forest and XGBoost:
#'   \describe{
#'     \item{\code{"full_sample"}}{Tune hyperparameters using full sample.}
#'     \item{\code{"fold"}}{Tuning is performed on the estimation part of the cross-fitting split.
#'                          Roughly \eqn{F} times more computationally intensive.}
#'     \item{\code{"no"}}{No tuning performed; uses default settings.}
#'   }
#' @param multinomial Optional character specifying multiclass handling approach:
#'   one of \code{c("one-vs-one", "one-vs-rest", "multiclass")}.
#' @param parallel Optional logical. If \code{TRUE}, enables parallelization for
#'   the One-Vs-One (OvO) multiclass routine. Requires \code{foreach} and
#'   \code{doParallel}.
#'
#' @return A list object that can be passed as input to \code{\link{ensemble}}.
#'
#' @details
#' Supported methods include:
#'
#' \itemize{
#'   \item \code{"mean"}: Mean difference.
#'
#'   \item \code{"ols"}: Ordinary Least Squares.
#'
#'   \item \code{"ridge"}: Ridge regression via \code{glmnet}. Uses 10-fold
#'   cross-validation (over 100 log-spaced values by default) to select the
#'   regularization strength \eqn{\lambda}. Users may supply a custom
#'   \eqn{\lambda} sequence, choose the CV loss, and set the number of folds.
#'
#'   \item \code{"plasso"}: Post-Lasso estimator via \code{plasso} built on top of the 
#'   \code{glmnet} package, with tuning settings identical to Ridge.
#'
#'   \item \code{"lasso"}: Standard Lasso via \code{glmnet}, with tuning settings
#'   identical to Ridge.
#'
#'   \item \code{"rlasso"}: Lasso via \code{hdm} with a theory-driven,
#'   data-dependent penalty robust to heteroskedastic and non-Gaussian errors.
#'   By default, \code{rlasso()} includes sets the penalty with theoretical choices 
#'   (\eqn{c = 1.1}, \eqn{\gamma = 0.1 / \log(n)}).
#'
#'   \item \code{"forest_grf"}: Regression forest via \code{grf}, with defaults
#'   \code{num.trees = 2000}, \code{min.node.size = 5},
#'   \code{sample.fraction = 0.5}, and \code{honesty = TRUE}.
#'
#'   \itemize{
#'     \item If \code{tuning = "full_sample"}, tuning is performed on the
#'     full sample \eqn{(X,Y)} over \code{sample.fraction}, \code{mtry},
#'     \code{min.node.size}, \code{honesty.fraction}, \code{honesty.prune.leaves},
#'     \code{alpha}, and \code{imbalance.penalty}.
#'
#'     \item If \code{tuning = "fold"}, tuning is performed on the
#'     estimation part of the cross-fitting split (\eqn{F–1} folds), which is roughly
#'     \eqn{F} times more computationally demanding.
#'   }
#'
#'   \item \code{"xgboost"}: Gradient boosting via \code{xgboost}, using 100 
#'   boosting rounds by default. These hyperparameters are fixed to unsure smoother 
#'   extraction: \code{reg_alpha = 0}, \code{subsample = 1}, \code{max_delta_step = 0}, 
#'   \code{base_score = 0}. Supports the same tuning logic as regression
#'   forests. The Hyperband-like tuning routine and tunable hyperparameters are
#'   described in \code{?tune_xgb_hyperband}.
#'
#'   \item \code{"knn"}: k-Nearest Neighbors via \code{FastKNN}, with
#'   \eqn{k = 10} neighbors by default.
#'
#'   \item \code{"forest_drf"}: Distributional random forest via \code{drf}, with
#'   defaults \code{min.node.size = 15}, \code{num.trees = 2000},
#'   \code{splitting.rule = "FourierMMD"}.
#'   
#'   \item \code{"glm"}: Binary regression models using \code{glm}. Supports logistic
#'   (logit) and probit link functions. To use probit, specify 
#'   \code{family = binomial(link = "probit")}. Default is logit.
#'
#'   \item \code{"logit"}: Logistic regression via \code{glmnet}. Uses
#'   \code{family = "binomial"} for binary outcomes and
#'   \code{family = "multinomial"} for multiclass outcomes.
#'
#'   \item \code{"logit_nnet"}: Logistic regression via \code{nnet::multinom()},
#'   supports binary and multiclass outcomes.
#'
#'   \item \code{"nb_gaussian"}: Gaussian Naive Bayes via
#'   \code{naivebayes::gaussian_naive_bayes()}.
#'
#'   \item \code{"nb_bernoulli"}: Bernoulli Naive Bayes via
#'   \code{naivebayes::naive_bayes()} with Bernoulli likelihood.
#'
#'   \item \code{"xgboost_prop"}: Gradient boosting via \code{xgboost} for
#'   classification or propensity score estimation. Uses
#'   \code{objective = "binary:logistic"} for binary outcomes and
#'   \code{objective = "multi:softprob"} for multiclass outcomes. Defaults to
#'   100 boosting rounds unless \code{nrounds} is specified.
#'
#'   \item \code{"svm"}: Support Vector Machine via \code{e1071::svm()}.
#'   Defaults to \code{kernel = "radial"}, \code{type = "C-classification"},
#'   and enables probability estimates (\code{probability = TRUE}).
#'   
#'   \item \code{"prob_forest"}: Probability forest via
#'   \code{grf::probability_forest()}. Grows 2000 trees by default.
#'
#'   \item \code{"ranger"}: Random forest classifier via \code{ranger}.
#'   Grows 500 trees by default.
#'
#'   \item \code{"knn_prop"}: k-Nearest Neighbors classifier via
#'   \code{kknn::train.kknn()}.
#' }
#'
#' @export
#'
#' @examples
#' # Create list of methods for ensemble
#' methods <- list(
#'   "ols" = create_method("ols"),
#'   "knn" = create_method("knn", arguments = list("k" = 3)),
#'   "forest_grf" = create_method("forest_grf"),
#'   "logit_ovo" = create_method("logit", multinomial = "one-vs-one", parallel = TRUE),
#'   "prob_forest" = create_method("prob_forest", multinomial = "multiclass")
#' )
create_method <- function(method = c("mean", "ols", "ridge", "plasso", "forest_grf", 
                                     "lasso", "knn", "forest_drf", "xgboost", "rlasso", 
                                     "glm", "logit", "logit_nnet", "nb_gaussian", "nb_bernoulli",
                                     "xgboost_prop", "svm", "prob_forest", "ranger", "knn_prop"
                                     ),
                          x_select = NULL,
                          arguments = list(),
                          tuning = "no",
                          multinomial = NULL,
                          parallel = FALSE) {
  # Sanity checks
  method <- match.arg(method)
  tuning <- match.arg(tuning, choices = c("no", "full_sample", "fold"))

  if (!(is.null(arguments) || is.list(arguments))) {
    stop("Provide either NULL or a list for arguments.")
  }

  if (!(is.null(x_select) || is.logical(x_select))) {
    stop("Provide either NULL or logical for x_select.")
  }

  if (!is.logical(parallel) || length(parallel) != 1) {
    stop("parallel must be a single logical value (TRUE or FALSE)")
  }

  if (!is.null(multinomial)) {
    multinomial <- match.arg(multinomial, choices = c("one-vs-one", "one-vs-rest", "multiclass"))
    if (multinomial == "multiclass" && method %in% c("svm")) {
      stop(
        "SVM does not support multiclass estimation.",
        "Please set multinomial to either 'one-vs-one' or 'one-vs-rest'."
      )
    }
  }

  return(list(
    method = method,
    x_select = x_select,
    arguments = arguments,
    tuning = tuning,
    multinomial = multinomial,
    parallel = parallel
  ))
}


#' Method creation for ensemble (internal version)
#'
#' Creates the methods to be actually used in the subsequent \code{\link{ensemble}} model,
#' allowing for either a simple method list or nuisance-parameter-specific method lists.
#'
#' @param methods Either:
#' \enumerate{
#'   \item{A list of methods (like in \code{create_method}) that will be used for all nuisance parameters, or}
#'   \item{A named list where each name corresponds to a nuisance parameter and contains a list of methods for that parameter}
#' }
#' @param NuPa Character vector of nuisance parameters to be estimated
#' @param K Number of unique treatment statuses (default is 2, binary treatment)
#'
#' @return Nested list object organized by nuisance parameter that will be passed as input to \code{\link{ensemble}}
#'
#' @keywords internal
process_methods <- function(methods, NuPa, K) {
  available_NuPa <- c("Y.hat", "Y.hat.d", "Y.hat.z", "D.hat", "D.hat.z", "Z.hat")
  multiclass_method <- c( "glm", "logit", "logit_nnet", "nb_gaussian", "nb_bernoulli", "xgboost_prop", "svm", "prob_forest", "ranger", "knn_prop")
  base_method <- c("mean", "ols", "ridge", "plasso", "forest_grf", "lasso", "knn", "forest_drf", "xgboost", "rlasso")

  # A function to extract method types from a methods list
  get_method_types <- function(method_list) {
    sapply(method_list, function(x) x$method)
  }

  # Check if a "one for all" methods list
  is_simple_list <- !any(names(methods) %in% available_NuPa)

  if (is_simple_list) {
    # Case 1: "One for all" methods list - replicate for all NuPa after filtering
    result <- stats::setNames(vector("list", length(NuPa)), NuPa)

    filter_multiclass <- NuPa[(NuPa %in% c("Y.hat", "Y.hat.d", "Y.hat.z"))]
    filter_base <- NuPa[NuPa %in% c("D.hat", "D.hat.z") & K > 2]

    # If any methods need to be removed
    method_types <- get_method_types(methods)
    multiclass_in_method <- any(method_types %in% multiclass_method)
    base_in_method <- any(method_types %in% base_method)

    if (length(filter_multiclass) > 0 && multiclass_in_method) {
      message(
        "Multiclass methods incompatible with NuPa: ",
        paste(filter_multiclass, collapse = ", "), ", they will be omitted."
      )
    }
    if (length(filter_base) > 0 && base_in_method) {
      message(
        "Base methods incompatible with NuPa: ",
        paste(filter_base, collapse = ", "), ", they will be omitted."
      )
    }

    # Apply filters
    for (nupa in NuPa) {
      if (nupa %in% filter_multiclass) {
        result[[nupa]] <- methods[!method_types %in% multiclass_method]
      } else if (nupa %in% filter_base) {
        result[[nupa]] <- methods[!method_types %in% base_method]
      } else {
        result[[nupa]] <- methods
      }
    }
  } else {
    # Case 2: NuPa-specific methods list
    result <- stats::setNames(vector("list", length(NuPa)), NuPa)
    multiclass_nupas <- base_nupas <- character(0)

    for (nupa in intersect(names(methods), NuPa)) {
      method_types <- get_method_types(methods[[nupa]])

      # Multiclass methods to remove
      if ((nupa %in% c("Y.hat", "Y.hat.d", "Y.hat.z"))) {
        if (any(method_types %in% multiclass_method)) {
          multiclass_nupas <- c(multiclass_nupas, nupa)
        }
      }
      # Base methods to remove
      if (nupa %in% c("D.hat", "D.hat.z") && K > 2) {
        if (any(method_types %in% base_method)) {
          base_nupas <- c(base_method, nupa)
        }
      }
    }

    if (length(multiclass_nupas) > 0) {
      message(
        "Multiclass methods incompatible with NuPa: ",
        paste(unique(multiclass_nupas), collapse = ", "), ", they will be omitted."
      )
    }
    if (length(base_nupas) > 0) {
      message(
        "Base methods incompatible with NuPa: ",
        paste(unique(base_nupas), collapse = ", "), ", they will be omitted."
      )
    }

    # Apply filters
    for (nupa in intersect(names(methods), NuPa)) {
      method_types <- get_method_types(methods[[nupa]])

      if ((nupa %in% c("Y.hat", "Y.hat.d", "Y.hat.z"))) {
        result[[nupa]] <- methods[[nupa]][!method_types %in% multiclass_method]
      } else if (nupa %in% c("D.hat", "D.hat.z") && K > 2) {
        result[[nupa]] <- methods[[nupa]][!method_types %in% base_method]
      } else {
        result[[nupa]] <- methods[[nupa]]
      }
    }
  }

  empty_nupas <- names(which(lengths(result) == 0))
  if (length(empty_nupas) > 0) {
    stop(paste(
      "The following nuisance parameters have empty methods lists:",
      paste(empty_nupas, collapse = ", "), "\nPlease provide at least one valid method for each NuPa"
    ))
  }

  return(result)
}


#' Add or Replace Nuisance Parameters in a `NuisanceParameters` Object
#'
#' This function allows to add or replace nuisance parameters and their corresponding models
#' from one `NuisanceParameters` object to another.
#'
#' @param np A `NuisanceParameters` object (the original object to be modified)
#' @param np_new A `NuisanceParameters` object (the donor object)
#' @param NuPa Character vector specifying which nuisance parameters to add/replace.
#'             If NULL (default), all available parameters from `np_new` will be used.
#' @param replace Logical indicating whether to replace existing parameters (TRUE)
#'                or only fill empty slots (FALSE, default)
#'
#' @return A modified `NuisanceParameters` object with updated nuisance parameters and models
#'
#' @export
#'
#' @examples
#' \dontrun{
#' set.seed(123)
#'
#' N <- 100
#' Y <- rnorm(N)
#' X <- matrix(rnorm(N * 10), ncol = 10)
#' D <- Z <- rbinom(N, 1, 0.5)
#'
#' methods1 <- list(
#'   "ols" = create_method("ols"),
#'   "plasso" = create_method("plasso")
#' )
#'
#' methods2 <- list(
#'   "forest_grf" = create_method("forest_grf"),
#'   "xgboost" = create_method("xgboost")
#' )
#'
#' np1 <- nuisance_parameters(
#'   NuPa = c("Y.hat", "D.hat"), X = X, Y = Y, D = D,
#'   methods = methods1, cf = 5, stacking = "short", storeModels = "Memory"
#' )
#'
#' np2 <- nuisance_parameters(
#'   NuPa = c("D.hat", "Z.hat"), X = X, D = D, Z = Z,
#'   cf_mat = np1$numbers$cf_mat, methods = methods2, cf = 5,
#'   stacking = "short", storeModels = "Memory"
#' )
#'
#' np3 <- add_nupa(np = np1, np_new = np2, NuPa = c("D.hat", "Z.hat"), replace = TRUE)
#' }
add_nupa <- function(np, np_new, NuPa = NULL, replace = FALSE) {
  # Check if both objects are of class NuisanceParameters
  if (!inherits(np, "NuisanceParameters") || !inherits(np_new, "NuisanceParameters")) {
    stop("Both np and np_new must be of class 'NuisanceParameters' (created by nuisance_parameters function)")
  }

  # Check if cross-fitting splits (cf_mat) are identical
  if (!identical(np$numbers$cf_mat, np_new$numbers$cf_mat)) {
    stop(
      "Cross-fitting splits are not identical. ",
      "Make sure to use the same cf_mat in both NuisanceParameters objects"
    )
  }

  # If NuPa not specified, use all available parameters from np_new
  if (is.null(NuPa)) {
    NuPa <- names(np_new$nuisance_parameters)[sapply(
      np_new$nuisance_parameters,
      function(x) {
        !is.character(x) ||
          !grepl("not specified", x)
      }
    )]
  }

  # Verify specified NuPa exist in np_new
  invalid_params <- setdiff(NuPa, names(np_new$nuisance_parameters))
  if (length(invalid_params) > 0) {
    stop(
      "The following parameters are not present in np_new: ",
      paste(invalid_params, collapse = ", ")
    )
  }

  ## Core
  # Process each specified parameter
  for (param in NuPa) {
    # Check if parameter is actually estimated in np_new
    if (is.character(np_new$nuisance_parameters[[param]]) &&
      grepl("not specified", np_new$nuisance_parameters[[param]])) {
      warning("Parameter '", param, "' is not specified in np_new and won't be added")
      next
    }

    # Check if we should replace or only fill empty slots
    if (replace ||
      (is.character(np$nuisance_parameters[[param]]) &&
        grepl("not specified", np$nuisance_parameters[[param]]))) {
      # Add the nuisance parameter and the methods
      np$nuisance_parameters[[param]] <- np_new$nuisance_parameters[[param]]
      np$numbers$methods[[param]] <- np_new$numbers$methods[[param]]

      if (np$numbers$cv == 1) {
        # Combine the ensemble weights
        w_ens <- as.matrix(np[["numbers"]][["ens_weights"]])
        w_new_ens <- as.matrix(np_new[["numbers"]][["ens_weights"]])

        all_methods <- union(rownames(w_ens), rownames(w_new_ens))
        all_nupas <- union(colnames(w_ens), colnames(w_new_ens))

        ens_weights <- matrix(0,
          nrow = length(all_methods), ncol = length(all_nupas),
          dimnames = list(all_methods, all_nupas)
        )

        ens_weights[rownames(w_ens), colnames(w_ens)] <- w_ens
        ens_weights[, NuPa] <- 0
        ens_weights[rownames(w_new_ens), NuPa] <- w_new_ens[, NuPa, drop = FALSE]

        ens_weights <- as.data.frame(ens_weights)
        class(ens_weights) <- c("ens_weights_short", "data.frame")
        np[["numbers"]][["ens_weights"]] <- ens_weights
      } else {
        w_ens <- np[["numbers"]][["ens_weights"]]
        w_new_ens <- np_new[["numbers"]][["ens_weights"]]
        ens_weights <- w_ens

        # Completely replace the NuPa elements with those from np_new
        for (nupa in NuPa) {
          if (nupa %in% names(w_new_ens)) {
            ens_weights[[nupa]] <- w_new_ens[[nupa]]
          }
        }

        # Add any new NuPa elements from np_new that weren't in the original
        new_nupas <- setdiff(names(w_new_ens), names(w_ens))
        new_nupas_to_add <- intersect(new_nupas, NuPa)

        for (nupa in new_nupas_to_add) {
          ens_weights[[nupa]] <- w_new_ens[[nupa]]
        }

        class(ens_weights) <- c("ens_weights_stand", "list")
        np[["numbers"]][["ens_weights"]] <- ens_weights
      }

      # Add the corresponding model if it exists
      model_name <- paste0(param, "_m")
      if (model_name %in% names(np_new$models)) {
        np$models[[model_name]] <- np_new$models[[model_name]]
      }
    }
  }

  return(np)
}


#' One-Hot Encoding (numeric output)
#'
#' Convert categorical vector to one-hot encoded matrix.
#'
#' @param Y A factor or character vector.
#'
#' @return A numeric matrix with one column per unique category.
#'
#' @keywords internal
#'
one_hot <- function(Y) {
  uY <- unique(Y)
  out <- matrix(0L, nrow = length(Y), ncol = length(uY))
  out[cbind(seq_along(Y), match(Y, uY))] <- 1L
  colnames(out) <- uY
  
  return(out)
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
#' @param seed Random seed for reproducibility. Default is 123.
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
tune_xgb_hyperband <- function(X, Y,
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

  for (budget in budgets) {
    scores <- numeric(length(configs))

    for (i in seq_along(configs)) {
      params <- configs[[i]]
      cv <- tryCatch(
        xgboost::xgb.cv(
          data = as.matrix(X),
          label = Y,
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
  
  # Important for standard stacking 
  best_params$base_score <- NULL
  
  return(list("params" = best_params, "best_score" = best_score))
}


#' Tune learner hyperparameters
#'
#' Tunes hyperparameters for specified learners in a \code{methods} list
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
  
  # Process each method in the list
  for (i in seq_along(methods)) {
    mtd <- methods[[i]]
    
    if (mtd$tuning == "no" || !identical(type, mtd$tuning)) {
      next
    }
    
    if (is.null(mtd$x_select)) {
      X_sub <- X
    } else {
      X_sub <- X[, mtd$x_select, drop = FALSE]
    }
    
    # grf's random forest tuning
    if (identical(mtd$method, "forest_grf")) {
      tuned_forest <- grf::regression_forest(
        X = X_sub,
        Y = Y,
        tune.parameters = "all",
        num.trees = 50
      )
      
      # Merge with existing arguments
      tuned_params <- tuned_forest[["tunable.params"]]
      methods[[i]][["arguments"]] <- utils::modifyList(methods[[i]][["arguments"]], tuned_params)
    }
    
    # XGBoost tuning using hyperband
    else if (identical(mtd$method, "xgboost")) {
      tuned_xgb <- tune_xgb_hyperband(X = X_sub, Y = Y)
      methods[[i]][["arguments"]] <- tuned_xgb[["params"]]
    }
  }
  
  return(methods)
}
