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
#' @keywords internal
#'
#' @examples
#' D <- factor(c("A", "B", "A", "C", "B", "C"))
#' d_mat <- prep_indicator_mat(D)
#' head(d_mat)
#'
prep_indicator_mat <- function(x) {
  # Input validation
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
#' @description
#' \code{\link{prep_cf_mat}} creates a matrix of binary cross-fitting fold
#' indicators (N x # cross-folds)
#'
#' @param N Number of observations
#' @param cf Number of cross-fitting folds
#' @param cluster Optional vector of cluster variable if cross-fitting should account
#'                for clusters within the data.
#' @param d_mat Optional logical matrix of treatment indicators (N x T+1 with T
#'              being the number of treatments). For example created by \code{\link{prep_indicator_mat}}.
#'              If specified, cross-fitting folds will preserve the treatment ratios from full sample.
#'              However, if cluster vector is provided, \code{d_mat} is ignored due to computational
#'              constraints and randomization as well as feasibility issues.
#'
#' @return Logical matrix of cross-fitting folds (N x # folds).
#'
#' @keywords internal
#'
#' @examples
#' N <- 1000
#' d_mat <- prep_indicator_mat(sample(3, N, replace = TRUE))
#' cf_mat <- prep_cf_mat(N, cf = 5, d_mat = d_mat)
#' head(cf_mat)
#'
prep_cf_mat <- function(N, 
                        cf, 
                        cluster = NULL, 
                        d_mat = NULL) {
  # Check if both cluster vector and treatment matrix are provided
  if (!is.null(cluster) & !is.null(d_mat)) {
    warning("You provided both a cluster vector and a treatment matrix. 
            This is not feasible due to computational constraints and randomization issues. 
            Thus, only cluster vector is considered.")
    d_mat <- NULL
  }

  # Only one fold (i.e. no cross-fitting)
  if (cf == 1) {
    cf_mat <- matrix(rep(1, N), ncol = 1)

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

    if (any(cf_mat_balance < (1 / cf) * 0.5)) 
      stop("There is a high imbalance in the cluster sizes. This poses a problem for the cross-fitting procedure. 
           Either choose a smaller number of folds or do not specify a cluster vector.")
  }

  colnames(cf_mat) <- sprintf("CF %d", 1:cf)

  return(cf_mat)
}


#' Non-negative least squares function for estimation of ensemble weights
#'
#' @description
#' \code{\link{nnls_weights}} estimates ensemble weights for a prediction problem
#' on the basis of the non-negative least squares algorithm that puts a positive
#' constraint on the regression coefficients.
#'
#' @param X A matrix where each column represents a different predictor variable
#' and each row represents an observation
#' @param Y A numeric vector of actual target values
#'
#' @return A numeric vector of the non-negative least weights.
#'
#' @keywords internal
#'
#' @examples
#' X <- matrix(c(0.2, 0.3, 0.4, 0.5, 0.6, 0.7), ncol = 2)
#' Y <- c(0.25, 0.45, 0.65)
#' \donttest{
#' nnls_w <- nnls_weights(X, Y)
#' }
#'
nnls_weights <- function(X, Y) {
  
  nnls_result <- nnls::nnls(as.matrix(X), Y)
  nnls_w <- nnls_result$x

  # In case of perfectly agreeing predictions, nnls provides only zeros.
  # In this scenario: uniform weights that add up to 1
  if (sum(nnls_w) == 0) {
    nnls_w <- nnls_w + 1 / length(nnls_w)
  }

  # Normalize weights
  nnls_w <- nnls_w / sum(nnls_w)

  if (!is.null(colnames(X))) names(nnls_w) <- colnames(X)

  return(nnls_w)
}


#' Adds an intercept to a matrix
#'
#' @description
#' \code{\link{add_intercept}} adds an intercept to a matrix.
#'
#' @param mat Any matrix.
#'
#' @return Matrix with intercept.
#'
#' @keywords internal
#'
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


#' Aggregate 3D array
#'
#' @description
#' This function computes the weighted sum along the third dimension of a
#' three dimensional array.
#'
#' @param a A 3D array
#' @param w A numeric vector.
#' The length of w should match the third dimension of a.
#'
#' @return A matrix of the same dimension as the first two dimensions of a.
#'
#' @keywords internal
#'
agg_array <- function(a, w) {
  return(apply(a, c(1, 2), function(x) sum(x * w)))
}


#' Prepare fit_cv matrix for cross-fitted ensemble predictions
#'
#' @description
#' Creates a matrix or array to store cross-fitted ensemble predictions.
#'
#' @param method List of methods built via \code{\link{create_method}} to be used in ensemble
#' @param N Number of observations
#' @param Y Outcome vector (used to determine number of classes for multinomial outcomes)
#'
#' @return A matrix of NAs with dimensions \code{N} x \code{length(method)} for binary outcomes,
#'         or a 3D array with dimensions \code{N} x \code{K} x \code{length(method)} for multinomial
#'         outcomes (where K is number of unique classes in Y)
#'
#' @keywords internal
#'
make_fit_cv <- function(method, N, Y) {
  is_mult <- !is.null(method[[1]]$multinomial)
  n_methods <- length(method)
  method_names <- names(method)
  
  if (is_mult) {
    n_classes <- length(unique(Y))
    fit_cv <- array(NA, dim = c(N, n_classes, n_methods))
    class_names <- if (!is.null(levels(Y))) levels(Y) else seq_len(n_classes) - 1
    dimnames(fit_cv) <- list(NULL, class_names, method_names)
  } else {
    fit_cv <- matrix(NA, nrow = N, ncol = n_methods)
    colnames(fit_cv) <- if (!is.null(method_names)) method_names else sprintf("method%s", seq_len(n_methods))
  }
  
  return(fit_cv)
}


#' Short- or Standard-stacking message
#'
#' @description
#' Prints a message into the console whether short- or standard-stacking is used.
#'
#' @param cv Cross-validation value for ensemble weights.
#'
#' @return Message in console. No object is returned.
#'
#' @keywords internal
#'
which_stacking <- function(cv = 1) {
  if (cv == 1) {
    message("Short-stacking is used.")
  } else if (cv > 1) {
    message("Standard-stacking is used.")
  }
}


#' Update progress bar
#'
#' @description
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
#' 
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


#' Format ensemble weights for short- or standard-stacking
#'
#' @param ens_weights List of ensemble weights to be formatted
#' @param cv Integer indicating number of cross-validation folds
#' 
#' @return A data.frame (\code{cv = 1}) or list of data.frames (\code{cv > 1}) with formatted weights
#' 
#' @keywords internal
#' 
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
#' @description
#' Creates a consistent color palette for visualizing ensemble method weights.
#' Returns a predefined set of colors for up to 12 methods, and generates
#' additional colors using \code{colorRampPalette} when more methods are needed.
#'
#' @param n Number of colors to return (number of distinct methods)
#'
#' @return A character vector of hex color codes of length `n`
#'
#' @keywords internal
#'
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
#' @description
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
#'
#' @method plot ens_weights_stand
#' 
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
    ggplot2::scale_y_continuous(limits = c(0, 1.01), expand = c(0, 0)) +
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
#' @description
#' Plots single set of ensemble weights for short stacking,
#' showing the weight assigned to each learner.
#'
#' @param x An object of class \code{ens_weights_short} containing ensemble weights
#' @param base_size Base font size for plot elements (default: 12)
#' @param ... Additional arguments passed to ggplot
#'
#' @return A ggplot object showing learner weights as a stacked bar plot
#'
#' @export
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
#' @method plot ens_weights_short
#' 
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
#' @description
#' Creates a customized progress bar to track progress.
#'
#' @param NuPa Character vector specifying the nuisance parameters to estimate.
#'             Currently supported options:
#'             \code{c("Y.hat", "Y.hat.d", "Y.hat.z", "D.hat", "D.hat.z", "Z.hat")}
#' @param n_d Integer, number of unique treatment values.
#' @param n_z Integer, number of unique instrument values.
#' @param cf_folds Integer, number of cross-fitting folds.
#' @param cv_folds Integer, number of cross-validation folds.
#' @param models List of methods to use for \code{\link{ensemble}} estimation.
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
#' 
setup_pb <- function(NuPa, n_d, n_z, cf_folds, cv_folds, models) {
  count_nested_lists <- function(models, n_d, n_z) {
    multipliers <- list("Y.hat.d" = n_d, "Y.hat.z" = n_z, "D.hat.z" = n_z)

    # Calculate weighted counts
    sum(sapply(names(models), function(key) {
      base_count <- length(models[[key]])
      multiplier <- ifelse(key %in% names(multipliers), multipliers[[key]], 1)
      base_count * multiplier
    }))
  }

  # Final ticks: cf folds × models × 2 (i.e. fitting+prediction) × cv folds
  total_ticks <- cf_folds * count_nested_lists(models, n_d, n_z) * 2 * if (cv_folds > 1) (cv_folds + 1) else 1 # && length(models) > 1

  pb <- progress::progress_bar$new(
    format = "[:bar] :percent | :current/:total | :nuisance | cf =:pb_cf, cv =:pb_cv | :task :model",
    total = total_ticks, clear = TRUE, width = 80, force = FALSE
  )
  
  return(pb)
}


#' Method creation for ensemble
#'
#' @description
#' Creates the methods to be used in the subsequent \code{\link{ensemble}} model.
#'
#' @param method Choose method from available options. See details for supported methods.
#' @param x_select Optional logical vector of length equal to the number of
#'   columns of the covariate matrix indicating which variables should be used by
#'   this method. E.g., tree-based methods usually should not be provided with
#'   interactions that Lasso is using.
#' @param arguments Optional list containing additional arguments that should be
#'   passed to the underlying method, or a character string specifying tuning
#'   mode. See details.
#' @param multinomial Optional character specifying multiclass handling approach.
#'   One of \code{c("one-vs-one", "one-vs-rest", "multiclass")}.
#' @param parallel Optional logical. If TRUE, attempts parallelization for the One-Vs-One (OvO) 
#'   multiclass routine. Requires \code{foreach} and \code{doParallel} packages to be installed.
#' @param name Optional string naming the method.
#'
#' @return List object that can be passed as input to \code{\link{ensemble}}
#'
#' @details
#' Supported methods: \code{"mean"}, \code{"ols"}, \code{"ridge"}, \code{"plasso"},
#' \code{"forest_grf"}, \code{"lasso"}, \code{"knn"}, \code{"forest_drf"},
#' \code{"xgboost"}, \code{"rlasso"}, \code{"logit"}, \code{"logit_nnet"},
#' \code{"nb_gaussian"}, \code{"nb_bernoulli"}, \code{"xgboost_prop"}, \code{"svm"},
#' \code{"prob_forest"}, \code{"ranger"}, \code{"knn_prop"}.
#'
#' For the \code{arguments} parameter, if a character string is provided, it must
#' be one of \code{c("tune_on_fold", "tune_full_sample")} to specify tuning mode.
#'
#' @export
#'
#' @examples
#' # Create list of methods for ensemble
#' method <- list(
#'   "ols" = create_method("ols"),
#'   "knn" = create_method("knn", arguments = list("k" = 3)),
#'   "forest_grf" = create_method("forest_grf", arguments = c("tune_full_sample")),
#'   "logit_ovo" = create_method("logit", multinomial = "one-vs-one", parallel = TRUE),
#'   "prob_forest" = create_method("prob_forest", multinomial = "multiclass")
#' )
#'
create_method <- function(
    method = c(
      "mean", "ols", "ridge", "plasso", "forest_grf", "lasso", "knn", "forest_drf", 
      "xgboost", "rlasso", "logit", "logit_nnet", "nb_gaussian", "nb_bernoulli", 
      "xgboost_prop", "svm", "prob_forest", "ranger", "knn_prop"
    ),
    x_select = NULL,
    multinomial = NULL,
    parallel = FALSE,
    arguments = list(),
    name = NULL) {
  
  # Sanity checks
  method <- match.arg(method)
  
  # Check if other inputs are valid
  if (!(is.null(arguments) || is.list(arguments) || is.character(arguments))) {
    stop("Provide either NULL, a list, or a character string for arguments.")
  }
  
  if (is.character(arguments) && length(arguments) == 1) {
    if (!arguments %in% c("tune_on_fold", "tune_full_sample")) {
      stop("If arguments is a character string, it must be one of: 'tune_on_fold', 'tune_full_sample'")
    }
  }
  
  if (!(is.null(x_select) || is.logical(x_select))) {
    stop("Provide either NULL or logical for x_select.")
  }
  
  if (!((is.character(name) && length(name) == 1) || is.null(name))) {
    stop("Provide single string to name method.")
  }
  
  if (!is.logical(parallel) || length(parallel) != 1) {
    stop("parallel must be a single logical value (TRUE or FALSE)")
  }
  
  if (!is.null(multinomial)) {
    multinomial <- match.arg(multinomial, choices = c("one-vs-one", "one-vs-rest", "multiclass"))
    if (multinomial == "multiclass" && method %in% c("svm")) {
      stop("SVM does not support multiclass estimation.", 
           "Please set multinomial to either 'one-vs-one' or 'one-vs-rest'.")
    }
  }
  
  # Return method specification
  return(list(
    method = method, 
    multinomial = multinomial, 
    parallel = parallel,
    arguments = arguments, 
    x_select = x_select, 
    name = name
  ))
}


#' Method creation for ensemble (internal version)
#'
#' @description
#' Creates the methods to be actually used in the subsequent \code{\link{ensemble}} model,
#' allowing for either a simple method list or nuisance-parameter-specific method lists.
#'
#' @param method Either:
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
create_final_method <- function(method, NuPa, K) {
  available_NuPa <- c("Y.hat", "Y.hat.d", "Y.hat.z", "D.hat", "D.hat.z", "Z.hat")
  multiclass_method <- c("logit", "logit_nnet", "nb_gaussian", "nb_bernoulli", "xgboost_prop", "svm", "prob_forest", "ranger", "knn_prop")
  base_method <- c("mean", "ols", "ridge", "plasso", "forest_grf", "lasso", "knn", "forest_drf", "xgboost", "rlasso")

  # A function to extract method types from a method list
  get_method_types <- function(method_list) {
    sapply(method_list, function(x) x$method)
  }

  # Check if a "one for all" method list
  is_simple_list <- !any(names(method) %in% available_NuPa)

  if (is_simple_list) {
    # Case 1: "One for all" method list - replicate for all NuPa after filtering
    result <- stats::setNames(vector("list", length(NuPa)), NuPa)

    filter_multiclass <- NuPa[(NuPa %in% c("Y.hat", "Y.hat.d", "Y.hat.z"))]
    filter_base <- NuPa[NuPa %in% c("D.hat", "D.hat.z") & K > 2]

    # If any methods need to be removed
    method_types <- get_method_types(method)
    multiclass_in_method <- any(method_types %in% multiclass_method)
    base_in_method <- any(method_types %in% base_method)

    if (length(filter_multiclass) > 0 && multiclass_in_method) {
      message("Multiclass methods incompatible with NuPa: ", 
              paste(filter_multiclass, collapse = ", "), ", they will be omitted.")
    }
    if (length(filter_base) > 0 && base_in_method) {
      message("Base methods incompatible with NuPa: ", 
              paste(filter_base, collapse = ", "), ", they will be omitted.")
    }

    # Apply filters
    for (nupa in NuPa) {
      if (nupa %in% filter_multiclass) {
        result[[nupa]] <- method[!method_types %in% multiclass_method]
      } else if (nupa %in% filter_base) {
        result[[nupa]] <- method[!method_types %in% base_method]
      } else {
        result[[nupa]] <- method
      }
    }
  } else {
    # Case 2: NuPa-specific method list
    result <- stats::setNames(vector("list", length(NuPa)), NuPa)
    multiclass_nupas <- base_nupas <- character(0)

    for (nupa in intersect(names(method), NuPa)) {
      method_types <- get_method_types(method[[nupa]])

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
      message("Multiclass methods incompatible with NuPa: ", 
              paste(unique(multiclass_nupas), collapse = ", "), ", they will be omitted.")
    }
    if (length(base_nupas) > 0) {
      message("Base methods incompatible with NuPa: ", 
              paste(unique(base_nupas), collapse = ", "), ", they will be omitted.")
    }

    # Apply filters
    for (nupa in intersect(names(method), NuPa)) {
      method_types <- get_method_types(method[[nupa]])

      if ((nupa %in% c("Y.hat", "Y.hat.d", "Y.hat.z"))) {
        result[[nupa]] <- method[[nupa]][!method_types %in% multiclass_method]
      } else if (nupa %in% c("D.hat", "D.hat.z") && K > 2) {
        result[[nupa]] <- method[[nupa]][!method_types %in% base_method]
      } else {
        result[[nupa]] <- method[[nupa]]
      }
    }
  }

  empty_nupas <- names(which(lengths(result) == 0))
  if (length(empty_nupas) > 0) {
    stop(paste("The following nuisance parameters have empty method lists:", 
               paste(empty_nupas, collapse = ", "), "\nPlease provide at least one valid method for each NuPa"))
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
#' method1 <- list(
#'   "ols" = create_method("ols"),
#'   "plasso" = create_method("plasso")
#' )
#' 
#' method2 <- list(
#'   "forest_grf" = create_method("forest_grf"),
#'   "xgboost" = create_method("xgboost")
#' )
#' 
#' np1 <- nuisance_parameters(
#'   NuPa = c("Y.hat", "D.hat"), X = X, Y = Y, D = D,
#'   method = method1, cf = 5, stacking = "short", storeModels = "Memory"
#' )
#' 
#' np2 <- nuisance_parameters(
#'   NuPa = c("D.hat", "Z.hat"), X = X, D = D, Z = Z, 
#'   cf_mat = np1$numbers$cf_mat, method = method2, cf = 5, 
#'   stacking = "short", storeModels = "Memory"
#' )
#' 
#' np3 <- add_nupa(np = np1, np_new = np2, NuPa = c("D.hat", "Z.hat"), replace = TRUE)
#' }
#' 
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
      # Add the nuisance parameter and the method
      np$nuisance_parameters[[param]] <- np_new$nuisance_parameters[[param]]
      np$numbers$method[[param]] <- np_new$numbers$method[[param]]
      
      if (np$numbers$cv == 1) {
        # Combine the ensemble weights
        w_ens <- as.matrix(np[["numbers"]][["ens_weights"]])
        w_new_ens <- as.matrix(np_new[["numbers"]][["ens_weights"]])
        
        all_methods <- union(rownames(w_ens), rownames(w_new_ens))
        all_nupas <- union(colnames(w_ens), colnames(w_new_ens))
        
        ens_weights <- matrix(0, nrow = length(all_methods), ncol = length(all_nupas),
                              dimnames = list(all_methods, all_nupas))
        
        ens_weights[rownames(w_ens), colnames(w_ens)] <- w_ens
        ens_weights[, NuPa] <- 0
        ens_weights[rownames(w_new_ens), NuPa] <- w_new_ens[, NuPa, drop = FALSE]
        
        ens_weights <- as.data.frame(ens_weights)
        class(ens_weights) <- c("ens_weights_short", "data.frame")
        np[["numbers"]][["ens_weights"]] <- ens_weights
      } 
      else {
        
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
#' @param x A factor or character vector.
#' 
#' @return A numeric matrix with one column per unique category.
#' 
#' @keywords internal
#' 
one_hot <- function(x) {
  ux <- unique(x)
  out <- matrix(0L, nrow = length(x), ncol = length(ux))
  out[cbind(seq_along(x), match(x, ux))] <- 1L
  colnames(out) <- ux
  out
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
#' Tunable hyperparameters include learning rate (eta), tree depth, regularization terms,
#' and subsampling parameters. The search space is designed for general-purpose optimization
#' with log-uniform sampling for continuous parameters.
#'
#' @keywords internal
#' 
tune_xgb_hyperband <- function(X, Y, max_rounds = 100, n_configs = 50,
                               eta_downfactor = 3, nfold = 5,
                               metrics = "rmse", seed = 123) {
  set.seed(seed)

  sample_params <- function() {
    list(
      eta = 10^runif(1, log10(1e-4), log10(0.3)),
      max_depth = sample(2:10, 1),
      min_child_weight = sample(c(1, 2, 5, 10), 1),
      colsample_bytree = runif(1, 0.5, 1),
      colsample_bylevel = runif(1, 0.5, 1),
      # tree_method = "hist",
      lambda = 10^runif(1, log10(1e-3), log10(10)),
      subsample = 1,
      alpha = 0, 
      max_delta_step = 0, 
      base_score = 0,
      objective = "reg:squarederror"
    )
  }

  configs <- replicate(n_configs, sample_params(), simplify = FALSE)
  budgets <- c(10, 30, max_rounds)

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
      scores[i] <- if (!is.null(cv)) min(cv$evaluation_log[[paste0("test_", metrics, "_mean")]]) else NA
    }

    # Remove configs with NA scores
    valid_idx <- which(!is.na(scores))
    if (length(valid_idx) == 0) break # nothing left

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

  return(list("params" = best_params, "best_score" = best_score))
}


#' Tune learner hyperparameters
#'
#' Automatically tunes hyperparameters for specified learners in a method list
#'
#' @param method List of methods for ensemble estimation.
#' @param type Tuning type, either "\code{tune_full_sample}" or "\code{tune_on_fold}".
#' @param X Covariate matrix.
#' @param Y Numeric vector containing the outcome variable.
#' 
#' @return Modified method list with tuned parameters.
#'
#' @keywords internal
#' 
tune_learners <- function(type = c("tune_full_sample", "tune_on_fold"),
                          method,
                          X, Y) {
  type <- match.arg(type)

  # Process each method in the list
  for (i in seq_along(method)) {
    current_method <- method[[i]]

    if (!is.null(current_method$arguments) && identical(current_method$arguments, type)) {
      
      if (is.null(current_method$x_select)) {
        X_sub <- X
      } else {
        X_sub <- X[, current_method$x_select, drop = FALSE]
      }
      
      # GRF forest tuning
      if (identical(current_method$method, "forest_grf")) {
        tuned_forest <- grf::regression_forest(
          X = X_sub,
          Y = Y,
          tune.parameters = "all",
          num.trees = 100
        )

        method[[i]][["arguments"]] <- tuned_forest[["tuning.output"]][["params"]]
      }

      # XGBoost tuning using hyperband
      else if (identical(current_method$method, "xgboost")) {
        tuned_xgb <- tune_xgb_hyperband(X = X_sub, Y = Y)

        method[[i]][["arguments"]] <- tuned_xgb[["params"]]
      }
    }
  }

  return(method)
}
