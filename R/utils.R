#' One-hot encoding of categorical variables
#'
#' @description
#' Creates a logical matrix of binary indicators (N x C) where each column 
#' represents a unique category from the input vector. Useful for encoding 
#' treatments, instruments, or other categorical variables.
#'
#' @param x Categorical vector (treatment, instrument, etc.). 
#'   Provide as factor to control category ordering, otherwise arranged in 
#'   ascending order or alphabetically.
#'
#' @return Logical matrix (N x C) where C is the number of unique categories in `x`.
#'   Each column is a binary indicator for one category.
#'
#' @importFrom stats model.matrix
#' @export
#'
#' @examples
#' D = factor(c("A", "B", "A", "C", "B", "C"))
#' d_mat = prep_indicator_mat(D)
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
#' \code{\link{check_cluster_compatibility}} checks if the cross-fitting procedure is feasible
#' given the cluster vector and desired number of cross-fitting folds.
#'
#' @param cl A vector representing the cluster assignments
#' @param cf The desired number of folds for cross-validation
#'
#' @return Does not return anything, only throws an error if checks fail.
#'
#' @export
#'
check_cluster_compatibility = function(cl, cf) {

  num_clusters = length(unique(cl))

  if (num_clusters < cf) {
    stop(
      "The number of clusters is less than the desired number of folds. Either choose a smaller number of folds or do not specify a cluster vector."
    )
  }

  cluster_shares = table(cl) / length(cl)
  max_cluster_share = max(cluster_shares)

  if (max_cluster_share > (1 / cf) * 0.9) {
    stop(
      "There is a high imbalance in the cluster sizes. This poses a problem for the cross-fitting procedure. Either choose a smaller number of folds or do not specify a cluster vector."
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
#' for clusters within the data.
#' @param d_mat Optional logical matrix of treatment indicators (N x T+1 with T
#' being the number of treatments).
#' For example created by \code{\link{prep_w_mat}}.
#' If specified, cross-fitting folds will preserve the treatment ratios from full sample.
#' However, if cluster vector is provided, d_mat is ignored due to computational
#' constraints and randomization as well as feasibility issues.
#'
#' @return Logical matrix of cross-fitting folds (N x # folds).
#'
#' @importFrom stats model.matrix quantile
#'
#' @export
#'
#' @examples
#' N = 1000
#' d_mat = prep_w_mat(sample(3, N, replace = TRUE))
#' cf_mat = prep_cf_mat(N, cf = 5, d_mat = d_mat)
#' head(cf_mat)
#'
prep_cf_mat = function(N, cf, cluster = NULL, d_mat = NULL) {

  # check if both cluster vector and treatment matrix are provided
  if (!is.null(cluster) & !is.null(d_mat)) {
    warning("You provided both a cluster vector and a treatment matrix. This is not feasible due to computational constraints and randomization issues. Thus, only cluster vector is considered.")
    d_mat = NULL
  }

  # only one fold (i.e. no cross-fitting)
  if (cf == 1) {

    cf_mat = matrix(rep(1, N), ncol = 1)

  # neither treatment matrix nor cluster vector provided
  } else if (is.null(d_mat) & is.null(cluster)) {

    rnd_id = sample(1:N, N)
    fold = factor(as.numeric(cut(rnd_id, breaks = stats::quantile(rnd_id, probs = seq(0, 1, length = cf + 1)), include.lowest = TRUE)))
    cf_mat = (stats::model.matrix(~ 0 + fold) == 1)

  # treatment matrix but no cluster vector provided
  } else if (!is.null(d_mat) & is.null(cluster)) {

    cf_mat = matrix(NA, nrow = N, ncol = cf)
    nw = colSums(d_mat)

    for (i in 1:ncol(d_mat)) {

      cf_mat_w = matrix(FALSE, nrow = nw[i], ncol = cf)
      rnd_id = sample(1:nw[i], nw[i])
      fold = as.numeric(cut(rnd_id, breaks = stats::quantile(rnd_id, probs = seq(0, 1, length = cf + 1)), include.lowest = TRUE))

      for (j in 1:cf) {

        cf_mat_w[fold == j, j] = TRUE

      }

      cf_mat[d_mat[, i], ] = cf_mat_w

    }

  # no treatment matrix but cluster vector provided
  } else if (is.null(d_mat) & !is.null(cluster)) {

    check_cluster_compatibility(cluster, cf)

    rnd_id = sample(1:length(unique(cluster)), length(unique(cluster)))
    fold = as.numeric(cut(rnd_id, breaks = stats::quantile(rnd_id, probs = seq(0, 1, length = cf + 1)), include.lowest = TRUE))
    fold = factor(fold[match(cluster, unique(cluster))])
    cf_mat = (stats::model.matrix(~ 0 + fold) == 1)

    cf_mat_balance = colSums(cf_mat) / nrow(cf_mat)

    if (any(cf_mat_balance < (1 / cf) * 0.5)) stop("There is a high imbalance in the cluster sizes. This poses a problem for the cross-fitting procedure. Either choose a smaller number of folds or do not specify a cluster vector.")

  }

  colnames(cf_mat) = sprintf("CF %d", 1:cf)

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
#' @importFrom nnls nnls
#'
#' @export
#'
#' @examples
#' X = matrix(c(0.2, 0.3, 0.4, 0.5, 0.6, 0.7), ncol = 2)
#' Y = c(0.25, 0.45, 0.65)
#' \donttest{nnls_w = nnls_weights(X, Y)}
#'
nnls_weights = function(X, Y) {
  nnls_result = nnls::nnls(X, Y)
  nnls_w = nnls_result$x

  # in case of perfectly agreeing predictions, nnls provides only zeros
  # in this scenario: uniform weights that add up to 1
  if (sum(nnls_w) == 0) {
    nnls_w = nnls_w + 1 / length(nnls_w)
  }

  # normalize weights
  nnls_w = nnls_w / sum(nnls_w)

  # assign names to weights
  if(!is.null(colnames(X))) names(nnls_w) = colnames(X)

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
add_intercept = function(mat) {
  if (is.null(dim(mat))) {
    mat = as.matrix(mat, ncol = 1)
  }
  if (all(mat[, 1] == 1)) {
    return(mat)
  } else {
    mat = cbind(rep(1, nrow(mat)), mat)
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
agg_array = function(a, w) {
  return(apply(a, c(1, 2), function(x) sum(x * w)))
}


#' Prepare fit_cv matrix
#'
#' @description
#' This function creates a matrix to store cross-fitted ensemble predictions.
#'
#' @param ml List of methods built via \code{\link{create_method}} to be used in
#' ensemble model
#' @param N Number of observations.
#' @param learner Vector of characters indicating whether to use S or T learner
#' or both.
#'
#' @return A matrix of NAs with dimensions \code{N} x \code{length(ml)}
#'
#' @keywords internal
#'
make_fit_cv = function(ml, N, learner = c("t", "s", "both")) {

  learner = match.arg(learner)

  fit_cv = matrix(NA, N, length(ml))
  colnames(fit_cv) = sprintf("method%s", seq(1:length(ml)))
  for (i in 1:length(ml)) {
    if (!is.null(names(ml))) colnames(fit_cv)[i] = names(ml)[i]
  }

  if(learner == "t") {
    colnames(fit_cv) = paste0("t-", colnames(fit_cv))
  } else if (learner == "s") {
    colnames(fit_cv) = paste0("s-", colnames(fit_cv))
  } else if (learner == "both") {
    ml_names = colnames(fit_cv)
    fit_cv = cbind(fit_cv, fit_cv)
    colnames(fit_cv) = c(paste0("t-", ml_names), paste0("s-", ml_names))
  }

  return(fit_cv)

}


#' Short- or Standard-stacking message
#'
#' @description
#' This function prints a message into the console whether short- or standard-
#' stacking is used.
#'
#' @param cv Cross-validation value for ensemble weights.
#'
#' @return Message in console. No object is returned.
#'
#' @keywords internal
#'
which_stacking = function(cv = 1) {
  if (cv == 1) {
    message("Short-stacking is used.")
  } else if (cv > 1) {
    message("Standard-stacking is used")
  }
}


#' Update progress bar
#' 
#' @description
#' This function updates progress bar with standardized formatting
#' 
#' @param pb Progress bar object from `progress` package
#' @param pb_np Current nuisance parameter (character)
#' @param pb_cf Current cross-fitting fold number (integer)
#' @param pb_cv Current cross-validation fold number (integer)
#' @param task Task description (character)
#' @param method Method name (character)
#' @keywords internal
update_progress <- function(pb, pb_np, pb_cf, pb_cv, task, method) {
  if (is.null(pb)) return(invisible(NULL))
  
  # Center-align formatter with padding
  format_center <- function(x, width) {
    x <- substr(as.character(x), 1, width)
    pad <- width - nchar(x)
    paste0(strrep(" ", ceiling(pad / 2)), x, strrep(" ", floor(pad / 2)))
  }
  
  pb$tick(tokens = list(
    nuisance = format_center(pb_np, 12),
    pb_cf    = format_center(pb_cf, 2),
    pb_cv    = format_center(pb_cv, 2), 
    task     = format_center(task, 7),
    model    = format_center(method, 10)
  ))
}


#' Format ensemble weights list for standard stacking
#'
#' @description 
#' Internal helper function that processes ensemble weights into a standardized 
#' format for cross-validated stacking (when cv > 1).
#'
#' @param ens_weights List of ensemble weights to be formatted
#' 
#' @return A list of data frames with formatted weights, where each element 
#'         corresponds to a NuPa (Y.hat.z0/Y.hat.z1) and contains fold-specific weights
#' 
#' @keywords internal
format_weights <- function(ens_weights) {
  weight_list <- list()
  for (NuPa in names(ens_weights)) {
    weight_matrix <- do.call(cbind, lapply(ens_weights[[NuPa]], unlist))
    colnames(weight_matrix) <- paste0("fold", seq_len(ncol(weight_matrix)))
    weight_list[[NuPa]] <- as.data.frame(weight_matrix)
  }
  return(weight_list)
}


#' Visualize standard-stacked ensemble weights
#'
#' @description
#' Plots cross-validated ensemble weights for standard stacking (multiple folds),
#' showing weight distributions across learners for each fold.
#'
#' @param x Object containing cross-validated ensemble weights
#' @param ncols Number of columns for facet grid (default: 3)
#' @param base_size Base font size for plot elements (default: 12)
#' @param ... Additional arguments passed to plotting functions
#'
#' @return A ggplot object showing weights by fold and learner
#'
#' @export
#'
#' @method plot ens_weights_stand
plot.ens_weights_stand <- function(x, 
                                   ncols = 3, 
                                   base_size = 12, 
                                   ...) {
  
  palette <- c("#FB8072", "#80B1D3",  "#FFED6F", "#BEBADA", 
               "#8DD3C7", "#FDB462", "#B3DE69", "#BC80BD")

  # Prepare plotting function (returns plot without legend)
  create_weight_plot <- function(model_name, model_data) {
    
    df <- as.data.frame(t(model_data))
    df$fold <- rownames(df)
    melted <- reshape2::melt(df, id.vars = "fold")
    
    ggplot2::ggplot(melted, ggplot2::aes(x = fold, y = value, fill = variable)) +
      ggplot2::geom_col(position = "stack", width = 0.7) +
      ggplot2::scale_fill_manual(values = palette) +
      ggplot2::scale_y_continuous(limits = c(0, 1.01), expand = c(0, 0)) +
      ggplot2::labs(title = model_name, x = NULL, y = NULL) +
      ggplot2::theme_minimal(base_size = base_size * 0.75) +
      ggplot2::theme(
        axis.text.x = ggplot2::element_text(angle = 45, hjust = 1, size = base_size * 0.6),
        plot.margin = ggplot2::unit(c(2, 2, 2, 2), "pt"), legend.position = "none"
      )
  }
  
  # Generate all plots (without legends)
  plot_list <- mapply(create_weight_plot, names(x), x, SIMPLIFY = FALSE)
  
  # Extract legend from first plot
  tmp_plot <- create_weight_plot(names(x)[1], x[[1]]) + 
    ggplot2::labs(fill = "Learners") +
    ggplot2::theme(
      legend.position = "bottom",
      legend.title = ggplot2::element_text(size = base_size * 0.8),
      legend.text = ggplot2::element_text(size = base_size * 0.7))
  
  n_plots <- length(plot_list)
  n_rows <- ceiling(n_plots / ncols)
  
  # Prevent empty slate printing 
  pdf(file = nullfile())
  legend <- gtable::gtable_filter(ggplot2::ggplotGrob(tmp_plot), "guide-box")
  plots <- gridExtra::arrangeGrob(grobs = plot_list, ncol = ncols, nrow = n_rows, padding = ggplot2::unit(0, "line"))
  invisible(dev.off())
  
  gridExtra::grid.arrange(plots, legend, nrow = 2, heights = c(0.85, 0.15))
}


#' Visualize short-stacked ensemble weights
#' 
#' @description 
#' Plots single set of ensemble weights for short stacking (no cross-validation),
#' showing the final weight assigned to each learner.
#'
#' @param x Object containing final ensemble weights
#' @param base_size Base font size for plot elements (default: 12)
#' @param ... Additional arguments passed to ggplot
#'
#' @return A ggplot object showing final learner weights
#'
#' @export
#'
#' @method plot ens_weights_short
plot.ens_weights_short <- function(x, 
                                   base_size = 12, 
                                   ...) {

  palette <- c("#FB8072", "#80B1D3",  "#FFED6F", "#BEBADA", 
               "#8DD3C7", "#FDB462", "#B3DE69", "#BC80BD")
  
  df <- as.data.frame(x)
  df$learner <- rownames(df)
  df_long <- reshape2::melt(df, id.vars = "learner")
  
  ggplot2::ggplot(df_long, ggplot2::aes(x = variable, y = value, fill = learner)) +
    ggplot2::geom_col(position = "stack", width = 0.7) +
    ggplot2::scale_fill_manual(values = palette) +
    ggplot2::labs(x = "Nuisance Parameter", y = "Weight", fill = "Learner") +
    ggplot2::theme_minimal(base_size = base_size) +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))
}
