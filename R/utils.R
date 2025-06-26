#' One-hot encoding of treatment vector
#'
#' @description
#' \code{\link{prep_w_mat}} creates a matrix of binary treatment indicators
#' (n x T+1 with T being the number of treatments)
#' where each column contains a binary indicator for one treatment (i.e. one-hot
#' encoding).
#'
#' @param w Treatment vector.
#' Provide as factor to control ordering of the treatments,
#' otherwise arranged in ascending order or alphabetically.
#'
#' @return Logical matrix of treatment indicators (n x T+1).
#'
#' @importFrom stats model.matrix
#'
#' @export
#'
#' @examples
#' w = factor(c("A", "B", "A", "C", "B", "C"))
#' w_mat = prep_w_mat(w)
#' head(w_mat)
#'
prep_w_mat = function(w) {
  # Checks
  if (length(unique(w)) <= 1) stop("Need at least two values in treatment vector")
  if (!is.factor(w)) w = factor(w, levels = sort(unique(w)))

  # Create one-hot matrix for each category
  w_mat = stats::model.matrix(~ 0 + w)
  colnames(w_mat) = gsub("w", "", colnames(w_mat))
  return(w_mat == 1)
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
#' indicators (n x # cross-folds)
#'
#' @param n Number of observations
#' @param cf Number of cross-fitting folds
#' @param cl Optional vector of cluster variable if cross-fitting should account
#' for clusters within the data.
#' @param w_mat Optional logical matrix of treatment indicators (n x T+1 with T
#' being the number of treatments).
#' For example created by \code{\link{prep_w_mat}}.
#' If specified, cross-fitting folds will preserve the treatment ratios from full sample.
#' However, if cluster vector is provided, w_mat is ignored due to computational
#' constraints and randomization as well as feasibility issues.
#'
#' @return Logical matrix of cross-fitting folds (n x # folds).
#'
#' @importFrom stats model.matrix quantile
#'
#' @export
#'
#' @examples
#' n = 1000
#' w_mat = prep_w_mat(sample(3, n, replace = TRUE))
#' cf_mat = prep_cf_mat(n, cf = 5, w_mat = w_mat)
#' head(cf_mat)
#'
prep_cf_mat = function(n, cf, cl = NULL, w_mat = NULL) {

  # check if both cluster vector and treatment matrix are provided
  if (!is.null(cl) & !is.null(w_mat)) {
    warning("You provided both a cluster vector and a treatment matrix. This is not feasible due to computational constraints and randomization issues. Thus, only cluster vector is considered.")
    w_mat = NULL
  }

  # only one fold (i.e. no cross-fitting)
  if (cf == 1) {

    cf_mat = matrix(rep(1, n), ncol = 1)

  # neither treatment matrix nor cluster vector provided
  } else if (is.null(w_mat) & is.null(cl)) {

    rnd_id = sample(1:n, n)
    fold = factor(as.numeric(cut(rnd_id, breaks = stats::quantile(rnd_id, probs = seq(0, 1, length = cf + 1)), include.lowest = TRUE)))
    cf_mat = (stats::model.matrix(~ 0 + fold) == 1)

  # treatment matrix but no cluster vector provided
  } else if (!is.null(w_mat) & is.null(cl)) {

    cf_mat = matrix(NA, nrow = n, ncol = cf)
    nw = colSums(w_mat)

    for (i in 1:ncol(w_mat)) {

      cf_mat_w = matrix(FALSE, nrow = nw[i], ncol = cf)
      rnd_id = sample(1:nw[i], nw[i])
      fold = as.numeric(cut(rnd_id, breaks = stats::quantile(rnd_id, probs = seq(0, 1, length = cf + 1)), include.lowest = TRUE))

      for (j in 1:cf) {

        cf_mat_w[fold == j, j] = TRUE

      }

      cf_mat[w_mat[, i], ] = cf_mat_w

    }

  # no treatment matrix but cluster vector provided
  } else if (is.null(w_mat) & !is.null(cl)) {

    check_cluster_compatibility(cl, cf)

    rnd_id = sample(1:length(unique(cl)), length(unique(cl)))
    fold = as.numeric(cut(rnd_id, breaks = stats::quantile(rnd_id, probs = seq(0, 1, length = cf + 1)), include.lowest = TRUE))
    fold = factor(fold[match(cl, unique(cl))])
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
#' @param y A numeric vector of actual target values
#'
#' @return A numeric vector of the non-negative least weights.
#'
#' @importFrom nnls nnls
#'
#' @export
#'
#' @examples
#' X = matrix(c(0.2, 0.3, 0.4, 0.5, 0.6, 0.7), ncol = 2)
#' y = c(0.25, 0.45, 0.65)
#' \donttest{nnls_w = nnls_weights(X, y)}
#'
nnls_weights = function(X, y) {
  nnls_result = nnls::nnls(X, y)
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
#' @param n Number of observations.
#' @param learner Vector of characters indicating whether to use S or T learner
#' or both.
#'
#' @return A matrix of NAs with dimensions \code{n} x \code{length(ml)}
#'
#' @keywords internal
#'
make_fit_cv = function(ml, n, learner = c("t", "s", "both")) {

  learner = match.arg(learner)

  fit_cv = matrix(NA, n, length(ml))
  colnames(fit_cv) = sprintf("method%s", seq(1:length(ml)))
  for (i in 1:length(ml)) {
    if (!is.null(ml[[i]]$name)) colnames(fit_cv)[i] = ml[[i]]$name
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


#' Visualize ensemble weights
#'
#' @description
#' This function plots the non-negative ensemble learner weights for both either
#' a short-stacked or standard-stacked ensemble learner.
#'
#' @param x Object containing predictions by ensemble learner (fit_cv)
#' and the cross-validated ensemble weights.
#' @param methods Vector of method names. In order of initial \code{\link{create_method}}
#' specification.
#' @param legend_pos Legend position.
#' @param legend_size Font size of legend.
#' @param ... Pass generic \code{\link[base]{plot}} options.
#' @return Plot object.
#'
#' @export
#'
#' @method plot ens.learner
#'
plot.ens.learner = function(x,
                            methods = NULL,
                            legend_pos = c("center", "bottomright", "bottom", "bottomleft", "left", "topleft", "top", "topright", "right"),
                            legend_size = 1,
                            ...) {

  legend_pos = match.arg(legend_pos)

  if (is.null(names(x[[1]]))) {
    w = matrix(x$nnls_w, nrow = 1)
    colnames(w) = names(x$nnls_w)
    xlab = ""
    x_labels = NULL
  } else {
    w = do.call(rbind, lapply(x, function(l) l$nnls_w))
    xlab = "Cross-Fitting Fold"
    x_labels = 1:nrow(w)
  }

  if(!is.null(methods)) colnames(w) = methods
  w_t = t(w)

  colors = grDevices::rainbow(ncol(w))
  names(colors) = colnames(w)

  graphics::barplot(w_t, beside = FALSE, col = colors,
                    names.arg = x_labels,
                    xlab = xlab, ylab = "Ensemble Weight", main = "Cross-Validated Ensemble Weights",
                    ...)
  graphics::legend(x = legend_pos, title = "Method", xpd = NA, legend = names(colors), col = colors,
                   ncol = 1, cex = legend_size, lty = 1)

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
  
  # Simplified center-align without ANSI interference
  format_center <- function(x, width) {
    if (getOption("knitr.in.progress", FALSE)) {
      return(x)  # Skip fancy formatting in notebooks
    }
    x <- substr(as.character(x), 1, width)
    pad <- width - nchar(x)
    paste0(strrep(" ", ceiling(pad / 2)), x, strrep(" ", floor(pad / 2)))
  }
  
  # Force plain text in notebooks
  if (getOption("knitr.in.progress", FALSE)) {
    pb$format <- gsub("\\033\\[[0-9;]*[mGKH]", "", pb$format)
  }
  
  pb$tick(tokens = list(
    nuisance = format_center(pb_np, 11),
    pb_cf    = format_center(pb_cf, 2),
    pb_cv    = format_center(pb_cv, 2), 
    task     = format_center(task, 7),
    model    = format_center(method, 10)
  ))
}