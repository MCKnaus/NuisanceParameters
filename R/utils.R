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

  if (max_cluster_share > (1/cf) * 0.9) {
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
#' @param w_mat Optional logical matrix of treatment indicators (n x T+1 with T
#' being the number of treatments).
#' For example created by \code{\link{prep_w_mat}}.
#' If specified, cross-fitting folds will preserve the treatment ratios from full sample.
#' @param cl Optional vector of cluster variable if cross-fitting should account for clusters.
#' If w_mat is provided, cluster vector is ignored due to computational and feasibility constraints
#' for the cross-fitting procedure to work.
#'
#' @return Logical matrix of cross-fitting folds (n x # folds).
#'
#' @importFrom stats model.matrix quantile
#'
#' @export
#'
prep_cf_mat = function(n, cf, w_mat = NULL, cl = NULL) {

  # only one fold (i.e. no cross-fitting)
  if (cf == 1) {

    cf_mat = matrix(rep(1,n), ncol = 1)

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

        cf_mat_w[fold == j,j] = TRUE

      }

      cf_mat[w_mat[,i],] = cf_mat_w

    }

  # no treatment matrix but cluster vector provided
  } else if (is.null(w_mat) & !is.null(cl)) {

    check_cluster_compatibility(cl, cf)

    rnd_id = sample(1:length(unique(cl)), length(unique(cl)))
    fold = as.numeric(cut(rnd_id, breaks = stats::quantile(rnd_id, probs = seq(0, 1, length = cf + 1)), include.lowest = TRUE))
    fold = factor(fold[match(cl, unique(cl))])
    cf_mat = (stats::model.matrix(~ 0 + fold) == 1)

    cf_mat_balance = colSums(cf_mat) / nrow(cf_mat)

    if (any(cf_mat_balance < (1/cf) * 0.5)) stop("There is a high imbalance in the cluster sizes. This poses a problem for the cross-fitting procedure. Either choose a smaller number of folds or do not specify a cluster vector.")

  } else if (!is.null(w_mat) & !is.null(cl)) {

    warning("You cannot provide both a treatment matrix and a cluster vector due to computational and feasibility constraints of the cross-fitting procedure.")

    cf_mat = matrix(NA, nrow = n, ncol = cf)
    nw = colSums(w_mat)

    for (i in 1:ncol(w_mat)) {

      cf_mat_w = matrix(FALSE, nrow = nw[i], ncol = cf)
      rnd_id = sample(1:nw[i], nw[i])
      fold = as.numeric(cut(rnd_id, breaks = stats::quantile(rnd_id, probs = seq(0, 1, length = cf + 1)), include.lowest = TRUE))

      for (j in 1:cf) {

        cf_mat_w[fold == j,j] = TRUE

      }

      cf_mat[w_mat[,i],] = cf_mat_w

    }

  }

  colnames(cf_mat) = sprintf("CF %d",1:cf)

  return(cf_mat)

}


#' Non-linear least weights function
#'
#' @description
#' \code{\link{nnls_weights}} calculates the non-linear least weights on the
#' basis of the non-negative least squares algorithm.
#'
#' @param X A matrix where each column represents a different predictor variable
#' and each row represents an observation
#' @param y A numeric vector of actual target values
#'
#' @return A numeric vector of non-linear least weights.
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

  return(nnls_w)
}


#' Adds an intercept to a matrix
#'
#' @description
#' \code{\link{add_intercept}} adds an intercept to a matrix.
#'
#' @param mat Any matrix (with column names).
#'
#' @return Matrix with intercept.
#'
#' @keywords internal
#'
add_intercept <- function(mat) {
  if (is.null(dim(mat))) {
    mat = as.matrix(mat,ncol = 1)
    colnames(mat) = "Var1"
  }

  if (all(mat[,1] == 1)) {
    colnames(mat)[1] == "(Intercept)"
    return(mat)
  } else {
    mat = cbind(rep(1,nrow(mat)),mat)
    colnames(mat)[1] = "(Intercept)"
    return(mat)
  }
}
