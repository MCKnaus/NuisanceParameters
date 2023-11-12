#' Simulated example data
#'
#' @description
#' Simulated data from a DGP with covariates matrix X, treatment vector w
#' and target vector y.
#'
#' @docType data
#'
#' @usage data(nuisance_data)
#'
#' @keywords datasets
#'
#' @examples
#' # load data
#' data(nuisance_data)
#' # extract target, treatment and covariates from data
#' # y = nuisance_data[, "y", drop = FALSE]
#' # w = nuisance_data[, "w", drop = FALSE]
#' # X = nuisance_data[, 3:ncol(nuisance_data)]
#'
"nuisance_data"
