#' Simulated example data
#'
#' @description
#' Simulated data based on the DGP of \code{\link[DoubleML]{make_plr_CCDDHNR2018}}.
#' Matrix of dimension 5000 X 12, i.e. 5000 observations and 12 variables of
#' which the first one refers to the target ("y"), the second one to the
#' treatment ("w") and the remaining 10 to the covariates ("x1" - "x10").
#' The latter 10 variables can be used for nuisance parameter estimation using
#' either "w" or "y" as target.
#' The data comes as \code{\link[base]{matrix}} class.
#'
#' @docType data
#'
#' @references https://cran.r-project.org/package=DoubleML
#'
#' @usage data(nuisance_data)
#'
#' @keywords datasets
#'
#' @examples
#' # load data
#' data(nuisance_data)
#' # extract target, treatment and covariates from data
#' y = nuisance_data[, "y", drop = FALSE]
#' w = nuisance_data[, "w", drop = FALSE]
#' X = nuisance_data[, 3:ncol(nuisance_data)]
#'
"nuisance_data"
