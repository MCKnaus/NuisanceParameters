#' Creates design matrix
#'
#' @description
#' \code{\link{design_matrix}} creates interactions and polynomials for a proper
#' design matrix that can be used as input for prediction functions.
#'
#' @param data Matrix with the main variables
#' @param int Vector of strings with the variables to be interacted
#' @param int_d Degree of interactions created. Default is 2.
#' @param poly Vector of strings with the variables for which polynomials should be created
#' @param poly_d Degree of polynomials to be created. Default is 2.
#' @param log Vector of strings with the variables for which the logged versions should be added
#'
#' @return Matrix including the main variables and the newly generated ones.
#'
#' @importFrom matrixStats colMins
#' @importFrom stats model.matrix as.formula
#'
#' @export
#'
design_matrix = function(data, int = NULL, int_d = 2, poly = NULL, poly_d = 2, log = NULL) {

  # interactions
  if (!is.null(int)) {
    int = paste0(" + (", paste0(int, collapse = " + "), ")^", toString(int_d))
  }

  # polynomials
  if (!is.null(poly)) {
    poly = paste0(" + poly(", paste0(poly, collapse = paste0(",", toString(poly_d), ", raw = TRUE) + poly(")), ",", toString(poly_d), ", raw = TRUE)",
                   "- (", paste0(paste0(poly, collapse = " + ")),")")
  }

  # logs
  if (!is.null(log)) {
    # check whether some variables cannot be logged because not positive
    ind_neg = matrixStats::colMins(data[, log]) <= 0
    if (sum(ind_neg) > 0) {
      cat("\n Following variables not modified to be logged because of non-positive values:", paste(colnames(data[, log])[ind_neg]), "\n")
      log = log[!ind_neg]
    }
    log = paste0(" + log(", paste0(log, collapse = ") + log("), ")")
  }

  # combine the three parts
  fmla = stats::as.formula(paste(" ~ 0", int, poly, log))

  # generate matrix
  data = stats::model.matrix(fmla, data = as.data.frame(data))

  # adjust variable names
  colnames(data) = gsub("poly\\(", "", colnames(data))
  colnames(data) = gsub(paste0(", ", toString(poly_d), ", raw = TRUE)") ,"", colnames(data))
  colnames(data) = gsub("log\\(","ln_", colnames(data))
  colnames(data) = gsub("\\)", "", colnames(data))

  return(data)

}

#' Data screening
#'
#' @description
#' \code{\link{data_screen}} takes a matrix of data and cleans it for the sake
#' of an improved input for subsequent prediction algorithms. \cr
#' It removes \cr
#' 1. Variables without variation \cr
#' 2. Dummy variables where one group is nearly empty (optionally in one of both treatment groups) \cr
#' 3. Redundant (i.e. highly correlated) variables \cr
#'
#' @param data Matrix the variables to be screened
#' @param treat Optional binary treatment vector if screening should be done within treatment groups
#' @param bin_cut Cut-off fraction under which nearly empty binary variables should be removed. Default 0.01.
#' @param corr_cut Cut-off above which highly correlated variables should be removed. Default 0.99.
#' @param quiet If False, details about the removed variables are shown at each step
#'
#' @return Screened matrix.
#'
#' @importFrom matrixStats colSds
#' @importFrom Matrix colMeans colSums
#' @importFrom stats cor
#'
#' @export
#'
data_screen = function(data, treat = NULL, bin_cut = 0.01, corr_cut = 0.99, quiet = TRUE) {


  ### eliminate variables with no variation ###

  # identify the names
  nm_del = colnames(data)[matrixStats::colSds(data) == 0]
  if (isFALSE(quiet)) {
    cat("\n\n Variables with no variation:", nm_del, "\n\n")
  }
  # remove identified variables
  if (identical(nm_del, character(0)) == FALSE) data = data[, !colnames(data) %in% nm_del]


  ### remove dummy variables lower than threshold in one of the two treatment groups ###

  # identify dummies
  bin = apply(data,2,function(x) { all(x %in% 0:1) })

  # calculate means of all variables and check whether they are potentially close to 0 or 1
  if (is.null(treat)) {
    mean = Matrix::colMeans(data)
    bel_cut = (mean < bin_cut | mean > (1 - bin_cut))
  } else {
    mean1 = Matrix::colMeans(data[treat == 1, ])
    mean0 = Matrix::colMeans(data[treat == 0, ])
    bel_cut = (mean1 < bin_cut | mean1 > (1 - bin_cut) | mean0 < bin_cut | mean0 > (1 - bin_cut))
  }

  nm_del = colnames(data)[bin & bel_cut]
  if (isFALSE(quiet)) {
    cat("\n\n Dummy variables close to 0 or 1:", nm_del, "\n\n")
  }
  # remove identified variables
  if (identical(nm_del, character(0)) == FALSE) data = data[, !colnames(data) %in% nm_del]


  ### remove all redundant variables ###

  # calculate correlation matrix and consider only upper diagonal
  cor_mat = (abs(cor(data)) > corr_cut)
  cor_mat[lower.tri(cor_mat, diag = TRUE)] = FALSE

  nm_del = colnames(cor_mat)[Matrix::colSums(cor_mat) > 0]
  if (isFALSE(quiet)) {
    cat("\n\n Variables (nearly) perfectly correlated:", nm_del, "\n\n")
  }
  # remove identified variables
  if (identical(nm_del, character(0)) == FALSE) data = data[, Matrix::colSums(cor_mat) == 0]

  return(data)

}
