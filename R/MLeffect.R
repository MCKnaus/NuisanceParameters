#' Double ML estimators with outcome smoothers
#'
#' Existing Double ML implementations are too general to easily extract smoother matrices
#' required to be compatible with the get_forest_weights() method. This motivates yet 
#' another Double ML implementation.
#'
#' @param Y Numeric vector containing the outcome variable.
#' @param W Optional binary treatment variable.
#' @param X Covariate matrix with N rows and p columns.
#' @param Z Optional binary instrumental variable.
#' @param estimators String (vector) indicating which estimators should be run.
#' Current menu: c("PLR","PLR_IV","AIPW_ATE","Wald_AIPW")
#' @param smoother Indicate which smoother to be used for nuisance parameter estimation.
#' Currently only available option \code{"honest_forest"} from the \pkg{grf} package.
#' @param n_cf_folds Number of cross-fitting folds. Default is 5.
#' @param n_reps Number of repetitions of cross-fitting. Default is 1.
#' @param ... Options to be passed to smoothers.
#' 
#' @return A list with three entries:
#' - \code{results}: a list storing the results, influence functions, and score functions of each estimator
#' - \code{NuPa.hat}: a list storing the estimated nuisance parameters and the outcome smoother matrices
#' 
#' @examples
#' \donttest{
#' # Sample from DGP borrowed from grf documentation
#' n = 200
#' p = 5
#' X = matrix(rbinom(n * p, 1, 0.5), n, p)
#' Z = rbinom(n, 1, 0.5)
#' Q = rbinom(n, 1, 0.5)
#' W = Q * Z
#' tau =  X[, 1] / 2
#' Y = rowSums(X[, 1:3]) + tau * W + Q + rnorm(n)
#' 
#' # Run outcome regression and extract smoother matrix
#' # Run DML and look at results
#' dml = dml_with_smoother(Y,W,X,Z)
#' results_dml = summary(dml)
#' plot(dml)
#' 
#' # Get weights
#' omega_dml = get_outcome_weights(dml)
#' 
#' # Observe that they perfectly replicate the original estimates
#' all.equal(as.numeric(omega_dml$omega %*% Y), 
#'           as.numeric(as.numeric(results_dml[,1])))
#'
#' # The weights can then be passed to the cobalt package for example.
#' }
#' 
#' @references 
#' Chernozhukov, V., Chetverikov, W., Demirer, M., Duflo, E., Hansen, C., Newey, W., & Robins, J. (2018). 
#' Double/debiased machine learning for treatment and structural parameters. The Econometrics Journal, 21(1), C1-C68.
#'     
#' Knaus, M. C. (2024). Treatment effect estimators as weighted outcomes, \url{https://arxiv.org/abs/2411.11559}.
#'      
#' @export
#' 
MLeffect = function(Y,D,X,Z=NULL,
                    NuPa.hat,
                    estimators = c("PLR","PLR_IV","AIPW_ATE","Wald_AIPW"),
                    ...) {
  
  ## Sanity checks
  supported_estimators = c("PLR","PLR_IV","AIPW_ATE","Wald_AIPW")
  not_supported = estimators[!estimators %in% supported_estimators]
  if (length(not_supported) > 0) {
    stop(paste("Error: The following specified estimators are not supported:", 
               paste(not_supported, collapse = ", ")))}
  
  supported_NuPa = c("Y.hat","Y.hat.d","Y.hat.z","D.hat","D.hat.z","Z.hat")
  not_supported = names(NuPa.hat)[!names(NuPa.hat) %in% supported_NuPa]
  if (length(not_supported) > 0) {
    stop(paste("Error: The following passed nuisance parameters are not supported:", 
               paste(not_supported, collapse = ", ")))}

  if (any(c("PLR_IV", "Wald_AIPW") %in% estimators) && is.null(Z)) {
    stop("Error: Z cannot be NULL when using either 'PLR_IV' or 'Wald_AIPW' as an estimator.")
  }
  
  # Extract useful information
  N = length(Y)
  n_estimators = length(estimators)
  n_reps = 1
  
  # Intialize empty
  dml_plr = dml_PLR_IV = dml_AIPW_ATE = dml_Wald_AIPW = 
    "This estimator was not run."
  
  # Run the specified DML estimators
  if ("PLR" %in% estimators) {
    D.hat = NuPa.hat$D.hat
    Y.hat = NuPa.hat$Y.hat
    psi.a = -(D - D.hat)^2
    psi.b = (Y - Y.hat) * (D - D.hat)
    dml_plr = dml_inference(psi.a,psi.b)
  }
  if ("PLR_IV" %in% estimators) {
    D.hat = NuPa.hat$D.hat
    Y.hat = NuPa.hat$Y.hat
    Z.hat = NuPa.hat$Z.hat
    psi.a = -(D - D.hat) * (Z - Z.hat)
    psi.b = (Y - Y.hat) * (Z - Z.hat) 
    dml_PLR_IV = dml_inference(psi.a,psi.b)
  }
  if ("AIPW_ATE" %in% estimators) {
    D.hat = NuPa.hat$D.hat
    Y.hat.d0 = NuPa.hat$Y.hat.d[,1]
    Y.hat.d1 = NuPa.hat$Y.hat.d[,2]
    psi.a = matrix(-1,N)
    psi.b = Y.hat.d1 - Y.hat.d0 + D * (Y - Y.hat.d1) / D.hat - (1 - D) * (Y - Y.hat.d0) / (1-D.hat)
    dml_AIPW_ATE = dml_inference(psi.a,psi.b)
  }
  if ("Wald_AIPW" %in% estimators) {
    Z.hat = NuPa.hat$Z.hat
    D.hat.z0 = NuPa.hat$D.hat.z[,1]
    D.hat.z1 = NuPa.hat$D.hat.z[,2]
    Y.hat.z0 = NuPa.hat$Y.hat.z[,1]
    Y.hat.z1 = NuPa.hat$Y.hat.z[,2]
    psi.a = -( D.hat.z1 - D.hat.z0 + Z * (D - D.hat.z1) / Z.hat - (1 - Z) * (D - D.hat.z0) / (1-Z.hat) )
    psi.b = Y.hat.z1 - Y.hat.z0 + Z * (Y - Y.hat.z1) / Z.hat - (1 - Z) * (Y - Y.hat.z0) / (1-Z.hat)
    dml_Wald_AIPW = dml_inference(psi.a,psi.b)
  }
  
  list_results = list(
    "PLR" = dml_plr,
    "PLR_IV" = dml_PLR_IV,
    "AIPW_ATE" = dml_AIPW_ATE,
    "Wald_AIPW" = dml_Wald_AIPW )
  
  list_data = list(
    "Y" = Y,
    "D" = D,
    "Z" = Z,
    "X" = X )
  
  list_nums = list(
    "N" = N,
    "n_estimators" = n_estimators
  )
  
  output = list("results" = list_results,
                "NuPa.hat" = NuPa.hat,
                "data" = list_data,
                "numbers" = list_nums)
  
  class(output) = c("dml_with_smoother")
  
  return(output)
}


#' Results and influence functions for linear double ML estimators.
#'
#' @param psi.a psi.a component of linear Neyman orthogonal score
#' @param psi.b psi.a component of linear Neyman orthogonal score
#'
#' @return List of three components:
#' - \code{TaPa} Matrix storing estimate, standard error, t-value and p-value of target parameter of each repetition in the rows
#' - \code{IF} Matrix storing the influence function of the target parameter of each repetition in the columns
#' - \code{score} Array of dimension n_reps x N x 3 where the matrix in the first dimension stores the score in the first, 
#' psi.a in the second, and psi.b in the third column
#'
#' @keywords internal
#' @noRd
#'
dml_inference = function(psi.a, psi.b) {
  N = length(psi.a)  # Use length() instead of nrow() since inputs are vectors
  n_reps = 1
  
  # Initialize output containers
  TaPa = matrix(NA, n_reps, 4)
  colnames(TaPa) = c("Estimate", "SE", "t", "p")
  IF = matrix(NA, N, n_reps)
  score = array(NA, c(n_reps, N, 3))
  dimnames(score) = list(NULL, NULL, third_dim_names = c("psi", "psi.a", "psi.b"))
  
  theta = -sum(psi.b) / sum(psi.a)
  psi = theta * psi.a + psi.b
  IF[,1] = -psi / mean(psi.a)
  se = sqrt(var(IF[,1]) / N)
  t = theta / se
  p = 2 * pt(abs(t), N, lower.tail = FALSE)
  TaPa[1,] = c(theta, se, t, p)
  score[1,,] = cbind(psi, psi.a, psi.b)
  
  list("TaPa" = TaPa, "IF" = IF, "score" = score)
}