#' Double ML estimators (binary treatments only)
#' 
#' @description
#' For test use only, copied from \code{OutcomeWeights} package
#'
#' @param Y Numeric vector containing the outcome variable.
#' @param D Optional binary treatment variable.
#' @param X Covariate matrix with N rows and p columns.
#' @param Z Optional binary instrumental variable.
#' @param NuPa.hat a list of estimated nuisance parameters from \code{\link{nuisance_parameters}}
#' @param estimators String (vector) indicating which estimators should be run.
#' Current menu: c("PLR","PLR_IV","AIPW_ATE","Wald_AIPW")
#' 
#' @return An object of class `dml_with_smoother` containing:
#' \itemize{
#'   \item \code{results} - A list of estimated target parameters (PLR, PLR_IV, AIPW_ATE, Wald_AIPW)
#'   \item \code{NuPa.hat} - A list of input nuisance parameters
#'   \item \code{data} - Input data (Y, D, Z, X)
#'   \item \code{numbers} - Sample information (N, n_estimators)
#' }
#' @export
#' 
MLeffect = function(Y,D,X,Z=NULL,
                    NuPa.hat,
                    estimators = c("PLR","PLR_IV","AIPW_ATE","Wald_AIPW")
                    ) {
  
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
  se = sqrt(stats::var(IF[,1]) / N)
  t = theta / se
  p = 2 * stats::pt(abs(t), N, lower.tail = FALSE)
  TaPa[1,] = c(theta, se, t, p)
  score[1,,] = cbind(psi, psi.a, psi.b)
  
  list("TaPa" = TaPa, "IF" = IF, "score" = score)
}