#' Nuisance Parameter Estimation
#'
#' Estimates nuisance parameters using ensemble methods.
#'
#' @param NuPa Character vector specifying the nuisance parameters to estimate.
#'             Currently supported options: 
#'             \code{c("Y.hat", "Y.hat.d", "Y.hat.z", "D.hat", "D.hat.z", "Z.hat")}
#' @param method List of methods to use for \code{\link{ensemble}} estimation.
#'               Methods can be created using \code{\link{create_method}}.
#' @param Y Numeric vector containing the outcome variable.
#' @param D Numeric vector containing the treatment variable.
#' @param Z Numeric vector containing the instrument variable.
#' @param X Covariate matrix.
#' @param cf Number of cross-fitting folds (default: 5).
#' @param cl Optional vector of cluster variables if cross-fitting should account
#'           for clusters within the data.
#' @param stratify Logical. If \code{TRUE}, cross-fitting folds will preserve the
#'                 treatment ratios from the full sample. Note: If cluster vector
#'                 is provided, stratification is ignored due to computational
#'                 constraints and randomization feasibility issues.
#' @param stacking Stacking option: Either "short" for short stacking or an integer >1
#'                 indicating the number of cross-validation folds within each
#'                 cross-fitting fold.
#' @param storeModels Character vector indicating whether to save individual
#'                    models for future processing (default: "No").
#' @param path Optional path to save the \code{\link{ensemble}} objects for
#'             later processing (saved as a list of models).
#' @param quiet Logical. If \code{FALSE}, progress output is printed to the console.
#'
#' @return A list containing:
#' \itemize{
#'   \item \code{nuisance_parameters}: Requested nuisance parameter estimates
#'   \item \code{models} (optional): \code{\link{ensemble}} objects (individual models)
#'    \item \code{numbers}: Additional objects used for downstream processing in
#'                          smoother matrices and outcome weights extraction
#' }
#'
#' @export
nuisance_parameters = function(NuPa = c("Y.hat","Y.hat.d","Y.hat.z","D.hat","D.hat.z","Z.hat"),
                               X, Y, D, Z,
                               method,
                               cluster = NULL,
                               stratify = FALSE,
                               cf, stacking,
                               storeModels = c("No", "Memory", "Disk"),
                               path = NULL,
                               quiet = TRUE) {
  
  
  ## Sanity checks
  storeModels = match.arg(storeModels)

  # Define path if not provided
  if (is.null(path)) path = getwd()
  
  supported_NuPas = c("Y.hat","Y.hat.d","Y.hat.z","D.hat","D.hat.z","Z.hat")
  not_supported = NuPa[!NuPa %in% supported_NuPas]
  
  if (length(not_supported) > 0) {
    stop(paste("Error: The following nuisance parameters specified in NuPa are not supported:", 
               paste(not_supported, collapse = ", ")))}
  if (is.null(Y) & any(NuPa %in% c("Y.hat","Y.hat.d","Y.hat.z"))) {
    stop("Please specify Y if at least one of c(\"Y.hat\",\"Y.hat.d\",\"Y.hat.z\") is specified in NuPa")}
  if (is.null(D) & any(NuPa %in% c("Y.hat.d","D.hat","D.hat.z"))) {
    stop("Please specify D if at least one of c(\"Y.hat.d\",\"D.hat\",\"D.hat.z\") is specified in NuPa")}
  if (is.null(Z) & any(NuPa %in% c("Y.hat.z","D.hat.z","Z.hat"))) {
    stop("Please specify Z if at least one of c(\"Y.hat.z\",\"D.hat.z\",\"Z.hat\") is specified in NuPa")}
  
  if (stacking == "short") { cv <- 1L
  } else if (is.numeric(stacking) && stacking > 1 && stacking == round(stacking)) {
    cv <- as.integer(stacking)
  } else { stop("`stacking` must be either 'short' or an integer > 1") }
  
  if (length(method) == 1) {cv <- 1
    message("With a single learner provided, estimation defaults to cross-fitted predictions.")
  }
  

  ## Preps
  N = nrow(X)
  d_mat = prep_indicator_mat(D)
  z_mat = prep_indicator_mat(Z)
  cf_mat = prep_cf_mat(N, cf = cf, d_mat = if (stratify) d_mat else NULL, cluster = cluster)

  if (isFALSE(quiet)) which_stacking(cv)
  ens_weights <- list()
  
  ## Initialize objects
  # Initialize nuisance parameters for all
  Y.hat = Y.hat.d = Y.hat.z = D.hat = D.hat.z = Z.hat = 
    "This nuisance parameter was not specified and is therefore empty."
  
  Y.hat_m = Y.hat.d_m = Y.hat.z_m = 
    "This model was not specified and is therefore empty."
  
  # Initialize nuisance parameters and estimated model lists to be filled
  if ("Y.hat" %in% NuPa) {Y.hat = matrix(NA, nrow(d_mat), ncol(d_mat)-1); Y.hat_m <- vector("list", ncol(d_mat)-1)}
  if ("Y.hat.d" %in% NuPa) {Y.hat.d = matrix(NA, nrow(d_mat), ncol(d_mat)); Y.hat.d_m <- vector("list", ncol(d_mat))}
  if ("Y.hat.z" %in% NuPa) {Y.hat.z = matrix(NA, nrow(z_mat), ncol(z_mat)); Y.hat.z_m <- vector("list", ncol(z_mat))}
  if ("D.hat" %in% NuPa) D.hat = matrix(NA, nrow(d_mat), ncol(d_mat)-1)
  if ("D.hat.z" %in% NuPa) D.hat.z = matrix(NA, nrow(z_mat), ncol(z_mat))
  if ("Z.hat" %in% NuPa) Z.hat = matrix(NA, nrow(z_mat), ncol(z_mat)-1)
  
  if ("Y.hat.d" %in% NuPa) colnames(Y.hat.d) = colnames(d_mat)
  if ("Y.hat.z" %in% NuPa) colnames(Y.hat.z) = colnames(z_mat)
  if ("D.hat.z" %in% NuPa) colnames(D.hat.z) = colnames(z_mat)
  
  
  ## Progress bar
  setup_progress_bar <- function(NuPa, n_d, n_z, cf_folds, cv_folds, models) {
    
    # Base ticks for each nuisance parameter
    base_ticks <- c("Y.hat.d"=n_d,"Y.hat.z"=n_z,"D.hat.z"=n_z,"Y.hat"=1,"D.hat"=1,"Z.hat"= 1)
    total_ticks <- sum(base_ticks[names(base_ticks) %in% NuPa])
    
    # Final ticks: cf folds × models × fitting/prediction × cv folds
    total_ticks <- total_ticks * cf_folds * length(models) * 2 * if (cv_folds > 1 && length(models) > 1) (cv_folds + 1) else 1
    
    pb <- progress::progress_bar$new(
      format = "[:bar] :percent | :current/:total | :nuisance | cf =:pb_cf, cv =:pb_cv | :task :model",
      total = total_ticks, clear = TRUE, width = 80, force = TRUE)
    return(pb)
    }
  
  pb <- if (isFALSE(quiet)) {setup_progress_bar(NuPa, ncol(d_mat), ncol(z_mat), ncol(cf_mat), cv, method)} else {NULL}


  ######
  # Outcome NuPa
  if ("Y.hat.d" %in% NuPa) {
    
    for (i in seq_len(ncol(d_mat))) {
      pb_np <- paste0("Y.hat.d, d=", i-1)
      
      temp <- nuisance_cf(
        method = method, Y = Y, X = X, cf_mat = cf_mat, cv = cv, subset = d_mat[, i],
        storeModels = "Memory", path = NULL, quiet = quiet, pb = pb, pb_np = pb_np)
      
      Y.hat.d[, i] <- temp$np
      Y.hat.d_m[[i]] <- temp$models
      ens_weights[[paste0("Y.hat.d", i-1)]] <- temp$ens_weights
    }
  }

  if ("Y.hat.z" %in% NuPa) {
    
    for (i in seq_len(ncol(z_mat))) {
      pb_np <- paste0("Y.hat.z, z=", i-1)
      
      temp <- nuisance_cf(
        method = method, Y = Y, X = X, cf_mat = cf_mat, cv = cv, subset = z_mat[, i],
        storeModels = "Memory", path = NULL, quiet = quiet, pb = pb, pb_np = pb_np)
      
      Y.hat.z[, i] <- temp$np
      Y.hat.z_m[[i]] <- temp$models
      ens_weights[[paste0("Y.hat.z", i-1)]] <- temp$ens_weights
    }
  }
  
  if ("Y.hat" %in% NuPa) {
    pb_np <- "Y.hat"
    
    temp <- nuisance_cf(
      method = method, Y = Y, X = X, cf_mat = cf_mat, cv = cv, subset = NULL,
      storeModels = "Memory", path = NULL, quiet = quiet, pb = pb, pb_np = pb_np)
    
    Y.hat <- temp$np
    Y.hat_m <- temp$models
    ens_weights[["Y.hat"]] <- temp$ens_weights
  }
  
  
  ######
  # Treatment NuPa
  if ("D.hat.z" %in% NuPa) {
    
    for (i in seq_len(ncol(z_mat))) {
      pb_np <- paste0("D.hat.z, z=", i-1)
      
      temp <- nuisance_cf(
        method = method, Y = D, X = X, cf_mat = cf_mat, cv = cv, subset = z_mat[, i],
        storeModels = "No", path = NULL, quiet = quiet, pb = pb, pb_np = pb_np)
      
      D.hat.z[, i] <- temp$np
      ens_weights[[paste0("D.hat.z", i-1)]] <- temp$ens_weights
    }
  }
  
  if ("D.hat" %in% NuPa) {
    pb_np <- "D.hat"
    
    temp <- nuisance_cf(
      method = method, Y = D, X = X, cf_mat = cf_mat, cv = cv, subset = NULL, 
      storeModels = "No", path = NULL, quiet = quiet, pb = pb, pb_np = pb_np)
    
    D.hat <- temp$np
    ens_weights[["D.hat"]] <- temp$ens_weights
  }
  

  ######
  # Instrument NuPa
  if ("Z.hat" %in% NuPa) {
    pb_np <- "Z.hat"
    
    temp <- nuisance_cf(
      method = method, Y = Z, X = X, cf_mat = cf_mat, cv = cv, subset = NULL, 
      storeModels = "No", path = NULL, quiet = quiet, pb = pb, pb_np = pb_np)
    
    Z.hat <- temp$np
    ens_weights[["Z.hat"]] <- temp$ens_weights
  }
  
  # Combine collected ens_weights into a matrix
  if (cv == 1) {
    ens_weights_mat <- as.data.frame(lapply(ens_weights, `length<-`, max(lengths(ens_weights))))
    ens_weights_mat[is.na(ens_weights_mat)] <- 0
    class(ens_weights_mat) <- c("ens_weights_short", "data.frame")
  } else {
    ens_weights_mat <- format_weights(ens_weights)
    class(ens_weights_mat) = c("ens_weights_stand", "data.frame")
  }

  np_list = list("Y.hat"=Y.hat,"Y.hat.d"=Y.hat.d,"Y.hat.z"=Y.hat.z, 
                 "D.hat"=D.hat,"D.hat.z"=D.hat.z,"Z.hat"=Z.hat)
  
  models_list = list("Y.hat_m"=Y.hat_m,"Y.hat.d_m"=Y.hat.d_m,"Y.hat.z_m"=Y.hat.z_m)
  
  nums_list = list("N" = N,"cv" = cv,"cf_mat" = cf_mat, "d_mat" = d_mat, "z_mat" = z_mat,
                   "method" = method, "X" = X, "Y" = Y, "ens_weights" = ens_weights_mat)
  
  if (storeModels == "Disk") {
    models_path <- file.path(path, "nuisance_models.rds")
    disk_data <- list("models" = models_list, "numbers" = nums_list)
    saveRDS(disk_data, file = models_path)
    if (!quiet) message("Models saved to: ", models_path)
  }
  
  return(list(
    nuisance_parameters = np_list,
    models = if (storeModels == "Memory") models_list else NULL,
    numbers = nums_list
  ))
  
}


#' Cross-fitting of nuisance parameters
#'
#' @description
#' \code{\link{nuisance_cf}} makes a cross-fitted ensemble prediction using
#' \code{\link{ensemble}}.
#'
#' @param method List of methods to be used in \code{\link{ensemble}} estimation.
#' Methods can be created by \code{\link{create_method}}.
#' @param Y Vector of variable to be predicted.
#' @param X Matrix of covariates.
#' @param cf_mat Logical matrix with k columns of indicators representing the different folds
#' (for example created by \code{\link{prep_cf_mat}}).
#' @param cv Number of cross-validation folds when estimating ensemble model.
#' Default value is 1 which then evaluates to a short-stacking procedure which
#' is computationally less demanding than standard stacking.
#' @param subset Optional logical vector if only subset of data should be used for prediction.
#' @param storeModels Vector of characters indicating where to save individual 
#' models for future processing (default "No").
#' @param path Path to save fit_cv matrix and non-negative least square weights to.
#' Optionally, weights as well.
#' @param quiet If FALSE, progress output is printed into console.
#'
#' @return Returns a vector of length N containing nuisance parameters.
#'
#' @keywords internal
#'
nuisance_cf = function(method, Y, X, cf_mat,
                       cv = 5,
                       subset = NULL,
                       storeModels = c("No", "Memory", "Disk"),
                       path,
                       quiet = TRUE,
                       pb = NULL, pb_np = NULL) {


  ## Checks
  if (is.numeric(cf_mat)) {
    if (!all(cf_mat %in% 0:1) | !is.matrix(cf_mat)) stop("Please provide cf_mat as binary indicator matrix. E.g. use function prep_cf_mat")
    if (nrow(cf_mat) != length(Y)) stop("cf_mat indicator matrix nrows different from # of obs.")
    if (length(Y) != sum(cf_mat)) stop("cf_mat indicator matrix does not sum to number of observations.")
  }
  
  if (is.null(subset)) subset = rep(TRUE, length(Y))

  np = rep(NA, length(Y))
  models <- NULL
  
  
  ### Hyperparameter tuning
  if (!is.null(method$forest_grf) && identical(method$forest_grf$arguments, list("tune_full"))) {
    tuning <- grf::regression_forest(X = X, Y = Y, tune.parameters = "all")
    method$forest_grf$arguments <- as.list(tuning$tuning.output$params)
  }
  

  ### Short-Stacking ###

  if (cv == 1) {
      
      ens <- ensemble_short(
        method = method, X = X, Y = Y, subset = subset, cf_mat = cf_mat,
        storeModels = storeModels, path = path, quiet = quiet, pb = pb, pb_np = pb_np)
      
      nnls_w <- if (length(method) == 1) 1 else nnls_weights(ens$fit_cv[subset, ], Y[subset])
      np <- predict(ens, w = nnls_w)
      
      fit_sub <- list("ens_object" = ens, "nnls_w" = nnls_w)
      class(fit_sub) <- "ens.learner"
      
      ens_weights <- nnls_w
      models <- fit_sub
  

  ### Standard-Stacking ###

  } else if (cv > 1) {

    fit_sub = ens_weights = list()

    for (i in 1:ncol(cf_mat)) {

      fold = cf_mat[, i]
      X_tr = X[!fold & subset, ]
      Y_tr = Y[!fold & subset]
      X_te = X[fold, ]
      method_fold = method
      
      ### Hyperparameter tuning
      if (!is.null(method_fold$forest_grf) && identical(method_fold$forest_grf$arguments, list("tune_fold"))) {
        tuning <- grf::regression_forest(X = X_tr, Y = Y_tr, tune.parameters = "all")
        method_fold$forest_grf$arguments <- as.list(tuning$tuning.output$params)
      }

      ens = ensemble(method = method_fold, X = X_tr, Y = Y_tr, nfolds = cv, 
                     quiet = quiet, pb = pb, pb_np = pb_np, pb_cf = i
                     )
      
      ens_pred = predict(ens, method_fold, X = X_tr, Y = Y_tr, Xnew = X_te, quiet = quiet, pb = pb, pb_np = pb_np, pb_cf = i)
      np[fold] = ens_pred$np

      fit_sub[[i]] = list("fit_cv" = ens_pred$fit_cv, "nnls_w" = ens$nnls_weights, "ens_object" = ens)
      ens_weights[[i]] <- ens$nnls_weights
    }

    class(fit_sub) = "ens.learner"
    models <- fit_sub

  }
  
  output <- list(
    "np" = np,
    "models" = models,
    "ens_weights" = ens_weights
  )

  return(output)

}