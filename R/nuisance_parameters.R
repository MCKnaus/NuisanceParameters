#' Nuisance parameter estimation
#'
#' @param NuPa String vector specifying the nuisance parameters to be estimated.
#' Currently supported: \code{c("Y.hat","Yw.hat","Yz.hat","W.hat","Wz.hat","Z.hat")} 
#' @param ml List of methods to be used in \code{\link{ensemble}} estimation.
#' Methods can be created by \code{\link{create_method}}.
#' @param y Numerical vector containing the outcome variable.
#' @param w Numerical vector containing the treatment variable.
#' @param z Numerical vector containing the instrument variable.
#' @param w_mat Logical matrix of treatment indicators (n x T+1). For example created by \code{\link{prep_w_mat}}.
#' @param z_mat Logical matrix of instrument indicators (n x T+1). For example created by \code{\link{prep_w_mat}}.
#' @param x Covariate matrix.
#' @param cf_mat Logical matrix with k columns of indicators representing the different folds
#' (for example created by \code{\link{prep_cf_mat}}).
#' @param cl Optional vector of cluster variable if cross-fitting should account
#' for clusters within the data.
#' @param cv Number of cross-validation folds when estimating ensemble (default 5).
#' @param learner Vector of characters indicating whether to use S or T learner
#' or a combination of the two.
#' @param storeModels Vector of characters indicating where to save individual 
#' models for future processing (default "No").
#' @param path Optional path to save the \code{\link{ensemble}} objects for later processing.
#' Saved as Ensemble_Yi where i is the number of the treatment in multiple treatment settings.
#' @param quiet If FALSE, progress output is printed into console.
#'
#' @return A list containing:
#' \itemize{
#'   \item{predictions:  }{Requested nuisance parameters}
#'   \item{models (optional):  }{\code{\link{ensemble}} objects (individual models)}
#' }
#' 
#' @export
nuisance_parameters = function(NuPa = c("Y.hat","Yw.hat","Yz.hat","W.hat","Wz.hat","Z.hat"),
                               ml,
                               x, y, w, z,
                               w_mat, z_mat, cf_mat,
                               cl = NULL,
                               cv = 5,
                               learner = c("t", "s", "both"),
                               storeModels = c("No", "Memory", "Disk"),
                               path = NULL,
                               quiet = TRUE) {
  
  
  ## Sanity checks
  learner = match.arg(learner)
  storeModels = match.arg(storeModels)
  
  supported_NuPas = c("Y.hat","Yw.hat","Yz.hat","W.hat","Wz.hat","Z.hat")
  not_supported = NuPa[!NuPa %in% supported_NuPas]
  if (length(not_supported) > 0) {
    stop(paste("Error: The following nuisance parameters specified in NuPa are not supported:", 
               paste(not_supported, collapse = ", ")))}
  if (is.null(y) & any(NuPa %in% c("Y.hat","Yw.hat","Yz.hat"))) {
    stop("Please specify Y if at least one of c(\"Y.hat\",\"Yw.hat\",\"Yz.hat\") is specified in NuPa")}
  if (is.null(w) & any(NuPa %in% c("Yw.hat","W.hat","Wz.hat"))) {
    stop("Please specify W if at least one of c(\"Yw.hat\",\"W.hat\",\"Wz.hat\") is specified in NuPa")}
  if (is.null(z) & any(NuPa %in% c("Yz.hat","Wz.hat","Z.hat"))) {
    stop("Please specify Z if at least one of c(\"Yz.hat\",\"Wz.hat\",\"Z.hat\") is specified in NuPa")}
  
  
  ## Preps
  n = nrow(x)
  
  # Define path if not provided
  if (is.null(path)) path = getwd()
  path_temp = paste0(path, "/Ensemble_Y", 1:ncol(w_mat))
  
  if((learner %in% c("s", "both")) & (ncol(w_mat) > 2)) stop("S-Learner cannot be combined with more than two treatments.")
  
  # Should be adjusted to accommodate different stacking for nuisance parameters
  if (isFALSE(quiet)) which_stacking(cv)
  
  
  ## Initialize objects
  # Initialize nuisance parameters for all
  Y.hat = Yw.hat = Yz.hat = W.hat = Wz.hat = Z.hat = 
    "This nuisance parameter was not specified and is therefore empty."
  
  # Initialize nuisance parameters and estimated model lists to be filled
  if ("Y.hat" %in% NuPa) {Y.hat = matrix(NA, nrow(w_mat), ncol(w_mat)-1); Y.hat_ml <- vector("list", ncol(w_mat)-1)}
  if ("Yw.hat" %in% NuPa) {Yw.hat = matrix(NA, nrow(w_mat), ncol(w_mat)); Yw.hat_ml <- vector("list", ncol(w_mat))}
  if ("Yz.hat" %in% NuPa) {Yz.hat = matrix(NA, nrow(z_mat), ncol(z_mat)); Yz.hat_ml <- vector("list", ncol(z_mat))}
  if ("W.hat" %in% NuPa) W.hat = matrix(NA, nrow(w_mat), ncol(w_mat)-1)
  if ("Wz.hat" %in% NuPa) Wz.hat = matrix(NA, nrow(z_mat), ncol(z_mat))
  if ("Z.hat" %in% NuPa) Z.hat = matrix(NA, nrow(z_mat), ncol(z_mat)-1)
  
  if ("Yw.hat" %in% NuPa) colnames(Yw.hat) = colnames(w_mat)
  if ("Yz.hat" %in% NuPa) colnames(Yz.hat) = colnames(z_mat)
  if ("Wz.hat" %in% NuPa) colnames(Wz.hat) = colnames(z_mat)
  
  
  ## Experimental: progress printing
  setup_progress_bar <- function(NuPa, n_w, n_z, cf_folds, cv_folds, models) {
    total_ticks <- 0
    
    # Calculate total steps
    if ("Yw.hat" %in% NuPa) total_ticks <- total_ticks + n_w
    if ("Yz.hat" %in% NuPa) total_ticks <- total_ticks + n_z
    if ("Wz.hat" %in% NuPa) total_ticks <- total_ticks + n_z
    if ("Y.hat" %in% NuPa) total_ticks <- total_ticks + 1
    if ("W.hat" %in% NuPa) total_ticks <- total_ticks + 1
    if ("Z.hat" %in% NuPa) total_ticks <- total_ticks + 1
    
    # Multiply by folds and models (fitting + prediction)
    total_ticks <- total_ticks * cf_folds * length(models) * 2 * if (cv_folds > 1 && length(models) > 1) (cv_folds+1) else 1
    
    pb <- progress::progress_bar$new(
      format = "[:bar] :percent | :current/:total | :nuisance | cf =:pb_cf, cv =:pb_cv | :task :model",
      total = total_ticks,
      clear = TRUE, width = 80, force = TRUE
      )
    
    return(pb)
  }
  
  # Initialize progress bar (only if quiet = FALSE)
  pb <- if (isFALSE(quiet)) {setup_progress_bar(NuPa, ncol(w_mat), ncol(z_mat), ncol(cf_mat), cv, ml)} else {NULL}
  
  
  ######
  # Outcome NuPa
  if ("Yw.hat" %in% NuPa) {
    
    for (i in seq_len(ncol(w_mat))) {
      pb_np <- paste0("Yw.hat, w=", i-1)
      
      temp <- nuisance_cf(
        ml = ml, y = y, x = x, cf_mat = cf_mat, learner = learner, cv = cv,
        subset = w_mat[, i], storeModels = storeModels,
        path = path_temp[i], quiet = quiet, pb = pb, pb_np = pb_np)
      
      Yw.hat[, i] <- temp$np
      Yw.hat_ml[[i]] <- temp$models
    }
  }
  
  if ("Yz.hat" %in% NuPa) {
    
    for (i in seq_len(ncol(z_mat))) {
      pb_np <- paste0("Yz.hat, z=", i-1)
      
      temp <- nuisance_cf(
        ml = ml, y = y, x = x, cf_mat = cf_mat, learner = learner, cv = cv,
        subset = z_mat[, i], storeModels = storeModels,
        path = path_temp[i], quiet = quiet, pb = pb, pb_np = pb_np)
      
      Yz.hat[, i] <- temp$np
      Yz.hat_ml[[i]] <- temp$models
    }
  }
  
  if ("Y.hat" %in% NuPa) {
    pb_np <- "Y.hat"
    
    temp <- nuisance_cf(
      ml = ml, y = y, x = x, cf_mat = cf_mat, learner = learner, cv = cv,
      subset = NULL, storeModels = storeModels,
      path = path_temp[1], quiet = quiet, pb = pb, pb_np = pb_np)
    
    Y.hat <- temp$np
    Y.hat_ml <- temp$models
  }
  
  
  ######
  # Treatment NuPa
  if ("Wz.hat" %in% NuPa) {
    
    for (i in seq_len(ncol(z_mat))) {
      pb_np <- paste0("Wz.hat, z=", i-1)
      
      temp <- nuisance_cf(
        ml = ml, y = w_mat[, i], x = x, cf_mat = cf_mat, learner = "t", cv = cv,
        subset = z_mat[, i], storeModels = "No",
        path = path_temp, quiet = quiet, pb = pb, pb_np = pb_np)
      
      Wz.hat[, i] <- temp$np
    }
  }
  
  if ("W.hat" %in% NuPa) {
    pb_np <- "W.hat"
    
    temp <- nuisance_cf(
      ml = ml, y = W, x = x, cf_mat = cf_mat, learner = "t", cv = cv,
      subset = NULL, storeModels = "No", 
      path = path_temp, quiet = quiet, pb = pb, pb_np = pb_np)
    
    W.hat <- temp$np
  }
  
  
  ######
  # Instrument NuPa
  if ("Z.hat" %in% NuPa) {
    pb_np <- "Z.hat"
    
    temp <- nuisance_cf(
      ml = ml, y = Z, x = x, cf_mat = cf_mat, learner = "t", cv = cv,
      subset = NULL, storeModels = "No", 
      path = path_temp, quiet = quiet, pb = pb, pb_np = pb_np)
    
    Z.hat <- temp$np
  }
  
  
  np_list = list("Y.hat"=Y.hat,"Yw.hat"=Yw.hat,"Yz.hat"=Yz.hat, 
                 "W.hat"=W.hat,"Wz.hat"=Wz.hat,"Z.hat"=Z.hat)
  models_list = list("Y.hat_ml"=Y.hat_ml,"Yw.hat_ml"=Yw.hat_ml,"Yz.hat_ml"=Yz.hat_ml)
  
  list("nuisance_parameters"=np_list, "models"=models_list)
  
}



#' Cross-fitting of nuisance parameters
#'
#' @description
#' \code{\link{nuisance_cf}} makes a cross-fitted ensemble prediction using
#' \code{\link{ensemble}}.
#'
#' @param ml List of methods to be used in \code{\link{ensemble}} estimation.
#' Methods can be created by \code{\link{create_method}}.
#' @param y Vector of variable to be predicted.
#' @param x Matrix of covariates.
#' @param cf_mat Logical matrix with k columns of indicators representing the different folds
#' (for example created by \code{\link{prep_cf_mat}}).
#' @param learner Vector of characters indicating whether to use S or T learner
#' or a combination of the two.
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
#' @return Returns a vector of length n containing nuisance parameters.
#'
#' @keywords internal
#'
nuisance_cf = function(ml, y, x, cf_mat,
                       learner = c("t", "s", "both"),
                       cv = 5,
                       subset = NULL,
                       storeModels = c("No", "Memory", "Disk"),
                       path,
                       quiet = TRUE,
                       pb = NULL, pb_np = NULL) {


  ### Parameter Configuration ###

  learner = match.arg(learner)



  ### Checks ###

  if (is.numeric(cf_mat)) {
    if (!all(cf_mat %in% 0:1) | !is.matrix(cf_mat)) stop("Please provide cf_mat as binary indicator matrix. E.g. use function prep_cf_mat")
    if (nrow(cf_mat) != length(y)) stop("cf_mat indicator matrix nrows different from # of obs.")
    if (length(y) != sum(cf_mat)) stop("cf_mat indicator matrix does not sum to number of observations.")
  }

  if (storeModels == "Disk" && is.null(path)) {
    stop("Provide a valid 'path' if storeModels is 'Disk' to save ensemble objects for later processing. Otherwise, set storeModels = 'No'.")
  }
  
  if (is.null(subset)) subset = rep(TRUE, length(y))

  np = rep(NA, length(y))
  models <- NULL
  

  ### Short-Stacking ###

  if (cv == 1) {
      
      ens <- ensemble_short(
        ml = ml, x = x, y = y, subset = subset,
        cf_mat = cf_mat, learner = learner,
        storeModels = storeModels, path = path, quiet = quiet, pb = pb, pb_np = pb_np
      )
      
      nnls_w <- nnls_weights(ens$fit_cv[subset, ], y[subset])
      np <- predict(ens, w = nnls_w)
      
      fit_sub <- list("ens_object" = ens, "nnls_w" = nnls_w)
      class(fit_sub) <- "ens.learner"
      
      if (storeModels == "Disk") {
        saveRDS(fit_sub, paste0(path, ".rds"))
      }
      
      if (storeModels == "Memory") {
        models <- fit_sub
      }
  

  ### Standard-Stacking ###

  } else if (cv > 1) {

    fit_sub = list()

    for (i in 1:ncol(cf_mat)) {

      fold = cf_mat[, i]
      x_tr = x[!fold & subset, ]
      y_tr = y[!fold & subset]
      x_te = x[fold, ]

      ens = ensemble(ml = ml, x = x_tr, y = y_tr, nfolds = cv, 
                     quiet = quiet, pb = pb, pb_np = pb_np, pb_cf = i
                     )
      
      ens_pred = predict(ens, ml, x = x_tr, y = y_tr, xnew = x_te, quiet = quiet, pb = pb, pb_np = pb_np, pb_cf = i)
      np[fold] = ens_pred$np

      fit_sub[[i]] = list("fit_cv" = ens_pred$fit_cv, "nnls_w" = ens$nnls_weights, "ens_object" = ens)

    }

    class(fit_sub) = "ens.learner"

    if (storeModels == "Disk") {
      saveRDS(fit_sub, paste0(path, ".rds"))
    }
    
    if (storeModels == "Memory") {
      models <- fit_sub
    }

  }
  
  output <- list(
    "np" = np,
    "models" = models
  )

  return(output)

}