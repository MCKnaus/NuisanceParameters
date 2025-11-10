#' Nuisance Parameter Estimation
#'
#' Estimates nuisance parameters using ensemble methods.
#'
#' @param NuPa Character vector specifying the nuisance parameters to estimate.
#'             Currently supported options:
#'             \code{c("Y.hat", "Y.hat.d", "Y.hat.z", "D.hat", "D.hat.z", "Z.hat")}
#' @param methods List of methods to use for \code{\link{ensemble}} estimation.
#'               Methods can be created using \code{\link{create_method}}.
#' @param Y A numeric vector indicating the outcome variable.
#' @param D A numeric vector indicating the treatment variable.
#' @param Z A numeric vector indicating the instrumental variable.
#' @param X A numeric matrix of covariates, with rows as observations and columns as covariates; no intercept included.
#' @param cf Number of cross-fitting folds
#' @param cf_mat Optional logical matrix of indicators representing the different
#'               cross-fitting folds, possibly from already estimated \code{NuisanceParameters} object.
#' @param cluster Optional vector of cluster variables if cross-fitting should
#'                account for clusters within the data.
#' @param stratify Logical. If \code{TRUE}, cross-fitting folds will preserve the
#'                 treatment ratios from the full sample. Note: If cluster vector
#'                 is provided, stratification is ignored due to computational
#'                 constraints and randomization feasibility issues.
#' @param stacking Stacking option: Either "short" for short stacking or an integer >1
#'                 indicating the number of cross-validation folds within each
#'                 cross-fitting fold.
#' @param store_models Character vector indicating whether to save individual
#'                    models for future processing (default: "no").
#'                    Supported options: \code{c("no", "memory", "disk")}
#' @param ensemble_type Method for calculating ensemble weights:
#'   \describe{
#'     \item{\code{"nnls"}}{Non-negative least squares; weights sum to 1 (default)}
#'     \item{\code{"bfgs"}}{BFGS optimization (for multivalued treatments only; falls to \code{nnls} otherwise)}
#'     \item{\code{"singlebest"}}{Weight of 1 on the learner with the lowest RMSE}
#'     \item{\code{"ols"}}{Ordinary least squares regression weights (falls to \code{nnls} for multivalued treatments)}
#'     \item{\code{"average"}}{Equal weights for all learners}
#'   }
#' @param path Optional path to save the \code{\link{ensemble}} objects for
#'             later processing (saved as a list of models).
#' @param quiet Logical. If \code{FALSE}, progress output is printed to the console.
#'
#' @return A list containing:
#' \itemize{
#'   \item \code{nuisance_parameters}: Requested nuisance parameter estimates
#'   \item \code{models} (optional): \code{\link{ensemble}} objects (individual models)
#'   \item \code{numbers}: Additional objects used for downstream processing in
#'                         smoother matrices and outcome weights extraction
#' }
#' 
#' @examples
#' \donttest{
#'   if (requireNamespace("hdm", quietly = TRUE)) {
#'     data(pension, package = "hdm")
#'     
#'     set.seed(123)
#'     D <- pension$p401
#'     Y <- pension$net_tfa
#'     X <- model.matrix(~ 0 + age + db + educ + fsize + hown + inc + male + 
#'                       marr + pira + twoearn, data = pension)
#'     
#'     methods = list(
#'       "ols" = create_method("ols"),
#'       "forest_grf" = create_method("forest_grf"),
#'       "xgboost" = create_method("xgboost")
#'     )
#'     
#'     np <- nuisance_parameters(
#'       NuPa = c("Y.hat", "D.hat"),
#'       X = X, Y = Y, D = D,
#'       methods = methods, cf = 2, stacking = "short"
#'     )
#'     
#'     print(head(np$nuisance_parameters[["Y.hat"]]))
#'     plot(np$numbers$ens_weights)
#'   }
#' }
#'
#' @export
nuisance_parameters <- function(NuPa = c("Y.hat", "Y.hat.d", "Y.hat.z", "D.hat", "D.hat.z", "Z.hat"),
                                X, 
                                Y = NULL, 
                                D = NULL, 
                                Z = NULL,
                                methods,
                                stacking,
                                cluster = NULL,
                                stratify = FALSE,
                                cf,
                                cf_mat = NULL,
                                ensemble_type = "nnls",
                                store_models = c("no", "memory", "disk"),
                                path = NULL,
                                quiet = TRUE
                                ) {
  ## Checks
  store_models <- match.arg(store_models)
  ensemble_type <- match.arg(ensemble_type, choices = c('nnls', 'bfgs', 'singlebest', 'ols', 'average'))
  
  if (cf < 1) stop("Number of cross-fitting folds must be a positive integer.")
  if (identical(stacking, "short") && cf == 1) stop("Short-stacking with cf = 1 is not implemented.")
  if (is.null(path)) path <- getwd()

  supported_NuPas <- c("Y.hat", "Y.hat.d", "Y.hat.z", "D.hat", "D.hat.z", "Z.hat")
  not_supported <- NuPa[!NuPa %in% supported_NuPas]

  if (length(not_supported) > 0) {
    stop(paste("Error: The following nuisance parameters specified in NuPa are not supported:",
               paste(not_supported, collapse = ", ")))}
  if (is.null(Y) & any(NuPa %in% c("Y.hat", "Y.hat.d", "Y.hat.z"))) {
    stop("Please specify Y if at least one of c(\"Y.hat\",\"Y.hat.d\",\"Y.hat.z\") is specified in NuPa")}
  if (is.null(D) & any(NuPa %in% c("Y.hat.d", "D.hat", "D.hat.z"))) {
    stop("Please specify D if at least one of c(\"Y.hat.d\",\"D.hat\",\"D.hat.z\") is specified in NuPa")}
  if (is.null(Z) & any(NuPa %in% c("Y.hat.z", "D.hat.z", "Z.hat"))) {
    stop("Please specify Z if at least one of c(\"Y.hat.z\",\"D.hat.z\",\"Z.hat\") is specified in NuPa")}

  if (stacking == "short") {
    cv <- 1L
  } else if (is.numeric(stacking) && stacking > 1 && stacking == round(stacking)) {
    cv <- as.integer(stacking)
  } else {
    stop("`stacking` must be either 'short' or an integer > 1")
  }

  if (length(methods) == 1) {
    cv <- 1
    message("With a single learner provided, estimation defaults to cross-fitted predictions.")
  }

  
  ## Preps
  N <- nrow(X)
  K <- length(unique(D))
  M <- length(unique(Y))
  d_mat <- if (!is.null(D)) prep_indicator_mat(D) else NULL
  z_mat <- if (!is.null(Z)) prep_indicator_mat(Z) else NULL
  if (is.null(cf_mat)) {
    cf_mat <- prep_cf_matrix(N, cf = cf, d_mat = if (stratify) d_mat else NULL, cluster = cluster)
  }

  if (isFALSE(quiet)) which_stacking(cv)
  ens_w <- list()

  
  ## Initialize objects
  Y.hat <- Y.hat.d <- Y.hat.z <- D.hat <- D.hat.z <- Z.hat <-
    "This nuisance parameter was not specified and is therefore empty."

  Y.hat_m <- Y.hat.d_m <- Y.hat.z_m <- D.hat_m <- D.hat.z_m <- Z.hat_m <-
    "This model was not specified and is therefore empty."

  if ("Y.hat" %in% NuPa) { 
    Y.hat <- rep(NA, NROW(Y))
    Y.hat_m <- list()
  }
  if ("Y.hat.d" %in% NuPa) {
    Y.hat.d <- matrix(NA, nrow(d_mat), ncol(d_mat))
    Y.hat.d_m <- vector("list", ncol(d_mat))
  }
  if ("Y.hat.z" %in% NuPa) {
    Y.hat.z <- matrix(NA, nrow(z_mat), ncol(z_mat))
    Y.hat.z_m <- vector("list", ncol(z_mat))
  }
  if ("D.hat" %in% NuPa) {
    D.hat <- matrix(NA, nrow(d_mat), ncol(d_mat) - 1)
    D.hat_m <- list()
  }
  if ("D.hat.z" %in% NuPa) {
    D.hat.z <- if (K > 2) { array(NA, dim = c(nrow(z_mat), K, ncol(z_mat)))
    } else { matrix(NA, nrow(z_mat), ncol(z_mat), dimnames = dimnames(z_mat)) }
    D.hat.z_m <- vector("list", ncol(z_mat))
  }
  if ("Z.hat" %in% NuPa) {
    Z.hat <- matrix(NA, nrow(z_mat), ncol(z_mat) - 1)
    Z.hat_m <- list()
  }

  if ("Y.hat.d" %in% NuPa) colnames(Y.hat.d) <- colnames(d_mat)
  if ("Y.hat.z" %in% NuPa) colnames(Y.hat.z) <- colnames(z_mat)
  if ("D.hat.z" %in% NuPa) {
    if (K > 2) { dimnames(D.hat.z) <- list(NULL, paste0("D", seq(K) - 1), colnames(z_mat))
    } else { colnames(D.hat.z) <- colnames(z_mat) }
  }

  methods <- process_methods(methods = methods, NuPa = NuPa, K = K, M = M)

  # Create progress bar
  pb <- if (isFALSE(quiet)) {
    setup_pb(NuPa = NuPa, n_d = ncol(d_mat), n_z = ncol(z_mat), cf_folds = cf, cv_folds = cv, methods = methods)
  } else { NULL }


  ## Outcome NuPa
  if ("Y.hat.d" %in% NuPa) {
    nupa <- "Y.hat.d"

    for (i in seq_len(ncol(d_mat))) {
      pb_np <- paste0(nupa, i - 1)

      np_cf <- nuisance_cf(
        methods = methods[[nupa]], Y = Y, X = X, cf_mat = cf_mat, 
        cv = cv, subset = d_mat[, i], store_models = store_models, 
        quiet = quiet, pb = pb, pb_np = pb_np, ensemble_type = ensemble_type
      )

      Y.hat.d[, i] <- np_cf$np
      Y.hat.d_m[[i]] <- np_cf$models_nupa
      ens_w[[pb_np]] <- extract_w(np_cf)
    }
  }

  if ("Y.hat.z" %in% NuPa) {
    nupa <- "Y.hat.z"

    for (i in seq_len(ncol(z_mat))) {
      pb_np <- paste0(nupa, i - 1)

      np_cf <- nuisance_cf(
        methods = methods[[nupa]], Y = Y, X = X, cf_mat = cf_mat, 
        cv = cv, subset = z_mat[, i], store_models = store_models, 
        quiet = quiet, pb = pb, pb_np = pb_np, ensemble_type = ensemble_type
      )

      Y.hat.z[, i] <- np_cf$np
      Y.hat.z_m[[i]] <- np_cf$models_nupa
      ens_w[[pb_np]] <- extract_w(np_cf)
    }
  }

  if ("Y.hat" %in% NuPa) {
    nupa <- pb_np <- "Y.hat"

    np_cf <- nuisance_cf(
      methods = methods[[nupa]], Y = Y, X = X, cf_mat = cf_mat, 
      cv = cv, subset = NULL, store_models = store_models, 
      quiet = quiet, pb = pb, pb_np = pb_np, ensemble_type = ensemble_type 
    )

    Y.hat <- np_cf$np
    Y.hat_m <- np_cf$models_nupa
    ens_w[[pb_np]] <- extract_w(np_cf)
  }


  ## Treatment NuPa
  if ("D.hat.z" %in% NuPa) {
    nupa <- "D.hat.z"

    for (i in seq_len(ncol(z_mat))) {
      pb_np <- paste0(nupa, i - 1)

      np_cf <- nuisance_cf(
        methods = methods[[nupa]], Y = D, X = X, cf_mat = cf_mat, 
        cv = cv, subset = z_mat[, i], store_models = store_models, 
        quiet = quiet, pb = pb, pb_np = pb_np, ensemble_type = ensemble_type
      )

      if (K > 2) { D.hat.z[, , i] <- np_cf$np } else { D.hat.z[, i] <- np_cf$np }
      D.hat.z_m[[i]] <- np_cf$models_nupa
      ens_w[[pb_np]] <- extract_w(np_cf)
    }
  }

  if ("D.hat" %in% NuPa) {
    nupa <- pb_np <- "D.hat"

    np_cf <- nuisance_cf(
      methods = methods[[nupa]], Y = D, X = X, cf_mat = cf_mat, 
      cv = cv, subset = NULL, store_models = store_models, 
      quiet = quiet, pb = pb, pb_np = pb_np, ensemble_type = ensemble_type
    )

    D.hat <- np_cf$np
    D.hat_m <- np_cf$models_nupa
    ens_w[[pb_np]] <- extract_w(np_cf)
  }


  ## Instrument NuPa
  if ("Z.hat" %in% NuPa) {
    nupa <- pb_np <- "Z.hat"

    np_cf <- nuisance_cf(
      methods = methods[[nupa]], Y = Z, X = X, cf_mat = cf_mat, 
      cv = cv, subset = NULL, store_models = store_models, 
      quiet = quiet, pb = pb, pb_np = pb_np, ensemble_type = ensemble_type
    )

    Z.hat <- np_cf$np
    Z.hat_m <- np_cf$models_nupa
    ens_w[[pb_np]] <- extract_w(np_cf)
  }

  
  ## Structure the output
  ens_weights_mat <- format_weights(ens_w, cv = cv)

  np_list <- list(
    "Y.hat" = Y.hat, "Y.hat.d" = Y.hat.d, "Y.hat.z" = Y.hat.z,
    "D.hat" = D.hat, "D.hat.z" = D.hat.z, "Z.hat" = Z.hat
  )

  models_list <- list(
    "Y.hat_m" = Y.hat_m, "Y.hat.d_m" = Y.hat.d_m, "Y.hat.z_m" = Y.hat.z_m,
    "D.hat.z_m" = D.hat.z_m, "D.hat_m" = D.hat_m, "Z.hat_m" = Z.hat_m)

  nums_list <- list(
    "N" = N, "cv" = cv, "cf_mat" = cf_mat, "d_mat" = d_mat, "z_mat" = z_mat,
    "methods" = methods, "X" = X, "Y" = Y, "ens_weights" = ens_weights_mat
  )
  
  if (store_models == "disk") {
    models_path <- file.path(path, "nuisance_models.rds")
    saved_object <- list(models = models_list, numbers = nums_list)
    class(saved_object) <- c("NuisanceParameters")
    saveRDS(saved_object, file = models_path)
    if (!quiet) message("Models saved to: ", models_path)
  }

  output <- list(
    nuisance_parameters = np_list,
    models = if (store_models == "memory") models_list else NULL,
    numbers = nums_list
  )

  class(output) <- c("NuisanceParameters")
  return(output)
}


#' Cross-fitting of nuisance parameters
#'
#' @description
#' \code{\link{nuisance_cf}} makes a cross-fitted ensemble prediction using
#' \code{\link{ensemble}}.
#'
#' @param methods List of methods to be used in \code{\link{ensemble}} estimation.
#' @param Y Vector of variable to be predicted.
#' @param X Matrix of covariates.
#' @param cf_mat Logical matrix with k columns of indicators representing the different folds.
#' @param cv Number of cross-validation folds when estimating ensemble model (default: 1).
#' @param subset Optional logical vector if only subset of data should be used for prediction.
#' @param store_models Vector of characters indicating where to save individual models for future processing.
#' @param ensemble_type Method for calculating ensemble weights:
#'   \describe{
#'     \item{\code{"nnls"}}{Non-negative least squares; weights sum to 1 (default)}
#'     \item{\code{"bfgs"}}{BFGS optimization (for multivalued treatments only; falls back to \code{nnls} otherwise)}
#'     \item{\code{"singlebest"}}{Weight of 1 on the learner with the lowest RMSE}
#'     \item{\code{"ols"}}{Ordinary least squares regression weights (falls back to \code{nnls} for multivalued treatments)}
#'     \item{\code{"average"}}{Equal weights for all learners}
#'   }
#' @param quiet If FALSE, progress output is printed into console.
#' @param pb A progress bar object to track overall computation progress.
#' @param pb_np String indicating the current nuisance parameter being
#'              estimated/predicted (for progress bar updates).
#'
#' @return A list containing:
#' \itemize{
#'   \item \code{np}: A vector (or matrix for multinomial outcomes) of length N 
#'                    containing nuisance parameter estimates
#'   \item \code{models}: An object of class \code{ens.learner} containing fitted ensemble models
#'  }
#'
#' @keywords internal
nuisance_cf <- function(methods, 
                        Y, X, 
                        cf_mat,
                        cv = 5,
                        subset = NULL,
                        store_models = c("no", "memory", "disk"),
                        ensemble_type, 
                        quiet = TRUE,
                        pb = NULL, 
                        pb_np = NULL) {
  ## Checks
  if (is.null(subset)) subset <- rep(TRUE, length(Y))

  # Is multiclass propensity score being estimated?
  is_mult <- !is.null(methods[[1]]$multinomial)
  
  if (is_mult) { np <- matrix(NA, nrow = length(Y), ncol = length(unique(Y))) 
  } else { np <- rep(NA, length(Y)) }
  
  # Full-sample hyperparameter tuning
  mtd_tuned <- tune_learners(type = "full_sample", X = X[subset, ], Y = Y[subset], methods = methods)
  
  ### Short-Stacking ###
  if (cv == 1) {
    ens <- ensemble_short(
      methods = mtd_tuned, X = X, Y = Y, subset = subset, cf_mat = cf_mat, 
      store_models = store_models, quiet = quiet, pb = pb, pb_np = pb_np)
    
    ens_w <- ens_weights_maker(
      X = ens$cf_preds, Y = Y, subset = subset, is_mult = is_mult, ensemble_type = ensemble_type
      )
    np <- if (is_mult) agg_array(a = ens$cf_preds, w = ens_w) else predict(ens, w = ens_w)
    
    models_nupa <- list("ens_object" = ens, "ens_w" = ens_w)
    class(models_nupa) <- "ens.learner"

    ### Standard-Stacking ###
  } else if (cv > 1) {
    models_nupa <- list()

    for (i in 1:ncol(cf_mat)) {
      fold  <- cf_mat[, i]
      no_cf <- ncol(cf_mat) == 1
      
      X_tr <- X[if (no_cf) subset else (!fold & subset), ]
      Y_tr <- Y[if (no_cf) subset else (!fold & subset)]
      X_te <- X[if (no_cf) TRUE else fold, ]

      # On-the-fold hyperparameter tuning
      mtd_tuned_fold <- tune_learners(type = "fold", X = X_tr, Y = Y_tr, methods = mtd_tuned)

      ens <- ensemble(
        methods = mtd_tuned_fold, X = X_tr, Y = Y_tr, nfolds = cv,
        ensemble_type = ensemble_type, quiet = quiet, pb = pb, pb_np = pb_np, pb_cf = i)
      
      ens_pred <- predict(
        object = ens, methods = mtd_tuned_fold, X = X_tr, Y = Y_tr, Xnew = X_te, 
        quiet = quiet, pb = pb, pb_np = pb_np, pb_cf = i)
        np[fold] <- ens_pred$np
        
      ens_object <- list("ens_models" = ens$ens_models, "cf_preds" = ens_pred$cf_preds)
      class(ens_object) <- "ensemble"
        
      models_nupa[[i]] <- list(
        "ens_object" = ens_object, 
        "ens_w" = ens$ens_w
        )
    }
    class(models_nupa) <- "ens.learner"
  }
  output <- list(
    "np" = np,
    "models_nupa" = models_nupa
  )

  return(output)
}
