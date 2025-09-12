#' Create a Design Matrix 
#' 
#' @description
#' \code{design_matrix} constructs an expanded design matrix by creating interaction terms,
#' polynomial expansions, and logarithmic transformations of specified variables.
#'
#' @param data A matrix or data.frame containing the predictor variables. Must have
#'   column names for all variables.
#' @param int Character vector of variable names to include in interaction terms.
#'   Use "all" to include all variables. If NULL (default), no interactions are created.
#' @param int_d Integer specifying the maximum degree of interactions to create.
#'   For example, `int_d = 2` creates all pairwise interactions. Default is 2.
#' @param poly Character vector of variable names for which polynomial terms should
#'   be created. Use "all" to include all variables. Polynomials are only created for degrees greater than 1.
#' @param poly_d Integer specifying the degree of polynomial expansion. For example,
#'   `poly_d = 3` creates terms up to cubic. Default is 2.
#' @param log Character vector of variable names to be log-transformed. Variables with
#'   non-positive values will be skipped with a warning.
#'
#' @return A numeric matrix with the expanded design matrix. The column names indicate
#'   the transformations applied:
#'   \itemize{
#'     \item Main effects retain their original names
#'     \item Interaction terms use colons (e.g., `var1:var2`)
#'     \item Polynomial terms use underscores (e.g., `var_2` for quadratic terms)
#'     \item Log-transformed variables are prefixed with "ln_" (e.g., `ln_var`)
#'   }
#'   The matrix includes an intercept term if present in the original data.
#' 
#' @examples
#' data(mtcars)
#' design_matrix(mtcars[, c("mpg", "wt", "hp")], 
#'               int = c("mpg", "wt"), 
#'               poly = "hp",
#'               log = "wt")
#'
#' @export
#'
design_matrix = function(data, int = NULL, int_d = 2, poly = NULL, poly_d = 2, log = NULL) {
  
  # Input validation
  if (!is.data.frame(data) && !is.matrix(data)) {
    stop("data must be a data.frame or matrix")
  }
  
  if (int_d < 1 || poly_d < 1) {
    stop("int_d and poly_d must be positive integers")
  }
  
  if (is.null(colnames(data)) || any(colnames(data) == "")) {
    stop("All columns in data must have names")
  }
  
  data_df <- as.data.frame(data)
  all_vars <- colnames(data_df)
  
  # Validate variable names
  validate_vars <- function(vars, name) {
    if (!is.null(vars) && length(vars) > 0 && !identical(vars, "all")) {
      invalid <- setdiff(vars, all_vars)
      if (length(invalid) > 0) {
        stop("Variables in ", name, " not found in data: ", paste(invalid, collapse = ", "))
      }
    }
  }
  
  validate_vars(int, "int")
  validate_vars(poly, "poly")
  validate_vars(log, "log")
  
  if (identical(int, "all")) {
    int <- all_vars
  }
  
  if (identical(poly, "all")) {
    poly <- all_vars
  }
  
  # Build formula components
  main_vars <- setdiff(all_vars, poly)  # Exclude variables that will have polynomials
  terms <- list()
  
  if (length(main_vars) > 0) {
    terms <- c(terms, paste(main_vars, collapse = " + "))
  }
  
  # Interactions
  if (!is.null(int) && length(int) > 0 && int_d > 1) {
    if (length(int) < int_d) {
      warning("Fewer interaction variables than interaction degree. Reducing int_d.")
      int_d <- min(length(int), int_d)
    }
    
    if (length(int) == 1) {
      # For single variable, just use main effect
      terms <- c(terms, int)
    } else {
      terms <- c(terms, paste0("(", paste(int, collapse = " + "), ")^", int_d))
    }
  }
  
  # Polynomials - only create if degree > 1
  if (!is.null(poly) && length(poly) > 0 && poly_d > 1) {
    poly_terms <- sapply(poly, function(var) {
      paste0("poly(", var, ", ", poly_d, ", raw = TRUE)")
    })
    terms <- c(terms, poly_terms)
  }
  
  # Logs
  if (!is.null(log) && length(log) > 0) {
    can_log <- sapply(log, function(var) {
      all(data_df[[var]] > 0, na.rm = TRUE)
    })
    
    if (any(!can_log)) {
      warning("Variables not logged due to non-positive values: ", 
              paste(log[!can_log], collapse = ", "))
      log <- log[can_log]
    }
    
    if (length(log) > 0) {
      log_terms <- paste0("log(", log, ")")
      terms <- c(terms, log_terms)
    }
  }
  
  # Build formula
  if (length(terms) == 0) {
    fmla <- stats::as.formula("~ 1")
  } else {
    fmla <- stats::as.formula(paste("~ 0 +", paste(terms, collapse = " + ")))
  }
  
  result <- stats::model.matrix(fmla, data = data_df)
  
  colnames(result) <- gsub("poly\\(([^,]+), \\d+, raw = TRUE\\)\\.(\\d+)", "\\1_\\2", colnames(result))
  colnames(result) <- gsub("poly\\(([^,]+), \\d+, raw = TRUE\\)", "\\1", colnames(result))
  colnames(result) <- gsub("log\\(([^)]+)\\)", "ln_\\1", colnames(result))
  colnames(result) <- make.names(colnames(result), unique = TRUE)
  
  return(result)
}


#' Data screening
#'
#' @description
#' Takes a matrix or data frame of numeric data and cleans it to improve its suitability for subsequent prediction algorithms.
#'
#' @details
#' This function performs the following steps in order:
#' \enumerate{
#'   \item Removes variables with zero or near-zero variance (standard deviation = 0 or NA).
#'   \item Removes dummy (binary) variables where the prevalence of one category is below a threshold
#'         in the overall sample or, if a treatment vector is provided, within any treatment group.
#'   \item Removes redundant variables that are highly correlated with others (keeping the first one encountered).
#' }
#'
#' @param data A matrix or data frame of numeric variables to be screened.
#' @param treat An optional treatment binary vector. If provided, binary variable 
#'              screening is performed separately within each treatment group.
#' @param bin_cut Numeric cut-off fraction (0 to 0.5) under which a category of a binary variable
#'                is considered "nearly empty". Default is 0.01.
#' @param corr_cut Numeric cut-off (0 to 1) above which the absolute correlation between two
#'                 variables is considered too high. Default is 0.99.
#' @param quiet Logical. If `FALSE`, details about the removed variables are printed at each step.
#'              Default is `TRUE`.
#'
#' @return A screened matrix with low-variance, near-constant binary, and redundant variables removed.
#'
#' @examples
#' set.seed(123)
#' n <- 100
#' data <- data.frame(
#'   constant = rep(1, n),
#'   good_predictor = rnorm(n),
#'   rare_dummy = sample(0:1, n, replace = TRUE, prob = c(0.99, 0.01)),
#'   correlated_var = good_predictor * 1.01 + rnorm(n, sd = 0.001)
#' )
#' screened_data <- data_screen(data, quiet = FALSE)
#' 
#' @export
#' 
data_screen <- function(data, treat = NULL, bin_cut = 0.01, corr_cut = 0.99, quiet = TRUE) {
  
  if (!is.null(treat)) {
    if (length(treat) != nrow(data)) {
      stop("Length of 'treat' must equal the number of rows in 'data'.")
    }
  }
  if (bin_cut < 0 || bin_cut > 0.5) {
    stop("'bin_cut' must be a value between 0 and 0.5.")
  }
  if (corr_cut < 0 || corr_cut > 1) {
    stop("'corr_cut' must be a value between 0 and 1.")
  }
  
  original_rownames <- rownames(data)
  orig_n_col <- ncol(data)
  
  if (!is.matrix(data)) data <- as.matrix(data)
  if (mode(data) != "numeric") {
    stop("'data' must contain only numeric values.")
  }
  
  # Eliminate variables with no (or near-zero) variation
  sds <- apply(data, 2, function(x) stats::sd(x, na.rm = TRUE))
  zero_var_mask <- sds == 0 | is.na(sds)
  nm_del_step1 <- colnames(data)[zero_var_mask]
  
  if (!quiet) {
    if (length(nm_del_step1) > 0) {
      cat("\nStep 1: Removed", length(nm_del_step1), "variables with no variation:",
          paste(nm_del_step1, collapse = ", "), "\n")
    } else {
      cat("\nStep 1: No variables removed for having no variation.\n")
    }
  }
  if (length(nm_del_step1) > 0) {
    data <- data[, !zero_var_mask, drop = FALSE]
  }
  
  # Early return if no variables left
  if (ncol(data) <= 1) {
    if (!quiet) cat("Screening complete. Fewer than 2 variables remain.\n")
    rownames(data) <- original_rownames
    return(data)
  }
  
  # Remove nearly-empty binary variables
  is_binary <- apply(data, 2, function(x) {
    all(na.omit(x) %in% 0:1)
  })
  
  nm_del_step2 <- character(0)
  if (any(is_binary)) {
    if (is.null(treat)) {
      # Check overall prevalence
      means <- Matrix::colMeans(data[, is_binary, drop = FALSE], na.rm = TRUE)
      below_cut <- (means < bin_cut | means > (1 - bin_cut))
    } else {
      # Check prevalence within each treatment group
      below_cut <- logical(sum(is_binary))
      names(below_cut) <- colnames(data)[is_binary]
      treat_groups <- unique(treat)
      
      for (col in colnames(data)[is_binary]) {
        for (grp in treat_groups) {
          group_data <- data[treat == grp, col]
          group_mean <- mean(group_data, na.rm = TRUE)
          if (group_mean < bin_cut || group_mean > (1 - bin_cut)) {
            below_cut[col] <- TRUE
            break # One failing group is enough to remove the variable
          }
        }
      }
    }
    
    nm_del_step2 <- names(below_cut)[below_cut]
    if (!quiet) {
      if (length(nm_del_step2) > 0) {
        cat("\nStep 2: Removed", length(nm_del_step2), "binary variables with a category prevalence <", bin_cut, ":",
            paste(nm_del_step2, collapse = ", "), "\n")
      } else {
        cat("\nStep 2: No binary variables removed for low prevalence.\n")
      }
    }
    if (length(nm_del_step2) > 0) {
      data <- data[, !colnames(data) %in% nm_del_step2, drop = FALSE]
    }
  } else if (!quiet) {
    cat("\nStep 2: No binary variables found to screen.\n")
  }
  
  # Early return if no variables left
  if (ncol(data) <= 1) {
    if (!quiet) cat("Screening complete. Fewer than 2 variables remain.\n")
    rownames(data) <- original_rownames
    return(data)
  }
  
  # Remove highly correlated variables
  cor_mat <- stats::cor(data, use = "pairwise.complete.obs")
  high_cor <- (abs(cor_mat) > corr_cut)
  high_cor[lower.tri(high_cor, diag = TRUE)] <- FALSE
  
  vars_to_remove <- character(0)
  if (any(high_cor)) {
    for (i in seq_len(ncol(high_cor) - 1)) {
      for (j in (i + 1):ncol(high_cor)) {
        if (high_cor[i, j] && !(colnames(data)[j] %in% vars_to_remove)) {
          vars_to_remove <- c(vars_to_remove, colnames(data)[j])
        }
      }
    }
  }
  
  nm_del_step3 <- unique(vars_to_remove)
  if (!quiet) {
    if (length(nm_del_step3) > 0) {
      cat("\nStep 3: Removed", length(nm_del_step3), "variables due to high correlation (r >", corr_cut, "):",
          paste(nm_del_step3, collapse = ", "), "\n")
    } else {
      cat("\nStep 3: No variables removed for high correlation.\n")
    }
  }
  if (length(nm_del_step3) > 0) {
    data <- data[, !colnames(data) %in% nm_del_step3, drop = FALSE]
  }
  
  if (!quiet) {
    cat("\nScreening complete. Removed", orig_n_col - ncol(data), "of", orig_n_col, "variables.\n")
  }
  
  rownames(data) <- original_rownames
  return(data)
}
