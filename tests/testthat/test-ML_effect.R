test_that("dml_inference works correctly", {
  # Test with simple inputs
  psi.a <- rep(1, 100)
  psi.b <- rep(2, 100)
  result <- dml_inference(psi.a, psi.b)
  
  # Check structure
  expect_named(result, c("TaPa", "IF", "score"))
  expect_equal(dim(result$TaPa), c(1, 4))
  expect_equal(dim(result$IF), c(100, 1))
  expect_equal(dim(result$score), c(1, 100, 3))
  
  # Check calculations
  expect_true(result$TaPa[1, "Estimate"] == -2)
  expect_true(result$TaPa[1, "SE"] == 0)
  expect_true(result$TaPa[1, "p"] == 0)
  
  # Test with different inputs
  psi.a <- 1:100
  psi.b <- 100:1
  result <- dml_inference(psi.a, psi.b)
  expect_true(result$TaPa[1, "Estimate"] == -mean(psi.b)/mean(psi.a))
})

# Test MLeffect function
test_that("MLeffect works correctly with basic inputs", {
  # Create test data
  set.seed(123)
  N <- 100
  Y <- rnorm(N)
  D <- sample(0:1, N, replace = TRUE)
  X <- matrix(rnorm(N*3), ncol = 3)
  
  # Create nuisance parameters
  NuPa.hat <- list(
    Y.hat = rnorm(N),
    Y.hat.d = cbind(rnorm(N), rnorm(N)),
    D.hat = stats::runif(N, 0.1, 0.9)
  )
  
  # Test PLR estimator
  result <- MLeffect(Y, D, X, NuPa.hat = NuPa.hat, estimators = "PLR")
  
  # Check class and structure
  expect_s3_class(result, "dml_with_smoother")
  expect_named(result, c("results", "NuPa.hat", "data", "numbers"))
  
  expect_true("PLR" %in% names(result$results))
  expect_type(result$results$PLR, "list")
  expect_named(result$results$PLR, c("TaPa", "IF", "score"))
  expect_equal(dim(result$results$PLR$TaPa), c(1, 4))
  
  # For estimators that weren't run
  not_run_estimators <- c("PLR_IV", "AIPW_ATE", "Wald_AIPW")
  for(est in not_run_estimators) {
    expect_identical(result$results[[est]], "This estimator was not run.")
  }

  expect_equal(result$numbers$N, N)
  expect_equal(result$numbers$n_estimators, 1)
})

test_that("MLeffect works with IV estimators when Z is provided", {
  set.seed(123)
  N <- 100
  Y <- rnorm(N)
  D <- sample(0:1, N, replace = TRUE)
  Z <- sample(0:1, N, replace = TRUE)
  X <- matrix(rnorm(N*3), ncol = 3)
  
  NuPa.hat <- list(
    Y.hat = rnorm(N),
    Y.hat.z = cbind(rnorm(N), rnorm(N)),
    D.hat = stats::runif(N, 0.1, 0.9),
    D.hat.z = cbind(rnorm(N), rnorm(N)),
    Z.hat = stats::runif(N, 0.1, 0.9)
  )
  
  # Test PLR_IV estimator
  result <- MLeffect(Y, D, X, Z, NuPa.hat = NuPa.hat, estimators = "PLR_IV")
  expect_named(result$results$PLR_IV, c("TaPa", "IF", "score"))
  expect_equal(dim(result$results$PLR_IV$TaPa), c(1, 4))
  
  # Test Wald_AIPW estimator
  result <- MLeffect(Y, D, X, Z, NuPa.hat = NuPa.hat, estimators = "Wald_AIPW")
  expect_named(result$results$Wald_AIPW, c("TaPa", "IF", "score"))
  expect_equal(dim(result$results$Wald_AIPW$TaPa), c(1, 4))
})

test_that("MLeffect handles errors correctly", {
  set.seed(123)
  N <- 100
  Y <- rnorm(N)
  D <- sample(0:1, N, replace = TRUE)
  X <- matrix(rnorm(N*3), ncol = 3)
  
  NuPa.hat <- list(
    Y.hat = rnorm(N),
    D.hat = stats::runif(N, 0.1, 0.9)
  )
  
  # Test unsupported estimator
  expect_error(MLeffect(Y, D, X, NuPa.hat = NuPa.hat, estimators = "UNSUPPORTED"),
               "The following specified estimators are not supported")
  
  # Test unsupported nuisance parameter
  bad_NuPa <- NuPa.hat
  bad_NuPa$unsupported <- rnorm(N)
  expect_error(MLeffect(Y, D, X, NuPa.hat = bad_NuPa),
               "The following passed nuisance parameters are not supported")
  
  # Test IV estimators without Z
  expect_error(MLeffect(Y, D, X, NuPa.hat = NuPa.hat, estimators = "PLR_IV"),
               "Z cannot be NULL when using either 'PLR_IV' or 'Wald_AIPW'")
  
  # Test non-binary D
  expect_error(MLeffect(Y, D = rep(1, N), X, NuPa.hat = NuPa.hat, estimators = "PLR"),
               "Treatment variable 'D' must be binary")
})

test_that("MLeffect works with multiple estimators", {
  set.seed(123)
  N <- 100
  Y <- rnorm(N)
  D <- sample(0:1, N, replace = TRUE)
  Z <- sample(0:1, N, replace = TRUE)
  X <- matrix(rnorm(N*3), ncol = 3)
  
  NuPa.hat <- list(
    Y.hat = rnorm(N),
    Y.hat.d = cbind(rnorm(N), rnorm(N)),
    Y.hat.z = cbind(rnorm(N), rnorm(N)),
    D.hat = stats::runif(N, 0.1, 0.9),
    D.hat.z = cbind(rnorm(N), rnorm(N)),
    Z.hat = stats::runif(N, 0.1, 0.9)
  )
  
  estimators <- c("PLR", "PLR_IV", "AIPW_ATE", "Wald_AIPW")
  result <- MLeffect(Y, D, X, Z, NuPa.hat = NuPa.hat, estimators = estimators)
  
  expect_length(result$results, 4)
  expect_named(result$results, estimators)
  expect_equal(result$numbers$n_estimators, 4)
  
  # Check all estimators produced results
  for(est in estimators) {
    expect_true(is.list(result$results[[est]]))
    expect_equal(dim(result$results[[est]]$TaPa), c(1, 4))
  }
})