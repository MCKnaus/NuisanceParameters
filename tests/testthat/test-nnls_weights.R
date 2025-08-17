test_that("nnls_weights produces correct weights", {
  set.seed(42)
  n = 1000
  p = 5
  X = matrix(runif(n*p), ncol = p)
  true_weights = c(0.3, 0, 0, 0.1, 0.6)
  Y = X %*% true_weights + rnorm(n, sd = 0.1)

  estimated_weights = nnls_weights(X, Y)

  expect_equal(estimated_weights, true_weights, tolerance = 0.05)
})

test_that("nnls_weights produces correct weights when features are highly correlated", {
  library(mvtnorm)
  set.seed(42)
  n = 1000
  p = 8
  cov_mat = matrix(rep(0.9, p*p), ncol = p, nrow = p)
  cov_mat[cbind(1:p, 1:p)] = 1
  X = mvtnorm::rmvnorm(n = n, mean = rep(0, p), sigma = cov_mat)
  true_weights = c(0.05, 0.05, 0.1, 0, 0.3, 0.3, 0.2, 0)
  Y = X %*% true_weights + rnorm(n, sd = 0.1)

  estimated_weights = nnls_weights(X, Y)

  expect_equal(estimated_weights, true_weights, tolerance = 0.05)
})
