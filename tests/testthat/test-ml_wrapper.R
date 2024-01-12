test_that("smoother weights lead to correct prediction of fitted values", {

  library(mvtnorm)

  n = 200
  n_test = 50
  p = 12
  p_act = 4

  pi = c(seq(1, 0.1, -(1 / p_act)) * rep(c(1, -1), p_act / 2), rep(0, p - p_act))
  cov_mat = toeplitz(0.7^(0:(p - 1)))

  X = mvtnorm::rmvnorm(n = n, mean = rep(0, p), sigma = cov_mat)
  Xnew = mvtnorm::rmvnorm(n = n_test, mean = rep(0, p), sigma = cov_mat)
  y = X %*% pi + rnorm(n, 0, 1)

  # mean
  m = mean_fit(x = X, y = y)
  p = predict(m, xnew = Xnew)
  w = weights(m, xnew = Xnew)
  p_s = as.vector(w %*% y)
  expect_equal(p, p_s, tolerance = 1e-9)

  # ols
  m = ols_fit(x = X, y = y)
  p = predict(m, xnew = Xnew)
  w = weights(m, x = X, xnew = Xnew)
  p_s = as.vector(w %*% y)
  expect_equal(p, p_s, tolerance = 1e-9)

  # ridge
  m = ridge_fit(x = X, y = y)
  p = predict(m, xnew = Xnew, y = y, exact = TRUE)
  w = weights(m, x = X, y = y, xnew = Xnew)
  p_s = as.vector(w %*% y)
  expect_equal(p, p_s, tolerance = 0.001)

})
