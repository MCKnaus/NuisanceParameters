test_that("smoother weights lead to correct prediction of fitted values", {

  library(mvtnorm)

  n = 2000
  n_test = 25
  p = 12
  p_act = 4

  pi = c(seq(1, 0.1, -(1 / p_act)) * rep(c(1, -1), p_act / 2), rep(0, p - p_act))
  cov_mat = toeplitz(0.7^(0:(p - 1)))

  X = mvtnorm::rmvnorm(n = n, mean = rep(0, p), sigma = cov_mat)
  Xnew = mvtnorm::rmvnorm(n = n_test, mean = rep(0, p), sigma = cov_mat)
  y = X %*% pi + rnorm(n, 0, 1)

  w_rows = rep(1, n_test)

  # mean
  m = mean_fit(x = X, y = y)
  p = predict(m, xnew = Xnew)
  w = weights(m, xnew = Xnew)
  p_s = as.vector(w %*% y)
  expect_equal(Matrix::rowSums(w), w_rows, tolerance = 1e-9)
  expect_true(all(w >= 0 & w <= 1))
  expect_equal(p, p_s, tolerance = 1e-9)

  # ols
  m = ols_fit(x = X, y = y)
  p = predict(m, xnew = Xnew)
  w = weights(m, x = X, xnew = Xnew)
  p_s = as.vector(w %*% y)
  expect_equal(Matrix::rowSums(w), w_rows, tolerance = 1e-9)
  expect_equal(p, p_s, tolerance = 1e-9)

  # ridge
  m = ridge_fit(x = X, y = y)
  p = predict(m, xnew = Xnew)
  w = weights(m, x = X, y = y, xnew = Xnew)
  p_s = as.vector(w %*% y)
  expect_equal(Matrix::rowSums(w), w_rows, tolerance = 1e-9)
  expect_equal(p, p_s, tolerance = 0.001)

  # plasso
  m = plasso_fit(x = X, y = y)
  p = predict(m, xnew = Xnew)
  w = weights(m, xnew = Xnew)
  p_s = as.vector(w %*% y)
  expect_equal(Matrix::rowSums(w), w_rows, tolerance = 1e-9)
  expect_equal(p, p_s, tolerance = 1e-9)

  # random forest
  m = forest_grf_fit(x = X, y = y)
  p = predict(m, xnew = Xnew)
  w = weights(m, xnew = Xnew)
  p_s = as.vector(w %*% y)
  expect_equal(Matrix::rowSums(w), w_rows, tolerance = 1e-9)
  expect_true(all(w >= 0 & w <= 1))
  expect_equal(p, p_s, tolerance = 1e-9)

  # knn
  m = knn_fit(arguments = list("k" = 5))
  p = predict(m, x = X, y = y, xnew = Xnew)
  w = weights(m, x = X, xnew = Xnew)
  p_s = as.vector(w %*% y)
  expect_equal(Matrix::rowSums(w), w_rows, tolerance = 1e-9)
  expect_true(all(w >= 0 & w <= 1))
  expect_equal(p, p_s, tolerance = 1e-9)

  # distributional random forest
  m = forest_drf_fit(x = X, y = y)
  p = predict(m, xnew = Xnew)
  w = weights(m, xnew = Xnew)
  p_s = as.vector(w %*% y)
  expect_equal(Matrix::rowSums(w), w_rows, tolerance = 1e-9)
  expect_true(all(w >= 0 & w <= 1))
  expect_equal(p, p_s, tolerance = 1e-9)

})


test_that("knn mickey mouse check", {

  X = matrix(c(
    1, 2, 3,
    7, 7, 5,
    1, 3, 1,
    6, 8, 6,
    2, 2, 2,
    3, 4, 6,
    4, 7, 9
  ), ncol = 3, byrow = TRUE)
  y = c(1, 3, 1, 3, 1, 9, 3)

  Xnew = matrix(c(
    1, 3, 2,
    5, 7, 8
  ), ncol = 3, byrow = TRUE)

  k = 3
  w_exp = matrix(c(
    1, 0, 1, 0, 1, 0, 0,
    0, 1, 0, 1, 0, 0, 1),
    nrow = 2, byrow = TRUE
  ) / k

  m = knn_fit(arguments = list("k" = k))
  w = weights(m, x = X, xnew = Xnew)

  expect_identical(w_exp, w)

  p = predict(m, x = X, y = y, xnew = Xnew)
  expect_identical(p, c(1, 3))

})


test_that("knn weight matrix plausibility check", {

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

  k = 9
  m = knn_fit(arguments = list("k" = k))
  p = predict(m, x = X, y = y, xnew = Xnew)
  expect_true(all(p >= min(y) & p <= max(y)))

  w = weights(m, x = X, xnew = Xnew)
  w_vals = unique(as.vector(as.matrix(w[-1])))
  expect_identical(sort(w_vals), c(0, 1/k))

})


test_that("(rough) equality of plasso and ols prediction solution for clear cut dgp", {

  library(mvtnorm)

  n = 100000
  n_test = 10
  p = 6

  pi = c(3, -3, 3, 0, 0, 0)
  cov_mat = matrix(0.1, ncol = p, nrow = p)
  diag(cov_mat) = 1

  X = mvtnorm::rmvnorm(n = n, mean = rep(0, p), sigma = cov_mat)
  Xnew = mvtnorm::rmvnorm(n = n_test, mean = rep(0, p), sigma = cov_mat)
  y = X %*% pi + rnorm(n, 0, 0.5)


  m_ols = ols_fit(x = X, y = y)
  p_ols = predict(m_ols, xnew = Xnew)

  m_p = plasso_fit(x = X, y = y)
  p_p = predict(m_p, xnew = Xnew)

  expect_equal(p_ols, p_p, tolerance = 0.01)

})
