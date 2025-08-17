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
  Y = as.vector(X %*% pi + rnorm(n, 0, 1))

  w_rows = rep(1, n_test)

  # mean
  m = mean_fit(X = X, Y = Y)
  p = NuisanceParameters:::predict.mean_fit(m, Xnew = Xnew)
  w = NuisanceParameters:::weights.mean_fit(m, Xnew = Xnew)
  p_s = as.vector(w %*% Y)
  expect_equal(Matrix::rowSums(w), w_rows, tolerance = 1e-9)
  expect_true(all(w >= 0 & w <= 1))
  expect_equal(p, p_s, tolerance = 1e-9)
  
  # ols
  m = ols_fit(X = X, Y = Y)
  p = NuisanceParameters:::predict.ols_fit(m, Xnew = Xnew)
  w = NuisanceParameters:::weights.ols_fit(m, X = X, Xnew = Xnew)
  p_s = as.vector(w %*% Y)
  expect_equal(Matrix::rowSums(w), w_rows, tolerance = 1e-9)
  expect_equal(p, p_s, tolerance = 1e-9)
  
  # ridge
  m = ridge_fit(X = X, Y = Y)
  p = NuisanceParameters:::predict.ridge_fit(m, Xnew = Xnew)
  w = NuisanceParameters:::weights.ridge_fit(m, X = X, Y = Y, Xnew = Xnew)
  p_s = as.vector(w %*% Y)
  expect_equal(Matrix::rowSums(w), w_rows, tolerance = 1e-9)
  expect_equal(p, p_s, tolerance = 0.001)
  
  # plasso
  m = plasso_fit(X = X, Y = Y)
  p = NuisanceParameters:::predict.plasso_fit(m, Xnew = Xnew)
  w = NuisanceParameters:::weights.plasso_fit(m, Xnew = Xnew)
  p_s = as.vector(w %*% Y)
  expect_equal(Matrix::rowSums(w), w_rows, tolerance = 1e-9)
  expect_equal(p, p_s, tolerance = 1e-9)
  
  # rlasso
  m = rlasso_fit(X = X, Y = Y)
  p = NuisanceParameters:::predict.rlasso_fit(m, Xnew = Xnew)
  w = NuisanceParameters:::weights.rlasso_fit(m, Xnew = Xnew)
  p_s = as.vector(w %*% Y)
  expect_equal(Matrix::rowSums(w), w_rows, tolerance = 1e-9)
  expect_equal(p, p_s, tolerance = 1e-9)
  
  # random forest
  m = forest_grf_fit(X = X, Y = Y)
  p = NuisanceParameters:::predict.forest_grf_fit(m, Xnew = Xnew)
  w = NuisanceParameters:::weights.forest_grf_fit(m, Xnew = Xnew)
  p_s = as.vector(w %*% Y)
  expect_equal(Matrix::rowSums(w), w_rows, tolerance = 1e-9)
  expect_true(all(w >= 0 & w <= 1))
  expect_equal(p, p_s, tolerance = 1e-9)
  
  # knn
  m = knn_fit(arguments = list("k" = 5))
  p = NuisanceParameters:::predict.knn_fit(m, X = X, Y = Y, Xnew = Xnew)
  w = NuisanceParameters:::weights.knn_fit(m, X = X, Xnew = Xnew)
  p_s = as.vector(w %*% Y)
  expect_equal(Matrix::rowSums(w), w_rows, tolerance = 1e-9)
  expect_true(all(w >= 0 & w <= 1))
  expect_equal(p, p_s, tolerance = 1e-9)
  
  # distributional random forest
  m = forest_drf_fit(X = X, Y = Y)
  p = NuisanceParameters:::predict.forest_drf_fit(m, Xnew = Xnew)
  w = NuisanceParameters:::weights.forest_drf_fit(m, Xnew = Xnew)
  p_s = as.vector(w %*% Y)
  expect_equal(Matrix::rowSums(w), w_rows, tolerance = 1e-9)
  expect_true(all(w >= 0 & w <= 1))
  expect_equal(p, p_s, tolerance = 1e-9)
  
  ## XGBoost: does not achieve perfect equalty
  # # xgboost
  # m = xgboost_fit(X = X, Y = Y)
  # p = NuisanceParameters:::predict.xgboost_fit(m, Xnew = Xnew)
  # w = NuisanceParameters:::weights.xgboost_fit(m, X = X, Xnew = Xnew)
  # p_s = as.vector(w %*% Y)
  # expect_equal(Matrix::rowSums(w), w_rows, tolerance = 1e-7)
  # expect_true(all(w >= 0 & w <= 1))
  # expect_equal(p, p_s, tolerance = 1e-7)
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
  Y = c(1, 3, 1, 3, 1, 9, 3)

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
  w = NuisanceParameters:::weights.knn_fit(m, X = X, Xnew = Xnew)

  expect_identical(w_exp, w)

  p = NuisanceParameters:::predict.knn_fit(m, X = X, Y = Y, Xnew = Xnew)
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
  Y = X %*% pi + rnorm(n, 0, 1)

  k = 9
  m = knn_fit(arguments = list("k" = k))
  p = NuisanceParameters:::predict.knn_fit(m, X = X, Y = Y, Xnew = Xnew)
  expect_true(all(p >= min(Y) & p <= max(Y)))

  w = NuisanceParameters:::weights.knn_fit(m, X = X, Xnew = Xnew)
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
  Y = X %*% pi + rnorm(n, 0, 0.5)


  m_ols = ols_fit(X = X, Y = Y)
  p_ols = NuisanceParameters:::predict.ols_fit(m_ols, Xnew = Xnew)

  m_p = plasso_fit(X = X, Y = Y)
  p_p = NuisanceParameters:::predict.plasso_fit(m_p, Xnew = Xnew)

  expect_equal(p_ols, p_p, tolerance = 0.01)

})
