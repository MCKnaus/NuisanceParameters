test_that("check prediction and weights method of ensemble_core", {

  library(mvtnorm)

  n = 2000
  n_test = 25
  p = 12
  p_act = 4

  pi = c(seq(1, 0.1, -(1 / p_act)) * rep(c(1, -1), p_act / 2), rep(0, p - p_act))
  cov_mat = toeplitz(0.7^(0:(p - 1)))

  x_tr = mvtnorm::rmvnorm(n = n, mean = rep(0, p), sigma = cov_mat)
  x_te = mvtnorm::rmvnorm(n = n_test, mean = rep(0, p), sigma = cov_mat)
  y_tr = as.vector(x_tr %*% pi + rnorm(n, 0, 1))

  # first ensemble
  ml_a = list("ols" = create_method("ols"),
              "plasso" = create_method("plasso"),
              "mean" = create_method("mean"))

  ens_core_a = ensemble_core(ml_a, x_tr = x_tr, y_tr = y_tr)
  pred_a = predict(ens_core_a, ml_a, x_tr, y_tr, x_te)
  w_a = weights(ens_core_a, ml_a, x_tr, y_tr, x_te)
  expect_equal(dim(w_a), c(nrow(x_te), nrow(x_tr), length(ml_a)))

  pred_s_a = sapply(1:length(ml_a), function(i) w_a[, , i] %*% y_tr)
  expect_equal(pred_a, pred_s_a, tolerance = 1e-9)

  w_sums_a = sapply(1:length(ml_a), function(i) rowSums(w_a[, , i]))
  expect_equal(matrix(1, nrow = nrow(x_te), ncol = length(ml_a)), w_sums_a, tolerance = 1e-9)


  # second ensemble
  ml_b = list("ols" = create_method("ols"),
              "plasso" = create_method("plasso"),
              "ridge" = create_method("ridge"),
              "forest_grf" = create_method("forest_grf"),
              "knn_10" = create_method("knn", arguments = list("k" = 10)))

  ens_core_b = ensemble_core(ml_b, x_tr = x_tr, y_tr = y_tr)
  pred_b = predict(ens_core_b, ml_b, x_tr, y_tr, x_te)
  w_b = weights(ens_core_b, ml_b, x_tr, y_tr, x_te)
  expect_equal(dim(w_b), c(nrow(x_te), nrow(x_tr), length(ml_b)))

  pred_s_b = sapply(1:length(ml_b), function(i) w_b[, , i] %*% y_tr)
  expect_equal(pred_b, pred_s_b, tolerance = 1e-3)

  w_sums_b = sapply(1:length(ml_b), function(i) rowSums(w_b[, , i]))
  expect_equal(matrix(1, nrow = nrow(x_te), ncol = length(ml_b)), w_sums_b, tolerance = 1e-9)

})
