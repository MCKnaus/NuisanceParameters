test_that("check prediction and weights method of ensemble_core", {

  skip_if_not_installed("OutcomeWeights")
  
  n = 2000
  n_test = 25
  p = 12
  p_act = 4

  pi = c(seq(1, 0.1, -(1 / p_act)) * rep(c(1, -1), p_act / 2), rep(0, p - p_act))
  cov_mat = toeplitz(0.7^(0:(p - 1)))

  X_tr = mvtnorm::rmvnorm(n = n, mean = rep(0, p), sigma = cov_mat)
  X_te = mvtnorm::rmvnorm(n = n_test, mean = rep(0, p), sigma = cov_mat)
  Y_tr = as.vector(X_tr %*% pi + rnorm(n, 0, 1))

  # first ensemble
  methods_a = list("ols" = create_method("ols"),
              "plasso" = create_method("plasso"),
              "mean" = create_method("mean"))

  ens_core_a = ensemble_core(methods_a, X_tr = X_tr, Y_tr = Y_tr)
  pred_a = predict(ens_core_a, methods_a, X_tr, X_te)
  weights_a = OutcomeWeights:::get_smoother_weights.ensemble_core(ens_core_a, methods_a, X_tr, X_te)
  expect_equal(dim(weights_a), c(nrow(X_te), nrow(X_tr), length(methods_a)))

  pred_s_a = sapply(1:length(methods_a), function(i) weights_a[, , i] %*% Y_tr)
  expect_equal(pred_a, pred_s_a, tolerance = 1e-9)

  weights_sums_a = sapply(1:length(methods_a), function(i) rowSums(weights_a[, , i]))
  expect_equal(matrix(1, nrow = nrow(X_te), ncol = length(methods_a)), weights_sums_a, tolerance = 1e-9)


  # second ensemble
  methods_b = list("ols" = create_method("ols"),
              "plasso" = create_method("plasso"),
              "ridge" = create_method("ridge"),
              "forest_grf" = create_method("forest_grf"),
              "knn_10" = create_method("knn", arguments = list("k" = 10)))

  ens_core_b = ensemble_core(methods_b, X_tr = X_tr, Y_tr = Y_tr)
  pred_b = predict(ens_core_b, methods_b, X_tr, X_te)
  weights_b = OutcomeWeights:::get_smoother_weights.ensemble_core(ens_core_b, methods_b, X_tr, X_te)
  expect_equal(dim(weights_b), c(nrow(X_te), nrow(X_tr), length(methods_b)))

  pred_s_b = sapply(1:length(methods_b), function(i) weights_b[, , i] %*% Y_tr)
  expect_equal(pred_b, pred_s_b, tolerance = 1e-3)

  weights_sums_b = sapply(1:length(methods_b), function(i) rowSums(weights_b[, , i]))
  expect_equal(matrix(1, nrow = nrow(X_te), ncol = length(methods_b)), weights_sums_b, tolerance = 1e-9)

})
