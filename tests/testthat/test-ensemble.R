test_that("check predictions and smoother weights of ensemble function", {
  
  skip_if_not_installed("OutcomeWeights")
  
  n = 600
  n_test = 200
  p = 12
  p_act = 4
  ensemble_type <- "nnls"
  
  pi = c(seq(1, 0.1, -(1 / p_act)) * rep(c(1, -1), p_act / 2), rep(0, p - p_act))
  cov_mat = toeplitz(0.7^(0:(p - 1)))

  X = mvtnorm::rmvnorm(n = n, mean = rep(0, p), sigma = cov_mat)
  Xnew = mvtnorm::rmvnorm(n = n_test, mean = rep(0, p), sigma = cov_mat)
  Y = as.vector(X %*% pi + rnorm(n, 0, 1))

  methods = list("ols" = create_method("ols"),
                 "forest_grf" = create_method("forest_grf"),
                 "ridge" = create_method("ridge"))

  ens = ensemble(methods = methods, X = X, Y = Y, nfolds = 4, ensemble_type = ensemble_type)
  pred = predict(object = ens, methods = methods, X = X, Xnew = Xnew)
  weights_mat = OutcomeWeights:::get_smoother_weights.ensemble(object = ens, methods = methods, ens_w = ens$ens_w, X = X, Xnew = Xnew)

  pred_s = as.vector(weights_mat %*% Y)

  weights_rows = rep(1, n_test)
  weights_mat_rows = Matrix::rowSums(weights_mat)

  expect_equal(weights_rows, weights_mat_rows, tolerance = 1e-5)
  expect_equal(pred$np, pred_s, tolerance = 1e-3)

})
