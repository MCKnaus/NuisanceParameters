test_that("check predictions and smoother weights of ensemble function", {

  library(mvtnorm)

  n = 600
  n_test = 200
  p = 12
  p_act = 4

  pi = c(seq(1, 0.1, -(1 / p_act)) * rep(c(1, -1), p_act / 2), rep(0, p - p_act))
  cov_mat = toeplitz(0.7^(0:(p - 1)))

  X = mvtnorm::rmvnorm(n = n, mean = rep(0, p), sigma = cov_mat)
  Xnew = mvtnorm::rmvnorm(n = n_test, mean = rep(0, p), sigma = cov_mat)
  y = X %*% pi + rnorm(n, 0, 1)

  ml = list("ols" = create_method("ols"),
            "forest_grf" = create_method("forest_grf"),
            "ridge" = create_method("ridge"))

  ens = ensemble(ml, x = X, y = y, nfolds = 4)
  pred = predict(object = ens, ml = ml, x = X, y = y, xnew = Xnew)
  w_mat = weights(object = ens, ml = ml, x = X, y = y, xnew = Xnew)

  pred_s = as.vector(w_mat %*% y)

  w_rows = rep(1, n_test)
  w_mat_rows = Matrix::rowSums(w_mat)

  expect_equal(w_rows, w_mat_rows, tolerance = 1e-5)
  expect_equal(pred, pred_s, tolerance = 1e-5)

})
