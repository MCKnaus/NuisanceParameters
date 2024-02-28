test_that("check predictions and smoother weights of ensemble_short function", {

  library(mvtnorm)

  n = 200
  p = 12
  p_act = 4

  pi = c(seq(1, 0.1, -(1 / p_act)) * rep(c(1, -1), p_act / 2), rep(0, p - p_act))
  cov_mat = toeplitz(0.7^(0:(p - 1)))

  X = mvtnorm::rmvnorm(n = n, mean = rep(0, p), sigma = cov_mat)
  y = X %*% pi + rnorm(n, 0, 1)

  cf_mat = prep_cf_mat(n, cf = 5)

  ml = list("ols" = create_method("ols"),
            "forest_grf" = create_method("forest_grf"),
            "knn" = create_method("knn", arguments = list("k" = 3)))

  ens = ensemble_short(ml, x = X, y = y, cf_mat = cf_mat, storeModels = "Memory")
  pred = predict(ens)
  w_mat = weights(ens, ml = ml, x = X, y = y, cf_mat = cf_mat)
  pred_s = as.vector(w_mat %*% y)

  w_rows = rep(1, n)
  w_mat_rows = Matrix::rowSums(w_mat)

  expect_equal(w_rows, w_mat_rows, tolerance = 1e-5)
  expect_equal(pred, pred_s, tolerance = 1e-5)

})


test_that("check storage of ensemble models and smoother weights production", {

  library(mvtnorm)

  n = 300
  p = 6
  p_act = 4

  pi = c(seq(1, 0.1, -(1 / p_act)) * rep(c(1, -1), p_act / 2), rep(0, p - p_act))
  cov_mat = toeplitz(0.7^(0:(p - 1)))

  X = mvtnorm::rmvnorm(n = n, mean = rep(0, p), sigma = cov_mat)
  y = X %*% pi + rnorm(n, 0, 1)

  cf_mat = prep_cf_mat(n, cf = 3)

  ml = list("ols" = create_method("ols"),
            "forest_grf" = create_method("forest_grf"),
            "mean" = create_method("mean"))

  ens_mem = ensemble_short(ml, x = X, y = y, cf_mat = cf_mat, storeModels = "Memory")
  w_mat_mem = weights(ens_mem, ml = ml, x = X, y = y, cf_mat = cf_mat)

  ens_disk = ensemble_short(ml, x = X, y = y, cf_mat = cf_mat, storeModels = "Disk", path = "..")
  expect_true(is.character(ens_disk$ml))
  expect_true(file.exists(ens_disk$ml))
  w_mat_disk = weights(ens_mem, ml = ml, x = X, y = y, cf_mat = cf_mat)

  ens_no = ensemble_short(ml, x = X, y = y, cf_mat = cf_mat, storeModels = "No")
  expect_error(weights(ens_no, ml = ml, x = X, y = y, cf_mat = cf_mat), "Ensemble models were not saved after training.")

  expect_identical(w_mat_mem, w_mat_disk)

  file.remove(ens_disk$ml)

})


test_that("check t, s leaner and their combination", {

  library(mvtnorm)

  quiet = TRUE
  n = 2000
  p = 7
  p_act = 4

  pi = c(seq(1, 0.1, -(1 / p_act)) * rep(c(1, -1), floor(p_act / 2)), rep(0, p - p_act))
  cov_mat = toeplitz(0.7^(0:(p - 1)))

  x = mvtnorm::rmvnorm(n = n, mean = rep(0, p), sigma = cov_mat)
  y = x %*% pi + rnorm(n, 2, 3)

  cf_mat = prep_cf_mat(n, cf = 3)

  subset = sample(c(TRUE, FALSE), length(y), replace = TRUE)

  ml = list("ols" = create_method("ols"),
            "ridge" = create_method("ridge"),
            "forest_drf" = create_method("forest_drf"))

  ens_t = ensemble_short(ml, x = x, y = y, subset = subset, cf_mat = cf_mat, learner = "t", storeModels = "Memory")
  nnls_w_t = nnls_weights(ens_t$fit_cv[subset, ], y[subset])
  ens_s = ensemble_short(ml, x = x, y = y, subset = subset, cf_mat = cf_mat, learner = "s", storeModels = "Memory")
  nnls_w_s = nnls_weights(ens_s$fit_cv[subset, ], y[subset])
  ens_both = ensemble_short(ml, x = x, y = y, subset = subset, cf_mat = cf_mat, learner = "both", storeModels = "Memory")
  nnls_w_both = nnls_weights(ens_both$fit_cv[subset, ], y[subset])

  expect_identical(ncol(ens_t$fit_cv), length(ml))
  expect_identical(ncol(ens_s$fit_cv), length(ml))
  expect_identical(ncol(ens_both$fit_cv), as.integer(length(ml)*2))

  pred_t = predict(ens_t, nnls_w_t)
  w_mat_t = weights(ens_t, ml = ml, x = x, y = y, subset = subset, w = nnls_w_t, cf_mat = cf_mat)
  pred_s_t = as.vector(w_mat_t %*% y)

  pred_s = predict(ens_s, nnls_w_s)
  w_mat_s = weights(ens_s, ml = ml, x = x, y = y, subset = subset, w = nnls_w_s, cf_mat = cf_mat)
  pred_s_s = as.vector(w_mat_s %*% y)

  pred_both = predict(ens_both, nnls_w_both)
  w_mat_both = weights(ens_both, ml = ml, x = x, y = y, subset = subset, w = nnls_w_both, cf_mat = cf_mat)
  pred_s_both = as.vector(w_mat_both %*% y)

  w_rows = rep(1, n)

  expect_equal(w_rows, Matrix::rowSums(w_mat_t), tolerance = 1e-5)
  expect_equal(w_rows, Matrix::rowSums(w_mat_s), tolerance = 1e-5)
  expect_equal(w_rows, Matrix::rowSums(w_mat_both), tolerance = 1e-5)

  expect_equal(pred_t, pred_s_t, tolerance = 1e-5)
  expect_equal(pred_s, pred_s_s, tolerance = 1e-5)
  expect_equal(pred_both, pred_s_both, tolerance = 1e-5)

})
