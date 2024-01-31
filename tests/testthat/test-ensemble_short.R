test_that("check predictions and smoother weights of ensemble_short function", {

  library(mvtnorm)

  n = 1000
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
  p_act = 3

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
