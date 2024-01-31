test_that("check nuisance_cf", {

  library(mvtnorm)

  n = 2000
  p = 8
  p_act = 4

  pi = c(seq(1, 0.1, -(1 / p_act)) * rep(c(1, -1), p_act / 2), rep(0, p - p_act))
  cov_mat = toeplitz(0.7^(0:(p - 1)))

  X = mvtnorm::rmvnorm(n = n, mean = rep(0, p), sigma = cov_mat)
  y = X %*% pi + rnorm(n, 0, 1)

  cf_mat = prep_cf_mat(n, cf = 3)

  ml = list("ols" = create_method("ols"),
            "forest_drf" = create_method("forest_drf"),
            "mean" = create_method("mean"))

  path = paste0("../weights_test.rds")


  t = Sys.time()
  np_standard = nuisance_cf(ml = ml, y = y, x = X, cf_mat = cf_mat,
                            cv = 5, subset = NULL, weights = FALSE, path = path, quiet = TRUE)
  t_standard = (Sys.time() - t) %>% as.numeric(units = "secs")
  expect_false(file.exists(path))


  t = Sys.time()
  np_short = nuisance_cf(ml = ml, y = y, x = X, cf_mat = cf_mat,
                         cv = 1, subset = NULL, weights = FALSE, path = path, quiet = TRUE)
  t_short = (Sys.time() - t) %>% as.numeric(units = "secs")
  expect_false(file.exists(path))


  t = Sys.time()
  np_standard_w = nuisance_cf(ml = ml, y = y, x = X, cf_mat = cf_mat,
                              cv = 5, subset = NULL, weights = TRUE, path = path, quiet = TRUE)
  t_standard_w = (Sys.time() - t) %>% as.numeric(units = "secs")
  expect_true(file.exists(path))
  w = readRDS(path)
  fold = cf_mat[, 1]
  expect_identical(as.vector(w[fold, fold]), rep(0, sum(fold)^2))
  file.remove(path)


  t = Sys.time()
  np_short_w = nuisance_cf(ml = ml, y = y, x = X, cf_mat = cf_mat,
                           cv = 1, subset = NULL, weights = TRUE, path = path, quiet = TRUE)
  t_short_w = (Sys.time() - t) %>% as.numeric(units = "secs")
  expect_true(file.exists(path))
  w = readRDS(path)
  fold = cf_mat[, 1]
  expect_identical(as.vector(w[fold, fold]), rep(0, sum(fold)^2))
  file.remove(path)


  # check computational time
  expect_true(t_short < t_short_w)
  expect_true(t_standard < t_standard_w)
  expect_true(t_short < t_standard)
  expect_true(t_short_w < t_standard_w)


  # check that nuisance parameters are better than simple mean prediction
  naive = rep(mean(y), length(y))
  rmse_naive = sqrt(mean((y - naive)^2))
  rmse_standard = sqrt(mean((y - np_standard)^2))
  rmse_short = sqrt(mean((y - np_short)^2))
  rmse_standard_w = sqrt(mean((y - np_standard_w)^2))
  rmse_short_w = sqrt(mean((y - np_short_w)^2))

  expect_true(rmse_naive > rmse_standard)
  expect_true(rmse_naive > rmse_short)
  expect_true(rmse_naive > rmse_standard_w)
  expect_true(rmse_naive > rmse_short_w)

})
