test_that("check predictions and smoother weights of ensemble_short function", {

  library(mvtnorm)

  n = 200
  p = 12
  p_act = 4

  pi = c(seq(1, 0.1, -(1 / p_act)) * rep(c(1, -1), p_act / 2), rep(0, p - p_act))
  cov_mat = toeplitz(0.7^(0:(p - 1)))

  X = mvtnorm::rmvnorm(n = n, mean = rep(0, p), sigma = cov_mat)
  Y = X %*% pi + rnorm(n, 0, 1)

  cf_mat = prep_cf_mat(n, cf = 5)

  methods = list("ols" = create_method("ols"),
                "forest_grf" = create_method("forest_grf"),
                "knn" = create_method("knn", arguments = list("k" = 3)))

  ens = ensemble_short(methods, X = X, Y = Y, cf_mat = cf_mat, storeModels = "Memory")
  pred = predict(ens)
  weights_mat = weights(ens, methods = methods, X = X, Y = Y, cf_mat = cf_mat)
  pred_s = as.vector(weights_mat %*% Y)

  weights_rows = rep(1, n)
  weights_mat_rows = Matrix::rowSums(weights_mat)

  expect_equal(weights_rows, weights_mat_rows, tolerance = 1e-5)
  expect_equal(pred, pred_s, tolerance = 1e-5)

})

## Model saving funcitonality: moved to the upper level

# test_that("check storage of ensemble models and smoother weights production", {
# 
#   library(mvtnorm)
# 
#   n = 300
#   p = 6
#   p_act = 4
# 
#   pi = c(seq(1, 0.1, -(1 / p_act)) * rep(c(1, -1), p_act / 2), rep(0, p - p_act))
#   cov_mat = toeplitz(0.7^(0:(p - 1)))
# 
#   X = mvtnorm::rmvnorm(n = n, mean = rep(0, p), sigma = cov_mat)
#   Y = X %*% pi + rnorm(n, 0, 1)
# 
#   cf_mat = prep_cf_mat(n, cf = 3)
# 
#   methods = list("ols" = create_method("ols"),
#                 "forest_grf" = create_method("forest_grf"),
#                 "mean" = create_method("mean"))
# 
#   ens_mem = ensemble_short(methods, X = X, Y = Y, cf_mat = cf_mat, storeModels = "Memory")
#   weights_mat_mem = weights(ens_mem, methods = methods, X = X, Y = Y, cf_mat = cf_mat)
# 
#   ens_disk = ensemble_short(methods, X = X, Y = Y, cf_mat = cf_mat, storeModels = "Disk", path = "..")
#   expect_true(is.character(ens_disk$methods))
#   expect_true(file.exists(ens_disk$methods))
#   weights_mat_disk = weights(ens_mem, methods = methods, X = X, Y = Y, cf_mat = cf_mat)
# 
#   ens_no = ensemble_short(methods, X = X, Y = Y, cf_mat = cf_mat, storeModels = "No")
#   expect_error(weights(ens_no, methods = methods, X = X, Y = Y, cf_mat = cf_mat), "Ensemble models were not saved after training.")
# 
#   expect_identical(weights_mat_mem, weights_mat_disk)
# 
#   file.remove(ens_disk$methods)
# 
# })

## S-learner functionality: temporarily suspended

# test_that("check t, s leaner and their combination", {
# 
#   library(mvtnorm)
# 
#   quiet = TRUE
#   n = 2000
#   p = 7
#   p_act = 4
# 
#   pi = c(seq(1, 0.1, -(1 / p_act)) * rep(c(1, -1), floor(p_act / 2)), rep(0, p - p_act))
#   cov_mat = toeplitz(0.7^(0:(p - 1)))
# 
#   X = mvtnorm::rmvnorm(n = n, mean = rep(0, p), sigma = cov_mat)
#   Y = X %*% pi + rnorm(n, 2, 3)
# 
#   cf_mat = prep_cf_mat(n, cf = 3)
# 
#   subset = sample(c(TRUE, FALSE), length(Y), replace = TRUE)
# 
#   methods = list("ols" = create_method("ols"),
#                 "ridge" = create_method("ridge"),
#                 "forest_drf" = create_method("forest_drf"))
# 
#   ens_t = ensemble_short(methods, X = X, Y = Y, subset = subset, cf_mat = cf_mat, learner = "t", storeModels = "Memory")
#   nnls_w_t = nnls_weights(ens_t$fit_cv[subset, ], Y[subset])
#   ens_s = ensemble_short(methods, X = X, Y = Y, subset = subset, cf_mat = cf_mat, learner = "s", storeModels = "Memory")
#   nnls_w_s = nnls_weights(ens_s$fit_cv[subset, ], Y[subset])
#   ens_both = ensemble_short(methods, X = X, Y = Y, subset = subset, cf_mat = cf_mat, learner = "both", storeModels = "Memory")
#   nnls_w_both = nnls_weights(ens_both$fit_cv[subset, ], Y[subset])
# 
#   expect_identical(ncol(ens_t$fit_cv), length(methods))
#   expect_identical(ncol(ens_s$fit_cv), length(methods))
#   expect_identical(ncol(ens_both$fit_cv), as.integer(length(methods)*2))
# 
#   pred_t = predict(ens_t, nnls_w_t)
#   weights_mat_t = weights(ens_t, methods = methods, X = X, Y = Y, subset = subset, w = nnls_w_t, cf_mat = cf_mat)
#   pred_s_t = as.vector(weights_mat_t %*% Y)
# 
#   pred_s = predict(ens_s, nnls_w_s)
#   weights_mat_s = weights(ens_s, methods = methods, X = X, Y = Y, subset = subset, w = nnls_w_s, cf_mat = cf_mat)
#   pred_s_s = as.vector(weights_mat_s %*% Y)
# 
#   pred_both = predict(ens_both, nnls_w_both)
#   weights_mat_both = weights(ens_both, methods = methods, X = X, Y = Y, subset = subset, w = nnls_w_both, cf_mat = cf_mat)
#   pred_s_both = as.vector(weights_mat_both %*% Y)
# 
#   weights_rows = rep(1, n)
# 
#   expect_equal(weights_rows, Matrix::rowSums(weights_mat_t), tolerance = 1e-5)
#   expect_equal(weights_rows, Matrix::rowSums(weights_mat_s), tolerance = 1e-5)
#   expect_equal(weights_rows, Matrix::rowSums(weights_mat_both), tolerance = 1e-5)
# 
#   expect_equal(pred_t, pred_s_t, tolerance = 1e-5)
#   expect_equal(pred_s, pred_s_s, tolerance = 1e-5)
#   expect_equal(pred_both, pred_s_both, tolerance = 1e-5)
# 
# })
