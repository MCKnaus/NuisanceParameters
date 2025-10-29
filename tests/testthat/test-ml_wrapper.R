test_that("smoother weights lead to correct prediction of fitted values", {
  
  skip_if_not_installed("OutcomeWeights")
  
  n = 1000
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
  p = predict.mean_fit(m, Xnew = Xnew)
  w = OutcomeWeights:::weights.mean_fit(m, Xnew = Xnew)
  p_s = as.vector(w %*% Y)
  expect_equal(Matrix::rowSums(w), w_rows, tolerance = 1e-9)
  expect_true(all(w >= 0 & w <= 1))
  expect_equal(p, p_s, tolerance = 1e-9)
  
  # ols
  m = ols_fit(X = X, Y = Y)
  p = predict.ols_fit(m, Xnew = Xnew)
  w = OutcomeWeights:::weights.ols_fit(m, Xnew = Xnew)
  p_s = as.vector(w %*% Y)
  expect_equal(Matrix::rowSums(w), w_rows, tolerance = 1e-9)
  expect_equal(p, p_s, tolerance = 1e-9)
  
  # ridge
  m = ridge_fit(X = X, Y = Y)
  p = predict.ridge_fit(m, Xnew = Xnew)
  w = OutcomeWeights:::weights.ridge_fit(m, Xnew = Xnew)
  p_s = as.vector(w %*% Y)
  expect_equal(Matrix::rowSums(w), w_rows, tolerance = 1e-9)
  expect_equal(p, p_s, tolerance = 0.001)
  
  # plasso
  m = plasso_fit(X = X, Y = Y)
  p = predict.plasso_fit(m, Xnew = Xnew)
  w = OutcomeWeights:::weights.plasso_fit(m, Xnew = Xnew)
  p_s = as.vector(w %*% Y)
  expect_equal(Matrix::rowSums(w), w_rows, tolerance = 1e-9)
  expect_equal(p, p_s, tolerance = 1e-9)
  
  # rlasso
  m = rlasso_fit(X = X, Y = Y)
  p = predict.rlasso_fit(m, Xnew = Xnew)
  w = OutcomeWeights:::weights.rlasso_fit(m, Xnew = Xnew)
  p_s = as.vector(w %*% Y)
  expect_equal(Matrix::rowSums(w), w_rows, tolerance = 1e-9)
  expect_equal(p, p_s, tolerance = 1e-9)
  
  # random forest
  m = forest_grf_fit(X = X, Y = Y)
  p = predict.forest_grf_fit(m, Xnew = Xnew)
  w = OutcomeWeights:::weights.forest_grf_fit(m, Xnew = Xnew)
  p_s = as.vector(w %*% Y)
  expect_equal(Matrix::rowSums(w), w_rows, tolerance = 1e-9)
  expect_true(all(w >= 0 & w <= 1))
  expect_equal(p, p_s, tolerance = 1e-9)
  
  # knn
  m = knn_fit(X = X, Y = Y, arguments = list("k" = 3))
  p = predict.knn_fit(m, Xnew = Xnew)
  w = OutcomeWeights:::weights.knn_fit(m, X = X, Xnew = Xnew)
  p_s = as.vector(w %*% Y)
  expect_equal(Matrix::rowSums(w), w_rows, tolerance = 1e-9)
  expect_true(all(w >= 0 & w <= 1))
  expect_equal(p, p_s, tolerance = 1e-9)
  
  # distributional random forest
  m = forest_drf_fit(X = X, Y = Y)
  p = predict.forest_drf_fit(m, Xnew = Xnew)
  w = OutcomeWeights:::weights.forest_drf_fit(m, Xnew = Xnew)
  p_s = as.vector(w %*% Y)
  expect_equal(Matrix::rowSums(w), w_rows, tolerance = 1e-9)
  expect_true(all(w >= 0 & w <= 1))
  expect_equal(p, p_s, tolerance = 1e-9)
  
  # xgboost
  m = xgboost_fit(X = X, Y = Y)
  p = predict.xgboost_fit(m, Xnew = Xnew)
  w = OutcomeWeights:::weights.xgboost_fit(m, Xnew = Xnew)
  p_s = as.vector(w %*% Y)
  expect_equal(p, p_s, tolerance = 1e-6)
})


test_that("knn mickey mouse check", {
  
  skip_if_not_installed("OutcomeWeights")

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

  m = knn_fit(X = X, Y = Y, arguments = list("k" = k))
  w = OutcomeWeights:::weights.knn_fit(m, X = X, Xnew = Xnew)

  expect_identical(w_exp, w)

  p = predict.knn_fit(m, Xnew = Xnew)
  expect_identical(p, c(1, 3))

})


test_that("knn weight matrix plausibility check", {

  skip_if_not_installed("OutcomeWeights")
  
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
  m = knn_fit(X = X, Y = Y, arguments = list("k" = k))
  p = predict.knn_fit(m, Xnew = Xnew)
  expect_true(all(p >= min(Y) & p <= max(Y)))

  w = OutcomeWeights:::weights.knn_fit(m, X = X, Xnew = Xnew)
  w_vals = unique(as.vector(as.matrix(w[-1])))
  expect_identical(sort(w_vals), c(0, 1/k))

})


test_that("(rough) equality of plasso and ols prediction solution for clear cut dgp", {
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
  p_ols = predict.ols_fit(m_ols, Xnew = Xnew)

  m_p = plasso_fit(X = X, Y = Y)
  p_p = predict.plasso_fit(m_p, Xnew = Xnew)

  expect_equal(p_ols, p_p, tolerance = 0.01)

})


set.seed(42)
X <- matrix(rnorm(50 * 5), ncol = 5)
colnames(X) <- paste0("X", 1:5)
Y_reg <- rnorm(50)
Y_bin <- sample(0:1, 50, replace = TRUE)
Y_multi <- sample(1:3, 50, replace = TRUE)

test_that("mean_fit and predict.mean_fit work", {
  fit <- mean_fit(X, Y_reg)
  expect_s3_class(fit, "mean_fit")
  pred <- predict(fit, X)
  pred_no_Xnew <- predict(fit, Xnew = NULL)
  expect_equal(pred, pred_no_Xnew)
  expect_equal(length(pred), nrow(X))
  expect_true(all(pred == mean(Y_reg)))
})

test_that("ols_fit and predict.ols_fit work", {
  fit <- ols_fit(X, Y_reg)
  expect_s3_class(fit, "ols_fit")
  pred <- predict(fit, X)
  pred_no_Xnew <- predict(fit, Xnew = NULL)
  expect_equal(pred, pred_no_Xnew)
  expect_equal(length(pred), nrow(X))
})

test_that("ridge_fit and predict.ridge_fit work", {
  skip_if_not_installed("glmnet")
  fit <- ridge_fit(X, Y_reg)
  expect_s3_class(fit, "ridge_fit")
  pred <- predict.ridge_fit(object = fit, Xnew = X)
  pred_no_Xnew <- predict.ridge_fit(fit, Xnew = NULL)
  expect_equal(pred, pred_no_Xnew)
  expect_equal(length(pred), nrow(X))
})

test_that("plasso_fit and predict.plasso_fit work", {
  skip_if_not_installed("plasso")
  fit <- plasso_fit(X, Y_reg)
  expect_s3_class(fit, "plasso_fit")
  pred <- predict.plasso_fit(fit, X)
  pred_no_Xnew <- predict.plasso_fit(fit, Xnew = NULL)
  expect_equal(pred, pred_no_Xnew)
  expect_equal(length(pred), nrow(X))
})

test_that("rlasso_fit and predict.rlasso_fit work", {
  skip_if_not_installed("hdm")
  fit <- rlasso_fit(X, Y_reg)
  expect_s3_class(fit, "rlasso_fit")
  pred <- predict.rlasso_fit(fit, X)
  pred_no_Xnew <- predict.rlasso_fit(fit, Xnew = NULL)
  expect_equal(pred, pred_no_Xnew)
  expect_equal(length(pred), nrow(X))
})

test_that("forest_grf_fit and predict.forest_grf_fit work", {
  skip_if_not_installed("grf")
  fit <- forest_grf_fit(X, Y_reg)
  expect_s3_class(fit, "forest_grf_fit")
  pred <- predict.forest_grf_fit(fit, X)
  pred_no_Xnew <- predict.forest_grf_fit(fit, Xnew = NULL)
  expect_equal(pred, pred_no_Xnew)
  expect_equal(length(pred), nrow(X))
})

test_that("xgboost_fit and predict.xgboost_fit work", {
  skip_if_not_installed("xgboost")
  fit <- xgboost_fit(X, Y_reg)
  expect_s3_class(fit, "xgboost_fit")
  pred <- predict.xgboost_fit(fit, X)
  pred_no_Xnew <- predict.xgboost_fit(fit, Xnew = NULL)
  expect_equal(pred, pred_no_Xnew)
  expect_equal(length(pred), nrow(X))
})

test_that("lasso_fit and predict.lasso_fit work", {
  skip_if_not_installed("glmnet")
  fit <- lasso_fit(X, Y_reg)
  expect_s3_class(fit, "lasso_fit")
  pred <- predict.lasso_fit(fit, X)
  pred_no_Xnew <- predict.lasso_fit(fit, Xnew = NULL)
  expect_equal(pred, pred_no_Xnew)
  expect_equal(length(pred), nrow(X))
})

test_that("forest_drf_fit and predict.forest_drf_fit work", {
  skip_if_not_installed("drf")
  fit <- forest_drf_fit(X, Y_reg)
  expect_s3_class(fit, "forest_drf_fit")
  pred <- predict.forest_drf_fit(fit, X)
  pred_no_Xnew <- predict.forest_drf_fit(fit, Xnew = NULL)
  expect_equal(pred, pred_no_Xnew)
  expect_equal(length(pred), nrow(X))
})

test_that("knn_fit and predict.knn_fit work", {
  fit <- knn_fit(X = X, Y = Y_reg, arguments = list(k = 5))
  expect_s3_class(fit, "knn_fit")
  pred <- predict.knn_fit(fit, X)
  pred_no_Xnew <- predict.knn_fit(fit, Xnew = NULL)
  expect_equal(pred, pred_no_Xnew)
  expect_equal(length(pred), nrow(X))
})

test_that("logit_fit and predict.logit_fit work for binary", {
  skip_if_not_installed("glmnet")
  fit <- logit_fit(X, Y_bin)
  pred <- predict.logit_fit(fit, X, Y_bin, X)
  pred_no_Xnew <- predict.logit_fit(fit, Xnew = NULL)
  expect_equal(pred, pred_no_Xnew)
  expect_true(is.vector(pred))
  })

test_that("logit_fit and predict.logit_fit work for multiclass", {
  skip_if_not_installed("glmnet")
  fit <- logit_fit(X, Y_multi)
  pred <- predict.logit_fit(fit, X)
  pred_no_Xnew <- predict.logit_fit(fit, Xnew = NULL)
  expect_equal(pred, pred_no_Xnew)
  expect_true(all(dim(pred) == c(nrow(X), length(unique(Y_multi)), 1)))
})

test_that("logit_nnet_fit and predict.logit_nnet_fit work", {
  skip_if_not_installed("nnet")
  fit <- logit_nnet_fit(X, Y_multi)
  pred <- predict.logit_nnet_fit(fit, X)
  pred_no_Xnew <- predict.logit_nnet_fit(fit, Xnew = NULL)
  expect_equal(pred, pred_no_Xnew)
  expect_true(all(dim(as.matrix(pred)) == c(nrow(X), length(unique(Y_multi)))))
})

test_that("nb_gaussian_fit and predict.nb_gaussian_fit work", {
  skip_if_not_installed("naivebayes")
  fit <- nb_gaussian_fit(X, Y_bin)
  pred <- predict.nb_gaussian_fit(fit, X)
  pred_no_Xnew <- predict.nb_gaussian_fit(fit, Xnew = NULL)
  expect_equal(pred, pred_no_Xnew)
  expect_true(is.vector(pred) || is.matrix(pred))
})

test_that("nb_bernoulli_fit and predict.nb_bernoulli_fit work", {
  skip_if_not_installed("naivebayes")
  fit <- nb_bernoulli_fit(X, Y_bin)
  pred <- predict.nb_bernoulli_fit(fit, X)
  pred_no_Xnew <- predict.nb_bernoulli_fit(fit, Xnew = NULL)
  expect_equal(pred, pred_no_Xnew)
  expect_true(is.vector(pred) || is.matrix(pred))
})

test_that("xgboost_prop_fit and predict.xgboost_prop_fit work", {
  skip_if_not_installed("xgboost")
  fit <- xgboost_prop_fit(X, Y_bin)
  pred <- predict.xgboost_prop_fit(fit, X)
  pred_no_Xnew <- predict.xgboost_prop_fit(fit, Xnew = NULL)
  expect_equal(pred, pred_no_Xnew)
  expect_true(is.vector(pred) || is.matrix(pred))
})

test_that("svm_fit and predict.svm_fit work", {
  skip_if_not_installed("e1071")
  fit <- svm_fit(X, Y_bin)
  pred <- predict.svm_fit(object = fit, Xnew = X)
  pred_no_Xnew <- predict.svm_fit(fit, Xnew = NULL)
  expect_equal(pred, pred_no_Xnew)
  expect_true(is.vector(pred) || is.matrix(pred))
})

test_that("prob_forest_fit and predict.prob_forest_fit work", {
  skip_if_not_installed("grf")
  fit <- prob_forest_fit(X, Y_bin)
  pred <- predict.prob_forest_fit(fit, X)
  pred_no_Xnew <- predict.prob_forest_fit(fit, Xnew = NULL)
  expect_equal(pred, pred_no_Xnew)
  expect_true(is.vector(pred) || is.matrix(pred))
})

test_that("knn_prop_fit and predict.knn_prop_fit work", {
  skip_if_not_installed("kknn")
  fit <- knn_prop_fit(X = X, Y = Y_multi)
  pred <- predict.knn_prop_fit(object = fit, Xnew = X)
  pred_no_Xnew <- predict.knn_prop_fit(fit, Xnew = NULL)
  expect_equal(pred, pred_no_Xnew)
  expect_true(is.vector(pred) || is.matrix(pred))
})

test_that("ranger_fit and predict.ranger_fit work", {
  skip_if_not_installed("ranger")
  fit <- ranger_fit(X, Y_multi)
  pred <- predict.ranger_fit(fit, X)
  pred_no_Xnew <- predict.ranger_fit(fit, Xnew = NULL)
  expect_equal(pred, pred_no_Xnew)
  expect_true(is.matrix(pred) || is.data.frame(pred))
})

