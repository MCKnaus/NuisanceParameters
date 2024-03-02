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

  path = paste0(gsub("\\\\", "/", tempdir()), "/Ensemble_Y")
  path_rds = paste0(path, ".rds")
  path_weights = paste0(path, "_Weights.rds")
  unlink(paste0(dirname(path), "/*"))


  t = Sys.time()
  np_standard = nuisance_cf(ml = ml, y = y, x = X, cf_mat = cf_mat,
                            cv = 5, subset = NULL, weights = FALSE, path = path, quiet = TRUE)
  t_standard = (Sys.time() - t) %>% as.numeric(units = "secs")
  expect_true(file.exists(path_rds))
  expect_false(file.exists(path_weights))
  unlink(paste0(dirname(path), "/*"))


  t = Sys.time()
  np_short = nuisance_cf(ml = ml, y = y, x = X, cf_mat = cf_mat,
                         cv = 1, subset = NULL, weights = FALSE, path = path, quiet = TRUE)
  t_short = (Sys.time() - t) %>% as.numeric(units = "secs")
  expect_true(file.exists(path_rds))
  expect_false(file.exists(path_weights))
  unlink(paste0(dirname(path), "/*"))


  t = Sys.time()
  np_standard_w = nuisance_cf(ml = ml, y = y, x = X, cf_mat = cf_mat,
                              cv = 5, subset = NULL, weights = TRUE, path = path, quiet = TRUE)
  t_standard_w = (Sys.time() - t) %>% as.numeric(units = "secs")
  expect_true(file.exists(path_rds))
  w = readRDS(path_weights)
  fold = cf_mat[, 1]
  expect_identical(as.vector(w[fold, fold]), rep(0, sum(fold)^2))
  expect_equal(np_standard_w, as.vector(w %*% y), tolerance = 1e-5)
  unlink(paste0(dirname(path), "/*"))


  t = Sys.time()
  np_short_w = nuisance_cf(ml = ml, y = y, x = X, cf_mat = cf_mat,
                           cv = 1, subset = NULL, weights = TRUE, path = path, quiet = TRUE)
  t_short_w = (Sys.time() - t) %>% as.numeric(units = "secs")
  expect_true(file.exists(path_rds))
  w = readRDS(path_weights)
  fold = cf_mat[, 1]
  expect_identical(as.vector(w[fold, fold]), rep(0, sum(fold)^2))
  expect_equal(np_short_w, as.vector(w %*% y), tolerance = 1e-5)
  unlink(paste0(dirname(path), "/*"))


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


test_that("check nuisance_w", {

  library(mvtnorm)
  set.seed(234)

  n = 1000
  p = 6
  cov_mat = toeplitz(0.7^(0:(6 - 1)))

  X = mvtnorm::rmvnorm(n = n, mean = rep(2, p), sigma = cov_mat)
  x1 = X[, 1]; x2 = X[, 2]; x3 = X[, 3]; x4 = X[, 4]; x5 = X[, 5]; x6 = X[, 6];

  effect_class1 = sin(x1) + x2^2 + x3*x4
  effect_class2 = cos(x2) + x4^2 + x5*x6
  effect_class3 = x1^2 + exp(x3) + x2*x5
  effect_class1 = (effect_class1 - mean(effect_class1)) / sd(effect_class1)
  effect_class2 = (effect_class2 - mean(effect_class2)) / sd(effect_class2)
  effect_class3 = (effect_class3 - mean(effect_class3)) / sd(effect_class3)
  raw_scores = cbind(effect_class1, effect_class2, effect_class3)
  prob_classes = t(apply(raw_scores, 1, function(x) exp(x) / sum(exp(x))))

  w = apply(prob_classes, 1, function(p) sample(1:3, 1, prob = p))

  w_mat = prep_w_mat(w)
  cf_mat = prep_cf_mat(n, cf = 2, w_mat = w_mat)

  ml = list("ols" = create_method("ols", name = "OLS"),
            "forest_grf" = create_method("forest_grf", name = "GRF"),
            "ridge" = create_method("ridge", name = "Ridge"),
            "knn" = create_method("knn", name = "KNN"))

  path = gsub("\\\\", "/", tempdir())
  unlink(paste0(path, "/*"))


  ### Short-Stacking ###
  expect_message(np_e_short <- nuisance_e(ml = ml, w_mat = w_mat, x = X, cf_mat = cf_mat, cv = 1,
                                          path = path, quiet = FALSE), "Short-stacking is used.")


  files = paste0(path, "/Ensemble_W", 1:3, ".rds")
  files_w = paste0(path, "/Ensemble_W", 1:3, "_Weights.rds")

  # check if ensemble output is stored in files
  expect_true(all(file.exists(files)))
  expect_false(any(file.exists(files_w)))

  # check if ensemble output contains correct info
  ens = readRDS(files[1])
  expect_length(ens, 2)
  expect_named(ens$nnls_w)
  expect_equal(dim(ens$fit_cv), c(n, length(ml)))
  expect_true(all(substr(colnames(ens$fit_cv), 1, 1) == "t"))

  # check for correct dimension
  expect_identical(dim(np_e_short), dim(w_mat))

  # check if probabilities sum up to 1 for every obs
  expect_equal(rep(1, n), rowSums(np_e_short), tolerance = 1e-9)

  unlink(paste0(path, "/*"))


  ### Standard-Stacking ###
  cv = 3
  expect_message(np_e_standard <- nuisance_e(ml = ml, w_mat = w_mat, x = X, cf_mat = cf_mat, cv = cv,
                                             path = path, quiet = FALSE), "Standard-stacking is used.")

  files = paste0(path, "/Ensemble_W", 1:3, ".rds")
  files_w = paste0(path, "/Ensemble_W", 1:3, "_Weights.rds")

  # check if ensemble output is stored in files
  expect_true(all(file.exists(files)))
  expect_false(any(file.exists(files_w)))

  # check if ensemble output contains correct info
  ens = readRDS(files[1])
  expect_length(ens, ncol(cf_mat))
  expect_named(ens[[1]]$nnls_w)
  expect_equal(dim(ens[[1]]$fit_cv), c(sum(cf_mat[, 1]), length(ml)))
  expect_true(all(substr(colnames(ens[[1]]$fit_cv), 1, 1) == "t"))

  # check for correct dimension
  expect_identical(dim(np_e_standard), dim(w_mat))

  # check if probabilities sum up to 1 for every obs
  expect_equal(rep(1, n), rowSums(np_e_standard), tolerance = 1e-9)

  unlink(paste0(path, "/*"))

})


test_that("check nuisance_m", {

 library(mvtnorm)
  set.seed(79)

  n = 1000
  p = 8
  p_act = 4

  pi = c(seq(1, 0.1, -(1 / p_act)) * rep(c(1, -1), p_act / 2), rep(0, p - p_act))
  cov_mat = toeplitz(0.7^(0:(p - 1)))

  x = mvtnorm::rmvnorm(n = n, mean = rep(0, p), sigma = cov_mat)
  y = as.vector(x %*% pi + rnorm(n, 0, 1))

  w_mods = 3
  w = sample(1:w_mods, n, replace = TRUE)
  w_mat = prep_w_mat(w)

  cf = 4
  cf_mat = prep_cf_mat(n, cf = cf, w_mat = w_mat)

  ml = list("ridge" = create_method("ridge"),
            "forest_grf" = create_method("forest_grf"),
            "mean" = create_method("mean"))

  path = gsub("\\\\", "/", tempdir())



  ### Short-Stacking ###
  expect_message(np_m_short <- nuisance_m(ml = ml, y = y, w_mat = w_mat, x = x, cf_mat = cf_mat,
                                          cv = 1, learner = "t", path = path, quiet = FALSE), "Short-stacking is used.")


  files = paste0(path, "/Ensemble_Y", 1:w_mods, ".rds")
  files_w = paste0(path, "/Ensemble_Y", 1:w_mods, "_Weights.rds")

  # check if ensemble output is stored in files
  expect_true(all(file.exists(files)))
  expect_false(any(file.exists(files_w)))

  # check if ensemble output contains correct info
  ens = readRDS(files[1])
  expect_length(ens, 2)
  expect_named(ens$nnls_w)
  expect_equal(dim(ens$fit_cv), c(n, length(ml)))
  expect_true(all(substr(colnames(ens$fit_cv), 1, 1) == "t"))

  # check for correct dimension
  expect_identical(dim(np_m_short), dim(w_mat))

  unlink(paste0(path, "/*"))



  ### Short-Stacking with Smoother Weights ###
  expect_message(np_m_short_w <- nuisance_m(ml = ml, y = y, w_mat = w_mat, x = x, cf_mat = cf_mat,
                                            cv = 1, learner = "t", path = path, quiet = FALSE,
                                            weights = TRUE), "Short-stacking is used.")


  files = paste0(path, "/Ensemble_Y", 1:w_mods, ".rds")
  files_w = paste0(path, "/Ensemble_Y", 1:w_mods, "_Weights.rds")

  # check if ensemble output is stored in files
  expect_true(all(file.exists(files)))
  expect_true(all(file.exists(files_w)))

  # check if ensemble output contains correct info
  ens = readRDS(files[1])
  expect_length(ens, 2)
  expect_named(ens$nnls_w)
  expect_equal(dim(ens$fit_cv), c(n, length(ml)))
  expect_true(all(substr(colnames(ens$fit_cv), 1, 1) == "t"))

  # check for correct dimension
  expect_identical(dim(np_m_short_w), dim(w_mat))

  # check smoother weights
  w = readRDS(files_w[1])
  expect_equal(np_m_short_w[, 1], as.vector(w %*% y), tolerance = 1e-3)


  unlink(paste0(path, "/*"))


  ### Standard-Stacking ###
  expect_message(np_m_standard <- nuisance_m(ml = ml, y = y, w_mat = w_mat, x = x, cf_mat = cf_mat,
                                             cv = 3, path = path, quiet = FALSE), "Standard-stacking is used.")


  files = paste0(path, "/Ensemble_Y", 1:w_mods, ".rds")
  files_w = paste0(path, "/Ensemble_Y", 1:w_mods, "_Weights.rds")

  # check if ensemble output is stored in files
  expect_true(all(file.exists(files)))
  expect_false(any(file.exists(files_w)))

  # check if ensemble output contains correct info
  ens = readRDS(files[1])
  expect_length(ens, cf)
  expect_named(ens[[1]]$nnls_w)
  expect_equal(dim(ens[[1]]$fit_cv), c(sum(cf_mat[, 1]), length(ml)))
  expect_true(all(substr(colnames(ens[[1]]$fit_cv), 1, 1) == "t"))

  # check for correct dimension
  expect_identical(dim(np_m_standard), dim(w_mat))

  unlink(paste0(path, "/*"))


  ### Standard-Stacking with Smoother Weights ###
  expect_message(np_m_standard_w <- nuisance_m(ml = ml, y = y, w_mat = w_mat, x = x, cf_mat = cf_mat,
                                               cv = 3, path = path, quiet = FALSE,
                                               weights = TRUE), "Standard-stacking is used.")


  files = paste0(path, "/Ensemble_Y", 1:w_mods, ".rds")
  files_w = paste0(path, "/Ensemble_Y", 1:w_mods, "_Weights.rds")

  # check if ensemble output is stored in files
  expect_true(all(file.exists(files)))
  expect_true(all(file.exists(files_w)))

  # check if ensemble output contains correct info
  ens = readRDS(files[1])
  expect_length(ens, cf)
  expect_named(ens[[1]]$nnls_w)
  expect_equal(dim(ens[[1]]$fit_cv), c(sum(cf_mat[, 1]), length(ml)))
  expect_true(all(substr(colnames(ens[[1]]$fit_cv), 1, 1) == "t"))

  # check for correct dimension
  expect_identical(dim(np_m_standard_w), dim(w_mat))

  # check smoother weights
  w = readRDS(files_w[1])
  expect_equal(np_m_standard_w[, 1], as.vector(w %*% y), tolerance = 1e-3)

  unlink(paste0(path, "/*"))


  ### S-Learner ###

  expect_error(
    nuisance_m(ml = ml, y = y, w_mat = w_mat, x = x, cf_mat = cf_mat,
               cv = 3, path = path, quiet = FALSE, weights = TRUE,
               learner = "s"),
    "S-Learner cannot be combined")

  expect_error(
    nuisance_m(ml = ml, y = y, w_mat = w_mat, x = x, cf_mat = cf_mat,
               cv = 1, path = path, quiet = FALSE, weights = FALSE,
               learner = "s"),
    "S-Learner cannot be combined")

  expect_error(
    nuisance_m(ml = ml, y = y, w_mat = w_mat, x = x, cf_mat = cf_mat,
               cv = 5, path = path, quiet = FALSE, weights = TRUE,
               learner = "both"),
    "S-Learner cannot be combined")

  expect_error(
    nuisance_m(ml = ml, y = y, w_mat = w_mat, x = x, cf_mat = cf_mat,
               cv = 2, path = path, quiet = FALSE, weights = FALSE,
               learner = "both"),
    "S-Learner cannot be combined")

})


