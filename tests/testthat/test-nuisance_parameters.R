test_that("check nuisance_parameters for Y.hat estimation", {

  library(mvtnorm)

  n = 2000
  p = 8
  p_act = 4

  pi = c(seq(1, 0.1, -(1 / p_act)) * rep(c(1, -1), p_act / 2), rep(0, p - p_act))
  cov_mat = toeplitz(0.7^(0:(p - 1)))

  X = mvtnorm::rmvnorm(n = n, mean = rep(0, p), sigma = cov_mat)
  Y = X %*% pi + rnorm(n, 0, 1)
  
  cf = 3

  methods = list("ols" = create_method("ols"),
            "forest_drf" = create_method("forest_drf"),
            "mean" = create_method("mean"))

  path = paste0(gsub("\\\\", "/", tempdir()))
  unlink(paste0(dirname(path), "/*"))

  t = Sys.time()
  np_standard <- nuisance_parameters(NuPa="Y.hat", X, Y, methods = methods, cf = cf, stacking = 5, storeModels = "No", path = path)
  t_standard = (Sys.time() - t) %>% as.numeric(units = "secs")
  unlink(paste0(dirname(path), "/*"))


  t = Sys.time()
  np_short <- nuisance_parameters(NuPa="Y.hat", X, Y, methods = methods, cf = cf, stacking = "short", storeModels = "No", path = path)
  t_short = (Sys.time() - t) %>% as.numeric(units = "secs")
  unlink(paste0(dirname(path), "/*"))


  t = Sys.time()
  np_standard_w = nuisance_parameters(NuPa="Y.hat", X, Y, methods = methods, cf = cf, stacking = 5, storeModels = "Disk", path = path)
  t_standard_w = (Sys.time() - t) %>% as.numeric(units = "secs")
  expect_true(file.exists(path))
  w = get_outcome_weights(object = file.path(path, "nuisance_models.RDS"), NuPa = "Y.hat")$Y.hat
  fold = np_standard_w$numbers$cf_mat[, 1]
  expect_identical(as.vector(w[fold, fold]), rep(0, sum(fold)^2))
  expect_equal(np_standard_w[["nuisance_parameters"]][["Y.hat"]], as.vector(w %*% Y), tolerance = 1e-5)
  unlink(paste0(dirname(path), "/*"))


  t = Sys.time()
  np_short_w = nuisance_parameters(NuPa="Y.hat", X, Y, methods = methods, cf = cf, stacking = "short", storeModels = "Disk", path = path)
  t_short_w = (Sys.time() - t) %>% as.numeric(units = "secs")
  expect_true(file.exists(path))
  w = get_outcome_weights(object = file.path(path, "nuisance_models.RDS"), NuPa = "Y.hat")$Y.hat
  fold = np_short_w$numbers$cf_mat[, 1]
  expect_identical(as.vector(w[fold, fold]), rep(0, sum(fold)^2))
  expect_equal(np_short_w[["nuisance_parameters"]][["Y.hat"]], as.vector(w %*% Y), tolerance = 1e-5)
  unlink(paste0(dirname(path), "/*"))


  # check computational time
  expect_true(t_short < t_short_w)
  expect_true(t_standard < t_standard_w)
  expect_true(t_short < t_standard)
  expect_true(t_short_w < t_standard_w)


  # check that nuisance parameters are better than simple mean prediction
  naive = rep(mean(Y), length(Y))
  rmse_naive = sqrt(mean((Y - naive)^2))
  rmse_standard = sqrt(mean((Y - np_standard[["nuisance_parameters"]][["Y.hat"]])^2))
  rmse_short = sqrt(mean((Y - np_short[["nuisance_parameters"]][["Y.hat"]])^2))
  rmse_standard_w = sqrt(mean((Y - np_standard_w[["nuisance_parameters"]][["Y.hat"]])^2))
  rmse_short_w = sqrt(mean((Y - np_short_w[["nuisance_parameters"]][["Y.hat"]])^2))

  expect_true(rmse_naive > rmse_standard)
  expect_true(rmse_naive > rmse_short)
  expect_true(rmse_naive > rmse_standard_w)
  expect_true(rmse_naive > rmse_short_w)

})


test_that("check nuisance_parameters for D.hat", {

  library(mvtnorm)
  set.seed(234)

  n = 1000
  p = 6
  cov_mat = toeplitz(0.7^(0:(6 - 1)))

  X = mvtnorm::rmvnorm(n = n, mean = rep(2, p), sigma = cov_mat)
  colnames(X) <- paste0("X", 1:6)
  
  x1 = X[, 1]; x2 = X[, 2]; x3 = X[, 3]; x4 = X[, 4]; x5 = X[, 5]; x6 = X[, 6];

  effect_class1 = sin(x1) + x2^2 + x3*x4
  effect_class2 = cos(x2) + x4^2 + x5*x6
  effect_class3 = x1^2 + exp(x3) + x2*x5
  effect_class1 = (effect_class1 - mean(effect_class1)) / sd(effect_class1)
  effect_class2 = (effect_class2 - mean(effect_class2)) / sd(effect_class2)
  effect_class3 = (effect_class3 - mean(effect_class3)) / sd(effect_class3)
  raw_scores = cbind(effect_class1, effect_class2, effect_class3)
  prob_classes = t(apply(raw_scores, 1, function(X) exp(X) / sum(exp(X))))

  D = apply(prob_classes, 1, function(p) sample(1:3, 1, prob = p))

  methods = list(
    "logit" = create_method("logit", multinomial = "one-vs-one"),
    "ranger" = create_method("ranger", multinomial = "one-vs-rest"),
    "logit_nnet" = create_method("logit_nnet", multinomial = "multiclass"),
    "prob_forest" = create_method("prob_forest", multinomial = "multiclass")
  )

  path = gsub("\\\\", "/", tempdir())
  unlink(paste0(path, "/*"))


  ### Short-Stacking ###
  expect_message(np_e_short <- nuisance_parameters(NuPa="D.hat", X = X, D = D, methods = methods, cf = 5, stacking = "short", 
                                                   storeModels = "No", path = path, quiet = FALSE), "Short-stacking is used.")
  
  ## Saving models for non-outcome NuPas is disabled
  
  # files = paste0(path, "/Ensemble_W", 1:3, ".rds")
  # files_w = paste0(path, "/Ensemble_W", 1:3, "_Weights.rds")
  # 
  # # check if ensemble output is stored in files
  # expect_true(all(file.exists(files)))
  # expect_false(any(file.exists(files_w)))
  # 
  # # check if ensemble output contains correct info
  # ens = readRDS(files[1])
  # expect_length(ens, 2)
  # expect_named(ens$nnls_w)
  # expect_equal(dim(ens$cf_preds), c(n, length(methods)))
  # expect_true(all(substr(colnames(ens$cf_preds), 1, 1) == "t"))
  # unlink(paste0(path, "/*"))

  # check for correct dimension
  expect_identical(dim(np_e_short$nuisance_parameters$D.hat), dim(np_e_short$numbers$d_mat))

  # check if probabilities sum up to 1 for every obs
  expect_equal(rep(1, n), rowSums(np_e_short$nuisance_parameters$D.hat), tolerance = 1e-9)

  # check that nuisance parameters are better than simple mean prediction
  D_indicator <- model.matrix(~ factor(D) - 1)
  K <- ncol(D_indicator)  # Number of classes
  brier_naive <- (K - 1)/K
  brier_np_short <- mean((np_e_short$nuisance_parameters$D.hat - D_indicator)^2)
  expect_true(brier_naive > brier_np_short)


  ### Standard-Stacking ###
  expect_message(np_e_standard <- nuisance_parameters(NuPa="D.hat", X = X, D = D, methods = methods, cf = 5, stacking = 3, 
                                                      storeModels = "No", path = path, quiet = FALSE), "Standard-stacking is used.")

  ## Saving models for non-outcome NuPas is disabled
  
  # files = paste0(path, "/Ensemble_W", 1:3, ".rds")
  # files_w = paste0(path, "/Ensemble_W", 1:3, "_Weights.rds")
  # 
  # # check if ensemble output is stored in files
  # expect_true(all(file.exists(files)))
  # expect_false(any(file.exists(files_w)))
  # 
  # # check if ensemble output contains correct info
  # ens = readRDS(files[1])
  # expect_length(ens, ncol(cf_mat))
  # expect_named(ens[[1]]$nnls_w)
  # expect_equal(dim(ens[[1]]$cf_preds), c(sum(cf_mat[, 1]), length(methods)))
  # expect_true(all(substr(colnames(ens[[1]]$cf_preds), 1, 1) == "t"))
  # unlink(paste0(path, "/*"))
  

  # check for correct dimension
  expect_identical(dim(np_e_standard$nuisance_parameters$D.hat), dim(np_e_standard$numbers$d_mat))

  # check if probabilities sum up to 1 for every obs
  expect_equal(rep(1, n), rowSums(np_e_standard$nuisance_parameters$D.hat), tolerance = 1e-9)
  
  # check that nuisance parameters are better than simple mean prediction
  brier_np_stand <- mean((np_e_standard$nuisance_parameters$D.hat - D_indicator)^2)
  expect_true(brier_naive > brier_np_stand)

})


test_that("check nuisance_parameters for Y.hat.d estimation", {
  
  set.seed(79)

  n = 1000
  p = 8
  p_act = 4

  pi = c(seq(1, 0.1, -(1 / p_act)) * rep(c(1, -1), p_act / 2), rep(0, p - p_act))
  cov_mat = toeplitz(0.7^(0:(p - 1)))

  X = mvtnorm::rmvnorm(n = n, mean = rep(0, p), sigma = cov_mat)
  Y = as.vector(X %*% pi + rnorm(n, 0, 1))

  d_mods = 3
  D = sample(1:d_mods, n, replace = TRUE)


  methods = list("ridge" = create_method("ridge"),
                "forest_grf" = create_method("forest_grf"),
                "rlasso" = create_method("rlasso"),
                "xgboost" = create_method("xgboost"))

  path = gsub("\\\\", "/", tempdir())


  
  ### Short-Stacking ###
  expect_message(np_m_short <- nuisance_parameters(NuPa=c("Y.hat.d"), X = X, D = D, Y = Y, methods = methods, cf = 4, stacking = "short", 
                                                   storeModels = "Disk", stratify = TRUE, path = path, quiet = FALSE), "Short-stacking is used.")


  # check if ensemble output is stored in files
  expect_true(file.exists(path))

  # check if ensemble output contains correct info
  complete_file = readRDS(file.path(path, "nuisance_models.RDS")) # open saved object
  ens = complete_file[["models"]][["Y.hat.d_m"]][[1]]             # select saved model for Y.hat.d[[1]] only
  expect_length(ens, 2)
  expect_named(ens$nnls_w)
  expect_equal(dim(ens[["ens_object"]][["cf_preds"]]), c(n, length(methods)))

  # check for correct dimension
  expect_identical(dim(np_m_short$nuisance_parameters$Y.hat.d), dim(np_m_short$numbers$d_mat))
  
  unlink(paste0(path, "/*"))


  ### Short-Stacking with Smoother Weights ###
  expect_message(np_m_short_w <- nuisance_parameters(NuPa=c("Y.hat.d"), X = X, D = D, Y = Y, methods = methods, cf = 4, stacking = "short", 
                                                     storeModels = "Disk", stratify = TRUE, path = path, quiet = FALSE), "Short-stacking is used.")


  # check if ensemble output is stored in files
  expect_true(file.exists(path))
  
  # check if ensemble output contains correct info
  complete_file = readRDS(file.path(path, "nuisance_models.RDS")) # open saved object
  ens = complete_file[["models"]][["Y.hat.d_m"]][[1]]             # select saved model for Y.hat.d[[1]] only
  expect_length(ens, 2)
  expect_named(ens$nnls_w)
  expect_equal(dim(ens[["ens_object"]][["cf_preds"]]), c(n, length(methods)))

  # check for correct dimension
  expect_identical(dim(np_m_short_w$nuisance_parameters$Y.hat.d), dim(np_m_short_w$numbers$d_mat))

  # check smoother weights
  w = get_outcome_weights(object = file.path(path, "nuisance_models.RDS"), NuPa = "Y.hat.d")$Y.hat.d[[1]]
  expect_equal(np_m_short_w[["nuisance_parameters"]][["Y.hat.d"]][, 1], as.vector(w %*% Y), tolerance = 1e-3)

  unlink(paste0(path, "/*"))


  
  ### Standard-Stacking ###
  expect_message(np_m_standard <- nuisance_parameters(NuPa=c("Y.hat.d"), X = X, D = D, Y = Y, methods = methods, cf = 4, stacking = 3, 
                                                      storeModels = "Disk", stratify = TRUE, path = path, quiet = FALSE), "Standard-stacking is used.")
  
  # check if ensemble output is stored in files
  expect_true(file.exists(path))
  
  # check if ensemble output contains correct info
  complete_file = readRDS(file.path(path, "nuisance_models.RDS")) # open saved object
  ens = complete_file[["models"]][["Y.hat.d_m"]][[1]]             # select saved model for Y.hat.d[1] only
  expect_length(ens, 4) # cf = 4
  expect_named(ens[[1]]$nnls_w)

  # check for correct dimension
  expect_identical(dim(np_m_standard$nuisance_parameters$Y.hat.d), dim(np_m_standard$numbers$d_mat))
  
  unlink(paste0(path, "/*"))


  
  ### Standard-Stacking with Smoother Weights ###
  expect_message(np_m_standard_w <- nuisance_parameters(NuPa=c("Y.hat.d"), X = X, D = D, Y = Y, methods = methods, cf = 4, stacking = 3, 
                                                      storeModels = "Disk", stratify = TRUE, path = path, quiet = FALSE), "Standard-stacking is used.")
  
  # check if ensemble output is stored in files
  expect_true(file.exists(path))
  
  # check if ensemble output contains correct info
  complete_file = readRDS(file.path(path, "nuisance_models.RDS")) # open saved object
  ens = complete_file[["models"]][["Y.hat.d_m"]][[1]]             # select saved model for Y.hat.d[[1]] only
  expect_length(ens, 4) # cf = 4
  expect_named(ens[[1]]$nnls_w)

  # check for correct dimension
  expect_identical(dim(np_m_standard_w$nuisance_parameters$Y.hat.d), dim(np_m_standard_w$numbers$d_mat))

  # check smoother weights
  w = get_outcome_weights(object = file.path(path, "nuisance_models.RDS"), NuPa = "Y.hat.d")$Y.hat.d[[1]]
  expect_equal(np_m_standard_w[["nuisance_parameters"]][["Y.hat.d"]][, 1], as.vector(w %*% Y), tolerance = 1e-3)
  
  unlink(paste0(path, "/*"))
  
  

  ## S-learner's functionality is temporarily suspended
  
  # ### S-Learner ###
  # 
  # expect_error(
  #   nuisance_m(methods = methods, Y = Y, w_mat = w_mat, X = X, cf_mat = cf_mat,
  #              cv = 3, path = path, quiet = FALSE, weights = TRUE,
  #              learner = "s"),
  #   "S-Learner cannot be combined")
  # 
  # expect_error(
  #   nuisance_m(methods = methods, Y = Y, w_mat = w_mat, X = X, cf_mat = cf_mat,
  #              cv = 1, path = path, quiet = FALSE, weights = FALSE,
  #              learner = "s"),
  #   "S-Learner cannot be combined")
  # 
  # expect_error(
  #   nuisance_m(methods = methods, Y = Y, w_mat = w_mat, X = X, cf_mat = cf_mat,
  #              cv = 5, path = path, quiet = FALSE, weights = TRUE,
  #              learner = "both"),
  #   "S-Learner cannot be combined")
  # 
  # expect_error(
  #   nuisance_m(methods = methods, Y = Y, w_mat = w_mat, X = X, cf_mat = cf_mat,
  #              cv = 2, path = path, quiet = FALSE, weights = FALSE,
  #              learner = "both"),
  #   "S-Learner cannot be combined")

})


