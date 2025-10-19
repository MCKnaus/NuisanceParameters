library(testthat)

test_that("prep_indicator_mat works correctly", {
  # Nominal case
  x <- factor(c("A", "B", "A", "C", "B", "C"))
  mat <- prep_indicator_mat(x)
  expect_true(is.matrix(mat))
  expect_true(is.logical(mat))
  expect_equal(dim(mat), c(length(x), 3))
  expect_equal(colnames(mat), c("A", "B", "C"))
  expect_true(all(rowSums(mat) == 1))
  
  # Works with character input
  y <- c("dog", "cat", "dog", "mouse")
  mat2 <- prep_indicator_mat(y)
  expect_equal(dim(mat2), c(4, 3))
  expect_equal(sort(colnames(mat2)), sort(unique(y)))
  
  # Throws error if only one unique value
  expect_error(prep_indicator_mat(rep("a", 5)), "at least two unique")
})

test_that("check_cluster_compatibility detects issues", {
  set.seed(1)
  cluster <- rep(1:4, each = 5)
  expect_silent(check_cluster_compatibility(cluster, 2))
  expect_error(check_cluster_compatibility(cluster, 5), "less than the desired number of folds")
  
  # Highly imbalanced
  cluster <- c(rep(1, 17), 2, 3)
  expect_error(check_cluster_compatibility(cluster, 3), "high imbalance")
})

test_that("prep_cf_mat creates a valid fold matrix", {
  set.seed(123)
  N <- 50
  cf <- 3
  # No cluster or d_mat
  mat <- prep_cf_mat(N, cf)
  expect_true(is.matrix(mat))
  expect_equal(dim(mat), c(N, cf))
  expect_true(all(rowSums(mat) == 1))
  
  # d_mat preserves ratios
  d_vec <- sample(c("A", "B", "C"), N, replace = TRUE)
  dmat <- prep_indicator_mat(d_vec)
  mat2 <- prep_cf_mat(N, cf, d_mat = dmat)
  expect_equal(dim(mat2), c(N, cf))
  expect_true(all(rowSums(mat2) == 1))
})

test_that("ens_weights_maker works for continuous outcomes", {
  set.seed(123)
  X <- matrix(c(0.2, 0.3, 0.4,
                0.5, 0.6, 0.7), ncol = 2,
              dimnames = list(NULL, c("learner1", "learner2")))
  Y <- c(0.25, 0.45, 0.65)
  
  w <- ens_weights_maker(X = X, Y = Y)
  
  expect_type(w, "double")
  expect_length(w, ncol(X))
  expect_named(w, colnames(X))
  expect_true(all(w >= 0))
  expect_equal(sum(w), 1, tolerance = 1e-8)
})

test_that("ens_weights_maker handles subset argument", {
  set.seed(123)
  X <- matrix(runif(20), ncol = 2,
              dimnames = list(NULL, c("l1", "l2")))
  Y <- rnorm(10)
  
  w_full <- ens_weights_maker(X, Y)
  w_sub <- ens_weights_maker(X, Y, subset = rep(c(TRUE, FALSE), 5))
  
  expect_true(all(w_full >= 0))
  expect_true(all(w_sub >= 0))
  expect_equal(sum(w_full), 1, tolerance = 1e-8)
  expect_equal(sum(w_sub), 1, tolerance = 1e-8)
})

test_that("ens_weights_maker works for multinomial outcomes with stacked NNLS", {
  set.seed(123)
  N <- 5; K <- 3; M <- 2
  X <- array(runif(N * K * M), dim = c(N, K, M),
             dimnames = list(NULL, paste0("class", 1:K), paste0("learner", 1:M)))
  Y <- factor(sample(1:K, N, replace = TRUE))
  
  w <- ens_weights_maker(X = X, Y = Y, is_mult = TRUE, ensemble_type = "nnls")
  
  expect_type(w, "double")
  expect_length(w, M)
  expect_named(w, dimnames(X)[[3]])
  expect_true(all(w >= 0))
  expect_equal(sum(w), 1, tolerance = 1e-8)
})

test_that("ens_weights_maker works for multinomial outcomes with BFGS optimization", {
  skip_if_not_installed("stats")
  
  set.seed(123)
  N <- 6; K <- 3; M <- 2
  X <- array(runif(N * K * M), dim = c(N, K, M),
             dimnames = list(NULL, paste0("class", 1:K), paste0("learner", 1:M)))
  Y <- factor(sample(1:K, N, replace = TRUE))
  
  w <- ens_weights_maker(X = X, Y = Y, is_mult = TRUE, ensemble_type = "bfgs")
  
  expect_type(w, "double")
  expect_length(w, M)
  expect_named(w, dimnames(X)[[3]])
  expect_true(all(w >= 0))
  expect_equal(sum(w), 1, tolerance = 1e-8)
})

test_that("ens_weights_maker falls back to uniform weights if solution is degenerate", {
  # Construct degenerate case where NNLS returns zero weights
  X <- matrix(0, nrow = 4, ncol = 3,
              dimnames = list(NULL, c("a", "b", "c")))
  Y <- rep(0, 4)
  
  w <- ens_weights_maker(X, Y)
  
  expect_equal(as.numeric(w), rep(1/3, 3))
  expect_equal(sum(w), 1)
})

test_that("add_intercept works as expected", {
  m <- matrix(1:6, ncol = 2)
  m2 <- add_intercept(m)
  expect_equal(ncol(m2), 3)
  expect_true(all(m2[, 1] == 1))
  # If already has intercept
  m3 <- cbind(1, m)
  m4 <- add_intercept(m3)
  expect_equal(m3, m4)
})

test_that("one_hot returns correct matrix", {
  y <- c("dog", "cat", "dog", "mouse")
  mat <- one_hot(y)
  expect_equal(dim(mat), c(4, 3))
  expect_equal(colnames(mat), unique(y))
  # Each row only one 1
  expect_true(all(rowSums(mat) == 1))
  expect_equal(as.numeric(mat[1, "dog"]), 1)
  expect_equal(as.numeric(mat[2, "cat"]), 1)
})

test_that("create_method returns expected structure for defaults", {
  m <- create_method("ols")
  
  expect_type(m, "list")
  expect_equal(m$method, "ols")
  expect_null(m$multinomial)
  expect_false(m$parallel)
  expect_equal(m$arguments, list())
  expect_null(m$x_select)
})

test_that("create_method supports list and character arguments", {
  m1 <- create_method("knn", arguments = list(k = 3))
  m2 <- create_method("forest_grf", arguments = list("num.trees" = 500))
  
  expect_equal(m1$arguments$k, 3)
  expect_equal(m2$arguments$num.trees, 500)
})

test_that("create_method validates x_select", {
  m <- create_method("ridge", x_select = c(TRUE, FALSE, TRUE))
  expect_equal(m$x_select, c(TRUE, FALSE, TRUE))
  
  expect_error(
    create_method("ridge", x_select = 1:3),
    "Provide either NULL or logical for x_select."
  )
})

test_that("create_method validates parallel flag", {
  m <- create_method("plasso", parallel = TRUE)
  expect_true(m$parallel)
  
  expect_error(
    create_method("plasso", parallel = c(TRUE, FALSE)),
    "parallel must be a single logical value"
  )
})

test_that("create_method handles multinomial options", {
  m1 <- create_method("logit", multinomial = "one-vs-rest")
  m2 <- create_method("prob_forest", multinomial = "multiclass")
  
  expect_equal(m1$multinomial, "one-vs-rest")
  expect_equal(m2$multinomial, "multiclass")
})

test_that("create_method errors on unsupported multinomial with svm", {
  expect_error(
    create_method("svm", multinomial = "multiclass"),
    "SVM does not support multiclass estimation"
  )
})

# helper to quickly create fake methods
fake_method <- function(name) {
  list(method = name, arguments = list(), multinomial = NULL,
       parallel = FALSE, x_select = NULL)
}

test_that("process_methods replicates a simple methods list across NuPas", {
  methods <- list(fake_method("ols"), fake_method("ridge"))
  NuPa <- c("Y.hat", "D.hat")
  
  result <- process_methods(methods, NuPa, K = 2)
  
  expect_type(result, "list")
  expect_named(result, NuPa)
  expect_equal(as.numeric(lengths(result)), c(2, 2))
  expect_equal(result$Y.hat[[1]]$method, "ols")
  expect_equal(result$D.hat[[2]]$method, "ridge")
})

test_that("process_methods applies NuPa-specific methods", {
  methods <- list(
    Y.hat = list(fake_method("ols")),
    D.hat = list(fake_method("ridge"))
  )
  NuPa <- c("Y.hat", "D.hat")
  
  result <- process_methods(methods, NuPa, K = 2)
  
  expect_equal(result$Y.hat[[1]]$method, "ols")
  expect_equal(result$D.hat[[1]]$method, "ridge")
})

test_that("process_methods filters out multiclass methods for Y.hat NuPa", {
  methods <- list(fake_method("ols"), fake_method("svm")) # svm is multiclass
  NuPa <- "Y.hat"
  
  expect_message(
    result <- process_methods(methods, NuPa, K = 2),
    "Multiclass methods incompatible with NuPa"
  )
  
  expect_equal(length(result$Y.hat), 1)
  expect_equal(result$Y.hat[[1]]$method, "ols")
})

test_that("process_methods filters out base methods for D.hat when K > 2", {
  methods <- list(fake_method("ols"), fake_method("ridge"))
  NuPa <- "D.hat"
  
  expect_error(
    process_methods(methods, NuPa, K = 3),
    "The following nuisance parameters have empty methods lists"
  )
})

test_that("process_methods allows base methods for D.hat when K = 2", {
  methods <- list(fake_method("ols"), fake_method("ridge"))
  NuPa <- "D.hat"
  
  result <- process_methods(methods, NuPa, K = 2)
  
  expect_equal(length(result$D.hat), 2)
})

test_that("process_methods NuPa-specific filtering works", {
  methods <- list(
    Y.hat = list(fake_method("svm"), fake_method("ols")),
    D.hat = list(fake_method("ridge"), fake_method("logit"))
  )
  NuPa <- c("Y.hat", "D.hat")
  
  expect_message(
    result <- process_methods(methods, NuPa, K = 3),
    "Multiclass methods incompatible with NuPa"
  )
  expect_equal(result$Y.hat[[1]]$method, "ols")
})

test_that("process_methods errors if any NuPa has no valid methods", {
  methods <- list(fake_method("svm")) # invalid for Y.hat
  NuPa <- "Y.hat"
  
  expect_error(
    process_methods(methods, NuPa, K = 2),
    "The following nuisance parameters have empty methods lists"
  )
})
