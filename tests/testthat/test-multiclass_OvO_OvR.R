set.seed(123)
n <- 50
X <- data.frame(x1 = rnorm(n), x2 = rnorm(n))
Y_3class <- sample(1:3, size = n, replace = TRUE)
Y_0_2 <- sample(0:2, size = n, replace = TRUE)
Xnew <- data.frame(x1 = rnorm(10), x2 = rnorm(10))

# --- OvO Tests ---

test_that("ovo_fit returns correct number of classifiers for 3 classes", {
  ovo <- ovo_fit(X, Y_3class, method = "logit", parallel = FALSE, verbose = FALSE)
  expect_type(ovo, "list")
  expect_equal(length(ovo), choose(length(unique(Y_3class)), 2))
  expect_true(all(grepl("_", names(ovo))))
})

test_that("ovo_fit works when classes start at 0", {
  ovo <- ovo_fit(X, Y_0_2, method = "logit", parallel = FALSE, verbose = FALSE)
  expect_equal(length(ovo), choose(length(unique(Y_0_2)), 2))
})

test_that("predict.ovo_fit returns matrix with correct dimensions and probabilities", {
  ovo <- ovo_fit(X, Y_3class, method = "logit", parallel = FALSE, verbose = FALSE)
  preds <- predict.ovo_fit(ovo, X, Y_3class, Xnew = Xnew, method = "logit", parallel = FALSE, verbose = FALSE)
  expect_type(preds, "double")
  expect_true(is.matrix(preds) | is.data.frame(preds))
  expect_equal(dim(preds), c(nrow(Xnew), length(unique(Y_3class))))
  expect_equal(round(rowSums(preds), 5), rep(1, nrow(Xnew)))
  expect_true(all(preds >= 0 & preds <= 1))
})

test_that("predict.ovo_fit works with training data as Xnew=NULL", {
  ovo <- ovo_fit(X, Y_3class, method = "logit", parallel = FALSE, verbose = FALSE)
  preds <- predict.ovo_fit(ovo, X, Y_3class, Xnew = NULL, method = "logit", parallel = FALSE, verbose = FALSE)
  expect_equal(nrow(preds), nrow(X))
})

# --- OvR Tests ---

test_that("ovr_fit returns correct number of classifiers for 3 classes", {
  ovr <- ovr_fit(X, Y_3class, method = "logit")
  expect_type(ovr, "list")
  expect_equal(length(ovr), length(unique(Y_3class)))
})

test_that("ovr_fit works if minimum class is 0", {
  ovr <- ovr_fit(X, Y_0_2, method = "logit")
  expect_equal(length(ovr), length(unique(Y_0_2)))
})

test_that("predict.ovr_fit returns matrix with correct dimensions and probabilities", {
  ovr <- ovr_fit(X, Y_3class, method = "logit")
  preds <- predict.ovr_fit(ovr, X, Y_3class, Xnew = Xnew, method = "logit")
  expect_true(is.data.frame(preds) | is.matrix(preds))
  expect_equal(dim(preds), c(nrow(Xnew), length(unique(Y_3class))))
  expect_equal(round(rowSums(preds), 5), rep(1, nrow(Xnew)))
  expect_true(all(preds >= 0 & preds <= 1))
})

test_that("predict.ovr_fit works with training data as Xnew=NULL", {
  ovr <- ovr_fit(X, Y_3class, method = "logit")
  preds <- predict.ovr_fit(ovr, X, Y_3class, Xnew = NULL, method = "logit")
  expect_equal(nrow(preds), nrow(X))
})

# --- KL Divergence Test ---

test_that("kl_convergence returns numeric and finite values", {
  p <- c(0.2, 0.5, 0.3)
  q_matrix <- matrix(runif(9, 0, 1), nrow = 3)
  diag(q_matrix) <- NA  # ignore diagonal as in code
  val <- kl_convergence(p, q_matrix)
  expect_type(val, "double")
  expect_true(is.finite(val))
})

test_that("kl_convergence is zero when p matches q_matrix exactly for pairs", {
  p <- c(0.5, 0.3, 0.2)
  q_matrix <- matrix(NA, 3, 3)
  q_matrix[1, 2] <- p[1] / (p[1] + p[2])
  q_matrix[1, 3] <- p[1] / (p[1] + p[3])
  q_matrix[2, 3] <- p[2] / (p[2] + p[3])
  q_matrix[2, 1] <- 1 - q_matrix[1, 2]
  q_matrix[3, 1] <- 1 - q_matrix[1, 3]
  q_matrix[3, 2] <- 1 - q_matrix[2, 3]
  val <- kl_convergence(p, q_matrix)
  expect_true(abs(val) < 1e-6)
})

