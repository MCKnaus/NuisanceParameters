test_that("prep_w_mat function works correctly for mickey mouse example", {
  # Example treatment vector
  w = c("A", "B", "A", "C", "B", "C")

  # Call the function
  w_mat = prep_w_mat(w)

  # Expected one-hot matrix
  expected_mat = matrix(c(
    TRUE, FALSE, FALSE,
    FALSE, TRUE, FALSE,
    TRUE, FALSE, FALSE,
    FALSE, FALSE, TRUE,
    FALSE, TRUE, FALSE,
    FALSE, FALSE, TRUE
  ), ncol = 3, byrow = TRUE)
  colnames(expected_mat) = unique(w)
  row.names(expected_mat) = as.character(1:length(w))

  # Check the result
  expect_identical(w_mat, expected_mat)
})


test_that("prep_w_mat function works correctly for large n", {

  n = 10000
  w = sample(x = 1:4, size = n, replace = TRUE, prob = c(0.5, 0.1, 0.05, 0.35))

  w_mat = prep_w_mat(w)

  expect = as.vector(table(w) / n)
  names(expect) = as.character(1:4)

  # Check class frequencies
  expect_identical(expect, colMeans(w_mat))
  # Check one-hot encoding
  expect_true(all(Matrix::rowSums(w_mat) == 1))

})
