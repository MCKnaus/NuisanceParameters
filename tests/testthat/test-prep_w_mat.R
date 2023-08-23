test_that("prep_w_mat function works correctly", {
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
