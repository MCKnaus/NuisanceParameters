## Tests for design_matrix() function ----------------------------

test_that("design_matrix creates basic matrix correctly", {
  data <- matrix(rnorm(100), ncol = 2)
  colnames(data) <- c("x1", "x2")
  
  # Test basic functionality
  result <- design_matrix(data)
  expect_equal(ncol(result), 2) # Should return same columns if no transformations
  expect_equal(colnames(result), c("x1", "x2"))
  expect_equal(nrow(result), 50)
})

test_that("design_matrix handles interactions correctly", {
  data <- matrix(1:6, ncol = 3)
  colnames(data) <- c("x1", "x2", "x3")
  
  # Test interactions
  result <- design_matrix(data, int = c("x1", "x2", "x3"), int_d = 2)
  expect_equal(ncol(result), 6) # "x1","x2","x3","x1:x2","x1:x3","x2:x3"
  expect_true("x1:x2" %in% colnames(result))
  expect_true("x2:x3" %in% colnames(result))
  expect_true("x1:x3" %in% colnames(result))
  
  # Test higher order interactions
  result <- design_matrix(data, int = c("x1", "x2", "x3"), int_d = 3)
  expect_true("x1:x2:x3" %in% colnames(result))
})

test_that("design_matrix handles polynomials correctly", {
  data <- matrix(1:4, ncol = 1)
  colnames(data) <- c("x1")
  
  result <- design_matrix(data, poly = "x1", poly_d = 2)
  expect_equal(ncol(result), 2) # x1, x1^2
  expect_true(all(c("x11", "x12") %in% colnames(result)))
})

test_that("design_matrix handles logs correctly", {
  data <- matrix(1:4, ncol = 1)
  colnames(data) <- c("x1")
  
  result <- design_matrix(data, log = "x1")
  expect_equal(ncol(result), 2) # x1, ln_x1
  expect_true("ln_x1" %in% colnames(result))
  
  # Test non-positive values handling
  data[,1] <- c(-1, 0, 1, 2)
  
  result <- design_matrix(data, log = "x1")
  expect_equal(ncol(result), 1)
  expect_equal(colnames(result), "x1")
})

test_that("design_matrix combines all transformations correctly", {
  data <- matrix(1:9, ncol = 3)
  colnames(data) <- c("x1", "x2", "x3")
  
  result <- design_matrix(data, 
                          int = c("x1", "x2"), 
                          poly = "x3", 
                          log = "x1")
  
  # Check all transformation types are present
  expect_true(any(grepl(":", colnames(result)))) # Interaction
  expect_true(any(grepl("\\d$", colnames(result)))) # Polynomial
  expect_true(any(grepl("^ln_", colnames(result)))) # Log
})

## Tests for data_screen() function ----------------------------

test_that("data_screen removes invariant variables", {
  data <- matrix(c(1,1,1,1, 1,2,3,4), ncol = 2)
  colnames(data) <- c("const", "var")
  
  result <- data_screen(data, quiet = TRUE)
  expect_equal(ncol(result), 1)
  expect_equal(colnames(result), "var")
})

test_that("data_screen handles binary variables correctly", {
  data <- matrix(c(1,1,0,0, 1,1,1,1, 0,0,0,1), ncol = 3)
  colnames(data) <- c("balanced", "all1", "rare1")
  
  # Test without treatment groups
  result <- data_screen(data, bin_cut = 0.2, quiet = TRUE)
  expect_equal(ncol(result), 2)
  expect_true("balanced" %in% colnames(result)) # Should be kept with default cutoff
  
  # # Test with stricter cutoff
  # result <- data_screen(data, bin_cut = 0.3, quiet = TRUE)
  # expect_equal(ncol(result), 1)
})

## Suspended until further decision

# test_that("data_screen handles treatment groups correctly", {
#   data <- matrix(c(1,1,1,0,0,0, 1,0,1,0,1,0), ncol = 2)
#   colnames(data) <- c("group1", "group2")
#   treat <- c(1,1,1,0,0,0) # First 3 obs in treatment
#   
#   # Test with treatment groups
#   result <- data_screen(data, treat = treat, bin_cut = 0.4, quiet = TRUE)
#   expect_equal(ncol(result), 1) # group2 should be kept (50% in each group)
# })
# 
# test_that("data_screen removes correlated variables", {
#   data <- matrix(c(1,2,3,4, 2,4,6,8, 1,3,5,7), ncol = 3)
#   colnames(data) <- c("x1", "x2", "x3")
#   
#   result <- data_screen(data, corr_cut = 0.9, quiet = TRUE)
#   expect_equal(ncol(result), 2) # Should remove one correlated variable
#   expect_false(all(c("x1", "x2") %in% colnames(result))) # Either x1 or x2 removed
# })
# 
# test_that("data_screen verbose mode works", {
#   data <- matrix(c(1,1,1,1, 1,2,3,4), ncol = 2)
#   colnames(data) <- c("const", "var")
#   
#   expect_message(data_screen(data, quiet = FALSE), 
#                  "Variables with no variation")
# })
# 
# test_that("data_screen handles edge cases", {
#   # Empty input
#   expect_error(data_screen(matrix(ncol = 0, nrow = 0)))
#   
#   # Single column input
#   result <- data_screen(matrix(1:4, ncol = 1))
#   expect_equal(ncol(result), 1)
#   
#   # All invariant columns
#   data <- matrix(1, nrow = 4, ncol = 3)
#   expect_warning(data_screen(data))
# })
