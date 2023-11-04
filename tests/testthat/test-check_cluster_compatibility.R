test_that("check_cluster_compatibility function with small examples", {

  # Test case 1: Imbalances with different # folds
  cl1 = c(1, 1, 2, 2, 2, 3, 3, 3, 3)
  cf1a = 3
  cf1b = 2
  expect_error(check_cluster_compatibility(cl1, cf1a))
  expect_no_error(check_cluster_compatibility(cl1, cf1b))

  # Test case 2: Too few clusters, should throw an error
  cl2 = c(1, 1, 2, 2, 2)
  cf2 = 3
  expect_error(check_cluster_compatibility(cl2, cf2), "The number of clusters is less than the desired number of folds.")

  # Test case 3: High imbalance in cluster sizes, should throw an error
  cl3 = c(1, 1, 1, 2, 2, 2, 2, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3)
  expect_error(check_cluster_compatibility(cl3, 2),  "There is a high imbalance in the cluster sizes.")
  expect_error(check_cluster_compatibility(cl3, 3),  "There is a high imbalance in the cluster sizes.")

})

# Define unit tests
test_that("check_cluster_compatibility function with larger cluster vector", {

  set.seed(72)
  n = 10000

  # Dealing with imbalances
  cl1 = sample(x = 1:4, size = n, replace = TRUE, prob = c(0.5,0.1,0.05,0.35))
  expect_error(check_cluster_compatibility(cl1, 2), "There is a high imbalance in the cluster sizes.")
  expect_error(check_cluster_compatibility(cl1, 3), "There is a high imbalance in the cluster sizes.")

  cl2 = sample(x = 1:21, size = n, replace = TRUE, prob = c(0.3, rep(0.035,20)))
  expect_no_error(check_cluster_compatibility(cl2, 2))
  expect_error(check_cluster_compatibility(cl2, 3), "There is a high imbalance in the cluster sizes.")

  cl3 = sample(x = 1:25, size = n, replace = TRUE, prob = c(0.15,0.15,0.1,0.05,0.05,rep(0.025,20)))
  expect_no_error(check_cluster_compatibility(cl3, 2))
  expect_no_error(check_cluster_compatibility(cl3, 5))
  expect_error(check_cluster_compatibility(cl3, 10), "There is a high imbalance in the cluster sizes.")


  # Too few clusters
  cl4 = sample(x = 1:5, size = n, replace = TRUE, prob = c(0.15,0.15,0.15,0.3,0.25))
  expect_no_error(check_cluster_compatibility(cl4, 3))
  expect_error(check_cluster_compatibility(cl4, 10), "The number of clusters is less than the desired number of folds.")

})
