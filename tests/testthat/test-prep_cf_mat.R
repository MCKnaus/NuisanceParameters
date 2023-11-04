test_that("prep_cf_mat uses check_cluster_compatibility correctly", {

  set.seed(3)

  # Test case 1: Imbalances with different # folds
  cl1 = c(1, 1, 1, 2, 2, 2, 3, 3, 3, 3)
  expect_no_error(prep_cf_mat(length(cl1), cf = 2, cl = cl1))
  expect_error(prep_cf_mat(length(cl1), cf = 3, cl = cl1), "There is a high imbalance in the cluster sizes.")

  # Test case 2: Too few clusters, should throw an error
  cl2 = c(1, 1, 2, 2, 2)
  expect_error(prep_cf_mat(length(cl2), cf = 3, cl = cl2), "The number of clusters is less than the desired number of folds.")

  # Test case 3: High imbalance in cluster sizes, should throw an error
  cl3 = c(1, 1, 1, 2, 2, 2, 2, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3)
  expect_error(prep_cf_mat(length(cl3), cf = 2, cl = cl3),  "There is a high imbalance in the cluster sizes.")
  expect_error(prep_cf_mat(length(cl3), cf = 3, cl = cl3),  "There is a high imbalance in the cluster sizes.")

})


test_that("prep_cf_mat returns correct dimensions", {

  set.seed(50)

  n = 1000
  cf = 5

  # Test without treatment matrix
  cf_mat1 = prep_cf_mat(n, cf)
  expect_equal(dim(cf_mat1), c(n, cf))

  # Test with treatment matrix
  treat = 3
  w_mat = prep_w_mat(sample(1:treat, n, replace = TRUE))
  cf_mat2 = prep_cf_mat(n, cf, w_mat = w_mat)
  expect_equal(dim(cf_mat2), c(n, cf))

  # Test with clusters
  cl = sample(x = c("a","b","c","d","e","f","g","h","i","j","k","l"), size = n, replace = TRUE, prob = c(rep(0.1,8), rep(0.05,4)))
  cf_mat3 = prep_cf_mat(n, cf, cl = cl)
  expect_equal(dim(cf_mat3), c(n, cf))

})


test_that("prep_cf_mat returns correct fold proportions with treatment matrix", {

  n = 1000
  cf = 5
  treatments = c("a","b","c","d")
  treatment_probs = c(0.5,0.4,0.08,0.02)
  w_mat = prep_w_mat(sample(x = treatments, size = n, replace = TRUE, prob = treatment_probs))

  cf_mat = prep_cf_mat(n, cf, w_mat = w_mat)

  # test fold sizes
  expect_equal(unname(colMeans(cf_mat)), rep(1/cf, cf), tolerance = 0.04)

  # test treatment proportions in different folds
  for (i in 1:cf){

    expect_equal(unname(colMeans(w_mat[cf_mat[,i],])), treatment_probs, tolerance = 0.1)

  }

})


test_that("prep_cf_mat returns warning if both treatment matrix and cluster vector is provided", {

  n = 1000
  cf = 5
  w_mat = prep_w_mat(sample(x = 1:3, size = n, replace = TRUE))
  cl = sample(x = c("b","d","o","k"), size = n, replace = TRUE)

  expect_warning(prep_cf_mat(n, cf, w_mat = w_mat, cl = cl), "You cannot provide both a treatment matrix and a cluster vector")

})


test_that("prep_cf_mat keeps observations from one cluster in the same cross-fitting fold", {

  set.seed(870)

  n = 10000
  cf = 10

  cl = sample(x = 1:100, size = n, replace = TRUE, prob = 100:1)

  cf_mat = prep_cf_mat(n, cf, cl = cl)

  # test fold sizes
  expect_equal(unname(colMeans(cf_mat)), rep(1/cf, cf), tolerance = 0.1)


  # test whether obs of a cluster are in the same fold
  collection = rep(NA, max(cl))
  for (i in 1:max(cl)) {
    cluster_indices = which(cl == i)
    fold_count = sum(unique(cf_mat[cluster_indices, ]))
    collection[i] = fold_count
  }

  expect_true(all(collection))

})
