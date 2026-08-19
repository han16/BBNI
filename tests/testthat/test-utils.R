test_that("update_ancestor_matrix computes correct transitive closure", {
  # Chain: node 1 (root) -> node 2 -> node 3
  # incid_matrix[i, j] == 1 means node j is a parent of node i
  incid_matrix <- matrix(0, nrow = 3, ncol = 3)
  incid_matrix[2, 1] <- 1 # node 1 is a parent of node 2
  incid_matrix[3, 2] <- 1 # node 2 is a parent of node 3

  ances <- update_ancestor_matrix(incid_matrix)

  expect_equal(ances[2, 1], 1)
  expect_equal(ances[3, 2], 1)
  expect_equal(ances[3, 1], 1) # transitive: node 1 is an ancestor of node 3
  expect_equal(sum(diag(ances)), 0) # should be acyclic and have no self-ancestry
})

test_that("check_ances_matrix detects a directed cycle", {
  # 2-node cycle: 1 -> 2 -> 1
  incid_matrix <- matrix(c(0, 1, 1, 0), nrow = 2, byrow = TRUE)
  ances <- update_ancestor_matrix(incid_matrix)

  expect_gt(check_ances_matrix(ances), 0)
})

test_that("check_ances_matrix reports zero loops for a valid DAG", {
  incid_matrix <- matrix(0, nrow = 3, ncol = 3)
  incid_matrix[2, 1] <- 1
  incid_matrix[3, 2] <- 1
  ances <- update_ancestor_matrix(incid_matrix)

  expect_equal(check_ances_matrix(ances), 0)
})

test_that("GenerateNetwork produces valid DAGs with a maximum in-degree of 2", {
  # test different network sizes to ensure robustness
  for (n in c(5, 10, 20)) {
    set.seed(35) # ensure reproducibility
    net <- GenerateNetwork(num.node = n)
    # structural check: must be a square matrix of size n x n
    expect_equal(nrow(net), n, info = paste("Failed for n =", n))
    expect_equal(ncol(net), n, info = paste("Failed for n =", n))
    # DAG validity check: transitive closure must have zero on the diagonal
    incid <- (net > 0) * 1
    ances <- update_ancestor_matrix(incid)
    expect_equal(sum(diag(ances)), 0, info = paste("Network of size", n, "contains a cycle"))
    # In-degree constraint check: no node should have more than 2 parents
    in_degrees <- rowSums(incid)
    expect_true(all(in_degrees <= 2), info = paste("Network of size", n, "has a node with in-degree > 2"))
  }
})
