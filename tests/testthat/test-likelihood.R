test_that("Error_LLH handles an all-root-node network correctly", {
  # A network with no edges: every node is a root node, exercising the
  # branch of Error_LLH that should report no error/mismatch terms.
  num.node <- 4
  SampleSize <- 50
  GeneData <- matrix(rbinom(num.node * SampleSize, 1, 0.5),
    nrow = num.node, ncol = SampleSize
  )
  TRFUM <- matrix(0, nrow = num.node, ncol = num.node)
  prior_para <- matrix(3, nrow = num.node + 1, ncol = 2)
  prior_para[num.node + 1, 1] <- 2
  prior_para[num.node + 1, 2] <- 100

  result <- Error_LLH(
    TRFUM = TRFUM, GeneData = GeneData, SampleSize = SampleSize,
    num.node = num.node, prior_para = prior_para, penalty = 0.1
  )

  expect_true(is.na(result[[1]][1])) # ErrorFactor: NA, no non-root nodes
  expect_true(is.na(result[[2]][num.node + 1])) # mismatch: NA
  expect_true(is.na(result[[2]][num.node + 2])) # Perror: NA

  log_post_model <- result[[1]][length(result[[1]])]
  expect_true(is.finite(log_post_model))
})

test_that("Error_LLH precisely shifts mismatch calculation for timeseries vs independent mode", {
  num.node <- 2
  SampleSize <- 4
  # Deterministic GeneData matrix:
  GeneData <- matrix(c(
    1, 0, 1, 0,
    0, 1, 0, 1
  ), nrow = num.node, byrow = TRUE)
  # Node 1 is the parent of Node 2 (Rule 11 = Identity).
  TRFUM <- matrix(0, nrow = num.node, ncol = num.node)
  TRFUM[2, 1] <- 11
  # Standard Beta priors
  prior_para <- matrix(3, nrow = num.node + 1, ncol = 2)
  prior_para[num.node + 1, 1] <- 2
  prior_para[num.node + 1, 2] <- 100
  # Test Time-Series Mode
  result_ts <- Error_LLH(
    TRFUM = TRFUM, GeneData = GeneData, SampleSize = SampleSize,
    num.node = num.node, prior_para = prior_para, penalty = 0.1,
    timeseries = TRUE
  )
  # Expect 0 actual mismatches + 0.0001 pseudo-count
  mismatch_ts <- result_ts[[2]][num.node + 1]
  expect_equal(mismatch_ts, 0.0001)
  # Test Independent Mode
  result_static <- Error_LLH(
    TRFUM = TRFUM, GeneData = GeneData, SampleSize = SampleSize,
    num.node = num.node, prior_para = prior_para, penalty = 0.1,
    timeseries = FALSE
  )
  # Expect 4 actual mismatches + 0.0001 pseudo-count
  mismatch_static <- result_static[[2]][num.node + 1]
  expect_equal(mismatch_static, 4.0001)
})

test_that("penalty = 1 produces uniform prior (no edge-count structural penalty)", {
  num.node <- 3
  SampleSize <- 10
  set.seed(42)
  GeneData <- matrix(rbinom(num.node * SampleSize, 1, 0.5),
                     nrow = num.node, ncol = SampleSize
  )
  prior_para <- matrix(3, nrow = num.node + 1, ncol = 2)
  prior_para[num.node + 1, 1] <- 2
  prior_para[num.node + 1, 2] <- 100
  # topology w/ 1 edge: node 2 has parent node 1 (identity function)
  TRFUM_1edge <- matrix(0, nrow = num.node, ncol = num.node)
  TRFUM_1edge[2, 1] <- 11
  # topology w/ 2 edges: node 2 has parent 1, node 3 has parent 1
  TRFUM_2edge <- matrix(0, nrow = num.node, ncol = num.node)
  TRFUM_2edge[2, 1] <- 11
  TRFUM_2edge[3, 1] <- 11
  # helper to get log_post_model (last element of first returned vector)
  get_logpost <- function(TRFUM, penalty) {
    res <- Error_LLH(
      TRFUM = TRFUM, GeneData = GeneData, SampleSize = SampleSize,
      num.node = num.node, prior_para = prior_para, penalty = penalty
    )
    res[[1]][length(res[[1]])]
  }
  # under penalty = 1, the structural term is d(T) * log(1) = 0 for all topologies
  # under penalty = 0.5, the structural term is d(T) * log(0.5)
  # therefore, logpost(pen=0.5) - logpost(pen=1) should equal d(T) * log(0.5)
  expect_equal(get_logpost(TRFUM_1edge, 0.5) - get_logpost(TRFUM_1edge, 1), 1 * log(0.5))
  expect_equal(get_logpost(TRFUM_2edge, 0.5) - get_logpost(TRFUM_2edge, 1), 2 * log(0.5))
  # under penalty = 1, the difference in log-posterior between two topologies
  # is purely from likelihood. the structural term contributes 0 to the difference.
  # under penalty = 0.5, the difference gains an extra (d2 - d1) * log(0.5) term.
  diff_unif <- get_logpost(TRFUM_2edge, 1) - get_logpost(TRFUM_1edge, 1)
  diff_pen  <- get_logpost(TRFUM_2edge, 0.5) - get_logpost(TRFUM_1edge, 0.5)
  expect_equal(diff_pen - diff_unif, (2 - 1) * log(0.5))
  # under penalty = 1, two topologies with the SAME edge count but different
  # edges should differ only by likelihood (structural term is identical: 0).
  TRFUM_2edge_alt <- matrix(0, nrow = num.node, ncol = num.node)
  TRFUM_2edge_alt[3, 2] <- 11
  TRFUM_2edge_alt[3, 1] <- 11
  diff_same_edge_count <- get_logpost(TRFUM_2edge, 1) - get_logpost(TRFUM_2edge_alt, 1)
  # structural terms cancel, so difference is finite and purely likelihood-based
  expect_true(is.finite(diff_same_edge_count))
})
