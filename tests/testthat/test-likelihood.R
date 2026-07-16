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
