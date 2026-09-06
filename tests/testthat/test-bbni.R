# Helper to create small valid network and data
make_small_data <- function(n = 3, m = 6, timeseries = TRUE) {
  set.seed(123)
  net <- GenerateNetwork(num.node = n)
  data <- GenerateSample(
    trans_matrix = net,
    SampleSize = m,
    para = rep(0.5, n),
    error = matrix(0, nrow = n, ncol = m),
    timeseries = timeseries
  )
  list(net = net, data = data, n = n, m = m)
}

test_that("run_bbni executes & returns expected structure (timeseries)", {
  d <- make_small_data(3, 8, timeseries = TRUE)
  res <- run_bbni(
    GeneData = d$data,
    num.node = d$n,
    SampleSize = d$m,
    num_update = 2,
    verbose = FALSE
  )
  expect_type(res, "list")
  expect_named(res, c("networks", "log_posterior", "post_edge_prob", "burn_in"))
  expect_true(is.list(res$networks))
  expect_true(is.numeric(res$log_posterior))
  expect_true(is.matrix(res$post_edge_prob))
  expect_equal(dim(res$post_edge_prob), c(d$n, d$n))
  expect_equal(res$burn_in, 0.7)
})

test_that("run_bbni works w/ independent (non-timeseries) data", {
  d <- make_small_data(3, 8, timeseries = FALSE)
  res <- run_bbni(
    GeneData = d$data,
    num_update = 2,
    timeseries = FALSE
  )
  expect_true(is.matrix(res$post_edge_prob))
  expect_equal(dim(res$post_edge_prob), c(3, 3))
})

test_that("run_bbni handles custom priors, penalty, and prop.ratio", {
  d <- make_small_data(4, 6, timeseries = TRUE)
  prior <- matrix(3, nrow = 5, ncol = 2)
  prior[5, ] <- c(2, 100)
  res <- run_bbni(
    GeneData = d$data,
    prior_para = prior,
    penalty = 1,
    prop.ratio = 0.5,
    burn_in = 0.3,
    num_update = 2
  )
  expect_equal(res$burn_in, 0.3)
  expect_true(is.matrix(res$post_edge_prob))
})

test_that("run_bbni prints verbose output w/o error", {
  d <- make_small_data(3, 5, timeseries = TRUE)
  expect_output(
    res <- run_bbni(
      GeneData = d$data,
      num_update = 1,
      verbose = TRUE
    ),
    "Running BBNI MCMC Sampling"
  )
  expect_type(res, "list")
})

test_that("run_bbni validates penalty and prop.ratio", {
  d <- make_small_data(3, 5)
  expect_error(
    run_bbni(GeneData = d$data, penalty = 2),
    "penalty"
  )
  expect_error(
    run_bbni(GeneData = d$data, prop.ratio = -0.1),
    "prop.ratio"
  )
})

test_that("run_bbni retains rownames in post_edge_prob", {
  d <- make_small_data(3, 6)
  rownames(d$data) <- c("G1", "G2", "G3")
  res <- run_bbni(GeneData = d$data, num_update = 1)
  expect_equal(rownames(res$post_edge_prob), c("G1", "G2", "G3"))
  expect_equal(colnames(res$post_edge_prob), c("G1", "G2", "G3"))
})
