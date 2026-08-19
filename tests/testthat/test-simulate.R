test_that("GenerateSample applies noise exactly via bitwise XOR", {
  # Very small 2-node network: Node 1 is root, Node 2 is exactly Node 1 (Rule 11)
  net <- matrix(0, nrow = 2, ncol = 2)
  net[2, 1] <- 11
  sample_size <- 10
  para <- c(0.5, 0.5)
  # generate perfectly clean data (no noise)
  clean_error <- matrix(0, nrow = 2, ncol = sample_size)
  set.seed(38) # Set seed for reproducibility
  clean_data <- GenerateSample(
    trans_matrix = net, SampleSize = sample_size, para = para, error = clean_error, timeseries = FALSE
  )
  # generate data with an all-1s error matrix (should flip every single bit)
  flip_error <- matrix(1, nrow = 2, ncol = sample_size)
  set.seed(38) # Must match the seed above so root nodes are identical before flipping
  flipped_data <- GenerateSample(
    trans_matrix = net, SampleSize = sample_size, para = para, error = flip_error, timeseries = FALSE
  )
  # assert that the flipped data is exactly the bitwise XOR of the clean data and 1
  # bitwXor removes dimensions in R, so wrap in matrix() to match output shape
  expected_flipped <- matrix(bitwXor(as.integer(clean_data), 1L), nrow = nrow(clean_data), ncol = ncol(clean_data))
  # root nodes are generated w/ rbinom and do not have the error matrix applied to them
  expected_flipped[1, ] <- clean_data[1, ]
  expect_equal(flipped_data, expected_flipped)
})
