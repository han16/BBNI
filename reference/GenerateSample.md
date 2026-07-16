# Simulate Boolean Network Observation Dataset

Simulates a synthetic binary observation dataset (\\G\\) based on a
given directed acyclic graph topology with Boolean rules. It starts by
running independent Bernoulli trials on root nodes. The remaining
non-root nodes are calculated based on their assigned Boolean logic
functions and parent states. A pre-generated binary noise matrix is
applied via a bitwise XOR operation to occasionally flip the Boolean
outputs, injecting natural biological noise expected by the model.

## Usage

``` r
GenerateSample(
  trans_matrix,
  SampleSize = 50,
  num.node = nrow(trans_matrix),
  para = rep(0.5, nrow(trans_matrix)),
  error = matrix(0, nrow = nrow(trans_matrix), ncol = SampleSize),
  timeseries = TRUE
)
```

## Arguments

- trans_matrix:

  A square matrix combining the network topology \\T\\ and integer-coded
  Boolean logic functions \\F\\ assigned to each directed edge.

- SampleSize:

  An integer representing the total number of samples or time points to
  simulate. Defaults to 50 if not specified.

- num.node:

  An integer representing the total number of network nodes. Defaults to
  `nrow(trans_matrix)` if not specified.

- para:

  A numeric vector of baseline success probabilities (\\\theta_i\\) used
  to generate the expression states of root nodes via independent
  Bernoulli trials. Defaults to a `rep(0.5, nrow(trans_matrix))` if not
  specified.

- error:

  A pre-generated binary noise matrix applied to occasionally flip
  Boolean outputs, injecting natural noise. Defaults to a zero matrix of
  size `nrow(trans_matrix) x SampleSize` if not specified (no noise).

- timeseries:

  Logical. If TRUE, simulates time-lagged data where child state at time
  \\t\\ is determined by parent stated at \\t-1\\. If FALSE, simulates
  independent samples where child and parent states are resolved
  simultaneously. Defaults to TRUE.

## Value

A simulated binary gene expression matrix \\G\\, where rows represent
individual genes/nodes and columns represent samples.

## Examples

``` r
# \donttest{
# 1. Generate a 5-node network
set.seed(123)
num_nodes <- 5
sample_size <- 10
true_network <- GenerateNetwork(num.node = num_nodes)

# 2. Set baseline probabilities and simulate zero-noise error matrix
root_probs <- rep(0.5, num_nodes)
error_matrix <- matrix(0, nrow = num_nodes, ncol = sample_size)

# 3. Generate cross-sectional (independent mode) synthetic data
dummy_data <- GenerateSample(
  trans_matrix = true_network,
  SampleSize = sample_size,
  para = root_probs,
  error = error_matrix,
  timeseries = FALSE
)
print(dummy_data)
#>      [,1] [,2] [,3] [,4] [,5] [,6] [,7] [,8] [,9] [,10]
#> [1,]    1    0    0    0    1    0    0    1    0     1
#> [2,]    1    0    0    0    1    0    0    1    0     1
#> [3,]    1    0    0    0    1    0    0    1    0     1
#> [4,]    1    0    0    0    1    0    0    1    0     1
#> [5,]    0    1    1    1    0    1    1    0    1     0
# }
```
