# Simulate Time-Series Observation Dataset

Simulates a synthetic time-series observation dataset (\\G\\). It starts
by running independent Bernoulli trials on root nodes. The remaining
non-root nodes are calculated from time \\t-1\\ to time \\t\\ using
their assigned Boolean logic functions. A pre-generated binary noise
matrix is applied via a bitwise XOR operation to occasionally flip the
Boolean outputs, injecting natural biological noise expected by the
model.

## Usage

``` r
GenerateSample(trans_matrix, num.node, SampleSize, para, error)
```

## Arguments

- trans_matrix:

  A square matrix combining the network topology \\T\\ and integer-coded
  Boolean logic functions \\F\\ assigned to each directed edge.

- num.node:

  An integer representing the total number of network nodes.

- SampleSize:

  An integer representing the total number of time points to simulate.

- para:

  A numeric vector of baseline success probabilities (\\\theta_i\\) used
  to generate the expression states of root nodes via independent
  Bernoulli trials.

- error:

  A pre-generated binary noise matrix applied to occasionally flip
  Boolean outputs, injecting natural noise.

## Value

A simulated binary gene expression matrix \\G\\, where rows represent
individual genes/nodes and columns represent sequential points in time.

## Examples

``` r
# 1. Generate a 5-node network
set.seed(123)
num_nodes <- 5
sample_size <- 10
true_network <- GenerateNetwork(num.node = num_nodes)

# 2. Set baseline probabilities and simulate zero-noise error matrix
root_probs <- rep(0.5, num_nodes)
error_matrix <- matrix(0, nrow = num_nodes, ncol = sample_size)

# 3. Generate the synthetic time-series data
dummy_data <- GenerateSample(
  trans_matrix = true_network,
  num.node = num_nodes,
  SampleSize = sample_size,
  para = root_probs,
  error = error_matrix
)
print(dummy_data)
#>      [,1] [,2] [,3] [,4] [,5] [,6] [,7] [,8] [,9] [,10]
#> [1,]    1    0    0    0    1    0    0    1    0     1
#> [2,]    0    0    1    0    0    0    1    0    0     1
#> [3,]    0    1    0    0    0    1    0    0    1     0
#> [4,]    0    1    0    0    0    1    0    0    1     0
#> [5,]    0    1    1    1    1    1    1    1    1     1
```
