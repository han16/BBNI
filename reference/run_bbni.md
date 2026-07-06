# Execute Metropolis-within-Gibbs MCMC Sampler for Boolean Networks

Executes a Metropolis-within-Gibbs Markov chain Monte Carlo (MCMC)
algorithm to sample from the joint posterior distribution of Directed
Acyclic Graph (DAG) topologies (\\T\\) and Boolean logic transition
functions (\\F\\). The algorithm iterates through individual network
nodes and proposes parent set mutations (edge additions, removals, or
swaps) paired with transition function reassignments to one of 14
candidate Boolean rules. Proposed states transitions are strictly
verified to follow the DAG constraint and evaluated with a
Metropolis-Hastings acceptance threshold using log-posterior values.

## Usage

``` r
run_bbni(
  GeneData,
  num.node,
  SampleSize,
  prior_para,
  num_update,
  penalty,
  prop.ratio,
  verbose = FALSE
)
```

## Arguments

- GeneData:

  A binary empirical observation matrix of the binary expression data
  (\\G\\).

- num.node:

  An integer representing the total number of network nodes.

- SampleSize:

  An integer representing the total number of time points in the
  dataset.

- prior_para:

  A matrix of Beta prior hyperparameters \\\alpha\\ and \\\beta\\ for
  root node probabilities and the global noise parameter \\e\\.

- num_update:

  An integer representing the total number of MCMC iterations to
  perform.

- penalty:

  A numeric value representing the structural prior probability per edge
  used to penalize network complexity \\P(T)\\.

- prop.ratio:

  A numeric probability threshold used to decide whether to sample a
  move from the empirical proposal distribution or a uniform random
  distribution.

- verbose:

  Logical. If TRUE, prints verbose MCMC iteration progress to the
  console. Default is FALSE.

## Value

A list containing the full trajectory of the MCMC chain. Specifically,
`networks` (a list of sampled transition function matrices) and
`log_posterior` (a numeric vector of log-posterior scores for each
iteration). These represent samples drawn from the marginal posterior
distribution \\P(T,F\|G)\\ used for Bayesian model averaging.

## Examples

``` r
# \donttest{
  # 1. Define network parameters
  set.seed(235)
  num_nodes <- 10
  sample_size <- 50

  # 2. Generate true network and simulate data
  true_network <- GenerateNetwork(num.node = num_nodes)

  # Set up Beta priors for root-node probabilities and the noise rate
  prior_para <- matrix(3, nrow = num_nodes + 1, ncol = 2)
  prior_para[num_nodes + 1, 1] <- 2
  prior_para[num_nodes + 1, 2] <- 100

  # Simulate parameters
  para <- numeric(num_nodes + 1)
  for (i in 1:(num_nodes + 1)) {
    para[i] <- stats::rbeta(1, prior_para[i, 1], prior_para[i, 2])
  }
  para[num_nodes + 1] <- 0.1 # Fixed noise rate for simulation

  error_matrix <- matrix(stats::rbinom(num_nodes * sample_size, 1, para[num_nodes + 1]),
    nrow = num_nodes, ncol = sample_size
  )

  dummy_data <- GenerateSample(
    trans_matrix = true_network,
    num.node = num_nodes,
    SampleSize = sample_size,
    para = para,
    error = error_matrix
  )

  # 3. Run the MCMC sampler (silently)
  mcmc_results <- run_bbni(
    GeneData = dummy_data,
    num.node = num_nodes,
    SampleSize = sample_size,
    prior_para = prior_para,
    num_update = 100, # Scaled down for example speed
    penalty = 0.1,
    prop.ratio = 0.1
  )

  # 4. Inspect results
  tail(mcmc_results$log_posterior)
#> [1] -279.5854 -279.5854 -281.8880 -281.8880 -281.8880 -281.8880
# }
```
