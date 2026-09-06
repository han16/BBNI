# Posterior Summary of a BBNI MCMC Object

Prints a statistical summary of the MCMC run and returns a structured
list of posterior quantities for further usage in downstream code.

## Usage

``` r
# S3 method for class 'bbni'
summary(object, threshold = 0.5, n_top = 5, ...)
```

## Arguments

- object:

  A `bbni` object returned by
  [`run_bbni()`](https://han16.github.io/BBNI/reference/run_bbni.md).

- threshold:

  Numeric. Posterior probability threshold for counting strong edges.

- n_top:

  Integer. Number of highest-posterior edges to display and include in
  `top_edges`.

- ...:

  Additional arguments (ignored).

## Value

A list with components `num_nodes`, `num_update`, `burn_in`,
`final_log_posterior`, `n_strong_edges`, and `top_edges` (data frame
with columns `child`, `parent`, and `posterior`).

## Examples

``` r
set.seed(123)
net <- GenerateNetwork(5)
data <- GenerateSample(net, 100)
res <- run_bbni(data, num_update = 50)
summary(res, threshold = 0.6, n_top = 3)
#> BBNI MCMC summary
#>   Nodes:                 5
#>   Outer iterations:      50
#>   Burn-in ratio:         0.70
#>   Final log-posterior:   -87.922
#>   Edges with P > 0.60:    5
#>   Highest-posterior edges (child <- parent):
#>     N3 <- N1 (P = 1.00)
#>     N4 <- N1 (P = 1.00)
#>     N5 <- N1 (P = 1.00)
```
