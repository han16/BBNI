# Plot the Inferred Bayesian Boolean Network

Plot the Inferred Bayesian Boolean Network

## Usage

``` r
plot_bbni(results, threshold = 0.5, node_names = NULL, ...)
```

## Arguments

- results:

  The list returned by
  [`run_bbni()`](https://han16.github.io/BBNI/reference/run_bbni.md),
  containing `networks` and `log_posterior`.

- threshold:

  Numeric. The minimum posterior probability required to draw an edge.
  Defaults to 0.5.

- node_names:

  Character vector. Optional names for the nodes. Defaults to "N1",
  "N2", etc.

- ...:

  Additional graphical parameters passed to
  [`igraph::plot.igraph()`](https://r.igraph.org/reference/plot.igraph.html).

## Value

An invisible `igraph` object.
