# BBNI MCMC Output Class

[`run_bbni()`](https://han16.github.io/BBNI/reference/run_bbni.md)
returns a list of class `bbni` containing the sampled path of the
Metropolis-within-Gibbs chain. Classed objects have
[`print()`](https://rdrr.io/r/base/print.html),
[`summary()`](https://rdrr.io/r/base/summary.html), and
[`plot()`](https://rdrr.io/r/graphics/plot.default.html) methods; all
list components remain directly accessible via `$`.

## Details

Components:

- `networks`: list of MCMC sampled transition-function matrices.

- `log_posterior`: numeric vector of collapsed log-posterior values.

- `post_edge_prob`: matrix of marginal posterior edge probabilities;
  entry `[i, j]` is the probability of the directed edge `j -> i`.

- `burn_in`: the burn-in ratio used for posterior summarization.

Run metadata (`num.node`, `SampleSize`, `num_update`, `timeseries`) is
stored as attributes.

## See also

[`run_bbni()`](https://han16.github.io/BBNI/reference/run_bbni.md),
[`print.bbni()`](https://han16.github.io/BBNI/reference/print.bbni.md),
[`summary.bbni()`](https://han16.github.io/BBNI/reference/summary.bbni.md),
[`plot.bbni()`](https://han16.github.io/BBNI/reference/plot.bbni.md)
