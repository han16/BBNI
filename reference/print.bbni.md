# Print a BBNI MCMC Object

The default display method for `bbni` objects. It provides a one-line
overview instead of a statistical summary; use
[`summary()`](https://rdrr.io/r/base/summary.html) for posterior
statistics and [`plot()`](https://rdrr.io/r/graphics/plot.default.html)
for trace plots.

## Usage

``` r
# S3 method for class 'bbni'
print(x, ...)
```

## Arguments

- x:

  A `bbni` object returned by
  [`run_bbni()`](https://han16.github.io/BBNI/reference/run_bbni.md).

- ...:

  Additional arguments (ignored).

## Value

`x`, invisibly.

## Examples

``` r
set.seed(123)
net <- GenerateNetwork(5)
data <- GenerateSample(net, 100)
res <- run_bbni(data, num_update = 50)
print(res)
#> <bbni> Bayesian Boolean Network Inference MCMC output
#> 5 nodes | 50 outer iterations | time-series data
#> Use summary() for posterior edge statistics and plot() for the trace plot.
```
