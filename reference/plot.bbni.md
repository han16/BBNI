# Plot a BBNI MCMC Object

Produces a trace plot of the log-posterior over MCMC iterations. This
method delegates to
[`plot_trace()`](https://han16.github.io/BBNI/reference/plot_trace.md).
For network visualization, use
[`plot_bbni()`](https://han16.github.io/BBNI/reference/plot_bbni.md).

## Usage

``` r
# S3 method for class 'bbni'
plot(x, ...)
```

## Arguments

- x:

  A `bbni` object returned by
  [`run_bbni()`](https://han16.github.io/BBNI/reference/run_bbni.md).

- ...:

  Extra arguments passed to
  [`plot_trace()`](https://han16.github.io/BBNI/reference/plot_trace.md).

## Value

The plot object (invisibly).

## Examples

``` r
set.seed(123)
net <- GenerateNetwork(5)
data <- GenerateSample(net, 100)
res <- run_bbni(data, num_update = 50)
plot(res)
```
