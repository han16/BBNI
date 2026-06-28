# BBNI

BBNI is a Bayesian approach to Boolean gene regulatory network inference
from noisy gene expression data. The method is discussed in more detail
in [Han et al. (2014)](https://doi.org/10.1371/journal.pone.0115806).
Unlike methods that return a single best-fit network topology, such as
REVEAL and BFE, found in
[BoolNet](https://cran.r-project.org/package=BoolNet), BBNI uses Markov
chain Monte Carlo (MCMC) to sample from a joint posterior distribution
of network topologies and Boolean transition functions. BBNI
deliberately takes biological noise into account and allows for
summarization (such as Bayesian model averaging, or BMA) that stabilizes
results through posterior edge probabilities rather than a single point
estimate.

## Installation

BBNI is not currently available on CRAN. The development version can be
installed from GitHub:

``` r

# install.packages("devtools")
devtools::install_github("anson-li8/BBNI")
```

## Example

The following example provides a minimal check that the package loads
successfully and can generate simulated data:

``` r

library(BBNI)

set.seed(1)
true_network <- GenerateNetwork(num.node = 5)
dummy_data <- GenerateSample(
  trans_matrix = true_network,
  num.node = 5,
  SampleSize = 20,
  para = rep(0.5, 5),
  error = matrix(0, nrow = 5, ncol = 20)
)
dummy_data
#>      [,1] [,2] [,3] [,4] [,5] [,6] [,7] [,8] [,9] [,10] [,11] [,12] [,13] [,14]
#> [1,]    0    0    0    0    0    1    1    1    0     0     1     1     1     0
#> [2,]    0    0    0    0    0    0    0    1    0     0     0     0     1     0
#> [3,]    0    1    1    1    1    1    0    0    0     1     1     0     0     0
#> [4,]    0    0    0    0    0    1    1    0    0     0     1     1     0     0
#> [5,]    1    1    1    1    0    0    1    1    1     0     0     1     1     0
#>      [,15] [,16] [,17] [,18] [,19] [,20]
#> [1,]     1     1     1     1     0     0
#> [2,]     0     0     0     1     0     0
#> [3,]     1     0     0     0     0     1
#> [4,]     1     1     1     0     0     0
#> [5,]     0     0     1     1     1     1
```

For a complete demonstration running the MCMC sampler and evaluating
final convergence and overall network recovery, see the [Introduction to
BBNI](https://anson-li8.github.io/BBNI/articles/Introduction_to_BBNI.html)
vignette.

## Citation

``` R
To cite package 'BBNI' in publications use:

  Han S, Wong RKW, Lee TCM, Shen L, Li S-YR, Fan X (2014). A Full
  Bayesian Approach for Boolean Genetic Network Inference. PLoS ONE
  9(12): e115806. doi:10.1371/journal.pone.0115806
```

## License

BSD-3-Clause
