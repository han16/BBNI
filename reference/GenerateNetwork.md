# Generate Initial Network Topology

Randomly generates a valid directed acyclic graph (DAG) topology \\T\\
and assigns a corresponding Boolean transition function \\F\\ to each
node. The algorithm samples parent set configurations, keeping the
constraint that the maximum in-degree for any node is 2, and further
ensures the resulting structure does not contain directed cyclic loops.

## Usage

``` r
GenerateNetwork(num.node)
```

## Arguments

- num.node:

  An integer representing the total number of genes/nodes in the
  network.

## Value

A square transition function matrix combining the initial DAG topology
with randomly assigned Boolean logic functions. Elements with a value of
0 indicate no directed edge, while positive integers indicate the
presence of an edge and specify the defining Boolean function type
(1-14).

## Examples

``` r
# Generate a true network topology and Boolean rules for 5 nodes
set.seed(123)
true_network <- GenerateNetwork(num.node = 5)
print(true_network)
#>      [,1] [,2] [,3] [,4] [,5]
#> [1,]    0    0    0    0    0
#> [2,]    0    0   11    0    0
#> [3,]   11    0    0    0    0
#> [4,]   11    0    0    0    0
#> [5,]    2    0    0    2    0
```
