#' Calculate Transitive Closure of Adjacency Matrix
#'
#' Calculates the transitive closure of the network's directed adjacency matrix
#' to locate all ancestor-descendant relationships. This is used by the MCMC
#' algorithm to ensure proposed network topologies represent directed acyclic
#' graphs (DAGs) and to prevent illegal cyclic moves when adding or swapping parents.
#' @param incid_matrix A square binary adjacency matrix representing the current network topology.
#'
#' @returns A binary matrix of the same dimensions as `incid_matrix`. An entry of 1 at (i, k) indicates that node i is an ancestor of node k through one or more directed edges.
#' @noRd
update_ancestor_matrix <- function(incid_matrix) {
  # BLAS matrix-power squaring can speed this up even more
  n <- nrow(incid_matrix)
  R <- incid_matrix
  steps <- ceiling(log2(n))
  for (s in seq_len(steps)) {
    R <- (R %*% R + R) > 0
  }
  storage.mode(R) <- "double"
  return(R)
}

#' Check Ancestor Matrix for Cyclic Loops
#'
#' Ensures that the proposed network topology fulfills the directed acyclic graph
#' (DAG) constraint by checking the diagonal of the ancestor matrix for illegal
#' cyclic loops. This enforces the requirement that adding or swapping parents
#' cannot lead to a cyclic graph during MCMC topology updates.
#'
#' @param ances_matrix A binary transitive closure matrix of the network topology, where an entry of 1 at (i, j) means node i is an ancestor of node j.
#'
#' @return An integer count of the number of nodes that are ancestors of themselves. A value greater than 0 indicates the proposed topology has cyclic loops and must be rejected.
#' @noRd
check_ances_matrix <- function(ances_matrix) {
  # vectorized
  sum(diag(ances_matrix) == 1)
}
