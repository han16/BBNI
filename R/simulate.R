#' Generate Initial Network Topology
#'
#' Randomly generates a valid directed acyclic graph (DAG) topology
#' \eqn{T}{T} and assigns a corresponding Boolean transition function \eqn{F}{F} to each node.
#' The algorithm samples parent set configurations, keeping the constraint
#' that the maximum in-degree for any node is 2, and further ensures
#' the resulting structure does not contain directed cyclic loops.
#'
#' @param num.node An integer representing the total number of genes/nodes in the network.
#'
#' @return A square transition function matrix combining the initial DAG topology with randomly assigned Boolean logic functions. Elements with a value of 0 indicate no directed edge, while positive integers indicate the presence of an edge and specify the defining Boolean function type (1-14).
#'
#' @examples
#' # Generate a true network topology and Boolean rules for 5 nodes
#' set.seed(123)
#' true_network <- GenerateNetwork(num.node = 5)
#' print(true_network)
#'
#' @importFrom stats rbinom runif
#' @export
GenerateNetwork <- function(num.node) {
  loop <- 1
  while (loop != 0) {
    tent_incid_matrix <- matrix(0, nrow = num.node, ncol = num.node)
    tent_trans_matrix <- matrix(0, nrow = num.node, ncol = num.node)
    for (i in 1:num.node) {
      u1 <- runif(1)
      uu <- 4
      if (u1 > 1 / 10 && u1 <= uu / 10) {
        position <- sample(seq(num.node)[-i], 1)
        tent_incid_matrix[i, position] <- 1
        u2 <- runif(1)
        if (u2 > 0.5) {
          tent_trans_matrix[i, position] <- 11
        }
        if (u2 < 0.5) {
          tent_trans_matrix[i, position] <- 12
        }
      }
      if (u1 > uu / 10) {
        position <- sample(seq(num.node)[-i], 2)
        tent_incid_matrix[i, position[1]] <- 1
        tent_incid_matrix[i, position[2]] <- 1
        func <- sample(10, 1)
        tent_trans_matrix[i, position[1]] <- func
        tent_trans_matrix[i, position[2]] <- func
      }
    }
    tent_ances_matrix <- update_ancestor_matrix(tent_incid_matrix)
    loop <- check_ances_matrix(tent_ances_matrix)
  }
  return(tent_trans_matrix)
}
#' Simulate Time-Series Observation Dataset
#'
#' Simulates a synthetic time-series observation dataset (\eqn{G}{G}).
#' It starts by running independent Bernoulli trials on root nodes.
#' The remaining non-root nodes are calculated from time \eqn{t-1}{t-1} to time \eqn{t}{t} using their assigned Boolean logic functions.
#' A pre-generated binary noise matrix is applied via a bitwise XOR operation
#' to occasionally flip the Boolean outputs, injecting natural biological
#' noise expected by the model.
#'
#' @param trans_matrix A square matrix combining the network topology \eqn{T}{T} and integer-coded Boolean logic functions \eqn{F}{F} assigned to each directed edge.
#' @param SampleSize An integer representing the total number of time points to simulate. Defaults to 50 if not specified.
#' @param num.node An integer representing the total number of network nodes. Defaults to \code{nrow(trans_matrix)} if not specified.
#' @param para A numeric vector of baseline success probabilities (\eqn{\theta_i}{\theta_i}) used to generate the expression states of root nodes via independent Bernoulli trials. Defaults to a \code{rep(0.5, nrow(trans_matrix))} if not specified.
#' @param error A pre-generated binary noise matrix applied to occasionally flip Boolean outputs, injecting natural noise. Defaults to a zero matrix of size \code{nrow(trans_matrix) x SampleSize} if not specified (no noise).
#'
#' @return A simulated binary gene expression matrix \eqn{G}{G}, where rows represent individual genes/nodes and columns represent sequential points in time.
#'
#' @examples
#' # 1. Generate a 5-node network
#' set.seed(123)
#' num_nodes <- 5
#' sample_size <- 10
#' true_network <- GenerateNetwork(num.node = num_nodes)
#'
#' # 2. Set baseline probabilities and simulate zero-noise error matrix
#' root_probs <- rep(0.5, num_nodes)
#' error_matrix <- matrix(0, nrow = num_nodes, ncol = sample_size)
#'
#' # 3. Generate the synthetic time-series data
#' dummy_data <- GenerateSample(
#'   trans_matrix = true_network,
#'   num.node = num_nodes,
#'   SampleSize = sample_size,
#'   para = root_probs,
#'   error = error_matrix
#' )
#' print(dummy_data)
#'
#' @importFrom bitops bitXor bitAnd bitOr
#' @export
GenerateSample <- function(trans_matrix, SampleSize = 50, num.node = nrow(trans_matrix), para = rep(0.5, nrow(trans_matrix)), error = matrix(0, nrow = nrow(trans_matrix), ncol = SampleSize)) {
  node_ances <- matrix(nrow = num.node, ncol = 2)
  GeneData <- matrix(0, nrow = num.node, ncol = SampleSize)
  incid_matrix <- trans_matrix
  for (i in seq_len(nrow(trans_matrix))) {
    for (j in seq_len(ncol(trans_matrix))) {
      if (trans_matrix[i, j] > 0) {
        incid_matrix[i, j] <- 1
      }
    }
  }
  ances_matrix <- update_ancestor_matrix(incid_matrix)
  for (i in 1:num.node) {
    node_ances[i, 1] <- i
    node_ances[i, 2] <- sum(ances_matrix[i, ])
  }
  node_ances <- node_ances[order(node_ances[, 2]), ]
  for (i in seq_len(nrow(node_ances))) {
    if (node_ances[i, 2] == 0) {
      GeneData[node_ances[i, 1], ] <- rbinom(SampleSize, 1, prob = para[node_ances[i, 1]])
    }
    if (node_ances[i, 2] != 0) {
      parent <- numeric()
      ii <- 1
      for (j in seq_len(ncol(incid_matrix))) {
        if (incid_matrix[node_ances[i, 1], j] != 0) {
          parent[ii] <- j
          ii <- ii + 1
        }
      }
      func <- trans_matrix[node_ances[i, 1], parent[1]]
      if (func == 1) {
        GeneData[node_ances[i, 1], 2:SampleSize] <- bitXor(bitAnd(GeneData[parent[1], 1:(SampleSize - 1)], GeneData[parent[2], 1:(SampleSize - 1)]), error[node_ances[i, 1], 1:(SampleSize - 1)])
      }
      if (func == 2) {
        GeneData[node_ances[i, 1], 2:SampleSize] <- bitXor(1 - bitAnd(GeneData[parent[1], 1:(SampleSize - 1)], GeneData[parent[2], 1:(SampleSize - 1)]), error[node_ances[i, 1], 1:(SampleSize - 1)])
      }
      if (func == 3) {
        GeneData[node_ances[i, 1], 2:SampleSize] <- bitXor(bitOr(GeneData[parent[1], 1:(SampleSize - 1)], GeneData[parent[2], 1:(SampleSize - 1)]), error[node_ances[i, 1], 1:(SampleSize - 1)])
      }
      if (func == 4) {
        GeneData[node_ances[i, 1], 2:SampleSize] <- bitXor(1 - bitOr(GeneData[parent[1], 1:(SampleSize - 1)], GeneData[parent[2], 1:(SampleSize - 1)]), error[node_ances[i, 1], 1:(SampleSize - 1)])
      }
      if (func == 5) {
        GeneData[node_ances[i, 1], 2:SampleSize] <- bitXor(bitOr(1 - GeneData[parent[1], 1:(SampleSize - 1)], GeneData[parent[2], 1:(SampleSize - 1)]), error[node_ances[i, 1], 1:(SampleSize - 1)])
      }
      if (func == 6) {
        GeneData[node_ances[i, 1], 2:SampleSize] <- bitXor(bitOr(GeneData[parent[1], 1:(SampleSize - 1)], 1 - GeneData[parent[2], 1:(SampleSize - 1)]), error[node_ances[i, 1], 1:(SampleSize - 1)])
      }
      if (func == 7) {
        GeneData[node_ances[i, 1], 2:SampleSize] <- bitXor(bitAnd(1 - GeneData[parent[1], 1:(SampleSize - 1)], GeneData[parent[2], 1:(SampleSize - 1)]), error[node_ances[i, 1], 1:(SampleSize - 1)])
      }
      if (func == 8) {
        GeneData[node_ances[i, 1], 2:SampleSize] <- bitXor(bitAnd(GeneData[parent[1], 1:(SampleSize - 1)], 1 - GeneData[parent[2], 1:(SampleSize - 1)]), error[node_ances[i, 1], 1:(SampleSize - 1)])
      }
      if (func == 9) {
        GeneData[node_ances[i, 1], 2:SampleSize] <- bitXor(bitXor(GeneData[parent[1], 1:(SampleSize - 1)], GeneData[parent[2], 1:(SampleSize - 1)]), error[node_ances[i, 1], 1:(SampleSize - 1)])
      }
      if (func == 10) {
        GeneData[node_ances[i, 1], 2:SampleSize] <- bitXor(1 - bitXor(GeneData[parent[1], 1:(SampleSize - 1)], GeneData[parent[2], 1:(SampleSize - 1)]), error[node_ances[i, 1], 1:(SampleSize - 1)])
      }
      if (func == 11) {
        GeneData[node_ances[i, 1], 2:SampleSize] <- bitXor(GeneData[parent[1], 1:(SampleSize - 1)], error[node_ances[i, 1], 1:(SampleSize - 1)])
      }
      if (func == 12) {
        GeneData[node_ances[i, 1], 2:SampleSize] <- bitXor(1 - GeneData[parent[1], 1:(SampleSize - 1)], error[node_ances[i, 1], 1:(SampleSize - 1)])
      }
    }
  }
  return(GeneData)
}
