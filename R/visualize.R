#' Plot the Inferred Bayesian Boolean Network
#'
#' @param results The list returned by \code{run_bbni()}, containing \code{networks} and \code{log_posterior}.
#' @param threshold Numeric. The minimum posterior probability required to draw an edge. Defaults to 0.5.
#' @param node_names Character vector. Optional names for the nodes. Defaults to "N1", "N2", etc.
#' @param ... Additional graphical parameters passed to \code{igraph::plot.igraph()}.
#'
#' @return An invisible \code{igraph} object.
#' @export
#'
#' @importFrom igraph graph_from_adjacency_matrix plot.igraph
plot_bbni <- function(results, threshold = 0.5, node_names = NULL, ...) {
  # Make sure igraph is available
  if (!requireNamespace("igraph", quietly = TRUE)) {
    stop("Package \"igraph\" is required for plotting. Please install it.", call. = FALSE)
  }
  # Extract list of networks and calculate posterior probabilities (res)
  networks <- results$networks
  res <- Reduce(`+`, lapply(networks, function(m) (m > 0) * 1)) / length(networks)
  # Filter probability matrix based on threshold
  adj_matrix <- ifelse(res >= threshold, 1, 0)
  # Assign node names if provided, otherwise default to N1, N2...
  num_nodes <- nrow(res)
  if (is.null(node_names)) {
    colnames(adj_matrix) <- rownames(adj_matrix) <- paste0("N", 1:num_nodes)
  } else {
    colnames(adj_matrix) <- rownames(adj_matrix) <- node_names
  }
  # Create directed graph
  g <- igraph::graph_from_adjacency_matrix(adj_matrix, mode = "directed")
  # Plot with defaults
  igraph::plot.igraph(g,
                      vertex.size = 20,
                      vertex.color = "lightblue",
                      vertex.label.color = "black",
                      vertex.label.cex = 0.8,
                      edge.arrow.size = 0.5,
                      edge.color = "darkgray",
                      main = paste("Inferred Network (Threshold >", threshold, ")"),
                      ...)
  # Return graph object invisibly in case user wants to manipulate it further
  invisible(g)
}

#' Plot MCMC Trace for BBNI
#'
#' @param results The list returned by \code{run_bbni()}, containing \code{networks} and \code{log_posterior}.
#'
#' @return A base R trace plot.
#' @export
plot_trace <- function(results) {
  logpost <- results$log_posterior
  plot(logpost, type = "l", col = "darkblue",
       xlab = "Iteration", ylab = "Log-Posterior",
       main = "MCMC Trace Plot",
       lwd = 1.5)
}
