#' Plot the Inferred Bayesian Boolean Network
#'
#' @description
#' Visualizes the causal network structure inferred by the BBNI MCMC sampler.
#' This function takes the marginal posterior probability of each directed edge
#' and plots the network using the \code{igraph} package. Edges with posterior
#' probabilities below the specified threshold are omitted from the plot.
#'
#' @param results The list returned by \code{run_bbni()}, containing \code{networks} and \code{log_posterior}.
#' @param threshold Numeric. The minimum posterior probability required to draw an edge. Defaults to 0.5.
#' @param node_names Character vector. Optional names for the nodes. Defaults to "N1", "N2", etc.
#' @param true_network Optional square matrix representing the true network topology. If provided, edges will be color-coded to indicate true positives (along with displaying wrong function inferences), false positives, and false negatives. Purely for simulation purposes.
#' @param ... Additional graphical parameters passed to \code{igraph::plot.igraph()}.
#'
#' @return An invisible \code{igraph} object.
#'
#' @examples
#' # 1. Generate synthetic network and time-series data
#' set.seed(123)
#' true_network <- GenerateNetwork(num.node = 5)
#' dummy_data <- GenerateSample(true_network, SampleSize = 15)
#'
#' # 2. Run BBNI sampler
#' prior_para <- matrix(3, nrow = 6, ncol = 2)
#' prior_para[6, 1] <- 2
#' prior_para[6, 2] <- 100
#'
#' results <- run_bbni(dummy_data, prior_para = prior_para, num_update = 100)
#'
#' # 3. Plot inferred network
#' plot_bbni(results, true_network = true_network, threshold = 0.5)
#'
#' @importFrom igraph graph_from_adjacency_matrix plot.igraph
#' @export
plot_bbni <- function(results, threshold = 0.5, node_names = NULL, true_network = NULL, ...) {
  # Make sure igraph is available
  if (!requireNamespace("igraph", quietly = TRUE)) {
    stop("Package \"igraph\" is required for plotting. Please install it.", call. = FALSE)
  }
  # Extract list of networks and calculate posterior probabilities (res)
  if (!is.null(results$post_edge_prob)) {
    res <- results$post_edge_prob
  } else {
    # fallback by manually calculate
    networks <- results$networks
    res <- Reduce(`+`, lapply(networks, function(m) (m > 0) * 1)) / length(networks)
  }
  # Filter probability matrix based on threshold
  adj_matrix <- ifelse(res >= threshold, 1, 0)
  # Assign node names if provided, otherwise default to N1, N2...
  num_nodes <- nrow(res)
    if (is.null(node_names)) {
    colnames(adj_matrix) <- rownames(adj_matrix) <- paste0("N", 1:num_nodes)
    vsize <- 12   # short N-labels: keep current look
    lcex  <- 0.7
  } else {
    colnames(adj_matrix) <- rownames(adj_matrix) <- node_names
    vsize <- 20   # real gene names (e.g. CDC20), bigger circle + smaller text so labels fit
    lcex  <- 0.6
  }
  # get functions for inferred edges
  inf_func_matrix <- matrix(0, num_nodes, num_nodes)
  # filter out burn-in samples so function inference isn't polluted by early chain states
  total_samples <- length(results$networks)
  burn_in_ratio <- if (!is.null(results$burn_in)) results$burn_in else 0
  burn_in_steps <- floor(burn_in_ratio * (total_samples - 1))
  post_samples <- results$networks[(burn_in_steps + 2):total_samples]
  for (i in 1:num_nodes) {
    for (j in 1:num_nodes) {
      if (adj_matrix[i, j] == 1) {
        vals <- sapply(post_samples, `[` , i, j)
        vals <- vals[vals > 0]
        if (length(vals) > 0) {
          inf_func_matrix[i, j] <- as.integer(names(which.max(table(vals))))
        }
      }
    }
  }
  # if ground truth provided, color-code to show model performance
  if (!is.null(true_network)) {
    true_adj <- (true_network > 0) * 1
    combined_adj <- (adj_matrix | true_adj) * 1
    g <- igraph::graph_from_adjacency_matrix(t(combined_adj), mode = "directed")
    edges <- igraph::as_edgelist(g, names = FALSE)
    idx <- cbind(edges[, 2], edges[, 1])
    is_true <- true_adj[idx] == 1
    is_inf  <- adj_matrix[idx] == 1
    is_correct_func <- (true_network[idx] == inf_func_matrix[idx])
    e_colors <- rep("darkgray", igraph::ecount(g))
    e_colors[is_true & is_inf & is_correct_func] <- "darkgreen"   # True Positive
    e_colors[is_true & is_inf & !is_correct_func] <- "darkorange" # True Positive, but wrong function
    e_colors[!is_true & is_inf] <- "firebrick"  # False Positive
    e_colors[is_true & !is_inf] <- "gray70"     # False Negative (Missed)
    e_lty <- rep(1, igraph::ecount(g))
    e_lty[is_true & !is_inf] <- 2               # Dashed for missed edges
    igraph::E(g)$color <- e_colors
    igraph::E(g)$lty <- e_lty
    main_title <- paste("Inferred vs True Network (Threshold >=", threshold, ")")
  } else {
    # default plotting w/o ground truth
    g <- igraph::graph_from_adjacency_matrix(t(adj_matrix), mode = "directed")
    igraph::E(g)$color <- "darkgray"
    igraph::E(g)$lty <- 1
    main_title <- paste("Inferred Network (Threshold >=", threshold, ")")
  }
  # plot graph
    # Capture extra arguments (like layout)
  dots <- list(...)
  # plot graph
  do.call(igraph::plot.igraph, c(list(g,
    vertex.size = vsize,
    vertex.color = "lightblue",
    vertex.label.color = "black",
    vertex.label.cex = lcex,
    edge.arrow.size = 0.4,
    edge.curved = 0.1,
    main = main_title), dots))
  # add legend if ground truth was provided
  if (!is.null(true_network)) {
    legend("bottomleft",
           legend = c("Correct (TP)", "Correct (TP), Wrong Function", "Spurious (FP)", "Missed (FN)"),
           col = c("darkgreen", "darkorange", "firebrick", "gray70"),
           lty = c(1, 1, 1, 2), lwd = 2, bty = "n", cex = 0.8)
  }
  # return graph object invisibly in case user wants to manipulate it further
  invisible(g)
}

#' Plot a Single Boolean Network
#'
#' Visualizes a single directed network topology. For instance, a synthetic "true" network
#' generated by \code{GenerateNetwork()}.
#'
#' @param trans_matrix A square matrix representing the network topology and Boolean rules (e.g., the output of \code{GenerateNetwork()}). Elements greater than 0 indicate a directed edge.
#' @param node_names Character vector. Optional names for the nodes. Defaults to "N1", "N2", etc.
#' @param ... Additional graphical parameters passed to \code{igraph::plot.igraph()}.
#'
#' @return An invisible \code{igraph} object.
#'
#' @examples
#' # 1. Generate synthetic 5-node network topology
#' set.seed(123)
#' true_network <- GenerateNetwork(num.node = 5)
#'
#' # 2. Plot the generated network
#' plot_network(true_network)
#'
#' @importFrom igraph graph_from_adjacency_matrix plot.igraph
#' @export
plot_network <- function(trans_matrix, node_names = NULL, ...) {
  # Make sure igraph is available
  if (!requireNamespace("igraph", quietly = TRUE)) {
    stop("Package \"igraph\" is required for plotting. Please install it.", call. = FALSE)
  }
  adj_matrix <- (trans_matrix > 0) * 1
  num_nodes <- nrow(adj_matrix)
  if (is.null(node_names)) {
    colnames(adj_matrix) <- rownames(adj_matrix) <- paste0("N", 1:num_nodes)
  } else {
    colnames(adj_matrix) <- rownames(adj_matrix) <- node_names
  }
  g <- igraph::graph_from_adjacency_matrix(t(adj_matrix), mode = "directed")
  igraph::plot.igraph(g, vertex.size = 20, vertex.color = "lightblue",
                      vertex.label.color = "black", vertex.label.cex = 0.8,
                      edge.arrow.size = 0.5, edge.color = "darkgray",
                      main = "Network Topology", ...)
  invisible(g)
}

#' Plot MCMC Trace for BBNI
#'
#' @description
#' Generates a trace plot of the log-posterior values over iterations of the MCMC
#' to visually assess convergence and stability of the executed Markov chain.
#' Burn-in line is graphed to show when the data started to be utilized for
#' edge-probability calculations.
#'
#' @param results The list returned by \code{run_bbni()}, containing \code{networks} and \code{log_posterior}.
#' @param every An integer specifying the thinning interval (sampling frequency) for plotting. Default is 1, which plots all log-posterior values. Values greater than 1 plot every \code{every}-th iteration.
#'
#' @return A base R trace plot.
#'
#' @examples
#' # 1. Generate synthetic network and time-series data
#' set.seed(123)
#' true_network <- GenerateNetwork(num.node = 5)
#' dummy_data <- GenerateSample(true_network, SampleSize = 15)
#'
#' # 2. Run BBNI sampler
#' prior_para <- matrix(3, nrow = 6, ncol = 2)
#' prior_para[6, 1] <- 2
#' prior_para[6, 2] <- 100
#'
#' results <- run_bbni(dummy_data, prior_para = prior_para, num_update = 100)
#'
#' # 3. Visualize MCMC results
#' plot_trace(results)
#'
#' @importFrom graphics plot abline legend
#' @export
plot_trace <- function(results, every = 1) {
  logpost <- results$log_posterior
  if (every > 1) logpost <- logpost[seq(1, length(logpost), by = every)]
  plot(logpost, type = "l", col = "darkblue",
      xlab = if (every == 1) "Node-level update" else "Outer iteration", 
      ylab = "Log-Posterior",
      main = "MCMC Trace Plot",
      lwd = 1.5)
  # add vertical line for burn-in if present
  if (!is.null(results$burn_in)) {
    abline(v = results$burn_in * length(logpost), col = "red", lty = 2)
    legend("bottomright", legend = paste("Burn-in =", results$burn_in), col = "red", lty = 2, bty = "n")
  }
}
