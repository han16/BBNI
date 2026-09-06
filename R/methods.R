#' BBNI MCMC Output Class
#'
#' @description
#' `run_bbni()` returns a list of class `bbni` containing the sampled path of
#' the Metropolis-within-Gibbs chain. Classed objects have `print()`, `summary()`,
#' and `plot()` methods; all list components remain directly accessible via `$`.
#'
#' @details
#' Components:
#' - `networks`: list of MCMC sampled transition-function matrices.
#' - `log_posterior`: numeric vector of collapsed log-posterior values.
#' - `post_edge_prob`: matrix of marginal posterior edge probabilities; entry
#'   `[i, j]` is the probability of the directed edge `j -> i`.
#' - `burn_in`: the burn-in ratio used for posterior summarization.
#'
#' Run metadata (`num.node`, `SampleSize`, `num_update`, `timeseries`) is
#' stored as attributes.
#'
#' @name bbni_class
#' @seealso [run_bbni()], [print.bbni()], [summary.bbni()], [plot.bbni()]
NULL

#' Print a BBNI MCMC Object
#'
#' The default display method for `bbni` objects. It provides
#' a one-line overview instead of a statistical summary; use
#' `summary()` for posterior statistics and `plot()` for trace plots.
#'
#' @param x A `bbni` object returned by `run_bbni()`.
#' @param ... Additional arguments (ignored).
#' @return `x`, invisibly.
#' @method print bbni
#' @export
#' @examples
#' set.seed(123)
#' net <- GenerateNetwork(5)
#' data <- GenerateSample(net, 100)
#' res <- run_bbni(data, num_update = 50)
#' print(res)
print.bbni <- function(x, ...) {
  n <- nrow(x$post_edge_prob)
  iters <- attr(x, "num_update")
  if (is.null(iters)) iters <- (length(x$networks) - 1L) / n
  cat("<bbni> Bayesian Boolean Network Inference MCMC output\n")
  cat(sprintf("%d nodes | %d outer iterations | %s data\n",
              n, iters,
              if (isTRUE(attr(x, "timeseries"))) "time-series" else "independent"))
  cat("Use summary() for posterior edge statistics and plot() for the trace plot.\n")
  invisible(x)
}

#' Posterior Summary of a BBNI MCMC Object
#'
#' Prints a statistical summary of the MCMC run and returns a structured list
#' of posterior quantities for further usage in downstream code.
#'
#' @param object A `bbni` object returned by `run_bbni()`.
#' @param threshold Numeric. Posterior probability threshold for counting
#'   strong edges.
#' @param n_top Integer. Number of highest-posterior edges to display and
#'   include in `top_edges`.
#' @param ... Additional arguments (ignored).
#' @return A list with components `num_nodes`, `num_update`, `burn_in`,
#'   `final_log_posterior`, `n_strong_edges`, and `top_edges`
#'   (data frame with columns `child`, `parent`, and `posterior`).
#' @method summary bbni
#' @export
#' @examples
#' set.seed(123)
#' net <- GenerateNetwork(5)
#' data <- GenerateSample(net, 100)
#' res <- run_bbni(data, num_update = 50)
#' summary(res, threshold = 0.6, n_top = 3)
summary.bbni <- function(object, threshold = 0.5, n_top = 5, ...) {
  pep <- object$post_edge_prob
  n <- nrow(pep)
  iters <- attr(object, "num_update")
  if (is.null(iters)) iters <- (length(object$networks) - 1L) / n
  rn <- rownames(pep)
  if (is.null(rn)) rn <- paste0("N", seq_len(n))
  m <- pep
  diag(m) <- NA
  v <- as.vector(m)
  na_pos <- is.na(v)
  v[na_pos] <- -Inf
  ord <- order(v, decreasing = TRUE)
  ord <- ord[!na_pos[ord]]
  top <- ord[seq_len(min(n_top, length(ord)))]
  ii <- (top - 1L) %% n + 1L
  jj <- (top - 1L) %/% n + 1L
  cat("BBNI MCMC summary\n")
  cat(sprintf("  Nodes:                 %d\n", n))
  cat(sprintf("  Outer iterations:      %d\n", iters))
  cat(sprintf("  Burn-in ratio:         %.2f\n", object$burn_in))
  cat(sprintf("  Final log-posterior:   %.3f\n",
              object$log_posterior[length(object$log_posterior)]))
  cat(sprintf("  Edges with P > %.2f:    %d\n", threshold,
              sum(pep > threshold & row(pep) != col(pep))))
  if (length(top) > 0L) {
    cat("  Highest-posterior edges (child <- parent):\n")
    for (k in seq_along(top)) {
      cat(sprintf("    %s <- %s (P = %.2f)\n", rn[ii[k]], rn[jj[k]], m[ii[k], jj[k]]))
    }
  }
  invisible(list(
    num_nodes = n,
    num_update = iters,
    burn_in = object$burn_in,
    final_log_posterior = object$log_posterior[length(object$log_posterior)],
    n_strong_edges = sum(pep > threshold & row(pep) != col(pep)),
    top_edges = data.frame(
      child = rn[ii],
      parent = rn[jj],
      posterior = m[cbind(ii, jj)],
      stringsAsFactors = FALSE
    )
  ))
}

#' Plot a BBNI MCMC Object
#'
#' Produces a trace plot of the log-posterior over MCMC iterations. This method
#' delegates to [plot_trace()]. For network visualization, use [plot_bbni()].
#'
#' @param x A `bbni` object returned by `run_bbni()`.
#' @param ... Extra arguments passed to [plot_trace()].
#' @return The plot object (invisibly).
#' @method plot bbni
#' @export
#' @examples
#' set.seed(123)
#' net <- GenerateNetwork(5)
#' data <- GenerateSample(net, 100)
#' res <- run_bbni(data, num_update = 50)
#' plot(res)
plot.bbni <- function(x, ...) {
  plot_trace(x, ...)
}
