#' Calculate Collapsed Posterior Log-Probability
#'
#' Calculates the collapsed posterior log-probability of a proposed network
#' topology (\eqn{T}{T}) and its Boolean functions (\eqn{F}{F}) given the observed data (\eqn{G}{G}).
#' It parses the transition function matrix (TRFUM) into root and non-root nodes,
#' predicts each node's value based on its Boolean rule, counts mismatches, and
#' evaluates the Beta-Binomial collapsed likelihood. By integrating out Bernoulli
#' parameters \eqn{\theta = \{p_E, p_i\}}{\theta = {p_E, p_i}} under Beta priors, the function computes
#' \eqn{\log p(T, F \mid G)}{\log p(T, F | G)}, the model score for the MCMC algorithm.
#'
#' @param TRFUM A matrix combining the current network topology and the specific Boolean logic functions assigned to each node (integer codes 1-12).
#' @param GeneData A matrix of the observed binary gene expression time-series dataset (rows = nodes, columns = time points).
#' @param SampleSize An integer representing the total number of time points/samples in `GeneData`.
#' @param num.node An integer representing the total number of genes/nodes in the network.
#' @param prior_para A matrix of Beta distribution hyperparameters (\eqn{\alpha, \beta}{\alpha, \beta}) for root nodes and the global noise parameter \eqn{e}{e}.
#' @param penalty A numeric value representing the structural prior probability per edge used to penalize network complexity \eqn{P(T)}{P(T)}.
#'
#' @return A list of two vectors evaluating model fit. `results[[1]]` contains `[ErrorFactor, RootFactor, likelihood, post_para, log_post_model]`, where `log_post_model` is the collapsed posterior metric \eqn{\log p(T, F \mid G)}{\log p(T, F | G)}. `results[[2]]` contains `[para_sample, mismatch, Perror]`, providing point estimates of the root ON-probabilities, total mismatches, and estimated noise error rate \eqn{e}{e}.
#' @noRd
Error_LLH <- function(TRFUM, GeneData, SampleSize, num.node, prior_para, penalty, timeseries = TRUE) {
  if (timeseries) {
    idx_child <- 2:SampleSize
    idx_parent <- 1:(SampleSize - 1)
  } else {
    idx_child <- 1:SampleSize
    idx_parent <- 1:SampleSize
  }
  InOutPair <- list()
  ii <- 0
  root_node <- numeric()
  root <- 0
  for (i in seq_len(nrow(TRFUM))) { #  i row   find each input output combination
    if (sum(TRFUM[i, ]) > 0) {
      ii <- ii + 1
      p <- numeric()
      k <- 1
      jj <- numeric()
      for (j in seq_len(nrow(TRFUM))) { # j column
        if (TRFUM[i, j] != 0) {
          p[k] <- TRFUM[i, j]
          jj[k] <- j
          InOutPair[[ii]] <- c(i, jj, p[k])
          k <- k + 1
        }
      }
    }
    if (sum(TRFUM[i, ]) == 0) {
      root <- root + 1
      root_node[root] <- i
    }
  }
  mismatch <- 0
  if (length(InOutPair) > 0) {
    for (i in seq_along(InOutPair)) {
      if (InOutPair[[i]][length(InOutPair[[i]])] <= 10) {
        in_node <- c(InOutPair[[i]][2], InOutPair[[i]][3])
        out_node <- InOutPair[[i]][1]
        if (InOutPair[[i]][length(InOutPair[[i]])] == 1) {
          mismatch <- mismatch + sum(bitXor(GeneData[out_node, idx_child], bitAnd(GeneData[in_node[1], idx_parent], GeneData[in_node[2], idx_parent])))
        }
        if (InOutPair[[i]][length(InOutPair[[i]])] == 2) {
          mismatch <- mismatch + sum(bitXor(GeneData[out_node, idx_child], 1 - bitAnd(GeneData[in_node[1], idx_parent], GeneData[in_node[2], idx_parent])))
        }
        if (InOutPair[[i]][length(InOutPair[[i]])] == 3) {
          mismatch <- mismatch + sum(bitXor(GeneData[out_node, idx_child], bitOr(GeneData[in_node[1], idx_parent], GeneData[in_node[2], idx_parent])))
        }
        if (InOutPair[[i]][length(InOutPair[[i]])] == 4) {
          mismatch <- mismatch + sum(bitXor(GeneData[out_node, idx_child], 1 - bitOr(GeneData[in_node[1], idx_parent], GeneData[in_node[2], idx_parent])))
        }
        if (InOutPair[[i]][length(InOutPair[[i]])] == 5) {
          mismatch <- mismatch + sum(bitXor(GeneData[out_node, idx_child], bitOr(1 - GeneData[in_node[1], idx_parent], GeneData[in_node[2], idx_parent])))
        }
        if (InOutPair[[i]][length(InOutPair[[i]])] == 6) {
          mismatch <- mismatch + sum(bitXor(GeneData[out_node, idx_child], bitOr(GeneData[in_node[1], idx_parent], 1 - GeneData[in_node[2], idx_parent])))
        }
        if (InOutPair[[i]][length(InOutPair[[i]])] == 7) {
          mismatch <- mismatch + sum(bitXor(GeneData[out_node, idx_child], bitAnd(1 - GeneData[in_node[1], idx_parent], GeneData[in_node[2], idx_parent])))
        }
        if (InOutPair[[i]][length(InOutPair[[i]])] == 8) {
          mismatch <- mismatch + sum(bitXor(GeneData[out_node, idx_child], bitAnd(GeneData[in_node[1], idx_parent], 1 - GeneData[in_node[2], idx_parent])))
        }
        if (InOutPair[[i]][length(InOutPair[[i]])] == 9) {
          mismatch <- mismatch + sum(bitXor(GeneData[out_node, idx_child], bitXor(GeneData[in_node[1], idx_parent], GeneData[in_node[2], idx_parent])))
        }
        if (InOutPair[[i]][length(InOutPair[[i]])] == 10) {
          mismatch <- mismatch + sum(bitXor(GeneData[out_node, idx_child], 1 - bitXor(GeneData[in_node[1], idx_parent], GeneData[in_node[2], idx_parent])))
        }
      }
      if (InOutPair[[i]][length(InOutPair[[i]])] > 10) {
        in_node <- InOutPair[[i]][2]
        out_node <- InOutPair[[i]][1]
        if (InOutPair[[i]][length(InOutPair[[i]])] == 11) {
          mismatch <- mismatch + sum(bitXor(GeneData[out_node, idx_child], GeneData[in_node, idx_parent] ))
        }
        if (InOutPair[[i]][length(InOutPair[[i]])] == 12) {
          mismatch <- mismatch + sum(bitXor(GeneData[out_node, idx_child], 1 - GeneData[in_node, idx_parent] ))
        }
      }
    }
  }
  pseudo_count <- 0.0001
  mismatch <- mismatch + pseudo_count
  Perror <- mismatch / (ii * SampleSize + pseudo_count)
  ErrorFactor <- mismatch * log(Perror) + (ii * SampleSize - mismatch + pseudo_count) * log(1 - Perror)
  ErrorPrior <- (prior_para[num.node + 1, 1] - 1) * log(Perror) + (prior_para[num.node + 1, 2] - 1) * log(1 - Perror)
  RootFactor <- numeric()
  RootPrior <- numeric()
  succ_count <- numeric()
  succ_prob <- numeric()
  for (i in seq_along(root_node)) {
    succ_count[i] <- sum(GeneData[root_node[i], ]) + pseudo_count
    succ_prob[i] <- succ_count[i] / (SampleSize + pseudo_count)
    RootFactor[i] <- succ_count[i] * log(succ_prob[i]) + (SampleSize + pseudo_count - succ_count[i]) * log(1 - succ_prob[i])

    RootPrior[i] <- (prior_para[root_node[i], 1] - 1) * log(succ_prob[i]) + (prior_para[root_node[i], 2] - 1) * log(1 - succ_prob[i])
  }
  if (length(InOutPair) == 0) { # all nodes are root nodes
    likelihood <- sum(RootFactor)
    post_para <- sum(RootFactor) + sum(RootPrior)
    log_post_model <- 0
    for (i in seq_len(nrow(TRFUM))) {
      nume <- lbeta(prior_para[i, 1] + sum(GeneData[i, ]), prior_para[i, 2] + SampleSize - sum(GeneData[i, ]))
      deno <- lbeta(prior_para[i, 1], prior_para[i, 2])
      log_post_model <- log_post_model + nume - deno
    }
    log_post_model <- log_post_model + length(TRFUM[TRFUM > 0]) * log(penalty)
    ErrorFactor <- NA
    Perror <- NA
    mismatch <- NA
  }
  if (length(InOutPair) > 0) { # exist non root nodes
    likelihood <- ErrorFactor + sum(RootFactor)
    post_para <- ErrorFactor + sum(RootFactor) + sum(RootPrior) + ErrorPrior
    log_post_model <- 0
    for (i in seq_along(root_node)) {
      index <- root_node[i]
      nume <- lbeta(prior_para[index, 1] + sum(GeneData[index, ]), prior_para[index, 2] + SampleSize - sum(GeneData[index, ]))
      deno <- lbeta(prior_para[index, 1], prior_para[index, 2])
      log_post_model <- log_post_model + nume - deno
    }
    noise_nume <- lbeta(mismatch + prior_para[num.node + 1, 1], length(InOutPair) * SampleSize - mismatch + prior_para[num.node + 1, 2])
    noise_deno <- lbeta(prior_para[num.node + 1, 1], prior_para[num.node + 1, 2])
    log_post_model <- log_post_model + noise_nume - noise_deno
    log_post_model <- log_post_model + length(TRFUM[TRFUM > 0]) * log(penalty)
  }
  result <- list()
  result[[1]] <- c(ErrorFactor, RootFactor, likelihood, post_para, log_post_model)
  para_sample <- rep(NA, num.node)
  for (i in 1:num.node) {
    if (i %in% root_node) {
      for (j in seq_along(root_node)) {
        if (i == root_node[j]) {
          para_sample[i] <- succ_prob[j]
        }
      }
    }
  }
  result[[2]] <- c(para_sample, mismatch, Perror)
  return(result)
}

#' Evaluate Boolean Function 1 (AND)
#'
#' Computes a Bayesian Information Criterion (BIC) score for a candidate
#' Boolean function by counting state mismatches between predicted and
#' observed outputs. It compresses the raw data into a compact 3-gene
#' frequency count (truth table) to efficiently evaluate matching patterns.
#'
#' @param test.stat An integer frequency count vector storing how many times each specific combination of parent and child states occurs together in the data.
#' @param pseudo.count A numeric pseudocount added to avoid zero cells.
#' @param SampleSize An integer representing the total sample size.
#' @param threshold A numeric rejection threshold. If the error rate exceeds this, the function is rejected.
#'
#' @return A numeric BIC value, or `NULL` if the mismatch count exceeds the threshold.
#' @noRd
BF1 <- function(test.stat, pseudo.count, SampleSize, threshold) { # model g_k=g_i and g_j
  test.stat <- test.stat + pseudo.count
  false.count <- sum(test.stat[2], test.stat[4], test.stat[6], test.stat[7])
  error.estimate <- false.count / (SampleSize + pseudo.count * 8)
  BIC.value <- -2 * (false.count * log(error.estimate) + (SampleSize + pseudo.count * 8 - false.count) * log(1 - error.estimate)) + 2 * log(SampleSize)
  if (false.count <= threshold) {
    return(BIC.value)
  }
}

#' Evaluate Boolean Function 2 (NAND)
#'
#' Computes a Bayesian Information Criterion (BIC) score for a candidate
#' Boolean function by counting state mismatches between predicted and
#' observed outputs. It compresses the raw data into a compact 3-gene
#' frequency count (truth table) to efficiently evaluate matching patterns.
#'
#' @param test.stat An integer frequency count vector storing how many times each specific combination of parent and child states occurs together in the data.
#' @param pseudo.count A numeric pseudocount added to avoid zero cells.
#' @param SampleSize An integer representing the total sample size.
#' @param threshold A numeric rejection threshold. If the error rate exceeds this, the function is rejected.
#'
#' @return A numeric BIC value, or `NULL` if the mismatch count exceeds the threshold.
#' @noRd
BF2 <- function(test.stat, pseudo.count, SampleSize, threshold) { # model g_k=complement(g_i and g_j)
  test.stat <- test.stat + pseudo.count
  false.count <- sum(test.stat[1], test.stat[3], test.stat[5], test.stat[8])
  error.estimate <- false.count / (SampleSize + pseudo.count * 8)
  BIC.value <- -2 * (false.count * log(error.estimate) + (SampleSize + pseudo.count * 8 - false.count) * log(1 - error.estimate)) + 2 * log(SampleSize)
  if (false.count <= threshold) {
    return(BIC.value)
  }
}

#' Evaluate Boolean Function 3 (OR)
#'
#' Computes a Bayesian Information Criterion (BIC) score for a candidate
#' Boolean function by counting state mismatches between predicted and
#' observed outputs. It compresses the raw data into a compact 3-gene
#' frequency count (truth table) to efficiently evaluate matching patterns.
#'
#' @param test.stat An integer frequency count vector storing how many times each specific combination of parent and child states occurs together in the data.
#' @param pseudo.count A numeric pseudocount added to avoid zero cells.
#' @param SampleSize An integer representing the total sample size.
#' @param threshold A numeric rejection threshold. If the error rate exceeds this, the function is rejected.
#'
#' @return A numeric BIC value, or `NULL` if the mismatch count exceeds the threshold.
#' @noRd
BF3 <- function(test.stat, pseudo.count, SampleSize, threshold) { # model g_k=(g_i or g_j)
  test.stat <- test.stat + pseudo.count
  false.count <- sum(test.stat[2], test.stat[3], test.stat[5], test.stat[7])
  error.estimate <- false.count / (SampleSize + pseudo.count * 8)
  BIC.value <- -2 * (false.count * log(error.estimate) + (SampleSize + pseudo.count * 8 - false.count) * log(1 - error.estimate)) + 2 * log(SampleSize)
  if (false.count <= threshold) {
    return(BIC.value)
  }
}

#' Evaluate Boolean Function 4 (NOR)
#'
#' Computes a Bayesian Information Criterion (BIC) score for a candidate
#' Boolean function by counting state mismatches between predicted and
#' observed outputs. It compresses the raw data into a compact 3-gene
#' frequency count (truth table) to efficiently evaluate matching patterns.
#'
#' @param test.stat An integer frequency count vector storing how many times each specific combination of parent and child states occurs together in the data.
#' @param pseudo.count A numeric pseudocount added to avoid zero cells.
#' @param SampleSize An integer representing the total sample size.
#' @param threshold A numeric rejection threshold. If the error rate exceeds this, the function is rejected.
#'
#' @return A numeric BIC value, or `NULL` if the mismatch count exceeds the threshold.
#' @noRd
BF4 <- function(test.stat, pseudo.count, SampleSize, threshold) { # model g_k=complement(g_i or g_j)
  test.stat <- test.stat + pseudo.count
  false.count <- sum(test.stat[1], test.stat[4], test.stat[6], test.stat[8])
  error.estimate <- false.count / (SampleSize + pseudo.count * 8)
  BIC.value <- -2 * (false.count * log(error.estimate) + (SampleSize + pseudo.count * 8 - false.count) * log(1 - error.estimate)) + 2 * log(SampleSize)
  if (false.count <= threshold) {
    return(BIC.value)
  }
}

#' Evaluate Boolean Function 5 (OR-NOT)
#'
#' Computes a Bayesian Information Criterion (BIC) score for a candidate
#' Boolean function by counting state mismatches between predicted and
#' observed outputs. It compresses the raw data into a compact 3-gene
#' frequency count (truth table) to efficiently evaluate matching patterns.
#'
#' @param test.stat An integer frequency count vector storing how many times each specific combination of parent and child states occurs together in the data.
#' @param pseudo.count A numeric pseudocount added to avoid zero cells.
#' @param SampleSize An integer representing the total sample size.
#' @param threshold A numeric rejection threshold. If the error rate exceeds this, the function is rejected.
#'
#' @return A numeric BIC value, or `NULL` if the mismatch count exceeds the threshold.
#' @noRd
BF5 <- function(test.stat, pseudo.count, SampleSize, threshold) { # model g_k=complement(g_i) or g_j
  test.stat <- test.stat + pseudo.count
  false.count <- sum(test.stat[1], test.stat[3], test.stat[6], test.stat[7])
  error.estimate <- false.count / (SampleSize + pseudo.count * 8)
  BIC.value <- -2 * (false.count * log(error.estimate) + (SampleSize + pseudo.count * 8 - false.count) * log(1 - error.estimate)) + 2 * log(SampleSize)
  if (false.count <= threshold) {
    return(BIC.value)
  }
}

#' Evaluate Boolean Function 6 (NOT-OR)
#'
#' Computes a Bayesian Information Criterion (BIC) score for a candidate
#' Boolean function by counting state mismatches between predicted and
#' observed outputs. It compresses the raw data into a compact 3-gene
#' frequency count (truth table) to efficiently evaluate matching patterns.
#'
#' @param test.stat An integer frequency count vector storing how many times each specific combination of parent and child states occurs together in the data.
#' @param pseudo.count A numeric pseudocount added to avoid zero cells.
#' @param SampleSize An integer representing the total sample size.
#' @param threshold A numeric rejection threshold. If the error rate exceeds this, the function is rejected.
#'
#' @return A numeric BIC value, or `NULL` if the mismatch count exceeds the threshold.
#' @noRd
BF6 <- function(test.stat, pseudo.count, SampleSize, threshold) { # model g_k=g_i or complement(g_j)
  test.stat <- test.stat + pseudo.count
  false.count <- sum(test.stat[1], test.stat[4], test.stat[5], test.stat[7])
  error.estimate <- false.count / (SampleSize + pseudo.count * 8)
  BIC.value <- -2 * (false.count * log(error.estimate) + (SampleSize + pseudo.count * 8 - false.count) * log(1 - error.estimate)) + 2 * log(SampleSize)
  if (false.count <= threshold) {
    return(BIC.value)
  }
}

#' Evaluate Boolean Function 7 (AND-NOT)
#'
#' Computes a Bayesian Information Criterion (BIC) score for a candidate
#' Boolean function by counting state mismatches between predicted and
#' observed outputs. It compresses the raw data into a compact 3-gene
#' frequency count (truth table) to efficiently evaluate matching patterns.
#'
#' @param test.stat An integer frequency count vector storing how many times each specific combination of parent and child states occurs together in the data.
#' @param pseudo.count A numeric pseudocount added to avoid zero cells.
#' @param SampleSize An integer representing the total sample size.
#' @param threshold A numeric rejection threshold. If the error rate exceeds this, the function is rejected.
#'
#' @return A numeric BIC value, or `NULL` if the mismatch count exceeds the threshold.
#' @noRd
BF7 <- function(test.stat, pseudo.count, SampleSize, threshold) { # model g_k=complement(g_i) and g_j
  test.stat <- test.stat + pseudo.count
  false.count <- sum(test.stat[2], test.stat[3], test.stat[6], test.stat[8])
  error.estimate <- false.count / (SampleSize + pseudo.count * 8)
  BIC.value <- -2 * (false.count * log(error.estimate) + (SampleSize + pseudo.count * 8 - false.count) * log(1 - error.estimate)) + 2 * log(SampleSize)
  if (false.count <= threshold) {
    return(BIC.value)
  }
}

#' Evaluate Boolean Function 8 (NOT-AND)
#'
#' Computes a Bayesian Information Criterion (BIC) score for a candidate
#' Boolean function by counting state mismatches between predicted and
#' observed outputs. It compresses the raw data into a compact 3-gene
#' frequency count (truth table) to efficiently evaluate matching patterns.
#'
#' @param test.stat An integer frequency count vector storing how many times each specific combination of parent and child states occurs together in the data.
#' @param pseudo.count A numeric pseudocount added to avoid zero cells.
#' @param SampleSize An integer representing the total sample size.
#' @param threshold A numeric rejection threshold. If the error rate exceeds this, the function is rejected.
#'
#' @return A numeric BIC value, or `NULL` if the mismatch count exceeds the threshold.
#' @noRd
BF8 <- function(test.stat, pseudo.count, SampleSize, threshold) { # model g_k=g_i and complement(g_j)
  test.stat <- test.stat + pseudo.count
  false.count <- sum(test.stat[2], test.stat[4], test.stat[5], test.stat[8])
  error.estimate <- false.count / (SampleSize + pseudo.count * 8)
  BIC.value <- -2 * (false.count * log(error.estimate) + (SampleSize + pseudo.count * 8 - false.count) * log(1 - error.estimate)) + 2 * log(SampleSize)
  if (false.count <= threshold) {
    return(BIC.value)
  }
}

#' Evaluate Boolean Function 9 (XOR)
#'
#' Computes a Bayesian Information Criterion (BIC) score for a candidate
#' Boolean function by counting state mismatches between predicted and
#' observed outputs. It compresses the raw data into a compact 3-gene
#' frequency count (truth table) to efficiently evaluate matching patterns.
#'
#' @param test.stat An integer frequency count vector storing how many times each specific combination of parent and child states occurs together in the data.
#' @param pseudo.count A numeric pseudocount added to avoid zero cells.
#' @param SampleSize An integer representing the total sample size.
#' @param threshold A numeric rejection threshold. If the error rate exceeds this, the function is rejected.
#'
#' @return A numeric BIC value, or `NULL` if the mismatch count exceeds the threshold.
#' @noRd
BF9 <- function(test.stat, pseudo.count, SampleSize, threshold) { # model g_k=g_i xor g_j
  test.stat <- test.stat + pseudo.count
  false.count <- sum(test.stat[2], test.stat[3], test.stat[5], test.stat[8])
  error.estimate <- false.count / (SampleSize + pseudo.count * 8)
  BIC.value <- -2 * (false.count * log(error.estimate) + (SampleSize + pseudo.count * 8 - false.count) * log(1 - error.estimate)) + 2 * log(SampleSize)
  if (false.count <= threshold) {
    return(BIC.value)
  }
}

#' Evaluate Boolean Function 10 (NXOR)
#'
#' Computes a Bayesian Information Criterion (BIC) score for a candidate
#' Boolean function by counting state mismatches between predicted and
#' observed outputs. It compresses the raw data into a compact 3-gene
#' frequency count (truth table) to efficiently evaluate matching patterns.
#'
#' @param test.stat An integer frequency count vector storing how many times each specific combination of parent and child states occurs together in the data.
#' @param pseudo.count A numeric pseudocount added to avoid zero cells.
#' @param SampleSize An integer representing the total sample size.
#' @param threshold A numeric rejection threshold. If the error rate exceeds this, the function is rejected.
#'
#' @return A numeric BIC value, or `NULL` if the mismatch count exceeds the threshold.
#' @noRd
BF10 <- function(test.stat, pseudo.count, SampleSize, threshold) { # model g_k=complement(g_i xor g_j)
  test.stat <- test.stat + pseudo.count
  false.count <- sum(test.stat[1], test.stat[4], test.stat[6], test.stat[7])
  error.estimate <- false.count / (SampleSize + pseudo.count * 8)
  BIC.value <- -2 * (false.count * log(error.estimate) + (SampleSize + pseudo.count * 8 - false.count) * log(1 - error.estimate)) + 2 * log(SampleSize)
  if (false.count <= threshold) {
    return(BIC.value)
  }
}

#' Evaluate Boolean Pairwise Relation 11 (Identity First Parent)
#'
#' Computes a Bayesian Information Criterion (BIC) score for a candidate
#' Boolean function by counting state mismatches between predicted and
#' observed outputs. It compresses the raw data into a compact 3-gene
#' frequency count (truth table) to efficiently evaluate matching patterns.
#'
#' @param test.stat An integer frequency count vector storing how many times each specific combination of parent and child states occurs together in the data.
#' @param pseudo.count A numeric pseudocount added to avoid zero cells.
#' @param SampleSize An integer representing the total sample size.
#' @param threshold A numeric rejection threshold. If the error rate exceeds this, the function is rejected.
#'
#' @return A numeric BIC value, or `NULL` if the mismatch count exceeds the threshold.
#' @noRd
BF11 <- function(test.stat, pseudo.count, SampleSize, threshold) { # model g_k=g_i
  test.stat <- test.stat + pseudo.count
  false.count <- sum(test.stat[2], test.stat[4], test.stat[5], test.stat[7])
  error.estimate <- false.count / (SampleSize + pseudo.count * 8)
  BIC.value <- -2 * (false.count * log(error.estimate) + (SampleSize + pseudo.count * 8 - false.count) * log(1 - error.estimate)) + 2 * log(SampleSize)
  if (false.count <= threshold) {
    return(BIC.value)
  }
}

#' Evaluate Boolean Pairwise Relation 12 (Identity Second Parent)
#'
#' Computes a Bayesian Information Criterion (BIC) score for a candidate
#' Boolean function by counting state mismatches between predicted and
#' observed outputs. It compresses the raw data into a compact 3-gene
#' frequency count (truth table) to efficiently evaluate matching patterns.
#'
#' @param test.stat An integer frequency count vector storing how many times each specific combination of parent and child states occurs together in the data.
#' @param pseudo.count A numeric pseudocount added to avoid zero cells.
#' @param SampleSize An integer representing the total sample size.
#' @param threshold A numeric rejection threshold. If the error rate exceeds this, the function is rejected.
#'
#' @return A numeric BIC value, or `NULL` if the mismatch count exceeds the threshold.
#' @noRd
BF12 <- function(test.stat, pseudo.count, SampleSize, threshold) { # model g_k=g_j
  test.stat <- test.stat + pseudo.count
  false.count <- sum(test.stat[2], test.stat[3], test.stat[6], test.stat[7])
  error.estimate <- false.count / (SampleSize + pseudo.count * 8)
  BIC.value <- -2 * (false.count * log(error.estimate) + (SampleSize + pseudo.count * 8 - false.count) * log(1 - error.estimate)) + 2 * log(SampleSize)
  if (false.count <= threshold) {
    return(BIC.value)
  }
}

#' Evaluate Boolean Negation 13 (NOT First Parent)
#'
#' Computes a Bayesian Information Criterion (BIC) score for a candidate
#' Boolean function by counting state mismatches between predicted and
#' observed outputs. It compresses the raw data into a compact 3-gene
#' frequency count (truth table) to efficiently evaluate matching patterns.
#'
#' @param test.stat An integer frequency count vector storing how many times each specific combination of parent and child states occurs together in the data.
#' @param pseudo.count A numeric pseudocount added to avoid zero cells.
#' @param SampleSize An integer representing the total sample size.
#' @param threshold A numeric rejection threshold. If the error rate exceeds this, the function is rejected.
#'
#' @return A numeric BIC value, or `NULL` if the mismatch count exceeds the threshold.
#' @noRd
BF13 <- function(test.stat, pseudo.count, SampleSize, threshold) { # model g_k=complement(g_i)
  test.stat <- test.stat + pseudo.count
  false.count <- sum(test.stat[1], test.stat[3], test.stat[6], test.stat[8])
  error.estimate <- false.count / (SampleSize + pseudo.count * 8)
  BIC.value <- -2 * (false.count * log(error.estimate) + (SampleSize + pseudo.count * 8 - false.count) * log(1 - error.estimate)) + 2 * log(SampleSize)
  if (false.count <= threshold) {
    return(BIC.value)
  }
}

#' Evaluate Boolean Negation 14 (NOT Second Parent)
#'
#' Computes a Bayesian Information Criterion (BIC) score for a candidate
#' Boolean function by counting state mismatches between predicted and
#' observed outputs. It compresses the raw data into a compact 3-gene
#' frequency count (truth table) to efficiently evaluate matching patterns.
#'
#' @param test.stat An integer frequency count vector storing how many times each specific combination of parent and child states occurs together in the data.
#' @param pseudo.count A numeric pseudocount added to avoid zero cells.
#' @param SampleSize An integer representing the total sample size.
#' @param threshold A numeric rejection threshold. If the error rate exceeds this, the function is rejected.
#'
#' @return A numeric BIC value, or `NULL` if the mismatch count exceeds the threshold.
#' @noRd
BF14 <- function(test.stat, pseudo.count, SampleSize, threshold) { # model g_k=complement(g_j)
  test.stat <- test.stat + pseudo.count
  false.count <- sum(test.stat[1], test.stat[4], test.stat[5], test.stat[8])
  error.estimate <- false.count / (SampleSize + pseudo.count * 8)
  BIC.value <- -2 * (false.count * log(error.estimate) + (SampleSize + pseudo.count * 8 - false.count) * log(1 - error.estimate)) + 2 * log(SampleSize)
  if (false.count <= threshold) {
    return(BIC.value)
  }
}
