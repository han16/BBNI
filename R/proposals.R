#' Pre-compute Weighted Proposal Distribution
#'
#' Pre-computes a weighted proposal distribution for the MCMC algorithm by
#' evaluating all possible 1-input and 2-input Boolean logic functions using
#' Bayesian Information Criterion (BIC). The inverse BIC values of the potential
#' candidates are used as selection weights for the MCMC sampler.
#'
#' @param GeneData A matrix of observational gene expression data.
#' @param SampleSize An integer representing the total number of time points in the dataset.
#'
#' @return A list containing two matrices (`Candidate[[1]]` for 2-input triplets and `Candidate[[2]]` for 1-input pairs). These matrices store the candidate parent nodes, the best-fitting Boolean logic function, the raw miscounts, and the calculated proposal selection weights.
#' @noRd
ProposalConstruction <- function(GeneData, SampleSize, timeseries = TRUE) {
  gene.data <- GeneData
  error.prop <- 0.4
  pseudo.count <- 0.01
  SampleSize <- SampleSize + pseudo.count * 8
  threshold <- SampleSize * error.prop
  if (timeseries) {
    idx_parent <- 1:(SampleSize - 1)
    idx_child <- 2:SampleSize
    child_offset <- 1
  } else {
    idx_parent <- 1:SampleSize
    idx_child <- 1:SampleSize
    child_offset <- 0
  }
  candidate.prior <- list()
  kk <- 1
  for (i in seq_len(nrow(gene.data))) {
    for (j in seq_len(nrow(gene.data))) {
      if (j != i) {
        for (k in seq_len(nrow(gene.data))) {
          if (k != i && k != j) {
            test.result <- list()
            gene.triplet <- rbind(rbind(gene.data[i, idx_parent], gene.data[j, idx_parent]), gene.data[k, idx_child])
            c000 <- 0
            c001 <- 0
            c010 <- 0
            c011 <- 0
            c100 <- 0
            c101 <- 0
            c110 <- 0
            c111 <- 0
            for (ii in seq_len(ncol(gene.triplet))) { # counts in each cell
              if (gene.triplet[1, ii] == 0 && gene.triplet[2, ii] == 0 && gene.triplet[3, ii] == 0) {
                c000 <- c000 + 1
              }
              if (gene.triplet[1, ii] == 0 && gene.triplet[2, ii] == 0 && gene.triplet[3, ii] == 1) {
                c001 <- c001 + 1
              }
              if (gene.triplet[1, ii] == 0 && gene.triplet[2, ii] == 1 && gene.triplet[3, ii] == 0) {
                c010 <- c010 + 1
              }
              if (gene.triplet[1, ii] == 0 && gene.triplet[2, ii] == 1 && gene.triplet[3, ii] == 1) {
                c011 <- c011 + 1
              }
              if (gene.triplet[1, ii] == 1 && gene.triplet[2, ii] == 0 && gene.triplet[3, ii] == 0) {
                c100 <- c100 + 1
              }
              if (gene.triplet[1, ii] == 1 && gene.triplet[2, ii] == 0 && gene.triplet[3, ii] == 1) {
                c101 <- c101 + 1
              }
              if (gene.triplet[1, ii] == 1 && gene.triplet[2, ii] == 1 && gene.triplet[3, ii] == 0) {
                c110 <- c110 + 1
              }
              if (gene.triplet[1, ii] == 1 && gene.triplet[2, ii] == 1 && gene.triplet[3, ii] == 1) {
                c111 <- c111 + 1
              }
            }
            test.stat <- c(c000, c001, c010, c011, c100, c101, c110, c111)
            test.result[[1]] <- c(i, j, k, 1, BF1(test.stat, pseudo.count, SampleSize, threshold))
            test.result[[2]] <- c(i, j, k, 2, BF2(test.stat, pseudo.count, SampleSize, threshold))
            test.result[[3]] <- c(i, j, k, 3, BF3(test.stat, pseudo.count, SampleSize, threshold))
            test.result[[4]] <- c(i, j, k, 4, BF4(test.stat, pseudo.count, SampleSize, threshold))
            test.result[[5]] <- c(i, j, k, 5, BF5(test.stat, pseudo.count, SampleSize, threshold))
            test.result[[6]] <- c(i, j, k, 6, BF6(test.stat, pseudo.count, SampleSize, threshold))
            test.result[[7]] <- c(i, j, k, 7, BF7(test.stat, pseudo.count, SampleSize, threshold))
            test.result[[8]] <- c(i, j, k, 8, BF8(test.stat, pseudo.count, SampleSize, threshold))
            test.result[[9]] <- c(i, j, k, 9, BF9(test.stat, pseudo.count, SampleSize, threshold))
            test.result[[10]] <- c(i, j, k, 10, BF10(test.stat, pseudo.count, SampleSize, threshold))
            test.result[[11]] <- c(i, j, k, 11, BF11(test.stat, pseudo.count, SampleSize, threshold)) # model g_k=g_i
            test.result[[12]] <- c(i, j, k, 12, BF12(test.stat, pseudo.count, SampleSize, threshold)) # model g_k=g_j
            test.result[[13]] <- c(i, j, k, 13, BF13(test.stat, pseudo.count, SampleSize, threshold)) # model g_k=complement(g_i)
            test.result[[14]] <- c(i, j, k, 14, BF14(test.stat, pseudo.count, SampleSize, threshold)) # model g_k=complement(g_j)
            miscount <- numeric()
            jj <- 1
            for (ii in seq_along(test.result)) {
              if (length(test.result[[ii]]) == 5) {
                miscount[jj] <- test.result[[ii]][5]
                jj <- jj + 1
              }
            }

            if (length(miscount) > 0) {
              min.miscount <- min(miscount)
              for (ii in seq_along(test.result)) {
                if (length(test.result[[ii]]) == 5 && test.result[[ii]][5] == min.miscount) {
                  candidate.prior[[kk]] <- test.result[[ii]]
                  kk <- kk + 1
                }
              }
            }
          }
        }
      }
    }
  }
  candidate <- matrix(nrow = length(candidate.prior), ncol = length(candidate.prior[[1]]))
  for (i in seq_along(candidate.prior)) {
    candidate[i, ] <- candidate.prior[[i]]
  }
  order.candidate <- candidate[order(candidate[, 3]), ] # order by output variables
  num.node <- max(order.candidate[, 3])
  triplet <- list()
  j1 <- 1
  pairwise <- list()
  j2 <- 1
  for (j in seq_len(nrow(order.candidate))) {
    if (order.candidate[j, 4] <= 10) {
      triplet[[j1]] <- order.candidate[j, ]
      j1 <- j1 + 1
    }
    if (order.candidate[j, 4] > 10) {
      if (order.candidate[j, 4] == 11) {
        pairwise[[j2]] <- c(order.candidate[j, 1], order.candidate[j, 3], 11, order.candidate[j, 5])
        j2 <- j2 + 1
      }
      if (order.candidate[j, 4] == 12) {
        pairwise[[j2]] <- c(order.candidate[j, 2], order.candidate[j, 3], 11, order.candidate[j, 5]) # model 11, g_k=g_i
        j2 <- j2 + 1
      }
      if (order.candidate[j, 4] == 13) {
        pairwise[[j2]] <- c(order.candidate[j, 1], order.candidate[j, 3], 12, order.candidate[j, 5])
        j2 <- j2 + 1
      }
      if (order.candidate[j, 4] == 14) {
        pairwise[[j2]] <- c(order.candidate[j, 2], order.candidate[j, 3], 12, order.candidate[j, 5]) # model 12: complement relation
        j2 <- j2 + 1
      }
    }
  }
  candidate.triplet <- matrix(nrow = length(triplet), ncol = 5)
  weighted.triplet <- matrix(nrow = length(triplet), ncol = 6)
  candidate.pairwise <- matrix(nrow = length(pairwise), ncol = 4)
  for (i in seq_along(triplet)) {
    candidate.triplet[i, ] <- triplet[[i]]
    for (j in 1:5) {
      weighted.triplet[i, j] <- candidate.triplet[i, j]
    }
    weighted.triplet[i, 6] <- candidate.triplet[i, 5]
  }
  for (i in seq_along(pairwise)) {
    candidate.pairwise[i, ] <- pairwise[[i]]
  }
  unique.pairwise <- unique(candidate.pairwise)
  weighted.pairwise <- matrix(nrow = length(unique.pairwise), ncol = 5)
  weighted.pairwise <- cbind(unique.pairwise, unique.pairwise[, 4])
  constant <- 0.001
  for (i in 1:num.node) {
    if (i %in% candidate.triplet[, 3]) {
      for (j in seq_len(nrow(candidate.triplet))) {
        if (candidate.triplet[j, 3] == i) {
          trip.matrix <- candidate.triplet[candidate.triplet[, 3] == i, 1:5]
          trip.matrix <- data.matrix(trip.matrix)
          if (ncol(trip.matrix) == 1) {
            trip.matrix <- t(trip.matrix)
          }
          score <- 1 / trip.matrix[, 5]
          prop <- (score + constant) / (sum(score + constant))
          trip.matrix[, 5] <- prop
          weighted.triplet[weighted.triplet[, 3] == i, ] <- cbind(trip.matrix, weighted.triplet[weighted.triplet[, 3] == i, 6])
        }
      }
    }
  }

  for (i in 1:num.node) {
    if (i %in% unique.pairwise[, 2]) {
      for (j in seq_len(nrow(unique.pairwise))) {
        if (unique.pairwise[j, 2] == i) {
          pair.matrix <- unique.pairwise[unique.pairwise[, 2] == i, 1:4]
          pair.matrix <- data.matrix(pair.matrix)
          if (ncol(pair.matrix) == 1) {
            pair.matrix <- t(pair.matrix)
          }
          score <- 1 / (pair.matrix[, 4] + constant)
          prop <- score / (sum(score))
          pair.matrix[, 4] <- prop
          weighted.pairwise[weighted.pairwise[, 2] == i, ] <- cbind(pair.matrix, weighted.pairwise[weighted.pairwise[, 2] == i, 5])
        }
      }
    }
  }
  CandidateTriplet <- weighted.triplet
  CandidatePairwise <- weighted.pairwise
  Candidate <- list()
  Candidate[[1]] <- CandidateTriplet
  Candidate[[2]] <- CandidatePairwise
  return(Candidate)
}

#' Build Initial Network Topology and Boolean Transition Matrix
#'
#' Builds a single valid starting network topology and Boolean transition function
#' matrix for the MCMC algorithm. It assigns parents using a 50/50 probability
#' of one versus two parents per node, then selects the most probable parent-child
#' relationships using pre-computed proposal weights. The function strictly ensures
#' the directed acyclic graph (DAG) constraint is maintained.
#'
#' @param Candidate A list containing transition function selection weights for 2-input and 1-input node combinations.
#' @param num.node An integer representing the total number of network nodes.
#'
#' @return A square transition function matrix representing the initial state of the network \eqn{T^{(0)}}{T^(0)} and its corresponding Boolean rules \eqn{F^{(0)}}{F^(0)}, providing an optimal acyclic starting network for the MCMC.
#' @importFrom stats runif
#' @noRd
ConstructInitial <- function(Candidate, num.node) {
  prior.triplet <- Candidate[[1]]
  prior.pairwise <- Candidate[[2]]
  trans_func_matrix <- matrix(0, nrow = num.node, ncol = num.node)
  prop_trans_func_matrix <- matrix(0, nrow = num.node, ncol = num.node)
  prop_incid_matrix <- matrix(0, nrow = num.node, ncol = num.node)
  incid_matrix <- matrix(0, nrow = num.node, ncol = num.node)
  ratio <- 0.5
  for (i in 1:num.node) {
    prop_trans_func_matrix <- trans_func_matrix
    prop_incid_matrix <- incid_matrix
    pairwise.prior.set <- matrix()
    pairwise.prior.set <- prior.pairwise[prior.pairwise[, 2] == i, ] # priors for node update_order[k]
    pairwise.prior.set <- data.matrix(pairwise.prior.set)
    if (ncol(pairwise.prior.set) == 1) {
      pairwise.prior.set <- t(pairwise.prior.set)
    }

    triplet.prior.set <- matrix()
    triplet.prior.set <- prior.triplet[prior.triplet[, 3] == i, ]
    triplet.prior.set <- data.matrix(triplet.prior.set)
    if (ncol(triplet.prior.set) == 1) {
      triplet.prior.set <- t(triplet.prior.set)
    }
    prop_prob <- runif(1)
    aa <- 0
    if (prop_prob >= ratio) { # consider two parents
      if (nrow(triplet.prior.set) >= 1) {
        candidate.pare.set <- list()
        for (jj in seq_len(nrow(triplet.prior.set))) {
          if (length(triplet.prior.set[jj, 1:2]) == 2) {
            aa <- aa + 1
            candidate.pare.set[[aa]] <- triplet.prior.set[jj, ]
          }
        }
        if (length(candidate.pare.set) > 0) {
          score <- numeric()
          for (jj in seq_along(candidate.pare.set)) {
            score[jj] <- candidate.pare.set[[jj]][5]
          }
          bbb <- matrix(nrow = length(candidate.pare.set), ncol = 6)
          for (jjj in seq_along(candidate.pare.set)) {
            bbb[jjj, ] <- candidate.pare.set[[jjj]]
          }

          max.score.candidate <- bbb[bbb[, 5] == max(score), ]
          if (is.numeric(max.score.candidate)) {
            add_parent <- max.score.candidate
          }
          if (is.matrix(max.score.candidate) && nrow(max.score.candidate) > 1) {
            add_parent <- max.score.candidate[sample.int(nrow(max.score.candidate), 1), ]
          }

          add_two_parent <- c(add_parent[1], add_parent[2])
        }
      }
      if (aa > 0) {
        func_order <- add_parent[4]
        prop_trans_func_matrix[i, add_two_parent[1]] <- func_order
        prop_trans_func_matrix[i, add_two_parent[2]] <- func_order
        for (ii in seq_len(nrow(prop_trans_func_matrix))) {
          for (jj in seq_len(ncol(prop_trans_func_matrix))) {
            if (prop_trans_func_matrix[ii, jj] > 0) {
              prop_incid_matrix[ii, jj] <- 1
            }
          }
        }
        prop_ances_matrix <- update_ancestor_matrix(prop_incid_matrix)
        if (check_ances_matrix(prop_ances_matrix) == 0) {
          trans_func_matrix[i, ] <- prop_trans_func_matrix[i, ]
          incid_matrix[i, ] <- prop_incid_matrix[i, ]
        }
      }
    }
    if (prop_prob < ratio && nrow(pairwise.prior.set) > 0) { # consider one parent
      candidate.pare.set <- pairwise.prior.set
      if (is.numeric(candidate.pare.set)) {
        add_parent <- candidate.pare.set
      }
      if (is.matrix(candidate.pare.set)) {
        add_parent <- candidate.pare.set[candidate.pare.set[, 4] == max(candidate.pare.set[, 4]), ]
      }
      if (is.matrix(add_parent) && nrow(add_parent) > 1) {
        add_parent <- add_parent[sample.int(nrow(add_parent), 1), ]
      }

      add_one_parent <- add_parent[1]

      func_order <- add_parent[3]
      prop_trans_func_matrix[i, add_one_parent] <- func_order
      for (ii in seq_len(nrow(prop_trans_func_matrix))) {
        for (jj in seq_len(ncol(prop_trans_func_matrix))) {
          if (prop_trans_func_matrix[ii, jj] > 0) {
            prop_incid_matrix[ii, jj] <- 1
          }
        }
      }
      prop_ances_matrix <- update_ancestor_matrix(prop_incid_matrix)
      if (check_ances_matrix(prop_ances_matrix) == 0) {
        trans_func_matrix[i, ] <- prop_trans_func_matrix[i, ]
        incid_matrix[i, ] <- prop_incid_matrix[i, ]
      }
    }
  }
  return(trans_func_matrix)
}
