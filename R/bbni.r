#' Execute Metropolis-within-Gibbs MCMC Sampler for Boolean Networks
#'
#' Executes a Metropolis-within-Gibbs Markov chain Monte Carlo (MCMC) algorithm to sample
#' from the joint posterior distribution of Directed Acyclic Graph (DAG) topologies (\eqn{T}{T}) and
#' Boolean logic transition functions (\eqn{F}{F}). The algorithm iterates through individual
#' network nodes and proposes parent set mutations (edge additions, removals, or swaps)
#' paired with transition function reassignments to one of 14 candidate Boolean rules.
#' Proposed states transitions are strictly verified to follow the DAG constraint and
#' evaluated with a Metropolis-Hastings acceptance threshold using log-posterior values.
#'
#' @param GeneData A binary empirical observation matrix (\eqn{G}{G}), where rows represent individual network nodes (genes) and columns represent independent samples or sequential time points.
#' @param num.node An integer representing the total number of network nodes. Defaults to \code{nrow(GeneData)} if not specified.)
#' @param SampleSize An integer representing the total number of time points or independent samples in the dataset. Defaults to \code{ncol(GeneData)} if not specified.
#' @param prior_para A matrix (with dimensions `(num.node + 1) x 2`) of Beta prior hyperparameters \eqn{\alpha}{\alpha} and \eqn{\beta}{\beta} for root node probabilities and the global noise parameter \eqn{e}{e}. Defaults to a flat prior if not specified.
#' @param num_update An integer representing the total number of MCMC iterations to perform. Defaults to 4000 if not specified.
#' @param penalty A numeric value representing the structural prior probability per edge used to penalize network complexity \eqn{P(T)}{P(T)}. Defaults to 0.1 if not specified.
#' @param prop.ratio A numeric value between 0 and 1 representing the probability of choosing a uniform proposal distribution over an empirical proposal distribution at each iteration. Defaults to 0.5
#' @param verbose Logical. If TRUE, prints verbose MCMC iteration progress to the console. Default is FALSE.
#' @param timeseries Logical. If TRUE, the algorithm assumes a time-series dataset. If FALSE, the algorithm assumes independent samples. Default is TRUE.
#' @param burn_in A numeric value between 0 and 1 representing the proportion of initial MCMC samples to discard as burn-in. Defaults to 0.7 if not specified.
#'
#' @return A list containing the full trajectory of the MCMC chain. Specifically, `networks` (a list of sampled transition function matrices) and `log_posterior` (a numeric vector of log-posterior scores for each iteration). These represent samples drawn from the marginal posterior distribution \eqn{P(T,F|G)}{P(T,F|G)} used for Bayesian model averaging. Additionally, the `post_edge_prob` (matrix of marginal posterior edge probabilities) and `burn_in` ratio are returned in the list.
#'
#' @examples
#' \donttest{
#'   # 1. Define network parameters
#'   set.seed(235)
#'   num_nodes <- 10
#'   sample_size <- 50
#'
#'   # 2. Generate true network and simulate data
#'   true_network <- GenerateNetwork(num.node = num_nodes)
#'
#'   # Set up Beta priors for root-node probabilities and the noise rate
#'   prior_para <- matrix(3, nrow = num_nodes + 1, ncol = 2)
#'   prior_para[num_nodes + 1, 1] <- 2
#'   prior_para[num_nodes + 1, 2] <- 100
#'
#'   # Simulate parameters
#'   para <- numeric(num_nodes + 1)
#'   for (i in 1:(num_nodes + 1)) {
#'     para[i] <- stats::rbeta(1, prior_para[i, 1], prior_para[i, 2])
#'   }
#'   para[num_nodes + 1] <- 0.1 # Fixed noise rate for simulation
#'
#'   error_matrix <- matrix(stats::rbinom(num_nodes * sample_size, 1, para[num_nodes + 1]),
#'     nrow = num_nodes, ncol = sample_size
#'   )
#'
#'   dummy_data <- GenerateSample(
#'     trans_matrix = true_network,
#'     SampleSize = sample_size,
#'     para = para,
#'     error = error_matrix
#'   )
#'
#'   # 3. Run the MCMC sampler (silently)
#'   mcmc_results <- run_bbni(
#'     GeneData = dummy_data,
#'     prior_para = prior_para,
#'     num_update = 100, # Scaled down for example speed
#'     prop.ratio = 0.1
#'   )
#'
#'   # 4. Visualize results
#'   plot_bbni(mcmc_results)
#'   plot_trace(mcmc_results)
#' }
#'
#' @importFrom stats runif
#' @importFrom utils setTxtProgressBar txtProgressBar
#' @export
run_bbni <- function(GeneData, num.node = nrow(GeneData), SampleSize = ncol(GeneData), prior_para = NULL,
                     num_update = 4000, penalty = 0.1, prop.ratio = 0.5, verbose = FALSE, timeseries = TRUE, burn_in = 0.7) {
  # Generate flat priors if not provided by user
  if (is.null(prior_para)) {
    prior_para <- matrix(1, nrow = num.node + 1, ncol = 2)
  }
  Candidate <- ProposalConstruction(GeneData, SampleSize, timeseries) # create the proposal for generated data
  prior.triplet <- Candidate[[1]]
  prior.pairwise <- Candidate[[2]]

  trans_func_matrix <- ConstructInitial(Candidate, num.node) # use the randomly selected initial
  incid_matrix <- trans_func_matrix
  for (i in seq_len(nrow(incid_matrix))) {
    for (j in seq_len(ncol(incid_matrix))) {
      if (incid_matrix[i, j] > 0) {
        incid_matrix[i, j] <- 1
      }
    }
  }
  ances_matrix <- update_ancestor_matrix(incid_matrix)


  Incidence_Matrix <- list()
  Ancestor_Matrix <- list() # Matrix: Ancestor Matrix recording ancestor-offspring relations for the whole chain
  Trans_Func_Matrix <- list() # Matrix: Transition Function Matrix for the whole chain
  Sample_Matrix <- matrix(nrow = num.node * num_update + 1, ncol = num.node + 2)
  Sample_Matrix[1, ] <- Error_LLH(GeneData = GeneData, SampleSize = SampleSize, num.node = num.node, prior_para = prior_para, penalty = penalty, TRFUM = trans_func_matrix, timeseries = timeseries)[[2]]
  Incidence_Matrix[[1]] <- incid_matrix
  Ancestor_Matrix[[1]] <- ances_matrix
  Trans_Func_Matrix[[1]] <- trans_func_matrix
  num <- numeric()
  logpost <- numeric()
  all_logpost <- numeric()
  aa <- Error_LLH(GeneData = GeneData, SampleSize = SampleSize, num.node = num.node, prior_para = prior_para, penalty = penalty, TRFUM = Trans_Func_Matrix[[1]], timeseries = timeseries)
  all_logpost[1] <- aa[[1]][length(aa[[1]])]
  logpost[1] <- all_logpost[1]
  n <- 1
  num[1] <- 1
  iter <- 1
  jump_point <- numeric()
  jump_point[1] <- 1 # parameters for each chain
  if (verbose) {
    # initialize progress bar
    cat("Running BBNI MCMC Sampling...\n")
    pb <- txtProgressBar(min = 0, max = num_update, style = 3)
  }
  for (ii in 1:num_update) { # run ii full rounds, with each round of num.node times
    if (ii <= round(0.1 * num_update)) { # use adaptive ratio of using proposal information
      prop.ratio <- 0.1
    }
    if (ii > round(0.1 * num_update)) {
      prop.ratio <- 0.9
    }

    update_order <- sample.int(num.node, num.node, replace = FALSE)
    for (k in seq_along(update_order)) { #   consider the updating node g_k
      if (verbose) {
        # update progress bar
        setTxtProgressBar(pb, ii)
      }
      old <- n
      current_incid_matrix <- Incidence_Matrix[[iter]]
      current_ances_matrix <- Ancestor_Matrix[[iter]]
      current_trans_func_matrix <- Trans_Func_Matrix[[iter]]
      current_post <- all_logpost[iter]

      parent_of_update <- numeric()
      j <- 1 # find the parent for node update_order[k]
      for (i in seq_len(ncol(current_incid_matrix))) {
        if (current_incid_matrix[update_order[k], i] != 0) {
          parent_of_update[j] <- i
          j <- j + 1
        }
      }
      swap_candi <- numeric()
      legal_parent <- numeric()
      j1 <- 1 # legal_parent is used for adding parent action
      for (i in 1:num.node) {
        if (i != update_order[k] && current_ances_matrix[i, update_order[k]] != 1 && !(i %in% parent_of_update)) {
          legal_parent[j1] <- i
          j1 <- j1 + 1
        }
      }
      swap_candi <- legal_parent
      pairwise.prior.set <- matrix() # clear variables
      pairwise.prior.pare <- matrix()
      pairwise.prior.set <- prior.pairwise[prior.pairwise[, 2] == update_order[k], ] # pairwise proposal for node update_order[k]
      pairwise.prior.set <- data.matrix(pairwise.prior.set)
      if (ncol(pairwise.prior.set) == 1) {
        pairwise.prior.set <- t(pairwise.prior.set)
        pairwise.prior.pare <- t(data.matrix(pairwise.prior.set[, 1]))
      }
      if (ncol(pairwise.prior.set) > 1 && nrow(pairwise.prior.set) > 0) {
        pairwise.prior.pare <- pairwise.prior.set[, 1]
      }

      triplet.prior.set <- matrix()
      triplet.prior.pare <- matrix()
      triplet.prior.set <- prior.triplet[prior.triplet[, 3] == update_order[k], ] # triplet proposal for node update_order[k]
      triplet.prior.set <- data.matrix(triplet.prior.set)
      if (ncol(triplet.prior.set) > 1 && nrow(triplet.prior.set) > 0) {
        triplet.prior.pare <- triplet.prior.set[, 1:2]
      }
      if (ncol(triplet.prior.set) == 1) {
        triplet.prior.set <- t(triplet.prior.set)
        triplet.prior.pare <- t(data.matrix(triplet.prior.set[, 1:2]))
      }
      pairwise.ratio <- 1 / 2
      # case one: no parent
      if (length(parent_of_update) == 0) { #  there are only add-parent moves, every time only one move can be proposed
        num_legal_parent <- length(legal_parent)
        if (num_legal_parent > 1) {
          uu <- runif(1)
        }
        if (num_legal_parent == 0) { # if no candidate nodes are available for adding parent, then ignore this node
          uu <- 0
        }
        if (num_legal_parent == 1 || (num_legal_parent > 1 && uu >= pairwise.ratio)) {
          # proposal move 1: add one parent
          prop.legal.overlap <- intersect(pairwise.prior.pare, legal_parent)

          prop.prob <- runif(1)
          aa <- 0
          if (prop.prob >= prop.ratio && length(prop.legal.overlap) >= 1 && nrow(pairwise.prior.set) > 0) { # proposal information is used with probability 1-prop.ratio
            aa <- aa + 1
            candidate.pare.set <- pairwise.prior.set[pairwise.prior.set[, 1] %in% prop.legal.overlap, ]
            if (is.numeric(candidate.pare.set)) {
              add_parent <- candidate.pare.set
            }
            if (is.matrix(candidate.pare.set)) {
              add_parent <- candidate.pare.set[candidate.pare.set[, 4] == max(candidate.pare.set[, 4]), ]
            } # use the most likely one
            if (is.matrix(add_parent) && nrow(add_parent) > 1) {
              add_parent <- add_parent[sample.int(nrow(add_parent), 1), ]
            }

            add_one_parent <- add_parent[1]
          }
          if (prop.prob < prop.ratio || aa == 0) {
            sample_node <- sample.int(length(legal_parent), 1, replace = FALSE)
            add_one_parent <- legal_parent[sample_node]
          }
          prop_incid_matrix <- current_incid_matrix
          prop_incid_matrix[update_order[k], add_one_parent] <- 1
          prop_ances_matrix <- update_ancestor_matrix(prop_incid_matrix)
          prop_trans_func_matrix <- current_trans_func_matrix
          if (aa > 0) {
            prop_trans_func_matrix[update_order[k], add_one_parent] <- add_parent[3]
          } # function prior
          if (aa == 0) {
            func_order <- 10 + sample.int(2, 1)
            prop_trans_func_matrix[update_order[k], add_one_parent] <- func_order
          }
          add_one_prob <- 1 / length(legal_parent)
          prop_trans_func_prob <- 1 / 2 #  there are in all two types of boolean functions to choose for pairwise genes
          xxx <- Error_LLH(GeneData = GeneData, SampleSize = SampleSize, num.node = num.node, prior_para = prior_para, penalty = penalty, TRFUM = prop_trans_func_matrix, timeseries = timeseries)
          prop_sample <- xxx[[2]]
          prop_post <- xxx[[1]][length(xxx[[1]])]
          prop_sample_prob <- numeric()
          p2c_prob <- numeric() # calculate acceptance probability
          curr_sample_prob <- numeric()
          c2p_prob <- numeric()
          prop_sample_prob[1] <- prop_post
          prop_sample_prob[2] <- 0 #  this is the prior && likelihood
          curr_sample_prob[1] <- current_post[1]
          curr_sample_prob[2] <- 0
          prop_sample_prob[3] <- log(1 / 2)
          curr_sample_prob[3] <- log(1) #  this is P(R|T)

          if (num_legal_parent > 1) {
            p2c_prob[1] <- 1 / 3
          } # there are 3 possible moves, add one; remove one; swap one
          if (num_legal_parent == 1) {
            p2c_prob[1] <- 1
          }

          if (num_legal_parent > 1) {
            c2p_prob[1] <- 1 / 2 * add_one_prob
          } else {
            c2p_prob[1] <- add_one_prob
          }

          p2c_prob[2] <- 1
          c2p_prob[2] <- prop_trans_func_prob
          p2c_prob[3] <- 1
          c2p_prob[3] <- 1

          nume <- sum((prop_sample_prob)) + sum(log(p2c_prob))
          deno <- sum((curr_sample_prob)) + sum(log(c2p_prob))
          acce_prob <- exp(nume - deno)
          ratio <- runif(1)
          if (ratio <= acce_prob) {
            n <- n + 1
          }
        }
        # proposal move 2 add two parents
        if (num_legal_parent > 1 && uu < pairwise.ratio) { # add two parents one time
          prop.prob <- runif(1)
          aa <- 0
          if (prop.prob >= prop.ratio && nrow(triplet.prior.set) >= 1) {
            candidate.pare.set <- list()
            for (jj in seq_len(nrow(triplet.prior.set))) {
              if (length(intersect(triplet.prior.set[jj, 1:2], legal_parent)) == 2) {
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
          if (prop.prob < prop.ratio || aa == 0) {
            sample_two_parent <- sample.int(length(legal_parent), 2, replace = FALSE)
            add_two_parent <- c(legal_parent[sample_two_parent[1]], legal_parent[sample_two_parent[2]])
          }
          prop_incid_matrix <- current_incid_matrix
          prop_incid_matrix[update_order[k], add_two_parent[1]] <- 1
          prop_incid_matrix[update_order[k], add_two_parent[2]] <- 1
          prop_ances_matrix <- update_ancestor_matrix(prop_incid_matrix)
          prop_trans_func_matrix <- current_trans_func_matrix
          if (aa > 0) {
            prop_trans_func_matrix[update_order[k], add_two_parent[1]] <- add_parent[4]
            prop_trans_func_matrix[update_order[k], add_two_parent[2]] <- add_parent[4]
          }
          if (aa == 0) {
            func_order <- sample.int(10, 1)
            prop_trans_func_matrix[update_order[k], add_two_parent[1]] <- func_order
            prop_trans_func_matrix[update_order[k], add_two_parent[2]] <- func_order
          }
          sample_two_prob <- 1 / choose(num_legal_parent, 2)
          prop_trans_func_prob <- 1 / 10
          xxx <- Error_LLH(GeneData = GeneData, SampleSize = SampleSize, num.node = num.node, prior_para = prior_para, penalty = penalty, TRFUM = prop_trans_func_matrix, timeseries = timeseries)
          prop_sample <- xxx[[2]]
          prop_post <- xxx[[1]][length(xxx[[1]])]
          prop_sample_prob <- numeric()
          p2c_prob <- numeric() # calculate acceptance probability
          curr_sample_prob <- numeric()
          c2p_prob <- numeric()
          prop_sample_prob[1] <- prop_post
          prop_sample_prob[2] <- 0 #  this is the prior && likelihood
          curr_sample_prob[1] <- current_post[1]
          curr_sample_prob[2] <- 0
          prop_sample_prob[3] <- log(prop_trans_func_prob)
          curr_sample_prob[3] <- log(1) #  this is P(R|T)

          if (length(swap_candi) > 1) {
            p2c_prob[1] <- 1 / 4
          } #  this is Q(T_c|T_p) inverse move, because there are 4 moves after adding two parents, depending on the number of nodes to swap
          if (length(swap_candi) == 1) {
            p2c_prob[1] <- 1 / 3
          }
          if (length(swap_candi) == 0) {
            p2c_prob[1] <- 1 / 2
          }
          c2p_prob[1] <- 1 / 2 * sample_two_prob # this is Q(T_p|T_c) there are adding one parent, adding two parents, two moves
          p2c_prob[2] <- 1
          c2p_prob[2] <- prop_trans_func_prob
          p2c_prob[3] <- 1
          c2p_prob[3] <- 1

          nume <- sum((prop_sample_prob)) + sum(log(p2c_prob))
          deno <- sum((curr_sample_prob)) + sum(log(c2p_prob))
          acce_prob <- exp(nume - deno)
          ratio <- runif(1)
          if (ratio <= acce_prob) {
            n <- n + 1
          }
        }
      }
      # case 2: one parent
      if (length(parent_of_update) == 1) {
        uu <- runif(1)
        one.parent.ratio <- 1 / 6
        # proposal move 1: add one parent
        if (length(legal_parent) > 0 && uu < one.parent.ratio) {
          prop.prob <- runif(1)
          aa <- 0
          if (prop.prob >= prop.ratio && nrow(triplet.prior.set) >= 1 && ncol(triplet.prior.set) > 1) {
            candidate.pare.set <- list()
            for (jj in seq_len(nrow(triplet.prior.set))) {
              if (length(intersect(parent_of_update, triplet.prior.pare[jj, ])) == 1 && parent_of_update %in% triplet.prior.pare[jj, ] && setdiff(triplet.prior.pare[jj, ], parent_of_update) %in% legal_parent) {
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

              add_one_parent <- setdiff(add_parent[1:2], parent_of_update)
            }
          }
          if (prop.prob < prop.ratio || aa == 0) {
            sample_node <- sample.int(length(legal_parent), 1, replace = FALSE)
            add_one_parent <- legal_parent[sample_node]
          }
          prop_incid_matrix <- current_incid_matrix
          prop_incid_matrix[update_order[k], add_one_parent] <- 1
          prop_ances_matrix <- update_ancestor_matrix(prop_incid_matrix)
          prop_trans_func_matrix <- current_trans_func_matrix
          if (aa > 0) {
            prop_trans_func_matrix[update_order[k], add_one_parent] <- add_parent[4]
            prop_trans_func_matrix[update_order[k], parent_of_update] <- add_parent[4]
          }
          if (aa == 0) {
            func_order <- sample.int(10, 1)
            prop_trans_func_matrix[update_order[k], add_one_parent] <- func_order
            prop_trans_func_matrix[update_order[k], parent_of_update] <- func_order
          }
          add_one_prob <- 1 / length(legal_parent)
          prop_trans_func_prob <- 1 / 10 # there are 10 functions to choose
          xxx <- Error_LLH(GeneData = GeneData, SampleSize = SampleSize, num.node = num.node, prior_para = prior_para, penalty = penalty, TRFUM = prop_trans_func_matrix, timeseries = timeseries)
          prop_sample <- xxx[[2]]
          prop_post <- xxx[[1]][length(xxx[[1]])]
          prop_sample_prob <- numeric()
          p2c_prob <- numeric() # calculate acceptance probability
          curr_sample_prob <- numeric()
          c2p_prob <- numeric()
          prop_sample_prob[1] <- prop_post
          prop_sample_prob[2] <- 0 #  this is the prior && likelihood
          curr_sample_prob[1] <- current_post[1]
          curr_sample_prob[2] <- 0
          prop_sample_prob[3] <- log(prop_trans_func_prob)
          curr_sample_prob[3] <- log(1) #  this is P(R|T)

          if (length(swap_candi) > 1) {
            p2c_prob[1] <- 1 / 4
          } #  this is Q(T_c|T_p) inverse move, because there are 4 moves after adding one parents, depending on the number of nodes to swap
          if (length(swap_candi) == 1) {
            p2c_prob[1] <- 1 / 3
          }
          if (length(swap_candi) == 0) {
            p2c_prob[1] <- 1 / 2
          }

          if (length(swap_candi) > 0) {
            c2p_prob[1] <- 1 / 3 * add_one_prob
          } else {
            c2p_prob[1] <- 1 / 2 * add_one_prob
          }

          p2c_prob[2] <- 1
          c2p_prob[2] <- prop_trans_func_prob
          p2c_prob[3] <- 1
          c2p_prob[3] <- 1

          nume <- sum((prop_sample_prob)) + sum(log(p2c_prob))
          deno <- sum((curr_sample_prob)) + sum(log(c2p_prob))
          acce_prob <- exp(nume - deno)
          ratio <- runif(1)
          if (ratio <= acce_prob) {
            n <- n + 1
          }
        }
        # proposal move 2, swap one parent
        if (length(swap_candi) > 0 && uu >= 1 / 6 && uu < 2 / 6) {
          swap.prior.overlap <- intersect(swap_candi, pairwise.prior.pare)
          prop.prob <- runif(1)
          aa <- 0
          if (prop.prob > prop.ratio && length(swap.prior.overlap) > 0 && nrow(pairwise.prior.set) > 0) {
            aa <- aa + 1
            swap.candidate <- pairwise.prior.set[pairwise.prior.set[, 1] %in% swap.prior.overlap, ]
            if (is.numeric(swap.candidate)) { # i.e. swap.candidate has one row
              swap_parent <- swap.candidate
            }
            if (is.matrix(swap.candidate)) {
              swap_parent <- swap.candidate[swap.candidate[, 4] == max(swap.candidate[, 4]), ]
            } # multiple rows may have same weight

            if (is.matrix(swap_parent) && nrow(swap_parent) > 1) {
              swap_parent <- swap_parent[sample.int(nrow(swap_parent), 1), ]
            } # every candidate has equal chance
            swap_one_node <- swap_parent[1]
          }
          if (prop.prob < prop.ratio || aa == 0) {
            sample_node <- sample.int(length(swap_candi), 1, replace = FALSE)
            swap_one_node <- swap_candi[sample_node]
          }
          prop_incid_matrix <- current_incid_matrix
          prop_incid_matrix[update_order[k], parent_of_update] <- 0
          prop_incid_matrix[update_order[k], swap_one_node] <- 1
          prop_ances_matrix <- update_ancestor_matrix(prop_incid_matrix)
          prop_trans_func_matrix <- current_trans_func_matrix
          prop_trans_func_matrix[update_order[k], parent_of_update] <- 0
          if (aa > 0) {
            prop_trans_func_matrix[update_order[k], swap_one_node] <- swap_parent[3]
          }
          if (aa == 0) {
            func_order <- 10 + sample.int(2, 1)
            prop_trans_func_matrix[update_order[k], swap_one_node] <- func_order
          }
          swap_one_prob <- 1 / length(swap_candi)
          prop_trans_func_prob <- 1 / 2 #  new proposal transition function matrix
          xxx <- Error_LLH(GeneData = GeneData, SampleSize = SampleSize, num.node = num.node, prior_para = prior_para, penalty = penalty, TRFUM = prop_trans_func_matrix, timeseries = timeseries)
          prop_sample <- xxx[[2]]
          prop_post <- xxx[[1]][length(xxx[[1]])]
          prop_sample_prob <- numeric()
          p2c_prob <- numeric() # calculate acceptance probability
          curr_sample_prob <- numeric()
          c2p_prob <- numeric()
          prop_sample_prob[1] <- prop_post
          prop_sample_prob[2] <- 0 #  this is the prior && likelihood
          curr_sample_prob[1] <- current_post[1]
          curr_sample_prob[2] <- 0
          prop_sample_prob[3] <- log(prop_trans_func_prob)
          curr_sample_prob[3] <- log(1) #  this is P(R|T)

          if (length(legal_parent) > 0) {
            p2c_prob[1] <- 1 / 3
          } #  this is Q(T_c|T_p) inverse move, because there are 4 moves after adding one parents, depending on the number of nodes to swap
          if (length(legal_parent) == 0) {
            p2c_prob[1] <- 1 / 2
          }

          if (length(legal_parent) > 0) {
            c2p_prob[1] <- 1 / 3 * swap_one_prob
          } # this is Q(T_p|T_c) probability of which move and which node selected
          if (length(legal_parent) == 0) {
            c2p_prob[1] <- 1 / 2 * swap_one_prob
          }
          p2c_prob[2] <- 1
          c2p_prob[2] <- prop_trans_func_prob
          p2c_prob[3] <- 1
          c2p_prob[3] <- 1

          nume <- sum((prop_sample_prob)) + sum(log(p2c_prob))
          deno <- sum((curr_sample_prob)) + sum(log(c2p_prob))
          acce_prob <- exp(nume - deno)
          ratio <- runif(1)
          if (ratio <= acce_prob) {
            n <- n + 1
          }
        }
        if (uu >= 2 / 6 && uu < 3 / 6) {
          remove_one_node <- parent_of_update
          remove_one_prob <- 1
          prop_incid_matrix <- current_incid_matrix
          prop_incid_matrix[update_order[k], remove_one_node] <- 0 # remove parent move does not need function.
          prop_ances_matrix <- update_ancestor_matrix(prop_incid_matrix)
          prop_trans_func_matrix <- current_trans_func_matrix
          prop_trans_func_matrix[update_order[k], remove_one_node] <- 0
          prop_trans_func_prob <- 1
          xxx <- Error_LLH(GeneData = GeneData, SampleSize = SampleSize, num.node = num.node, prior_para = prior_para, penalty = penalty, TRFUM = prop_trans_func_matrix, timeseries = timeseries)
          prop_sample <- xxx[[2]]
          prop_post <- xxx[[1]][length(xxx[[1]])]
          prop_sample_prob <- numeric()
          p2c_prob <- numeric() # calculate acceptance probability
          curr_sample_prob <- numeric()
          c2p_prob <- numeric()
          prop_sample_prob[1] <- prop_post
          prop_sample_prob[2] <- 0 #  this is the prior && likelihood
          curr_sample_prob[1] <- current_post[1]
          curr_sample_prob[2] <- 0
          prop_sample_prob[3] <- log(prop_trans_func_prob)
          curr_sample_prob[3] <- log(1) #  this is P(R|T)

          if (length(legal_parent) == 0) {
            p2c_prob[1] <- 1
          }
          if (length(legal_parent) == 1) {
            p2c_prob[1] <- 1
          } #  this is Q(T_c|T_p) inverse move
          if (length(legal_parent) > 1) {
            p2c_prob[1] <- 1 / 2
          }

          if (length(legal_parent) > 0) {
            c2p_prob[1] <- 1 / 3 * remove_one_prob
          } # this is Q(T_p|T_c) probability of which move and which node selected
          if (length(legal_parent) == 0) {
            c2p_prob[1] <- 1 / 2 * remove_one_prob
          }
          p2c_prob[2] <- 1 / 2
          c2p_prob[2] <- prop_trans_func_prob
          p2c_prob[3] <- 1
          c2p_prob[3] <- 1

          nume <- sum((prop_sample_prob)) + sum(log(p2c_prob))
          deno <- sum((curr_sample_prob)) + sum(log(c2p_prob))
          acce_prob <- exp(nume - deno)
          ratio <- runif(1)
          if (ratio <= acce_prob) {
            n <- n + 1
          }
        }
        # additional move 1: reverse one arc in pairwise genes
        if (uu >= 3 / 6 && uu < 4 / 6) {
          prop_incid_matrix <- current_incid_matrix
          prop_incid_matrix[update_order[k], parent_of_update] <- 0
          prop_incid_matrix[parent_of_update, update_order[k]] <- 1
          parent_parent <- numeric()
          j <- 0
          for (i in seq_len(nrow(current_incid_matrix))) {
            if (current_incid_matrix[parent_of_update, i] > 0) {
              j <- j + 1
              parent_parent[j] <- i
            }
          }
          prop_ances_matrix <- update_ancestor_matrix(prop_incid_matrix)
          if (check_ances_matrix(prop_ances_matrix) == 0 && length(parent_parent) == 0) { # make sure no loops && parent_of_update no parents
            prop_trans_func_matrix <- current_trans_func_matrix
            func_order <- current_trans_func_matrix[update_order[k], parent_of_update]
            prop_trans_func_matrix[update_order[k], parent_of_update] <- 0
            prop_trans_func_matrix[parent_of_update, update_order[k]] <- func_order
            xxx <- Error_LLH(GeneData = GeneData, SampleSize = SampleSize, num.node = num.node, prior_para = prior_para, penalty = penalty, TRFUM = prop_trans_func_matrix, timeseries = timeseries)
            prop_sample <- xxx[[2]]
            prop_post <- xxx[[1]][length(xxx[[1]])]

            nume <- sum((prop_post))
            deno <- sum((current_post))
            acce_prob <- exp(nume - deno)
            ratio <- runif(1)
            if (ratio <= acce_prob) {
              n <- n + 1
            }
          }
        }
        # additional move 2: reverse two arcs simultanously
        if (uu >= 4 / 6 && uu < 5 / 6) {
          parent_parent <- numeric()
          j1 <- 0
          children <- numeric()
          j <- 0 # find the children of parent_of_update
          for (i in seq_len(ncol(current_incid_matrix))) {
            if (i != update_order[k] && current_incid_matrix[i, parent_of_update] > 0) {
              j <- j + 1
              children[j] <- i
            }
            if (current_incid_matrix[parent_of_update, i] > 0) {
              j1 <- j1 + 1
              parent_parent[j1] <- i
            }
          }
          if (length(children) == 1 && length(parent_parent) == 0) { # parent_of_update has two exact 2 children and no parents
            parent_child <- numeric()
            j <- 0 # find parent of children to determine the trans_func_matrix
            for (i in 1:num.node) {
              if (current_incid_matrix[children, i] > 0) {
                j <- j + 1
                parent_child[j] <- i
              }
            }

            prop_incid_matrix <- current_incid_matrix
            prop_incid_matrix[update_order[k], parent_of_update] <- 0
            prop_incid_matrix[children, parent_of_update] <- 0
            prop_incid_matrix[parent_of_update, update_order[k]] <- 1
            prop_incid_matrix[parent_of_update, children] <- 1
            prop_ances_matrix <- update_ancestor_matrix(prop_incid_matrix)
            if (check_ances_matrix(prop_ances_matrix) == 0) { # ensure no directed cycles
              prop_trans_func_matrix <- current_trans_func_matrix
              prop_trans_func_matrix[update_order[k], parent_of_update] <- 0
              prop_trans_func_matrix[children, parent_of_update] <- 0
              func_order <- sample.int(10, 1)
              prop_trans_func_matrix[parent_of_update, update_order[k]] <- func_order
              prop_trans_func_matrix[parent_of_update, children] <- func_order
              if (length(parent_child) == 2) {
                prop_trans_func_matrix[children, setdiff(parent_child, parent_of_update)] <- 10 + sample.int(2, 1)
              }
              xxx <- Error_LLH(GeneData = GeneData, SampleSize = SampleSize, num.node = num.node, prior_para = prior_para, penalty = penalty, TRFUM = prop_trans_func_matrix, timeseries = timeseries)
              prop_sample <- xxx[[2]]
              prop_post <- xxx[[1]][length(xxx[[1]])]

              nume <- sum((prop_post))
              deno <- sum((current_post))
              acce_prob <- exp(nume - deno)
              ratio <- runif(1)
              if (ratio <= acce_prob) {
                n <- n + 1
              }
            }
          }
        }
        if (uu >= 5 / 6 && length(swap_candi) > 1) {
          legal.triplet.pare <- intersect(swap_candi, triplet.prior.pare)
          prop.prob <- runif(1)
          aa <- 0
          if (prop.prob > prop.ratio) {
            swap.candidate <- list()
            if (length(intersect(parent_of_update, legal.triplet.pare)) == 0 && nrow(triplet.prior.set) > 0) {
              for (jj in seq_len(nrow(triplet.prior.pare))) {
                if (length(intersect(parent_of_update, triplet.prior.pare[jj, ])) == 0 && length(intersect(triplet.prior.pare[jj, ], swap_candi)) == 2) {
                  aa <- aa + 1
                  swap.candidate[[aa]] <- triplet.prior.set[jj, ]
                }
              }
            }
            if (length(swap.candidate) > 0) {
              score <- numeric()
              for (jj in seq_along(swap.candidate)) {
                score[jj] <- swap.candidate[[jj]][5]
              }
              bbb <- matrix(nrow = length(swap.candidate), ncol = 6) # convert list to matrix
              for (jjj in seq_along(swap.candidate)) {
                bbb[jjj, ] <- swap.candidate[[jjj]]
              }
              max.score.candidate <- bbb[bbb[, 5] == max(score), ]
              if (is.numeric(max.score.candidate)) {
                swap_parent <- max.score.candidate
              }
              if (is.matrix(max.score.candidate) && nrow(max.score.candidate) > 1) {
                swap_parent <- max.score.candidate[sample.int(nrow(max.score.candidate), 1), ]
              }
              swap_two_parent <- swap_parent[1:2]
            }
          }
          if (prop.prob < prop.ratio || aa == 0) {
            sample_two_node <- sample.int(length(swap_candi), 2, replace = FALSE)
            swap_two_parent <- c(swap_candi[sample_two_node[1]], swap_candi[sample_two_node[2]])
          }
          prop_incid_matrix <- current_incid_matrix
          prop_incid_matrix[update_order[k], parent_of_update] <- 0
          prop_incid_matrix[update_order[k], swap_two_parent[1]] <- 1
          prop_incid_matrix[update_order[k], swap_two_parent[2]] <- 1
          prop_ances_matrix <- update_ancestor_matrix(prop_incid_matrix)
          prop_trans_func_matrix <- current_trans_func_matrix
          prop_trans_func_matrix[update_order[k], parent_of_update] <- 0
          if (aa > 0) {
            prop_trans_func_matrix[update_order[k], swap_two_parent[1]] <- swap_parent[4]
            prop_trans_func_matrix[update_order[k], swap_two_parent[2]] <- swap_parent[4]
          }
          if (aa == 0) {
            func_order <- sample.int(10, 1)
            prop_trans_func_matrix[update_order[k], swap_two_parent[1]] <- func_order
            prop_trans_func_matrix[update_order[k], swap_two_parent[2]] <- func_order
          }
          swap_two_prob <- 1 / choose(length(swap_candi), 2)
          prop_trans_func_prob <- 1 / 10
          xxx <- Error_LLH(GeneData = GeneData, SampleSize = SampleSize, num.node = num.node, prior_para = prior_para, penalty = penalty, TRFUM = prop_trans_func_matrix, timeseries = timeseries)
          prop_sample <- xxx[[2]]
          prop_post <- xxx[[1]][length(xxx[[1]])]
          prop_sample_prob <- numeric()
          p2c_prob <- numeric() # calculate acceptance probability
          curr_sample_prob <- numeric()
          c2p_prob <- numeric()

          prop_sample_prob[1] <- prop_post
          prop_sample_prob[2] <- 0 #  this is the prior && likelihood
          curr_sample_prob[1] <- current_post[1]
          curr_sample_prob[2] <- 0
          prop_sample_prob[3] <- log(prop_trans_func_prob)
          curr_sample_prob[3] <- log(1) #  this is P(R|T)

          if (length(swap_candi) > 1) {
            p2c_prob[1] <- 1 / 4
          } #  this is Q(T_c|T_p) inverse move, because there are 4 moves after adding one parents, depending on the number of nodes to swap
          if (length(swap_candi) == 1) {
            p2c_prob[1] <- 1 / 3
          }

          if (length(swap_candi) > 1) {
            c2p_prob[1] <- 1 / 4 * swap_two_prob
          } # this is Q(T_p|T_c) probability of which move and which node selected
          if (length(swap_candi) == 1) {
            c2p_prob[1] <- 1 / 3 * swap_two_prob
          }
          p2c_prob[2] <- 1 / 10
          c2p_prob[2] <- prop_trans_func_prob
          p2c_prob[3] <- 1
          c2p_prob[3] <- 1

          nume <- sum((prop_sample_prob)) + sum(log(p2c_prob))
          deno <- sum((curr_sample_prob)) + sum(log(c2p_prob))
          acce_prob <- exp(nume - deno)
          ratio <- runif(1)
          if (ratio <= acce_prob) {
            n <- n + 1
          }
        }
      }
      # case3: two parents
      if (length(parent_of_update) == 2) {
        uu <- runif(1)
        two.parent.ratio <- 1 / 7
        # proposal move 1: remove one parent
        if (uu < two.parent.ratio) {
          pare.prior.overlap <- intersect(parent_of_update, pairwise.prior.pare)
          prop.prob <- runif(1)
          aa <- 0
          if (prop.prob >= prop.ratio && length(pare.prior.overlap) == 1 && nrow(pairwise.prior.set) > 0) {
            aa <- aa + 1
            remove_one_node <- setdiff(parent_of_update, pairwise.prior.pare)
            remove_pare <- pairwise.prior.set[pairwise.prior.set[, 1] == setdiff(parent_of_update, remove_one_node), ]
            if (is.numeric(remove_pare)) {
              func_order <- remove_pare[3]
            }
            if (is.matrix(remove_pare)) {
              Remove_Pare <- remove_pare[remove_pare[, 4] == max(remove_pare[, 4]), ]
              if (is.matrix(Remove_Pare) && nrow(Remove_Pare) > 1) {
                Remove_Pare <- Remove_Pare[sample.int(nrow(Remove_Pare), 1), ]
              } # each candidate has equal chance
              func_order <- Remove_Pare[3]
            }
          }
          if (prop.prob < prop.ratio || aa == 0) {
            sample_node <- sample.int(length(parent_of_update), 1, replace = FALSE)
            remove_one_node <- parent_of_update[sample_node]
          }
          remove_one_prob <- 1 / length(parent_of_update)
          prop_incid_matrix <- current_incid_matrix
          prop_incid_matrix[update_order[k], remove_one_node] <- 0
          prop_ances_matrix <- update_ancestor_matrix(prop_incid_matrix)
          prop_trans_func_matrix <- current_trans_func_matrix
          prop_trans_func_matrix[update_order[k], remove_one_node] <- 0
          if (aa > 0) {
            prop_trans_func_matrix[update_order[k], setdiff(parent_of_update, remove_one_node)] <- func_order
          } # by the definition
          if (aa == 0) {
            func_order <- 10 + sample.int(2, 1)
          }
          prop_trans_func_matrix[update_order[k], setdiff(parent_of_update, remove_one_node)] <- func_order
          prop_trans_func_prob <- 1 / 2
          xxx <- Error_LLH(GeneData = GeneData, SampleSize = SampleSize, num.node = num.node, prior_para = prior_para, penalty = penalty, TRFUM = prop_trans_func_matrix, timeseries = timeseries)
          prop_sample <- xxx[[2]]
          prop_post <- xxx[[1]][length(xxx[[1]])]
          prop_sample_prob <- numeric()
          p2c_prob <- numeric() # calculate acceptance probability
          curr_sample_prob <- numeric()
          c2p_prob <- numeric()

          prop_sample_prob[1] <- prop_post
          prop_sample_prob[2] <- 0 #  this is the prior && likelihood
          curr_sample_prob[1] <- current_post[1]
          curr_sample_prob[2] <- 0
          prop_sample_prob[3] <- log(prop_trans_func_prob)
          curr_sample_prob[3] <- log(1) #  this is P(R|T)

          if (length(swap_candi) > 0) {
            p2c_prob[1] <- 1 / 3
          } #  this is Q(T_c|T_p) inverse move, because there are 4 moves after adding one parents, depending on the number of nodes to swap
          if (length(swap_candi) == 0) {
            p2c_prob[1] <- 1 / 2
          }

          if (length(swap_candi) == 0) {
            c2p_prob[1] <- 1 / 2
          }
          if (length(swap_candi) > 1) {
            c2p_prob[1] <- 1 / 4 * remove_one_prob
          } # this is Q(T_p|T_c) probability of which move and which node selected
          if (length(swap_candi) == 1) {
            c2p_prob[1] <- 1 / 3 * remove_one_prob
          }

          p2c_prob[2] <- 1 / 10
          c2p_prob[2] <- prop_trans_func_prob
          p2c_prob[3] <- 1
          c2p_prob[3] <- 1

          nume <- sum((prop_sample_prob)) + sum(log(p2c_prob))
          deno <- sum((curr_sample_prob)) + sum(log(c2p_prob))
          acce_prob <- exp(nume - deno)
          ratio <- runif(1)
          if (ratio <= acce_prob) {
            n <- n + 1
          }
        }
        # proposal move 2; remove two parents
        if (uu >= two.parent.ratio && uu < 2 / 7) {
          remove_node <- parent_of_update
          remove_prob <- 1
          prop_incid_matrix <- current_incid_matrix
          prop_incid_matrix[update_order[k], remove_node[1]] <- 0
          prop_incid_matrix[update_order[k], remove_node[2]] <- 0
          prop_ances_matrix <- update_ancestor_matrix(prop_incid_matrix)
          prop_trans_func_matrix <- current_trans_func_matrix
          prop_trans_func_matrix[update_order[k], remove_node[1]] <- 0
          prop_trans_func_matrix[update_order[k], remove_node[2]] <- 0
          prop_trans_func_prob <- 1
          xxx <- Error_LLH(GeneData = GeneData, SampleSize = SampleSize, num.node = num.node, prior_para = prior_para, penalty = penalty, TRFUM = prop_trans_func_matrix, timeseries = timeseries)
          prop_sample <- xxx[[2]]
          prop_post <- xxx[[1]][length(xxx[[1]])]
          prop_sample_prob <- numeric()
          p2c_prob <- numeric() # calculate acceptance probability
          curr_sample_prob <- numeric()
          c2p_prob <- numeric()
          prop_sample_prob[1] <- prop_post
          prop_sample_prob[2] <- 0 #  this is the prior && likelihood
          curr_sample_prob[1] <- current_post[1]
          curr_sample_prob[2] <- 0
          prop_sample_prob[3] <- log(prop_trans_func_prob)
          curr_sample_prob[3] <- log(1) #  this is P(R|T)

          p2c_prob[1] <- 1 / 2 #  this is Q(T_c|T_p) inverse move
          c2p_prob[1] <- 1 / 4 * remove_prob # this is Q(T_p|T_c) probability of which move and which node selected
          p2c_prob[2] <- 1 / 10
          c2p_prob[2] <- prop_trans_func_prob
          p2c_prob[3] <- 1
          c2p_prob[3] <- 1

          nume <- sum((prop_sample_prob)) + sum(log(p2c_prob))
          deno <- sum((curr_sample_prob)) + sum(log(c2p_prob))
          acce_prob <- exp(nume - deno)
          ratio <- runif(1)
          if (ratio <= acce_prob) {
            n <- n + 1
          }
        }
        # proposal move 3: swap one parent
        if (uu >= 2 / 7 && uu < 3 / 7 && length(swap_candi) > 0) {
          legal.triplet.pare <- intersect(swap_candi, triplet.prior.pare)
          prop.prob <- runif(1)
          aa <- 0
          if (prop.prob >= prop.ratio) {
            swap.candidate <- list()
            if (length(intersect(parent_of_update, legal.triplet.pare)) == 1 && nrow(triplet.prior.set) > 0) {
              for (jj in seq_len(nrow(triplet.prior.pare))) {
                if (length(intersect(parent_of_update, triplet.prior.pare[jj, ])) == 1) {
                  aa <- aa + 1
                  swap.candidate[[aa]] <- triplet.prior.set[jj, ]
                }
              }
            }
            if (length(swap.candidate) > 0) {
              score <- numeric()
              for (jj in seq_along(swap.candidate)) {
                score[jj] <- swap.candidate[[jj]][5]
              }
              bbb <- matrix(nrow = length(swap.candidate), ncol = 6)
              for (jjj in seq_along(swap.candidate)) {
                bbb[jjj, ] <- swap.candidate[[jjj]]
              }

              max.score.candidate <- bbb[bbb[, 5] == max(score), ]
              if (is.numeric(max.score.candidate)) {
                swap_parent <- max.score.candidate
              }
              if (is.matrix(max.score.candidate) && nrow(max.score.candidate) > 1) {
                swap_parent <- max.score.candidate[sample.int(nrow(max.score.candidate), 1), ]
              }
              swap_one_node <- setdiff(swap_parent[1:2], parent_of_update)
              sample_one_parent <- setdiff(parent_of_update, swap_parent[1:2])
            }
          }
          if (prop.prob < prop.ratio || aa == 0) {
            sample_node <- sample.int(length(swap_candi), 1, replace = FALSE)
            sample_parent <- sample.int(2, 1, replace = FALSE)
            swap_one_node <- swap_candi[sample_node]
            sample_one_parent <- parent_of_update[sample_parent]
          }
          prop_incid_matrix <- current_incid_matrix
          prop_incid_matrix[update_order[k], sample_one_parent] <- 0
          prop_incid_matrix[update_order[k], swap_one_node] <- 1
          prop_ances_matrix <- update_ancestor_matrix(prop_incid_matrix)
          prop_trans_func_matrix <- current_trans_func_matrix
          prop_trans_func_matrix[update_order[k], sample_one_parent] <- 0
          if (aa > 0) {
            prop_trans_func_matrix[update_order[k], setdiff(parent_of_update, sample_one_parent)] <- swap_parent[4]
            prop_trans_func_matrix[update_order[k], swap_one_node] <- swap_parent[4]
          }
          if (aa == 0) {
            func_order <- sample.int(10, 1)
            prop_trans_func_matrix[update_order[k], setdiff(parent_of_update, sample_one_parent)] <- func_order
            prop_trans_func_matrix[update_order[k], swap_one_node] <- func_order
          }
          swap_one_prob <- 1 / length(swap_candi) * 1 / 2
          prop_trans_func_prob <- 1 / 10
          xxx <- Error_LLH(GeneData = GeneData, SampleSize = SampleSize, num.node = num.node, prior_para = prior_para, penalty = penalty, TRFUM = prop_trans_func_matrix, timeseries = timeseries)
          prop_sample <- xxx[[2]]
          prop_post <- xxx[[1]][length(xxx[[1]])]
          prop_sample_prob <- numeric()
          p2c_prob <- numeric() # calculate acceptance probability
          curr_sample_prob <- numeric()
          c2p_prob <- numeric()
          prop_sample_prob[1] <- prop_post
          prop_sample_prob[2] <- 0 #  this is the prior && likelihood
          curr_sample_prob[1] <- current_post[1]
          curr_sample_prob[2] <- 0
          prop_sample_prob[3] <- log(prop_trans_func_prob)
          curr_sample_prob[3] <- log(1) #  this is P(R|T)

          if (length(swap_candi) == 0) {
            p2c_prob[1] <- 1
          } # since we use log(p2c_prob)
          if (length(swap_candi) > 1) {
            p2c_prob[1] <- 1 / 4
          } #  this is Q(T_c|T_p) inverse move, because there are 4 moves after adding one parents, depending on the number of nodes to swap
          if (length(swap_candi) == 1) {
            p2c_prob[1] <- 1 / 3
          }

          if (length(swap_candi) == 0) {
            c2p_prob[1] <- 1
          }
          if (length(swap_candi) > 1) {
            c2p_prob[1] <- 1 / 4 * swap_one_prob
          } # this is Q(T_p|T_c) probability of which move and which node selected
          if (length(swap_candi) == 1) {
            c2p_prob[1] <- 1 / 3 * swap_one_prob
          }
          p2c_prob[2] <- 1 / 10
          c2p_prob[2] <- prop_trans_func_prob
          p2c_prob[3] <- 1
          c2p_prob[3] <- 1

          nume <- sum((prop_sample_prob)) + sum(log(p2c_prob))
          deno <- sum((curr_sample_prob)) + sum(log(c2p_prob))
          acce_prob <- exp(nume - deno)
          ratio <- runif(1)
          if (ratio <= acce_prob) {
            n <- n + 1
          }
        }
        # proposal move 4, swap two parents
        if (uu >= 3 / 7 && uu < 4 / 7 && length(swap_candi) > 1) {
          legal.triplet.pare <- intersect(swap_candi, triplet.prior.pare)
          prop.prob <- runif(1)
          aa <- 0
          if (prop.prob > prop.ratio) {
            swap.candidate <- list()
            if (length(intersect(parent_of_update, legal.triplet.pare)) == 0 && nrow(triplet.prior.set) > 0) {
              for (jj in seq_len(nrow(triplet.prior.pare))) {
                if (length(intersect(parent_of_update, triplet.prior.pare[jj, ])) == 0 && length(intersect(triplet.prior.pare[jj, ], swap_candi)) == 2) {
                  aa <- aa + 1
                  swap.candidate[[aa]] <- triplet.prior.set[jj, ]
                }
              }
            }
            if (length(swap.candidate) > 0) {
              score <- numeric()
              for (jj in seq_along(swap.candidate)) {
                score[jj] <- swap.candidate[[jj]][5]
              }
              bbb <- matrix(nrow = length(swap.candidate), ncol = 6) # convert list to matrix
              for (jjj in seq_along(swap.candidate)) {
                bbb[jjj, ] <- swap.candidate[[jjj]]
              }
              max.score.candidate <- bbb[bbb[, 5] == max(score), ]
              if (is.numeric(max.score.candidate)) {
                swap_parent <- max.score.candidate
              }
              if (is.matrix(max.score.candidate) && nrow(max.score.candidate) > 1) {
                swap_parent <- max.score.candidate[sample.int(nrow(max.score.candidate), 1), ]
              }
              swap_two_parent <- swap_parent[1:2]
            }
          }
          if (prop.prob < prop.ratio || aa == 0) {
            sample_two_node <- sample.int(length(swap_candi), 2, replace = FALSE)
            swap_two_parent <- c(swap_candi[sample_two_node[1]], swap_candi[sample_two_node[2]])
          }
          prop_incid_matrix <- current_incid_matrix
          prop_incid_matrix[update_order[k], parent_of_update[1]] <- 0
          prop_incid_matrix[update_order[k], parent_of_update[2]] <- 0
          prop_incid_matrix[update_order[k], swap_two_parent[1]] <- 1
          prop_incid_matrix[update_order[k], swap_two_parent[2]] <- 1
          prop_ances_matrix <- update_ancestor_matrix(prop_incid_matrix)
          prop_trans_func_matrix <- current_trans_func_matrix
          prop_trans_func_matrix[update_order[k], parent_of_update[1]] <- 0
          prop_trans_func_matrix[update_order[k], parent_of_update[2]] <- 0
          if (aa > 0) {
            prop_trans_func_matrix[update_order[k], swap_two_parent[1]] <- swap_parent[4]
            prop_trans_func_matrix[update_order[k], swap_two_parent[2]] <- swap_parent[4]
          }
          if (aa == 0) {
            func_order <- sample.int(10, 1)
            prop_trans_func_matrix[update_order[k], swap_two_parent[1]] <- func_order
            prop_trans_func_matrix[update_order[k], swap_two_parent[2]] <- func_order
          }
          swap_two_prob <- 1 / choose(length(swap_candi), 2)
          prop_trans_func_prob <- 1 / 10
          xxx <- Error_LLH(GeneData = GeneData, SampleSize = SampleSize, num.node = num.node, prior_para = prior_para, penalty = penalty, TRFUM = prop_trans_func_matrix, timeseries = timeseries)
          prop_sample <- xxx[[2]]
          prop_post <- xxx[[1]][length(xxx[[1]])]
          prop_sample_prob <- numeric()
          p2c_prob <- numeric() # calculate acceptance probability
          curr_sample_prob <- numeric()
          c2p_prob <- numeric()

          prop_sample_prob[1] <- prop_post
          prop_sample_prob[2] <- 0 #  this is the prior && likelihood
          curr_sample_prob[1] <- current_post[1]
          curr_sample_prob[2] <- 0
          prop_sample_prob[3] <- log(prop_trans_func_prob)
          curr_sample_prob[3] <- log(1) #  this is P(R|T)

          if (length(swap_candi) > 1) {
            p2c_prob[1] <- 1 / 4
          } #  this is Q(T_c|T_p) inverse move, because there are 4 moves after adding one parents, depending on the number of nodes to swap
          if (length(swap_candi) == 1) {
            p2c_prob[1] <- 1 / 3
          }

          if (length(swap_candi) > 1) {
            c2p_prob[1] <- 1 / 4 * swap_two_prob
          } # this is Q(T_p|T_c) probability of which move and which node selected
          if (length(swap_candi) == 1) {
            c2p_prob[1] <- 1 / 3 * swap_two_prob
          }
          p2c_prob[2] <- 1 / 10
          c2p_prob[2] <- prop_trans_func_prob
          p2c_prob[3] <- 1
          c2p_prob[3] <- 1

          nume <- sum((prop_sample_prob)) + sum(log(p2c_prob))
          deno <- sum((curr_sample_prob)) + sum(log(c2p_prob))
          acce_prob <- exp(nume - deno)
          ratio <- runif(1)
          if (ratio <= acce_prob) {
            n <- n + 1
          }
        }
        #  additional move 1: reorder the input and output variables  for XOR relation
        if (uu >= 4 / 7 && uu < 5 / 7 && current_trans_func_matrix[update_order[k], parent_of_update[1]] == 9) {
          prop_incid_matrix <- current_incid_matrix
          sample.child <- parent_of_update[sample.int(2, 1)]
          for (i in 1:2) {
            prop_incid_matrix[update_order[k], parent_of_update[i]] <- 0
          }
          parent_parent <- numeric()
          mm <- 0
          for (i in 1:num.node) {
            if (current_incid_matrix[sample.child, i] > 0) {
              mm <- mm + 1
              parent_parent[mm] <- i
            }
          }
          if (length(parent_parent) == 0) { # selected child should have no parents before this operation
            prop_incid_matrix[sample.child, update_order[k]] <- 1
            prop_incid_matrix[sample.child, setdiff(parent_of_update, sample.child)] <- 1
            prop_ances_matrix <- update_ancestor_matrix(prop_incid_matrix)
            if (check_ances_matrix(prop_ances_matrix) == 0) { # make sure no directed cycles
              prop_trans_func_matrix <- current_trans_func_matrix
              func_order <- current_trans_func_matrix[update_order[k], parent_of_update[1]]
              for (i in 1:2) {
                prop_trans_func_matrix[update_order[k], parent_of_update[i]] <- 0
              }
              prop_trans_func_matrix[sample.child, update_order[k]] <- func_order
              prop_trans_func_matrix[sample.child, setdiff(parent_of_update, sample.child)] <- func_order
              xxx <- Error_LLH(GeneData = GeneData, SampleSize = SampleSize, num.node = num.node, prior_para = prior_para, penalty = penalty, TRFUM = prop_trans_func_matrix, timeseries = timeseries)
              prop_sample <- xxx[[2]]
              prop_post <- xxx[[1]][length(xxx[[1]])]

              nume <- sum((prop_post))
              deno <- sum((current_post))
              acce_prob <- exp(nume - deno)
              ratio <- runif(1)
              if (ratio <= acce_prob) {
                n <- n + 1
              }
            }
          }
        }
        # additional move 2: reverse one arc among gene triplet
        if (uu >= 5 / 7 && uu < 6 / 7) {
          reverse_pare <- parent_of_update[sample.int(2, 1)]
          remain_pare <- setdiff(parent_of_update, reverse_pare)
          reverse_node_pare <- numeric()
          j <- 1 # find the already existing parent for reverse_pare
          for (i in seq_len(ncol(current_incid_matrix))) {
            if (current_incid_matrix[reverse_pare, i] != 0) {
              reverse_node_pare[j] <- i
              j <- j + 1
            }
          }
          if (length(reverse_node_pare) == 0) { # currently no parent
            prop_incid_matrix <- current_incid_matrix
            prop_incid_matrix[update_order[k], reverse_pare] <- 0
            prop_incid_matrix[reverse_pare, update_order[k]] <- 1
            prop_ances_matrix <- update_ancestor_matrix(prop_incid_matrix)
            if (check_ances_matrix(prop_ances_matrix) == 0) { # make sure no loops
              prop_trans_func_matrix <- current_trans_func_matrix
              prop_trans_func_matrix[update_order[k], reverse_pare] <- 0
              prop_trans_func_matrix[update_order[k], remain_pare] <- 10 + sample.int(2, 1)
              prop_trans_func_matrix[reverse_pare, update_order[k]] <- 10 + sample.int(2, 1)
              xxx <- Error_LLH(GeneData = GeneData, SampleSize = SampleSize, num.node = num.node, prior_para = prior_para, penalty = penalty, TRFUM = prop_trans_func_matrix, timeseries = timeseries)
              prop_sample <- xxx[[2]]
              prop_post <- xxx[[1]][length(xxx[[1]])]

              nume <- sum((prop_post))
              deno <- sum((current_post))
              acce_prob <- exp(nume - deno)
              ratio <- runif(1)
              if (ratio <= acce_prob) {
                n <- n + 1
              }
            }
          }
          if (length(reverse_node_pare) == 1) { # currently one parent
            prop_incid_matrix <- current_incid_matrix
            prop_incid_matrix[update_order[k], reverse_pare] <- 0
            prop_incid_matrix[reverse_pare, update_order[k]] <- 1
            prop_ances_matrix <- update_ancestor_matrix(prop_incid_matrix)
            if (check_ances_matrix(prop_ances_matrix) == 0) { # make sure no loops
              prop_trans_func_matrix <- current_trans_func_matrix
              prop_trans_func_matrix[update_order[k], reverse_pare] <- 0
              prop_trans_func_matrix[update_order[k], remain_pare] <- 10 + sample.int(2, 1)
              func_order <- sample.int(10, 1)
              prop_trans_func_matrix[reverse_pare, update_order[k]] <- func_order
              prop_trans_func_matrix[reverse_pare, reverse_node_pare] <- func_order
              xxx <- Error_LLH(GeneData = GeneData, SampleSize = SampleSize, num.node = num.node, prior_para = prior_para, penalty = penalty, TRFUM = prop_trans_func_matrix, timeseries = timeseries)
              prop_sample <- xxx[[2]]
              prop_post <- xxx[[1]][length(xxx[[1]])]

              nume <- sum((prop_post))
              deno <- sum((current_post))
              acce_prob <- exp(nume - deno)
              ratio <- runif(1)
              if (ratio <= acce_prob) {
                n <- n + 1
              }
            }
          }
        }
        # additional move 3: only change boolean functions among gene triplets  note: the network may remain the same
        if (uu >= 6 / 7) {
          prop.prob <- runif(1)
          aa <- 0
          if (prop.prob >= prop.ratio && nrow(triplet.prior.set) >= 1 && ncol(triplet.prior.set) > 1) {
            candidate.pare.set <- list()
            for (jj in seq_len(nrow(triplet.prior.set))) {
              if (length(intersect(parent_of_update, triplet.prior.pare[jj, ])) == 2) {
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
            }
          }
          prop_incid_matrix <- current_incid_matrix
          prop_ances_matrix <- update_ancestor_matrix(prop_incid_matrix)
          prop_trans_func_matrix <- current_trans_func_matrix
          if (aa > 0) {
            for (i in 1:2) {
              prop_trans_func_matrix[update_order[k], parent_of_update[i]] <- add_parent[4]
            }
          }
          if (prop.prob < prop.ratio || aa == 0) {
            func_order <- sample.int(10, 1)
            for (i in 1:2) {
              prop_trans_func_matrix[update_order[k], parent_of_update[i]] <- func_order
            }
          }
          xxx <- Error_LLH(GeneData = GeneData, SampleSize = SampleSize, num.node = num.node, prior_para = prior_para, penalty = penalty, TRFUM = prop_trans_func_matrix, timeseries = timeseries)
          prop_sample <- xxx[[2]]
          prop_post <- xxx[[1]][length(xxx[[1]])]
          nume <- sum(prop_post)
          deno <- sum(current_post)
          acce_prob <- exp(nume - deno)
          ratio <- runif(1)
          if (acce_prob != 1 && ratio <= acce_prob) {
            n <- n + 1
          }
        }
      }
      if (old == n) { # no update
        iter <- iter + 1
        xxx <- Error_LLH(GeneData = GeneData, SampleSize = SampleSize, num.node = num.node, prior_para = prior_para, penalty = penalty, TRFUM = current_trans_func_matrix, timeseries = timeseries)
        current_sample <- xxx[[2]]
        Sample_Matrix[iter, ] <- current_sample
        all_logpost[iter] <- xxx[[1]][length(xxx[[1]])]
        Incidence_Matrix[[iter]] <- current_incid_matrix
        Ancestor_Matrix[[iter]] <- current_ances_matrix
        Trans_Func_Matrix[[iter]] <- current_trans_func_matrix
      }
      if (old != n) { # when either T or F changes,  recount the number
        iter <- iter + 1
        jump_point[n] <- iter

        logpost[n] <- sum(prop_post)
        Sample_Matrix[iter, ] <- prop_sample
        all_logpost[iter] <- sum(prop_post)
        Incidence_Matrix[[iter]] <- prop_incid_matrix
        Ancestor_Matrix[[iter]] <- prop_ances_matrix
        Trans_Func_Matrix[[iter]] <- prop_trans_func_matrix
      }
    } # end of updating
  } # end of num_update
  # calculate burn-in
  burn_in_steps <- floor(burn_in * num_update * num.node)
  # post-burn-in samples
  post_samples <- Trans_Func_Matrix[(burn_in_steps + 2):(num_update * num.node + 1)]
  # calculate the posterior probability of each edge
  post_edge_prob <- Reduce(`+`, lapply(post_samples, function(m) (m > 0) * 1)) / length(post_samples)
  rownames(post_edge_prob) <- rownames(GeneData)
  colnames(post_edge_prob) <- rownames(GeneData)
  if (verbose) {
    close(pb)
    strong_edges <- sum(post_edge_prob > 0.5)
    # print a clean summary block
    cat("\n")
    cat("=========================================\n")
    cat("          BBNI Analysis Summary          \n")
    cat("=========================================\n")
    cat(sprintf("Nodes Analyzed:          %d\n", num.node))
    cat(sprintf("Samples Processed:       %d\n", SampleSize))
    cat(sprintf("MCMC Iterations:         %d\n", num_update))
    cat(sprintf("Burn-in ratio:           %.2f\n", burn_in))
    cat(sprintf("Final Log-Posterior:     %.3f\n", all_logpost[length(all_logpost)]))
    cat(sprintf("Strong Edges (P > 0.5):  %d\n", strong_edges))
    cat("=========================================\n")
  }
  return(list(
    networks = Trans_Func_Matrix,
    log_posterior = all_logpost,
    post_edge_prob = post_edge_prob,
    burn_in = burn_in
  ))
}
