##  Auxiliary functions

#' Count HMM Parameters
#'
#' Count the number of parameters in an HMM larger than a small number epsilon.
#' Auxiliary function used in \code{\link{partition.bic}} and
#' \code{\link{cluster.bic}}.
#'
#' @param hmm An \code{hmm} object (see \code{\link[seqHMM]{build_hmm}}).
#' @param eps A threshold epsilon for counting parameters.
#' @return The number of parameters larger than \code{eps}.
#' @seealso Used in \code{\link{partition.bic}} and \code{\link{cluster.bic}}.
#' @export
count.parameters <- function(hmm, eps = 1e-3) {
  d <- (hmm$n_states*(hmm$n_symbols-1) - sum(hmm$emission_probs < eps)) +
    (hmm$n_states*(hmm$n_states-1) - sum(hmm$transition_probs < eps)) +
    (length(hmm$initial_probs)-1)  # subtract the parameters that are smaller
  return(d)
}

#' Smooth Probabilities
#'
#' Smooth a vector of probabilities using absolute discounting. Auxiliary
#' function used in \code{\link{smooth.hmm}}.
#'
#' @param probs A vector of raw probabilities.
#' @param smoothing Smoothing parameter for absolute discounting.
#' @return A vector of smoothed probabilities.
#' @seealso Used in \code{\link{smooth.hmm}}.
#' @export
smooth.probabilities <- function(probs, smoothing = 1e-4) {
  is.zero <- (probs == 0)
  if(any(is.zero)) {
    # Only smooth if there are parameters equal to zero
    n.zeros <- sum(is.zero)
    n.nonzeros <- sum(!is.zero)
    probs[is.zero] <- n.nonzeros*smoothing/n.zeros  # divide smoothing amount
                                                    # equally
    probs[!is.zero] <- probs[!is.zero] - smoothing  # subtract smoothing amount
  }
  return(probs)
}

#' Smooth HMM Parameters
#'
#' Smooth the parameters of an HMM using absolute discounting given a threshold
#' epsilon. Auxiliary function used in \code{\link{select.seeds}},
#' \code{\link{assign.clusters}}, and \code{\link{hmm.clust}}.
#'
#' @param hmm A raw \code{hmm} object (see \code{\link[seqHMM]{build_hmm}}).
#' @param smoothing Smoothing parameter for absolute discounting in
#'   \code{\link{smooth.probabilities}}.
#' @return An \code{hmm} object with smoothed probabilities.
#' @seealso Used in \code{\link{select.seeds}}, \code{\link{assign.clusters}},
#'   and main function for the DBHC algorithm \code{\link{hmm.clust}}.
#' @export
smooth.hmm <- function(hmm, smoothing = 1e-04) {
  # Smooth all probabilities
  temp.trans <- t(apply(hmm$transition_probs, 1, smooth.probabilities,
                        smoothing = smoothing))
  temp.emiss <- t(apply(hmm$emission_probs, 1, smooth.probabilities,
                        smoothing = smoothing))
  temp.init <- smooth.probabilities(hmm$initial_probs, smoothing = smoothing)

  # Construct the new hmm object
  new.hmm <- build_hmm(hmm$observations,  transition_probs = temp.trans,
                      emission_probs = temp.emiss,
                      initial_probs = temp.init)
  return(new.hmm)
}

#' Get HMM Log Likelihood
#'
#' Get the log likelihood of an HMM object and check if it is feasible (i.e.,
#' contains no illegal emissions). Auxiliary function used in
#' \code{\link{partition.bic}}.
#'
#' @param hmm An \code{hmm} object (see \code{\link[seqHMM]{build_hmm}}).
#' @return The log likelihood of the \code{hmm} object, print warning if model
#'   is infeasible (i.e., if the log likelihood is evaluated for a sequence that
#'   contains emissions that are assigned probability 0 in the \code{hmm}
#'   object).
#' @seealso Used in \code{\link{partition.bic}}.
#' @export
model.ll <- function(hmm) {
  ll <- logLik(hmm)
  if(is.nan(ll)) {
    # Model is infeasible
    warning("Evaluated NaN ll.")
    ll <- -Inf  # set ll to minus Inf
  }
  return(ll)
}

#' Partition BIC
#'
#' Compute the BIC of a partition given a threshold epsilon for counting
#' parameters. Auxiliary function used in \code{\link{hmm.clust}}.
#'
#' @param partition A list object with the partition of HMMs, a mixture of HMMs.
#' @param eps A threshold epsilon for counting parameters in
#'   \code{\link{count.parameters}}.
#' @return The BIC of the partition.
#' @seealso Used in main function for the DBHC algorithm
#'   \code{\link{hmm.clust}}.
#' @export
partition.bic <- function(partition, eps = 1e-03) {
  K <- length(partition)  # number of clusters
  d <- sum(sapply(partition, count.parameters, eps = eps))  # total number of
                                                            # parameters
  n <- sum(sapply(partition, function(hmm) hmm$n_sequences))  # total number of
                                                              # sequences
  ll.sum <- sum(sapply(partition, model.ll))  # sum log likelihoods
  bic <- (-2)*ll.sum + (K + d)*log(n)  # compute BIC
  return(bic)
}

#' HMM BIC
#'
#' Compute the BIC of a single HMM given a threshold epsilon for counting
#' parameters. Auxiliary function used in \code{\link{size.search}}.
#'
#' @inheritParams count.parameters
#' @return The BIC of \code{hmm}.
#' @seealso Used in \code{\link{size.search}}.
#' @export
cluster.bic <- function(hmm, eps = 1e-03) {
  d <- count.parameters(hmm)
  bic <- (-2)*logLik(hmm) + d*log(hmm$n_sequences)
  return(bic)
}

#' Sequence-to-HMM Likelihood
#'
#' Compute the sequence-to-HMM likelihood of an HMM evaluated for a single
#' sequence and check if the sequence contains emissions that are not possible
#' according to the HMM. Auxiliary function used in \code{\link{select.seeds}}
#' and \code{\link{assign.clusters}}.
#'
#' @param hmm An \code{hmm} object (see \code{\link[seqHMM]{build_hmm}})
#'   containing a single sequence.
#' @return The log likelihood of the sequence contained in \code{hmm}, value
#'   will be set to minus infinity if the sequence contains illegal emissions.
#' @seealso Used in \code{\link{select.seeds}} and
#'   \code{\link{assign.clusters}}.
#' @export
seq2hmm.ll <- function(hmm) {
  if(nrow(hmm$observations) > 1) {
    stop("HMM cannot contain more than one sequence")
  }
  temp.ll <- logLik(hmm)  # get ll of sequence
  if(is.nan(temp.ll)) {
    # Model is infeasible
    ll <- -Inf  # set ll to minus Inf
  } else {
    ll <- temp.ll  # set ll to actual
  }
  return(ll)
}

#' Seed Selection Procedure
#'
#' Seed selection procedure of the DBHC algorithm, also invokes size search
#' algorithm for seed in \link{size.search}. Used in \code{\link{hmm.clust}}.
#'
#' @param sequences An \code{stslist} object (see
#'   \code{\link[TraMineR]{seqdef}}) of sequences with discrete observations.
#' @param log_space Logical, parameter provided to
#'   \code{\link[seqHMM]{fit_model}} for whether to use optimization in log
#'   space or not.
#' @param K The number of seeds to select, equal to the number of clusters in a
#'   partition.
#' @param seed.size Seed size, the number of sequences to be selected for a
#'   seed.
#' @param init.size The number of HMM states in an initial HMM.
#' @param print Logical, whether to print intermediate steps or not.
#' @param smoothing Smoothing parameter for absolute discounting in
#'   \code{\link{smooth.probabilities}}.
#' @return A partition as a list object with HMMs for the selected seeds.
#' @seealso Used in main function for the DBHC algorithm
#'   \code{\link{hmm.clust}}.
#' @export
select.seeds <- function(sequences, log_space = FALSE, K, seed.size = 3,
                         init.size = 2, print = FALSE, smoothing = 1e-04) {
  # Error handling
  if(is.null(sequences)) {
    stop("Need sequences for selecting a seed.")
  }

  partition <- list()
  taken.ids <- NULL
  for(i in 1:K) {
    # Select a seed for each cluster
    if(i == 1) {
      # Randomly select first sequence of first seed
      ids <- rownames(sequences)
      id.selected <- ids[sample(1:length(ids), 1)]
    } else {
      # Select first sequences for seeds two and larger
      ll.lowest <- Inf
      # Get ids of available sequences
      available.ids <- rownames(sequences)[!(rownames(sequences) %in% taken.ids)]
      # Initialize a list with best lls for each sequence
      highest.lls <- matrix(-Inf, nrow = length(available.ids))
      for(k in 1:length(partition)) {
        # Loop over all cluster models to find worst represented sequence
        # Smooth hmm of current cluster
        hmm <- smooth.hmm(build_hmm(sequences[rownames(sequences) %in% available.ids,],
                         transition_probs = partition[[k]]$transition_probs,
                         emission_probs = partition[[k]]$emission_probs,
                         initial_probs = partition[[k]]$initial_probs),
                         smoothing = smoothing)
        ll.list <- logLik(hmm, partials = T)  # get likelihoods for sequences
        ll.list[is.na(ll.list)] <- -Inf  # set lls of sequences for which the
                                         # model is infeasible to minus Inf
        # Update list of best lls
        highest.lls[ll.list > highest.lls] <- ll.list[ll.list > highest.lls]
      }
      ll.lowest <- min(highest.lls)  # find lowest ll
      hmm.id <- available.ids[which(highest.lls == ll.lowest)]
      if(length(hmm.id) > 1) {
        # More than one sequence worst represented
        # Randomly select one sequence from worst represented
        id.selected <- hmm.id[sample(1:length(hmm.id), 1)]
      } else {
        id.selected <- hmm.id  # select worst represented sequence
      }
    }
    selected <- id.selected  # selected in this iteration
    taken.ids <- c(taken.ids, selected)  # taken ids in the entire search

    if(seed.size > 1) {
      # Build initial hmm for selected sequence
      init.hmm <- build_hmm(sequences[id.selected,], n_states = init.size)
      init.fit <- fit_model(init.hmm,
                            log_space = log_space,
                            control_em = list(restart = list(times = 10)))
      temp.hmm <- init.fit$model

      other.hmms <- list()
      other.ids <- rownames(sequences)[!(rownames(sequences) %in% taken.ids)]
      for(j in 1:length(other.ids)) {
        # Construct hmm object for all other sequences using the initial hmm
        other.hmms <- c(other.hmms,
                        list(smooth.hmm(build_hmm(sequences[rownames(sequences) == other.ids[j],],
                                                  transition_probs = temp.hmm$transition_probs,
                                                  emission_probs = temp.hmm$emission_probs,
                                                  initial_probs = temp.hmm$initial_probs),
                                        smoothing = smoothing)))
      }
      # Get lls for all hmm objects and select the (seed.size-1) best represented
      other.lls <- data.frame(other.ids,sapply(other.hmms, seq2hmm.ll))
      sorted.indices <- sort(other.lls[,2], decreasing = T, index.return = T)$ix
      to.add <- as.character(other.lls[sorted.indices[1:(seed.size-1)],1])
      selected <- c(selected, to.add)  # update selected ids for current seed
      taken.ids <- c(taken.ids, to.add)  # update taken ids among all seeds
    }
    # Invoke size search algorithm on seed
    search.obj <- size.search(sequences[rownames(sequences) %in% selected,],
                              log_space = log_space,
                              print = print)
    partition <- c(partition, list(search.obj$hmm))  # update partition
  }
  return(partition)
}

#' Cluster Assignment
#'
#' Assign sequences to cluster models that give the highest sequence-to-hmm
#' likelihood. Used in \code{\link{hmm.clust}}.
#'
#' @param partition A list object with the partition, a mixture of HMMs. Each
#'   element in the list is an \code{hmm} object (see
#'   \code{\link[seqHMM]{build_hmm}}).
#' @param memberships A matrix with cluster memberships for each sequence.
#' @param sequences An \code{stslist} object (see
#'   \code{\link[TraMineR]{seqdef}}) of sequences with discrete observations.
#' @param smoothing Smoothing parameter for absolute discounting in
#'   \code{\link{smooth.probabilities}}.
#' @return The updated matrix with cluster memberships for each sequence.
#' @seealso Used in main function for the DBHC algorithm
#'   \code{\link{hmm.clust}}.
#' @export
assign.clusters <- function(partition, memberships, sequences, smoothing = 1e-04) {
  new.memberships <- memberships
  anyInf <- F  # indicator for infeasible models
  for(id in rownames(sequences)) {
    # Find model with highest ll for each sequence
    highest.k <- 0  # initialize cluster with highest ll with temp value
    highest.ll <- -Inf  # initialize highest ll at minus Inf
    for(k in 1:length(partition)) {
      # Evaluate ll for each cluster model in the partition
      # Build HMM object and smooth parameters
      temp.hmm <- smooth.hmm(build_hmm(sequences[id,],
                                       transition_probs = partition[[k]]$transition_probs,
                                       emission_probs = partition[[k]]$emission_probs,
                                       initial_probs = partition[[k]]$initial_probs),
                             smoothing = smoothing)
      curr.ll <- seq2hmm.ll(temp.hmm)  # get ll
      if(curr.ll > highest.ll) {
        # Update variables for best cluster model
        highest.ll <- curr.ll
        highest.k <- k
      }
    }
    if(highest.k == 0) {
      # None of the cluster models is feasible for current sequence
      anyInf <- T
      highest.k <- sample(1:length(partition), 1)  # randomly assign sequence
    }
    # Update membership matrix
    new.memberships[memberships[,"id"] == id,"cluster"] <- highest.k
  }
  if(anyInf) {
    # Print warning if any of HMMs is infeasible
    print(paste0("WARNING: Infeasible HMM at K=", length(partition)))
  }
  return(new.memberships)
}


#' Size Search Algorithm
#'
#' The size search algorithm finds the optimal number of HMM states for a set of
#' sequences and returns both the optimal \code{hmm} object and the
#' corresponding number of hidden states. Used in \code{\link{select.seeds}}.
#'
#' @param sequences An \code{stslist} object (see
#'   \code{\link[TraMineR]{seqdef}}) of sequences with discrete observations.
#' @param log_space Logical, parameter provided to
#'   \code{\link[seqHMM]{fit_model}} for whether to use optimization in log
#'   space or not.
#' @param print Logical, whether to print intermediate steps or not.
#' @return A list with the optimal number of HMM states and the optimal
#'   \code{hmm} object.
#' @seealso Used in the DBHC seed selection procedure in
#'   \code{\link{select.seeds}}.
#' @export
size.search <- function(sequences, log_space = FALSE, print = FALSE) {
  H <- 2  # initialize size at 2
  # Build initial hmm and estimate parameters
  hmm.init <- build_hmm(sequences, n_states = H)
  curr.fit <- fit_model(hmm.init,
                        log_space = log_space,
                        control_em = list(restart = list(times = 10)))
  curr.hmm <- curr.fit$model
  opt.size.found <- F

  if(print) {
    # Print progress
    print("---------")
    print("Start new size search: ")
  }
  while(!opt.size.found) {
    # Search until BIC does not increase anymore
    if(print) {
      # Print progress
      print(paste0(H, paste0(": ", cluster.bic(curr.hmm))))
      print(paste0("LogLik: ", logLik(curr.hmm)))
    }
    H <- H + 1  # increment number of states
    prev.hmm <- curr.hmm # save previous hmm

    # Build hmm with current size and estimate parameters
    hmm.init <- build_hmm(sequences, n_states = H)
    curr.fit <- fit_model(hmm.init,
                          log_space = log_space,
                          control_em = list(restart = list(times = 10)))
    curr.hmm <- curr.fit$model  # extract model

    if(cluster.bic(curr.hmm) > cluster.bic(prev.hmm)) {
      # BIC did not increase
      opt.size.found <- T
      if(print) {
        # Print progress
        print(paste0(H, paste0(": ", cluster.bic(curr.hmm))))
        print(paste0("LogLik: ", logLik(curr.hmm)))
        print(paste0("Size found: ", prev.hmm$n_states))
      }
    }
  }
  # Return list with number of states and hmm object for optimal model
  return(list(size = prev.hmm$n_states, hmm = prev.hmm))
}


## Main HMM Clustering algorithm

#' DBHC Algorithm
#'
#' Implementation of the DBHC algorithm, an HMM clustering algorithm that finds
#' a mixture of discrete-output HMMs. The algorithm uses heuristics based on BIC
#' to search for the optimal number of hidden states in each HMM and the optimal
#' number of clusters.
#'
#' @param sequences An \code{stslist} object (see
#'   \code{\link[TraMineR]{seqdef}}) of sequences with discrete observations or
#'   a \code{data.frame}.
#' @param id A vector with ids that identify the sequences in \code{sequences}.
#' @param smoothing Smoothing parameter for absolute discounting in
#'   \code{\link{smooth.probabilities}}.
#' @param eps A threshold epsilon for counting parameters in
#'   \code{\link{count.parameters}}.
#' @param init.size The number of HMM states in an initial HMM.
#' @param alphabet The alphabet of output labels, if not provided alphabet is
#'   taken from \code{stslist} object (see \code{\link[TraMineR]{seqdef}}).
#' @param K.max Maximum number of clusters, if not provided algorithm searches
#'   for the optimal number itself.
#' @param log_space Logical, parameter provided to
#'   \code{\link[seqHMM]{fit_model}} for whether to use optimization in log
#'   space or not.
#' @param print Logical, whether to print intermediate steps or not.
#' @param seed.size Seed size, the number of sequences to be selected for a seed
#' @return A list with components: \describe{ \item{\code{sequences}}{An
#'   \code{stslist} object of sequences with discrete observations.}
#'   \item{\code{id}}{A vector with ids that identify the sequences in
#'   \code{sequences}.} \item{\code{cluster}}{A vector with found cluster
#'   memberships for the sequences.} \item{\code{partition}}{A list object with
#'   the partition, a mixture of HMMs. Each element in the list is an \code{hmm}
#'   object.} \item{\code{memberships}}{A matrix with cluster memberships for
#'   each sequence.} \item{\code{n.clusters}}{Numerical, the found number of
#'   clusters.} \item{\code{sizes}}{A vector with the number of HMM states for
#'   each cluster model.} \item{\code{bic}}{A vector with the BICs for each
#'   cluster model.} }
#' @examples
#' ## Simulated data
#' library(seqHMM)
#' output.labels <-  c("H", "T")
#'
#' # HMM 1
#' states.1 <- c("A", "B", "C")
#' transitions.1 <- matrix(c(0.8,0.1,0.1,0.1,0.8,0.1,0.1,0.1,0.8), nrow = 3)
#' rownames(transitions.1) <- states.1
#' colnames(transitions.1) <- states.1
#' emissions.1 <- matrix(c(0.5,0.75,0.25,0.5,0.25,0.75), nrow = 3)
#' rownames(emissions.1) <- states.1
#' colnames(emissions.1) <- output.labels
#' initials.1 <- c(1/3,1/3,1/3)
#'
#' # HMM 2
#' states.2 <- c("A", "B")
#' transitions.2 <- matrix(c(0.75,0.25,0.25,0.75), nrow = 2)
#' rownames(transitions.2) <- states.2
#' colnames(transitions.2) <- states.2
#' emissions.2 <- matrix(c(0.8,0.6,0.2,0.4), nrow = 2)
#' rownames(emissions.2) <- states.2
#' colnames(emissions.2) <- output.labels
#' initials.2 <- c(0.5,0.5)
#'
#' # Simulate
#' hmm.sim.1 <- simulate_hmm(n_sequences = 100,
#'                           initial_probs = initials.1,
#'                           transition_probs = transitions.1,
#'                           emission_probs = emissions.1,
#'                           sequence_length = 25)
#' hmm.sim.2 <- simulate_hmm(n_sequences = 100,
#'                           initial_probs = initials.2,
#'                           transition_probs = transitions.2,
#'                           emission_probs = emissions.2,
#'                           sequence_length = 25)
#' sequences <- rbind(hmm.sim.1$observations, hmm.sim.2$observations)
#' n <- nrow(sequences)
#'
#' # Clustering algorithm
#' id <- paste0("K-", 1:n)
#' rownames(sequences) <- id
#' sequences <- sequences[sample(1:n, n),]
#' res <- hmm.clust(sequences, id)
#'
#'
#' #############################################################################
#'
#' ## Swiss Household Data
#' data("biofam", package = "TraMineR")
#'
#' # Clustering algorithm
#' new.alphabet <- c("P", "L", "M", "LM", "C", "LC", "LMC", "D")
#' sequences <- seqdef(biofam[,10:25], alphabet = 0:7, states = new.alphabet)
## Code below not tested because it takes a lot of time to run #################
#' \dontrun{res <- hmm.clust(sequences)
#'
#' # Heatmaps
#' cluster <- 1  # display heatmaps for cluster 1
#' transition.heatmap(res$partition[[cluster]]$transition_probs,
#'                    res$partition[[cluster]]$initial_probs)
#' emission.heatmap(res$partition[[cluster]]$emission_probs)
#' }
#'
#' ## A smaller example, which takes less time to run
#' subset <- sequences[sample(1:nrow(sequences), 20, replace = FALSE),]
#'
#' # Clustering algorithm
#' res <- hmm.clust(subset, K.max = 3)
#'
#' # Number of clusters
#' print(res$n.clusters)
#'
#' # Table of cluster memberships
#' table(res$memberships[,"cluster"])
#'
#' # BIC for each number of clusters
#' print(res$bic)
#'
#' # Heatmaps
#' cluster <- 1  # display heatmaps for cluster 1
#' transition.heatmap(res$partition[[cluster]]$transition_probs,
#'                    res$partition[[cluster]]$initial_probs)
#' emission.heatmap(res$partition[[cluster]]$emission_probs)
#' @export
hmm.clust <- function(sequences, id = NULL, smoothing = 1e-04, eps = 1e-03,
                      init.size  = 2, alphabet = NULL, K.max = NULL,
                       log_space = FALSE, print = FALSE, seed.size = 3) {
  # Error handling
  if(!("stslist" %in% class(sequences))) {
    # User did not provide stslist object
    if(class(sequences) == "data.frame") {
      # Convert data.frame to stslist
      sequences <- seqdef(sequences)
    } else {
      # User did not provide valid object
      stop(paste0("sequences must be a sequence data object of type stslist or",
                  " must be a data.frame"))
    }
  }
  if(is.null(id)) {
    # Take rownames data as id
    id <- rownames(sequences)
  } else {
    if(length(id) != nrow(sequences)) {
      stop("id and sequences need to have equal lengths.")
    }
    if(length(unique(id)) != length(id)) {
      stop("sequences must have unique ids.")
    }
  }
  if(eps < 0 | eps >= 1) {
    stop("eps must be valid threshold for counting probabilities: 0 <= eps < 1")
  }
  if(smoothing > 1e-01) {
    stop("smoothing parameter must be a small number: smoothing <= 1e-01")
  }
  if(init.size < 2) {
    stop("initial hmm must have at least two hidden states")
  }
  if(!is.null(alphabet)) {
    alphabet(sequences) <- alphabet
  }

  # Initialization
  n <- nrow(sequences)
  K <- 1  # initialize number of clusters at 1
  done <- F
  bic.vec <- NULL
  curr.memberships <- cbind(id = id, cluster = rep(1, n))  # initialize matrix
                                                           # with memberships
  rownames(sequences) <- id  # set row names of stslist object to ids

  # Seed selection procedure for first partition
  temp.partition <- select.seeds(sequences = sequences,
                                 log_space = log_space,
                                 K = K,
                                 seed.size = seed.size,
                                 print = print,
                                 smoothing = smoothing,
                                 init.size = init.size)
  # Create hmm object for all sequences and smooth hmm parameters
  temp.hmm <- smooth.hmm(
            build_hmm(sequences,
                      transition_probs = temp.partition[[1]]$transition_probs,
                      emission_probs = temp.partition[[1]]$emission_probs,
                      initial_probs = temp.partition[[1]]$initial_probs),
            smoothing = smoothing)
  curr.partition <- list(temp.hmm)  # initialize partition list
  curr.bic <- partition.bic(curr.partition, eps = eps)  # store BIC
  bic.vec <- c(bic.vec, curr.bic)

  while(!done) {
    # Keep expanding partition until BIC does not increase anymore
    if(print) {
      print(paste0("Number of clusters: ", K))
      print(paste0("BIC: ", curr.bic))
    }
    K <- K + 1  # increment number of clusters
    prev.memberships <- curr.memberships  # store memberships
    prev.partition <- curr.partition  # store partition

    # Select seeds for all clusters
    curr.partition <- select.seeds(sequences = sequences,
                                   log_space = log_space,
                                   K = K,
                                   seed.size = seed.size,
                                   print = print,
                                   smoothing = smoothing,
                                   init.size = init.size)

    # Distribute sequences over cluster models
    curr.memberships <- assign.clusters(curr.partition, curr.memberships, sequences)
    curr.clusters <- curr.memberships[,"cluster"]

    changes <- T
    while(changes) {
      # Keep alternating until no sequence changes
      prev.clusters <- curr.clusters

      # Re-estimate and smooth hmms for all clusters
      for(k in 1:length(curr.partition)) {
        temp.hmm <- smooth.hmm(
          build_hmm(sequences[curr.memberships[,"cluster"] == k,],
                    transition_probs = curr.partition[[k]]$transition_probs,
                    emission_probs = curr.partition[[k]]$emission_probs,
                    initial_probs = curr.partition[[k]]$initial_probs),
          smoothing = smoothing)
        fit.hmm <- fit_model(temp.hmm, log_space = log_space)
        curr.partition[[k]] <- smooth.hmm(fit.hmm$model,
                                          smoothing = smoothing)  # smooth hmm parameters
      }

      # Re-assign clusters to cluster models
      curr.memberships <- assign.clusters(curr.partition, curr.memberships, sequences)
      curr.clusters <- curr.memberships[,"cluster"]
      for(k in 1:length(curr.partition)) {
        # Update partition list according to cluster assignment and smooth hmm
        # parameters of all cluster models
        temp.hmm <- curr.partition[[k]]
        curr.partition[[k]] <- smooth.hmm(build_hmm(sequences[curr.memberships[,"cluster"] == k,],
                              transition_probs = temp.hmm$transition_probs,
                              emission_probs = temp.hmm$emission_probs,
                              initial_probs = temp.hmm$initial_probs),
                              smoothing = smoothing)
      }

      if(sum(prev.clusters == curr.clusters) == length(curr.clusters)) {
        # Check if cluster membership changed for any sequence
        changes <- F
      }
    }

    prev.bic <- curr.bic  # store previous BIC
    curr.bic <- partition.bic(curr.partition, eps = eps)  # get current BIC
    bic.vec <- c(bic.vec, curr.bic)
    if(is.null(K.max)) {
      if(curr.bic >= prev.bic) {
        # Optimal partition found
        done <- T
        if(print) {
          print(paste0("Number of clusters found: ", length(prev.partition)))
        }
      }
    } else {
      if(K == K.max) {
        # Max number of clusters reached
        done <- T
        prev.partition <- curr.partition
        prev.memberships <- curr.memberships
        if(print) {
          print(paste0("Maximum number of clusters reached: ", K.max))
          print(paste0("BIC: ", curr.bic))
        }
      }
    }
  }

  final.partition <- prev.partition  # store final partition
  n.states <- sapply(final.partition, function(hmm) hmm$n_states) # get number of
                                                                  # states
  names(bic.vec) <- c(1:length(final.partition))

  # Return list with final partition, data, and diagnostic variables
  return(list(sequences = sequences,
              id = id,
              cluster = prev.memberships[,"cluster"],
              partition = final.partition,
              memberships = prev.memberships,
              n.clusters = length(final.partition),
              sizes = n.states,
              bic = bic.vec))
}
