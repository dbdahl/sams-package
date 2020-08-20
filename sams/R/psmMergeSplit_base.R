#' Base Functionality for the psmMergeSplit Function
#'
#' Merge-split proposals for conjugate "Chinese Restaurant Process" (CRP)
#' mixture models using sequentially-allocated elements. Allocation is performed
#' with weights derived from a previously-calculated pairwise similarity matrix.
#'
#' @param partition A numeric vector of cluster labels representing the current
#'   partition.
#' @param psm A matrix of previously-calculated pairwise similarity
#'   probabilities for each pair of data indices.
#' @param logPosteriorPredictiveDensity A function taking an index \eqn{i} (as a
#'   numeric vector of length one) and a subset of integers \eqn{subset}, and
#'   returning the natural logarithm of \eqn{p( y_i | y_subset )}, i.e., that
#'   item's contribution to the log integrated likelihood given a subset of the
#'   other items. The default value "turns off" the likelihood, resulting in
#'   prior simulation (rather than posterior simulation).
#' @param mass A specification of the mass (concentration) parameter in the CRP
#'   prior. Must be greater than the \code{-discount} argument.
#' @param discount A numeric value on the interval [0,1) corresponding to the
#'   discount parameter in the two-parameter CRP prior.
#' @param nUpdates An integer giving the number of merge-split proposals before
#'   returning. This has the effect of thinning the Markov chain.
#' @param selectionWeights A matrix or data frame whose first two columns are
#'   the unique pairs of data indices, along with a column of weights
#'   representing how likely each pair is to be selected at the beginning of
#'   each merge-split update.
#'
#' @return \describe{ \item{partition}{A numeric vector giving the updated
#'   partition encoded using cluster labels.} \item{accept}{The acceptance rate
#'   of the Metropolis-Hastings proposals, i.e. the number of accepted proposals
#'   divided by \code{nUpdates}.} }
#'
#' @seealso \code{\link{psmMergeSplit}}
#'
#' @import stats

psmMergeSplit_base <- function(partition,
                               psm,
                               logPosteriorPredictiveDensity = function(i, subset) 0.0,
                               mass = 1.0,
                               discount = 0.0,
                               nUpdates = 1L,
                               selectionWeights = NULL) {
  # check function arguments
  if (!is.function(logPosteriorPredictiveDensity)) {
    stop("Function argument 'logPosteriorPredictiveDensity' must be of type 'closure'")
  }
  if (discount < 0 | discount >= 1) {
    stop("Function argument 'discount' must be on the interval [0,1).")
  }
  if (mass <= -discount) {
    stop("Function argument 'mass' must be strictly greater than the negative of the function argument 'discount'.")
  }
  if (!is.integer(nUpdates)) {
    nUpdates <- as.integer(nUpdates)
    if (nUpdates < 1) {
      stop("Function argument 'nUpdates' must be a positive integer.")
    }
  }
  if (!isCanonical(partition)) {
    partition <- asCanonical(partition)
  }
  nItems <- length(partition)
  # define sampling mechanism for (i,j) pair
  if (is.null(selectionWeights)) {
    samplePair <- function() sample(nItems, 2, replace = FALSE)
  } else {
    samplePair <- function() {
      as.integer(selectionWeights[sample(nrow(selectionWeights), 1, prob = selectionWeights[,3]), 1:2])
    }
  }

  psmWeights <- function(k, set) {
    wt <- c(0,0)
    clusterSizes <- sapply(set, length)
    t <- sum(clusterSizes)
    num <- sapply(1:2, function(x) sum(psm[set[[x]], k]))
    wts <- num / sum(num)
    wts
  }

  mkLogPriorRatio <- function(d) {
    if (d == 0) {
      function(doSplit) {
        if (doSplit) {
          log(mass) + lfactorial(n_si_split - 1) + lfactorial(n_sj_split - 1) - lfactorial(n_si -1)
        } else {
          lfactorial(n_si_merge - 1) - log(mass) - lfactorial(n_si - 1) - lfactorial(n_sj - 1)
        }
      }
    } else {
      function(doSplit) {
        if (doSplit) {
          log(mass + d*q) + lgamma(n_si_split - d) + lgamma(n_sj_split - d) -
            lgamma(1-d) - lgamma(n_si - d)
        } else {
          lgamma(n_si_merge - d) + lgamma(1-d) -
            log(mass + d*(q-1)) - lgamma(n_si - d) - lgamma(n_sj - d)
        }
      }
    }
  }

  logPriorRatio <- mkLogPriorRatio(discount)

  accept <- 0
  for (u in 1:nUpdates) {
    q <- length(unique(partition))
    ijPair <- samplePair()
    clusterForI <- clusterWithItem(ijPair[1], partition)
    clusterForJ <- clusterWithItem(ijPair[2], partition)
    doSplit <- clusterForI$which == clusterForJ$which
    # propose a split of i and j belong to the same cluster
    if (doSplit) {
      s <- which(partition == clusterForI$which)
      clusterForJ$which <- max(unique(partition)) + 1
      s <- s[!s %in% ijPair]
      s_i <- ijPair[1]
      s_j <- ijPair[2]
      n_s <- length(s)

      # randomly and uniformly allocate remaining items with i or j
      if (n_s > 0) {
        if (n_s > 1) {
          permuteS <- sample(s)
        } else {
          permuteS <- s
        }
        q_k <- numeric(n_s)
        for (k in 1:n_s) {
          wts <- psmWeights(permuteS[k], list(s_i, s_j))
          if (sum(wts) == 0) wts <- c(1,1)
          chooseThisCluster <- sample(1:2, 1, prob = wts)
          if (chooseThisCluster == 1) {
            s_i <- sort(c(s_i, permuteS[k]))
          } else {
            s_j <- sort(c(s_j, permuteS[k]))
          }
          q_k[k] <- wts[chooseThisCluster]/sum(wts)
        }
      } else {
        q_k <- 1
      }

      # get proposed state
      proposedPartition <- partition
      proposedPartition[s_i] <- clusterForI$which
      proposedPartition[s_j] <- clusterForJ$which

      # calculations for MH ratio on log scale
      si_split <- clusterWithItem(ijPair[1], proposedPartition)$cluster
      n_si_split <- length(si_split)
      ik_split <- sapply(1:n_si_split, function(i) logPosteriorPredictiveDensity(si_split[i], si_split[0:(i-1)]))

      sj_split <- clusterWithItem(ijPair[2], proposedPartition)$cluster
      n_sj_split <- length(sj_split)
      jk_split <- sapply(1:n_sj_split, function(j) logPosteriorPredictiveDensity(sj_split[j], sj_split[0:(j-1)]))

      si <- clusterWithItem(ijPair[1], partition)$cluster
      n_si <- length(si)
      ik <- sapply(1:n_si, function(i) logPosteriorPredictiveDensity(si[i], si[0:(i-1)]))

      pRatio <- logPriorRatio(doSplit)
      lRatio <- sum(ik_split) + sum(jk_split) - sum(ik)
      qRatio <- -sum(log(q_k))

      mhRatio <- qRatio + pRatio + lRatio
      if (log(runif(1)) < mhRatio) {
        partition <- asCanonical(proposedPartition)
        accept <- accept + 1
      }
    } else { # propose a merge if i != j
      # propose a merge if i and j are in different components
      s <- union(clusterForI$cluster, clusterForJ$cluster)
      proposedPartition <- partition
      proposedPartition[s] <- partition[ijPair[1]]
      # imaginary split using permuted indices for transition probabilities
      s_split <- s[!s %in% ijPair] # remove i and j from s
      s_i <- ijPair[1] # create singleton sets for i and j
      s_j <- ijPair[2]
      n_s <- length(s_split)

      # randomly and uniformly allocate remaining items with i or j
      if (n_s > 0) {
        if (n_s > 1) {
          permuteS <- sample(s_split)
        } else {
          permuteS <- s_split
        }
        q_k <- numeric(n_s)
        for (k in 1:n_s) {
          wts <- psmWeights(permuteS[k], list(s_i, s_j))
          if (sum(wts) == 0) wts <- c(1,1)
          chooseThisCluster <- sample(1:2, 1, prob = wts)
          if (chooseThisCluster == 1) {
            s_i <- c(s_i, permuteS[k])
          } else {
            s_j <- c(s_j, permuteS[k])
          }
          q_k[k] <- wts[chooseThisCluster]/sum(wts)
        }
      } else {
        q_k <- 1
      }


      # mh ratio calculations on log scale
      si_merge <- sort(union(clusterForI$cluster, clusterForJ$cluster))
      n_si_merge <- length(si_merge)
      ik_merge <- sapply(1:n_si_merge, function(i) logPosteriorPredictiveDensity(si_merge[i], si_merge[0:(i-1)]))

      si <- clusterForI$cluster
      n_si <- length(si)
      ik <- sapply(1:n_si, function(i) logPosteriorPredictiveDensity(si[i], si[0:(i-1)]))

      sj <- clusterForJ$cluster
      n_sj <- length(sj)
      jk <- sapply(1:n_sj, function(j) logPosteriorPredictiveDensity(sj[j], sj[0:(j-1)]))

      pRatio <- logPriorRatio(doSplit)
      lRatio <- sum(ik_merge) - sum(ik) - sum(jk)
      qRatio <- sum(log(q_k))

      mhRatio <- pRatio + lRatio + qRatio
      if (log(runif(1)) < mhRatio) {
        partition <- asCanonical(proposedPartition)
        accept <- accept + 1
      }
    }
  }
  list(partition = partition, accept = accept/nUpdates)
}
