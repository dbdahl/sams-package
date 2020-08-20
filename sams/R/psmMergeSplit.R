#' Merge-Split Sampling for a Partition Based on Sequential Allocation Informed
#' by Pairwise Similarities
#'
#' Merge-split proposals for conjugate "Chinese Restaurant Process" (CRP)
#' mixture models using sequentially-allocated elements. Allocation is performed
#' with weights derived from a previously-calculated pairwise similarity matrix,
#' and optionally complemented with "restricted Gibbs" scans as discussed in
#' Jain & Neal (2004).
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
#' @param t A non-negative integer indicating the number of restricted Gibbs
#'   scans to perform for each merge/split proposal.
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
#' @references Jain, S., & Neal, R. M. (2004). A split-merge Markov chain Monte
#' Carlo procedure for the Dirichlet process mixture model. \emph{Journal of
#' computational and Graphical Statistics}, 13(1), 158-182.
#'
#' @import stats
#'
#' @export
#'
#' @examples
#'
#' # Neal (2000) model and data
#' nealData <- c(-1.48, -1.40, -1.16, -1.08, -1.02, 0.14, 0.51, 0.53, 0.78)
#' mkLogPosteriorPredictiveDensity <- function(data = nealData,
#'                                             sigma2 = 0.1^2,
#'                                             mu0 = 0,
#'                                             sigma02 = 1) {
#'   function(i, subset) {
#'     posteriorVariance <- 1 / ( 1/sigma02 + length(subset)/sigma2 )
#'     posteriorMean <- posteriorVariance * ( mu0/sigma02 + sum(data[subset])/sigma2 )
#'     posteriorPredictiveSD <- sqrt(posteriorVariance + sigma2)
#'     dnorm(data[i], posteriorMean, posteriorPredictiveSD, log=TRUE)
#'   }
#' }
#'
#' logPostPredict <- mkLogPosteriorPredictiveDensity()
#'
#' nSamples <- 1100L
#' nBurn <- 100
#' partitions <- matrix(0, nrow=nSamples, ncol=length(nealData))
#'
#' # initial draws to inform similarity matrix
#' for ( i in 2:nBurn ) {
#'   partitions[i,] <- nealAlgorithm3(partitions[i-1,],
#'                                    logPostPredict,
#'                                    mass = 1,
#'                                    nUpdates = 1)
#' }
#'
#' # Generate pairwise similarity matrix from initial draws
#' psm.mat <- psm(partitions[1:nBurn,])
#'
#' accept <- 0
#' for ( i in (nBurn+1):nSamples ) {
#'   ms <- psmMergeSplit(partitions[i-1,],
#'                       psm.mat,
#'                       logPostPredict,
#'                       t = 1,
#'                       mass = 1.0,
#'                       nUpdates = 1)
#'   partitions[i,] <- ms$partition
#'   accept <- accept + ms$accept
#' }
#'
#' accept / (nSamples - nBurn) # post burn-in M-H acceptance rate
#' nSubsets <- apply(partitions, 1, function(x) length(unique(x)))
#' mean(nSubsets)
#' sum(acf(nSubsets)$acf)-1   # Autocorrelation time
#'
#' entropy <- apply(partitions, 1, partitionEntropy)
#' plot.ts(entropy)
#'

psmMergeSplit <- function(partition,
                          psm,
                          logPosteriorPredictiveDensity = function(i, subset) 0.0,
                          t = 1,
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
  if (t < 1) {
    if (t < 0) {
      stop("Function argument 't' must be a non-negative integer.")
    }
    return(psmMergeSplit_base(partition,
                              psm,
                              logPosteriorPredictiveDensity,
                              mass,
                              discount,
                              nUpdates,
                              selectionWeights))
  }
  if (!isCanonical(partition)) {
    partition <- asCanonical(partition)
  }
  # define sampling mechanism for (i,j) pair
  nItems <- length(partition)

  if (is.null(selectionWeights)) {
    samplePair <- function() sample(nItems, 2, replace = FALSE)
  } else {
    samplePair <- function() {
      as.integer(selectionWeights[sample(nrow(selectionWeights), 1,
                                         prob = selectionWeights[,3]), 1:2])
    }
  }

  # function for calculating restricted Gibbs probabilities
  restrictedGibbs <- function(k, rgPartition) {
    containsK <- rgPartition[k]
    rgPartitionNoK <- rgPartition
    rgPartitionNoK[k] <- NA
    clusters <- unique(rgPartition[-k])
    clusters <- sort(clusters[!is.na(clusters)])
    clusterSizes <- sapply(clusters, function(x) sum(x == rgPartitionNoK, na.rm = TRUE))
    wts <- clusterSizes*exp(sapply(clusters, function(x) {
      logPosteriorPredictiveDensity(k, which(rgPartitionNoK == x))}))
    names(wts) <- clusters
    wts
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
  for (updates in 1:nUpdates) {
    q <- length(unique(partition))
    # randomly and uniformly select an (i,j) pair of indices
    ijPair <- sample(nItems, 2, replace = FALSE)
    clusterForI <- partition[ijPair[1]]
    clusterForJ <- partition[ijPair[2]]
    doSplit <- clusterForI == clusterForJ

    # propose a split if i and j belong to same cluster
    if (doSplit) {
      s <- sort(which(partition == clusterForI))
      s <- s[!s %in% ijPair]
      s_i <- ijPair[1]
      s_j <- ijPair[2]
      n_s <- length(s)
      doRestricted <- n_s > 0
      clusterOptions <- sort(c(clusterForI, max(partition) + 1))

      # t restricted Gibbs scans if there are other elements
      if (doRestricted) {
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
        }

        # prepare restricted Gibbs partition
        rgPartition <- rep(NA, nItems)
        rgPartition[s_i] <- clusterOptions[1]
        rgPartition[s_j] <- clusterOptions[2]

        for (i in 1:t) {
          for (k in s) {
            wts <- restrictedGibbs(k, rgPartition)
            if (sum(wts) == 0) wts <- c(1,1)
            clusterForK <- sample(clusterOptions, 1, prob = wts)
            rgPartition[k] <- clusterForK
          }
        }

        # get launch state after an additional restricted Gibbs scan
        launchPartition <- partition
        launchPartition[which(rgPartition == clusterForI)] <- clusterForI
        launchPartition[which(rgPartition == clusterOptions[2])] <- clusterOptions[2]

        # keep track of transition probabilities
        q_k <- numeric(n_s)
        for (k in 1:n_s) {
          wts <- restrictedGibbs(s[k], rgPartition)
          if (sum(wts) == 0) wts <- c(1,1)
          whichClusterForK <- sample(1:2, 1, prob = wts)
          rgPartition[s[k]] <- clusterOptions[whichClusterForK]
          q_k[k] <- wts[whichClusterForK]/sum(wts)
        }

        # get proposal state after t+1 restricted Gibbs scans
        proposedPartition <- launchPartition
        proposedPartition[which(rgPartition == clusterForI)] <- clusterForI
        proposedPartition[which(rgPartition == clusterOptions[2])] <- clusterOptions[2]
      } else {
        proposedPartition <- partition
        proposedPartition[s_i] <- clusterForI
        proposedPartition[s_j] <- clusterOptions[2]
        q_k <- 1 # only one way to split if i and j are singleton sets
      }

      # calculations for MH ratio on log scale
      si_split <- which(proposedPartition == clusterOptions[1])
      n_si_split <- length(si_split)
      ik_split <- sapply(1:n_si_split,
                         function(i) logPosteriorPredictiveDensity(si_split[i], si_split[0:(i-1)]))

      sj_split <- which(proposedPartition == clusterOptions[2])
      n_sj_split <- length(sj_split)
      jk_split <- sapply(1:n_sj_split,
                         function(j) logPosteriorPredictiveDensity(sj_split[j], sj_split[0:(j-1)]))

      si <- which(partition == clusterForI)
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
    } else {
      # propose a merged state if i and j are in different clusters
      inClusterForI <- which(partition == clusterForI)
      inClusterForJ <- which(partition == clusterForJ)
      s <- sort(union(inClusterForI, inClusterForJ))
      s <- s[!s %in% ijPair]
      s_i <- ijPair[1]
      s_j <- ijPair[2]
      n_s <- length(s)
      doRestricted <- n_s > 0
      clusterOptions <- sort(c(clusterForI, clusterForJ))

      if (doRestricted) {
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
        }

        # prepare restricted Gibbs partition
        rgPartition <- rep(NA, nItems)
        rgPartition[s_i] <- clusterOptions[1]
        rgPartition[s_j] <- clusterOptions[2]

        # do t restricted Gibbs scans, keeping track of transition
        # probabilities on the last scan
        launchPartition <- partition
        if (t > 1) {
          for (i in 1:(t-1)) {
            for (k in s) {
              wts <- restrictedGibbs(k, rgPartition)
              if (sum(wts) == 0) wts <- c(1,1)
              clusterForK <- sample(clusterOptions, 1, prob = wts)
              rgPartition[k] <- clusterForK
            }
          }

          # transition probabilities
          q_k <- numeric(n_s)
          for (k in 1:n_s) {
            wts <- restrictedGibbs(s[k], rgPartition)
            if (sum(wts) == 0) wts <- c(1,1)
            whichClusterForK <- sample(1:2, 1, prob = wts)
            rgPartition[s[k]] <- clusterOptions[whichClusterForK]
            q_k[k] <- wts[whichClusterForK]/sum(wts)
          }

          # establish hypothetical launch state
          launchPartition[which(rgPartition == clusterForI)] <- clusterForI
          launchPartition[which(rgPartition == clusterOptions[2])] <- clusterForJ
        } else {
          # transition probabilities
          q_k <- numeric(n_s)
          for (k in 1:n_s) {
            wts <- restrictedGibbs(s[k], rgPartition)
            if (sum(wts) == 0) wts <- c(1,1)
            whichClusterForK <- sample(1:2, 1, prob = wts)
            rgPartition[s[k]] <- clusterOptions[whichClusterForK]
            q_k[k] <- wts[whichClusterForK]/sum(wts)
          }

          # establish hypothetical launch state
          launchPartition[which(rgPartition == clusterForI)] <- clusterForI
          launchPartition[which(rgPartition == clusterOptions[2])] <- clusterForJ
        }
      } else {
        q_k <- 1
      }

      # only one way to merge two clusters
      proposedPartition <- partition
      proposedPartition[inClusterForJ] <- clusterForI

      # mh ratio calculations on log scale
      si_merge <- sort(union(inClusterForI, inClusterForJ))
      n_si_merge <- length(si_merge)
      ik_merge <- sapply(1:n_si_merge,
                         function(i) logPosteriorPredictiveDensity(si_merge[i], si_merge[0:(i-1)]))

      si <- inClusterForI
      n_si <- length(si)
      ik <- sapply(1:n_si, function(i) logPosteriorPredictiveDensity(si[i], si[0:(i-1)]))

      sj <- inClusterForJ
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
  list(partition = asCanonical(partition), accept = accept/nUpdates)
}

