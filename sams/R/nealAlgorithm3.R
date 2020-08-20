#' Conjugate Gibbs Sampler for a Partition
#'
#' Algorithm 3 from Neal (2000) to update the state of a partition based on the
#' "Chinese Restaurant Process" (CRP) prior and a user-supplied log posterior
#' predictive density function, with additional functionality for the two
#' parameter CRP prior.
#'
#' @param partition A numeric vector of cluster labels representing the current
#'   partition.
#' @param logPosteriorPredictiveDensity A function taking an index \eqn{i} (as a
#'   numeric vector of length one) and a subset of integers \eqn{subset}, and
#'   returning the natural logarithm of \eqn{p( y_i | y_subset )}, i.e., that
#'   item's contribution to the log integrated likelihood given a subset of the
#'   other items. The default value "turns off" the likelihood, resulting in
#'   prior simulation (rather than posterior simulation).
#' @param mass A specification of the mass (concentration) parameter in the CRP
#'   prior. Must be greater than the \code{-discount} argument.
#' @param discount A numeric value on the interval [0,1) corresponding to the
#'   discount parameter in the two parameter CRP prior.  Set to zero for the
#'   usual, one parameter CRP prior.
#' @param nUpdates An integer giving the number of Gibbs scans before returning.
#'   This has the effect of thinning the Markov chain.
#'
#' @return A numeric vector giving the updated partition encoded using cluster
#'   labels.
#'
#' @references Neal, R. M. (2000). Markov chain sampling methods for Dirichlet
#' process mixture models. \emph{Journal of computational and graphical
#' statistics}, 9(2), 249-265.
#'
#' @export
#'
#' @examples
#'
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
#' nSamples <- 1000L
#' partitions <- matrix(0, nrow = nSamples, ncol = length(nealData))
#' for (i in 2:nSamples) {
#'   partitions[i,] <- nealAlgorithm3(partitions[i-1,], logPostPredict, mass = 1.0, nUpdates = 1)
#' }
#'
#' # convergence and mixing diagnostics
#' nSubsets <- apply(partitions, 1, function(x) length(unique(x)))
#' mean(nSubsets)
#' sum(acf(nSubsets)$acf) - 1   # Autocorrelation time
#'
#' entropy <- apply(partitions, 1, partitionEntropy)
#' plot.ts(entropy)
#'

nealAlgorithm3 <- function(partition,
                           logPosteriorPredictiveDensity = function(i, subset) 0.0,
                           mass = 1.0,
                           discount = 0.0,
                           nUpdates = 1L) {
  # check function arguments
  if (!is.function(logPosteriorPredictiveDensity)) {
    stop("Function argument 'logPosteriorPredictiveDensity' must be of type 'closure'")
  }
  if (discount < 0 | discount >= 1) {
    stop("Function argument 'discount' must be on the interval [0,1).")
  }
  if (mass <= -discount) {
    stop("Function argument 'mass' must be greater than the negative of the function argument 'discount'.")
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

  mkAllocationWeights <- function(d) {
    if (d == 0) {
      function(partition) {
        partition.noi <- partition
        partition.noi[i] <- NA
        clusters <- unique(partition[-i])
        clusterSizes <- sapply(clusters, function(x) sum(x == partition.noi, na.rm = TRUE))
        logWts <- sapply(c(clusters, max(clusters) + 1),
                         function(x) logPosteriorPredictiveDensity(i, which(partition.noi == x)))
        wts <- c(clusterSizes, mass)*exp(logWts)
        wts
      }
    } else {
      function(partition) {
        partition.noi <- partition
        partition.noi[i] <- NA
        clusters <- unique(partition[-i])
        clusterSizes <- sapply(clusters, function(x) sum(x == partition.noi, na.rm = TRUE))
        q <- length(clusterSizes)
        logWts <- sapply(c(clusters, max(clusters) + 1),
                         function(x) logPosteriorPredictiveDensity(i, which(partition.noi == x)))
        wts <- c(clusterSizes - d, mass + d*q)*exp(logWts)
        wts
      }
    }
  }

  allocationWeights <- mkAllocationWeights(discount)

  # do a full scan of the indices 'nUpdates' times
  nItems <- length(partition)
  for (u in 1:nUpdates) {
    for (i in 1:nItems) {
      clusters <- unique(partition[-i])
      wts <- allocationWeights(partition)
      partition[i] <- sample(c(clusters, max(clusters) + 1), 1, prob = wts)
      partition <- asCanonical(partition)
    }
  }

  # output cluster labels in canonical form
  partition
}
