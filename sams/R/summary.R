#' Calculate the Entropy of a Set Partition
#'
#' @param partition A numeric vector representing a partition of the integers
#'   \eqn{1}, ..., \eqn{n} using cluster labels
#'
#' @return Calculated partition entropy as a numeric vector of length one
#'
#' @export
#'
#' @examples
#'
#' p <- c(0,0,0,1,1,2) # n = 6, 3 unique clusters
#' partitionEntropy(p)

partitionEntropy <- function(partition) {
  if (!is.numeric(partition)) {
    stop("Argument 'partition' must be a numeric vector of cluster labels.")
  }
  x <- table(partition)/length(partition)
  sum(-x*log(x))
}

# ----------------------------------------------------------------------------------------------- #

#' Calculate the Number of Items in the Largest Cluster of a Set Partition
#'
#' @param partition A numeric vector representing a partition of the integers
#'   \eqn{1}, ..., \eqn{n} using cluster labels
#'
#' @return The number of items in the largest cluster of the given partition as
#'   a numeric vector of length one.
#'
#' @export
#' @examples
#'
#' p <- c(0,1,1,1,1,1,2)
#' sizeOfLargestCluster(p)

sizeOfLargestCluster <- function(partition) {
  if (!is.numeric(partition)) {
    stop("Argument 'partition' must be a numeric vector of cluster labels.")
  }
  unname(sort(table(partition), decreasing = TRUE))[1]
}

# ----------------------------------------------------------------------------------------------- #

#' Count the Number of Clusters in a Set Partition
#'
#' @param partition A numeric vector representing a partition of the integers
#'   \eqn{1}, ..., \eqn{n} using cluster labels
#'
#' @return The number of clusters in the given set partition as a numeric vector
#'   of length one.
#'
#' @export
#'
#' @examples
#'
#' p <- c(0,1,1,2,3,2,4,4,2)
#' nClusters(p)

nClusters <- function(partition) {
  if (!is.numeric(partition)) {
    stop("Argument 'partition' must be a numeric vector of cluster labels.")
  }
  length(unique(partition))
}

# ----------------------------------------------------------------------------------------------- #

#' Compute the Proportion of Items in Each Cluster for All Partitions
#'
#' @param partitions A matrix, with each row representing a set partition of the
#'   integers \eqn{1}, ..., \eqn{n} as cluster labels
#'
#' @return A matrix whose columns represent the cumulative proportion of the
#'   data that correspond to that cluster.
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
#' nSamples <- 500L
#' partitions <- matrix(0, nrow=nSamples, ncol=length(nealData))
#' for ( i in 2:nSamples ) {
#'   partitions[i,] <- nealAlgorithm3(partitions[i-1,], logPostPredict, mass = 1.0, nUpdates = 2)
#' }
#' clusterProportions(partitions)

clusterProportions <- function(partitions) {
  if (is.null(dim(partitions))) {
    if (!is.numeric(partitions)) {
      stop("Each row of 'partitions' must be a numeric vector of cluster labels.")
    }
      clusterSizes <- cumsum(table(partitions))
  } else {
    if (!is.matrix(partitions)) {
      if (!is.data.frame(partitions)) {
        stop("Argument 'partitions' must be a matrix or data.frame.")
      } else {
        if (any(apply(partitions, 1, is.numeric)) == 0) {
          stop("Each row of 'partitions' must be a numeric vector of cluster labels.")
        }
        partitions <- as.matrix(partitions)
      }
    }
    nSubsets <- apply(partitions, 1, function(x) nClusters(x))
    clusterSizes <- matrix(0, ncol = max(nSubsets), nrow = nrow(partitions))
    for (i in 1:nrow(partitions)) {
      insert <- sort(table(partitions[i,]), decreasing = TRUE)
      clusterSizes[i,] <- c(insert, rep(0, max(nSubsets) - length(insert)))
    }
    clusterSizes <- t(sapply(1:nrow(clusterSizes),
                             function(k) cumsum(clusterSizes[k,])/ncol(partitions)))
  }
  clusterSizes
}

# ----------------------------------------------------------------------------------------------- #

#' Plot Traces of Cluster Sizes
#' @import graphics
#'
#' @param partitions A matrix, with each row a numeric vector cluster labels
#' @param plot.cols A character vector of valid color names, whose length
#'   represents the maximum number of stacked traces to be plotted
#' @param plot.title A character string to be used as the main title on the
#'   trace plot
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
#' nSamples <- 500L
#' partitions <- matrix(0, nrow=nSamples, ncol=length(nealData))
#' for ( i in 2:nSamples ) {
#'   partitions[i,] <- nealAlgorithm3(partitions[i-1,], logPostPredict, mass = 1.0, nUpdates = 2)
#' }
#'
#' clusterTrace(partitions, plot.title = "Neal (2000) Data")

clusterTrace <- function(partitions,
                         plot.cols = rep("black", ncol(partitions)),
                         plot.title = "") {
  nSubsets <- apply(partitions, 1, function(x) nClusters(x))
  clusterSizes <- clusterProportions(partitions)
  plot(clusterSizes[,1], type = 'l', ylim = c(0,1), main = plot.title,
       ylab = "Proportion of Items in Cluster", col = plot.cols[1])
  nLines <- min(length(plot.cols), ncol(clusterSizes))
  for (i in 2:nLines) {
    lines(clusterSizes[,i], col = plot.cols[i])
  }
}

# ------------------------------------------------------------------------------------------------ #

#' Compute the Posterior Pairwise Similarity for All Pairs of Items
#'
#' @param partitions A matrix, with each row a numeric vector cluster labels
#'
#' @return A symmetric matrix of pairwise similarities based on the partitions
#'   given.
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
#' nSamples <- 500L
#' partitions <- matrix(0, nrow=nSamples, ncol=length(nealData))
#' for ( i in 2:nSamples ) {
#'   partitions[i,] <- nealAlgorithm3(partitions[i-1,], logPostPredict, mass = 1.0, nUpdates = 2)
#' }
#'
#' psm(partitions)

psm <- function(partitions) {
  if (!is.matrix(partitions)) partitions <- matrix(partitions, nrow = 1)
  n <- ncol(partitions)
  pairs <- expand.grid(1:n, 1:n)
  out.vec <- apply(pairs, 1, function(x) mean(partitions[,x[1]] == partitions[,x[2]]))
  matrix(out.vec, byrow = TRUE, ncol = n, nrow = n)
}
