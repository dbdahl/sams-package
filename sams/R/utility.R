#' Check if a Vector of Cluster Labels is in Canonical Form
#'
#' @param partition A numeric vector representing a partition of the integers
#'   \eqn{1}, ..., \eqn{n} using cluster labels
#'
#' @return Logical, indicating whether \code{partition} is in canonical form.

isCanonical <- function(partition) {
  u <- unique(partition)
  if ( min(u) != 0L ) return(FALSE)
  if ( max(u) != length(u)-1 ) return(FALSE)
  !any(diff(u) < 0)
}

# ----------------------------------------------------------------------------------------------- #

#' Coerce a Vector of Cluster Labels to Canonical Form
#'
#' @param partition A numeric vector representing a set partition of the
#'   integers \eqn{1}, ..., \eqn{n} using cluster labels
#'
#' @return A numeric vector representing \code{partition}, but now in canonical
#'   form.

asCanonical <- function(partition) {
  if (!is.numeric(partition)) {
    stop("Argument 'partition' must be a numeric vector of cluster labels.")
  }
  temp <- integer(length(partition))
  i <- 0
  for (s in unique(partition)) {
    temp[which(partition == s)] <- i
    i <- i + 1
  }
  temp
}

# ----------------------------------------------------------------------------------------------- #

#' Identify Which Cluster Contains a Given Item
#'
#' @param i Item index as an integer vector of length one
#' @param partition Set partition of the integers \eqn{1}, ..., \eqn{n}
#'   represented as either a numeric vector of cluster labels, or a list
#'   containing subsets of these integers
#'
#' @return A list consisting of \describe{ \item{which}{An integer representing
#'   which cluster \code{i} belongs to} \item{cluster}{The subset of indices
#'   that correspond to the same cluster as \code{i}} }

clusterWithItem <- function(i, partition) UseMethod("clusterWithItem", partition)

clusterWithItem.numeric <- function(i, partition) {
  list(which = partition[i], cluster = which(partition == partition[i]))
}

clusterWithItem.list <- function(i, partition) {
  index <- which(sapply(1:length(partition), function(x) i %in% partition[[x]]))
  list(which = index, cluster = partition[[index]])
}

# ----------------------------------------------------------------------------------------------- #

#' Create a New Cluster with Given Item
#'
#' @param i Item index as an integer vector of length one
#' @param partition Set partition of the integers \eqn{1}, ..., \eqn{n}
#'   represented as either a numeric vector of cluster labels, or a list
#'   containing subsets of these integers
#' @return Updated partition with a new cluster.

createNewCluster <- function(i, partition) UseMethod("createNewCluster", partition)

createNewCluster.numeric <- function(i, partition) {
  partition[i] <- max(unique(partition))+1
  asCanonical(partition)
}

createNewCluster.list <- function(i, partition) {
  removeFrom <- clusterWithItem(i, partition)$which
  partition[[removeFrom]] <- partition[[removeFrom]][partition[[removeFrom]] != i]
  partition[[length(partition) + 1]] <- i
  partition[which(lapply(partition, length) > 0)]
}

# ----------------------------------------------------------------------------------------------- #

#' Join Item to an Existing Cluster
#'
#' @param i Item index as an integer vector of length one
#' @param join Label or index of cluster that \code{i} must join
#' @param partition Set partition of the integers \eqn{1}, ..., \eqn{n}
#'   represented as either a numeric vector of cluster labels, or a list
#'   containing subsets of these integers
#'
#' @return Updated partition.

joinExistingCluster <- function(i, join, partition) UseMethod("joinExistingCluster", partition)

joinExistingCluster.numeric <- function(i, join, partition) {
  partition[i] <- join
  partition
}

joinExistingCluster.list <- function(i, join, partition) {
  removeFrom <- clusterWithItem(i, partition)$which
  partition[[removeFrom]] <- partition[[removeFrom]][partition[[removeFrom]] != i]
  partition[[join]] <- sort(c(partition[[join]], i))
  partition[which(lapply(partition, length) > 0)]
}

# ----------------------------------------------------------------------------------------------- #

#' Coerce a Set Partition in List Structure to Numeric Vectors of Cluster Label
#'
#' @param partition A list representing a set partition of the integers \eqn{1},
#'   ..., \eqn{n}
#' @return A numeric vector representing the set partition using cluster labels.

asClusterLabels <- function(partition) {
  if (is.numeric(partition)) return(asCanonical(partition))
  nItems <- length(unlist(partition))
  partition <- sapply(1:nItems, function(x) clusterWithItem(x, partition)$which)
  asCanonical(partition)
}

# ----------------------------------------------------------------------------------------------- #

#' Coerce a Set Partition as Numeric Vectors of Cluster Labels to a List
#' Structure
#'
#' @param partition A numeric vector representing a partition of the integers
#'   \eqn{1}, ..., \eqn{n} using cluster labels
#' @return The set partition in a list structure.
#'

asSetPartition <- function(partition) {
  if (!isCanonical(partition)) partition <- asCanonical(partition)
  lapply(unique(partition), function(x) which(partition == x))
}

# ----------------------------------------------------------------------------------------------- #

#' Get theta Parameters from a Numeric Vector of Cluster Labels and Unique phi
#' Values
#'
#' @param partition A numeric vector representing a partition of the integers
#'   \eqn{1}, ..., \eqn{n} using cluster labels
#' @param phi A list of unique model parameters whose length must equal the
#'   number of unique cluster labels in \code{partition}
#'
#' @return A numeric vector of model parameters \eqn{theta_1}, ...,
#'   \eqn{theta_n}.

getThetas <- function(partition, phi) {
  if (!isCanonical(partition)) {
    partition <- asCanonical(partition)
  }
  if (length(phi) != length(unique(partition))) {
    stop("Length of 'phi' must equal number of unique cluster labels in 'partition'.")
  }
  thetas <- lapply(1:length(partition), function(x) {
    phi[[partition[x] + 1]]
  })
  thetas
}

# ----------------------------------------------------------------------------------------------- #

#' Enumerate Transformed Weights for Choosing i and j Non-Uniformly
#'
#' @param m A square matrix of pairwise similarilities between items
#' @param fn A function that maps pairwise similarities. Default is the identity
#'   function.
#' @param eps A numeric value close to 0 to give some nonzero weight to pairs of
#'   items with 0 or 1 pairwise similarity
#'
#' @export

transformedWeights <- function(m, fn = function(x) x, eps = 1e-12) {
  if (!is.function(fn)) {
    stop("Function argument 'fn' must be of type closure.")
  }
  dims <- dim(m)
  if (dims[1] != dims[2]) {
    stop("Matrix m must be square.")
  }
  n <- dims[1]
  indices <- expand.grid(1:n, 1:n)
  indices$keep <- as.vector(upper.tri(m))
  indices$wt <- as.vector(m)
  colnames(indices) <- c("i", "j", "keep", "wt")
  indices$wt <- sapply(indices$wt, function(i) max(eps, fn(i)))
  return(subset(indices, indices$keep == TRUE, select = c(1,2,4)))
}

# ------------------------------------------------------------------------------------------------ #

#' Compute the Pochhammer Symbol (Rising Factorials) With Increment
#'
#' @param x Non-negative numeric value
#' @param y Non-negative real value representing increment parameter for
#'   Pochhammer function. If \code{NULL}, there is no increment (i.e.
#'   \code{y=1}).
#' @param n Non-negative integer representing subscript in Pochhammer symbol
#' @param log Logical value indicating whether to return results on log scale
#'
#' @return A numeric value indicating the result of Pochhammer function.
#'
#' @export
#'
#' @examples
#' # effect of increment parameter
#' poch(5, y = NULL, n = 3, log = FALSE)
#' poch(5, y = 1, n = 3, log = FALSE)
#' poch(5, y = 1:4, n = 3, log = FALSE)
#'
#' # increment being NULL is equivalent to ratio of gamma functions
#' a <- 7
#' b <- 3
#' out1 <- poch(a, y = NULL, n = b, log = FALSE)
#' out2 <- gamma(a + b) / gamma(a)
#'
poch <- function(x, y = NULL, n = 1, log = FALSE) {
  if (!log) {
    if (is.null(y)) {
      return(gamma(x + n) / gamma(x))
    } else {
      iter <- 0
      acc <- 1
      while(iter < n) {
        acc <- (x + iter*y) * acc
        iter <- iter + 1
      }
      return(acc)
    }
  } else {
    if (is.null(y)) {
      return(lgamma(x + n) - lgamma(x))
    } else {
      iter <- 0
      acc <- 0
      while(iter < n) {
        acc <- log(x + iter*y) + acc
        iter <- iter + 1
      }
      return(acc)
    }
  }
}

# ------------------------------------------------------------------------------------------------ #

#' Compute Probability Mass of a Partition Under the Two Parameter Chinese
#' Restaurant Process (CRP)
#'
#' @param partition A numeric vector of cluster labels, or a matrix whose rows
#'   are numeric vectors of cluster labels
#' @param mass A numeric value indicating the mass parameter in the CRP, which
#'   must be greater than the \code{-discount} argument
#' @param discount A numeric value on the interval [0,1), indicating the
#'   discount parameter of the two parameter CRP
#' @param log A logical value indicating whether results should be returned on
#'   the log scale
#'
#' @return A numeric vector of probabilities, or log probabilities if \code{log
#'   = TRUE}.
#'
#' @export
#'
#' @examples
#'
#' partitions <- matrix(c(0,0,0,0,0,
#'                        0,0,0,0,1,
#'                        0,0,0,1,2,
#'                        0,0,1,2,3,
#'                        0,1,2,3,4), ncol = 5, nrow = 5, byrow = TRUE)
#'
#' # discount = 0 shows higher probability for lower quantity of components
#' dCRP(partitions, mass = 1, discount = 0, log = FALSE)
#'
#' # discount = 0.5 shows higher probability for higher quantity of components
#' dCRP(partitions, mass = 1, discount = 0.5, log = FALSE)

dCRP <- function(partition, mass = 1.0, discount = 0.0, log = FALSE) {
  if (discount < 0 | discount >= 1) {
    stop("Function argument 'discount' must be on the interval [0,1).")
  }
  if (mass <= -discount) {
    stop("Function argument 'mass' must be strictly greater than then negative of the function argument 'discount'.")
  }
  if (!is.logical(log)) {
    stop("Function argument 'log' must be either TRUE or FALSE.")
  }

  if (!is.matrix(partition)) {
    partition <- matrix(partition, nrow = 1)
  }
  nItems <- ncol(partition)
  nSamples <- nrow(partition)

  mkFunction <- function(m, d, log) {
    function(p) {
      sizes <- sapply(unique(p), function(x) sum(p == x))
      q <- length(sizes)
      if (d == 0) {
        if (log) {
          lgamma(m) + q*log(m) + sum(lgamma(sizes)) - lgamma(m + nItems)
        } else {
          gamma(m) * m^q * prod(gamma(sizes)) / gamma(m + nItems)
        }
      } else {
        if (log) {
          poch(m, y = d, n = q, log = TRUE) +
            sum(poch(1-d, y = NULL, n = sizes - 1, log = TRUE)) -
            poch(m, y = NULL, n = nItems, log = TRUE)
        } else {
          poch(m, y = d, n = q, log = FALSE) *
            prod(poch(1-d, y = NULL, n = sizes - 1, log = FALSE)) /
            poch(m, y = NULL, n = nItems, log = FALSE)
        }
      }
    }
  }

  dCRP.once <- mkFunction(mass, discount, log)
  apply(partition, 1, function(z) dCRP.once(z))
}
