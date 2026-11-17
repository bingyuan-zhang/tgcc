#' Upper end of the lambda sequence
#'
#' @description
#' Finds `lambdaMax`, the smallest value of lambda at which all samples are
#' fused into a single cluster, by a doubling search followed by bisection.
#' The number of clusters at a given lambda is that of the exact solution of
#' the tree-structured problem at this lambda (a single dynamic programming
#' step without cluster fusion), using the same tree and edge weights as
#' [tgCC()].
#'
#' The search starts at `lambdaInit`. If all samples are fused there, lambda is
#' halved until they are not; otherwise lambda is doubled until they are. The
#' interval obtained is then bisected until its relative width is below
#' `relTol`, and the upper end is returned.
#'
#' Some samples, such as outliers with very small edge weights, may never be
#' fused. The doubling therefore stops when the number of clusters has not
#' decreased in `maxStall` consecutive doublings. In this case, `lambdaMax` is
#' the smallest value at which the number of clusters reaches its final value,
#' and a warning is issued.
#'
#' @param data matrix with samples as rows and features as columns.
#' @param bandwidth,useNorm,depthThresh,probThresh,isNaive,mst,famstControl
#'   as in [tgCC()]; use the same values as in the subsequent call of [tgCC()].
#'   For `mst = "famst"`, set `seed` in `famstControl` so that both calls use
#'   the same tree.
#' @param lambdaInit starting value of the search. Default is `nrow(data)`.
#' @param relTol relative width of the final bisection interval. Default is
#'   `0.01`.
#' @param maxStall number of consecutive doublings without a decrease in the
#'   number of clusters after which the search stops. Default is `10`.
#' @return A list containing
#'   \item{lambdaMax}{the upper end of the lambda sequence;}
#'   \item{numClusters}{number of clusters at `lambdaMax` (1 unless some
#'     samples are never fused);}
#'   \item{search}{data frame of the evaluated values of lambda and the
#'     corresponding numbers of clusters.}
#'
#' @examples
#' set.seed(1)
#' dl <- tgcc:::make.mixgaussian(n = 400)
#' res <- lambdaMax(dl$data)
#' res$lambdaMax
#' lambdaSeq <- seq(res$lambdaMax / 100, res$lambdaMax, length.out = 100)
#' tgccFit <- tgCC(dl$data, lambdaSeq)
#'
#' @export
lambdaMax <- function(data,
  bandwidth = NULL,
  useNorm = TRUE,
  depthThresh = 10,
  probThresh = 0.1,
  isNaive = FALSE,
  mst = c("exact", "famst"),
  famstControl = list(nNeighbors = 15, numRandomEdges = 3, metric = "euclidean"),
  lambdaInit = nrow(data),
  relTol = 0.01,
  maxStall = 10) {

  mst <- match.arg(mst)
  stopifnot(lambdaInit > 0, relTol > 0, maxStall >= 1)

  # Tree and edge weights, as in tgCC().
  params <- initParams(data, bandwidth, useNorm, isNaive, mst, famstControl)$params
  edgeWeights <- params$EdgeWeights
  outlierIndex <- (params$Partitions < depthThresh)
  edgeWeights[outlierIndex] <-
    pmax(edgeWeights[outlierIndex], stats::quantile(edgeWeights[outlierIndex], probThresh))

  lambdas <- numeric(0)
  counts <- integer(0)
  numClustersAt <- function(lambda) {
    theta <- matrix(0, nrow(data), ncol(data))
    for (j in seq_len(ncol(data))) {
      theta[, j] <- computeTheta(data[, j], lambda, params$Vertices, params$Types,
        params$Parents, params$NodeWeights, edgeWeights, params$Children)
    }
    k <- nrow(unique(theta))
    lambdas <<- c(lambdas, lambda)
    counts <<- c(counts, k)
    k
  }

  # Doubling search for an interval (lo, hi] with the final number of clusters at hi.
  lo <- 0
  hi <- lambdaInit
  kHi <- numClustersAt(hi)
  stall <- 0
  while (kHi > 1 && stall < maxStall) {
    lo <- hi
    hi <- 2 * hi
    kNew <- numClustersAt(hi)
    stall <- if (kNew < kHi) 0 else stall + 1
    kHi <- kNew
  }
  kTarget <- kHi
  if (kTarget > 1) {
    warning(sprintf(paste0("The number of clusters did not decrease in %d consecutive ",
      "doublings; %d clusters remain at lambda = %g."), maxStall, kTarget, hi))
    # Smallest evaluated lambda at which the final number of clusters is reached.
    hi <- min(lambdas[counts == kTarget])
    below <- lambdas[lambdas < hi]
    lo <- if (length(below)) max(below) else 0
  }

  # If lambdaInit is already beyond lambdaMax, halve until the samples are not all fused.
  if (lo == 0) {
    for (i in 1:60) {
      mid <- hi / 2
      if (numClustersAt(mid) <= kTarget) hi <- mid else { lo <- mid; break }
    }
    if (lo == 0)   # e.g. all rows of data are identical
      return(list(lambdaMax = hi, numClusters = kTarget,
        search = data.frame(lambda = lambdas, numClusters = counts)))
  }

  # Bisection on (lo, hi].
  while ((hi - lo) / hi > relTol) {
    mid <- (lo + hi) / 2
    if (numClustersAt(mid) <= kTarget) hi <- mid else lo <- mid
  }

  list(
    lambdaMax = hi,
    numClusters = kTarget,
    search = data.frame(lambda = lambdas, numClusters = counts)
  )
}
