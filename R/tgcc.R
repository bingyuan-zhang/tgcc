#' Tree-Guided Convex Clustering
#'
#' @description
#' Perform the tree-guided convex clustering using dynamic programming method.
#'
#' @param data matrix with samples as rows and features as columns.
#' @param lambdaSeq tuning parameters sequence.
#' @param bandwidth bandwidth for the Gaussian kernel between rows. Default is heuristic assigned.
#' @param useNorm logical, whether to use normalized distance. Default is `TRUE`.
#' @param depthThresh depth threshold target. Default is `10`.
#' @param probThresh probability threshold for edge weight. Default is `0.1`.
#' @param isNaive use naive method in searching the exact MST. Default is `FALSE`.
#' @param mst method used to construct the tree: `"exact"` for the exact
#'   Euclidean MST computed by `mlpack::emst`, or `"famst"` for the
#'   approximate MST computed by [famst()]. Default is `"exact"`.
#' @param famstControl list of arguments passed to [famst()] when
#'   `mst = "famst"`: `nNeighbors` (default `15`), `numRandomEdges`
#'   (default `3`), `metric` (default `"euclidean"`) and `seed`.
#' @param keepTheta logical, whether to return the fitted matrix for every
#'   value of lambda (`theta`). Each element is an n x p matrix, so set this
#'   to `FALSE` for large data when only the cluster labels are needed;
#'   [clusterLabel()] does not use `theta`. Default is `TRUE`.
#' @return A list containing the clustered results.
#'
#' @useDynLib tgcc, .registration=TRUE
#' @importFrom stats quantile
#'
#' @examples
#'
#' library(tgcc)
#' dl <- tgcc:::make.mixgaussian(n = 400)
#' data <- dl$data
#' label <- dl$label
#'
#' lambdaSeq <- seq(1, 2000, length.out = 100)
#' tgccFit <- tgCC( data = data, lambdaSeq = lambdaSeq)
#' makeDendrogram(tgccFit, label)
#'
#' # approximate MST by FAMST
#' tgccFit2 <- tgCC(data = data, lambdaSeq = lambdaSeq, mst = "famst")
#' table(clusterLabel(tgccFit2, numClusters = 3), label)
#'
#' @export
tgCC <- function(
  data,
  lambdaSeq,
  bandwidth = NULL,
  useNorm = TRUE,
  depthThresh = 10,
  probThresh = 0.1,
  isNaive = FALSE,
  mst = c("exact", "famst"),
  famstControl = list(),
  keepTheta = TRUE) {

  mst <- match.arg(mst)
  if (!is.list(famstControl))
    stop("`famstControl` must be a list.")

  # Initialize the parameters
  init <- initParams(data, bandwidth, useNorm, isNaive, mst, famstControl)
  params <- init$params

  # Extract parameters
  nodeTypes <- params$Types
  partitionSizes <- params$Partitions
  edgeWeights <- params$EdgeWeights
  nodeWeights <- params$NodeWeights
  vertices <- params$Vertices
  parents <- params$Parents
  children <- params$Children

  # Set a threshold for edge weights of outliers
  outlierIndex <- (partitionSizes < depthThresh)
  edgeWeights[outlierIndex] <-
    pmax(edgeWeights[outlierIndex], stats::quantile(edgeWeights[outlierIndex], probThresh))

  # Optimization using dynamic programming
  coordinates <- pointerList <- list()
  updatedData <- data
  totalFitTime <- 0

  for (lambda in lambdaSeq) {
    numSamples <- nrow(updatedData)
    numFeatures <- ncol(updatedData)
    thetaMatrix <- matrix(0, nrow = numSamples, ncol = numFeatures)
    startTime <- proc.time()

    # Dynamic Programming
    for (featureIdx in seq_len(numFeatures)) {
      thetaMatrix[, featureIdx] <-
        computeTheta(
          updatedData[, featureIdx],
          lambda,
          vertices,
          nodeTypes,
          parents,
          nodeWeights,
          edgeWeights,
          children
        )
    }

    # Update parameters after clustering step
    updatedParams <-
      updatenew(
        thetaMatrix,
        updatedData,
        vertices,
        nodeTypes,
        parents,
        nodeWeights,
        edgeWeights,
        children
      )
    endTime <- proc.time()
    totalFitTime <- totalFitTime + (endTime - startTime)[3]

    # Update data and other parameters
    updatedData  <- updatedParams$newInput
    nodeTypes <- updatedParams$newTypes
    edgeWeights <- updatedParams$newEdgeWeights
    nodeWeights <- updatedParams$newNodeWeights
    vertices  <- updatedParams$newVertices
    children  <- updatedParams$newChildrenList
    parents  <- updatedParams$newParents

    coordinates[[length(coordinates) + 1]] <- updatedData
    pointerList[[length(pointerList) + 1]] <- updatedParams$pointer

    # Break condition
    if (length(vertices) == 1)
      break
  }

  # Obtain the matrices corresponding to different lambda values
  thetaList <- NULL
  if (keepTheta) {
    thetaList <- list()
    sampleIndices <- seq_len(nrow(data))
    for (i in seq_along(coordinates)) {
      sampleIndices <- pointerList[[i]][sampleIndices] + 1
      thetaList[[i]] <- coordinates[[i]][sampleIndices,]
    }
  }

  list(
    data = data,
    theta = thetaList,
    coordinates = coordinates,
    pointer = pointerList,
    lambdaSeq = lambdaSeq,
    tgccTime = totalFitTime,
    initTime = init$initTime,
    mstTime = init$mstTime,
    mst = mst
  )
}
