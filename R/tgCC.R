#' Tree-Guided Convex Clustering
#'
#' @param data matrix with samples as rows and features as columns.
#' @param lambdaSeq sequence of tuning parameters.
#' @param bandwidth bandwidth for the Gaussian kernel between rows.
#' @param useNorm logical, whether to use normalized distance. Default is TRUE.
#' @param depthThresh depth threshold target. Default is 10.
#' @param probThresh probability threshold for edge weight. Default is 0.1.
#' @return A list containing the clustered results and timings.
#' @export
tgCC <- function(data, lambdaSeq, bandwidth = NULL, useNorm = TRUE, depthThresh = 10, probThresh = 0.1) {

  # Initialize the parameters
  init <- initParams(data, bandwidth, useNorm)
  mstTime <- init$mstTime
  initTime <- init$initTime
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
  edgeWeights[outlierIndex] <- pmax(edgeWeights[outlierIndex], stats::quantile(edgeWeights[outlierIndex], probThresh))

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
      thetaMatrix[, featureIdx] <- computeTheta(updatedData[, featureIdx], lambda, vertices, nodeTypes, parents, nodeWeights, edgeWeights, children)
    }

    # Update parameters after clustering step
    updatedParams <- updatenew(thetaMatrix, updatedData, vertices, nodeTypes, parents, nodeWeights, edgeWeights, children)
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
    if (length(vertices) == 1) break
  }

  # Obtain the matrices corresponding to different lambda values
  thetaList <- list()
  sampleIndices <- seq_len(nrow(data))
  for (i in seq_along(coordinates)) {
    sampleIndices <- pointerList[[i]][sampleIndices] + 1
    thetaList[[i]] <- coordinates[[i]][sampleIndices, ]
  }

  list(
    data = data,
    theta = thetaList,
    coordinates = coordinates,
    pointer = pointerList,
    lambdaSeq = lambdaSeq,
    tgccTime = totalFitTime,
    initTime = initTime,
    mstTime = mstTime
  )
}
