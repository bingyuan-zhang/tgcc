#' Sparse Tree-Guided Convex Clustering
#'
#' @description
#' Perform the sparse tree-guided convex clustering.
#'
#' @param data matrix with samples as rows and features as columns.
#' @param lambdaSeq tuning parameters sequence for samples.
#' @param gammaSeq tuning parameters sequence for features.
#' @param bandwidth bandwidth for the Gaussian kernel between rows. Default is heuristic assigned.
#' @param useNorm logical, whether to use normalized distance. Default is `TRUE`.
#' @param depthThresh depth threshold target. Default is `10`.
#' @param probThresh probability threshold for edge weight. Default is `0.1`.
#' @param maxIter maximum iteration for the Dykstra algorithm.
#' @param threshold stop condition for the Dykstra algorithm.
#' @param isNaive use naive method in searching MST. Default is `FALSE`.
#' @return A list containing the clustered results.
#'
#' @export
spTGCC <- function(
  data,
  lambdaSeq,
  gammaSeq,
  bandwidth = NULL,
  useNorm = TRUE,
  depthThresh = 10,
  probThresh = 0.1,
  maxIter = 500,
  threshold = 1e-5,
  isNaive = FALSE) {

  # Initialize the parameters
  init <- initParams(data, bandwidth, useNorm, isNaive)
  params <- init$params

  nodeTypes <- params$Types
  edgeWeights <- params$EdgeWeights
  nodeWeights <- params$NodeWeights
  partitionSizes <- params$Partitions
  vertices  <- params$Vertices
  parents  <- params$Parents
  children  <- params$Children

  # Set a threshold for edge weights of outliers
  ind <- (partitionSizes < depthThresh)
  edgeWeights[ind] <-
    pmax(edgeWeights[ind], quantile(edgeWeights[ind], probThresh))

  # Optimization using dynamic programming
  coordinates <- pointerList <- leftFeatureList <- list()
  updatedData <- data
  leftFeature <- 1:ncol(data)

  for(k in seq_along(lambdaSeq)) {

    numSamples <- nrow(updatedData)
    numFeatures <- ncol(updatedData)
    lambda <- lambdaSeq[k]
    gamma <- gammaSeq[k]
    updatedDataBefore <- updatedData

    # Dynamic Programming
    if (numSamples == 1) {
      updatedData[, 1] <- computeTheta(
        updatedData,
        gamma,
        vertices,
        nodeTypes,
        parents,
        nodeWeights,
        edgeWeights,
        children)
    } else if (numFeatures != 1 & numSamples != 1) {
      updatedData <- updateUSC(
        updatedData,
        updatedDataBefore,
        lambda,
        gamma,
        maxIter,
        threshold,
        vertices,
        nodeTypes,
        parents,
        nodeWeights,
        edgeWeights,
        children
      )
    }

    # Select non-zero features
    selectedFeature <- which(colSums(updatedData) != 0)
    updatedData <- updatedData[ ,selectedFeature]
    updatedDataBefore <- updatedDataBefore[,selectedFeature]

    if(!is.matrix(updatedData)) {
      updatedData <- matrix(updatedData, length(updatedData), ncol=1)
      updatedDataBefore <- matrix(updatedDataBefore, length(updatedDataBefore), ncol=1)
    }

    # Clustering step
    updatedParams <- updatenew(updatedData, updatedDataBefore,
                        vertices, nodeTypes, parents,
                        nodeWeights, edgeWeights, children)

    updatedData <- updatedParams$newInput
    nodeTypes <- updatedParams$newTypes
    edgeWeights <- updatedParams$newEdgeWeights
    nodeWeights <- updatedParams$newNodeWeights
    vertices  <- updatedParams$newVertices
    children  <- updatedParams$newChildrenList
    parents  <- updatedParams$newParents

    leftFeature <- leftFeature[selectedFeature]
    coordinates[[length(coordinates) + 1]] <- updatedData
    pointerList[[length(pointerList) + 1]] <- updatedParams$pointer
    leftFeatureList[[length(leftFeatureList) + 1]] <- leftFeature

    # ---- break condition ----
    if (length(vertices) == 1) break
  }

  # Obtain the matrices corresponding to different lambda values
  thetaList <- list()
  theta <- matrix(0, nrow(data), ncol(data))
  sampleIndices <- seq_len(nrow(data))
  for (i in seq_along(coordinates)) {
    sampleIndices <- pointerList[[i]][sampleIndices] + 1
    featureIndices <- leftFeatureList[[i]]
    theta[, featureIndices] <- coordinates[[i]][sampleIndices, ]
    theta[, -featureIndices] <- NA
    thetaList[[i]] <- theta
  }

  list(
    data = data,
    theta = thetaList,
    coordinates = coordinates,
    pointer = pointerList,
    leftFeature = leftFeatureList,
    lambdaSeq = lambdaSeq,
    gammaSeq = gammaSeq
  )
}

