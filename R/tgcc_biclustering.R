#' Tree-Guided Convex Bi-Clustering
#'
#' @description
#' Perform the tree-guided convex biclustering.
#'
#' @param data matrix with samples as rows and features as columns.
#' @param lambdaSeq tuning parameters sequence for samples.
#' @param gammaSeq tuning parameters sequence for features.
#' @param bandwidth_s bandwidth_s for the Gaussian kernel between samples.
#' @param bandwidth_f bandwidth_f for the Gaussian kernel between features.
#' @param useNorm logical, whether to use normalized distance. Default is `TRUE`.
#' @param depthThresh_s depth threshold target for the tree structure of samples. Default is `10`..
#' @param depthThresh_f depth threshold target for the tree structure of features. Default is `10`..
#' @param probThresh probability threshold for edge weight. Default is `0.1`.
#' @param maxIter maximum iteration of the Dykstra algorithm.
#' @param threshold stop condition of the Dykstra algorithm.
#' @param isNaive use naive method in searching MST. Default is `FALSE`.
#' @export
#'
#'
#' @export
biTGCC <- function (
  data,
  lambdaSeq,
  gammaSeq,
  bandwidth_s = NULL,
  bandwidth_f = NULL,
  useNorm = TRUE,
  depthThresh_s = 10,
  depthThresh_f = 10,
  probThresh = 0.1,
  maxIter = 500,
  threshold = 1e-5,
  isNaive = FALSE) {

  # initialize parameters
  init <- initParamsBicc(data, bandwidth_s, bandwidth_f, useNorm, isNaive)
  paramsSamples <- init$paramsSamples
  paramsFeatures <- init$paramsFeatures

  nodeTypes_s <- paramsSamples$Types
  edgeWeights_s <- paramsSamples$EdgeWeights
  nodeWeights_s <- paramsSamples$NodeWeights
  partitionSizes_s <- paramsSamples$Partitions
  vertices_s <- paramsSamples$Vertices
  parents_s <- paramsSamples$Parents
  children_s <- paramsSamples$Children

  nodeTypes_f <- paramsFeatures$Types
  edgeWeights_f <- paramsFeatures$EdgeWeights
  nodeWeights_f <- paramsFeatures$NodeWeights
  partitionSizes_f <- paramsFeatures$Partitions
  vertices_f <- paramsFeatures$Vertices
  parents_f <- paramsFeatures$Parents
  children_f <- paramsFeatures$Children

  # threshold for outliers
  ind_s <- (partitionSizes_s < depthThresh_s)
  ind_f <- (partitionSizes_f < depthThresh_f)
  edgeWeights_s[ind_s] <-
    pmax(edgeWeights_s[ind_s], quantile(edgeWeights_s[ind_s], probThresh))
  edgeWeights_f[ind_f] <-
    pmax(edgeWeights_f[ind_f], quantile(edgeWeights_f[ind_f], probThresh))

  # Optimization using dynamic programming
  coordinates <- pointerList_s <- pointerList_f <- list()
  updatedData <- data

  for(k in seq_along(lambdaSeq)) {

    lambda <- lambdaSeq[k]
    gamma <- gammaSeq[k]
    numSamples <- nrow(updatedData)
    numFeatures <- ncol(updatedData)
    updatedDataBefore <- updatedData

    # Dynamic Programming
    if (numSamples == 1) {
      updatedData[1, ] <- computeTheta (
        updatedData,
        gamma,
        vertices_f,
        nodeTypes_f,
        parents_f,
        nodeWeights_f,
        edgeWeights_f,
        children_f)

    } else if (numFeatures == 1) {
      updatedData[, 1] <- computeTheta (
        updatedData,
        lambda,
        vertices_s,
        nodeTypes_s,
        parents_s,
        nodeWeights_s,
        edgeWeights_s,
        children_s)

    } else if (numFeatures != 1 & numSamples != 1) {
      updatedData <- updateUBC(
        updatedData, updatedDataBefore,
        lambda, gamma,
        maxIter, threshold,
        vertices_s, nodeTypes_s, parents_s,
        nodeWeights_s, edgeWeights_s, children_s,
        vertices_f, nodeTypes_f, parents_f,
        nodeWeights_f, edgeWeights_f, children_f)
    }

    # clustering step
    updatedParams_s <-
      updatenew(
        updatedData,
        updatedDataBefore,
        vertices_s,
        nodeTypes_s,
        parents_s,
        nodeWeights_s,
        edgeWeights_s,
        children_s
      )

    updatedData_s <- updatedParams_s$newInput
    nodeTypes_s <- updatedParams_s$newTypes
    edgeWeights_s <- updatedParams_s$newEdgeWeights
    nodeWeights_s <- updatedParams_s$newNodeWeights
    vertices_s <- updatedParams_s$newVertices
    children_s <- updatedParams_s$newChildrenList
    parents_s <- updatedParams_s$newParents

    # break condition
    if (length(vertices_s) == 1) {
      updatedData <- updatedData_s
      coordinates[[k]] <- updatedData_s
      pointerList_s[[k]] <- updatedParams_s$pointer
      pointerList_f[[k]] <- updatedParams_s$pointer
      break
    }

    updatedParams_f <-
      updatenew(
        t(updatedData),
        t(updatedData_s),
        vertices_f,
        nodeTypes_f,
        parents_f,
        nodeWeights_f,
        edgeWeights_f,
        children_f
      )
    updatedData_f <- t(updatedParams_f$newInput)
    nodeTypes_f <- updatedParams_f$newTypes
    edgeWeights_f <- updatedParams_f$newEdgeWeights
    nodeWeights_f <- updatedParams_f$newNodeWeights
    vertices_f <- updatedParams_f$newVertices
    children_f <- updatedParams_f$newChildrenList
    parents_f <- updatedParams_f$newParents

    # break condition
    if (length(vertices_f) == 1) {
      updatedData <- updatedData_f
      coordinates[[k]] <- updatedData_f
      pointerList_s[[k]] <- updatedParams_s$pointer
      pointerList_f[[k]] <- updatedParams_f$pointer
      break
    }

    updatedData <- updatedData_f
    coordinates[[k]] <- updatedData_f
    pointerList_s[[k]] <- updatedParams_s$pointer
    pointerList_f[[k]] <- updatedParams_f$pointer
  }

  # Obtain the matrices corresponding to different lambda values
  thetaList <- list()
  theta <- matrix(0, nrow(data), ncol(data))
  ind_s <- seq_len(nrow(data))
  ind_f <- seq_len(ncol(data))
  for (i in seq_along(coordinates)) {
    ind_s <- pointerList_s[[i]][ind_s] + 1
    ind_f <- pointerList_f[[i]][ind_f] + 1
    theta <- coordinates[[i]][ind_s, ind_f]
    thetaList[[i]] <- theta
  }

  return(
    list(
      data = data,
      theta = thetaList,
      coordinates = coordinates,
      pointer = pointerList_s,
      pointer_f = pointerList_f,
      lambdaSeq = lambdaSeq,
      gammaSeq = gammaSeq
    )
  )
}








