# DEBUG code
updateRootTheta <- function(nodeValue, nodeWeight, boundary) {
  # 复制 boundary
  b <- boundary

  # 计算 theta 的初始值
  theta <- (nodeWeight * nodeValue - b$intercepts[1]) / (b$slopes[1] + nodeWeight)

  # 迭代更新 theta
  while (b$boundaryPoints[1] < theta) {
    b$intercepts <- b$intercepts[-1]
    b$slopes <- b$slopes[-1]
    b$boundaryPoints <- b$boundaryPoints[-1]
    theta <- (nodeWeight * nodeValue - b$intercepts[1]) / (b$slopes[1] + nodeWeight)
    if (length(b$slopes) == 1) break
  }

  return(list(boundary = b, theta = theta))
}

updateBoundariesForLeaf <- function(nodeValue, nodeWeight, lambda, boundary) {
  b <- boundary

  lb <- nodeValue - lambda / nodeWeight
  ub <- nodeValue + lambda / nodeWeight

  b$boundaryPoints <- c(lb, b$boundaryPoints, ub)
  b$intercepts <- c(-lambda, b$intercepts, -nodeWeight * nodeValue, lambda)
  b$slopes <- c(0, b$slopes, nodeWeight, 0)

  return(list(boundary = b, lb = lb, ub = ub))
}

updateBoundaries <-
  function(
    nodeValue,
    nodeWeight,
    lambda,
    boundary) {
  # 复制 boundary

  b <- boundary
  lb <- (-lambda + nodeWeight * nodeValue - b$intercepts[1]) / (b$slopes[1] + nodeWeight)

  while (b$boundaryPoints[1] < lb) {
    b$intercepts <- b$intercepts[-1]
    b$slopes <- b$slopes[-1]
    b$boundaryPoints <- b$boundaryPoints[-1]
    lb <- (-lambda + nodeWeight * nodeValue - b$intercepts[1]) / (b$slopes[1] + nodeWeight)

    if (length(b$slopes) == 1) break
  }

  # 计算初始上界
  if (length(b$slopes) == 1) {
    ub <- (lambda + nodeWeight * nodeValue - b$intercepts[1]) / (b$slopes[1] + nodeWeight)
  } else {
    ub <- (lambda + nodeWeight * nodeValue - b$intercepts[length(b$intercepts)]) / (b$slopes[length(b$slopes)] + nodeWeight)

    # 通过迭代更新上界
    while (b$boundaryPoints[length(b$boundaryPoints)] > ub) {
      b$intercepts <- b$intercepts[-length(b$intercepts)]
      b$slopes <- b$slopes[-length(b$slopes)]
      b$boundaryPoints <- b$boundaryPoints[-length(b$boundaryPoints)]
      ub <- (lambda + nodeWeight * nodeValue - b$intercepts[length(b$intercepts)]) / (b$slopes[length(b$slopes)] + nodeWeight)
      if (length(b$slopes) == 1) break
    }
  }

  # 更新 boundaryPoints
  b$boundaryPoints <- c(lb, b$boundaryPoints, ub)

  # 更新 intercepts
  for (i in 1:length(b$intercept)){
    b$intercepts[i] <- b$intercepts[i] - nodeWeight * nodeValue
  }
  b$intercepts <- c(-lambda, b$intercepts, lambda)

  # 更新 slopes
  for (i in 1:length(b$slopes)){
    b$slopes[i] <- b$slopes[i] + nodeWeight
  }
  b$slopes <- c(0, b$slopes, 0)

  return(list(boundary = b, lb = lb, ub = ub))
}

mergeQueues <- function(children, Boundaries) {

  numChildren <- length(children)

  bpLengths <- rep(0, numChildren)
  mergedBoundaryPoints <- c()

  for (i in 1:numChildren) {
    child <- children[i]
    bpLengths[i] <- length(Boundaries[[child]]$boundaryPoints)
    mergedBoundaryPoints <-
      c(mergedBoundaryPoints,
        Boundaries[[child]]$boundaryPoints)
  }

  totalLength <- length(mergedBoundaryPoints)

  # 获得 boundaryPoints 的排序顺序
  bpSort <- order(mergedBoundaryPoints)
  ranks <- integer(totalLength)
  ranks[bpSort] <- seq_along(bpSort)

  # 排序 mergedBoundaryPoints
  mergedBoundaryPoints <- sort(mergedBoundaryPoints)

  # 合并队列
  interceptsToSum <- matrix(0, nrow = numChildren, ncol = totalLength + 1)
  slopesToSum <- matrix(0, nrow = numChildren, ncol = totalLength + 1)

  start1 <- 1
  for (i in 1:numChildren) {
    child <- children[i]
    start2 <- 1
    for (j in 1:bpLengths[i]) {
      # 填充 intercepts
      intercept_idx <- start2:ranks[start1 + j - 1]
      interceptsToSum[i, intercept_idx] <- Boundaries[[child]]$intercepts[j]

      # 填充 slopes
      slopesToSum[i, intercept_idx] <- Boundaries[[child]]$slopes[j]

      start2 <- ranks[start1 + j - 1] + 1
    }

    # 填充剩余部分
    interceptsToSum[i, start2:(totalLength + 1)] <- Boundaries[[child]]$intercepts[bpLengths[i] + 1]
    slopesToSum[i, start2:(totalLength + 1)] <- Boundaries[[child]]$slopes[bpLengths[i] + 1]

    start1 <- start1 + bpLengths[i]
  }

  newIntercepts <- rep(0.0, totalLength + 1)
  newSlopes <- rep(0.0, totalLength + 1)

  # 更新 intercepts 和 slopes
  for (j in 1:(totalLength + 1)) {
    for (i in 1:numChildren) {
      newIntercepts[j] <- newIntercepts[j] + interceptsToSum[i, j]
      newSlopes[j] <- newSlopes[j] + slopesToSum[i, j]
    }
  }

  # 创建新的 Boundary
  b <- list(
    boundaryPoints = mergedBoundaryPoints,
    intercepts = newIntercepts,
    slopes = newSlopes
  )

  return(b)
}

computeTheta <-
  function(
    nodeValues,
    lambda,
    vertices,
    nodeTypes,
    parents,
    nodeWeights,
    edgeWeights,
    childrenList) {

  numNodes <- length(nodeValues)

  # 初始化 nodeBoundaries、lowerBounds、upperBounds
  nodeBoundaries <- vector("list", numNodes)
  lowerBounds <- numeric(numNodes)
  upperBounds <- numeric(numNodes)

  # 初始化 theta
  theta <- rep(0, numNodes)

  # 前向遍历
  for (i in 1:numNodes) {

    vertex <- vertices[i] + 1
    nodeType <- nodeTypes[vertex]  # R索引从1开始
    nodeWeight <- nodeWeights[vertex]
    edgeWeight <- edgeWeights[vertex]
    nodeValue <- nodeValues[vertex]

    if (nodeType == 0) {  # leaves
      res <- updateBoundariesForLeaf(
        nodeValue = nodeValue,
        nodeWeight = nodeWeight,
        lambda = edgeWeight * lambda,
        boundary = nodeBoundaries[[vertex]]
      )
      nodeBoundaries[[vertex]] <- res$boundary
      lowerBounds[i] <- res$lb
      upperBounds[i] <- res$ub

    } else if (nodeType == 1) {  # single child
      child <- childrenList[[vertex]][1] + 1
      nodeBoundaries[[vertex]] <- nodeBoundaries[[child]]

      # if (i == 6816) {
      #
      #   # debugList2 <<-  list(
      #   #     nodeValue = nodeValue,
      #   #     nodeWeight = nodeWeight,
      #   #     lambda = edgeWeight * lambda,
      #   #     boundary = nodeBoundaries[[vertex]])
      #   #
      #   # saveRDS(
      #   #   debugList2,
      #   #   file = "./simulation/runtime_by_samples/debug2.rds")
      #   browser()
      # }

      res <- updateBoundaries(
        nodeValue = nodeValue,
        nodeWeight = nodeWeight,
        lambda = edgeWeight * lambda,
        boundary = nodeBoundaries[[vertex]]
      )
      nodeBoundaries[[vertex]] <- res$boundary
      lowerBounds[i] <- res$lb
      upperBounds[i] <- res$ub

    } else if (nodeType == 2) {  # multiple children
      children <- childrenList[[vertex]] + 1
      updatedBoundary <- mergeQueues(children, nodeBoundaries)

      res <- updateBoundaries(
        nodeValue = nodeValue,
        nodeWeight = nodeWeight,
        lambda = edgeWeight * lambda,
        boundary = updatedBoundary
      )
      nodeBoundaries[[vertex]] <- res$boundary
      lowerBounds[i] <- res$lb
      upperBounds[i] <- res$ub

    } else {  # root
      if (nodeType == 3) {  # root with single child
        child <- childrenList[[vertex]][1] + 1
        res <- updateRootTheta(
          nodeValue = nodeValue,
          nodeWeight = nodeWeight,
          boundary = nodeBoundaries[[child]]
        )
        nodeBoundaries[[child]] <- res$boundary
        theta[vertex] <- res$theta

      } else {  # root with multiple children
        children <- childrenList[[vertex]] + 1
        updatedBoundary <- mergeQueues(children, nodeBoundaries)

        res <- updateRootTheta(
          nodeValue = nodeValue,
          nodeWeight = nodeWeight,
          boundary = updatedBoundary
        )
        nodeBoundaries[[vertex]] <- res$boundary
        theta[vertex] <- res$theta
      }
    }
  }

  # 反向遍历
  for (i in (numNodes - 1):1) {  # R中的索引从1开始
    vertex <- vertices[i] + 1
    parent <- parents[vertex] + 1

    if (theta[parent] < lowerBounds[i]) {
      theta[vertex] <- lowerBounds[i]
    } else if (theta[parent] > upperBounds[i]) {
      theta[vertex] <- upperBounds[i]
    } else {
      theta[vertex] <- theta[parent]
    }
  }

  return(theta)
}



################################################################################
# 示例数据
# set.seed(1)
# dl <- tgcc:::make.mixgaussian(n = 1e3)
# data <- dl$data
# lamFrom <- 1
# lamTo <- norm(dl$data, type="2")^2
# lenNum <- 10
# lambdaSeq <- seq(lamFrom, lamTo, length.out = lenNum)
#
# bandwidth = NULL
# useNorm = TRUE
# depthThresh = 10
# probThresh = 0.1
# isNaive = FALSE
#
# # Initialize the parameters
# init <- tgcc:::initParams(data, bandwidth, useNorm, isNaive)
# params <- init$params
#
# # Extract parameters
# nodeTypes <- params$Types
# partitionSizes <- params$Partitions
# edgeWeights <- params$EdgeWeights
# nodeWeights <- params$NodeWeights
# vertices <- params$Vertices
# parents <- params$Parents
# childrenList <- params$Children
#
# # Set a threshold for edge weights of outliers
# outlierIndex <- (partitionSizes < depthThresh)
# edgeWeights[outlierIndex] <-
#   pmax(edgeWeights[outlierIndex], stats::quantile(edgeWeights[outlierIndex], probThresh))
#
# # Optimization using dynamic programming
# coordinates <- pointerList <- list()
# updatedData <- data
#
# lambda <- lambdaSeq[1]
# numSamples <- nrow(updatedData)
# numFeatures <- ncol(updatedData)
# thetaMatrix <- matrix(0, nrow = numSamples, ncol = numFeatures)
#
# # 调用 computeTheta
# theta <- computeTheta(
#   updatedData[, 1],
#   lambda,
#   vertices,
#   nodeTypes,
#   parents,
#   nodeWeights,
#   edgeWeights,
#   childrenList
# )
#
# theta0 <- tgcc:::computeTheta(
#   updatedData[, 1],
#   lambda,
#   vertices,
#   nodeTypes,
#   parents,
#   nodeWeights,
#   edgeWeights,
#   childrenList
# )
#
# all(theta == theta0)
#
