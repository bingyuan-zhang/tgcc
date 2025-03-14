## find the loss function
findloss <- function(data, theta, lam, parents, nodeWeights) {
  penalty <- 0
  n <- length(parents)
  for (node in 1:n) {
    penalty = penalty + sum(abs(theta[node, ] - theta[parents[node] + 1, ])) * nodeWeights[node]
  }
  penalty = penalty * lam
  loss = 0.5 * norm(data - theta, type = "f") ^ 2 + penalty
  loss
}

tgccLoss <- function(
    data,
    lamseq,
    bandwidth = NULL,
    useNorm = TRUE,
    depthThresh = 10,
    probThresh = 0.1) {

  # inputs: data (n*p) and lamseq (K).
  input <- data
  k <- length(lamseq)

  # parameters for tree-guided L_1 convex clustering
  init <- tgcc:::initParams(data, bandwidth, useNorm, isNaive = FALSE)
  params <- init$params
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
    pmax(edgeWeights[outlierIndex],
      stats::quantile(edgeWeights[outlierIndex], probThresh))
  params$EdgeWeights <- edgeWeights

  # Naive method Losses
  Losstrue = rep(0, k)
  Losstgcc = rep(0, k)
  n = nrow(input)
  p = ncol(input)
  theta = matrix(0, nrow = n, ncol = p)
  Theta = list()
  Label = list()

  for (i in 1:k) {

    # Compute Loss True
    for (j in 1:p) {
      theta[, j] <- tgcc:::computeTheta(
        input[, j],
        lamseq[i],
        vertices,
        nodeTypes,
        parents,
        nodeWeights,
        edgeWeights,
        children)
    }

    par <- tgcc:::updatenew(
      theta,
      input,
      vertices,
      nodeTypes,
      parents,
      nodeWeights,
      edgeWeights,
      children)

    Label[[i]] <- par$pointer + 1
    Theta[[i]] <- theta
    Losstrue[i] <- findloss(
      data,
      theta,
      lamseq[i],
      parents,
      edgeWeights)
  }

  # TGCC method losses
  idx = 1:nrow(data) - 1
  clustnum = rep(0, k)
  Theta2 = list()
  Label2 = list()

  updatedData <- input
  for (i in 1:k) {
    n = nrow(updatedData)
    p = ncol(updatedData)
    clustnum[i] = n
    theta = matrix(0, nrow = n, ncol = p)

    for (j in 1:p) {
      theta[, j] <- tgcc:::computeTheta(
        updatedData[, j],
        lamseq[i],
        vertices,
        nodeTypes,
        parents,
        nodeWeights,
        edgeWeights,
        children)
    }
    Theta2[[i]] = theta[idx + 1, ]
    Losstgcc[i] = findloss(
      data,
      Theta2[[i]],
      lamseq[i],
      params$Parents,
      params$EdgeWeights)

    # update parameters
    updatedParams <-
      tgcc:::updatenew(
        theta,
        updatedData,
        vertices,
        nodeTypes,
        parents,
        nodeWeights,
        edgeWeights,
        children
      )
    updatedData  <- updatedParams$newInput
    nodeTypes <- updatedParams$newTypes
    edgeWeights <- updatedParams$newEdgeWeights
    nodeWeights <- updatedParams$newNodeWeights
    vertices  <- updatedParams$newVertices
    children  <- updatedParams$newChildrenList
    parents  <- updatedParams$newParents

    idx = updatedParams$pointer[idx + 1]
    Label2[[i]] = idx + 1
    if (length(vertices) == 1) break
  }
  # TGCC stops when all clusters are merged,
  # we change the lambda value to assign the loss of tgcc.
  if (length(vertices) == 1) {
    if (i != k) {
      for (j in (i + 1):k) {
        clustnum[j] = 1
        Losstgcc[j] = findloss(
          data,
          updatedData[idx + 1, ],
          lamseq[j],
          params$Parents,
          params$EdgeWeights)
      }
    }
  }

  list(
    ltrue = Losstrue,
    ltgcc = Losstgcc,
    clustnum = clustnum,
    Thetatrue = Theta,
    Thetatgcc = Theta2,
    Labeltrue = Label,
    Labeltgcc = Label2
  )
}
