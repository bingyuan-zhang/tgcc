#' @importFrom mlpack emst
#' @importFrom stats median
initParams <- function(data,
  bandwidth = NULL,
  normalize = TRUE,
  isNaive = TRUE,
  mst = "exact",
  famstControl = list(nNeighbors = 15, numRandomEdges = 3, metric = "euclidean")) {

  # Prepare the minimum spanning tree for samples:
  # a three-column edge list (0-based from, 0-based to, edge length).
  tStart <- proc.time()
  if (mst == "exact") {
    resMST <- mlpack::emst(data, naive = isNaive)$output
  } else {
    resFAMST <- do.call(famst, c(list(data = data), famstControl))
    if (resFAMST$numComponents != 1)
      stop("FAMST did not return a spanning tree; try a larger `nNeighbors`.")
    resMST <- cbind(resFAMST$mst$u, resFAMST$mst$v, resFAMST$mst$w)
  }
  tEnd <- proc.time()
  mstTime <- (tEnd - tStart)[3]

  # Normalize distances between samples.
  if (normalize) {
    distSq <- resMST[, 3] ^ 2
    resMST[, 3] <- sqrt(distSq / mean(distSq))
  }

  # Use the median heuristic as bandwidth if not provided.
  if (is.null(bandwidth))
    bandwidth <- stats::median(resMST[, 3] ^ 2)

  # Initialize parameters.
  tStart <- proc.time()
  params <- initPrepare(resMST[, 1], resMST[, 2], resMST[, 3], bandwidth)
  tEnd <- proc.time()
  initTime <- (tEnd - tStart)[3]

  list(params = params,
    mstTime = mstTime,
    initTime = initTime)
}

initParamsBicc <- function(data,
  bandwidthSamples = NULL,
  bandwidthFeatures = NULL,
  normalize = TRUE,
  isNaive = TRUE) {
  # Prepare minimum spanning tree for samples and features.
  resMSTSamples <- mlpack::emst(data, naive = isNaive)$output
  resMSTFeatures <- mlpack::emst(t(data), naive = isNaive)$output

  # Normalize distances.
  if (normalize) {
    distSqSamples <- resMSTSamples[, 3] ^ 2
    resMSTSamples[, 3] <- sqrt(distSqSamples / mean(distSqSamples))
    distSqFeatures <- resMSTFeatures[, 3] ^ 2
    resMSTFeatures[, 3] <- sqrt(distSqFeatures / mean(distSqFeatures))
  }

  # Use median heuristic for bandwidth if not provided.
  if (is.null(bandwidthSamples))
    bandwidthSamples <- stats::median(resMSTSamples[, 3] ^ 2)
  if (is.null(bandwidthFeatures))
    bandwidthFeatures <- stats::median(resMSTFeatures[, 3] ^ 2)

  # Initialize parameters for samples and features.
  paramsSamples <- initPrepare(resMSTSamples[, 1],
    resMSTSamples[, 2],
    resMSTSamples[, 3],
    bandwidthSamples)
  paramsFeatures <- initPrepare(resMSTFeatures[, 1],
    resMSTFeatures[, 2],
    resMSTFeatures[, 3],
    bandwidthFeatures)

  list(paramsSamples = paramsSamples, paramsFeatures = paramsFeatures)
}
