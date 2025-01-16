checkLamTGCC <- function(
    data,
  bandwidth = NULL,
  stopat = -Inf,
  step_size = 1000,
  min_step = 50,
  normalize = TRUE,
  depthThresh = 10,
  probThresh = 0.1,
  maxlam = NULL) {

  lam_vec <- c()
  num_vec <- c()
  tmp_lam <- step_size
  tmp_step <- step_size
  cnt <- 0
  tmp <- TRUE

  if(is.null(maxlam)) maxlam = sum(data^2)/2

  initpars = tgcc:::initParams(
    data, bandwidth, normalize)
  init = initpars$params

  while (tmp == TRUE && tmp_lam < maxlam) {
    cnt <- cnt + 1
    estlabel = dpLam(
      data,
      tmp_lam,
      init,
      depthThresh,
      probThresh)

    tmp_num <- length(unique(estlabel))
    lam_vec <- c(lam_vec, tmp_lam)
    num_vec <- c(num_vec, tmp_num)
    if(tmp_num < stopat) break

    if(tmp_num == 1){
      tmp_step <- tmp_step/2
      tmp_lam <- tmp_lam - tmp_step
      if(tmp_step < min_step){
        tmp = FALSE
        break
      }
    } else {
      if(min(num_vec[1:cnt]) == 1){
        tmp_step <- tmp_step/2
      }
      tmp_lam <- tmp_lam + tmp_step
    }
  }

  list(
    lam_max = tmp_lam + tmp_step,
    lam_vec = lam_vec,
    num_vec = num_vec
  )
}

dpLam <- function(
    data,
  lambda,
  init,
  depthThresh,
  probThresh) {

  nodeTypes <- init$Types
  partitionSizes <- init$Partitions
  edgeWeights <- init$EdgeWeights
  nodeWeights <- init$NodeWeights
  vertices <- init$Vertices
  parents <- init$Parents
  children <- init$Children

  # threshold for outlier
  tstart = proc.time()
  outlierIndex <- (partitionSizes < depthThresh)
  edgeWeights[outlierIndex] <-
    pmax(edgeWeights[outlierIndex], stats::quantile(edgeWeights[outlierIndex], probThresh))

  input = data
  n = nrow(input)
  p = ncol(input)
  thetaMatrix = matrix(0, nrow = n, ncol = p)

  for(j in 1:p){
    thetaMatrix[,j]  = tgcc:::computeTheta(
      input[,j],
      lambda,
      vertices,
      nodeTypes,
      parents,
      nodeWeights,
      edgeWeights,
      children)
  }

  uniquetheta <- unique(thetaMatrix)
  label = rep(0, n)
  for(j in 1:nrow(uniquetheta)) {
    idx = rep(TRUE, n)
    for(p in 1:ncol(uniquetheta)){
      idx = (idx & uniquetheta[j,p] == thetaMatrix[,p])
    }
    label[idx] = j
  }

  return(label)
}


