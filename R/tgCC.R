#' tree-guided convex clustering
#'
#' @param data matrix with samples as rows and features as columns.
#' @param lamseq tuning parameters sequence.
#' @param bw bandwidth for the Gaussian kernel between rows.
#' @param norm to use normalized distance if true.
#' @param dth threshold target with in depth_thres.
#' @param prob threshold of edge weight.
#' @export
tgCC <- function(
    data,
    lamseq,
    bw = NULL,
    norm = TRUE,
    dth = 10,
    prob = 0.1) {

  # initialize the parameters.
  init <- init_params(data, bw, norm)
  mst_time <- init$mst_time
  init_time <- init$init_time
  params <- init$params

  # parameters for tree-guided L-1 convex clustering.
  Tp <- params$Types
  Pt <- params$Partitions
  WE <- params$EdgeWeights
  WV <- params$NodeWeights
  V  <- params$Vertices
  P  <- params$Parents
  C  <- params$Children

  # set a threshold for edge weights of outliers.
  ind <- (Pt < dth)
  WE[ind] <- pmax(WE[ind], quantile(WE[ind], prob))

  # optimization using dynamic programming.
  Coord <- Pointer <- list()
  U <- data
  K <- length(lamseq)
  fit_time <- 0

  for(k in 1:K) {
    n <- nrow(U)
    p <- ncol(U)
    lam <- lamseq[k]
    Theta <- matrix(0, nrow = n, ncol = p)

    t_start = proc.time()

    # Dynamic Programming
    for (j in 1:p) {
      Theta[, j] = computeTheta (U[, j], lam, V, Tp, P, WV, WE, C)
    }

    # Cluster step
    params = updatenew(Theta, U, V, Tp, P, WV, WE, C)

    t_end = proc.time()
    fit_time = fit_time + t_end[3] - t_start[3]

    U  <- params$newInput
    Tp <- params$newTypes
    WE <- params$newEdgeWeights
    WV <- params$newNodeWeights
    V  <- params$newVertices
    C  <- params$newChildrenList
    P  <- params$newParents

    Coord[[k]] = U
    Pointer[[k]] = params$pointer

    # break condition
    if (length(V) == 1) break
  }

  # obtain the matrices corresponding to different lambda values.
  Theta <- list()
  id <- 1:nrow(data)
  for(i in seq_along(Coord)){
    id <- Pointer[[i]][id] + 1;
    Theta[[i]] <- Coord[[i]][id,]
  }

  list(
    data = data,
    theta = Theta,
    coord = Coord,
    pointer = Pointer,
    lamseq = lamseq,
    tgcc.time = fit_time,
    init.time = init_time,
    mst.time = mst_time
  )

}
