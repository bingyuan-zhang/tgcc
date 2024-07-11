#' sparse tree-guided convex clustering
#'
#' @param data matrix with samples as rows and features as columns.
#' @param lamseq tuning parameters sequence for rows.
#' @param bw bandwidth for the Gaussian kernel between rows.
#' @param norm to use normalized distance if true.
#' @param dth threshold target of rows with in dth.
#' @param prob threshold of edge weight.
#' @param max_iter maximum iteration of Dykstra algorithm.
#' @param thres the stop condition of Dykstra algorithm.
#' @export
spTGCC <- function(
  data,
  lamseq,
  gamseq,
  bw = NULL,
  norm = TRUE,
  dth = 50,
  prob = 0.1,
  max_iter = 500,
  thres = 1e-5) {

  # find the initialization parameters
  init <- init_params(data, bw, norm)
  init_params <- init$params

  Tp <- init_params$Type
  WE <- init_params$WeightE
  WV <- init_params$WeightV
  Pt <- init_params$Partition
  V  <- init_params$Vertice
  P  <- init_params$Parent
  C  <- init_params$Children

  # threshold for outliers
  ind <- (Pt < dth)
  WE[ind] <- pmax(WE[ind], quantile(WE[ind], prob))

  # Optimization
  # Coord: Updated Parameter matrix
  # Pointer (s: sample, f: feature) : pointers that record the merge events
  # Fleft: the nonzero features of original column index with each lambda
  Coord <- Pointer <- Fleft <- list()
  U <- data
  fleft <- 1:ncol(data)

  for(k in seq_along(lamseq)) {
    n <- nrow(U)
    p <- ncol(U)
    lam <- lamseq[k]
    gam <- gamseq[k]
    U0 <- U

    # Dynamic Programming
    if (n == 1) {
      U[,1] <- dp (U, gam, V, Tp, P, WV, WE, C)
    } else if (p != 1 & n != 1) {
      U <- updateUSC(U, U0, lam, gam, max_iter, thres, V, Tp, P, WV, WE, C)
    }

    # Record non zero features
    select.features <- which(colSums(U) != 0)
    U <- U[,select.features]
    U0 <- U0[,select.features]

    # Clustering step
    params <- updatenew(U, U0, V, Tp, Pa, WV, WE, C)
    U <- params$newinput
    V <- params$newV
    P <- params$newP
    C <- params$newC
    Tp <- params$newT
    WE <- params$newWE
    WV <- params$newWV

    fleft <- fleft[select.features]
    Coord[[k]] <- U
    Pointer[[k]] <- params$pointer
    Fleft[[k]] <- fleft

    # ---- break condition ----
    if (length(V) == 1) break

  }

  return(list(
    data = data,
    coord = Coord,
    pointer = Pointer,
    fleft = Fleft,
    lamseq = lamseq[1:k],
    gamseq = gamseq[1:k]
  ))
}

#' bi-tree-guided convex clustering
#'
#' @param data matrix with samples as rows and features as columns.
#' @param lamseq tuning parameters sequence for rows.
#' @param gamseq tuning parameters sequence for columns.
#' @param bw_s bandwidth_s for the Gaussian kernel between rows.
#' @param bw_f bandwidth_f for the Gaussian kernel between columns.
#' @param norm to use normalized distance if true.
#' @param dth_s threshold target of rows with in dth_s.
#' @param dth_f threshold target of columns with in dth_f.
#' @param prob threshold of edge weight.
#' @param max_iter maximum iteration of Dykstra algorithm.
#' @param thres the stop condition of Dykstra algorithm.
#' @export
biTGCC <- function (
  data,
  lamseq,
  gamseq,
  bw_s = NULL,
  bw_f = NULL,
  norm = TRUE,
  dth_s = 10,
  dth_f = 10,
  prob = 0.1,
  max_iter = 500,
  thres = 1e-5) {

  # initialize parameters
  init <- init_params_bicc(data, bw_s, bw_f, norm)
  init_params_s <- init$params_s
  init_params_f <- init$params_s

  Tp_s <- init_params_s$Type
  WE_s <- init_params_s$WeightE
  WV_s <- init_params_s$WeightV
  Pt_s <- init_params_s$Partition
  V_s <- init_params_s$Vertice
  P_s <- init_params_s$Parent
  C_s <- init_params_s$Children

  Tp_f <- init_f$Type
  WE_f <- init_f$WeightE
  WV_f <- init_f$WeightV
  Pt_f <- init_f$Partition
  V_f <- init_f$Vertice
  P_f <- init_f$Parent
  C_f <- init_f$Children

  # threshold for outliers
  ind_s <- (Pt_s < dth_s)
  ind_f <- (Pt_f < dth_f)
  WE_s[ind_s] <- pmax(WE_s[ind_s], quantile(WE_s[ind_s], prob))
  WE_f[ind_f] <- pmax(WE_f[ind_f], quantile(WE_f[ind_f], prob))

  # Optimization
  # Coord: Updated Parameter matrix
  # Pointer (s: sample, f: feature) : pointers that record the merge events
  Coord <- list()
  Pointer_s <- list()
  Pointer_f <- list()
  U <- data

  for(k in 1:seq_along(lamseq)) {
    lam <- lamseq[k]
    gam <- gamseq[k]
    n <- nrow(U)
    p <- ncol(U)
    U0 <- U

    # dynamic programming
    if (n == 1) {
      U[1,] = dp (U, gam, V_f, Tp_f, P_f, WV_f, WE_f, C_f)
    } else if (p == 1) {
      U[,1] = dp (U, lam, V_s, Tp_s, P_s, WV_s, WE_s, C_s)
    } else if (p != 1 & n != 1) {
      U = updateUBC(U, U0, lam, gam, max_iter, thres,
        V_s, Tp_s, P_s, WV_s, WE_s, C_s,
        V_f, Tp_f, P_f, WV_f, WE_f, C_f)
    }

    # clustering step
    params_s <- updatenew(U, U0, V_s, Tp_s, P_s, WV_s, WE_s, C_s)

    U_s <- params_s$newinput
    V_s <- params_s$newV
    C_s <- params_s$newC
    P_s <- params_s$newP
    WE_s <- params_s$newWE
    WV_s <- params_s$newWV
    Tp_s <- params_s$newT

    # break condition
    if (length(V_s) == 1) {
      U <- U_s
      Coord[[k]] <- U_s
      Pointer_s[[k]] <- params_s$pointer
      Pointer_f[[k]] <- params_s$pointer
      break
    }

    params_f <- updatenew(t(U), t(U_s), V_f, Tp_f, P_f, WV_f, WE_f, C_f)
    U_f <- t(params_f$newinput)
    V_f <- params_f$newV
    P_f <- params_f$newP
    C_f <- params_f$newC
    Tp_f <- params_f$newT
    WE_f <- params_f$newWE
    WV_f <- params_f$newWV

    # break condition
    if (length(V_f) == 1) {
      U <- U_f
      Coord[[k]] <- U_f
      Pointer_s[[k]] <- params_s$pointer
      Pointer_f[[k]] <- params_f$pointer
      break
    }

    U <- U_f
    Coord[[k]] <- U_f
    Pointer_s[[k]] <- params_s$pointer
    Pointer_f[[k]] <- params_f$pointer
  }

  return(
    list(
      data = data,
      coord = Coord,
      pointer_s = Pointer_s,
      pointer_f = Pointer_f,
      lamseq = lamseq[1:k],
      gamseq = gamseq[1:k]
    )
  )
}








