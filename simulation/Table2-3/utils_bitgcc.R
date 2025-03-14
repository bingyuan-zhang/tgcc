library(cvxbiclustr)
library(tgcc)
#===============================================================================
# Rewrite the cobra function for comparison
#===============================================================================
# for BICC, the convergence condition is sqrt(\|U.new-Y.new\|_F^2)
gkn_weights2 <- function(X,
  phi = 0.5,
  k_row = 5,
  k_col = 5) {
  p <- nrow(X)
  n <- ncol(X)

  ## Construct Gaussian kernel weights
  w_row <- cvxbiclustr:::kernel_weights(t(X), phi / median(dist(t(X)) ^
      2))
  w_col <- cvxbiclustr:::kernel_weights(X, phi / median(dist(X) ^ 2))
  ## Assign weights to k-nearest neighbors
  w_row <- cvxbiclustr:::knn_weights(w_row, k_row, p)
  w_col <- cvxbiclustr:::knn_weights(w_col, k_col, n)

  ## Construct edge-incidence matrices
  E_row <- cvxbiclustr:::create_edge_incidence(w_row, p)
  E_col <- cvxbiclustr:::create_edge_incidence(w_col, n)

  ## Get connectivity information
  nRowComp <- length(cvxbiclustr:::find_clusters(cvxbiclustr:::weights_graph(w = w_row, p))$size)
  nColComp <- length(cvxbiclustr:::find_clusters(cvxbiclustr:::weights_graph(w = w_col, n))$size)

  list(
    w_row = w_row@x,
    w_col = w_col@x,
    E_row = E_row,
    E_col = E_col,
    nRowComp = nRowComp,
    nColComp = nColComp
  )
}

find_cobra_clust <- function(U) {
  U <- round(U, 5)
  clust.lab <- rep(1, nrow(U))
  unique.clust <- unique(U)
  nclust <- nrow(unique.clust)
  for (i in 1:nclust) {
    is.cur.clust <- sapply(1:nrow(U), function(j)
      identical(U[j, ], unique.clust[i, ]))
    clust.lab[is.cur.clust] <- i
  }
  clust.lab
}

find_cobra_clustseq <- function(Useq, k) {
  nseq = length(Useq)

  for (i in 1:nseq) {
    U <- round(Useq[[i]], 5)
    clust.num <- nrow(unique(U))
    if (clust.num <= k)
      break
  }

  clust.lab <- rep(1, nrow(U))
  unique.clust <- unique(U)
  nclust <- nrow(unique.clust)
  for (i in 1:nclust) {
    is.cur.clust <- sapply(1:nrow(U), function(j)
      identical(U[j, ], unique.clust[i, ]))
    clust.lab[is.cur.clust] <- i
  }
  clust.lab
}

#===============================================================================
# Generate checkerboard data
#===============================================================================
generate_bicc_datalist <- function(){
  rep = 10
  n_nums = c(100, 200, 300, 400, 500)
  gdata_list <- list()

  for(n in n_nums)for(i in 1:rep){
    gdata.cur = tgcc:::make.checkerboard(n)
    gdata_list[[paste(n)]][[paste(i)]] = gdata.cur
  }
  gdata_list
}

#===============================================================================
# Perform evaluation for bitgcc
#===============================================================================
evaluate_bitgcc <- function (gdata_list) {
  # parameters
  lamseq <- gamseq <- seq(1, 200, length.out=50)
  thres <- 1e-5
  max_iter <- 1000
  rep <- 10

  # run simulation evaluation
  res.cobra <- res.bitgcc <- list()
  ARI.mat <- AC.mat <- Time.mat <- array(0, dim = c(2,5,rep))
  for(i in 1:5) for(j in 1:rep){
    cat(i,"_", j,"\n")

    gdata.cur = gdata_list[[i]][[j]]
    data = gdata.cur$data
    n = nrow(data)
    p = ncol(data)
    label = gdata.cur$label

    # cobra
    tic <- proc.time()
    invisible(capture.output({
      k <- 5
      wts <- gkn_weights2(t(data), phi=1, k_row=k, k_col=k)
      sol <-
        cobra(t(data),
          wts$E_row,
          wts$E_col,
          wts$w_row,
          wts$w_col,
          gamma = lamseq,
          tol = thres,
          max_iter = max_iter)
    }))
    toc <- proc.time()
    templab <- find_cobra_clustseq(sol$U, 4) # true cluster number 4
    ac1 <- 1-mclust::classError(templab, label)$errorRate
    ari1 <- mclust::adjustedRandIndex(templab, label)
    t1 <- (toc - tic)[3]

    # bitgcc
    tic <- proc.time()
    tgccFit <- tgcc::biTGCC(
      data = data,
      lambdaSeq = lamseq, gammaSeq = gamseq,
      threshold = thres*n*p, # The same condition for biTGCC and cvxbiclustr
      maxIter = max_iter)
    toc <- proc.time()
    templab <- tgcc::clusterLabel(tgccFit, 4)
    ac2 <- 1-mclust::classError(templab, label)$errorRate
    ari2 <- mclust::adjustedRandIndex(templab, label)
    t2 <- (toc - tic)[3]

    # store results
    ARI.mat[1,i,j] <- ari1
    AC.mat[1,i,j] <- ac1
    Time.mat[1,i,j] <- t1

    ARI.mat[2,i,j] <- ari2
    AC.mat[2,i,j] <- ac2
    Time.mat[2,i,j] <- t2
  }

  list(ARI = ARI.mat, AC = AC.mat, Time = Time.mat)
}

