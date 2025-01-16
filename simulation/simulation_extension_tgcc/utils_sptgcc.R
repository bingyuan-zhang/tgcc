library(scvxclustr)
library(tgcc)
#===============================================================================
# Prepare function for comparison
#===============================================================================
gkn_weights_spcc <- function(X, phi=1, k=5) {
  # default setting in scvxclustr
  n <- nrow(X)
  p <- ncol(X)
  w <- scvxclustr::dist_weight( t(X) / sqrt(p), phi, dist.type = "euclidean", p = 2)
  w <- knn_weights(w, k, n)
  w
}

find_clustseq <- function(Useq, k){
  nseq = length(Useq)

  for(i in 1:nseq) {
    U <- round(Useq[[i]], 3)
    clust.num <- nrow(unique(U))
    if(clust.num <= k) break
  }
  #cat("i is", i, "\n")
  clust.lab <- rep(1, nrow(U))
  unique.clust <- unique(U)
  nclust <- nrow(unique.clust)
  for(i in 1:nclust){
    is.cur.clust <- sapply(1:nrow(U), function(j) identical(U[j,],unique.clust[i,]))
    clust.lab[is.cur.clust] <- i
  }
  clust.lab
}

#===============================================================================
# Generate fourspherical data
#===============================================================================
generate_spcc_datalist <- function(){
  rep = 10
  n_nums = c(100, 200, 300, 400, 500)
  gdata_list <- list()

  for(n in n_nums)for(i in 1:rep){
    gdata.cur = tgcc:::make.fourspherical(n)
    gdata_list[[paste(n)]][[paste(i)]] = gdata.cur
  }
  gdata_list
}

#===============================================================================
# Perform evaluation for bitgcc
#===============================================================================
evaluate_sptgcc <- function (gdata_list) {
  # Note, we have to set the maxIter big enough
  # Because the spcc method requires a large iteration number
  max_iter <- 10000

  # parameters
  gamma <- 1
  gamseq <- rep(gamma, 20)
  thres <- 1e-05
  rep <- 10

  # evaluation
  res.spcc <- res.sptgcc <- list()
  ARI.mat <- AC.mat <- Time.mat <- array(0, dim = c(2,5,rep))

  for(i in 1:5) for(j in 1:rep){
    cat(i,"_", j,"\n")

    gdata.cur = gdata_list[[i]][[j]]
    data = gdata.cur$data
    label = gdata.cur$label
    # This is a bit wired here (+ 3 +)
    n <<- nrow(data) # do not change
    # To avoid the bug in scvxclustr pacakge
    p = ncol(data)

    lamseq <- seq(1, n, length.out=20)

    # spcc
    tic <- proc.time()
    Gamma2.weight <- rep(1, p)
    w <- gkn_weights_spcc(data)
    nu <- AMA_step_size(w, n) /2
    spccFit <-
      scvxclust(
        X = data,
        w = w,
        Gamma1 = lamseq,
        Gamma2 = gamma, # only take a single input
        Gamma2_weight = Gamma2.weight,
        nu = nu,
        max_iter = max_iter,
        tol_abs = thres
      )
    toc <- proc.time()
    templab <- find_clustseq(spccFit$U, 4)
    ac1 <- 1-mclust::classError(templab, label)$errorRate
    ari1 <- mclust::adjustedRandIndex(templab, label)
    t1 <- (toc - tic)[3]

    # sptgcc
    tic <- proc.time()
    tgccFit <- tgcc::spTGCC(
      data = data,
      lambdaSeq = lamseq,
      gammaSeq = gamseq,
      threshold = thres,
      maxIter = max_iter
    )
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

