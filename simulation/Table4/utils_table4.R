library(dbscan)
library(kernlab)
library(CCMMR)
source("./simulation/standard_clustering_method/Kmeans.R")
source("./simulation/standard_clustering_method/nn.R")
source("./simulation/standard_clustering_method/SC.R")
Rcpp::sourceCpp("./simulation/standard_clustering_method/qpp.cpp")
source("./simulation/Table1/utils.R") # cpaint related function

run_slc <- function(data, label, n_rep = 3) {
  K <- length(unique(label))
  t <- 0
  for(i in 1:n_rep) {
    print(paste("slc", i))
    tic <- proc.time()
    hclust.fit <- hclust(dist(data), method = "single")
    toc <- proc.time()
    t <- t + (toc - tic)[3]
  }
  estLabel <- cutree(hclust.fit, k=K)
  ac <- 1 - mclust::classError(estLabel, label)$errorRate
  ari <- mclust::adjustedRandIndex(estLabel, label)
  t <- t / n_rep
  list(ac = ac, ari = ari, rt = t)
}

run_clc <- function(data, label, n_rep = 3) {
  K <- length(unique(label))
  t <- 0
  for(i in 1:n_rep) {
    print(paste("clc", i))
    tic <- proc.time()
    hclust.fit <- hclust(dist(data), method = "complete")
    toc <- proc.time()
    t <- t + (toc - tic)[3]
  }
  estLabel <- cutree(hclust.fit, k=K)
  ac <- 1 - mclust::classError(estLabel, label)$errorRate
  ari <- mclust::adjustedRandIndex(estLabel, label)
  t <- t / n_rep
  list(ac = ac, ari = ari, rt = t)
}

run_sc <- function(data, label, n_rep = 3) {
  K <- length(unique(label))
  t <- 0
  gamma_sc <- c(0.5,1,2,5,10,20,50)
  AC <- ARI <- rep(0, length(gamma_sc))
  for(i in seq_along(gamma_sc)) {
    for(j in 1:n_rep) {
      print(paste("sc:", "gamma", i, "repeat", j))
      tic <- proc.time()
      sp_ng <- SC(
        data,
        K,
        gamma = gamma_sc[i],
        neighsize = ceiling(log(nrow(data))),
        nstart = 100
      )
      toc <- proc.time()
      t <- t + (toc - tic)[3]
    }
    AC[i] <- 1 - mclust::classError(sp_ng, label)$errorRate
    ARI[i] <- mclust::adjustedRandIndex(sp_ng, label)
  }
  choose <- which.max(AC)
  ac <- AC[choose]
  ari <- ARI[choose]
  t <- t / n_rep / length(gamma_sc)
  list(ac = ac, ari = ari, rt = t)
}

run_km <- function(data, label, n_rep = 3) {
  K <- length(unique(label))
  t <- 0
  for(i in 1:n_rep) {
    print(paste("KM", i))
    tic <- proc.time()
    res.kmeans <- kmeans_pp(data, k = K, nstart = 100, iter.max = 10000)
    toc <- proc.time()
    t <- t + (toc - tic)[3]
  }
  ac <- 1 - mclust::classError(res.kmeans$cluster, label)$errorRate
  ari <- mclust::adjustedRandIndex(res.kmeans$cluster, label)
  t <- t / n_rep
  list(ac = ac, ari = ari, rt = t)
}

run_dbscan <- function(data, label, n_rep = 3) {
  K <- length(unique(label))
  t <- 0
  eps_scale <- c(0.2,0.4,0.6,0.8,1,2,4)
  minP <- c(5,10)
  AC <- ARI <- c()
  for(i in seq_along(eps_scale)) {
    for(j in seq_along(minP)) {
      for(k in 1:n_rep) {
        tic <- proc.time()
        print(paste("dbscan:", "eps", i, "minP", j, "rep", k))
        res.dbscan <-
          dbscan::dbscan(
            data,
            eps = eps_scale[i]/mean(dist(data)),
            minPts = minP[j])
        toc <- proc.time()
        t <- t + (toc - tic)[3]
      }
      AC <- c(AC, 1 - mclust::classError(res.dbscan$cluster, label)$errorRate)
      ARI <- c(ARI, mclust::adjustedRandIndex(res.dbscan$cluster, label))
    }
  }
  choose <- which.max(AC)
  ac <- AC[choose]
  ari <- ARI[choose]
  t <- t / n_rep / length(eps_scale) / length(minP)
  list(ac = ac, ari = ari, rt = t)
}

run_carp <- function(data, label, n_rep = 3) {
  K <- length(unique(label))
  t <- 0
  gammas <- c(0.5,1,2,5,10,20,50)
  AC <- ARI <- c()
  for(i in seq_along(gammas)) {
    for(j in 1:n_rep) {
      print(paste("carp:", "gamma", i, "repeat", j))
      tic <- proc.time()
      weightfunc <- clustRviz::sparse_rbf_kernel_weights(phi = 1, k = 10)
      res <- weightfunc(data)
      kap <- mean(-log(res$weight_mat[res$weight_mat!=0]))
      carp.fit <- clustRviz::CARP(data,
        weights = clustRviz::sparse_rbf_kernel_weights(
          phi = 1/gammas[i]/kap, k = 10))
      estlabel <- cutree(carp.fit$dendrogram, k = K)
      toc <- proc.time()
      t <- t + (toc - tic)[3]
    }
    AC <- c(AC, 1-mclust::classError(estlabel, label)$errorRate)
    ARI <- c(ARI, mclust::adjustedRandIndex(estlabel, label))
  }
  choose <- which.max(AC)
  ac <- AC[choose]
  ari <- ARI[choose]
  t <- t / n_rep / length(gammas)
  list(ac = ac, ari = ari, rt = t)
}

run_cpaint <- function(data, label, n_rep = 3) {
  K <- length(unique(label))
  t <- 0
  lam_max <- dpcc::find_lambda(data)
  Lam <- seq(0.01*lam_max, lam_max, length.out = 100)
  for(i in 1:n_rep) {
    print(paste("cpaint", i))
    tic <- proc.time()
    cpaint.fit <- dpcc::cpaint(data, Lam)
    toc <- proc.time()
    t <- t + (toc - tic)[3]
  }
  res <- cpaint_best_clust(cpaint.fit, label)
  ac <- res$ac
  ari <- res$ari
  t <- t / n_rep
  list(ac = ac, ari = ari, rt = t)
}

run_ccmm <- function(data, label, n_rep = 3) {
  K <- length(unique(label))
  t <- 0
  gammas <- c(0.5,1,2,5,10,20,50)
  lambdas <- seq(1, 2*nrow(data), length.out = 1000)
  AC <- ARI <- c()
  for(i in seq_along(gammas)) {
      for(j in 1:n_rep) {
        tic <- proc.time()
        print(paste("ccmm:", "gamma", i, "rep", j))
        W <- CCMMR::sparse_weights(data, k = 10, phi = 1/gammas[i])
        ccmm.fit <- CCMMR::convex_clusterpath(data, W, lambdas, save_clusterpath = FALSE)
        toc <- proc.time()
        t <- t + (toc - tic)[3]
      }
    estLabel <- clusters(ccmm.fit, K)
    AC <- c(AC, 1 - mclust::classError(estLabel, label)$errorRate)
    ARI <- c(ARI, mclust::adjustedRandIndex(estLabel, label))
  }
  choose <- which.max(AC)
  ac <- AC[choose]
  ari <- ARI[choose]
  t <- t / n_rep / length(gammas)
  list(ac = ac, ari = ari, rt = t)
}

run_tgcc <- function(data, label, n_rep = 3) {
  K <- length(unique(label))
  t <- 0
  gammas <- c(1,2,5,10,20,50,100)
  lamseq = seq(1, 2*nrow(data), length.out = 1000)
  AC <- ARI <- c()
  for(i in seq_along(gammas)) {
    for(j in 1:n_rep) {
      tic <- proc.time()
      print(paste("tgcc:", "gamma", i, "rep", j))
      tgcc.fit <-
        tgcc::tgCC(data = data,
          lambdaSeq = lamseq,
          bandwidth = gammas[i])
      toc <- proc.time()
      t <- t + (toc - tic)[3]
    }
    estLabel <- tgcc::clusterLabel(tgcc.fit, numClusters = K)
    AC <- c(AC, 1 - mclust::classError(estLabel, label)$errorRate)
    ARI <- c(ARI, mclust::adjustedRandIndex(estLabel, label))
  }
  choose <- which.max(AC)
  ac <- AC[choose]
  ari <- ARI[choose]
  t <- t / n_rep / length(gammas)
  list(ac = ac, ari = ari, rt = t)
}


# dl <- tgcc:::make.mixgaussian(400)
# data <- dl$data
# label <- dl$label
#
# slc.res <- run_slc(data, label)
# clc.res <- run_clc(data, label)
# sc.res <- run_sc(data, label)
# km.res <- run_km(data, label)
# dbscan.res <- run_dbscan(data, label)
# cpaint.res <- run_cpaint(data, label)
# carp.res <- run_carp(data, label)
# ccmm.res <- run_ccmm(data, label)
# tgcc.res <- run_tgcc(data, label)
