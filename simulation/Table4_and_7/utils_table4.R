library(dbscan)
library(kernlab)
library(CCMMR)
source("./simulation/standard_clustering_method/Kmeans.R")
source("./simulation/standard_clustering_method/nn.R")
source("./simulation/standard_clustering_method/SC.R")
Rcpp::sourceCpp("./simulation/standard_clustering_method/qpp.cpp")

cpaint_best_clust <- function(cpaint.fit, label){
  n <- nrow(cpaint.fit$dim1)
  ac <- ari <- rep(0, n)
  for(i in 1:n) {
    X <- t(rbind(cpaint.fit$dim1[i,], cpaint.fit$dim2[i, ]))
    dist_mat <- dist(X)
    estlabel <- cutree(hclust(dist_mat, method = 'single'), h = 1e-5)
    ac[i] <- 1- mclust::classError(estlabel, label)$errorRate
    ari[i] <- mclust::adjustedRandIndex(estlabel, label)
  }
  idx <- which.max(ac)
  list(ac = ac[idx], ari = ari[idx])
}

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
  skip <- 0
  AC <- ARI <- c()
  for (i in seq_along(gammas)) {
    for (j in 1:n_rep) {
      tic <- proc.time()
      print(paste("ccmm:", "gamma", i, "rep", j))
      W <- CCMMR::sparse_weights(
        data, k = ceiling(log(nrow(data))), phi = 1/gammas[i])
      ccmmr.fit <- CCMMR::convex_clustering(
        data, W, target_low = K, target_high = 10*K)
      toc <- proc.time()
      t <- t + (toc - tic)[3]
    }
    tryCatch({
      estLabel <- clusters(
        ccmmr.fit, ccmmr.fit$num_clusters[length(ccmmr.fit$num_clusters)])
      AC <- c(AC, 1 - mclust::classError(estLabel, label)$errorRate)
      ARI <- c(ARI, mclust::adjustedRandIndex(estLabel, label))
      },
      error = function(e) {
        skip <<- skip + 1
        message("current case passed")
      },
      slient = TRUE)
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
  lamseq = seq(1, 2*nrow(data), length.out = 100)
  AC <- ARI <- c()
  for(i in seq_along(gammas)) {
    for(j in 1:n_rep) {
      gc()
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


run_methods <- function(
  data,
  label,
  methods = c("slc",
    "clc",
    "sc",
    "km",
    "dbscan",
    "cpaint",
    "carp",
    "ccmm",
    "tgcc")) {

  df <- data.frame()

  if ("slc" %in% methods) {
    res <- run_slc(data, label)
    df <- rbind(df,
      data.frame(
        method = rep("slc", 3),
        index = c("ac", "ari", "rt"),
        val = c(res$ac, res$ari, res$rt)
      ))
  }

  if ("clc" %in% methods) {
    res <- run_clc(data, label)
    df <- rbind(df,
      data.frame(
        method = rep("clc", 3),
        index = c("ac", "ari", "rt"),
        val = c(res$ac, res$ari, res$rt)
      ))
  }

  if ("sc" %in% methods) {
    res <- run_sc(data, label)
    df <- rbind(df,
      data.frame(
        method = rep("sc", 3),
        index = c("ac", "ari", "rt"),
        val = c(res$ac, res$ari, res$rt)
      ))
  }

  if ("km" %in% methods) {
    res <- run_km(data, label)
    df <- rbind(df,
      data.frame(
        method = rep("km", 3),
        index = c("ac", "ari", "rt"),
        val = c(res$ac, res$ari, res$rt)
      ))
  }

  if ("dbscan" %in% methods) {
    res <- run_dbscan(data, label)
    df <- rbind(df,
      data.frame(
        method = rep("dbscan", 3),
        index = c("ac", "ari", "rt"),
        val = c(res$ac, res$ari, res$rt)
      ))
  }

  if ("cpaint" %in% methods) {
    res <- run_cpaint(data, label)
    df <- rbind(df,
      data.frame(
        method = rep("cpaint", 3),
        index = c("ac", "ari", "rt"),
        val = c(res$ac, res$ari, res$rt)
      ))
  }

  if ("carp" %in% methods) {
    res <- run_carp(data, label)
    df <- rbind(df,
      data.frame(
        method = rep("carp", 3),
        index = c("ac", "ari", "rt"),
        val = c(res$ac, res$ari, res$rt)
      ))
  }

  if ("ccmm" %in% methods) {
    res <- run_ccmm(data, label)
    df <- rbind(df,
      data.frame(
        method = rep("ccmm", 3),
        index = c("ac", "ari", "rt"),
        val = c(res$ac, res$ari, res$rt)
      ))
  }

  if ("tgcc" %in% methods) {
    res <- run_tgcc(data, label)
    df <- rbind(df,
      data.frame(
        method = rep("tgcc", 3),
        index = c("ac", "ari", "rt"),
        val = c(res$ac, res$ari, res$rt)
      ))
  }

  return(df)
}

