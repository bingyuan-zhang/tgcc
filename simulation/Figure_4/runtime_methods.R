library(CCMMR)
library(tgcc)
library(dpcc)
library(clustRviz)
library(dbscan)
library(kernlab)
library(RSpectra)
source("./simulation/standard_clustering_method/Kmeans.R")
source("./simulation/standard_clustering_method/nn.R")
source("./simulation/standard_clustering_method/SC.R")
source("./simulation/Figure_4/modified_func.R")
Rcpp::sourceCpp("./simulation/standard_clustering_method/qpp.cpp")

set.seed(2025)
n_size = c(1e2, 2e2, 5e2, 1e3, 2e3, 5e3, 1e4, 2e4, 5e4, 1e5, 2e5, 5e5, 1e6)
data_list  = list()
for(i in 1:length(n_size)){
  data_list[[i]] = tgcc:::make.mixgaussian(n_size[i])$data
}
LamMax <- 2*n_size
n_rep <- 5

# CARP up to 5000
n_size[6]
to <- 6
rt_carp <- data.frame(
  method = rep("CARP", to),
  sample = n_size[1:to],
  rtime = rep(0, to)
)
for(i in 1:to) {
  data <- data_list[[i]]
  t <- 0
  for(j in 1:n_rep) {
    print(c(i,j))
    w <- sparse_rbf_kernel_weights(data)
    carp.fit <- CARP_modified(data)
    t <- t + carp.fit$fit_time
  }
  rt_carp$rtime[i] <- t/n_rep
}
saveRDS(rt_carp, file = "./simulation/Figure_4/rt_carp.rds")


# SLC up to 50000
n_size[8]
to <- 8
rt_slc <- data.frame(
  method = rep("SLC", to),
  sample = n_size[1:to],
  rtime = rep(0, to)
)
for(i in 1:to) {
  data <- data_list[[i]]
  t <- 0
  for(j in 1:n_rep) {
    print(c(i,j))
    tic <- proc.time()
    res <- hclust(d = dist(data), method = "single")
    toc <- proc.time()
    t <- t + (toc - tic)[3]
  }
  rt_slc$rtime[i] = t/n_rep
}
saveRDS(rt_slc, file = "./simulation/Figure_4/rt_slc.rds")


# CPAINT up to 1e6
n_size[13]
to <- 13
rt_cpaint <- data.frame(
  method = rep("CPAINT", to),
  sample = n_size[1:to],
  rtime = rep(0, to)
)

for(i in 1:to) {
  data <- data_list[[i]]
  t <- 0
  for(j in 1:n_rep) {
    print(c(i,j))
    tic <- proc.time()
    Lam <- seq(1, nrow(data)/2, length.out = 100)
    res <- cpaint(data, Lam)
    toc <- proc.time()
    t <- t + (toc - tic)[3]
  }
  rt_cpaint$rtime[i] = t/n_rep
}
saveRDS(rt_cpaint, file = "./simulation/Figure_4/rt_cpaint.rds")


# CCMM up to 2e5
n_size[11]
to <- 11
rt_ccmm <- data.frame(
  method = rep("CCMM", to),
  sample = n_size[1:to],
  rtime = rep(0, to)
)
rt_knn <- data.frame(
  method = rep("KNN", to),
  sample = n_size[1:to],
  rtime = rep(0, to)
)
for(i in 1:to) {
  data <- data_list[[i]]
  t <- t_knn <- 0
  for(j in 1:n_rep) {
    print(c(i,j))
    k <- ceiling(log(nrow(data)))
    tic_knn <-proc.time()
    W <- sparse_weights(data, k = k,  1)
    toc_knn <-proc.time()
    t_knn <- t_knn + (toc_knn - tic_knn)[3]
    lambdas <- seq(1, LamMax[i], length.out = 100)
    ccmm.fit <- convex_clusterpath_modified(data, W, lambdas, save_clusterpath = FALSE)
    t <- t + ccmm.fit$elapsed_time
  }
  rt_knn$rtime[i] = t_knn/n_rep
  rt_ccmm$rtime[i] = t/n_rep
}
saveRDS(rt_ccmm, file = "./simulation/Figure_4/rt_ccmm.rds")
saveRDS(rt_knn, file = "./simulation/Figure_4/rt_knn.rds")

# TGCC 1e6
n_size[13]
to <- 13
rt_tgcc <- data.frame(
  method = rep("TGCC", to),
  sample = n_size[1:to],
  rtime = rep(0, to)
)
rt_mst <- data.frame(
  method = rep("MST", to),
  sample = n_size[1:to],
  rtime = rep(0, to)
)
for(i in 1:to) {
  data <- data_list[[i]]
  t1 <- t2 <- 0
  for(j in 1:n_rep) {
    print(c(i,j))
    lambdas <- seq(1, LamMax[i], length.out = 100)
    tgcc.fit <- tgcc::tgCC(data, lambdaSeq = lambdas, bandwidth = 10)
    t1 <- t1 + tgcc.fit$tgccTime
    t2 <- t2 + tgcc.fit$mstTime
  }
  rt_tgcc$rtime[i] = t1/n_rep
  rt_mst$rtime[i] = t2/n_rep
}
saveRDS(rt_tgcc, file = "./simulation/Figure_4/rt_tgcc.rds")
saveRDS(rt_mst, file = "./simulation/Figure_4/rt_mst.rds")


# SC upto 2e5
n_size[8]
to <- 8
rt_sc <- data.frame(
  method = rep("SC", to),
  sample = n_size[1:to],
  rtime = rep(0, to)
)
for(i in 1:to) {
  data <- data_list[[i]]
  t <- 0
  for(j in 1:n_rep) {
    print(c(i,j))
    tic <- proc.time()
    sp_ng = SC(
      data,
      K=3,
      gamma = 1,
      neighsize = ceiling(log(nrow(data))),
      nstart = 1
    )
    toc <- proc.time()
    t <- t + (toc-tic)[3]
  }
  rt_sc$rtime[i] = t/n_rep
}
saveRDS(rt_sc, file = "./simulation/Figure_4/rt_sc.rds")

# DBSCAN to 5e5
n_size[12]
to <- 12
rt_dbscan <- data.frame(
  method = rep("DBS", to),
  sample = n_size[1:to],
  rtime = rep(0, to)
)
for(i in 1:to) {
  data <- data_list[[i]]
  t <- 0
  for(j in 1:n_rep) {
    print(c(i,j))
    tic <- proc.time()
    res.dbscan <-
      dbscan::dbscan(
        data,
        eps = 0.5,
        minPts = 5)
    toc <- proc.time()
    t <- t + (toc-tic)[3]
  }
  rt_dbscan$rtime[i] = t/n_rep
}
saveRDS(rt_dbscan, file = "./simulation/Figure_4/rt_dbscan.rds")


# K-means up to 1e6
n_size[13]
to <- 13
rt_km <- data.frame(
  method = rep("KM", to),
  sample = n_size[1:to],
  rtime = rep(0, to)
)
for(i in 1:to) {
  data <- data_list[[i]]
  t <- 0
  for(j in 1:n_rep) {
    print(c(i,j))
    tic <- proc.time()
    res <- kmeans_pp(data, k = 3, nstart = 1, iter.max = 10000)
    toc <- proc.time()
    t <- t + (toc-tic)[3]
  }
  rt_km$rtime[i] = t/n_rep
}
saveRDS(rt_km, file = "./simulation/Figure_4/rt_km.rds")



