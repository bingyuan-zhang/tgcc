library(clustRviz)
library(CCMMR)

model <- c("GM1", "GM2", "TM", "TC")
num_clust <- c(3,3,2,2)
n_rep <- 50
data_list <- list()

generate.model <- function(n, model = "GM1") {
  if(model == "GM1") return(tgcc:::make.mixgaussian(n))
  if(model == "GM2") return(tgcc:::make.threeclust(n))
  if(model == "TM") return(tgcc:::make.twomoons(n))
  if(model == "TC") return(tgcc:::make.twocircles(n))
}
set.seed(2025)
for(i in seq_along(model)){
  rep_list <- list()
  for(j in 1:n_rep){
    rep_list[[j]] = generate.model(n = 400, model = model[i])
  }
  data_list[[ model[i] ]] = rep_list
}

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

# CC methods
df_carp <- df_ccmm <- df_tgcc <- df_cpaint <- data.frame()
gammas_knn = c(0.5,1,2,5,10,20,50)
gammas_tree = c(1,2,5,10,20,50,100)
for(i in seq_along(gammas_knn)){
  AC = matrix(0, nrow = 4, ncol = n_rep)
  ARI = matrix(0, nrow = 4, ncol = n_rep)
  for(j in seq_along(model)){
    for(k in 1:n_rep){
      print(c("i =", i, "j =", j, "k =", k))
      # assign data/label
      dl <- data_list[[ model[j] ]][[k]]
      data <- dl$data
      label <- dl$label

      # CARP method
      tryCatch({
        weightfunc <- clustRviz::sparse_rbf_kernel_weights(phi = 1, k = 10)
        res <- weightfunc(data)
        kap <- mean(-log(res$weight_mat[res$weight_mat!=0]))
        carp.fit <- clustRviz::CARP(data,
          weights = clustRviz::sparse_rbf_kernel_weights(
            phi = 1/gammas_knn[i]/kap,
            k = 10))
        estlabel <- cutree(carp.fit$dendrogram, k = num_clust[j])
        er <- mclust::classError(estlabel, label)$errorRate
        ari <- mclust::adjustedRandIndex(estlabel, label)
        AC[1,k] = 1-er; ARI[1,k] = ari
        },
        error = function(e){
          message("pass the current replicate for CARP")
        },
        silent = TRUE
      )


      # TGCC method
      tryCatch({
        lam <- tgcc:::checkLamTGCC(
          data,
          bandwidth = gammas_tree[i],
          stopat = 3,
          step_size = nrow(data) / 2
        )
        lamseq = seq(1, lam$lam_max, length.out = 10000)
        tgcc.fit <-
          tgcc::tgCC(data = data,
            lambdaSeq = lamseq,
            bandwidth = gammas_tree[i])
        estlabel <- tgcc::clusterLabel(tgcc.fit, numClusters = num_clust[j])
        er <- mclust::classError(estlabel, label)$errorRate
        ari <- mclust::adjustedRandIndex(estlabel, label)
        AC[2,k] = 1-er; ARI[2,k] = ari
      },
        error = function(e){
          message("pass the current replicates for TGCC")
        },
        silent = TRUE
      )

      # CCMM method
      tryCatch({
        W <- CCMMR::sparse_weights(data, k = 10, phi = 1/gammas_knn[i])
        ccmmr.fit <-
          CCMMR::convex_clustering(data,
            W,
            target_low = num_clust[j],
            target_high = 10 * num_clust[j])
        estlabel <- clusters(
          ccmmr.fit, ccmmr.fit$num_clusters[length(ccmmr.fit$num_clusters)])
        er <- mclust::classError(estlabel, label)$errorRate
        ari <- mclust::adjustedRandIndex(estlabel, label)
        AC[3,k] <- 1-er; ARI[3,k] <- ari
        },
        error = function(e){
          message("pass the current replicate for CCMM")
        },
        silent = TRUE
      )

      # CPAINT method
      tryCatch({
        lam_max <- dpcc::find_lambda(data)
        Lam <- seq(0.01*lam_max, lam_max, length.out = 100)
        cpaint.fit <- dpcc::cpaint(data, Lam)
        res <- cpaint_best_clust(cpaint.fit, label)
        AC[4,k] <- res$ac; ARI[4,k] <- res$ari
      },
        error = function(e){
          message("pass the current replicate for CPAINT")
        },
        silent = TRUE
      )
    }

      carp_new <- data.frame(
        method = "CARP",
        AC = median(AC[1,]),
        ARI = median(ARI[1,]),
        gamma = gammas_knn[i],
        model = model[j]
      )

      tgcc_new = data.frame(
        method = "TGCC",
        AC = median(AC[2,]),
        ARI = median(ARI[2,]),
        gamma = gammas_tree[i],
        model = model[j]
      )

      ccmm_new = data.frame(
        method = "CCMM",
        AC = median(AC[3,]),
        ARI = median(ARI[3,]),
        gamma = gammas_knn[i],
        model = model[j]
      )

      cpaint_new = data.frame(
        method = "CPAINT",
        AC = median(AC[4,]),
        ARI = median(ARI[4,]),
        gamma = "Not Applicable",
        model = model[j]
      )

      df_carp <- rbind(df_carp, carp_new)
      df_tgcc <- rbind(df_tgcc, tgcc_new)
      df_ccmm <- rbind(df_ccmm, ccmm_new)
      df_cpaint <- rbind(df_cpaint, cpaint_new)
    }
}

library(dbscan)
library(kernlab)
library(Rspectra)
source("./simulation/standard_clustering_method/Kmeans.R")
source("./simulation/standard_clustering_method/nn.R")
source("./simulation/standard_clustering_method/SC.R")
Rcpp::sourceCpp("./simulation/standard_clustering_method/qpp.cpp")

# standard methods
AC.spectral <- array(0, dim = c(4,7,50))
ARI.spectral <- array(0, dim = c(4,7,50))
AC.dbscan <- array(0, dim = c(4,2,7,50))
ARI.dbscan <- array(0, dim = c(4,2,7,50))
AC.km <- array(0, dim = c(4,50))
ARI.km <- array(0, dim = c(4,50))
for(i in seq_along(model)){
  for(j in 1:n_rep){
    print(c(i, j))
    #i <- 1; j <- 1;
    dl <- data_list[[ model[i] ]][[j]]
    data <- dl$data
    label <- dl$label
    K <- length(unique(label))

    # Spectral
    gamma_sc <- c(0.5,1,2,5,10,20,50)
    for (k in 1:7) {
      print(c("SC", "gamma_choice = ", k))
      tryCatch({
        sp_ng = SC(
          data,
          K,
          gamma = gamma_sc[k],
          neighsize = ceiling(log(nrow(data))),
          nstart = 100
        )
        AC.spectral[i,k,j] <- 1 - mclust::classError(sp_ng, label)$errorRate
        ARI.spectral[i,k,j] <- mclust::adjustedRandIndex(sp_ng, label)
      },
        error = function(e){},
        silent = TRUE
      )
    }

    # DBSCAN
    eps_scale <- c(0.2,0.4,0.6,0.8,1,2,4)
    minP <- c(5,10)
    for(k in 1:7) for(l in 1:2) {
        print(c("DBSCAN", k, l))
        tryCatch({
          res.dbscan <-
            dbscan::dbscan(
              data,
              eps = eps_scale[k]/mean(dist(data)),
              minPts = minP[l])
          AC.dbscan[i,l,k,j] <-
            1 - mclust::classError(res.dbscan$cluster, label)$errorRate
          ARI.dbscan[i,l,k,j] <-
            mclust::adjustedRandIndex(res.dbscan$cluster, label)
        },
          error = function(e){},
          silent = TRUE
        )
    }

    # KM
    set.seed(2025)
    res.kmeans <- kmeans_pp(data, k = K, nstart = 100, iter.max = 10000)
    AC.km[i,j] = 1 - mclust::classError(res.kmeans$cluster, label)$errorRate
    ARI.km[i,j] = mclust::adjustedRandIndex(res.kmeans$cluster, label)
  }
}

# Hierarchical clustering methods
AC.slc <- ARI.slc <- AC.clc <- ARI.clc <- array(0, dim = c(4, 50))
for(i in seq_along(model)){
  for(j in 1:n_rep){
    print(c(i, j))

    dl <- data_list[[ model[i] ]][[j]]
    data <- dl$data
    label <- dl$label
    K <- length(unique(label))

    # SLC
    hclust.fit <- hclust(dist(data), method = "single")
    estLabel <- cutree(hclust.fit, k=K)
    AC.slc[i,j] = 1 - mclust::classError(estLabel, label)$errorRate
    ARI.slc[i,j] = mclust::adjustedRandIndex(estLabel, label)

    # CLC
    hclust.fit <- hclust(dist(data), method = "complete")
    estLabel <- cutree(hclust.fit, k=K)
    AC.clc[i,j] = 1 - mclust::classError(estLabel, label)$errorRate
    ARI.clc[i,j] = mclust::adjustedRandIndex(estLabel, label)
  }
}

# Finish the table
ac_mat <- matrix(0, nrow = 4, ncol = 9)
rownames(ac_mat) <- c("GM1", "GM2", "TM", "TC")
colnames(ac_mat) <-
  c("SLC", "CL", "SC", "DBS", "KM", "CPAINT", "CARP", "CCMM", "TGCC")


ac_mat[,1] <- apply(AC.slc, 1, median)
ac_mat[,2] <- apply(AC.clc, 1, median)
ac_mat[,3] <- lapply(1:4, function(x)
  apply(AC.spectral[x, , ], 1, median) %>% max) %>% unlist
ac_mat[,4] <- lapply(1:4,
  function(x) {
    best <- 0;
    for (i in 1:2) for (j in 1:7) best <- pmax(median(AC.dbscan[x, i, j,]), best)
    best}) %>% unlist
ac_mat[,5] <- apply(AC.km, 1, median)
ac_mat[,6] <- lapply(1:4,
  function(x) max(df_cpaint$AC[df_cpaint$model == model[x]]) ) %>% unlist
ac_mat[,7] <- lapply(1:4,
  function(x) max(df_carp$AC[df_carp$model == model[x]]) ) %>% unlist
ac_mat[,8] <- lapply(1:4,
  function(x) max(df_ccmm$AC[df_ccmm$model == model[x]]) ) %>% unlist
ac_mat[,9] <- lapply(1:4,
  function(x) max(df_tgcc$AC[df_tgcc$model == model[x]]) ) %>% unlist
saveRDS(ac_mat, file = "./simulation/Table1/ac_mat.rds")
