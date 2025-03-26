library(clustRviz)
library(CCMMR)
library(dplyr)

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
    rep_list[[j]] = generate.model(n = 5000, model = model[i])
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
df_ccmm <- df_tgcc <- df_cpaint <- data.frame()
gammas_knn = c(0.5,1,2,5,10,20,50)
gammas_tree = c(1,2,5,10,20,50,100)
for(i in seq_along(gammas_knn)){
  AC = matrix(0, nrow = 3, ncol = n_rep)
  ARI = matrix(0, nrow = 3, ncol = n_rep)
  for(j in seq_along(model)){

    for(k in 1:n_rep){
      print(c("i =", i, "j =", j, "k =", k))
      # assign data/label
      dl <- data_list[[ model[j] ]][[k]]
      data <- dl$data
      label <- dl$label

      # TGCC method
      message("compute TGCC")
      tryCatch({
        lamseq = seq(1, nrow(data), length.out = 10000)
        tgcc.fit <-
          tgcc::tgCC(data = data,
            lambdaSeq = lamseq,
            bandwidth = gammas_tree[i])
        estlabel <- tgcc::clusterLabel(tgcc.fit, numClusters = num_clust[j])
        er <- mclust::classError(estlabel, label)$errorRate
        ari <- mclust::adjustedRandIndex(estlabel, label)
        AC[1,k] = 1-er; ARI[1,k] = ari
      },
        error = function(e){
          message("pass the current replicates for TGCC")
        },
        silent = TRUE
      )

      message("compute CCMM")
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
        AC[2,k] <- 1-er; ARI[2,k] <- ari
      },
        error = function(e){
          message("pass the current replicate for CCMM")
        },
        silent = TRUE
      )

      message("compute CPAINT")
      # CPAINT method
      tryCatch({
        lam_max <- dpcc::find_lambda(data)
        Lam <- seq(0.01*lam_max, lam_max, length.out = 100)
        cpaint.fit <- dpcc::cpaint(data, Lam)
        res <- cpaint_best_clust(cpaint.fit, label)
        AC[3,k] <- res$ac; ARI[3,k] <- res$ari
      },
        error = function(e){
          message("pass the current replicate for CPAINT")
        },
        silent = TRUE
      )
    }

    tgcc_new = data.frame(
      method = "TGCC",
      AC = median(AC[1,]),
      ARI = median(ARI[1,]),
      gamma = gammas_tree[i],
      model = model[j]
    )

    ccmm_new = data.frame(
      method = "CCMM",
      AC = median(AC[2,]),
      ARI = median(ARI[2,]),
      gamma = gammas_knn[i],
      model = model[j]
    )

    cpaint_new = data.frame(
      method = "CPAINT",
      AC = median(AC[3,]),
      ARI = median(ARI[3,]),
      gamma = "Not Applicable",
      model = model[j]
    )

    df_tgcc <- rbind(df_tgcc, tgcc_new)
    df_ccmm <- rbind(df_ccmm, ccmm_new)
    df_cpaint <- rbind(df_cpaint, cpaint_new)
  }
}

# Finish the table
ac_mat <- matrix(0, nrow = 4, ncol = 3)
rownames(ac_mat) <- paste(c("GM1", "GM2", "TM", "TC"), "AC")
colnames(ac_mat) <- c("CPAINT", "CCMM", "TGCC")

ac_mat[,1] <- lapply(1:4,
  function(x) max(df_cpaint$AC[df_cpaint$model == model[x]]) ) %>% unlist
ac_mat[,2] <- lapply(1:4,
  function(x) max(df_ccmm$AC[df_ccmm$model == model[x]]) ) %>% unlist
ac_mat[,3] <- lapply(1:4,
  function(x) max(df_tgcc$AC[df_tgcc$model == model[x]]) ) %>% unlist

ari_mat <- matrix(0, nrow = 4, ncol = 3)
rownames(ari_mat) <- paste(c("GM1", "GM2", "TM", "TC"), "ARI")
colnames(ari_mat) <- c("CPAINT", "CCMM", "TGCC")

ari_mat[,1] <- lapply(1:4,
  function(x) max(df_cpaint$ARI[df_cpaint$model == model[x]]) ) %>% unlist
ari_mat[,2] <- lapply(1:4,
  function(x) max(df_ccmm$ARI[df_ccmm$model == model[x]]) ) %>% unlist
ari_mat[,3] <- lapply(1:4,
  function(x) max(df_tgcc$ARI[df_tgcc$model == model[x]]) ) %>% unlist

saveRDS(ac_mat, file = "./simulation/Table5-6/ac_mat_n5000.rds")
saveRDS(ari_mat, file = "./simulation/Table5-6/ari_mat_n5000.rds")
