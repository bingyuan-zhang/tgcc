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


