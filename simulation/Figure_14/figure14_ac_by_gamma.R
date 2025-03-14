library(tgcc)
library(CCMMR)
library(clustRviz)

# Generate data
lab_true <- c(3,3,2,2)
n_rep <- 5
data_list <- list()
set.seed(2025)
for(r in 1:n_rep){
  data_list[["gm1"]][[r]] = tgcc:::make.mixgaussian(n=400)
  data_list[["gm2"]][[r]] = tgcc:::make.threeclust(n=400)
  data_list[["tm"]][[r]] = tgcc:::make.twomoons(n=400)
  data_list[["tc"]][[r]] = tgcc:::make.twocircles(n=400)
}

# Calculate AC with different gamma
gammas <- c(1,2,5,10,20,50,100)
n_gam <- 7
n_mod <- 4

df_tgcc <- data.frame()
for(i in 1:n_mod){
  ac <- matrix(0, nrow = n_rep, ncol = n_gam)
  for(j in 1:n_rep){

    dl <- data_list[[i]][[j]]
    data <- scale(dl$data)
    label <- dl$label
    lamseq <- seq(1, 2000, length.out = 10000)
    for(k in 1:n_gam) {
      print(c(i, j, k))
      tgcc.fit <- tgcc::tgCC(data, lambdaSeq = lamseq, bandwidth = gammas[k])
      estlabel <- tgcc::clusterLabel(tgcc.fit, numClusters = lab_true[i])
      ac[j,k] <- 1 - mclust::classError(estlabel, label)$errorRate
    }
  }

  ac_median = rep(0, n_gam);
  for(k in 1:n_gam) ac_median[k] = median(ac[,k])

  df_tgcc <- rbind(
    df_tgcc,
    data.frame(
      gammas = as.factor(gammas),
      ac = ac_median,
      method = rep("TGCC", n_gam),
      model = rep(i, n_gam)
    )
  )
}
saveRDS(df_tgcc, file = "./simulation/Figure_14/tgcc_gamma.rds")


df_ccmm <- data.frame()
for(i in 1:n_mod){
  ac <- matrix(0, nrow = n_rep, ncol = n_gam)
  for(j in 1:n_rep){
    dl <- data_list[[i]][[j]]
    data <- scale(dl$data)
    label <- dl$label
    for(k in 1:n_gam) {
      print(c(i, j, k))
      W <- CCMMR::sparse_weights(data, k = 10, phi = 1/gammas[k])
      ccmmr.fit <- convex_clustering(data, W, target_low = lab_true[i], target_high = 10*lab_true[i])
      estlabel <- clusters(ccmmr.fit, ccmmr.fit$num_clusters[length(ccmmr.fit$num_clusters)])
      ac[j,k] <- 1 - mclust::classError(estlabel, label)$errorRate
    }
  }

  ac_median = rep(0, n_gam);
  for(l in 1:n_gam) ac_median[l] = median(ac[,l])

  df_ccmm <- rbind(
    df_ccmm,
    data.frame(
      gammas = as.factor(gammas),
      ac = ac_median,
      method = rep("CCMM", n_gam),
      model = rep(i, n_gam)
    )
  )
}
saveRDS(df_ccmm, file = "./simulation/Figure_14/ccmm_gamma.rds")

df_carp <- data.frame()
for(i in 1:n_mod){
  ac <- matrix(0, nrow = n_rep, ncol = n_gam)

  for(j in 1:n_rep){
    dl <- data_list[[i]][[j]]
    data <- scale(dl$data)
    label <- dl$label
    for(k in 1:n_gam) {
      print(c(i, j, k))

      tryCatch({
        weightfunc <- sparse_rbf_kernel_weights(phi = 1)
        res <- weightfunc(data)
        kap <- mean(-log(res$weight_mat[res$weight_mat!=0]))
        carp.fit <-
          CARP(data, weights = sparse_rbf_kernel_weights(phi = 1/gammas[k]/kap))
        estlabel <- cutree(carp.fit$dendrogram, k = lab_true[i])

        1 - mclust::classError(estlabel, label)$errorRate
        ac[j,k] <- 1 - mclust::classError(estlabel, label)$errorRate
      },
      error = function(e){ message("pass the current example for CARP")},
      silent = TRUE)
    }
  }

  ac_median = rep(0, n_gam);
  for(k in 1:n_gam) ac_median[k] = median(ac[,k])

  df_carp <- rbind(
    df_carp,
    data.frame(
      gammas = as.factor(gammas),
      ac = ac_median,
      method = rep("CARP", n_gam),
      model = rep(i, n_gam)
    )
  )
}
saveRDS(df_carp, file = "./simulation/Figure_14/carp_gamma.rds")

color <- c("royalblue3", "orange", "seagreen3")
df <- rbind(df_tgcc, df_ccmm,  df_carp)
line <- 1:3

df_subset <- df[df$model == 4,]
plot(x = NULL, y = NULL,
  xlim = c(1,7), ylim = c(0,1),
  ylab = "", xlab = "",
  cex.axis = 1.5, cex.lab = 1.5, cex.main = 1.5,
  xaxt = "n")
axis(1, at = 1:7, labels = gammas, cex.axis = 1.5)
methods <- c("TGCC", "CCMM", "CARP")
for (i in seq_along(methods)) {
  cur_df <- df_subset[df_subset$method == methods[i],]
  lines(cur_df$gammas, cur_df$ac, col = color[i], lty = line[i], lwd = 2)
}
legend(
  "bottomright",
  legend = methods,
  col = color[1:length(methods)],
  lty = line[1:length(methods)],
  lwd = 2, cex = 2, bty = "n"
)
