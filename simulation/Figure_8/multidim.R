library(tgcc)

set.seed(2025)
n_dim <- 2:20
n_samp <- 1e5
lamseq <- seq(50, n_samp/2, length.out = 100)
n_rep <- 5
rt_mst <- rt_tgcc <- array(NA, dim = c(length(n_dim), n_rep))

for(i in seq_along(n_dim)) {
  p <- n_dim[i]
  for(j in 1:n_rep){
    print(paste("dimension =", p, "repeat =", j))
    dataset <- tgcc:::make.multidim.threeclust(n_samp, p)
    res.fit <- tgcc::tgCC(dataset$data, lamseq, bandwidth = 10)
    rt_mst[i,j] <- res.fit$mstTime
    rt_tgcc[i,j] <- res.fit$tgccTime
  }
}
saveRDS(rt_mst, file = "./simulation/Figure_8/rt_mst.rds")
saveRDS(rt_tgcc, file = "./simulation/Figure_8/rt_tgcc.rds")


rt_mat <- rbind(rowMeans(rt_tgcc), rowMeans(rt_mst))
barplot(rt_mat,
  beside = FALSE,
  col = c("black", "lightgray"),
  names.arg = 2:20,
  # legend.text = c("TGCC", "MST"),
  args.legend = list(x = "topleft", bty = "n"),
  border = "black")

