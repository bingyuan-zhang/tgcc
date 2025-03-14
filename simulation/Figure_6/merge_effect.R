library(tgcc)
source("./simulation/Figure_6/tgcc_loss.R")

set.seed(2025)
n_samples <- 5000
n_rep <- 100
lamseq <- seq(50, 6000, length.out = 1000)
loss_true <- loss_tgcc <- rdiff <- array(NA, dim = c(n_rep, 1000))

for(r in 1:n_rep){
  print(paste("r =", r))
  dataset <- tgcc:::make.mixgaussian(n_samples)
  res.fit <- tgccLoss(dataset$data, lamseq, bandwidth = 10)
  # max((res.fit$ltgcc - res.fit$ltrue)/res.fit$ltrue)
  loss_true[r,] <- res.fit$ltrue
  loss_tgcc[r,] <- res.fit$ltgcc
  rdiff[r,] <- (res.fit$ltgcc - res.fit$ltrue)/res.fit$ltrue
}

saveRDS(rdiff, "./simulation/Figure_6/effect_merge_rf.rds")
rdiff <- readRDS("./simulation/Figure_6/effect_merge_rf.rds")

# width = 800, height = 600
plot(NULL, type = "l", xlim = c(0, 6000), ylim = c(0, 1), xlab = "", ylab = "", cex.axis = 1.5)
for(r in 1:n_rep){
  lines(x = lamseq, y = 100*rdiff[r,], type = "l", lty = 1)
}


