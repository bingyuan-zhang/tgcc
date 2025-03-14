library(tgcc)
source("./simulation/Figure_5/tgcc_loss.R")

set.seed(2025)
n_samp <- seq(5e2, 5e5, length.out = 20)
n_lamb <- c(20, 50, 100)
len_samp <- length(n_samp)
len_lamb <- length(n_lamb)
n_rep <- 5
rt_naive <- rt_tgcc <- array(NA, dim = c(len_samp, len_lamb, n_rep))

for(i in 1:len_samp) {
  n <- n_samp[i]
  for(j in 1:len_lamb){
    len_l <- n_lamb[j]
    lamseq <- seq(50, n/2, length.out = len_l)
    for(k in 1:n_rep){
      print(paste("i =", i, "j =", j, "k =", k))
      dataset <- tgcc:::make.mixgaussian(n)
      res.fit <- tgcc_compare_with_naive(dataset$data, lamseq, bandwidth = 10)
      rt_naive[i,j,k] <- res.fit$runtimeNaive
      rt_tgcc[i,j,k] <- res.fit$runtimeTGCC
    }
  }
}

saveRDS(rt_naive, file = "./simulation/Figure_5/rt_naive.rds")
saveRDS(rt_tgcc, file = "./simulation/Figure_5/rt_tgcc.rds")

# large sample size for tgcc
set.seed(2025)
n_samp <- seq(5e5, 1e6, length.out = 5)
n_lamb <- c(20, 50, 100)
len_samp <- length(n_samp)
len_lamb <- length(n_lamb)
n_rep <- 5
rt_tgcc_large <- array(NA, dim = c(len_samp, len_lamb, n_rep))

for(i in 1:len_samp) {
  n <- n_samp[i]
  for(j in 1:len_lamb){
    len_l <- n_lamb[j]
    lamseq <- seq(50, n/2, length.out = len_l)
    for(k in 1:n_rep){
      print(paste("i =", i, "j =", j, "k =", k))
      dataset <- tgcc:::make.mixgaussian(n)
      res.fit <- tgcc::tgCC(dataset$data, lamseq, bandwidth = 10)
      rt_tgcc_large[i,j,k] <- res.fit$tgccTime
    }
  }
}
saveRDS(rt_tgcc_large, file = "./simulation/Figure_5/rt_tgcc_large.rds")


# plot out
rt1 <- readRDS("./simulation/Figure_5/rt_naive.rds")
rt2 <- readRDS("./simulation/Figure_5/rt_tgcc.rds")
rt3 <- readRDS("./simulation/Figure_5/rt_tgcc_large.rds")
rt_naive_mean <- rt_tgcc_mean <- matrix(0,20,3)
rt_tgcc_large_mean <- matrix(0,5,3)
for(j in 1:3) {
  for(i in 1:20) {
    rt_naive_mean[i,j] <- mean(rt1[i,j,])
    rt_tgcc_mean[i,j] <- mean(rt2[i,j,])
  }
  for(i in 1:5) {
    rt_tgcc_large_mean[i,j] <- mean(rt3[i,j,])
  }
}
rt_tgcc_mean <- rbind(rt_tgcc_mean, rt_tgcc_large_mean)

# width = 800, height = 600
plot(
  NULL,
  type = "l",
  xlim = c(3, 6),
  ylim = c(0, 200),
  xlab = "",
  ylab = "",
  cex.axis = 1.5,
  xaxt = "n"
)
axis(
  1,
  at = 3:6,
  labels = c(
    expression(10 ^ 3),
    expression(10 ^ 4),
    expression(10 ^ 5),
    expression(10 ^ 6)
  ), cex.axis = 1.5)
for(i in 1:3) {
  lines(
    x = log10(seq(5e2, 5e5, length.out = 20)),
    y = rt_naive_mean[, i],
    type = "l",
    lty = i,
    col = "blue"
  )
}
for(i in 1:3) {
  lines(
    x = log10(c(
      seq(5e2, 5e5, length.out = 20), seq(5e5, 1e6, length.out = 5)
    )),
    y = rt_tgcc_mean[, i],
    type = "l",
    lty = i,
    col = "red"
  )
}
legend(
  "topleft",
  legend = c(
    "Naive with 20 values",
    "Naive with 50 values",
    "Naive with 100 values",
    "TGCC with 20 values",
    "TGCC with 50 values",
    "TGCC with 100 values"
    ),
  lty = c(1:3, 1:3),
  col = c(rep("blue", 3), rep("red", 3)),
  lwd = ,
  border = FALSE,
  bty = "o")
