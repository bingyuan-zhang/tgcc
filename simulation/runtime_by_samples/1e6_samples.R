# Use MixGaussian Example
library(tgcc)
repeatT <- 10
t_lam1_bandNull <- t_lam1_band100 <- rep(0, repeatT)
t_lam100_bandNull <- t_lam100_band100 <- rep(0, repeatT)

################################################################################
set.seed(123)
for(i in 1:repeatT) {
  # Sample generation
  dl <- tgcc:::make.mixgaussian(n = 1e6)
  data <- dl$data
  lamTo <- norm(dl$data, type="2")^2
  lenNum <- 100

  lamFrom <- 1
  lambdaSeq <- seq(lamFrom, lamTo, length.out = lenNum)

  tic <- proc.time()
  tgccFit <- tgcc::tgCC(
    data = data,
    lambdaSeq = lambdaSeq,
    isNaive = FALSE)
  toc <- proc.time()
  t_lam1_bandNull[i] <- (toc - tic)[3]

  tic <- proc.time()
  tgccFit <- tgcc::tgCC(
    data = data,
    lambdaSeq = lambdaSeq,
    bandwidth = 100, # bandwidth = 100
    isNaive = FALSE)
  toc <- proc.time()
  t_lam1_band100[i] <- (toc - tic)[3]

  # different start for lambda

  lamFrom <- 100
  lambdaSeq <- seq(lamFrom, lamTo, length.out = lenNum)

  tic <- proc.time()
  tgccFit <- tgcc::tgCC(
    data = data,
    lambdaSeq = lambdaSeq,
    isNaive = FALSE)
  toc <- proc.time()
  t_lam100_bandNull[i] <- (toc - tic)[3]

  tic <- proc.time()
  tgccFit <- tgcc::tgCC(
    data = data,
    lambdaSeq = lambdaSeq,
    bandwidth = 100, # bandwidth = 100
    isNaive = FALSE)
  toc <- proc.time()
  t_lam100_band100[i] <- (toc - tic)[3]

  print(i)
}

t_lam100_bandNull
t_lam100_band100
t_lam1_band100
t_lam1_bandNull

