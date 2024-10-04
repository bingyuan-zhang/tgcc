make.twomoons <- function(n, a = 2, diff = -0.7, noise=0.5){
  m1 = n/2
  m2 = n/2
  x1 <- runif(m1, min = 0, max = pi)
  x2 <- runif(m2, min = pi/2, max = (1+1/2)*pi)
  y1 <- a*sin(x1) + diff/2
  y2 <- a*cos(x2) - diff/2
  data <- matrix(NA, nrow = n, ncol = 2)
  data[,1] <- c(x1, x2) + rnorm(n, sd = noise/2)
  data[,2] <- c(y1, y2) + rnorm(n, sd = noise/2)
  target <- c(rep(1,m1),rep(2,m2))
  return(list(data = data, label = target))
}

make.mixgaussian <- function(n) {
  n1 = round(n/3)
  n2 = n1
  n3 = n-n1-n2
  p = 2
  mean1 = t(matrix(c(1,2.5), nrow = p, ncol = n1))
  mean2 = t(matrix(c(2.5,-1.8), nrow = p, ncol = n2))
  mean3 = t(matrix(c(-2.5,-2), nrow = p, ncol = n3))
  noise1 = matrix(rnorm(n1*p, sd = 1), nrow = n1, ncol = p)
  noise2 = matrix(rnorm(n2*p, sd = 1), nrow = n2, ncol = p)
  noise3 = matrix(rnorm(n3*p, sd = 1), nrow = n3, ncol = p)
  data = rbind(noise1 + mean1, noise2 + mean2, noise3 + mean3)
  label = c(rep(1,n1), rep(2,n2), rep(3,n3))
  return(list(data = data, label = label))
}

make.threeclust <- function (n) {
  n1 = round(n/3)
  n2 = n1
  n3 = n-n1-n2
  mean1 = c(1.3,3.5)
  mean2 = c(2,-2)
  mean3 = c(-1.2,4)
  sigma = matrix(c(1,0.9,0.9,1.2), ncol = 2)
  data = rbind(rmvnorm(n1, mean1, sigma), rmvnorm(n2, mean2, sigma), rmvnorm(n3, mean3, sigma))
  label = c(rep(1,n1), rep(2,n2), rep(3,n3))

  return(list(data = data, label = label))
}

make.twocircles <-function (n) {
  n1 <- n*9/20
  n2 <- n*1/20
  t1 <- 2 * pi * runif(n1)
  t2 <- 2 * pi * runif(n1)
  t3 <- 2 * pi * runif(n2)
  t4 <- 2 * pi * runif(n2)
  y4<-x4<-y3<-x3<-y2<-x2<-y1<-x1<-NULL
  for(t in t1) x1 <- c(x1,cos(t)*runif(1,0.8,0.9))
  for(t in t1) y1 <- c(y1,sin(t)*runif(1,0.8,0.9))
  X1 <- t(as.matrix(rbind(x1, y1)))
  for(t in t2) x2 <- c(x2,cos(t)*runif(1,0.3,0.5))
  for(t in t2) y2 <- c(y2,sin(t)*runif(1,0.3,0.5))
  X2 <- t(as.matrix(rbind(x2, y2)))
  for(t in t3) x3 <- c(x3,cos(t)*runif(1,0.6,0.8))
  for(t in t3) y3 <- c(y3,sin(t)*runif(1,0.6,0.8))
  X3 <- t(as.matrix(rbind(x3, y3)))
  for(t in t4) x4 <- c(x4,cos(t)*runif(1,0.4,0.6))
  for(t in t4) y4 <- c(y4,sin(t)*runif(1,0.4,0.6))
  X4 <- t(as.matrix(rbind(x4, y4)))
  data <- as.matrix(rbind(X1, X2, X3, X4))
  label <-c(rep(1,n1),rep(2,n1),rep(1,n2),rep(2,n2))
  return(list(data = data, label = label))
}

make.sphericalshells <- function (n) {
  t1 <- 2 * pi * runif(n/2)
  t2 <- pi * runif(n/2)
  t3 <- 2 * pi * runif(n/2)
  t4 <- pi * runif(n/2)
  z.shell1 <- sin(t2)
  x.shell1 <- cos(t2)*sin(t1)
  y.shell1 <- cos(t2)*cos(t1)
  z.shell2 <- sin(t4)*1.4
  x.shell2 <- cos(t4)*sin(t3)*1.4
  y.shell2 <- cos(t4)*cos(t3)*1.4
  data = rbind(cbind(x.shell1, y.shell1, z.shell1), cbind(x.shell2, y.shell2, z.shell2))
  label = c(rep(1, n/2), rep(2, n/2))
  return(list(data = data, label = label))
}

make.chainlink <- function (n) {
  t1 <- 2 * pi * runif(n/2)
  x1 <- sin(t1)
  y1 <- cos(t1)
  z1 <- 0
  t2 <- 2 * pi * runif(n/2)
  x2 <- 0
  y2 <- sin(t2) + 1.3
  z2 <- cos(t2)
  data <- rbind(cbind(x1, y1, z1), cbind(x2, y2, z2))
  data <- data + matrix(rnorm(n*3, sd = 0.1), nrow = n)
  label <- c(rep(1, n/2), rep(2, n/2))
  return(list(data = data, label = label))
}

make.atom = function (n) {
  phi <- runif(n/2, max = 2 * pi)
  r <- runif(n/2)^(1/3) * 0.3
  cos_theta <- runif(n/2, min = -1, max = 1)
  x1 <- r * sqrt(1-cos_theta^2) * cos(phi)
  y1 <- r * sqrt(1-cos_theta^2) * sin(phi)
  z1 <- r * cos_theta
  t3 <- 2 * pi * runif(n/2)
  t4 <- 2 * pi * runif(n/2)
  z2 <- sin(t4)*1.5
  x2 <- cos(t4)*sin(t3)*1.5
  y2 <- cos(t4)*cos(t3)*1.5
  data = rbind(cbind(x1, y1, z1), cbind(x2, y2, z2))
  label = c(rep(1, n/2), rep(2, n/2))
  return(list(data = data, label = label))
}

make.tworolls <- function (n) {
  t1 <- runif(n/2,min=1.5,max=5)
  t2 <- runif(n/2,min=0,max=0.5)
  x1 <- sqrt(t1)*cos(2*pi*sqrt(t1))
  y1 <- sqrt(t1)*sin(2*pi*sqrt(t1))
  z1 <- t2
  t3 <- runif(n/2,min=1.5,max=5)
  t4 <- runif(n/2,min=0,max=0.5)
  x2 <- 0.7*sqrt(t3)*cos(2*pi*sqrt(t3))
  y2 <- 0.7*sqrt(t3)*sin(2*pi*sqrt(t3))
  z2 <- t4
  data <- rbind(cbind(x1,y1,z1), cbind(x2,y2,z2))
  data <- matrix(rnorm(n*3, sd = 0.25^2), ncol = 3) + data
  label <- c(rep(1, n/2), rep(2, n/2))
  return(list(data = data, label = label))
}

make.checkboard <- function (n, isGroundTruth = FALSE) {
  # 4 * 4 = 8 check board
  candidates <- seq(from = -8, to = 8, by = 0.7)
  centers <- sample(candidates, 16, replace = FALSE)

  row_cluster_id <-
    sort(sample(
      1:4,
      size = n,
      replace = TRUE,
      prob = exp(1:4 / 4)
    ))
  column_cluster_id <-
    sort(sample(
      1:4,
      size = n,
      replace = TRUE,
      prob = exp(1:4 / 4)
    ))
  sigma <- ifelse(isGroundTruth, 0.01, 3)

  X <- matrix(0, nrow = n, ncol = n)
  for (row_cluster in 1:4) {
    rowid <- row_cluster == row_cluster_id
    for (col_cluster in 1:4) {
      colid <- col_cluster == column_cluster_id
      X[rowid, colid] = rnorm(n = sum(rowid) * sum(colid),
        mean = centers[(row_cluster - 1) * 4 + col_cluster],
        sd = sigma)
    }
  }

  return(list(data = X, label = row_cluster_id))
}

make.fourspherical <- function (n, isGroundTruth = FALSE){
  mu = 2
  # generate informative features from 4 groups
  X <- matrix(0, nrow = n, ncol = 20)
  candidates <- sample(1:4, n, replace = TRUE)
  num_samples = table(candidates)
  sigma = diag(x = 1, nrow=20, ncol=20)
  X[candidates==1] <- mvtnorm::rmvnorm(n = num_samples[1], mean = rep(mu, 20), sigma)
  X[candidates==2] <- mvtnorm::rmvnorm(n = num_samples[2], mean = rep(-mu, 20), sigma)
  X[candidates==3] <- mvtnorm::rmvnorm(n = num_samples[3], mean = c(rep(-mu, 10), rep(mu, 10)), sigma)
  X[candidates==4] <- mvtnorm::rmvnorm(n = num_samples[4], mean = c(rep(mu, 10), rep(-mu, 10)), sigma)

  # generate noises
  X_noise <- matrix(rnorm(n=80*n, mean=0, sd=1), nrow=n)
  X_NA <- matrix(NA, nrow=n, ncol=80)
  if (isGroundTruth) {
    X <- cbind(X, X_NA)
  } else {
    X <- cbind(X, X_noise)
  }
  list(data=X, label=candidates)
}
