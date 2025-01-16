library(tgcc)

nsamples <- c(1000, 5000)
n_len <- length(nsamples)
len_lam <- c(100, 500)
n_rep <- 100
m <- length(len_lam)
time_mst <- rep(0,n_len)
lam_max_mat <- matrix(0,n_rep, n_len)

list_loss = list(array(NA, dim = c(n_len, n_rep, 100)),  array(NA, dim = c(n_len, n_rep, 500)), array(NA, dim = c(n_len, n_rep, 1000)) )
list_loss_tgcc = list(array(NA, dim = c(n_len, n_rep, 100)),  array(NA, dim = c(n_len, n_rep, 500)), array(NA, dim = c(n_len, n_rep, 1000)) )
list_rdiff = list(array(NA, dim = c(n_len, n_rep, 100)),  array(NA, dim = c(n_len, n_rep, 500)), array(NA, dim = c(n_len, n_rep, 1000)) )
list_nclus = list(array(NA, dim = c(n_len, n_rep, 100)),  array(NA, dim = c(n_len, n_rep, 500)), array(NA, dim = c(n_len, n_rep, 1000)) )
list_nclus_tgcc = list(array(NA, dim = c(n_len, n_rep, 100)),  array(NA, dim = c(n_len, n_rep, 500)), array(NA, dim = c(n_len, n_rep, 1000)) )
list_lamsq = list(array(NA, dim = c(n_len, n_rep, 100)),  array(NA, dim = c(n_len, n_rep, 500)), array(NA, dim = c(n_len, n_rep, 1000)) )
list_acc = list(array(NA, dim = c(n_len, n_rep, 100)),  array(NA, dim = c(n_len, n_rep, 500)), array(NA, dim = c(n_len, n_rep, 1000)) )
gamma <- 0.05
set.seed(1)
#for(i in 1:n_len){
i <- 1
n <- nsamples[i]
print(paste("Start: n =", n))
for(r in 1:n_rep){
  data = dataset <- make.sinusoidals(m = n/2, n = n/2, diff = - 0.7, noise = 0.5)$data
  #generate_data("twomoons", n = 1000)$data
  #init = init_pars(data, gamma = gamma)
  #res_lam <- check_lam(data = data, step_size = 1000, ini_len = 200 + 10000/2, min_step = 10, init = init)
  #lam_max_mat[r, i] <- res_lam$lam_max
  for(j in 1:m){
    if(i == 1){
      lamseq = seq(50, 1000, length.out = len_lam[j])
    }else{
      lamseq = seq(500, 3000, length.out = len_lam[j])
    }

    res.fit = tgccloss(data, gamma = gamma, lamseq, delta = NULL, prob = 0.01)
    list_loss[[j]][i,r,] <- res.fit$ltrue
    list_loss_tgcc[[j]][i,r,] <- res.fit$ltgcc
    list_rdiff[[j]][i,r,] <- (res.fit$ltgcc - res.fit$ltrue)/res.fit$ltrue
    list_nclus[[j]][i,r,] <- res.fit$clustnum
    list_lamsq[[j]][i,r,] <- lamseq
    tmp <- length(res.fit$Lab2)
    for(t in 1:len_lam[j]){
      if(t > tmp){
        res.fit$Lab2[[t]] <- rep(1, n)
      }
      list_nclus[[j]][i,r,t] <- max((res.fit$Lab[[t]]))
      list_nclus_tgcc[[j]][i,r,t] <- max((res.fit$Lab2[[t]]))
      list_acc[[j]][i,r,t] <- 1-mclust::classError(classification = res.fit$Lab2[[t]], class = res.fit$Lab[[t]])$errorRate
    }
  }
}
#}

res_all <- list(rdiff = list_rdiff, nclus = list_nclus, nclus_tgcc = list_nclus_tgcc, lamsq = list_lamsq, acc = list_acc)
saveRDS(res_all, "effect_merge_clip.rds")

#res_all <- readRDS("sample_effect_merge.rds")


par(mfrow = c(3, 1))
#n = 1000
lam_max <- max(res_all$lamsq[[1]][1,,])
col_vec <- rainbow(n_rep)
plot(NULL, type = "l", xlim = c(0,1000), ylim = c(0,1), col = col_vec[1], main = "T = 100")
for(r in 1:n_rep){
  lines(x = res_all$lamsq[[1]][1,r,], y = 100*res_all$rdiff[[1]][1,r,], type = "l", col = col_vec[r], lty = 1)
}

#n = 1000
lam_max <- max(res_all$lamsq[[1]][1,,])
col_vec <- rainbow(n_rep)
plot(NULL, type = "l", xlim = c(0,1000), ylim = c(0,1), col = col_vec[1], main = "T = 100")
for(r in 1:n_rep){
  #tmp <- min(which(res_all$nclus[[1]][1,r,]==1)) - 1
  lines(x = res_all$lamsq[[1]][1,r,], y = res_all$acc[[1]][1,r,], type = "l", col = col_vec[r], lty = 1)
}

#n = 1000
lam_max <- max(res_all$lamsq[[1]][1,,])
col_vec <- rainbow(n_rep)
plot(NULL, type = "l", xlim = c(0,1000), ylim = c(-5,5), col = col_vec[1], main = "T = 100")
for(r in 1:n_rep){
  #tmp <- min(which(res_all$nclus[[1]][1,r,]==1)) - 1
  lines(x = res_all$lamsq[[1]][1,r,], y = res_all$nclus[[1]][1,r,]-res_all$nclus_tgcc[[1]][1,r,], type = "l", col = col_vec[r], lty = 1)
  #lines(x = res_all$lamsq[[1]][1,r,], y = res_all$nclus[[1]][1,r,], type = "l", col = col_vec[r], lty = 2)
}

#n = 1000
lam_max <- max(res_all$lamsq[[1]][1,,])
col_vec <- rainbow(n_rep)
plot(NULL, type = "l", xlim = c(0,1000), ylim = c(0,1), col = col_vec[1], main = "T = 100")
for(r in 1:n_rep){
  #tmp <- min(which(res_all$nclus[[1]][1,r,]==1)) - 1
  lines(x = res_all$lamsq[[1]][1,r,], y = res_all$acc[[1]][1,r,], type = "l", col = col_vec[r], lty = 1)
}


lam_max <- max(list_lamsq[[2]][1,,])
col_vec <- rainbow(n_rep)
plot(NULL, type = "l", xlim = c(0,1000), ylim = c(0,1), col = col_vec[1], main = "T = 500")
for(r in 1:n_rep){
  lines(x = list_lamsq[[2]][1,r,], y = 100*list_rdiff[[2]][1,r,], type = "l", col = col_vec[r], lty = 1)
}

lam_max <- max(res_all$lamsq[[2]][1,,])
col_vec <- rainbow(n_rep)
plot(NULL, type = "l", xlim = c(0,1000), ylim = c(0,1), col = col_vec[1], main = "T = 500")
for(r in 1:n_rep){
  #tmp <- min(which(res_all$nclus[[1]][1,r,]==1)) - 1
  lines(x = res_all$lamsq[[2]][1,r,], y = res_all$acc[[2]][1,r,], type = "l", col = col_vec[r], lty = 1)
}


lam_max <- max(list_lamsq[[3]][1,,])
col_vec <- rainbow(n_rep)
plot(NULL, type = "l", xlim = c(0,800), ylim = c(0,1), col = col_vec[1], main = "T = 1000")
for(r in 1:n_rep){
  lines(x = list_lamsq[[3]][1,r,], y = 100*list_rdiff[[3]][1,r,], type = "l", col = col_vec[r], lty = 1)
}

#n = 5000
lam_max <- max(list_lamsq[[1]][2,,])
col_vec <- rainbow(n_rep)
plot(NULL, type = "l", xlim = c(0,lam_max), ylim = c(0,1), col = col_vec[1])
for(r in 1:n_rep){
  for(j in 1:m){
    lines(x = list_lamsq[[j]][2,r,], y = 100*list_rdiff[[j]][2,r,], type = "l", xlim = c(0,lam_max), col = col_vec[r], lty = j)
  }
}



