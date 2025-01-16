library(CCMMR)

checkLamCCMM <- function (
  data,
  k = 10,
  gamma = 1,
  stopat = -Inf,
  step_size = 1000,
  min_step = 50,
  maxlam = NULL) {

  lam_vec <- c()
  num_vec <- c()
  tmp_lam <- step_size
  tmp_step <- step_size
  cnt <- 0
  tmp <- TRUE

  if(is.null(maxlam)){
    maxlam = sum(data^2)/2
  }

  W <- CCMMR::sparse_weights(data, k, phi = 1/gamma)

  while (tmp == TRUE && tmp_lam < maxlam) {
    cnt <- cnt + 1
    ccmm.fit <- CCMMR::convex_clusterpath(data, W, tmp_lam, save_clusterpath = FALSE)

    tmp_num <- ccmm.fit$num_clusters
    lam_vec <- c(lam_vec, tmp_lam)
    num_vec <- c(num_vec, tmp_num)
    if(tmp_num < stopat) break

    if(tmp_num == 1){
      tmp_step <- tmp_step/2
      tmp_lam <- tmp_lam - tmp_step
      if(tmp_step < min_step){
        tmp = FALSE
        break
      }
    }else{
      if( min(num_vec[1:cnt]) == 1){
        tmp_step <- tmp_step/2
      }
      tmp_lam <- tmp_lam + tmp_step
    }
  }

  list(
    lam_max = tmp_lam + tmp_step,
    lam_vec = lam_vec,
    num_vec = num_vec
  )
}




