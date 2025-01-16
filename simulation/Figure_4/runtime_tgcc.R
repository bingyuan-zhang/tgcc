### Find Maximum Lambda for TGCC / CCMM
LamMax_CCMM <- LamMax_TGCC <- rep(0, 10)
for(i in 1:10) {
  print(i)
  data <- data_list[[i]]

  # check the maximum lambda
  res_lam_ccmm <-
    checkLamCCMM(
      data,
      k = ceiling(log(nrow(data))),
      gamma = 1,
      stopat = 3,
      step_size = nrow(data) / 2
    )
  res_lam_tgcc <-
    tgcc:::checkLamTGCC(
      data,
      bandwidth = 10,
      stopat = 3,
      step_size = nrow(data) / 2
    )

  LamMax_CCMM[i] <- res_lam_ccmm$lam_max
  LamMax_TGCC[i] <- res_lam_tgcc$lam_max
}

res <- list(
  max_lambda_ccmm = LamMax_CCMM,
  max_lambda_tgcc = LamMax_TGCC
)

saveRDS(res, file = "./simulation/Figure_4/max_lambda_ccmm_tgcc.rds")
