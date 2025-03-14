spRes <- readRDS("./simulation/Table2-3/spRes.rds")
biRes <- readRDS("./simulation/Table2-3/biRes.rds")

biDF <- data.frame()
for(i in 1:5) {
  cobra_df <- data.frame(
    "method" = "COBRA",
    "AC_mean" = mean(biRes$AC[1,i,]),
    "AC_sd" = sd(biRes$AC[1,i,]),
    "ARI_mean" = mean(biRes$ARI[1,i,]),
    "ARI_sd" = sd(biRes$ARI[1,i,]),
    "rt_mean" = mean(biRes$Time[1,i,]),
    "rt_sd" = sd(biRes$Time[1,i,])
  )
  bitgcc_df <- data.frame(
    "method" = "biTGCC",
    "AC_mean" = mean(biRes$AC[2,i,]),
    "AC_sd" = sd(biRes$AC[2,i,]),
    "ARI_mean" = mean(biRes$ARI[2,i,]),
    "ARI_sd" = sd(biRes$ARI[2,i,]),
    "rt_mean" = mean(biRes$Time[2,i,]),
    "rt_sd" = sd(biRes$Time[2,i,])
  )
  biDF <- rbind(biDF, cobra_df, bitgcc_df)
}

spDF <- data.frame()
for(i in 1:5) {
  spcc_df <- data.frame(
    "method" = "SPCC",
    "AC_mean" = mean(spRes$AC[1,i,]),
    "AC_sd" = sd(spRes$AC[1,i,]),
    "ARI_mean" = mean(spRes$ARI[1,i,]),
    "ARI_sd" = sd(spRes$ARI[1,i,]),
    "rt_mean" = mean(spRes$Time[1,i,]),
    "rt_sd" = sd(spRes$Time[1,i,])
  )
  sptgcc_df <- data.frame(
    "method" = "spTGCC",
    "AC_mean" = mean(spRes$AC[2,i,]),
    "AC_sd" = sd(spRes$AC[2,i,]),
    "ARI_mean" = mean(spRes$ARI[2,i,]),
    "ARI_sd" = sd(spRes$ARI[2,i,]),
    "rt_mean" = mean(spRes$Time[2,i,]),
    "rt_sd" = sd(spRes$Time[2,i,])
  )
  spDF <- rbind(spDF, spcc_df, sptgcc_df)
}


biDF; spDF
