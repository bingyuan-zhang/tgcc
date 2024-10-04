# # prepare the dataframe corresponding to the clusterpath
getPathDF <- function (tgccFit) {

  K <- length(tgccFit$coordinates)
  startCood <- tgccFit$data
  n <- nrow(startCood)
  p <- ncol(startCood)
  idx <- 1:n
  df <- data.frame()

  for(k in 1:K){
    endCoord <- matrix(tgccFit$coordinates[[k]], ncol = p)
    idxNew <- tgccFit$pointer[[k]] + 1
    startNew <- data.frame(
      x = startCood[, 1],
      y = startCood[, 2],
      groups = paste0(k,"_", 1:nrow(startCood)),
      steps = k
    )
    endNew <- data.frame(
      x = endCoord[idxNew, 1],
      y = endCoord[idxNew, 2],
      groups = paste0(k,"_", 1:nrow(startCood)),
      steps = k
    )

    df <- rbind(df, rbind(startNew, endNew))
    startCood <- endCoord
  }

  df
}



