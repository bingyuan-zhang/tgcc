#' Obtain cluster labels
#'
#' @description
#' obtain the estimated cluster labels with desired number of clusters
#'
#' @param tgccFit output from `tgcc`, `biTGCC`, and `spTGCC`.
#' @param numClusters the desired number of clusters.
#' @return A vector indicates the estimated clusters of samples.
#'
#' @importFrom stats as.dendrogram
#' @importFrom dendextend cutree
#'
#' @export
clusterLabel <- function (tgccFit, numClusters = 2) {
  data <- tgccFit$data
  pointer <- tgccFit$pointer
  lamseq <- tgccFit$lambdaSeq

  n = length(pointer)
  m = length(unique(pointer[[n]]))
  if (m != 1) {
    # force all clusters to be merged.
    pointer[[n+1]] = rep(0, m)
    lamseq = c(lamseq, lamseq[n] * 2)
  }

  if(numClusters <= 50) {
    I = 1
    tmpshow = length(pointer[[I]])
    len_max = length(pointer) - 1
    while (tmpshow >= 100 & I < len_max) {
      I = I + 1
      tmpshow = length(pointer[[I]])
    }
    if (I != 1) {
      p = 1:nrow(data)
      for (i in 1:(I)) {
        pnew = pointer[[i]] + 1
        p = pnew[p]
      }
    } else {
      p = pointer[[1]] + 1
    }
    n = length(pointer)
    tgccObj = list(pointer = pointer[(I+1):n], lambdaSeq = lamseq[(I+1):n])
    pars = parDendrogram(tgccObj)
    hcTgcc = list(merge = pars$m, height = pars$h, order = pars$iorder)
    class(hcTgcc) = "hclust"
    dendTgcc = stats::as.dendrogram(hcTgcc)
    dendLabel = dendextend::cutree(dendTgcc, k  = numClusters)
    estlabel = dendLabel[p]
  } else {
    I = 1
    tmpshow = length(pointer[[I]])
    len_max = length(pointer) - 1
    while (tmpshow > numClusters & I < len_max) {
      I = I + 1
      tmpshow = length(pointer[[I]])
    }
    if (I != 1) {
      p = 1:nrow(data)
      for (i in 1:(I-2)) {
        pnew = pointer[[i]] + 1
        p = pnew[p]
      }
    } else {
      p = pointer[[1]] + 1
    }
    estlabel = p
  }

  estlabel
}

