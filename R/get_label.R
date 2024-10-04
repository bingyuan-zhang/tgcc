#' Obtain cluster labels
#'
#' @description
#' obtain the estimated cluster labels with desired number of clusters
#'
#' @param tgccFit output from `tgcc`, `biTGCC`, and `spTGCC`.
#' @param numClusters the desired number of clusters.
#' @return A vector indicates the estimated clusters of samples.
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

################################################################################
# The rest are not checked
#
# # get the mode in a label sequence
# getmode <- function(v) {
#   uniqv <- unique(v)
#   uniqv[which.max(tabulate(match(v, uniqv)))]
# }
#
# # prepare the parameters in the heatmap.
# # data/label are the input/true labels
# # tgcc.fit.pointer/lamseqnew are obtained from the result of the TGCC function
# # show is the intermediate cluster number
# # k is the given cluster number of TGCC
# pars_heat = function(tgcc.fit, label, show = 100, k = 2, showlog = FALSE){
#
#   data = tgcc.fit$input
#   tgcc.fit.pointer = tgcc.fit$pointer
#   lamseqnew = tgcc.fit$lamseq
#
#   I = 1
#   tmpshow = length(tgcc.fit.pointer[[I]])
#   while (tmpshow >= show) {
#     I = I + 1
#     tmpshow = length(tgcc.fit.pointer[[I]])
#   }
#   if (I != 1) {
#     p = 1:nrow(data)
#     for (i in 1:(I)) {
#       pnew = tgcc.fit.pointer[[i]] + 1
#       p = pnew[p]
#     }
#   } else {
#     p = tgcc.fit.pointer[[1]] + 1
#   }
#   datanew = matrix(nrow = tmpshow, ncol = ncol(data))
#   labelnew = rep(0, tmpshow)
#
#   for (i in 1:length(unique(p))) {
#     if (sum(p == i) != 1) {
#       datanew[i,] = colSums(data[p == i,]) / sum(p == i)
#       labelnew[i] = getmode(label[p == i])
#     } else {
#       datanew[i,] = data[p == i, ]
#       labelnew[i] = label[p == i]
#     }
#   }
#   n = length(tgcc.fit.pointer)
#
#   pars2 = par_dendrogram(tgcc.fit.pointer[(I+1):n], lamseqnew[(I+1):n])
#
#   if(showlog == FALSE) {
#     showheight = pars2$h
#   } else {
#     showheight = log(pars2$h, base = showlog)
#   }
#
#   h.tgcc2 = list(merge = pars2$m, height = showheight, order = pars2$iorder)
#   class(h.tgcc2) = "hclust"
#   dend.tgcc2 = stats::as.dendrogram(h.tgcc2)
#   labelnew = labelnew[pars2$iorder]
#   datanew = datanew[pars2$iorder,]
#   dend.label = dendextend::cutree(dend.tgcc2, k  = k)
#   est.label = dend.label[p]
#
#   return(list(
#     datanew = datanew,
#     dend = dend.tgcc2,
#     rowcolor = labelnew,
#     estlabel = est.label
#   ))
#
# }
