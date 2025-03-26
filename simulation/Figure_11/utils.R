getmode <- function(v) {
  uniqv <- unique(v)
  uniqv[which.max(tabulate(match(v, uniqv)))]
}

# prepare the parameters in the heatmap.
# data/label are the input/true labels
# tgcc.fit.pointer/lamseqnew are obtained from the result of the TGCC function
# show is the intermediate cluster number
# k is the given cluster number of TGCC
parsHeatmap <- function(tgcc.fit, label, show = 100, k = 2, showlog = FALSE){
  data <- tgcc.fit$data
  tgcc.fit.pointer <- tgcc.fit$pointer
  lamseqnew <- tgcc.fit$lambdaSeq

  I = 1
  tmpshow = length(tgcc.fit.pointer[[I]])
  while (tmpshow >= show) {
    I = I + 1
    tmpshow = length(tgcc.fit.pointer[[I]])
  }
  if (I != 1) {
    p = 1:nrow(data)
    for (i in 1:(I)) {
      pnew = tgcc.fit.pointer[[i]] + 1
      p = pnew[p]
    }
  } else {
    p = tgcc.fit.pointer[[1]] + 1
  }
  datanew = matrix(nrow = tmpshow, ncol = ncol(data))
  labelnew = rep(0, tmpshow)
  for (i in 1:length(unique(p))) {
    if (sum(p == i) != 1) {
      datanew[i,] = colSums(data[p == i,]) / sum(p == i)
      labelnew[i] = getmode(label[p == i])
    } else {
      datanew[i,] = data[p == i, ]
      labelnew[i] = label[p == i]
    }
  }
  n = length(tgcc.fit.pointer)

  tgcc.fit.new <- list(
    pointer = tgcc.fit.pointer[(I+1):n],
    lambdaSeq = lamseqnew[(I+1):n])

  pars2 <- tgcc:::parDendrogram(tgcc.fit.new)
  if(showlog == FALSE) {
    showheight = pars2$h
  } else {
    showheight = log(pars2$h, base = showlog)
  }

  h.tgcc2 = list(merge = pars2$m, height = showheight, order = pars2$iorder)
  class(h.tgcc2) = "hclust"
  dend.tgcc2 = stats::as.dendrogram(h.tgcc2)
  labelnew = labelnew[pars2$iorder]
  datanew = datanew[pars2$iorder,]
  dend.label = dendextend::cutree(dend.tgcc2, k  = k)
  est.label = dend.label[p]

  list(
    datanew = datanew,
    dend = dend.tgcc2,
    rowcolor = labelnew,
    estlabel = est.label
  )

}
