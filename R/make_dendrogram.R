# prepare parameters of a dendrogram (m, h, iorder)
parDendrogram <- function(tgccFit) {

  pointer <- tgccFit$pointer
  lamseq <- tgccFit$lambdaSeq
  # check if all samples were merged into one final cluster
  # if not, force the final merge and assign an arbitrary lambda value for the dendrogram
  len_pt = length(pointer)
  cluster_num = length(unique(pointer[[len_pt]]))
  if (cluster_num != 1) {
    # message = paste("clusters are unmerged with current lambda sequence")
    # warning(message)
    pointer[[len_pt + 1]] = rep(0, cluster_num)
    lamseq = c(lamseq, lamseq[len_pt] * 1.5)
  }

  n = length(pointer[[1]])
  m = c()
  l = c()
  l.new = c()
  h = c()
  l = -(1:n)
  j = 0

  for(p in 1:length(pointer)){
    cur_pointer = pointer[[p]] + 1 # start from 1 not 0
    K = length(unique(cur_pointer))
    lam = lamseq[p]

    for(k in 1:K){
      ind = which(cur_pointer == k)
      I = l[ind]

      if(length(I) > 1) {

        if(sum(I > 0) == 0) {
          m.new = c(I[1], I[2]) # new cluster
          h = c(h, lam)
          j = j + 1

          if (length(I) > 2 ){
            for (t in 3:length(ind)) {
              newrow = c(I[t], j)
              m.new = rbind(m.new, newrow)
              h = c(h, lam)
              j = j + 1
            }
          }
          l.new[k] = j

        } else {
          to = which(I == max(I))
          tt = I[to]
          m.new = c()
          for (t in I[-to]) {
            newrow = c(t, tt)
            m.new = rbind(m.new, newrow)
            h = c(h, lam)
            tt = j + 1
            j = j + 1
          }
          l.new[k] = j
        }
        m = rbind(m, m.new)
      } else {
        l.new[k] = l[ind]
      }

    }

    l = l.new
    l.new = c()
  }

  iia = m[,1]; iib = m[,2]
  iorder = rep(0, n)
  iorder[1] <- iia[n - 1]
  iorder[2] <- iib[n - 1]
  loc <- 2
  for (i in (n - 2):1){
    for (j in 1:loc){
      if (iorder[j] == i) {
        iorder[j] <- iia[i]
        if (j == loc) {
          loc <- loc + 1
          iorder[loc] <- iib[i]
        } else {
          loc <- loc + 1
          for (k in loc:(j + 2)){
            iorder[k] <- iorder[k - 1]
          }
          iorder[j + 1] <- iib[i]
        }
        break
      }
    }
  }

  list(m = m, h = h, iorder = -iorder)
}


getmode <- function(v) {
  uniqv <- unique(v)
  uniqv[which.max(tabulate(match(v, uniqv)))]
}

# prepare the parameters in the heatmap for very large dataset.
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

  pars2 <- parDendrogram(tgcc.fit.new)
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

# Generate color hues for dendrogram
#' @importFrom grDevices hcl
ggColorHue <- function(numColors, lightness = 65) {
  hues <- seq(15, 375, length=numColors + 1)
  grDevices::hcl(h=hues, l=lightness, c=100)[1:numColors]
}

#' Plot the dendrogram
#'
#' @description
#' Plot the dendrogram based on the TGCC result.
#'
#' @param tgccFit output from `tgcc`, `biTGCC`, and `spTGCC`.
#' @param labels the underlying ground truth cluster labels.
#' @return A dendrogram of hierarchical clustering indicates how samples are merged into clusters.
#'
#' @export
makeDendrogram <- function(tgccFit, labels) {
  params = parDendrogram(tgccFit)
  hcTgcc = list(merge = params$m, height = params$h, order = params$iorder)
  class(hcTgcc) = "hclust"

  # Define color bar
  ggColors <- ggColorHue(length(unique(labels)))
  barColors <- ggColors[labels[params$iorder]]

  dendrogram <- stats::as.dendrogram(hcTgcc)
  dendextend::labels(dendrogram) <- rep("", length(labels(dendrogram)))

  plot(dendrogram)

  dendextend::colored_bars(barColors, text_shift = NULL, y_shift = -5)
}


