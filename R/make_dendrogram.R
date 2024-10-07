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

  dendextend::colored_bars(barColors, text_shift = NULL, y_shift = -50)
}


