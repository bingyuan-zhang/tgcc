SC <- function(data, K, gamma, neighsize = NULL, nstart = 100, iter.max = 1000) {
  n = nrow(data)
  #gamma = norm(data, type = "F")^2/(n*n-n)*gamma
  if(is.null(neighsize)){
    # Fully connected graph
    kfunc = kernlab::rbfdot(sigma = 1/gamma)
    kmat <- kernlab::kernelMatrix(kernel=kfunc, x=data)
    kmat = as(kmat, "dgCMatrix")
  } else {
    # k Nearest Neighbor graph
    nn = annoy_search(data, neighsize)
    kmat = diag(1, n)
    kmat = as(kmat, "dgCMatrix")
    for(i in 1:n){
      nd = exp(-nn$dist[i,-1]^2/gamma)
      kmat[i, nn$idx[i,-1]] = nd
      kmat[nn$idx[i,-1], i] = nd
    }
  }
  
  # Perform spectral clustering
  #-------------------
  n <- nrow(data)
  d_vec <- Matrix::colSums(kmat)
  Dhinv = as(diag(1/sqrt(d_vec)), "dgCMatrix")
  #Ng algorithm
  #--------------------
  Ln <- diag(n) - Dhinv%*%kmat%*%Dhinv
  # Fast computation of SC method
  res.eigen <- RSpectra::eigs_sym(Ln, k = K, which = "LM", sigma = 0)
  V.nlap <- (res.eigen$vector[,1:K])
  u_mat <- t(scale(t(V.nlap), center = FALSE))
  sp_ng <- stats::kmeans(u_mat, centers = K, nstart = nstart, iter.max = iter.max)$cluster
  
  sp_ng
}

# ----------------------
# dl = generate_data(model = "GM1", n = 5000)
# data = dl$data
# label = dl$label
# K = 3
# 
# k = ceiling(log(nrow(data)))
# system.time({sp_ng = SC(data, K, gamma = 1, neighsize = k, nstart = 1)})
# mclust::classError(sp_ng, label)$errorRate
