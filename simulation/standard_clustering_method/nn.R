library(RcppAnnoy)

annoy_search <- function(X, k, search_k = 100 * k, n_trees = 50) {
  nr <- nrow(X)
  nc <- ncol(X)
  ann <- methods::new(RcppAnnoy::AnnoyEuclidean, nc)
  for (i in 1:nr) {
    ann$addItem(i - 1, X[i, ])
  }
  ann$build(n_trees)
  
  nr <- nrow(X)
  idx <- matrix(nrow = nr, ncol = k+1)
  dist <- matrix(nrow = nr, ncol = k+1)
  nstars <- 50
  for (i in 1:nr) {
    res <- ann$getNNsByVectorList(X[i, ], k+1, search_k, TRUE)
    idx[i, ] <- res$item
    dist[i, ] <- res$distance
  }
  list(idx = (idx + 1), dist = dist)
}