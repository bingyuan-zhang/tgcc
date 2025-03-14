#K-means++
#-------------------------
kmeans_pp <- function(x, k, iter.max = 10, nstart = 1, ...){
  n <- nrow(x) # number of data points
  centers <- numeric(k) # IDs of centers
  distances <- matrix(numeric(n * (k - 1)), ncol = k - 1) 
  # distances[i, j]: The distance between x[i,] and x[centers[j],]
  res.best <- list(tot.withinss = Inf) # the best result among <nstart> iterations
  for (rep in 1:nstart) {
    ini_centers <- qpp_arma(x, k, r = 2)$centers
    ## Perform k-means with the obtained centers
    res <- stats::kmeans(x, ini_centers, iter.max = iter.max, nstart = 1, ...)
    res$initial.centers <- x[centers, ]
    ## Store the best result
    if (res$tot.withinss < res.best$tot.withinss) {
      res.best <- res
    }
  }
  res.best
}

#-------------------------
# dl = generate_data(model = "GM1", n = 10000)
# data = dl$data
# label = dl$label
# K = 3
# 
# system.time({
#   res.kmeans <- kmeans_pp(data, k = K, nstart = 100, iter.max = 10000)
# })
# 
# system.time({
#   res.kmeans <- kmeans_pp(data, k = K, nstart = 1, iter.max = 10000)
# })
# 
# res.kmeans$tot.withinss
# mclust::classError(res.kmeans$cluster, label)$errorRate

