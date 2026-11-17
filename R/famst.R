#' Fast approximate minimum spanning tree (FAMST)
#'
#' @description
#' Constructs an approximate minimum spanning tree (MST) of the rows of `data`
#' with the FAMST algorithm of Almansoori and Telek (2025). This is an R/C++
#' re-implementation of the authors' Python code
#' (<https://github.com/mahmood-almansoori/FAMST>), with the approximate
#' k-nearest-neighbour search provided by the `rnndescent` package (an R port
#' of PyNNDescent). Unlike the Python code, the approximate kNN search is
#' always used, also for small or low-dimensional data.
#'
#' The algorithm has four stages:
#' 1. an approximate k-nearest-neighbour (kNN) graph by nearest-neighbour descent;
#' 2. the connected components of the (symmetrised) kNN graph;
#' 3. random bridge edges between every pair of components, refined by
#'    hill climbing over the local neighbourhoods;
#' 4. Kruskal's algorithm on the kNN and bridge edges.
#'
#' @param data matrix with samples as rows and features as columns.
#' @param nNeighbors number of nearest neighbours requested from the kNN search
#'   (the point itself included, so each point keeps `nNeighbors - 1`
#'   neighbours). Default is `15`.
#' @param numRandomEdges number of bridge edges kept between each pair of
#'   components of the kNN graph. Default is `3`.
#' @param metric distance used for the kNN search and the tree edges, one of
#'   `"euclidean"`, `"manhattan"` and `"cosine"`. Default is `"euclidean"`.
#' @param seed optional integer. If given, `set.seed(seed)` is called before
#'   the kNN search, which makes the tree reproducible; note that this
#'   changes the global random number generator state.
#' @return A list with components
#'   \item{mst}{data frame of the tree edges with 0-based endpoints `u`, `v`
#'     and edge length `w`;}
#'   \item{weight}{total length of the tree;}
#'   \item{numComponents}{number of connected components of the result
#'     (1 for a spanning tree);}
#'   \item{nKnnComponents}{number of connected components of the kNN graph;}
#'   \item{timing}{elapsed seconds of the kNN search, the bridge-edge
#'     refinement and the MST extraction.}
#'
#' @references
#' Almansoori, M. K. M. and Telek, M. (2025). FAMST: Fast approximate minimum
#' spanning tree construction for large-scale and high-dimensional data.
#' arXiv:2507.14261.
#'
#' @examples
#' set.seed(1)
#' x <- matrix(rnorm(2000), ncol = 2)
#' tree <- famst(x, seed = 1)
#' head(tree$mst)
#'
#' @export
famst <- function(data,
  nNeighbors = 15,
  numRandomEdges = 3,
  metric = c("euclidean", "manhattan", "cosine"),
  seed = NULL) {

  metric <- match.arg(metric)
  data <- as.matrix(data)
  n <- nrow(data)
  nNeighbors <- as.integer(nNeighbors)
  numRandomEdges <- as.integer(numRandomEdges)
  if (n < 2)
    stop("`data` must have at least two rows.")
  if (nNeighbors < 2 || nNeighbors > n)
    stop("`nNeighbors` must be between 2 and nrow(data).")
  if (numRandomEdges < 1)
    stop("`numRandomEdges` must be at least 1.")
  if (!is.null(seed))
    set.seed(seed)

  # Stage 1: approximate kNN graph.
  tStart <- proc.time()[["elapsed"]]
  knn <- famstKnn(data, nNeighbors, metric)
  knnTime <- proc.time()[["elapsed"]] - tStart

  # Stages 2-3: components of the kNN graph and refined bridge edges.
  tStart <- proc.time()[["elapsed"]]
  neighbors <- matrix(as.integer(knn$neighbors), nrow(knn$neighbors))
  cc <- famstComponents(neighbors, n)
  compId <- cc$comp_id
  components <- famstComponentList(compId, cc$n_components)
  bridge <- famstBridgeEdges(data, components, numRandomEdges, metric)
  bridgeEdges <- bridge$edges
  repeat {
    refined <- famstRefineBridgeEdges(data, bridgeEdges, bridge$comps,
      cc$adj_offsets, cc$adj_indices, compId, metric)
    bridgeEdges <- refined$edges
    if (!refined$changed) break
  }
  refineTime <- proc.time()[["elapsed"]] - tStart

  # Stage 4: Kruskal on the kNN edges and the bridge edges.
  tStart <- proc.time()[["elapsed"]]
  k <- ncol(knn$neighbors)
  uAll <- rep(0:(n - 1), times = k)
  vAll <- as.integer(knn$neighbors)
  wAll <- as.numeric(knn$distances)
  if (length(bridgeEdges) > 0) {
    bm <- do.call(rbind, bridgeEdges)
    uAll <- c(uAll, bm[, 1])
    vAll <- c(vAll, bm[, 2])
    wAll <- c(wAll, bm[, 3])
  }
  km <- famstKruskal(as.integer(uAll), as.integer(vAll), as.numeric(wAll), n)
  mstTime <- proc.time()[["elapsed"]] - tStart

  list(
    mst = data.frame(u = km$u, v = km$v, w = km$w),
    weight = km$weight,
    numComponents = km$num_components,
    nKnnComponents = cc$n_components,
    timing = list(knn = knnTime, refinement = refineTime, mst = mstTime)
  )
}

# Distance between rows u and v (0-based) of `data`.
famstDist <- function(data, u, v, metric) {
  a <- data[u + 1L, ]
  b <- data[v + 1L, ]
  if (metric == "manhattan") {
    sum(abs(a - b))
  } else if (metric == "cosine") {
    den <- sqrt(sum(a * a)) * sqrt(sum(b * b))
    if (den == 0) 0 else 1 - sum(a * b) / den
  } else {
    d <- a - b
    sqrt(sum(d * d))
  }
}

# Approximate kNN graph. rnnd_knn() returns 1-based indices with the point
# itself in the first column; drop it and convert to 0-based indices.
#' @importFrom rnndescent rnnd_knn
famstKnn <- function(data, nNeighbors, metric) {
  res <- rnndescent::rnnd_knn(data, k = nNeighbors, metric = metric,
    verbose = FALSE)
  list(neighbors = res$idx[, -1, drop = FALSE] - 1L,
    distances = res$dist[, -1, drop = FALSE])
}

# List of 0-based node indices of each component.
famstComponentList <- function(compId, nComponents) {
  cc <- split(0:(length(compId) - 1L), compId)[as.character(0:(nComponents - 1L))]
  lapply(cc, as.integer)
}

# Candidate bridge edges: for every pair of components, draw
# numRandomEdges^2 random pairs and keep the numRandomEdges shortest.
famstBridgeEdges <- function(data, components, numRandomEdges, metric) {
  nc <- length(components)
  edges <- list()
  comps <- list()
  idx <- 0L
  if (nc < 2)
    return(list(edges = edges, comps = comps))
  m <- numRandomEdges ^ 2
  for (i in seq_len(nc - 1)) {
    for (j in (i + 1):nc) {
      ci <- sort(components[[i]])
      cj <- sort(components[[j]])
      us <- if (length(ci) == 1) rep(ci, m) else sample(ci, m, replace = TRUE)
      vs <- if (length(cj) == 1) rep(cj, m) else sample(cj, m, replace = TRUE)
      ds <- vapply(seq_len(m),
        function(t) famstDist(data, us[t], vs[t], metric), numeric(1))
      for (t in order(ds)[seq_len(numRandomEdges)]) {
        idx <- idx + 1L
        edges[[idx]] <- c(us[t], vs[t], ds[t])
        comps[[idx]] <- c(i - 1L, j - 1L)
      }
    }
  }
  list(edges = edges, comps = comps)
}

# One hill-climbing pass over the bridge edges: move each endpoint to a
# neighbour in the same component if that shortens the edge.
famstRefineBridgeEdges <- function(data, bridgeEdges, bridgeComps,
  offsets, indices, compId, metric) {
  neighborsOf <- function(u) {
    s <- offsets[u + 1L]
    e <- offsets[u + 2L]
    if (e > s) indices[(s + 1L):e] else integer(0)
  }
  refined <- bridgeEdges
  changed <- FALSE
  for (k in seq_along(bridgeEdges)) {
    e <- bridgeEdges[[k]]
    u <- e[1]; v <- e[2]; d <- e[3]
    ci <- bridgeComps[[k]][1]
    cj <- bridgeComps[[k]][2]
    bestU <- u; bestV <- v; bestD <- d
    for (nu in neighborsOf(u)) {
      if (compId[nu + 1L] == ci && nu != v) {
        dc <- famstDist(data, nu, v, metric)
        if (dc < bestD) { bestD <- dc; bestU <- nu }
      }
    }
    for (nv in neighborsOf(v)) {
      if (compId[nv + 1L] == cj && nv != u) {
        dc <- famstDist(data, bestU, nv, metric)
        if (dc < bestD) { bestD <- dc; bestV <- nv }
      }
    }
    refined[[k]] <- c(bestU, bestV, bestD)
    if (bestU != u || bestV != v) changed <- TRUE
  }
  list(edges = refined, changed = changed)
}
