init_params <- function (
    data,
    bw = NULL,
    norm = TRUE) {

  # prepare the minimum spanning tree for samples.
  t_start <- proc.time()
  res_mst <- mlpack::emst(data, naive = (ncol(data) > 10))$output
  t_end <- proc.time()
  mst_time <- (t_end - t_start)[3]

  # normalize distances between samples.
  if (norm) {
    dist_sq <- res_mst[,3]^2
    res_mst[,3] <- sqrt(dist_sq / mean(dist_sq))
  }

  # use the median heuristic as bandwidth if not provided.
  if(is.null(bw)) bw <- median(res_mst[,3]^2)

  # initialize parameters.
  t_start <- proc.time()
  params <- init_prepare(res_mst[, 1], res_mst[, 2], res_mst[, 3], bw)
  t_end <- proc.time()
  init_time <- (t_end - t_start)[3]

  list(params = params, mst_time = mst_time, init_time = init_time)
}


init_params_bicc <- function (
    data,
    bw_s = NULL,
    bw_f = NULL,
    norm = TRUE) {

  # prepare minimum spanning tree for samples and features.
  res_mst_s <- mlpack::emst(data, naive = (ncol(data) > 10))$output
  res_mst_f <- mlpack::emst(t(data), naive = (nrow(data) > 10))$output

  # normalize distances.
  if(norm == TRUE) {
    dist_sq_s <- res_mst_s[, 3]^2
    res_mst_s[, 3] <- sqrt(dist_sq_s / mean(dist_sq_s))
    dist_sq_f <- res_mst_f[, 3]^2
    res_mst_f[, 3] <- sqrt(dist_sq_f / mean(dist_sq_f))
  }

  # use median heuristic for bandwidth if not provided.
  if (is.null(bw_s)) bw_s <- median(res_mst_s[, 3]^2)
  if (is.null(bw_f)) bw_f <- median(res_mst_f[, 3]^2)

  # initialize parameters for samples and features.
  params_s <- init_prepare(res_mst_s[, 1], res_mst_s[, 2], res_mst_s[, 3], bw_s)
  params_f <- init_prepare(res_mst_f[, 1], res_mst_f[, 2], res_mst_f[, 3], bw_f)

  list(params_s = params_s, params_f = params_f)
}
