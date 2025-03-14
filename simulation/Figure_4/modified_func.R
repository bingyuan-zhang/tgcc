# We modify the CCMMR::convex_clusterpath function by eliminating the post -
#   fitting step, as it is not essential for runtime comparison.

convex_clusterpath_modified <- function(
  X,
  W,
  lambdas,
  tau = 1e-3,
  center = TRUE,
  scale = TRUE,
  eps_conv = 1e-6,
  burnin_iter = 25,
  max_iter_conv = 5000,
  save_clusterpath = FALSE,
  target_losses = NULL,
  save_losses = FALSE,
  save_convergence_norms = FALSE)
{

  # Input checks
  CCMMR:::.check_array(X, 2, "X")
  CCMMR:::.check_weights(W)
  CCMMR:::.check_lambdas(lambdas)
  CCMMR:::.check_scalar(tau, TRUE, "tau", upper_bound = 1)
  CCMMR:::.check_boolean(center, "center")
  CCMMR:::.check_boolean(scale, "scale")
  CCMMR:::.check_scalar(eps_conv, TRUE, "eps_conv", upper_bound = 1)
  CCMMR:::.check_int(burnin_iter, FALSE, "burnin_iter")
  CCMMR:::.check_int(max_iter_conv, FALSE, "max_iter_conv")
  CCMMR:::.check_boolean(save_clusterpath, "save_clusterpath")
  CCMMR:::.check_boolean(save_losses, "save_losses")
  CCMMR:::.check_boolean(save_convergence_norms, "save_convergence_norms")

  # Check the vector of target losses
  if (!is.null(target_losses)) {
    CCMMR:::.check_array(target_losses, 1, "target_losses")

    if (length(lambdas) != length(target_losses)) {
      message = "target_losses should have the same length as lambdas"
      stop(message)
    }

    use_target = TRUE
  } else {
    target_losses = rep(-1, length(lambdas))
    use_target = FALSE
  }

  # Preliminaries
  n = nrow(X)

  # Set the means of each column of X to zero
  if (center) {
    X_ = X - matrix(apply(X, 2, mean), byrow = TRUE, ncol = ncol(X),
      nrow = nrow(X))
  } else {
    X_ = X
  }

  # Transpose X
  X_ = t(X_)

  # Separate the weights into keys and values
  W_idx = t(W$keys) - 1
  W_val = W$values

  # Compute fusion threshold
  eps_fusions = CCMMR:::.fusion_threshold(X_, tau)

  t_start = Sys.time()
  clust = CCMMR:::.convex_clusterpath(X_, W_idx, W_val, lambdas, target_losses,
    eps_conv, eps_fusions, scale, save_clusterpath,
    use_target, save_losses,
    save_convergence_norms, burnin_iter,
    max_iter_conv)
  elapsed_time = difftime(Sys.time(), t_start, units = "secs")

  # Construct result
  result = list()
  result$info = data.frame(
    clust$info_d[1, ],
    clust$info_i[2, ],
    clust$info_d[2, ],
    clust$info_i[1, ]
  )
  names(result$info) = c("lambda", "clusters", "loss", "iterations")

  # Merge table and height vector
  result$merge = t(clust$merge)
  result$height = clust$height

  # # Determine order of the observations for a dendrogram
  # # Start with an entry in a hashmap for each observation
  # D = list()
  # for (i in 1:n) {
  #   D[as.character(-i)] = i
  # }
  #
  # # Work through the merge table to make sure that everything that is
  # # merged is next to each other
  # if (nrow(result$merge) >= 1) {
  #   for (i in 1:nrow(result$merge)) {
  #     D[[as.character(i)]] = c(D[[as.character(result$merge[i, 1])]],
  #       D[[as.character(result$merge[i, 2])]])
  #     D[as.character(result$merge[i, 1])] = NULL
  #     D[as.character(result$merge[i, 2])] = NULL
  #   }
  # }
  #
  # # Finally, create a vector with the order of the observations
  # result$order = c()
  # keys = names(D)
  # for (key in keys) {
  #   result$order = c(result$order, D[[key]])
  # }
  #
  # Add elapsed time
  result$elapsed_time = elapsed_time
  #
  # # Add clusterpath coordinates
  # if (save_clusterpath) {
  #   result$coordinates = t(clust$clusterpath)
  # }
  #
  # # Add lambdas
  # result$lambdas = result$info$lambda
  #
  # # Add fusion threshold
  # result$eps_fusions = eps_fusions
  #
  # # Add vector of possible cluster counts
  # result$num_clusters = unique(result$info$clusters)
  #
  # # Add the number of observations
  # result$n = nrow(X)
  #
  # # Give the result a class
  # class(result) = "cvxclust"
  #
  # # Add losses
  # if (save_losses) {
  #   result$losses = clust$losses
  # }
  #
  # # Add convergence norms
  # if (save_convergence_norms) {
  #   result$convergence_norms = clust$convergence_norms
  # }

  return(result)
}

library(rlang)
#
# X <- data_list[[4]]
# weights = sparse_rbf_kernel_weights(k = "auto",
#   phi = "auto",
#   dist.method = "euclidean",
#   p = 2)
# labels = rownames(X)
# X.center = TRUE
# X.scale = FALSE
# back_track = FALSE
# exact = FALSE
# norm = 2
# t = 1.05
# npcs = min(4L, NCOL(X), NROW(X))
# dendrogram.scale = NULL
# impute_func = function(X) {if(anyNA(X)) missForest(X)$ximp else X}
# status = (interactive() && (clustRviz_logger_level() %in% c("MESSAGE", "WARNING", "ERROR")))
CARP_modified <- function(X,
  weights = sparse_rbf_kernel_weights(k = "auto",
    phi = "auto",
    dist.method = "euclidean",
    p = 2),
  labels = rownames(X),
  X.center = TRUE,
  X.scale = FALSE,
  back_track = FALSE,
  exact = FALSE,
  norm = 2,
  t = 1.05,
  npcs = min(4L, NCOL(X), NROW(X)),
  dendrogram.scale = NULL,
  impute_func = function(X) {if(anyNA(X)) missForest(X)$ximp else X},
  status = (interactive() && (clustRviz_logger_level() %in% c("MESSAGE", "WARNING", "ERROR")))) {

  tic <- Sys.time()

  ####################
  ##
  ## Input validation
  ##
  ####################

  if (!is.matrix(X)) {
    crv_warning(sQuote("X"), " should be a matrix, not a " , class(X)[1],
      ". Converting with as.matrix().")
    X <- as.matrix(X)
  }

  if (!is.numeric(X)) {
    crv_error(sQuote("X"), " must be numeric.")
  }

  # Missing data mask: M_{ij} = 1 means we see X_{ij};
  M <- 1 - is.na(X)

  # Impute missing values in X
  # By default, we use the "Missing Forest" function from the missForest package
  # though other imputers can be supplied by the user.
  X.orig <- X

  if(anyNA(X)) {
    X <- impute_func(X)
  }

  ## Check that imputation was successful.
  if (anyNA(X)) {
    crv_error("Imputation failed. Missing values found in ", sQuote("X"), " even after imputation.")
  }

  if (!all(is.finite(X))) {
    crv_error("All elements of ", sQuote("X"), " must be finite.")
  }

  if (!clustRviz:::is_logical_scalar(X.center)) {
    crv_error(sQuote("X.center"), "must be either ", sQuote("TRUE"), " or ", sQuote("FALSE."))
  }

  if (!clustRviz:::is_logical_scalar(X.scale)) {
    crv_error(sQuote("X.scale"), "must be either ", sQuote("TRUE"), " or ", sQuote("FALSE."))
  }

  if (!clustRviz:::is_logical_scalar(back_track)) {
    crv_error(sQuote("back_track"), "must be either ", sQuote("TRUE"), " or ", sQuote("FALSE."))
  }

  if (!clustRviz:::is_logical_scalar(exact)) {
    crv_error(sQuote("exact"), "must be either ", sQuote("TRUE"), " or ", sQuote("FALSE."))
  }

  l1 <- (norm == 1)

  if ( (!clustRviz:::is_numeric_scalar(t)) || (t <= 1) ) {
    crv_error(sQuote("t"), " must be a scalar greater than 1.")
  }

  if (!is.null(dendrogram.scale)) {
    if (dendrogram.scale %not.in% c("original", "log")) {
      crv_error("If not NULL, ", sQuote("dendrogram.scale"), " must be either ", sQuote("original"), " or ", sQuote("log."))
    }
  }

  if ( (!clustRviz:::is_integer_scalar(npcs)) || (npcs < 2) || (npcs > NCOL(X)) || (npcs > NROW(X)) ){
    crv_error(sQuote("npcs"), " must be an integer scalar between 2 and ", sQuote("min(dim(X))."))
  }

  ## Get row (observation) labels
  if (is.null(labels)) {
    labels <- paste0("Obs", seq_len(NROW(X)))
  }

  if ( length(labels) != NROW(X) ){
    crv_error(sQuote("labels"), " must be of length ", sQuote("NROW(X)."))
  }

  rownames(X.orig) <- rownames(X) <- labels <- make.unique(as.character(labels), sep="_")

  n <- NROW(X)
  p <- NCOL(X)

  # Center and scale X
  if (X.center | X.scale) {
    X <- scale(X, center = X.center, scale = X.scale)
  }

  scale_vector  <- attr(X, "scaled:scale", exact=TRUE)  %||% rep(1, p)
  center_vector <- attr(X, "scaled:center", exact=TRUE) %||% rep(0, p)

  clustRviz:::crv_message("Pre-computing weights and edge sets")

  # Calculate clustering weights
  if (is.function(weights)) {
    # Usual case, `weights` is a function which calculates the weight matrix
    weight_result <- weights(X)

    if (is.matrix(weight_result)) {
      weight_matrix <- weight_result
      weight_type   <- UserFunction()
    } else {
      weight_matrix <- weight_result$weight_mat
      weight_type   <- weight_result$type
    }
  } else if (is.matrix(weights)) {

    if (!is_square(weights)) {
      crv_error(sQuote("weights"), " must be a square matrix.")
    }

    if (NROW(weights) != NROW(X)) {
      crv_error(sQuote("NROW(weights)"), " must be equal to ", sQuote("NROW(X)."))
    }

    weight_matrix <- weights
    weight_type   <- UserMatrix()
  } else {
    crv_error(sQuote("CARP"), " does not know how to handle ", sQuote("weights"),
      " of class ", class(weights)[1], ".")
  }

  if (any(weight_matrix < 0) || anyNA(weight_matrix)) {
    crv_error("All fusion weights must be positive or zero.")
  }

  if (!clustRviz:::is_connected_adj_mat(weight_matrix != 0)) {
    crv_error("Weights do not imply a connected graph. Clustering will not succeed.")
  }

  weight_matrix_ut <- weight_matrix * upper.tri(weight_matrix);

  edge_list <- which(weight_matrix_ut != 0, arr.ind = TRUE)
  edge_list <- edge_list[order(edge_list[, 1], edge_list[, 2]), ]
  cardE <- NROW(edge_list)
  D <- matrix(0, ncol = n, nrow = cardE)
  D[cbind(seq_len(cardE), edge_list[,1])] <-  1
  D[cbind(seq_len(cardE), edge_list[,2])] <- -1

  weight_vec <- clustRviz:::weight_mat_to_vec(weight_matrix)

  clustRviz:::crv_message("Computing Convex Clustering [CARP] Path")
  tic_inner <- proc.time()

  carp.sol.path <- clustRviz:::CARPcpp(X = X,
    M = M,
    D = D,
    t = t,
    epsilon = clustRviz:::.clustRvizOptionsEnv[["epsilon"]],
    weights = weight_vec[weight_vec != 0],
    rho = clustRviz:::.clustRvizOptionsEnv[["rho"]],
    thresh = clustRviz:::.clustRvizOptionsEnv[["stopping_threshold"]],
    max_iter = clustRviz:::.clustRvizOptionsEnv[["max_iter"]],
    max_inner_iter = clustRviz:::.clustRvizOptionsEnv[["max_inner_iter"]],
    burn_in = clustRviz:::.clustRvizOptionsEnv[["burn_in"]],
    viz_max_inner_iter = clustRviz:::.clustRvizOptionsEnv[["viz_max_inner_iter"]],
    viz_initial_step = clustRviz:::.clustRvizOptionsEnv[["viz_initial_step"]],
    viz_small_step = clustRviz:::.clustRvizOptionsEnv[["viz_small_step"]],
    keep = clustRviz:::.clustRvizOptionsEnv[["keep"]],
    l1 = l1,
    show_progress = status,
    back_track = back_track,
    exact = exact)

  toc_inner <- proc.time()

  ## FIXME - Convert gamma.path to a single column matrix instead of a vector
  ##         RcppArmadillo returns a arma::vec as a n-by-1 matrix
  ##         RcppEigen returns an Eigen::VectorXd as a n-length vector
  ##         Something downstream cares about the difference, so just change
  ##         the type here for now
  carp.sol.path$gamma_path <- matrix(carp.sol.path$gamma_path, ncol=1)

  # crv_message("Post-processing")
  #
  # post_processing_results <- ConvexClusteringPostProcess(X = X,
  #   edge_matrix      = edge_list,
  #   gamma_path       = carp.sol.path$gamma_path,
  #   u_path           = carp.sol.path$u_path,
  #   v_path           = carp.sol.path$v_path,
  #   v_zero_indices   = carp.sol.path$v_zero_inds,
  #   labels           = labels,
  #   dendrogram_scale = dendrogram.scale,
  #   npcs             = npcs,
  #   smooth_U         = TRUE)

  carp.fit <- list(
    X = X.orig,
    M = M,
    D = D,
    # U = post_processing_results$U,
    # dendrogram = post_processing_results$dendrogram,
    # rotation_matrix = post_processing_results$rotation_matrix,
    # cluster_membership = post_processing_results$membership_info,
    # n = n,
    # p = p,
    # weights = weight_matrix,
    # weight_type = weight_type,
    # back_track = back_track,
    # exact = exact,
    # norm = norm,
    # t = t,
    # X.center = X.center,
    # center_vector = center_vector,
    # X.scale = X.scale,
    # scale_vector = scale_vector,
    time = Sys.time() - tic,
    fit_time = (toc_inner - tic_inner)[3]
  )

  # if (.clustRvizOptionsEnv[["keep_debug_info"]]) {
  #   carp.fit[["debug"]] <- list(path = carp.sol.path,
  #     row  = post_processing_results$debug)
  # }
  #
  # class(carp.fit) <- "CARP"

  return(carp.fit)
}
