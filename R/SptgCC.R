#' #' sparse tree-guided convex clustering
#' #'
#' #' @param data matrix with samples as rows and features as columns.
#' #' @param lamseq tuning parameters sequence for rows.
#' #' @param bw bandwidth for the Gaussian kernel between rows.
#' #' @param norm to use normalized distance if true.
#' #' @param dth threshold target of rows with in dth.
#' #' @param prob threshold of edge weight.
#' #' @param max_iter maximum iteration of Dykstra algorithm.
#' #' @param thres the stop condition of Dykstra algorithm.
#' #' @export
#' spTGCC <- function(
#'   data,
#'   lamseq,
#'   gamseq,
#'   bw = NULL,
#'   norm = TRUE,
#'   dth = 50,
#'   prob = 0.1,
#'   max_iter = 500,
#'   thres = 1e-5) {
#'
#'   # find the initialization parameters
#'   init <- init_params(data, bw, norm)
#'   init_params <- init$params
#'
#'   Tp <- init_params$Type
#'   WE <- init_params$WeightE
#'   WV <- init_params$WeightV
#'   Pt <- init_params$Partition
#'   V  <- init_params$Vertice
#'   P  <- init_params$Parent
#'   C  <- init_params$Children
#'
#'   # threshold for outliers
#'   ind <- (Pt < dth)
#'   WE[ind] <- pmax(WE[ind], quantile(WE[ind], prob))
#'
#'   # Optimization
#'   # Coord: Updated Parameter matrix
#'   # Pointer (s: sample, f: feature) : pointers that record the merge events
#'   # Fleft: the nonzero features of original column index with each lambda
#'   Coord <- Pointer <- Fleft <- list()
#'   U <- data
#'   fleft <- 1:ncol(data)
#'
#'   for(k in seq_along(lamseq)) {
#'     n <- nrow(U)
#'     p <- ncol(U)
#'     lam <- lamseq[k]
#'     gam <- gamseq[k]
#'     U0 <- U
#'
#'     # Dynamic Programming
#'     if (n == 1) {
#'       U[,1] <- dp (U, gam, V, Tp, P, WV, WE, C)
#'     } else if (p != 1 & n != 1) {
#'       U <- updateUSC(U, U0, lam, gam, max_iter, thres, V, Tp, P, WV, WE, C)
#'     }
#'
#'     # Record non zero features
#'     select.features <- which(colSums(U) != 0)
#'     U <- U[,select.features]
#'     U0 <- U0[,select.features]
#'
#'     # Clustering step
#'     params <- updatenew(U, U0, V, Tp, Pa, WV, WE, C)
#'     U <- params$newinput
#'     V <- params$newV
#'     P <- params$newP
#'     C <- params$newC
#'     Tp <- params$newT
#'     WE <- params$newWE
#'     WV <- params$newWV
#'
#'     fleft <- fleft[select.features]
#'     Coord[[k]] <- U
#'     Pointer[[k]] <- params$pointer
#'     Fleft[[k]] <- fleft
#'
#'     # ---- break condition ----
#'     if (length(V) == 1) break
#'
#'   }
#'
#'   return(list(
#'     data = data,
#'     coord = Coord,
#'     pointer = Pointer,
#'     fleft = Fleft,
#'     lamseq = lamseq[1:k],
#'     gamseq = gamseq[1:k]
#'   ))
#' }
#'
