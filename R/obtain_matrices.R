#' #' obtain the matrices for biTGCC
#' #'
#' #' @param bitgcc.fit
#' #' @export
#' obtain_matrices_bitgcc <- function(bitgcc.fit){
#'   data <- bitgcc.fit$data
#'   coord <- bitgcc.fit$coord
#'   pointer_s <- bitgcc.fit$pointer_s
#'   pointer_f <- bitgcc.fit$pointer_f
#'
#'   Theta <- list()
#'   id_s <- 1:nrow(data)
#'   id_f <- 1:ncol(data)
#'
#'   for(i in seq_along(coord)){
#'     id_s <- pointer_s[[i]][id_s] + 1
#'     id_f <- pointer_f[[i]][id_f] + 1
#'     Theta[[i]] <- coord[[i]][id_s, id_f]
#'   }
#'   Theta[[length(coord) + 1]] <- matrix(mean(data), nrow = nrow(data), ncol = ncol(data))
#'   Theta
#' }
#'
#' #' obtain matrices for spTGCC
#' #' @param sptgcc.fit
#' #' @export
#' obtain_matrices_sptgcc <- function(sptgcc.fit){
#'
#'   data <- sptgcc.fit$data
#'   fleft <- sptgcc.fit$fleft
#'   coord <- sptgcc.fit$coord
#'   pointer <- sptgcc.fit$pointer
#'
#'   Theta <- list()
#'   mat <- matrix(0, nrow(data), ncol(data))
#'
#'   id_s <- 1:nrow(data)
#'
#'   for(i in seq_along(coord)){
#'     id_s <- pointer[[i]][id_s] + 1
#'     id_f <- fleft[[i]]
#'     mat[,id_f] <- coord[[i]][id_s, ]
#'     mat[,-id_f] <- NA
#'     Theta[[i]] <- mat
#'   }
#'   Theta
#' }
#'
#'
#'
#'
