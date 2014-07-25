


#' @title design matrix parent population
#' @param sspop
#' @param par.list
#' @param hybrides
#' @export
popstr.matrix <- function (names, sspop, par.list) {
  names <- as.character(names)
  sspop <- as.character(sspop)
  parents <- sort(unique(unlist(par.list)))
  out <- matrix(0, length(sspop), length(parents))
  colnames(out) <- parents
  for (i in 1:length(sspop)) {
    if (sspop[i] != "P") {
      inx <- unlist(par.list[unlist(strsplit(sspop[i], split=""))])
    } else {
      inx <- match(names[i], parents)
    }
    out[i, inx] <- 1
  }
  out
}



# #' @title design matrix parent population
# #' @param sspop
# #' @param par.list
# #' @param hybrides
# #' @export
# popstr.matrix <- function (sspop, par.list, hybrides) {
#   if (is.null(sspop)) {
#     stop ("'sspop' is missing")
#   }
#   sspop <- as.character(sspop)
#   if(missing(hybrides)) {
#     if (is.null(par.list)) {
#       stop ("'par.list' is missing")
#     }
#     parents <- sort(unique(unlist(par.list)))
#     out <- matrix(0, length(sspop), length(parents))
#     colnames(out) <- parents
#     for (i in 1:length(sspop)) {
#       inx <- unlist(par.list[unlist(strsplit(sspop[i], split=""))])
#       out[i, inx] <- 1
#     }
#     return(out)
#   } else {
#     out <- matrix(0, length(sspop), length(hybrides))
#     colnames(out) <- hybrides
#     for (i in 1:length(sspop)) {
#       inx <- match(unlist(strsplit(sspop[i], split="")), hybrides)
#       out[i, inx] <- 1
#     }
#     return(out)
#   }
# }
