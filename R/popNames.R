#' @title Sub-sub-population names
#' @description ...
#' 
#' @param names names of individuals
#' 
#' @export
sspop.names <- function(names) substring(names, 1, 2)


#' @title Sub-population names
#' @description ...
#' 
#' @param names names of individuals
#' 
#' @export
spop.names <- function(names) {
  sspop <- sspop.names(names)
  revsspop <- strrev(sspop)
  spop <- sapply(1:length(sspop), 
                 function(i, x, y) paste(sort(c(x[i], y[i])), collapse = ""), 
                 x = sspop, y = revsspop)
  spop
}


strrev <- function(x) sapply(lapply(strsplit(x, NULL), rev), paste, collapse = "")
