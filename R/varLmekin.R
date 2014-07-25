


#' @title Variace of lmekin
#' @param x model
#' @export
varLmekin <- function (x, ...) {
  random <- x$vcoef
  gname <- names(random)
  nrow <- sapply(random, function(x) if (is.matrix(x)) nrow(x) else length(x))
  maxcol <- max(sapply(random, function(x) if (is.matrix(x)) 1+ncol(x) else 2))
  temp1 <- matrix(NA, nrow=sum(nrow)+1, ncol=maxcol)
  indx <- 0
  for (term in random) {
    if (is.matrix(term)) {
      k <- nrow(term)
      nc <- ncol(term)
      for (j in 1:k) {
        temp1[j+indx, 1] <- sqrt(term[j, j])
        temp1[j+indx, 2] <- term[j, j]
        if (nc > j) {
          indx2 <- (j+1):nc
          temp1[j+indx, 1+indx2] <- term[j, indx2]
        }
      }
    }
    else {
      k <- length(term)
      temp1[1:k+indx, 1] <- sqrt(term)
      temp1[1:k+indx, 2] <- term
    }
    indx <- indx+k
  }
  temp1[, 1] <- temp1[, 1]*x$sigma
  temp1[, 2] <- temp1[, 2]*x$sigma^2
  temp1[nrow(temp1), ] <- c(x$sigma, x$sigma^2)
  colnames(temp1) <- c("SD", "Var")
  rownames(temp1)  <- c(gname, "residuals")
  temp1
}
